"""Leave-observed-out real census test with a physically sanitized fitting file.

Assignment uses declared E0/campaign covariates, never E1 fate/DBH. Evaluation
owns the truth table; fit and predict receive only visible CSV and covariates.
The first origin reconstruction calls BOTH old read_core/fit/reconstruct intact.
The tagged-endpoint projection uses their sampled unknown-fate kernel, retaining
its cohort RMS diameter limitation rather than substituting a better estimator.
"""
import copy
import datetime as dt
import math
from collections import Counter
from pathlib import Path
import numpy as np
from h03_common import (ROOT, RUN, PROTOCOL, frozen, empirical, demographic,
                        rng, write_csv, write_json, TraceRNG)

# Frozen fitting code ignores ancillary AGB/measurement-ID fields, but they
# could reveal a masked diameter to a future consumer. Permit identity/taxon
# metadata only; endpoint date remains a separate registered query covariate.
VISIBLE_IDENTITY_FIELDS={"stem.id","tree.id","tag","stem.tag","sp"}


def secure(old, new):
    if not old or old["df.status"] != "alive" or new["tree.id"] != old["tree.id"] or new["sp"] != old["sp"]:
        return False
    if frozen.fate(new) is None:
        return False
    tokens = set(new["codes"].split(";"))
    if new["df.status"] == "alive":
        return new["status"] == "A"
    return "D" in tokens and not ({"Miss", "MT", "S"} & tokens)


def covariates(e0, e1):
    result = []
    denominator, missing = np.zeros(64), np.zeros(64)
    for key, old in e0.items():
        if old["df.status"] == "alive":
            c = empirical.cell(old)
            denominator[c] += 1
            if frozen.fate(e1.get(key)) is None:
                missing[c] += 1
    fraction = np.divide(missing, denominator, out=np.zeros(64), where=denominator>0)
    for key in sorted(e1, key=int):
        new, old = e1[key], e0.get(key)
        if secure(old, new):
            interval = (dt.date.fromisoformat(new["exact.date"])-dt.date.fromisoformat(old["exact.date"])).days/365.25
            if interval > 0:
                c = empirical.cell(old)
                result.append(dict(stem_id=key, cell=c, years=interval,
                                   campaign_month=int(new["exact.date"][5:7]),
                                   observed_unresolved_fraction=float(fraction[c])))
    return result


def mask_probability(cov, family):
    if set(cov)!={"stem_id","cell","years","campaign_month","observed_unresolved_fraction"}:
        raise ValueError("Mask assignment accepts only registered pre-outcome covariates")
    c = cov["cell"]
    if family == 0:
        return 0.15
    if family == 1:
        return 0.08+0.08*(c//16)+0.035*(c%4)
    if family == 2:
        return [0.10,0.16,0.25,0.55][(c//4)%4]
    return float(np.clip((0.10+0.65*cov["observed_unresolved_fraction"])
                        *(1+0.1*(cov["campaign_month"]%3-1)), 0.05, 0.8))


def fit_visible(sources, path, withheld):
    rows = __import__("h03_common").read_csv(path)
    for row in rows:
        if row["stem.id"] in withheld:
            if any(row[field] != "NA" for field in PROTOCOL["calibration"]["heldout_fields"]):
                raise ValueError("Held-out truth leaked into fit input")
            if row["exact.date"] != "NA":
                raise ValueError("A held-out dated fate remains available")
            if any(value!="NA" for field,value in row.items() if field not in VISIBLE_IDENTITY_FIELDS):
                raise ValueError("Ancillary heldout measurement-derived field leaked")
    # Old parser accepts absolute provenance-bearing paths; no HFD01 write.
    visible_sources = copy.deepcopy(sources)
    visible_sources["HF253-E1"]["local_cache_path"] = str(path)
    core = empirical.read_core(visible_sources)
    return core, empirical.fit(core)


def endpoint_kernel(model,cov,trace):
    cells = np.array([x["cell"] for x in cov])
    years = np.array([x["years"] for x in cov])
    hazard=np.array(trace.trace["choice"][:64]).reshape(64)
    growth=np.clip(trace.trace["normal"][0],0,0.7)
    shifts=trace.trace["normal"][1]+((cells//4)%4==3)*trace.trace["normal"][2]
    p=np.exp(-hazard[cells]*np.exp(shifts)*years)
    d=np.clip(model["diameter0"][cells]+growth[cells]*years+trace.trace["normal"][4][cells],1,300)
    return p,d


def predict(core, model, cov, family, replicate):
    probabilities, sizes, alive_draws, origin = [], [], [], []
    for draw in range(PROTOCOL["calibration"]["draws"]):
        trace = TraceRNG(rng("calibration", family, replicate, draw, 0))
        state = demographic.reconstruct(core, model, trace)
        origin.append([int(state["n"].sum()), float(np.sum(state["n"]*state["d"]**2)*math.pi/40000)])
        # The same kernel is projected to known actual endpoint dates; no true E1
        # size or fate is provided. All same-cell diameters remain cohort RMS.
        p,d=endpoint_kernel(model,cov,trace)
        probabilities.append(p)
        sizes.append(d)
        alive_draws.append(rng("calibration", family, replicate, draw, 1).binomial(1,p))
    return np.array(probabilities), np.array(sizes), np.array(alive_draws), origin


def metrics(mask, replicate, label, truth, probability, size, alive_draw):
    p = probability.mean(axis=0)
    fate = truth["fate"]
    live = (fate == 1) & np.isfinite(truth["dbh"])
    fields = dict(mask="M"+str(mask), replicate=replicate, subgroup=label, n=len(fate),
                  deaths=int((fate==0).sum()), brier=float(np.mean((p-fate)**2)),
                  log_score=float(-np.mean(fate*np.log(np.clip(p,1e-9,1-1e-9))
                                          +(1-fate)*np.log(np.clip(1-p,1e-9,1-1e-9)))))
    if label=="ALL":
        fields.update(taxon_group="ALL",E0_size_class="ALL",stratum="ALL")
    else:
        c=int(label.split("_")[1])
        fields.update(taxon_group=["hemlock","other_canopy","shrub","other"][c//16],
                      E0_size_class=["1-3cm","3-10cm","10-30cm","30+cm"][c%4],
                      stratum=["west","central","east","wet_candidate_proxy"][(c//4)%4])
    estimate = size.mean(axis=0)
    error = estimate[live]-truth["dbh"][live]
    fields.update(dbh_n=int(live.sum()), dbh_bias=float(error.mean()) if len(error) else None,
                  dbh_mae=float(np.abs(error).mean()) if len(error) else None,
                  dbh_rmse=float(np.sqrt((error**2).mean())) if len(error) else None)
    growth_error=error/truth["years"][live]
    fields.update(growth_bias=float(growth_error.mean()) if len(error) else None,
                  growth_mae=float(np.abs(growth_error).mean()) if len(error) else None,
                  growth_rmse=float(np.sqrt((growth_error**2).mean())) if len(error) else None)
    for coverage in PROTOCOL["calibration"]["intervals"]:
        lo, hi = np.quantile(size[:,live],[(1-coverage)/2,(1+coverage)/2],axis=0)
        fields["coverage_"+str(int(coverage*100))] = float(np.mean(
            (truth["dbh"][live]>=lo)&(truth["dbh"][live]<=hi))) if live.any() else None
        fields["width_"+str(int(coverage*100))] = float((hi-lo).mean()) if live.any() else None
        fields["growth_width_"+str(int(coverage*100))] = float(((hi-lo)/truth["years"][live]).mean()) if live.any() else None
    observed_ba = float(np.sum(truth["dbh"][live]**2)*math.pi/40000)
    ba_eligible=np.isfinite(truth["dbh"])|(fate==0)
    # No partial observed BA total is compared to a full predicted total.
    # Missing-alive DBH is excluded ONLY at evaluation, not given to estimator.
    predicted_ba = np.sum(alive_draw[:,ba_eligible]*size[:,ba_eligible]**2,axis=1)*math.pi/40000
    observed_juvenile = int(np.sum((truth["dbh"][live]>=1)&(truth["dbh"][live]<10)))
    predicted_juvenile=(alive_draw[:,ba_eligible]*(size[:,ba_eligible]<10)).sum(axis=1)
    fields.update(living_truth=int(fate.sum()), living_prediction=float(alive_draw.sum(axis=1).mean()),
                  living_error=float(alive_draw.sum(axis=1).mean()-fate.sum()),
                  ba_truth=observed_ba, ba_prediction=float(predicted_ba.mean()),
                  ba_error=float(predicted_ba.mean()-observed_ba),
                  aggregate_size_truth_excluded=int((~ba_eligible).sum()),
                  juvenile_truth=observed_juvenile,
                  juvenile_error=float(predicted_juvenile.mean()-observed_juvenile))
    for coverage in PROTOCOL["calibration"]["intervals"]:
        for name,draws,actual in [("living",alive_draw.sum(axis=1),fate.sum()),
                                 ("ba",predicted_ba,observed_ba),
                                 ("juvenile",predicted_juvenile,observed_juvenile)]:
            lo,hi=np.quantile(draws,[(1-coverage)/2,(1+coverage)/2])
            suffix=str(int(coverage*100))
            fields[name+"_lower_"+suffix]=float(lo)
            fields[name+"_upper_"+suffix]=float(hi)
            fields[name+"_covered_"+suffix]=bool(lo<=actual<=hi)
            fields[name+"_width_"+suffix]=float(hi-lo)
    return fields


def run(sources, e0, e1):
    cfg = PROTOCOL["calibration"]
    pool = covariates(e0,e1)
    RUN.mkdir(parents=True,exist_ok=True)
    pool_rows = []
    for x in pool:
        new = e1[x["stem_id"]]
        coordinate = all(empirical.number(new[field]) is not None for field in ("gx","gy","quadrat"))
        pool_rows.append(dict(stem_id=x["stem_id"], cell=x["cell"],
                              fate_eligible=True, dbh_eligible=frozen.fate(new)==1 and empirical.number(new["dbh"]) is not None,
                              coordinate_eligible=coordinate, provenance="OBSERVED",
                              source_id="HF253-E1", secure_death="D corroboration; Miss/MT/S excluded"))
    write_csv("reconstruction_validation/eligible_truth_pool.csv",pool_rows)
    groups, curves, runs, coordinate_rows, compositions = [], [], [], [], []
    assignments = []
    for family in range(4):
        for replicate in range(cfg["replicates"]):
            selection = rng("mask",family,replicate).random(len(pool)) < np.array([mask_probability(x,family) for x in pool])
            held = [x for x, take in zip(pool,selection) if take]
            ids = {x["stem_id"] for x in held}
            assignments.extend(dict(mask="M"+str(family),replicate=replicate,stem_id=x["stem_id"]) for x in held)
            visible = []
            for key,row in e1.items():
                row = dict(row)
                if key in ids:
                    row={field:value if field in VISIBLE_IDENTITY_FIELDS else "NA" for field,value in row.items()}
                visible.append(row)
            path = RUN / ("visible_M%d_%d.csv" % (family,replicate))
            write_csv(str(path),visible)
            core, model = fit_visible(sources,path,ids)
            probabilities, sizes, alive, origins = predict(core,model,held,family,replicate)
            truth = dict(fate=np.array([frozen.fate(e1[x["stem_id"]]) for x in held]),
                         dbh=np.array([empirical.number(e1[x["stem_id"]]["dbh"]) if frozen.fate(e1[x["stem_id"]])==1 else np.nan for x in held],float),
                         years=np.array([x["years"] for x in held]))
            cells = np.array([x["cell"] for x in held])
            labels = [("ALL",np.ones(len(held),bool))]
            labels += [("cell_"+str(c),cells==c) for c in np.unique(cells)]
            for label,take in labels:
                subtruth={k:v[take] for k,v in truth.items()}
                groups.append(metrics(family,replicate,label,subtruth,probabilities[:,take],sizes[:,take],alive[:,take]))
            p = probabilities.mean(axis=0)
            for low in np.arange(0,1,.1):
                take = (p>=low)&(p<low+.1+1e-12)
                if take.any():
                    curves.append(dict(mask="M"+str(family),replicate=replicate,bin_lower=float(low),
                                       n=int(take.sum()),predicted=float(p[take].mean()),observed=float(truth["fate"][take].mean())))
            truth_tax = np.array([truth["fate"][cells//16==t].sum() for t in range(4)],float)
            predicted_tax = np.array([alive[:,cells//16==t].sum(axis=1).mean() for t in range(4)])
            for tax in range(4):
                compositions.append(dict(mask="M"+str(family),replicate=replicate,taxon_group=tax,
                                          truth_share=float(truth_tax[tax]/max(truth_tax.sum(),1)),
                                          predicted_share=float(predicted_tax[tax]/max(predicted_tax.sum(),1))))
            for field in ("gx","gy","quadrat"):
                actual = [x for x in held if empirical.number(e1[x["stem_id"]][field]) is not None]
                errors = [float(e0[x["stem_id"]][field])-float(e1[x["stem_id"]][field]) for x in actual
                          if empirical.number(e0[x["stem_id"]][field]) is not None]
                coordinate_rows.append(dict(mask="M"+str(family),replicate=replicate,field=field,
                                            observed_n=len(actual),recovered_n=len(errors),
                                            mae=float(np.abs(errors).mean()) if errors else None,
                                            provenance="RECOVERED_FROM_PRIOR_OBSERVATION"))
            runs.append(dict(mask="M"+str(family),replicate=replicate,heldout=len(ids),
                             sanitized_visible_sha256=frozen.sha(path),
                             fit_survival_n=core["summary"]["survival_training_count"],
                             fit_growth_n=core["summary"]["growth_training_count"],
                             fitted_pooled_hazard=model["pooled_hazard"],
                             origin_count_median=float(np.median(np.array(origins)[:,0])),
                             origin_ba_median=float(np.median(np.array(origins)[:,1])),
                             fit_interface="Frozen read_core -> fit -> reconstruct; truth fields physically removed"))
            print("Calibration M%d/%d heldout %d" % (family,replicate,len(ids)),flush=True)
    write_csv("reconstruction_validation/mask_assignments.csv",assignments)
    write_csv("reconstruction_validation/calibration.csv",groups)
    write_csv("reconstruction_validation/fate_curve.csv",curves)
    write_csv("reconstruction_validation/coordinate_recovery.csv",coordinate_rows)
    write_csv("reconstruction_validation/taxon_composition.csv",compositions)
    write_json("reconstruction_validation/summary.json",dict(eligible=len(pool),runs=runs,
               alignment=cfg["target_alignment"],primary_truth="Secure stem-level fate, not whole-tree death",
               no_growth_outlier_truth_removal=True,masking_config_sha256=__import__("h03_common").PROTOCOL_HASH))
    return groups
