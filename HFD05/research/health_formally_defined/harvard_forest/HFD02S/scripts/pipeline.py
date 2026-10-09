#!/usr/bin/env python3
"""Replayable finite empirical-synthetic demonstration. Outputs stay isolated."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for directory in ["scripts","models","validation","ablation"]:
    sys.path.insert(0,str(ROOT/directory))
from contracts import load, sha, canonical, verify_frozen_inputs, Seeds, quantity
from empirical import read_core, fit
from demography import reconstruct, parameters, simulate, response, mass
from controls import hostile_tests, require_design, require_draws, require_reference
from worlds import validate_worlds
from representations import factor_checks, paired_scores
from external import external_check


def freeze_description():
    sources = {str(p.relative_to(ROOT)):sha(p) for p in sorted(ROOT.rglob("*.py"))
               if ".runtime" not in p.parts and ".runs" not in p.parts}
    configs = {name:sha(ROOT/name) for name in ["config/design.json","config/requirements.txt","config/future_protocol.json",
                                              "provenance/random_seed_registry.json",
                                              "provenance/provenance_schema.json",
                                              "reconstruction/observed_core_schema.json"]}
    return dict(schema="HFD02S-pre-outcome-freeze-v1",sources=sources,configs=configs,
                design_canonical_sha256=hashlib.sha256(canonical(load("config/design.json"))).hexdigest(),
                truth_isolation="validation/estimator.py does not import worlds; strict ObservationPacket",
                frozen_before_primary_results=True,no_E2_inputs=True)


def write_json(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n")


def write_csv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:
        raise ValueError("Empty required output: "+str(path))
    with path.open("w",newline="") as handle:
        writer = csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def progress(message):
    print(message,flush=True)


def summarize(bank,lawful,model,cfg):
    # Axis: state,param,model,scenario,path,annual time,metric.
    states,params,variants,scenarios,paths,times,_ = bank.shape
    outer = states*params*variants
    b = bank.reshape(outer,scenarios,paths,times,6)
    law = lawful.reshape(outer,scenarios,paths)
    structure = b[:,:,:,:,1]/model["baseline_ba"]
    juvenile = b[:,:,:,:,2]/model["baseline_juveniles"]
    reserve = np.divide(b[:,:,:,:,2],b[:,:,:,0:1,2],out=np.zeros_like(juvenile),
                        where=b[:,:,:,0:1,2]>0)
    rows,capacity,attrition,convergence = [],[],[],[]
    semantic = cfg["semantics"]
    reference_probability = None
    semantic_means = []
    for alpha in semantic["alpha_grid"]:
        for beta in semantic["beta_grid"]:
            realizes = (structure>=alpha)&(juvenile>=beta)
            longest = np.zeros((outer,scenarios,paths),int)
            current = longest.copy()
            present = realizes[:,:,0,0]
            for year in range(times):
                current = np.where(realizes[:,:,:,year],0,current+1)
                longest = np.maximum(longest,current)
                if year not in cfg["horizons_years"]:
                    continue
                for tau in semantic["recovery_years_grid"]:
                    viable = law&realizes[:,:,:,year]&(longest<=tau)
                    p = viable.mean(axis=2)
                    semantic_means.append(float(p.mean()))
                    if (alpha,beta,tau)==(.75,.4,2):
                        if year==10:
                            reference_probability = p.copy()
                        for s,scenario in enumerate(cfg["scenarios"]):
                            capacity.append(dict(scenario=scenario["id"],horizon=year,
                                                 viable_mass_mean=float(p[:,s].mean()),
                                                 viable_mass_low=float(np.quantile(p[:,s],.05)),
                                                 viable_mass_high=float(np.quantile(p[:,s],.95)),
                                                 failure_mass_mean=float(1-p[:,s].mean()),
                                                 present_fraction=float(present[:,s].mean())))
                            unlawful = float((~law[:,s]).mean())
                            attrition.append(dict(scenario=scenario["id"],horizon=year,candidate_mass=1.,
                                                  unlawful_mass=unlawful,
                                                  lawful_failure_mass=float((law[:,s]&~viable[:,s]).mean()),
                                                  viable_mass=float(p[:,s].mean())))
                    for theta in semantic["theta_grid"]:
                        for query,rho in [("Q1",None)]+[("Q2",x) for x in semantic["reserve_ratio_grid"]]:
                            rp = p if rho is None else (viable&(reserve[:,:,:,year]>=rho)).mean(axis=2)
                            adequate = (p>=theta)&(rp>=theta)
                            health = present&adequate
                            for s,scenario in enumerate(cfg["scenarios"]):
                                clustered = health[:,s].reshape(states,params,variants).mean(axis=(1,2))
                                rows.append(dict(scenario=scenario["id"],horizon=year,alpha=alpha,beta=beta,tau=tau,
                                                 theta=theta,query=query,reserve_ratio=rho,
                                                 present_fraction=float(present[:,s].mean()),
                                                 adequate_fraction=float(adequate[:,s].mean()),
                                                 health_ensemble_fraction=float(health[:,s].mean()),
                                                 state_cluster_mc_se=float(clustered.std(ddof=1)/np.sqrt(states)),
                                                 viable_mass_mean=float(p[:,s].mean()),
                                                 reserve_viable_mass_mean=float(rp[:,s].mean())))
                            if (alpha,beta,tau,theta,query,rho)==(.75,.4,2,.75,"Q2",1):
                                for prefix in cfg["ensemble"]["mc_prefixes"]:
                                    pp = viable[:,:,:prefix].mean(axis=2)
                                    rr = (viable[:,:,:prefix]&(reserve[:,:,:prefix,year]>=rho)).mean(axis=2)
                                    hp = present&(pp>=theta)&(rr>=theta)
                                    for s,scenario in enumerate(cfg["scenarios"]):
                                        convergence.append(dict(scenario=scenario["id"],horizon=year,prefix=prefix,
                                                                max_viable_mass_difference=float(np.max(abs(pp[:,s]-p[:,s]))),
                                                                mean_abs_viable_mass_difference=float(np.mean(abs(pp[:,s]-p[:,s]))),
                                                                health_disagreement_with_full=float(np.mean(hp[:,s]!=health[:,s]))))
    uncertainty = []
    q = reference_probability.reshape(states,params,variants,scenarios)
    # Descriptive conditional contrasts, not additive independent Sobol effects.
    effects = [
        ("imputation_and_measurement_joint",np.var(q.mean(axis=(1,2,3)),ddof=1)),
        ("parameter",np.mean(np.var(q.mean(axis=(2,3)),axis=1,ddof=1))),
        ("model",np.mean(np.var(q.mean(axis=3),axis=2,ddof=1))),
        ("scenario",np.mean(np.var(q,axis=3,ddof=1))),
        ("semantic",np.var(semantic_means,ddof=1))]
    for name,value in effects:
        uncertainty.append(dict(source=name,conditional_variance=float(value),
                                method="Nested conditional variance/finite-grid sensitivity; nonadditive, designs dependent",
                                metric="unconditional viable mass"))
    return rows,capacity,attrition,convergence,uncertainty


def execute():
    cfg = load("config/design.json")
    frozen = load("config/freeze.json")
    if freeze_description()!=frozen:
        raise ValueError("Code/config no longer equals pre-outcome freeze")
    require_design(cfg,frozen["design_canonical_sha256"])
    sources = verify_frozen_inputs()
    seeds = Seeds()
    controls = hostile_tests(seeds)
    factors = factor_checks()
    progress("observed core and E0→E1 fit (no E2)")
    core = read_core(sources); model = fit(core)
    e = cfg["ensemble"]
    run_id = hashlib.sha256(canonical(frozen)).hexdigest()[:16]
    out = ROOT/".runs"/run_id
    out.mkdir(parents=True,exist_ok=True)
    inputs = {k:{name:sources[k][name] for name in ["source_id","sha256","retrieval_url"]}
              for k in cfg["training_artifacts"]}
    write_json(out/"inputs.json",dict(input_manifest_sha256=cfg["hfd01_input_sha256"],artifacts=inputs,
                                    identity_ledger_sha256=core["identity_ledger_sha256"]))
    write_json(out/"observed.json",{k:core[k] for k in ["summary","raw_status","baseline_ba","baseline_juveniles","identity_ledger_sha256"]})
    write_csv(out/"missingness.csv",core["fields"])
    write_json(out/"fit.json",{k:x.tolist() if isinstance(x,np.ndarray) else x for k,x in model.items()})
    shape = (e["state_draws"],e["parameter_draws_per_state"],len(e["model_variants"]),
             len(cfg["scenarios"]),e["future_draws"],max(cfg["horizons_years"])+1,6)
    bank = np.empty(shape)
    lawful = np.empty(shape[:5],dtype=bool)
    states, state_rows = [],[]
    for i in range(e["state_draws"]):
        state = reconstruct(core,model,seeds.rng("reconstruct",i))
        states.append(state)
        c = state["components"]
        state_rows.append(dict(draw=i,living_stems=int(state["n"].sum()),
                               basal_area=float((state["n"]*state["d"]**2).sum()*np.pi/40000),
                               juvenile_support=int((state["n"]*(state["d"]<10)).sum()),
                               imputed_ba=c["imputed_ba"],measured_projected_ba=c["measured_projected_ba"],
                               wet_imputed_stems=c["wet_imputed_stems"],mnar_shift=c["mnar_shift"],
                               wet_shift=c["wet_shift"],association_e0_branches=c["association_branches"][0],
                               association_e1_branches=c["association_branches"][1]))
        for p in range(e["parameter_draws_per_state"]):
            param = parameters(model,seeds.rng("parameters",i,p))
            for m,variant in enumerate(e["model_variants"]):
                for s,scenario in enumerate(cfg["scenarios"]):
                    bank[i,p,m,s],lawful[i,p,m,s] = simulate(
                        state,param,scenario,max(cfg["horizons_years"]),e["future_draws"],
                        seeds.rng("future",i,p,m,s),variant)
        if (i+1)%16==0:
            progress("primary reconstruction/futures %d/%d"%(i+1,e["state_draws"]))
    require_draws(states)
    np.save(out/"primary_histories.npy",bank)
    np.save(out/"primary_lawful.npy",lawful)
    np.save(out/"latent_counts.npy",np.array([s["n"] for s in states]))
    np.save(out/"latent_diameters.npy",np.array([s["d"] for s in states]))
    write_csv(out/"reconstruction.csv",state_rows)
    progress("frozen R1/C1 and capacity/Health surfaces")
    phase,capacity,attrition,convergence,uncertainty = summarize(bank,lawful,model,cfg)
    write_csv(out/"health_phase.csv",phase)
    write_csv(out/"capacity.csv",capacity)
    write_csv(out/"attrition.csv",attrition)
    write_csv(out/"monte_carlo.csv",convergence)
    # Separating measurement from imputation remains a conditional diagnostic,
    # not an asserted independent variance decomposition.
    measured = np.array([s["components"]["measured_projected_ba"] for s in states])
    imputed = np.array([s["components"]["imputed_ba"] for s in states])
    for name,values in [("measurement_and_date_projection",measured),("imputation",imputed)]:
        uncertainty.append(dict(source=name,conditional_variance=float(np.var(values/model["baseline_ba"],ddof=1)),
                                method="Reconstruction component BA variance; includes coupled demographic projection; not separable causal attribution",
                                metric="initial structural ratio, NOT directly comparable to viable-mass variance"))
    write_csv(out/"uncertainty.csv",uncertainty)
    rich_rows = []
    for s,scenario in enumerate(cfg["scenarios"]):
        rr = response(bank[0,0,0,s],lawful[0,0,0,s],model,10,.75,.4,2,1)
        ledger = mass(rr["viable"],lawful[0,0,0,s],np.full(e["future_draws"],1/e["future_draws"]))
        if abs(sum(ledger[k] for k in ["unlawful_mass","lawful_failure_mass","viable_mass"])-1)>1e-12:
            raise ValueError("Unconditional mass conservation")
        for h in range(e["future_draws"]):
            rich_rows.append(dict(history_id="state000-param0-model0-"+scenario["id"]+"-path%03d"%h,
                                  unconditional_weight=1/e["future_draws"],outer_weight=1/(e["state_draws"]*2*2),
                                  lawful=bool(lawful[0,0,0,s,h]),viable=bool(rr["viable"][h]),
                                  final_structure=float(rr["structure"][h]),final_regeneration=float(rr["regeneration"][h]),
                                  regenerative_reserve=float(rr["reserve"][h]),
                                  hemlock_ba_fraction=float(bank[0,0,0,s,h,10,3]/max(bank[0,0,0,s,h,10,1],1e-9)),
                                  longest_excursion=int(rr["longest_excursion"][h]),
                                  provenance_class="SIMULATED",source_id=scenario["id"],method_id="cohort_demography_v1"))
    write_csv(out/"rich_capacity_fixture.csv",rich_rows)
    primary_digest = sha(out/"health_phase.csv")
    # No fitted input or above-primary computation reads external E2 aggregates.
    external = external_check(sources,out/"health_phase.csv")
    external["D0_five_year_primary_mean"] = {
        "cumulative_deaths":float(bank[:,:,:,0,:,5,4].mean()),
        "cumulative_threshold_entries":float(bank[:,:,:,0,:,5,5].mean())}
    write_json(out/"external.json",external)
    progress("truth-isolated synthetic reference and 32 worlds")
    scores, reference, anchor, packets = validate_worlds(core,model,seeds,progress)
    require_reference(0)
    for item in reference.pop("history_banks"):
        np.save(out/("reference_"+item["scenario"]+"_histories.npy"),item["bank"])
        np.save(out/("reference_"+item["scenario"]+"_lawful.npy"),item["lawful"])
    write_json(out/"reference_world.json",reference)
    write_csv(out/"world_scores.csv",scores)
    write_json(out/"observation_packets.json",[dict(counts=p.counts,diameters_cm=p.diameters_cm,
                                                     observation_probabilities=p.observation_probabilities) for p in packets])
    write_json(out/"public_anchor.json",{k:x.tolist() for k,x in anchor.items()})
    write_csv(out/"paired_ablation.csv",paired_scores(scores))
    write_json(out/"factor_surfaces.json",factors)
    write_json(out/"hostile_controls.json",controls)
    write_json(out/"seeds_used.json",seeds.used)
    summaries = dict(
        reconstructed_ba_90_interval=[float(x) for x in np.quantile([r["basal_area"] for r in state_rows],[.05,.5,.95])],
        reconstructed_stems_90_interval=[float(x) for x in np.quantile([r["living_stems"] for r in state_rows],[.05,.5,.95])],
        world_ba_coverage=float(np.mean([r["ba_covered"] for r in scores])),
        world_juvenile_coverage=float(np.mean([r["juvenile_covered"] for r in scores])),
        world_viable_mass_coverage=float(np.mean([r["viable_covered"] for r in scores])),
        mean_health_brier=float(np.mean([r["brier"] for r in scores])),
        true_health_positive_cases=sum(int(r["true_health"]) for r in scores),
        health_discrimination=("NOT_ESTABLISHED_DEGENERATE_TEST_SET"
                               if len(set(r["true_health"] for r in scores))<2 else "REPORT_SCORES_WITHOUT_ACCEPTANCE_GATE"),
        mean_absolute_viable_mass_error=float(np.mean([abs(r["estimated_viable_mass"]-r["true_viable_mass"]) for r in scores])),
        max_mc_health_disagreement=max(r["health_disagreement_with_full"] for r in convergence),
        all_candidate_lawful=bool(lawful.all()),phase_rows=len(phase),
        capacity_atom_census="All and only lawful C1 histories, unique outer/scenario/path identity; original weights retained",
        probability_wording=cfg["reporting"]["health_probability_wording"])
    write_json(out/"summary.json",summaries)
    # Every quantity family has an explicit class/source/method contract.
    families = dict(
        recorded_census=quantity("immutable tagged CSV fields","OBSERVED","HF253-E0/E1","raw_csv"),
        recovered_locations=quantity("stable recorded root identities","RECOVERED_FROM_PRIOR_OBSERVATION","HF253-E0","stationary_root"),
        observed_aggregates=quantity("BA/juvenile strata","DERIVED","HF253-E0/E1","cohort_aggregation"),
        latent_origin=quantity("256 count/diameter states","IMPUTED","HF253-E0/E1","MNAR_cohort_reconstruction"),
        model_parameters=quantity("interval transition likelihood","DERIVED","HF253-E0/E1","tempered_demographic_fit"),
        scenarios_and_thresholds=quantity("frozen grids","EXTERNALLY_CONSTRAINED","HFD02S-frozen-assumptions","declared_not_biologically_established"),
        response_histories=quantity("unconditional cohort futures","SIMULATED","D0-D3","demographic_kernel"),
        biological_mechanism_rationale=quantity("qualitative drought/hemlock pathways, NOT magnitudes","LITERATURE_LICENSED","HFD01-L04/L06/L13","primary_source_review"),
        synthetic_truth=quantity("fixed finite reference laws and latent worlds","SIMULATED","HF-SYNTH-REF-001","truth_harness"),
        health_surfaces=quantity("R1 AND Q1/Q2 over declared ensemble","DERIVED","HFD02S-primary-histories","frozen_semantic_evaluation"),
        external_e2=quantity("already exposed published aggregates","OBSERVED","HF-third-census-news","external_check_only"))
    write_json(out/"field_provenance.json",families)
    column_contract = {}
    for path in sorted(out.glob("*.csv")):
        with path.open(newline="") as handle:
            columns = next(csv.reader(handle))
        pc = "IMPUTED" if path.name=="reconstruction.csv" else "DERIVED"
        if path.name=="rich_capacity_fixture.csv":
            pc = "SIMULATED"
        column_contract[path.name] = {
            col:dict(provenance_class=("SIMULATED" if col.startswith("true_") else pc),
                     source_id="HFD02S-frozen-run-"+run_id,
                     method_id=path.stem+"_field_derivation",
                     value_location=path.name+"::"+col,
                     parents=["HF253-E0","HF253-E1","HFD02S-frozen-assumptions"])
            for col in columns}
    write_json(out/"column_provenance.json",column_contract)
    outputs = {p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file() and p.name!="run_manifest.json"}
    manifest = dict(schema="HFD02S-run-v1",run_id=run_id,freeze_sha256=sha(ROOT/"config/freeze.json"),
                    input_manifest_sha256=cfg["hfd01_input_sha256"],numpy_version=np.__version__,
                    seeds_used=len(seeds.used),outputs=outputs,summary=summaries,
                    terminal_statuses=dict(OPERATIONAL_HARVARD_FOREST_ASSESSMENT="NOT_AUTHORIZED",
                                           EMPIRICAL_SYNTHETIC_DEMONSTRATION="COMPLETE",
                                           RETROSPECTIVE_E2_VALIDATION="NOT_CLAIMED",
                                           E2_EXTERNAL_CONSISTENCY_CHECK="PERFORMED",
                                           FUTURE_PROSPECTIVE_PROTOCOL="FROZEN"))
    write_json(out/"run_manifest.json",manifest)
    progress("run complete: "+str(out))
    return out


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--describe-freeze",action="store_true")
    args = parser.parse_args()
    if args.describe_freeze:
        print(json.dumps(freeze_description(),indent=2,sort_keys=True))
    else:
        execute()
