"""Read-only frozen tagged records; selective demonstrative transition fit.

No fate from absence, no overwrite, no synonym/plant crosswalk invented.
The returned X is explicitly cohort-resolution, not fabricated wet-area trees.
"""
import csv
import datetime as dt
import math
import hashlib
from collections import Counter
import numpy as np
from contracts import HFD01, load, fate, recover_coordinate


def number(value):
    return None if value in {None, "", "NA"} else float(value)


def taxon(code):
    if code == "tsugca":
        return 0
    if code in {"acerru", "querru", "pinust", "fagugr", "betual", "betule", "betupa", "acersa"}:
        return 1
    if code in {"ilexve", "vaccco", "kalmla", "hamavi", "ilexmu", "nemomu", "lyonli", "alnuin", "vibuca", "vibunu", "vibual", "vibula"}:
        return 2
    return 3


def sector(row, old=None):
    if old is not None and sector(old) == 3:
        return 3  # Fixed E0 proxy membership, not a newly certified wet mask.
    source = row if number(row.get("gx")) is not None else old
    if source is None:
        return 1
    x = float(source["gx"])
    s = min(2, max(0, int(x / (700 / 3))))
    date = source.get("exact.date", "")
    if s == 1 and date[:4] in {"2012", "2013", "2014"} and date[5:7] in {"11", "12", "01", "02", "03"}:
        return 3
    return s


def cell(row, old=None):
    diameter = number(row.get("dbh"))
    if diameter is None and old:
        diameter = number(old.get("dbh"))
    diameter = diameter if diameter is not None else 1.5
    size = int(np.searchsorted([3, 10, 30], diameter, side="right"))
    return (taxon(row["sp"]) * 4 + sector(row, old)) * 4 + size


def read_core(sources):
    tables = []
    for key in ["HF253-E0", "HF253-E1"]:
        with (HFD01 / sources[key]["local_cache_path"]).open(newline="") as handle:
            tables.append({r["stem.id"]: r for r in csv.DictReader(handle)})
    e0, e1 = tables
    config = load("config/design.json")
    origin = dt.date.fromisoformat(config["origin"])
    summary = Counter()
    field_summary = []
    for key, table in zip(["E0", "E1"], tables):
        for field in ["tree.id", "sp", "dbh", "gx", "gy", "quadrat", "exact.date", "df.status"]:
            missing = sum(r[field] in {"", "NA"} for r in table.values())
            recovered = 0
            if key == "E1" and field in {"gx", "gy", "quadrat"}:
                recovered = sum(r[field] in {"", "NA"} and recover_coordinate(e0.get(k), r, field) is not None
                                for k, r in table.items())
            field_summary.append(dict(epoch=key, field=field, rows=len(table), observed=len(table)-missing,
                                      missing=missing, recovered=recovered, unresolved=missing-recovered))
    stats = {"n0": np.zeros(64, int), "ba0": np.zeros(64), "n1": np.zeros(64, int), "ba1": np.zeros(64),
             "unknown": np.zeros(64, int), "unknown_ba0": np.zeros(64), "unknown_dt": np.zeros(64),
             "tail_years": np.zeros(64), "alive": np.zeros(64, int), "dead": np.zeros(64, int),
             "entrants": np.zeros(64, int), "alive_missing_size": np.zeros(64, int),
             "alive_missing_size_tail": np.zeros(64)}
    durations = [[] for _ in range(64)]
    growth = [[] for _ in range(64)]
    raw_status = Counter(r["df.status"] for r in e1.values())
    conflicts = []
    for key, old in e0.items():
        if old["df.status"] != "alive":
            continue
        c = cell(old)
        d0 = float(old["dbh"])
        stats["n0"][c] += 1
        stats["ba0"][c] += math.pi * d0*d0 / 40000
        new = e1.get(key)
        identity = new is not None and new["tree.id"] == old["tree.id"]
        same_taxon = new is not None and new["sp"] == old["sp"]
        f = fate(new) if identity else None
        if new is None:
            summary["absent_e1_not_death"] += 1
        elif not identity:
            summary["association_conflict"] += 1
        elif not same_taxon:
            summary["taxon_code_change_training_excluded"] += 1
        if new is not None and not identity:
            conflicts.append(dict(old_cell=c, new_cell=cell(new, old),old_diameter=d0,
                                  new_diameter=number(new["dbh"]),new_fate=fate(new),
                                  tail=max(0,(origin-dt.date.fromisoformat(new["exact.date"])).days/365.25)
                                  if fate(new) is not None else 0,
                                  years=(origin-dt.date.fromisoformat(old["exact.date"])).days/365.25))
        if f is None and (new is None or identity):
            stats["unknown"][c] += 1
            stats["unknown_ba0"][c] += math.pi*d0*d0/40000
            stats["unknown_dt"][c] += (origin-dt.date.fromisoformat(old["exact.date"])).days / 365.25
            summary["unresolved_fate_not_assumed_dead"] += 1
        if f is None or not same_taxon:
            continue
        interval = (dt.date.fromisoformat(new["exact.date"])-dt.date.fromisoformat(old["exact.date"])).days/365.25
        if interval <= 0:
            summary["noncausal_measurement_excluded"] += 1
            continue
        durations[c].append((interval, f))
        stats["alive" if f else "dead"][c] += 1
        if f:
            d1 = number(new["dbh"])
            h0, h1 = number(old["hom"]), number(new["hom"])
            if d1 is not None and h0 is not None and h1 is not None and abs((h0 or 1.3)-(h1 or 1.3)) <= config["fit"]["same_pom_tolerance_m"]:
                g = (d1-d0)/interval
                lo, hi = config["fit"]["growth_cm_per_year_inclusion"]
                if lo <= g <= hi:
                    growth[c].append(g)
                else:
                    summary["growth_outlier_quarantined_not_overwritten"] += 1
    # Every dated living E1 measurement contributes observed structure, even
    # when its taxon recoding excludes it from transition fitting.
    for key, row in e1.items():
        old = e0.get(key)
        if old is not None and row["tree.id"] != old["tree.id"]:
            continue  # One unresolved stem, never two latent copies.
        if fate(row) == 1 and number(row["dbh"]) is not None:
            c = cell(row, old)
            diameter = float(row["dbh"])
            stats["n1"][c] += 1
            stats["ba1"][c] += math.pi*diameter*diameter/40000
            stats["tail_years"][c] += max(0, (origin-dt.date.fromisoformat(row["exact.date"])).days/365.25)
            if old is None:
                stats["entrants"][c] += 1
        elif fate(row) == 1:
            # Alive is recorded but DBH missing; latent size is reconstructed,
            # not the observation overwritten.
            summary["dated_alive_missing_dbh"] += 1
            stats["alive_missing_size"][cell(row, old)] += 1
            stats["alive_missing_size_tail"][cell(row, old)] += max(
                0,(origin-dt.date.fromisoformat(row["exact.date"])).days/365.25)
    summary.update(e0_rows=len(e0), e1_rows=len(e1), union_stem_ids=len(e0.keys()|e1.keys()),
                   growth_training_count=sum(map(len,growth)), survival_training_count=sum(map(len,durations)),
                   observed_new_id_alive_count=int(stats["entrants"].sum()))
    return dict(stats=stats, durations=durations, growth=growth, summary=dict(summary),
                fields=field_summary, raw_status=dict(raw_status),
                conflicts=conflicts,
                identity_ledger_sha256=hashlib.sha256("".join(
                    repr((key,e0.get(key,{}).get("tree.id"),e1.get(key,{}).get("tree.id"),
                          e0.get(key,{}).get("sp"),e1.get(key,{}).get("sp")))+"\n"
                    for key in sorted(e0.keys()|e1.keys())).encode()).hexdigest(),
                baseline_ba=float(stats["ba0"].sum()), baseline_juveniles=int(stats["n0"].reshape(16,4)[:,:2].sum()))


def fit(core):
    cfg = load("config/design.json")["fit"]
    stats = core["stats"]
    known = stats["alive"]+stats["dead"]
    dt_all = [x for cell in core["durations"] for x in cell]
    global_survival = stats["alive"].sum()/known.sum()
    pooled_h = -math.log(global_survival) / np.mean([x[0] for x in dt_all])
    grid = np.linspace(-math.log(cfg["annual_survival_bounds"][1]),
                       -math.log(cfg["annual_survival_bounds"][0]), 256)
    all_growth = np.array([x for cell in core["growth"] for x in cell])
    mean = np.zeros(64)
    sd = np.zeros(64)
    se = np.zeros(64)
    posterior = np.zeros((64,len(grid)))
    for c in range(64):
        strength = cfg["survival_beta_prior_strength"]
        a = math.exp(-pooled_h)*strength
        b = strength-a
        prior = -a*grid + (b-1)*np.log(-np.expm1(-grid))  # Beta on annual survival, Jacobian included.
        data = core["durations"][c]
        if data:
            exposure_live = sum(t for t,z in data if z)
            td = np.array([t for t,z in data if not z])
            ll = -grid*exposure_live
            if len(td):
                ll += np.log(-np.expm1(-grid[:,None]*td[None,:])).sum(axis=1)
            prior += ll * min(1, cfg["parameter_effective_sample_cap"]/len(data))
        p = np.exp(prior-prior.max())
        posterior[c] = p/p.sum()
        g = np.array(core["growth"][c])
        # Shrink sparse cells to pooled demonstrative growth; retain variation.
        n = len(g)
        mean[c] = np.clip(((g.sum() if n else 0)+32*all_growth.mean())/(n+32), *cfg["growth_mean_bounds"])
        sd[c] = max(0.02, float(g.std()) if n>1 else float(all_growth.std()))
        se[c] = sd[c]/math.sqrt(min(max(n,1),cfg["parameter_effective_sample_cap"]))
    recruits = stats["entrants"].astype(float)/cfg["recruitment_exposure_years"]
    # Unmeasured wet proxy has no fitted observed recruitment exposure.
    # Borrow species/size share under an explicitly wide latent component.
    for tax in range(4):
        cells = np.arange(tax*16,(tax+1)*16)
        wet = cells[(cells//4)%4==3]
        dry = cells[(cells//4)%4!=3]
        share = stats["n0"][wet].sum()/max(stats["n0"][dry].sum(),1)
        recruits[wet[0]] = recruits[dry].sum()*share
    diameter0 = np.sqrt(np.divide(stats["ba0"]*40000/math.pi,stats["n0"],out=np.zeros(64),where=stats["n0"]>0))
    diameter1 = np.sqrt(np.divide(stats["ba1"]*40000/math.pi,stats["n1"],out=diameter0.copy(),where=stats["n1"]>0))
    diameter0[diameter0==0] = np.tile([1.7,6,18,40],16)[diameter0==0]
    diameter1[diameter1==0] = diameter0[diameter1==0]
    model = dict(hazard_grid=grid,hazard_posterior=posterior,growth_mean=mean,growth_sd=sd,growth_se=se,
                 recruitment_mean=recruits,diameter0=diameter0,diameter1=diameter1,
                 pooled_hazard=pooled_h,baseline_ba=core["baseline_ba"],baseline_juveniles=core["baseline_juveniles"])
    return model
