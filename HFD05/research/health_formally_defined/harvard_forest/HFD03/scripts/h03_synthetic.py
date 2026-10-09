"""Full predeclared truth grid; truth census is persisted BEFORE estimator work."""
import itertools
import math
import sys
import numpy as np
from h03_common import ROOT, RUN, PROTOCOL, demographic, rng, write_csv, write_json
from estimator import ObservationPacket

sys.path.insert(0,str(ROOT/"synthetic_validation"))
import importlib.util
spec=importlib.util.spec_from_file_location("h03_observation_estimator",ROOT/"synthetic_validation/estimator.py")
est=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=est
spec.loader.exec_module(est)


def make_state(core,model,s,r,world):
    cfg=PROTOCOL["synthetic"]
    generator=rng("world",world,0)
    n=core["stats"]["n0"].copy()
    d=model["diameter0"].copy()
    juvenile=(d<10)
    target_juveniles=cfg["regeneration"][r]*model["baseline_juveniles"]*generator.lognormal(-.03**2/2,.03)
    n[juvenile]=np.rint(n[juvenile]*target_juveniles/max(n[juvenile].sum(),1)).astype(int)
    target_ba=cfg["structure"][s]*model["baseline_ba"]*generator.lognormal(-.03**2/2,.03)
    juvenile_ba=(n[juvenile]*d[juvenile]**2).sum()*math.pi/40000
    adult_ba=(n[~juvenile]*d[~juvenile]**2).sum()*math.pi/40000
    d[~juvenile]*=math.sqrt(max(target_ba-juvenile_ba,1)/max(adult_ba,1))
    d=np.clip(d,1,300)
    return dict(n=n,d=d)


def truth_simulation(state,params,model,world,suite):
    cfg=PROTOCOL["synthetic"]
    scenario=dict(hazard_factor=1,growth_factor=1,recruitment_factor=1,
                  hemlock_extra_hazard=.12 if suite==1 else 0)
    K=cfg["truth_paths"]
    generator=rng("truth_future",world,suite)
    if suite!=2:
        return demographic.simulate(state,params,scenario,cfg["horizon"],K,generator)
    # Independent path-level annual lognormal entry intensity, absent from
    # estimator. No changing world selection; same means, broad dispersion.
    class EntryRNG:
        def __getattr__(self,name):
            if name!="poisson":
                return getattr(generator,name)
            def call(lam,size=None):
                multiplier=generator.lognormal(-1.2**2/2,1.2,size=(size[0],1))
                return generator.poisson(np.broadcast_to(lam,size)*multiplier)
            return call
    return demographic.simulate(state,params,scenario,cfg["horizon"],K,EntryRNG())


def scores(rows):
    y=np.array([r["truth_health"] for r in rows],bool)
    p=np.array([r["estimate_health_probability"] for r in rows])
    predicted=p>=.5
    tp=int((y&predicted).sum());tn=int((~y&~predicted).sum())
    fp=int((~y&predicted).sum());fn=int((y&~predicted).sum())
    sensitivity=tp/max(int(y.sum()),1)
    specificity=tn/max(int((~y).sum()),1)
    near=[r for r in rows if abs(r["truth_mass"]-PROTOCOL["synthetic"]["theta"])<=PROTOCOL["synthetic"]["near_boundary"]]
    truth_mass=np.array([r["truth_mass"] for r in rows])
    estimate=np.array([r["estimate_mass"] for r in rows])
    lower=np.array([r["estimate_lower"] for r in rows]);upper=np.array([r["estimate_upper"] for r in rows])
    bins=[]
    for lo in (0,.2,.4,.6,.8):
        take=(p>=lo)&(p<lo+.2 if lo<.8 else p<=1)
        if take.any():
            bins.append(dict(lower=lo,n=int(take.sum()),predicted=float(p[take].mean()),truth=float(y[take].mean())))
    return dict(worlds=len(rows),positives=int(y.sum()),negatives=int((~y).sum()),
                sensitivity=sensitivity,specificity=specificity,balanced_accuracy=(sensitivity+specificity)/2,
                brier=float(np.mean((p-y)**2)),viable_mass_mae=float(np.abs(estimate-truth_mass).mean()),
                viable_mass_bias=float((estimate-truth_mass).mean()),nominal90_mass_coverage=float(
                    ((truth_mass>=lower)&(truth_mass<=upper)).mean()),
                confusion=dict(TP=tp,TN=tn,FP=fp,FN=fn),calibration=bins,
                threshold_adjacent_n=len(near),threshold_adjacent_mae=float(np.mean([
                    abs(r["estimate_mass"]-r["truth_mass"]) for r in near])) if near else None,
                unresolved_truth_cells=sum(r["truth_mc_status"]=="MC_UNRESOLVED" for r in rows),
                clustering="324 independent world states; suites paired on same state, not 972 independent worlds")


def run(core,model):
    cfg=PROTOCOL["synthetic"]
    worlds=[]
    truth_rows=[]
    grid=list(itertools.product(range(3),range(3),range(3),range(3),range(cfg["replicates_per_cell"])))
    for world,(s,r,p,recovery,replicate) in enumerate(grid):
        state=make_state(core,model,s,r,world)
        context=est.PublicContext(p,recovery)
        params=est.public_parameters(model,context,rng("world",world,1))
        generator=rng("observation",world)
        counts=generator.binomial(state["n"],cfg["observation_thinning"])
        diameters=np.clip(state["d"]+generator.normal(0,cfg["dbh_noise_sd_cm"],64),1,300)
        packet=ObservationPacket(tuple(map(int,counts)),tuple(float(d) if n else None for n,d in zip(counts,diameters)),
                                 (cfg["observation_thinning"],)*64)
        worlds.append((packet,context))
        for suite,label in enumerate(cfg["misspecifications"]):
            bank,lawful=truth_simulation(state,params,model,world,suite)
            response=demographic.response(bank,lawful,model,cfg["horizon"],cfg["alpha"],cfg["beta"],
                                          cfg["tau"],cfg["reserve_ratio"])
            mass=float(response["viable"].mean())
            reserve=float(response["reserve_viable"].mean())
            status,lo,hi=__import__("h03_common").classify(response["present"],
                [int(response["viable"].sum()),int(response["reserve_viable"].sum())],len(bank),cfg["theta"])
            truth_rows.append(dict(world=world,suite=label,structure_level=s,regeneration_level=r,
                                   pressure_level=p,recovery_level=recovery,replicate=replicate,
                                   truth_mass=mass,truth_reserve=reserve,truth_present=response["present"],
                                   truth_health=bool(response["present"] and mass>=cfg["theta"] and reserve>=cfg["theta"]),
                                   truth_mc_status=status,truth_lower=float(min(lo)),truth_upper=float(min(hi))))
        if world%24==0:
            print("Synthetic truth census world",world,flush=True)
    write_csv("synthetic_validation/truth_census.csv",truth_rows)
    census={label:dict(positive=sum(r["truth_health"] for r in truth_rows if r["suite"]==label),
                       negative=sum(not r["truth_health"] for r in truth_rows if r["suite"]==label))
            for label in cfg["misspecifications"]}
    write_json("synthetic_validation/pre_estimator_census.json",dict(
        id=cfg["id"],worlds=len(grid),all_worlds_retained=True,truth_census=census,
        truth_census_sha256=__import__("h03_common").frozen.sha(ROOT/"synthetic_validation/truth_census.csv"),
        boundary="Truth Health is finite 4096-path reference, with its own MC status; not exact infinite-law oracle"))
    anchor=dict(n=core["stats"]["n0"].copy(),d=model["diameter0"].copy())
    rows=[]
    for world,(packet,context) in enumerate(worlds):
        result=est.evaluate(packet,context,anchor,model,world)
        for row in truth_rows[world*3:(world+1)*3]:
            rows.append(dict(row,estimate_mass=result["mass"],estimate_lower=result["mass_lower"],
                             estimate_upper=result["mass_upper"],estimate_reserve=result["reserve_mass"],
                             estimate_health_probability=result["health_probability"]))
        if world%24==0:
            print("Observation-only estimator world",world,flush=True)
    write_csv("synthetic_validation/performance.csv",rows)
    results={label:scores([r for r in rows if r["suite"]==label]) for label in cfg["misspecifications"]}
    correct=results["correct"]
    established=(correct["positives"]>=50 and correct["negatives"]>=50 and
                 correct["sensitivity"]>=.8 and correct["specificity"]>=.8)
    status="ESTABLISHED" if established else ("PARTIAL" if correct["positives"] and correct["negatives"] else "NOT_ESTABLISHED")
    write_json("synthetic_validation/summary.json",dict(benchmark=cfg["id"],results=results,
               discrimination=status,correct_and_misspecified_never_pooled=True,
               truth_params_in_estimator=False,selection="Full frozen grid, no pruning",
               prior_experiment="Pressure/recovery labels are observed imposed design contexts; actual rates and state hidden"))
    return results
