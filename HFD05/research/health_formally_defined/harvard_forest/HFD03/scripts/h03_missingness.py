"""Five frozen reconstruction families propagated through unconditioned futures."""
import itertools
import math
import time
import numpy as np
from h03_common import (RUN, PROTOCOL, frozen, demographic, rng, reconstruct,
                        wilson, write_csv, write_json)


def query_grid(old_config):
    sem = old_config["semantics"]
    cells=[]
    for h,a,b,t in itertools.product(PROTOCOL["missingness"]["horizons"],
                                     sem["alpha_grid"],sem["beta_grid"],sem["recovery_years_grid"]):
        for q,r in [("Q1",0.5),("Q2",0.5),("Q2",1)]:
            for theta in sem["theta_grid"]:
                cells.append((h,a,b,t,q,r,theta))
    return cells


def summarize_bank(bank,lawful,model,cells):
    cached={}
    mass, reserve, present = [],[],[]
    for h,a,b,t,q,r,theta in cells:
        key=(h,a,b,t)
        if key not in cached:
            cached[key] = demographic.response(bank,lawful,model,h,a,b,t,0.5)
        response=cached[key]
        v=response["viable"]
        rv=v&(response["reserve"]>=r)
        mass.append(v.sum())
        reserve.append(rv.sum())
        present.append(response["present"])
    v=np.array(mass)
    r=np.array(reserve)
    p=np.array(present)
    theta=np.array([c[-1] for c in cells])
    q2=np.array([c[4]=="Q2" for c in cells])
    vl,vu=wilson(v,len(bank),0.95)
    ql,qu=wilson(v,len(bank),0.975)
    rl,ru=wilson(r,len(bank),0.975)
    lower=np.where(q2,np.minimum(ql,rl),vl)
    upper=np.where(q2,np.minimum(qu,ru),vu)
    status=np.where(~p | (upper<theta),0,np.where(lower>theta,1,2))
    point=p&(v/len(bank)>=theta)&(~q2|(r/len(bank)>=theta))
    return v/len(bank),r/len(bank),point,status


def run(core,model):
    cfg=PROTOCOL["missingness"]
    old=frozen.load("config/design.json")
    cells=query_grid(old)
    N=cfg["state_draws"]
    K=cfg["future_draws"]
    models=cfg["model_variants"]
    scenarios=old["scenarios"]
    unit_count=N*len(models)
    phase_rows,state_rows,provenance_rows=[],[],[]
    reference=None
    results={}
    for family in ("B0","B1","B2","B3","B4"):
        start=time.monotonic()
        masses=np.empty((len(scenarios),unit_count,len(cells)),np.float32)
        reserves=np.empty_like(masses)
        points=np.empty(masses.shape,bool)
        statuses=np.empty(masses.shape,np.uint8)
        states=[]
        for i in range(N):
            state=reconstruct(core,model,i,family)
            states.append(state)
            n=state["n"]
            ba=float((n*state["d"]**2).sum()*math.pi/40000)
            juv=int((n*(state["d"]<10)).sum())
            state_rows.append(dict(family=family,state=i,count=int(n.sum()),basal_area=ba,juveniles=juv,
                                   wet_count=int(n[(np.arange(64)//4)%4==3].sum())))
            provenance_rows.append(dict(family=family,state=i,
                observed_anchored_projected_ba=state["components"]["measured_projected_ba"],
                imputed_ba=state["components"]["imputed_ba"],
                observed_anchor_fraction=state["components"]["measured_projected_ba"]/max(ba,1e-9),
                directly_observed_origin_ba=0,origin_provenance="IMPUTED",
                location_evidence="OBSERVED or RECOVERED_FROM_PRIOR_OBSERVATION; no invented positions"))
            params=demographic.parameters(model,rng("future",i,0))
            for m,variant in enumerate(models):
                for s,scenario in enumerate(scenarios):
                    bank,lawful=demographic.simulate(state,params,scenario,20,K,rng("future",i,m,s+1),variant)
                    v,r,p,z=summarize_bank(bank,lawful,model,cells)
                    unit=i*len(models)+m
                    masses[s,unit]=v
                    reserves[s,unit]=r
                    points[s,unit]=p
                    statuses[s,unit]=z
            if i%16==0:
                print("Missingness",family,"state",i,flush=True)
        if reference is None:
            reference=(masses.copy(),points.copy())
        for s,scenario in enumerate(scenarios):
            for j,(h,a,b,t,q,r,theta) in enumerate(cells):
                v=masses[s,:,j]
                z=statuses[s,:,j]
                phase_rows.append(dict(family=family,scenario=scenario["id"],horizon=h,alpha=a,beta=b,tau=t,
                    query=q,rho=r,theta=theta,viable_mass=float(v.mean()),
                    viable_p05=float(np.quantile(v,.05)),viable_p95=float(np.quantile(v,.95)),
                    reserve_mass=float(reserves[s,:,j].mean()),point_health_fraction=float(points[s,:,j].mean()),
                    resolved_health_lower=float((z==1).mean()),resolved_health_upper=float((z!=0).mean()),
                    true_units=int((z==1).sum()),false_units=int((z==0).sum()),unresolved_units=int((z==2).sum()),
                    paired_health_disagreement=float((points[s,:,j]!=reference[1][s,:,j]).mean()),
                    paired_mass_change=float((v-reference[0][s,:,j]).mean()),
                    state_provenance_link="provenance/state_contributions.csv",
                    weight_interpretation="Equal finite state/model design weights, not ecological Health probability"))
        np.savez_compressed(RUN/("missingness_"+family+".npz"),mass=masses,reserve=reserves,
                            point=points,status=statuses)
        results[family]=dict(phase_cells=len(cells)*len(scenarios),conditional_units=unit_count,
                            point_disagreement=float(np.mean(points!=reference[1])),
                            mean_mass_change=float(np.mean(masses-reference[0])),
                            elapsed_seconds=time.monotonic()-start)
        print("Completed",family,results[family]["point_disagreement"],flush=True)
    # Do not commit clock-dependent execution durations into scientific results.
    for value in results.values():
        value.pop("elapsed_seconds")
    write_csv("missingness/phase_surface.csv",phase_rows)
    write_csv("missingness/state_distribution.csv",state_rows)
    write_csv("provenance/state_contributions.csv",provenance_rows)
    write_json("missingness/summary.json",dict(results=results,config=cfg,
        phase_metric="Disagreement across paired conditional Health cells; also paired continuous mass changes",
        unresolved_policy="Bracket from conditional TRUE/FALSE/MC_UNRESOLVED; point labels descriptive only"))
    return results
