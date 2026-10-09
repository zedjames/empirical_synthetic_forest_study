"""Maximum-first paired-prefix MC; conditional Health intervals and state blocks."""
import itertools
import numpy as np
from h03_common import ROOT, RUN, PROTOCOL, frozen, demographic, rng, reconstruct, wilson, classify, write_csv, write_json


def count_queries(bank,lawful,model):
    cfg=PROTOCOL["monte_carlo"]
    counts,present=[],[]
    for horizon in cfg["horizons"]:
        response=demographic.response(bank,lawful,model,horizon,cfg["alpha"],cfg["beta"],cfg["tau"],.5)
        counts.append([int(response["viable"].sum()),int(response["reserve_viable"].sum()),
                       int((response["viable"]&(response["reserve"]>=1)).sum())])
        present.append(response["present"])
    return np.array(counts,np.uint16),np.array(present,bool)


def query_counts(counts,query):
    return counts[:1] if query==0 else counts[[0,query]]


def run(core,model):
    cfg=PROTOCOL["monte_carlo"]
    scenarios=[s for s in frozen.load("config/design.json")["scenarios"] if s["id"] in cfg["scenarios"]]
    rows=[]
    queries=[("Q1",.5),("Q2",.5),("Q2",1)]
    maxK=max(cfg["K"])
    for i in range(cfg["inner_states"]):
        state=reconstruct(core,model,i)
        for p in range(cfg["outer_parameter_draws"]):
            params=demographic.parameters(model,rng("inner_future",i,p,0))
            for m,variant in enumerate(cfg["outer_model_variants"]):
                for s,scenario in enumerate(scenarios):
                    bank,lawful=demographic.simulate(state,params,scenario,20,maxK,rng("inner_future",i,p,m,s+1),variant)
                    for K in cfg["K"]:
                        counts,present=count_queries(bank[:K],lawful[:K],model)
                        for h,horizon in enumerate(cfg["horizons"]):
                            for q,(query,rho) in enumerate(queries):
                                selected=query_counts(counts[h],q)
                                for theta in cfg["theta"]:
                                    status,lower,upper=classify(present[h],selected,K,theta)
                                    rows.append(dict(state=i,parameter=p,model=variant,scenario=scenario["id"],
                                                     K=K,horizon=horizon,query=query,rho=rho,theta=theta,
                                                     present=bool(present[h]),viable_count=int(counts[h,0]),
                                                     reserve_count=int(counts[h,q]) if q else int(counts[h,1]),
                                                     viable_mass=float(counts[h,0]/K),reserve_mass=float(counts[h,q]/K) if q else float(counts[h,1]/K),
                                                     lower=float(min(lower)),upper=float(min(upper)),status=status,
                                                     max_half_width=float(np.max((upper-lower)/2)),
                                                     provenance_link="provenance/state_contributions.csv"))
        print("MC paired inner state",i,flush=True)
    write_csv("monte_carlo/inner_resolution.csv",rows)
    N=max(cfg["outer_N"])
    units=cfg["outer_parameter_draws"]*len(cfg["outer_model_variants"])
    counts=np.empty((N,units,len(scenarios),len(cfg["horizons"]),3),np.uint16)
    present=np.empty(counts.shape[:-1],bool)
    provenance=[]
    for i in range(N):
        state=reconstruct(core,model,i)
        total_ba=float((state["n"]*state["d"]**2).sum()*np.pi/40000)
        provenance.append(dict(family="MC_B0",state=i,count=int(state["n"].sum()),basal_area=total_ba,
            observed_anchored_projected_ba=state["components"]["measured_projected_ba"],
            imputed_ba=state["components"]["imputed_ba"],
            observed_anchor_fraction=state["components"]["measured_projected_ba"]/max(total_ba,1e-9),
            directly_observed_origin_ba=0,origin_provenance="IMPUTED"))
        for p in range(cfg["outer_parameter_draws"]):
            params=demographic.parameters(model,rng("outer_future",i,p,0))
            for m,variant in enumerate(cfg["outer_model_variants"]):
                for s,scenario in enumerate(scenarios):
                    bank,lawful=demographic.simulate(state,params,scenario,20,cfg["outer_K"],
                                                    rng("outer_future",i,p,m,s+1),variant)
                    counts[i,p*2+m,s],present[i,p*2+m,s]=count_queries(bank,lawful,model)
        if i%64==0:
            print("MC outer state",i,flush=True)
    # Compact sufficient census, not an enormous trajectory bank. This permits
    # independent recomputation of every reported prefix/status/state-block SE.
    target=ROOT/"monte_carlo/outer_sufficient_statistics.npz"
    np.savez_compressed(target,counts=counts,present=present)
    write_csv("provenance/mc_state_contributions.csv",provenance)
    outer=[]
    for n in cfg["outer_N"]:
        for s,scenario in enumerate(scenarios):
            for h,horizon in enumerate(cfg["horizons"]):
                for q,(query,rho) in enumerate(queries):
                    v=counts[:n,:,s,h,0]/cfg["outer_K"]
                    r=counts[:n,:,s,h,q]/cfg["outer_K"] if q else v
                    for theta in cfg["theta"]:
                        point=present[:n,:,s,h]&(v>=theta)&(r>=theta)
                        states=point.mean(axis=1)
                        se=float(states.std(ddof=1)/np.sqrt(n))
                        mass_se=float(v.mean(axis=1).std(ddof=1)/np.sqrt(n))
                        statuses=[]
                        for i,u in itertools.product(range(n),range(units)):
                            statuses.append(classify(present[i,u,s,h],query_counts(counts[i,u,s,h],q),
                                                     cfg["outer_K"],theta)[0])
                        z=np.array(statuses)
                        outer.append(dict(N_X=n,outer_units=n*units,K=cfg["outer_K"],scenario=scenario["id"],
                                          horizon=horizon,query=query,rho=rho,theta=theta,
                                          viable_mass=float(v.mean()),state_block_mass_se=mass_se,
                                          point_health_fraction=float(point.mean()),state_block_health_se=se,
                                          health_fraction_lower=float((z=="TRUE").mean()),
                                          health_fraction_upper=float((z!="FALSE").mean()),
                                          unresolved_fraction=float((z=="MC_UNRESOLVED").mean()),
                                          weighting="Equal finite state/parameter/model design; block SE conditional on assumed state sampler"))
    write_csv("monte_carlo/outer_resolution.csv",outer)
    maximum=[r for r in rows if r["K"]==maxK]
    resolved=sum(r["status"]!="MC_UNRESOLVED" for r in maximum)/len(maximum)
    halfwidth=max(r["max_half_width"] for r in maximum)
    outer_se=max(r["state_block_health_se"] for r in outer if r["N_X"]==N)
    precision="ADEQUATE" if resolved>=.9 and halfwidth<=.04 and outer_se<=.03 else (
        "MIXED" if resolved>0 else "INADEQUATE")
    write_json("monte_carlo/summary.json",dict(config=cfg,maximum_K=maxK,maximum_N=N,
                inner_resolved_fraction=resolved,inner_max_half_width=halfwidth,
                outer_max_state_block_health_se=outer_se,precision=precision,
                inner_conditional_cells=len(rows),outer_surface_rows=len(outer),
                optional_4096_extension="Not executed; no adaptive stopping",
                outer_count_sha256=frozen.sha(target),
                distinction="Inner K study on fixed 16-state panel; outer N study uses K64. No claim of N1024 with K1024 everywhere.",
                systematic_limits="MC precision is not ecological validity or missingness/model/semantic uncertainty"))
    return rows,outer
