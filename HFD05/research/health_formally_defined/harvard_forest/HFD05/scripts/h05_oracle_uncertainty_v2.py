"""Explicit second uncertainty implementation: retain undefined empty strata."""
import itertools,json,math
from collections import defaultdict,Counter
import numpy as np
import h05_common as h
import h05_oracle_engine as e
from h05_inference import wilson

def boolean(v):return v is True or v=="True"
def audit():
    h.guard("oracle");h.guard("oracle_uncertainty")
    raw=h.read_csv(h.ROOT/"oracle_attribution/draw_census.csv");old={(r["world"],r["draw"]):r for r in h.read_csv(h.FOUR/"context_ablation/estimator_draws.csv") if r["cohort"]=="boundary" and r["policy"]=="B0_full"};keys=set()
    for r in raw:
        nv,nr,K=map(int,[r["nv"],r["nr"],r["K"]]);assert 0<=nr<=nv<=K
        assert boolean(r["finite_health"])==(boolean(r["present"]) and (nv if r["query"]=="Q1" else nr)*4>=K*3)
        roles={"E00":("EstimatedState","EstimatedParameters"),"E10":("OracleState","EstimatedParameters"),"E01":("EstimatedState","OracleParameters"),"E11":("OracleState","OracleParameters")}
        assert (r["role_state"],r["role_parameters"])==roles[r["condition"]]
        if r["condition"]=="E00" and r["replicate"]=="0" and r["K"]=="256":
            o=old[(r["world"],r["draw"])];assert r["nv"]==o["nv"] and r["nr"]==o["nr"] and r["present"]==o["present"];keys.add((r["world"],r["draw"]))
    assert len(keys)==81*32
    dummy=dict(n=np.ones(64,int),d=np.ones(64));pp={k:np.ones(64)*.01 for k in ["hazard","growth","sd","recruit"]};rejections=0
    for sr,pr in [(e.OracleState,e.EstimatedParameters),(e.EstimatedState,e.OracleParameters),(e.OracleState,e.OracleParameters)]:
        try:e.predict("E00",e.pack_state(dummy,sr),e.pack_parameters(pp,pr),1,h.rng("fixture",13))
        except ValueError:rejections+=1
        else:raise AssertionError("Oracle role accepted by E00")
    contrasts=h.read_csv(h.ROOT/"oracle_attribution/interactions.csv")
    losses=defaultdict(dict)
    rows=h.read_csv(h.ROOT/"oracle_attribution/paired_world_losses.csv")
    for r in rows:losses[tuple(r[k] for k in ["world","panel","replicate","K","suite","query"])][r["condition"]]=float(r["mass_absolute_error"])
    for r in contrasts:
        L=losses[tuple(r[k] for k in ["world","panel","replicate","K","suite","query"])];assert set(L)=={"E00","E10","E01","E11"}
        assert math.isclose(float(r["interaction"]),L["E11"]-L["E10"]-L["E01"]+L["E00"],abs_tol=1e-15)
    h.write_json("oracle_attribution/independent_audit.json",dict(status="PASS",original_E00_draws_exact=2592,strict_role_hostile_rejections=rejections,
        oracle_draw_query_prefix_rows=len(raw),complete_factorial_cells=len(losses),all_worlds=81,source_sha256=h.sha(__file__),forced_additivity=False))
    print("PASS independent original oracle audit",len(raw),flush=True)

def uncertainty():
    h.guard("oracle_uncertainty_v2");from h04_analysis import score
    rows=h.read_csv(h.ROOT/"oracle_attribution/paired_world_losses.csv");parsed=[]
    for r in rows:
        r=dict(r)
        for key in ["health_fraction","estimated_viable","estimated_reserve","reference_viable","reference_reserve","viable_lower","viable_upper","reserve_lower","reserve_upper","margin","mass_absolute_error"]:r[key]=float(r[key])
        parsed.append(r)
    intervals=[];near=[];paired=[]
    for qi,query in enumerate(["Q1","Q2"]):
        for si,suite in enumerate(h.original.config()["boundary"]["misspecifications"]):
            for ki,(panel,K) in enumerate([("BASELINE32",256),("PREFIX8",256),("PREFIX8",1024),("PREFIX8",4096)]):
                base=[r for r in parsed if r["query"]==query and r["suite"]==suite and r["panel"]==panel and r["K"]==str(K) and r["replicate"]=="0"]
                for subi,subset in enumerate(["ALL_ROSTER","PRESENT_NEAR_0.1"]):
                    selected=base if subset=="ALL_ROSTER" else [r for r in base if boolean(r["reference_present"]) and abs(r["margin"])<=.1]
                    for condition in ["E00","E10","E01","E11"]:
                        chosen=[r for r in selected if r["condition"]==condition];result=score(chosen)
                        near.append(dict(query=query,suite=suite,panel=panel,K=K,subset=subset,condition=condition,**result))
                        for metric,n,success in [("sensitivity",result["positives"],result["TP"]),("specificity",result["negatives"],result["TN"])]:
                            lo,hi=wilson(success,n) if n else (None,None)
                            intervals.append(dict(query=query,suite=suite,panel=panel,K=K,subset=subset,condition=condition,metric=metric,n=n,lower=lo,upper=hi,scope="Pointwise design-conditional Wilson; not ecological population"))
                    byworld=defaultdict(dict)
                    for r in selected:byworld[r["world"]][r["condition"]]=r
                    values=[];families=[]
                    for world,v in sorted(byworld.items()):
                        L={k:r["mass_absolute_error"] for k,r in v.items()};values.append([L["E00"]-L["E10"],L["E00"]-L["E01"],L["E11"]-L["E10"]-L["E01"]+L["E00"]]);families.append(v["E00"]["family"])
                    if not values:
                        for metric in ["state_gain","parameter_gain","interaction"]:
                            paired.append(dict(query=query,suite=suite,panel=panel,K=K,subset=subset,metric=metric,n=0,estimate=None,lower=None,upper=None,scope="UNDEFINED_NO_PRESENT_NEAR_CASES; empty registered stratum retained"))
                        continue
                    data=np.array(values);g=h.rng("bootstrap",100,qi,si,ki,subi);draws=[]
                    for iteration in range(1000):
                        index=np.concatenate([g.choice(np.flatnonzero(np.array(families)==f),sum(v==f for v in families),replace=True) for f in sorted(set(families))]);draws.append(data[index].mean(axis=0))
                    estimates=data.mean(axis=0);bounds=np.quantile(draws,[.025,.975],axis=0)
                    for j,metric in enumerate(["state_gain","parameter_gain","interaction"]):paired.append(dict(query=query,suite=suite,panel=panel,K=K,subset=subset,metric=metric,n=len(data),estimate=float(estimates[j]),lower=float(bounds[0,j]),upper=float(bounds[1,j]),scope="Paired-world family-stratified bootstrap, conditional fixed constructed design"))
    h.write_csv("oracle_attribution/near_boundary_scores.csv",near);h.write_csv("uncertainty/oracle_classification_intervals.csv",intervals);h.write_csv("uncertainty/oracle_paired_contrasts.csv",paired)
    # The fixed seed panel changes one source at a time relative to the same
    # original8-draw/1024-prefix baseline; no forced variance decomposition.
    index={(r["world"],r["condition"],r["query"],r["suite"]):r for r in parsed if r["panel"]=="PREFIX8" and r["K"]=="1024"};sources=[]
    for r in parsed:
        if not r["panel"].endswith("ONLY8"):continue
        base=index[(r["world"],r["condition"],r["query"],r["suite"])];name="estimated_viable" if r["query"]=="Q1" else "estimated_reserve"
        sources.append(dict(world=r["world"],variation_axis=r["panel"],replicate=r["replicate"],condition=r["condition"],query=r["query"],suite=r["suite"],mass_delta=r[name]-base[name],health_probability_delta=r["health_fraction"]-base["health_fraction"],
            decision_changed=boolean(r["predicted_health"])!=boolean(base["predicted_health"]),scope="One-source seed contrast; other packet/inference/future seed identities held as declared"))
    h.write_csv("oracle_attribution/seed_source_contrasts.csv",sources)
    print("Oracle uncertainty and one-source contrasts generated",flush=True)
if __name__=="__main__":
    import sys
    {"audit":audit,"uncertainty":uncertainty}[sys.argv[1]]()

