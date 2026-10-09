"""Derived decision scores; no truth generation or estimator parameter changes."""
from collections import defaultdict
import itertools, json
import numpy as np
import h04_common as h
def boolean(value):return value is True or value=="True"
def signed_region(m):
    if m<-.1:return "OUTSIDE_BELOW"
    if m<-.05:return "CLEARLY_BELOW"
    if m<-.02:return "MODERATELY_BELOW"
    if m<0:return "IMMEDIATELY_BELOW"
    if m<=.02:return "IMMEDIATELY_ABOVE"
    if m<=.05:return "MODERATELY_ABOVE"
    if m<=.1:return "CLEARLY_ABOVE"
    return "OUTSIDE_ABOVE"
REGIONS=["OUTSIDE_BELOW","CLEARLY_BELOW","MODERATELY_BELOW","IMMEDIATELY_BELOW","IMMEDIATELY_ABOVE","MODERATELY_ABOVE","CLEARLY_ABOVE","OUTSIDE_ABOVE"]
def absolute_region(m):
    a=abs(m)
    return "0-.02" if a<=.02 else (".02-.05" if a<=.05 else (".05-.1" if a<=.1 else ">.1"))
def references():
    rows=[]
    for row in h.read_csv(h.ROOT/"boundary_validation/truth_prefixes.csv"):
        if int(row["K"])==65536:rows.append(dict(cohort="boundary",**row))
    for r in h.read_csv(h.THREE/"synthetic_validation/truth_census.csv"):
        nv=round(float(r["truth_mass"])*4096);nr=round(float(r["truth_reserve"])*4096);present=boolean(r["truth_present"])
        for q in ["Q1","Q2"]:
            counts=[nv] if q=="Q1" else [nv,nr];d=h.classify(present,counts,4096,.75)
            rows.append(dict(cohort="HFD03",world=r["world"],family="HFD03_factorial",suite=r["suite"],query=q,K=4096,nv=nv,nr=nr,present=present,
               FINITE_BANK_HEALTH=d["FINITE_BANK_HEALTH"],KERNEL_MC_STATUS=d["KERNEL_MC_STATUS"],margin=min(counts)/4096-.75,
               failure_class="REALIZATION" if not present else ("Q1_CAPACITY" if nv/4096<.75 else ("Q2_RESERVE" if q=="Q2" and nr/4096<.75 else "HEALTHY"))))
    return rows
def score(rows):
    certified=[r for r in rows if r["reference_status"]!="MC_UNRESOLVED"]
    y=np.array([boolean(r["reference_health"]) for r in certified],bool)
    pred=np.array([boolean(r["predicted_health"]) for r in certified],bool)
    p=np.array([r["health_fraction"] for r in certified])
    tp=int((y&pred).sum());tn=int((~y&~pred).sum());fp=int((~y&pred).sum());fn=int((y&~pred).sum())
    pos=int(y.sum());neg=int((~y).sum())
    sens=tp/pos if pos else None;spec=tn/neg if neg else None
    result=dict(n=len(rows),certified_n=len(certified),positives=pos,negatives=neg,TP=tp,TN=tn,FP=fp,FN=fn,
      sensitivity=sens,specificity=spec,FPR=fp/neg if neg else None,FNR=fn/pos if pos else None,
      balanced_accuracy=(sens+spec)/2 if sens is not None and spec is not None else None,
      brier=float(((p-y)**2).mean()) if len(y) else None,
      certified_error=(fp+fn)/len(y) if len(y) else None,
      reference_unresolved_fraction=(len(rows)-len(certified))/len(rows) if rows else None,
      estimator_unresolved_fraction=sum(r["estimate_status"]=="MC_UNRESOLVED" for r in rows)/len(rows) if rows else None)
    for name in ["viable","reserve"]:
        errors=np.array([r["estimated_"+name]-r["reference_"+name] for r in rows])
        result[name+"_bias"]=float(errors.mean()) if len(errors) else None
        result[name+"_MAE"]=float(np.abs(errors).mean()) if len(errors) else None
        result[name+"_coverage90"]=sum(r[name+"_lower"]<=r["reference_"+name]<=r[name+"_upper"] for r in rows)/len(rows) if rows else None
        result[name+"_width90"]=np.mean([r[name+"_upper"]-r[name+"_lower"] for r in rows]).item() if rows else None
    return result
def analyze():
    grouped=defaultdict(list)
    for r in h.read_csv(h.ROOT/"context_ablation/estimator_draws.csv"):
        grouped[(r["cohort"],r["world"],r["policy"])].append(r)
    refs=references();rows=[]
    for ref in refs:
        for policy in h.config()["estimator"]["contexts"]:
            draw=grouped[(ref["cohort"],str(ref["world"]),policy)]
            mass=np.array([int(r["nv"])/int(r["K"]) for r in draw]);reserve=np.array([int(r["nr"])/int(r["K"]) for r in draw])
            q=ref["query"];p=np.mean([boolean(r[q]) for r in draw]).item()
            hl=np.mean([r[q+"_status"]=="TRUE" for r in draw]).item();hu=np.mean([r[q+"_status"]!="FALSE" for r in draw]).item()
            status="TRUE" if hl>.5 else ("FALSE" if hu<.5 else "MC_UNRESOLVED")
            lo,hi=np.quantile(mass,[.05,.95]);rl,rh=np.quantile(reserve,[.05,.95])
            margin=float(ref["margin"])
            rows.append(dict(cohort=ref["cohort"],world=int(ref["world"]),family=ref["family"],suite=ref["suite"],query=q,policy=policy,
              reference_present=boolean(ref["present"]),reference_health=boolean(ref["FINITE_BANK_HEALTH"]),reference_status=ref["KERNEL_MC_STATUS"],
              reference_viable=int(ref["nv"])/int(ref["K"]),reference_reserve=int(ref["nr"])/int(ref["K"]),reference_K=int(ref["K"]),
              margin=margin,signed_region=signed_region(margin),absolute_region=absolute_region(margin),failure_class=ref["failure_class"],
              estimated_viable=float(mass.mean()),estimated_reserve=float(reserve.mean()),viable_lower=float(lo),viable_upper=float(hi),reserve_lower=float(rl),reserve_upper=float(rh),
              health_fraction=p,health_fraction_MC_lower=hl,health_fraction_MC_upper=hu,predicted_health=p>=.5,estimate_status=status,
              context_interpretation="KNOWN_CONDITIONAL_CONTEXT" if policy in ["B0_full","B3_label_only"] else "PREDICTION_UNDER_HIDDEN_REALIZED_CONTEXT",
              observation_interpretation="PRIOR_STATE_ONLY" if policy=="B3_label_only" else "ACTUAL_OBSERVATION_PACKET"))
    h.write_csv("boundary_validation/performance.csv",rows)
    summary=[];byregion=[];byabsolute=[];byfailure=[]
    cfg=h.config()["boundary"]
    for cohort,suite,q,policy in itertools.product(["boundary","HFD03"],cfg["misspecifications"],["Q1","Q2"],h.config()["estimator"]["contexts"]):
        selected=[r for r in rows if (r["cohort"],r["suite"],r["query"],r["policy"])==(cohort,suite,q,policy)]
        summary.append(dict(cohort=cohort,suite=suite,query=q,policy=policy,**score(selected)))
        for region in REGIONS:
            byregion.append(dict(cohort=cohort,suite=suite,query=q,policy=policy,region=region,**score([r for r in selected if r["signed_region"]==region])))
        for region in ["0-.02",".02-.05",".05-.1",">.1"]:
            byabsolute.append(dict(cohort=cohort,suite=suite,query=q,policy=policy,region=region,**score([r for r in selected if r["absolute_region"]==region])))
        for failure in ["REALIZATION","Q1_CAPACITY","Q2_RESERVE","HEALTHY"]:
            byfailure.append(dict(cohort=cohort,suite=suite,query=q,policy=policy,failure_class=failure,**score([r for r in selected if r["failure_class"]==failure])))
    baseline={(r["cohort"],r["suite"],r["query"]):r for r in summary if r["policy"]=="B0_full"}
    for row in summary:
        base=baseline[(row["cohort"],row["suite"],row["query"])]
        for metric in ["balanced_accuracy","viable_MAE","reserve_MAE","brier","estimator_unresolved_fraction"]:
            row["delta_vs_B0_"+metric]=row[metric]-base[metric] if row[metric] is not None and base[metric] is not None else None
    h.write_csv("context_ablation/comparison.csv",summary)
    h.write_csv("boundary_validation/signed_margin_scores.csv",byregion)
    h.write_csv("boundary_validation/absolute_margin_scores.csv",byabsolute)
    h.write_csv("boundary_validation/failure_scores.csv",byfailure)
    primary={}
    for q in ["Q1","Q2"]:
        selected=[r for r in rows if r["cohort"]=="boundary" and r["suite"]=="correct" and r["query"]==q and r["policy"]=="B0_full" and r["reference_present"] and abs(r["margin"])<=.1]
        s=score(selected)
        good=s["positives"]>=8 and s["negatives"]>=8 and s["sensitivity"] is not None and s["specificity"] is not None and min(s["sensitivity"],s["specificity"])>=.8
        s["status"]="ESTABLISHED" if good else ("PARTIAL" if s["positives"] and s["negatives"] else "NOT_ESTABLISHED");primary[q]=s
    h.write_json("boundary_validation/primary_discrimination.json",primary)
    print("Derived",len(rows),"world/suite/query/context comparisons",flush=True)
def old_decisions():
    import h03_missingness as old
    cells=old.query_grid(h.frozen.load("config/design.json"));records=[];results=[]
    for family in ["B0","B1","B2","B3","B4"]:
        bank=np.load(h.prior.RUN/("missingness_"+family+".npz"))
        nv=np.rint(bank["mass"]*256).astype(int);nr=np.rint(bank["reserve"]*256).astype(int)
        present=np.repeat(bank["point"][:,:,::5],5,axis=2)
        for theta in [0,.5,.75,.9,1]:
            for query in ["Q1","Q2"]:
                indices=[j for j,c in enumerate(cells) if c[-1]==theta and c[4]==query]
                v=nv[:,:,indices].ravel();r=nr[:,:,indices].ravel();p=present[:,:,indices].ravel();z=bank["status"][:,:,indices].ravel()
                code=(((p.astype(int)*257+v)*257+r)*3+z)
                unique,counts=np.unique(code,return_counts=True);aggregate=defaultdict(int)
                for value,count in zip(unique,counts):
                    value=int(value);old_status=value%3;value//=3;rr=value%257;value//=257;vv=value%257;pp=bool(value//257)
                    decision=h.classify(pp,[vv] if query=="Q1" else [vv,rr],256,theta)
                    aggregate[decision["KERNEL_MC_STATUS"]]+=int(count)
                    records.append(dict(family=family,query=query,theta=theta,K=256,nv=vv,nr=rr,present=pp,old_status=["FALSE","TRUE","MC_UNRESOLVED"][old_status],evaluations=int(count),
                      FINITE_BANK_HEALTH=decision["FINITE_BANK_HEALTH"],KERNEL_MC_STATUS=decision["KERNEL_MC_STATUS"]))
                results.append(dict(family=family,query=query,theta=theta,evaluations=len(v),TRUE=aggregate["TRUE"],FALSE=aggregate["FALSE"],MC_UNRESOLVED=aggregate["MC_UNRESOLVED"]))
    h.write_csv("decision_rules/hfd03_sufficient_census.csv",records);h.write_csv("decision_rules/versioned_resolution.csv",results)
    h.write_json("decision_rules/interpretation.json",dict(version="HFD04-derived-status-v1",source="Immutable HFD03 missingness banks",
        historical_unresolved_total=963976,historical_interior_unresolved=3984,theta0="Adequacy algebraic; Health present-dependent",
        theta1="Finite bank all success may be TRUE; continuous-law massone never certified by all-success sample",
        inclusive_finite_rule=True,strict_MC_rule=True,renormalized_successes=False))
    print("Versioned exact/kernel decisions",len(records),"sufficient census rows",flush=True)
if __name__=="__main__":
    import sys
    old_decisions() if sys.argv[1]=="decisions" else analyze()

