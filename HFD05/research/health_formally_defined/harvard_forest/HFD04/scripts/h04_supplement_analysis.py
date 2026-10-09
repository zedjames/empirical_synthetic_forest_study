"""Postprocessing of separately registered context-prior audit; no tuning."""
import itertools,json
from collections import Counter,defaultdict
import numpy as np
import h04_common as h
import h04_analysis as a

def run():
    config=json.loads((h.ROOT/"config/context_prior_audit.json").read_text())
    record=json.loads((h.ROOT/"provenance/context_prior_registration.json").read_text())
    for name,digest in record["sources"].items():
        if h.sha(h.ROOT/name)!=digest:raise ValueError("Supplement registration changed")
    draws=defaultdict(list)
    for r in h.read_csv(h.ROOT/"context_ablation/prior_reconciled_draws.csv"):
        draws[int(r["world"])].append(r)
    rows=[]
    for ref in a.references():
        if ref["cohort"]!="boundary":continue
        draw=draws[int(ref["world"])];q=ref["query"]
        m=np.array([int(r["nv"])/int(r["K"]) for r in draw]);v=np.array([int(r["nr"])/int(r["K"]) for r in draw])
        p=np.mean([a.boolean(r[q]) for r in draw]).item()
        lo=np.mean([r[q+"_status"]=="TRUE" for r in draw]).item();hi=np.mean([r[q+"_status"]!="FALSE" for r in draw]).item()
        ml,mh=np.quantile(m,[.05,.95]);vl,vh=np.quantile(v,[.05,.95]);margin=float(ref["margin"])
        rows.append(dict(cohort="boundary",world=int(ref["world"]),suite=ref["suite"],query=q,policy="B2_roster_prior_v2",
            reference_present=a.boolean(ref["present"]),reference_health=a.boolean(ref["FINITE_BANK_HEALTH"]),reference_status=ref["KERNEL_MC_STATUS"],
            reference_viable=int(ref["nv"])/int(ref["K"]),reference_reserve=int(ref["nr"])/int(ref["K"]),
            margin=margin,signed_region=a.signed_region(margin),absolute_region=a.absolute_region(margin),failure_class=ref["failure_class"],
            estimated_viable=float(m.mean()),estimated_reserve=float(v.mean()),viable_lower=float(ml),viable_upper=float(mh),reserve_lower=float(vl),reserve_upper=float(vh),
            health_fraction=p,health_fraction_MC_lower=lo,health_fraction_MC_upper=hi,predicted_health=p>=.5,
            estimate_status="TRUE" if lo>.5 else ("FALSE" if hi<.5 else "MC_UNRESOLVED"),
            context_interpretation="PREDICTION_UNDER_HIDDEN_REALIZED_CONTEXT",observation_interpretation="ACTUAL_OBSERVATION_PACKET"))
    h.write_csv("context_ablation/prior_reconciled_performance.csv",rows)
    summary=[];strata=[]
    for suite,q in itertools.product(h.config()["boundary"]["misspecifications"],["Q1","Q2"]):
        selected=[r for r in rows if r["suite"]==suite and r["query"]==q]
        summary.append(dict(cohort="boundary",suite=suite,query=q,policy="B2_roster_prior_v2",**a.score(selected)))
        for region in a.REGIONS:
            strata.append(dict(suite=suite,query=q,kind="signed_margin",region=region,**a.score([r for r in selected if r["signed_region"]==region])))
        for region in ["0-.02",".02-.05",".05-.1",">.1"]:
            strata.append(dict(suite=suite,query=q,kind="absolute_margin",region=region,**a.score([r for r in selected if r["absolute_region"]==region])))
        for region in ["REALIZATION","Q1_CAPACITY","Q2_RESERVE","HEALTHY"]:
            strata.append(dict(suite=suite,query=q,kind="failure_class",region=region,**a.score([r for r in selected if r["failure_class"]==region])))
    h.write_csv("context_ablation/prior_reconciled_comparison.csv",summary)
    h.write_csv("context_ablation/prior_reconciled_strata.csv",strata)
    changes=[]
    for r in h.read_csv(h.ROOT/"juvenile_boundary/propagation.csv"):
        if a.boolean(r["coarse_Health"])==a.boolean(r["fine_Health"]):continue
        nv,nr,fv,fr=[int(r[k]) for k in ["coarse_nv","coarse_nr","fine_nv","fine_nr"]]
        changes.append(dict(**r,viable_change=(fv-nv)/256,reserve_change=(fr-nr)/256,
            coarse_reserve_given_viable=nr/nv if nv else None,fine_reserve_given_viable=fr/fv if fv else None,
            explanation="CONTINUATION_DRIVEN" if nr==nv and fr==fv else "NOT_SEPARATELY_IDENTIFIED",
            original_change_class_is_causal=False))
    h.write_csv("juvenile_boundary/changed_case_audit.csv",changes)
    print("Supplement scores",len(rows),"cells;",len(changes),"representation changes",flush=True)

if __name__=="__main__":run()
