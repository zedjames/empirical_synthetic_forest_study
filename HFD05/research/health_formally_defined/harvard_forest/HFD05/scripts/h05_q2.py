"""Count-only mathematical correction. Frozen bank bytes are read, never written."""
import itertools,json,math
from collections import Counter,defaultdict
from fractions import Fraction
import numpy as np
import h05_common as h
from h05_inference import classify,wilson_array

def boolean(v):return v is True or v=="True"
def add(rows,family,identity,nv,nr,K,present,theta,emitted=None,multiplicity=1):
    nv,nr,K=int(nv),int(nr),int(K);present=boolean(present);theta=float(theta)
    old=h.original.classify(present,[nv,nr],K,theta);new=classify(present,nv,nr,K,theta,"Q2")
    if old["FINITE_BANK_HEALTH"]!=new["FINITE_BANK_HEALTH"]:raise ValueError("Equivalent Q2 finite-bank disagreement")
    rows.append(dict(family=family,identity=str(identity),nv=nv,nr=nr,K=K,present=present,theta=theta,multiplicity=int(multiplicity),
        original_point=old["FINITE_BANK_HEALTH"],simplified_point=new["FINITE_BANK_HEALTH"],finite_identical=True,
        old_status=old["KERNEL_MC_STATUS"],corrected_status=new["KERNEL_MC_STATUS"],historical_emitted_status=emitted,
        old_lower=min(old["lower"]),old_upper=min(old["upper"]),corrected_lower=new["lower"],corrected_upper=new["upper"],
        old_width=min(old["upper"])-min(old["lower"]),corrected_width=new["upper"]-new["lower"],
        old_margin=min(nv,nr)/K-theta,corrected_margin=nr/K-theta,
        point_change_cause="NONE_EXACT_EQUIVALENCE",status_change_cause="SINGLE_RESERVE_WILSON95" if old["KERNEL_MC_STATUS"]!=new["KERNEL_MC_STATUS"] else "NONE"))

def finite_tests():
    cases=0
    for K in range(1,33):
        for nv in range(K+1):
            for nr in range(nv+1):
                for theta,present in itertools.product([0,.5,.75,.9,1],[False,True]):
                    t=Fraction(str(theta));a=present and Fraction(nv,K)>=t and Fraction(nr,K)>=t
                    if a!=classify(present,nv,nr,K,theta,"Q2")["FINITE_BANK_HEALTH"]:raise ValueError("Finite theorem counterexample")
                    cases+=1
    weighted=0
    for weights in itertools.product([0,1,2],repeat=3):
        denominator=max(sum(weights),1);w=[Fraction(v,denominator) for v in weights]
        for V,R in itertools.product(range(8),repeat=2):
            mv=sum(w[i] for i in range(3) if V>>i&1);mr=sum(w[i] for i in range(3) if V>>i&1 and R>>i&1)
            for theta in [Fraction(0),Fraction(1,2),Fraction(3,4),Fraction(1)]:
                if not (0<=mr<=mv<=1) or ((mv>=theta and mr>=theta)!=(mr>=theta)):raise ValueError("Weighted measure theorem failure")
                weighted+=1
    return dict(integer_truth_tables=cases,weighted_truth_tables=weighted,all_pass=True,scope="Finite fixtures accompany general nonnegative measure proof, not replace it")

def historical_counts():
    rows=[]
    paths=[("HFD04_BOUNDARY","boundary_validation/truth_prefixes.csv"),("HFD04_ESTIMATOR","context_ablation/estimator_draws.csv"),("HFD04_ROSTER_PRIOR","context_ablation/prior_reconciled_draws.csv")]
    for family,name in paths:
        for i,r in enumerate(h.read_csv(h.FOUR/name)):
            if r.get("query","Q2")!="Q2":continue
            add(rows,family,i,r["nv"],r["nr"],r["K"],r["present"],.75,r.get("KERNEL_MC_STATUS",r.get("Q2_status")))
    for i,r in enumerate(h.read_csv(h.FOUR/"decision_rules/hfd03_sufficient_census.csv")):
        if r["query"]=="Q2":add(rows,"HFD03_FULL_MISSINGNESS_CENSUS",i,r["nv"],r["nr"],r["K"],r["present"],r["theta"],r["old_status"],r["evaluations"])
    for i,r in enumerate(h.read_csv(h.THREE/"synthetic_validation/truth_census.csv")):
        add(rows,"HFD03_SYNTHETIC_TRUTH",i,round(float(r["truth_mass"])*4096),round(float(r["truth_reserve"])*4096),4096,r["truth_present"],.75,r["truth_mc_status"])
    for i,r in enumerate(h.read_csv(h.FOUR/"juvenile_boundary/propagation.csv")):
        for label in ["coarse","fine"]:add(rows,"HFD04_REPRESENTATION_"+label,str(i),r[label+"_nv"],r[label+"_nr"],r["K"],r[label+"_present"],r["theta"],r[label+"_status"])
    for i,r in enumerate(h.read_csv(h.FOUR/"effective_count/propagation.csv")):
        if r["query"]=="Q2":add(rows,"HFD04_CAP",i,r["nv"],r["nr"],r["K"],r["present"],r["theta"],r["KERNEL_MC_STATUS"])
    # Direct subset and archived prefix verification: no count relation assumed.
    import hashlib
    run=h.FOUR/".runs"/h.sha(h.FOUR/"config/protocol.json")[:16]/"truth"
    banks=json.loads((h.FOUR/"provenance/truth_bank_manifest.json").read_text())
    prefixes=defaultdict(list)
    for record in h.read_csv(h.FOUR/"boundary_validation/truth_prefixes.csv"):
        if record["query"]=="Q2":prefixes[(record["world"],record["suite"])].append(record)
    for item in banks:
        with np.load(run/(str(item["world"])+"_"+item["suite"]+".npz")) as b:
            v,r,l=b["viable"],b["reserve"],b["lawful"]
            if np.any(r&~v) or np.any(v&~l):raise ValueError("Reserve support is not a subset")
            digest=hashlib.sha256(v.tobytes()+r.tobytes()+l.tobytes()).hexdigest()
            if digest!=item["raw_response_sha256"]:raise ValueError("Original response bank mutated")
            for record in prefixes[(str(item["world"]),item["suite"])]:
                K=int(record["K"])
                if int(v[:K].sum())!=int(record["nv"]) or int(r[:K].sum())!=int(record["nr"]):raise ValueError("Prefix count/archive disagreement")
    return rows

def primary():
    manifest=json.loads((h.OLD/"provenance/run_manifest.json").read_text());run=h.OLD/".runs"/manifest["run_id"]
    for name in ["primary_histories.npy","primary_lawful.npy"]:
        if h.sha(run/name)!=manifest["outputs"][name]["sha256"]:raise ValueError("Original primary archive changed")
    bank=np.load(run/"primary_histories.npy",mmap_mode="r");law=np.load(run/"primary_lawful.npy",mmap_mode="r")
    c=h.original.frozen.load("config/design.json");model=json.loads((run/"fit.json").read_text())
    b=bank.reshape(-1,4,64,21,6);law=law.reshape(-1,4,64)
    structure=b[:,:,:,:,1]/model["baseline_ba"];reg=b[:,:,:,:,2]/model["baseline_juveniles"]
    reserve=np.divide(b[:,:,:,:,2],b[:,:,:,0:1,2],out=np.zeros_like(reg),where=b[:,:,:,0:1,2]>0)
    original={(r["scenario"],int(r["horizon"]),float(r["alpha"]),float(r["beta"]),int(r["tau"]),float(r["theta"]),float(r["reserve_ratio"])):r
        for r in h.read_csv(h.OLD/"health/health_phase.csv") if r["query"]=="Q2"}
    summaries=[];sufficient=[]
    for alpha,beta in itertools.product(c["semantics"]["alpha_grid"],c["semantics"]["beta_grid"]):
        realizes=(structure>=alpha)&(reg>=beta);present=realizes[:,:,0,0];longest=np.zeros(law.shape,int);current=longest.copy()
        for year in range(21):
            current=np.where(realizes[:,:,:,year],0,current+1);longest=np.maximum(longest,current)
            if year not in c["horizons_years"]:continue
            for tau,rho in itertools.product(c["semantics"]["recovery_years_grid"],c["semantics"]["reserve_ratio_grid"]):
                viable=law&realizes[:,:,:,year]&(longest<=tau);nv=viable.sum(axis=2);nr=(viable&(reserve[:,:,:,year]>=rho)).sum(axis=2)
                if np.any(nr>nv):raise ValueError("Primary reserve support violation")
                for s,scenario in enumerate(c["scenarios"]):
                    code=present[:,s].astype(int)*65*65+nv[:,s]*65+nr[:,s]
                    values,multiplicities=np.unique(code,return_counts=True)
                    for value,count in zip(values,multiplicities):
                        value=int(value);rr=value%65;value//=65;vv=value%65;pp=bool(value//65)
                        sufficient.append(dict(scenario=scenario["id"],horizon=year,alpha=alpha,beta=beta,tau=tau,rho=rho,nv=vv,nr=rr,K=64,present=pp,multiplicity=int(count)))
                    for theta in c["semantics"]["theta_grid"]:
                        oldpoint=present[:,s]&(nv[:,s]/64>=theta)&(nr[:,s]/64>=theta);newpoint=present[:,s]&(nr[:,s]/64>=theta)
                        if not np.array_equal(oldpoint,newpoint):raise ValueError("Primary exact finite equivalence defect")
                        old=h.original.prior.wilson(np.stack([nv[:,s],nr[:,s]]),64,.975);new=wilson_array(nr[:,s],64)
                        def statuses(lo,hi):
                            return np.where(~present[:,s],0,np.where(theta==0,1,np.where(hi<theta,0,np.where(lo>theta,1,2))))
                        os=statuses(np.min(old[0],axis=0),np.min(old[1],axis=0));ns=statuses(new[0],new[1])
                        key=(scenario["id"],year,alpha,beta,tau,theta,rho);reported=float(original[key]["health_ensemble_fraction"])
                        if oldpoint.mean()!=reported:raise ValueError("Primary original headline not reproduced")
                        summaries.append(dict(scenario=scenario["id"],horizon=year,alpha=alpha,beta=beta,tau=tau,theta=theta,rho=rho,outer_units=len(nv),
                            original_health_fraction=reported,corrected_health_fraction=float(newpoint.mean()),finite_identical=True,
                            original_procedure_unresolved=int((os==2).sum()),corrected_unresolved=int((ns==2).sum()),
                            old_mean_width=float(np.mean(np.min(old[1],axis=0)-np.min(old[0],axis=0))),corrected_mean_width=float(np.mean(new[1]-new[0]))))
        print("Q2 primary semantic cell",alpha,beta,flush=True)
    h.write_csv("q2_correction/primary_sufficient_counts.csv",sufficient)
    h.write_csv("q2_correction/primary_surface_identity.csv",summaries)
    return dict(original_primary_query_cells=len(summaries),all_finite_identical=True,source_bank_sha256=manifest["outputs"]["primary_histories.npy"]["sha256"])

def boundary_rescore():
    from h04_analysis import score,REGIONS
    draws=defaultdict(list)
    for r in h.read_csv(h.FOUR/"context_ablation/estimator_draws.csv"):
        status=classify(boolean(r["present"]),int(r["nv"]),int(r["nr"]),int(r["K"]),.75,"Q2")["KERNEL_MC_STATUS"]
        draws[(r["cohort"],r["world"],r["policy"])].append(status)
    refs={(r["world"],r["suite"]):r for r in h.read_csv(h.FOUR/"boundary_validation/truth_prefixes.csv") if r["query"]=="Q2" and r["K"]=="65536"}
    rows=[]
    for r in h.read_csv(h.FOUR/"boundary_validation/performance.csv"):
        if r["cohort"]!="boundary" or r["query"]!="Q2":continue
        row=dict(r);ref=refs[(r["world"],r["suite"])]
        revised=classify(boolean(ref["present"]),int(ref["nv"]),int(ref["nr"]),65536,.75,"Q2")
        statuses=draws[(r["cohort"],r["world"],r["policy"])];lo=sum(s=="TRUE" for s in statuses)/32;hi=sum(s!="FALSE" for s in statuses)/32
        row.update(reference_status=revised["KERNEL_MC_STATUS"],estimate_status="TRUE" if lo>.5 else ("FALSE" if hi<.5 else "MC_UNRESOLVED"),
            health_fraction_MC_lower=lo,health_fraction_MC_upper=hi,old_reference_status=r["reference_status"],reference_procedure="HFD05_CORRECTED_INFERENCE")
        for name in ["health_fraction","estimated_viable","estimated_reserve","viable_lower","viable_upper","reserve_lower","reserve_upper","reference_viable","reference_reserve","margin"]:row[name]=float(row[name])
        rows.append(row)
    h.write_csv("q2_correction/boundary_rescoring.csv",rows)
    summaries=[];strata=[];intervals=[]
    for suite,policy in itertools.product(h.original.config()["boundary"]["misspecifications"],h.original.config()["estimator"]["contexts"]):
        selected=[r for r in rows if r["suite"]==suite and r["policy"]==policy]
        subsets=[("ALL_ROSTER",selected),("PRESENT_NEAR_0.1",[r for r in selected if boolean(r["reference_present"]) and abs(r["margin"])<=.1])]
        for subset,records in subsets:
            result=score(records);summaries.append(dict(suite=suite,policy=policy,query="Q2",subset=subset,**result))
            certified=[r for r in records if r["reference_status"]!="MC_UNRESOLVED"]
            positive=[r for r in certified if boolean(r["reference_health"])];negative=[r for r in certified if not boolean(r["reference_health"])]
            for name,sample,inverted in [("sensitivity",positive,False),("specificity",negative,True)]:
                n=sum(not boolean(r["predicted_health"]) if inverted else boolean(r["predicted_health"]) for r in sample)
                lo,hi=__import__("h05_inference").wilson(n,len(sample)) if sample else (None,None)
                intervals.append(dict(suite=suite,policy=policy,subset=subset,metric=name,n=len(sample),lower=lo,upper=hi,scope="POINTWISE_DESIGN_CONDITIONAL_WILSON; finite constructed worlds, not ecological population"))
            if positive and negative:
                g=h.rng("bootstrap",list(h.original.config()["boundary"]["misspecifications"]).index(suite),list(h.original.config()["estimator"]["contexts"]).index(policy),0 if subset=="ALL_ROSTER" else 1)
                ps=np.array([boolean(r["predicted_health"]) for r in positive]);ns=np.array([not boolean(r["predicted_health"]) for r in negative])
                samples=.5*(g.binomial(len(ps),ps.mean(),1000)/len(ps)+g.binomial(len(ns),ns.mean(),1000)/len(ns));lo,hi=np.quantile(samples,[.025,.975])
                intervals.append(dict(suite=suite,policy=policy,subset=subset,metric="balanced_accuracy",n=len(certified),lower=float(lo),upper=float(hi),scope="CLASS_STRATIFIED_PARAMETRIC_BOOTSTRAP; conditional on certified constructed-world outcomes"))
        for region in REGIONS:
            strata.append(dict(suite=suite,policy=policy,region=region,**score([r for r in selected if r["signed_region"]==region])))
    h.write_csv("q2_correction/boundary_score_summary.csv",summaries);h.write_csv("q2_correction/boundary_margin_strata.csv",strata)
    h.write_csv("q2_correction/design_conditional_intervals.csv",intervals)
    ambiguity=[]
    for row in rows:
        if row["suite"]=="correct" and row["policy"]=="B0_full" and row["reference_present"]=="True" and abs(row["margin"])<=.1 and row["reference_status"]=="MC_UNRESOLVED":
            ambiguity.append(dict(world=row["world"],predicted_health=row["predicted_health"],health_fraction=row["health_fraction"],reference_status=row["reference_status"],old_reference_status=row["old_reference_status"]))
    # Every unresolved assignment retained; not a chosen favorable truth roster.
    selected=[r for r in rows if r["suite"]=="correct" and r["policy"]=="B0_full" and boolean(r["reference_present"]) and abs(r["margin"])<=.1]
    uncertain=[r for r in selected if r["reference_status"]=="MC_UNRESOLVED"]
    if len(uncertain)>16:raise ValueError("Ambiguity enumeration needs separately registered larger design")
    evaluated=[]
    for assignment in itertools.product([False,True],repeat=len(uncertain)):
        copied=[dict(r) for r in selected];j=0
        for r in copied:
            if r["reference_status"]=="MC_UNRESOLVED":r["reference_health"]=assignment[j];r["reference_status"]="TRUE" if assignment[j] else "FALSE";j+=1
        result=score(copied);evaluated.append(result)
    h.write_json("q2_correction/reference_ambiguity.json",dict(original_unresolved=6,corrected_unresolved=len(uncertain),assignments=len(evaluated),
        uncertainty_is_reference_truth_not_estimator_sampling=True,balanced_accuracy_min=min(r["balanced_accuracy"] for r in evaluated),balanced_accuracy_max=max(r["balanced_accuracy"] for r in evaluated),
        uncertain_worlds=[r["world"] for r in uncertain],scope="All binary assignments of uncertified constructed-world truths; descriptive identification bracket"))

def run():
    h.guard("q2");tests=finite_tests();rows=historical_counts()
    h.write_csv("q2_correction/finite_bank_identity.csv",rows)
    summary=[]
    for family in sorted({r["family"] for r in rows}):
        selected=[r for r in rows if r["family"]==family];total=sum(r["multiplicity"] for r in selected)
        summary.append(dict(family=family,sufficient_cells=len(selected),evaluations=total,finite_disagreements=0,
            old_unresolved=sum(r["multiplicity"] for r in selected if r["old_status"]=="MC_UNRESOLVED"),
            corrected_unresolved=sum(r["multiplicity"] for r in selected if r["corrected_status"]=="MC_UNRESOLVED"),
            status_changes=sum(r["multiplicity"] for r in selected if r["old_status"]!=r["corrected_status"]),
            mean_old_width=sum(r["old_width"]*r["multiplicity"] for r in selected)/total,mean_corrected_width=sum(r["corrected_width"]*r["multiplicity"] for r in selected)/total))
    h.write_csv("q2_correction/resolution_comparison.csv",summary)
    primary_result=primary();boundary_rescore()
    h.write_json("mathematical_audit/q2_verification.json",dict(status="PASS",tests=tests,primary=primary_result,archived_count_cells=len(rows),raw_boundary_subset_banks=243,
        old_point_labels_preserved=True,corrected_intervals="SINGLE_RESERVE_WILSON95",health_definition_changed=False))
    print("PASS Q2 mathematical correction",len(rows),"archived count cells",flush=True)

if __name__=="__main__":run()
