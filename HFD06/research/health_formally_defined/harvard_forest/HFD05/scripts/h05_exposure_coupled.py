"""Controlled individual-exposure reconstruction with captured original nuisance draws."""
import math
import numpy as np
import h05_common as h
from h05_exposure import records

def reconstruct(core,model,groups,generator,survivor_rng,family="B0"):
    wrapped=h.original.prior.MechanismRNG(generator,family) if family!="B0" else generator
    trace=h.original.prior.TraceRNG(wrapped);old=h.original.demography.reconstruct(core,model,trace);t=trace.trace
    if [len(t.get(k,[])) for k in ["choice","gamma","normal","binomial","poisson"]]!=[64,1,7,3,1]:raise ValueError("Original reconstruction trace contract changed")
    hazard=np.array(t["choice"],float).reshape(64);growth=np.clip(t["normal"][0],0,.7);wet=(np.arange(64)//4)%4==3
    unknown_h=hazard*np.exp(float(t["normal"][1])+wet*float(t["normal"][2]));shift=np.zeros(64)
    if family=="B1":
        fc=h.original.prior.PROTOCOL["missingness"][family];c=np.arange(64)
        shift=np.array(fc["taxon_log_shift"])[c//16]+np.array(fc["size_log_shift"])[c%4]+np.array(fc["stratum_log_shift"])[(c//4)%4]
    old_counts=[t["binomial"][0],t["binomial"][1],t["binomial"][2]];unseen=t["poisson"][0]
    names=["measured_alive","unresolved_fate","alive_missing_size"];bases=[model["diameter1"],model["diameter0"],model["diameter0"]]
    noise=[t["normal"][3],t["normal"][4],t["normal"][6]];hazards=[hazard,unknown_h*np.exp(shift),hazard]
    old_squares=[];new_counts=[];new_squares=[];analytic=[]
    for name,base,error,hz,original_n in zip(names,bases,noise,hazards,old_counts):
        N=np.zeros(64,int);S=np.zeros(64);O=np.zeros(64);A=np.zeros(64)
        for cell,times in enumerate(groups[name]):
            mean_time=float(times.mean()) if len(times) else 0
            pooled_d=np.clip(base[cell]+growth[cell]*mean_time+error[cell],1,300);O[cell]=original_n[cell]*pooled_d**2
            p=np.exp(-hz[cell]*times);alive=survivor_rng.binomial(1,p).astype(bool)
            d=np.clip(base[cell]+growth[cell]*times+error[cell],1,300)
            N[cell]=int(alive.sum());S[cell]=float((d[alive]**2).sum());A[cell]=float(p.sum())
        old_squares.append(O);new_counts.append(N);new_squares.append(S);analytic.append(A)
    branch_n=old["n"]-sum(old_counts)-unseen
    branch_square=old["n"]*old["d"]**2-sum(old_squares)-unseen*1.7**2
    if np.any(branch_n<0) or np.any(branch_square < -1e-7):raise ValueError("Captured original component ledger inconsistent")
    branch_square=np.maximum(branch_square,0)
    n=sum(new_counts)+unseen+branch_n;d2=sum(new_squares)+unseen*1.7**2+branch_square
    d=np.sqrt(np.divide(d2,n,out=model["diameter0"]**2,where=n>0))
    components=dict(old["components"])
    components.update(measured_projected_ba=float(new_squares[0].sum()*math.pi/40000),
        imputed_ba=float((new_squares[1].sum()+new_squares[2].sum()+(unseen*1.7**2).sum()+branch_square.sum())*math.pi/40000),
        wet_imputed_stems=int((new_counts[1]+new_counts[2]+unseen)[wet].sum()))
    new=dict(n=n,d=d,components=components,ancestry=old["ancestry"])
    ledger=dict(unchanged_unseen_count=int(unseen.sum()),unchanged_association_count=int(branch_n.sum()),
        conditional_individual_survivor_mean=float(sum(v.sum() for v in analytic)),conditional_original_survivor_mean=float(sum(np.exp(-hz[cell]*(times.mean() if len(times) else 0))*len(times) for hz,name in zip(hazards,names) for cell,times in enumerate(groups[name]))),
        original_dated_survivors=int(sum(v.sum() for v in old_counts)),individual_dated_survivors=int(sum(v.sum() for v in new_counts)))
    if ledger["conditional_individual_survivor_mean"]<ledger["conditional_original_survivor_mean"]-1e-8:raise ValueError("Conditional Jensen gap negative")
    return old,new,ledger

def run():
    h.guard("exposure_coupling");_,e0,e1,core,model=h.original.load_data();_,groups=records(e0,e1,core);c=h.config("exposure_coupling_protocol.json")["design"];rows=[];strata=[]
    designs=[("HFD02S_PRIMARY","B0",c["primary_draws"])] + [("HFD03_MISSINGNESS",family,c["missingness_draws"]) for family in c["families"]]
    for di,(design,family,total) in enumerate(designs):
        for i in range(total):
            g=h.original.frozen.Seeds().rng("reconstruct",i) if design=="HFD02S_PRIMARY" else h.original.prior.rng("state",i)
            old,new,ledger=reconstruct(core,model,groups,g,h.rng("individual_state",di,i),family)
            def values(state,mask=np.ones(64,bool)):
                n,d=state["n"][mask],state["d"][mask]
                return dict(N=int(n.sum()),BA=float((n*d*d).sum()*math.pi/40000),J=int((n*(d<10)).sum()))
            a,b=values(old),values(new)
            rows.append(dict(design=design,family=family,state=i,**{"original_"+k:v for k,v in a.items()},**{"individual_"+k:v for k,v in b.items()},**{k+"_delta":b[k]-a[k] for k in a},**ledger,
                original_measured_ba=old["components"]["measured_projected_ba"],individual_measured_ba=new["components"]["measured_projected_ba"],original_imputed_ba=old["components"]["imputed_ba"],individual_imputed_ba=new["components"]["imputed_ba"]))
            for axis,labels in [("taxon",np.arange(64)//16),("E0_sector_proxy",(np.arange(64)//4)%4)]:
                for label in range(4):
                    a,b=values(old,labels==label),values(new,labels==label)
                    strata.append(dict(design=design,family=family,state=i,axis=axis,label=label,**{"original_"+k:v for k,v in a.items()},**{"individual_"+k:v for k,v in b.items()},**{k+"_delta":b[k]-a[k] for k in a}))
            if i%32==0:print("Controlled exposure",design,family,i,flush=True)
    h.write_csv("exposure_model/controlled_origin_comparison.csv",rows);h.write_csv("exposure_model/controlled_stratum_comparison.csv",strata)
    summary=[]
    for design,family,_ in designs:
        selected=[r for r in rows if r["design"]==design and r["family"]==family]
        for metric in ["N","BA","J"]:
            d=np.array([r[metric+"_delta"] for r in selected])
            summary.append(dict(design=design,family=family,metric=metric,n=len(d),mean_delta=float(d.mean()),q95_abs_delta=float(np.quantile(abs(d),.95)),max_abs_delta=float(abs(d).max()),
                scope="Other nuisance draws fixed; dated survivor Monte Carlo remains; local extrema are observed finite-bank diagnostics"))
    h.write_csv("exposure_model/controlled_summary.csv",summary)
if __name__=="__main__":run()
