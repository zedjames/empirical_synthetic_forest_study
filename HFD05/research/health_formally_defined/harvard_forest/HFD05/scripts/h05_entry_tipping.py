"""Controlled individual-exposure reconstruction with captured original nuisance draws."""
import math
import numpy as np
import h05_common as h
from h05_exposure import records
from h05_growth import simulate
from h05_inference import classify

def reconstruct(core,model,groups,generator,survivor_rng,family="B0",diagnostic_log_shift=0):
    wrapped=h.original.prior.MechanismRNG(generator,family) if family!="B0" else generator
    trace=h.original.prior.TraceRNG(wrapped);old=h.original.demography.reconstruct(core,model,trace);t=trace.trace
    if [len(t.get(k,[])) for k in ["choice","gamma","normal","binomial","poisson"]]!=[64,1,7,3,1]:raise ValueError("Original reconstruction trace contract changed")
    hazard=np.array(t["choice"],float).reshape(64);growth=np.clip(t["normal"][0],0,.7);wet=(np.arange(64)//4)%4==3
    unknown_h=hazard*np.exp(float(t["normal"][1])+wet*float(t["normal"][2])+diagnostic_log_shift);shift=np.zeros(64)
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
    h.guard("entry_tipping");_,e0,e1,core,model=h.original.load_data();_,groups=records(e0,e1,core);design=h.original.frozen.load("config/design.json");rows=[]
    for i in range(16):
        params=h.original.demography.parameters(model,h.original.prior.rng("future",i,0))
        for shift in [-.5,0,.5]:
            old,state,ledger=reconstruct(core,model,groups,h.original.prior.rng("state",i),h.rng("tipping_state",i),"B0",shift)
            for weight in [0,.5,1]:
                p={k:v.copy() for k,v in params.items()};p["recruit"]*=weight
                for si in [0,3]:
                    bank,law=simulate(state,p,design["scenarios"][si],5,1024,h.rng("tipping_future",i,si),"gamma_growth",True)
                    rr=h.original.response(bank,law,model);nv=int(rr["viable"].sum());nr=int(rr["reserve_viable"].sum())
                    for query in ["Q1","Q2"]:
                        n=nv if query=="Q1" else nr
                        for theta in [0,.5,.75,.9,1]:
                            d=classify(rr["present"],nv,nr,1024,theta,query)
                            rows.append(dict(state=i,unknown_log_hazard_shift=shift,entry_membership_multiplier=weight,scenario=design["scenarios"][si]["id"],horizon=5,query=query,theta=theta,K=1024,nv=nv,nr=nr,present=rr["present"],mass=n/1024,finite_health=d["FINITE_BANK_HEALTH"],MC_status=d["KERNEL_MC_STATUS"],latent_fate_and_entry_weights="DIAGNOSTIC_NOT_POSTERIOR",realization_stratum="PRESENT" if rr["present"] else "ABSENT",margin=n/1024-theta,N=int(state["n"].sum()),J=int((state["n"]*(state["d"]<10)).sum()),measured_projected_ba=state["components"]["measured_projected_ba"],imputed_ba=state["components"]["imputed_ba"]))
        print("Entry/fate diagnostic state",i,"of16",flush=True)
    h.write_csv("entry_missingness/diagnostic_tipping.csv",rows)
    h.write_json("entry_missingness/diagnostic_completion.json",dict(status="COMPLETE",states=16,entry_weights=[0,.5,1],unknown_hazard_shifts=[-.5,0,.5],all_ranges_retained=True,weights_are_posterior=False,source_sha256=h.sha(__file__)))

if __name__=="__main__":run()
