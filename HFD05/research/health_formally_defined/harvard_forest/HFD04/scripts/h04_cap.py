"""Isolated exact original fit/parameter/reconstruct code under four caps."""
import math
import numpy as np
import h04_common as h
import h03_calibration as cal
from h04_juvenile import visible_inputs
def run():
    h.confirm_guard()
    c=h.config()["effective_count"];sources,e0,e1,core,base=h.load_data();cfg=h.frozen.load("config/design.json")
    posterior_rows=[];state_rows=[];future_rows=[];calibration_rows=[]
    held,path=visible_inputs("M0");masked_core,_=cal.fit_visible(sources,path,held)
    cov=[x for x in cal.covariates(e0,e1) if x["stem_id"] in held]
    cells=np.array([x["cell"] for x in cov]);truth=np.array([h.frozen.fate(e1[x["stem_id"]]) for x in cov])
    dbh=np.array([h.empirical.number(e1[x["stem_id"]]["dbh"]) or np.nan for x in cov])
    live=(truth==1)&np.isfinite(dbh)
    for cap in c["caps"]:
        empirical,demographic=h.isolated_models(cap);model=empirical.fit(core)
        for cell,p in enumerate(model["hazard_posterior"]):
            cumulative=np.cumsum(p);grid=model["hazard_grid"]
            q=[float(grid[min(np.searchsorted(cumulative,x),255)]) for x in [.05,.5,.95]]
            posterior_rows.append(dict(cap=cap,cell=cell,n_fit=len(core["durations"][cell]),temper=min(1,cap/max(len(core["durations"][cell]),1)),
              hazard_mean=float((grid*p).sum()),hazard_p05=q[0],hazard_median=q[1],hazard_p95=q[2],hazard_width90=q[2]-q[0],
              growth_mean=float(model["growth_mean"][cell]),growth_sd=float(model["growth_sd"][cell]),growth_se=float(model["growth_se"][cell])))
        masked=empirical.fit(masked_core);prob=[];sizes=[]
        for draw in range(c["calibration_draws"]):
            trace=h.prior.TraceRNG(h.rng("cap_calibration",draw))
            demographic.reconstruct(masked_core,masked,trace)
            p,d=cal.endpoint_kernel(masked,cov,trace);prob.append(p);sizes.append(d)
        prob=np.array(prob);sizes=np.array(sizes)
        lo,hi=np.quantile(sizes[:,live],[.05,.95],axis=0)
        calibration_rows.append(dict(cap=cap,mask="M0",replicate=0,draws=c["calibration_draws"],n=len(cov),fate_brier=float(((prob.mean(axis=0)-truth)**2).mean()),
           DBH_MAE=float(np.abs(sizes[:,live].mean(axis=0)-dbh[live]).mean()),DBH_coverage90=float(((lo<=dbh[live])&(hi>=dbh[live])).mean()),
           J_prediction=float((prob*(sizes<10)).sum(axis=1).mean()),J_truth=int(((dbh>=1)&(dbh<10)&(truth==1)).sum()),visible_source_sha256=h.sha(path)))
        for i in range(c["states"]):
            state=demographic.reconstruct(core,model,h.rng("cap_state",i,0));params=demographic.parameters(model,h.rng("cap_state",i,1))
            state_rows.append(dict(cap=cap,state=i,N=int(state["n"].sum()),BA=float((state["n"]*state["d"]**2).sum()*math.pi/40000),
                                  J=int((state["n"]*(state["d"]<10)).sum())))
            for si,scenario in enumerate(cfg["scenarios"]):
                for vi,variant in enumerate(c["variants"]):
                    bank,lawful=demographic.simulate(state,params,scenario,max(c["horizons"]),c["paths"],h.rng("cap_future",i,si,vi),variant)
                    for horizon in c["horizons"]:
                        for rho in [.5,1]:
                            rr=demographic.response(bank,lawful,model,horizon,.75,.4,2,rho)
                            nv=int(rr["viable"].sum());nr=int(rr["reserve_viable"].sum())
                            for query in ["Q1","Q2"]:
                                for theta in [.5,.75,.9]:
                                    d=h.classify(rr["present"],[nv] if query=="Q1" else [nv,nr],c["paths"],theta)
                                    future_rows.append(dict(cap=cap,state=i,scenario=scenario["id"],variant=variant,horizon=horizon,rho=rho,query=query,theta=theta,
                                        present=rr["present"],nv=nv,nr=nr,K=c["paths"],FINITE_BANK_HEALTH=d["FINITE_BANK_HEALTH"],KERNEL_MC_STATUS=d["KERNEL_MC_STATUS"]))
            if i%8==0:print("Cap",cap,"state",i,flush=True)
    h.write_csv("effective_count/posteriors.csv",posterior_rows)
    h.write_csv("effective_count/origin_states.csv",state_rows)
    h.write_csv("effective_count/propagation.csv",future_rows)
    h.write_csv("effective_count/masking_calibration.csv",calibration_rows)
    h.write_json("effective_count/fit_contract.json",dict(cap_grid=c["caps"],reference=200,cap_selected=False,E2_used_in_fit=False,
       hazard_grid_points=256,prior="Beta annual survival32, Jacobian included; unequal observed durations; likelihood min(1,cap/n)",
       growth_se="max(.02, observedSD)/sqrt(min(max(n_growth,1),cap))",parameter_entry_shape="min(entrants+1,cap); positive pseudo-exposure1/6",shared_module_mutation=False))
if __name__=="__main__":run()
