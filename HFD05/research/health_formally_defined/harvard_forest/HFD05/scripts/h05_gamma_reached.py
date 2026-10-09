"""Predeclared impact panel for every source-detected reached primary identity."""
import numpy as np
import h05_common as h
from h05_growth import simulate,laws
from h05_inference import classify
def run():
    h.guard("gamma_reached");c=h.config("gamma_reached_protocol.json");_,_,_,core,model=h.original.load_data();seed=h.original.frozen.Seeds()
    detected=sorted({int(r["parameter_identity"]) for r in h.read_csv(h.ROOT/"gamma_growth/branch_census.csv") if r["family"]=="HFD02S_PRIMARY" and r["positive_shape_floor"]=="True"})
    if detected!=c["parameter_identities"]:raise ValueError("Reached branch roster changed")
    rows=[];clipping=[]
    for identity in detected:
        i,rep=divmod(identity,2);state=h.original.demography.reconstruct(core,model,seed.rng("reconstruct",i));params=h.original.demography.parameters(model,seed.rng("parameters",i,rep))
        for si,scenario in enumerate(h.original.frozen.load("config/design.json")["scenarios"]):
            outputs=[]
            for corrected in [False,True]:
                g=h.rng("gamma_impact",2,identity,si);bank,lawful=simulate(state,params,scenario,20,c["paths"],g,corrected=corrected)
                outputs.append((bank,lawful,g.bit_generator.state))
            if outputs[0][2]!=outputs[1][2]:raise ValueError("Declared Gamma-only coupling changed RNG consumption")
            for horizon in c["horizons"]:
                a=[]
                for bank,lawful,_ in outputs:
                    rr=h.original.demography.response(bank,lawful,model,horizon,.75,.4,2,1);nv=int(rr["viable"].sum());nr=int(rr["reserve_viable"].sum());d=classify(rr["present"],nv,nr,c["paths"],.75,"Q2")
                    a.append(dict(nv=nv,nr=nr,present=rr["present"],health=d["FINITE_BANK_HEALTH"],status=d["KERNEL_MC_STATUS"],mean_BA=float(bank[:,horizon,1].mean()),mean_J=float(bank[:,horizon,2].mean())))
                old,new=a
                rows.append(dict(parameter_identity=identity,state=i,replicate=rep,scenario=scenario["id"],horizon=horizon,K=c["paths"],
                    **{"original_"+k:v for k,v in old.items()},**{"corrected_"+k:v for k,v in new.items()},viable_delta=(new["nv"]-old["nv"])/c["paths"],reserve_delta=(new["nr"]-old["nr"])/c["paths"],health_changed=new["health"]!=old["health"],same_rng_terminal_state=True))
            print("Reached Gamma panel",identity,scenario["id"],flush=True)
    h.write_csv("gamma_growth/reached_paired_impact.csv",rows)
    h.write_json("gamma_growth/reached_impact_summary.json",dict(cells=len(rows),all_reached_primary_identities=True,max_abs_viable_delta=max(abs(r["viable_delta"]) for r in rows),max_abs_reserve_delta=max(abs(r["reserve_delta"]) for r in rows),health_changes=sum(r["health_changed"] for r in rows),
        conditional_moment_defect="REACHED_MATERIAL_MOMENTS",scope="Registered all-reached-primary-identity panel; not a full surface or ecological robustness guarantee"))
if __name__=="__main__":run()
