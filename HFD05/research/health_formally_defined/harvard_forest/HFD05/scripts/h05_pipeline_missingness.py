"""Full B0-B4 original semantic surfaces under R0-R4, no original overwrites."""
import json,hashlib,math
import numpy as np
import h05_common as h
from h05_exposure import records
from h05_exposure_coupled import reconstruct
from h05_growth import simulate
from h05_bank_counts import count,infer
from h05_pipeline_primary import require_oracle

def run():
    h.guard("pipeline_missingness");require_oracle();import h03_missingness as historical
    _,e0,e1,core,model=h.original.load_data();_,groups=records(e0,e1,core);design=h.original.frozen.load("config/design.json");cells=historical.query_grid(design)
    if design["semantics"]["theta_grid"]!=[0,.5,.75,.9,1]:raise ValueError("Original present-at-zero index contract changed")
    N=128;K=256;shape=(4,N*2,len(cells));outputs=[];manifest=[];unit_rows=[]
    folder=h.ROOT/".runs"/("missingness_"+h.sha(h.ROOT/"config/pipeline_protocol.json")[:16]);folder.mkdir(parents=True,exist_ok=True)
    archived={(r["family"],r["scenario"],int(r["horizon"]),float(r["alpha"]),float(r["beta"]),int(r["tau"]),r["query"],float(r["rho"]),float(r["theta"])):r for r in h.read_csv(h.THREE/"missingness/phase_surface.csv")}
    for fi,family in enumerate(["B0","B1","B2","B3","B4"]):
        original_path=h.original.prior.RUN/("missingness_"+family+".npz");source_sha=h.sha(original_path)
        with np.load(original_path) as bank:
            r0=dict(nv=np.rint(bank["mass"]*K).astype(np.uint16),nr=np.rint(bank["reserve"]*K).astype(np.uint16),present=np.repeat(bank["point"][:,:,::5],5,axis=2),finite=bank["point"].copy(),status=bank["status"].copy())
        r1=dict(r0);r1["finite"],r1["status"]=infer(r0["nv"],r0["nr"],r0["present"],K,cells)
        if not np.array_equal(r0["finite"],r1["finite"]):raise ValueError("Simplification changed original finite labels")
        arrays={"R0":r0,"R1":r1}
        for result in ["R2","R3","R4"]:arrays[result]={key:np.empty(shape,dtype=np.uint16 if key in ["nv","nr"] else (np.uint8 if key=="status" else bool)) for key in r0}
        for i in range(N):
            a,b,ledger=reconstruct(core,model,groups,h.original.prior.rng("state",i),h.rng("individual_state",fi+1,i),family)
            params=h.original.demography.parameters(model,h.original.prior.rng("future",i,0))
            for vi,variant in enumerate(["gamma_growth","normal_growth"]):
                unit=i*2+vi
                for si,scenario in enumerate(design["scenarios"]):
                    for result,initial,corrected in [("R2",a,True),("R3",b,False),("R4",b,True)]:
                        index=(si,unit)
                        if result=="R2" and vi==1:
                            for key in r0:arrays[result][key][index]=r1[key][index]
                            continue
                        if result=="R4" and vi==1:
                            for key in r0:arrays[result][key][index]=arrays["R3"][key][index]
                            continue
                        bank,lawful=simulate(initial,params,scenario,20,K,h.original.prior.rng("future",i,vi,si+1),variant,corrected)
                        nv,nr,p=count(bank,lawful,model,cells);point,status=infer(nv,nr,p,K,cells)
                        for key,value in [("nv",nv),("nr",nr),("present",p),("finite",point),("status",status)]:arrays[result][key][index]=value
                        unit_rows.append(dict(family=family,state=i,variant=variant,scenario=scenario["id"],result=result,raw_history_sha256=hashlib.sha256(bank.tobytes()+lawful.tobytes()).hexdigest(),K=K,
                            source_identity="HFD05_CORRECTED_GENERATOR",inference_identity="HFD05_CORRECTED_INFERENCE",weights="Unconditional equally weighted candidate histories"))
            if i%16==0:print("Full corrected missingness",family,i,"of128",flush=True)
        for result,values in arrays.items():
            if result in ["R2","R3","R4"]:
                path=folder/(family+"_"+result+".npz")
                if path.exists():raise ValueError("Corrected sufficient bank already exists; explicit resume/replay path required")
                np.savez_compressed(path,**values)
                manifest.append(dict(family=family,result=result,sha256=h.sha(path),bytes=path.stat().st_size,private_path=str(path.relative_to(h.ROOT)),original_bank_sha256=source_sha))
            for si,scenario in enumerate(design["scenarios"]):
                for j,(horizon,alpha,beta,tau,query,rho,theta) in enumerate(cells):
                    point=values["finite"][si,:,j];status=values["status"][si,:,j];nv=values["nv"][si,:,j].astype(int);nr=values["nr"][si,:,j].astype(int);ov=r0["nv"][si,:,j].astype(int);orr=r0["nr"][si,:,j].astype(int)
                    if result=="R0":
                        old=archived[(family,scenario["id"],horizon,alpha,beta,tau,query,rho,theta)]
                        if float(point.mean())!=float(old["point_health_fraction"]) or not math.isclose(float(nv.mean()/K),float(old["viable_mass"]),abs_tol=1e-8):raise ValueError("Original missingness surface not reproduced")
                    outputs.append(dict(family=family,result=result,scenario=scenario["id"],horizon=horizon,alpha=alpha,beta=beta,tau=tau,query=query,rho=rho,theta=theta,outer_units=N*2,K=K,
                        viable_mass=float(nv.mean()/K),reserve_mass=float(nr.mean()/K),finite_health_fraction=float(point.mean()),present_fraction=float(values["present"][si,:,j].mean()),
                        kernel_TRUE_units=int((status==1).sum()),kernel_FALSE_units=int((status==0).sum()),kernel_MC_UNRESOLVED_units=int((status==2).sum()),
                        finite_health_disagreement_vs_original=float(np.mean(point!=r0["finite"][si,:,j])),mean_viable_delta_vs_original=float((nv-ov).mean()/K),mean_reserve_delta_vs_original=float((nr-orr).mean()/K),
                        local_max_abs_viable_delta=float(abs(nv-ov).max()/K),local_q95_abs_viable_delta=float(np.quantile(abs(nv-ov),.95)/K),local_max_abs_reserve_delta=float(abs(nr-orr).max()/K),local_q95_abs_reserve_delta=float(np.quantile(abs(nr-orr),.95)/K),
                        historical_status_preserved=result=="R0",case_weights="Original equal finite state/model design weights; not ecological posterior"))
        if h.sha(original_path)!=source_sha:raise ValueError("Original missingness bank changed during replay")
        h.write_csv("corrected_pipeline/missingness_surface_partial.csv",outputs);h.write_csv("corrected_pipeline/missingness_unit_provenance_partial.csv",unit_rows)
        h.write_json("corrected_pipeline/missingness_manifest_partial.json",dict(completed_families=fi+1,banks=manifest,scope="Partial execution checkpoint, never terminal evidence"))
    h.write_csv("corrected_pipeline/missingness_surface.csv",outputs);h.write_csv("corrected_pipeline/missingness_unit_provenance.csv",unit_rows)
    h.write_json("corrected_pipeline/missingness_bank_manifest.json",dict(status="COMPLETE",families=5,variants=5,query_cells_per_family_variant=4*len(cells),rows=len(outputs),banks=manifest,all_original_finite_surfaces_reproduced=True,old_banks_overwritten=False))
    print("Full missingness factorial complete",len(outputs),flush=True)
if __name__=="__main__":run()
