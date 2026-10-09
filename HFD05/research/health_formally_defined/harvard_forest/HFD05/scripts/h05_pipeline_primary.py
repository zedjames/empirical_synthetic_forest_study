"""Complete original primary grid under separately addressed R0-R4 banks."""
import itertools,json,hashlib
from pathlib import Path
import numpy as np
import h05_common as h
from h05_growth import simulate
from h05_exposure import records
from h05_exposure_coupled import reconstruct
from h05_inference import wilson_array

def paths():
    manifest=json.loads((h.OLD/"provenance/run_manifest.json").read_text());old=h.OLD/".runs"/manifest["run_id"]
    for name in ["primary_histories.npy","primary_lawful.npy"]:
        if h.sha(old/name)!=manifest["outputs"][name]["sha256"]:raise ValueError("Original primary bank changed")
    run=h.ROOT/".runs"/("primary_"+h.sha(h.ROOT/"config/pipeline_protocol.json")[:16]);run.mkdir(parents=True,exist_ok=True)
    if run.resolve()==old.resolve() or h.ROOT.resolve() not in run.resolve().parents:raise ValueError("Unisolated corrected bank")
    return old,run
def require_oracle():
    path=h.ROOT/"oracle_attribution/analysis_summary.json"
    if not path.exists() or json.loads(path.read_text())["status"]!="EXECUTED_ORIGINAL_ROSTER":raise ValueError("Original oracle attribution must finish before corrected pipeline")
    h.guard("oracle")
def generate():
    h.guard("pipeline_primary");require_oracle();old,run=paths();original=np.load(old/"primary_histories.npy",mmap_mode="r");original_l=np.load(old/"primary_lawful.npy",mmap_mode="r")
    banks={};flags={}
    for identity in ["R2","R3","R4"]:
        path=run/(identity+"_histories.npy");lp=run/(identity+"_lawful.npy")
        if path.exists()!=lp.exists():raise ValueError("Partially allocated corrected pair requires explicit repair")
        banks[identity]=np.lib.format.open_memmap(path,mode="r+" if path.exists() else "w+",dtype=original.dtype,shape=original.shape)
        flags[identity]=np.lib.format.open_memmap(lp,mode="r+" if lp.exists() else "w+",dtype=original_l.dtype,shape=original_l.shape)
    progress=run/"progress.json";start=json.loads(progress.read_text())["next_state"] if progress.exists() else 0
    _,e0,e1,core,model=h.original.load_data();_,groups=records(e0,e1,core);seed=h.original.frozen.Seeds();design=h.original.frozen.load("config/design.json")
    for i in range(start,256):
        a,b,ledger=reconstruct(core,model,groups,seed.rng("reconstruct",i),h.rng("individual_state",0,i),"B0")
        for rep in range(2):
            params=h.original.demography.parameters(model,seed.rng("parameters",i,rep))
            for vi,variant in enumerate(["gamma_growth","normal_growth"]):
                for si,scenario in enumerate(design["scenarios"]):
                    index=(i,rep,vi,si)
                    if variant=="normal_growth":banks["R2"][index]=original[index];flags["R2"][index]=original_l[index]
                    else:
                        bank,lawful=simulate(a,params,scenario,20,64,seed.rng("future",i,rep,vi,si),variant,True)
                        banks["R2"][index]=bank;flags["R2"][index]=lawful
                    for identity,corrected in [("R3",False),("R4",True)]:
                        if identity=="R4" and variant=="normal_growth":bank,lawful=banks["R3"][index],flags["R3"][index]
                        else:bank,lawful=simulate(b,params,scenario,20,64,seed.rng("future",i,rep,vi,si),variant,corrected)
                        banks[identity][index]=bank;flags[identity][index]=lawful
        for array in list(banks.values())+list(flags.values()):array.flush()
        progress.write_text(json.dumps(dict(next_state=i+1,phase="pipeline_primary",source_sha256=h.sha(__file__),protocol_sha256=h.sha(h.ROOT/"config/pipeline_protocol.json")),sort_keys=True)+"\n")
        if i%8==0:print("Corrected full primary state",i,"of256",flush=True)
    h.write_json("corrected_pipeline/primary_bank_manifest.json",dict(identities={identity:{name:dict(sha256=h.sha(run/(identity+"_"+name+".npy")),bytes=(run/(identity+"_"+name+".npy")).stat().st_size,
        private_path=str((run/(identity+"_"+name+".npy")).relative_to(h.ROOT))) for name in ["histories","lawful"]} for identity in banks},
        original_source_archive="6b1fb4cff2149a58",old_banks_overwritten=False,all_states_retained=256,model_variants=2,parameter_replicates=2,scenarios=4,
        source_identity="HFD05_CORRECTED_GENERATOR",original_reference_identity="HFD02S_ORIGINAL_GENERATOR",controlled_exposure_source=h.sha(h.ROOT/"scripts/h05_exposure_coupled.py")))

def evaluate():
    h.guard("pipeline_primary");require_oracle();old,run=paths();manifest=json.loads((h.ROOT/"corrected_pipeline/primary_bank_manifest.json").read_text());design=h.original.frozen.load("config/design.json");model=json.loads((old/"fit.json").read_text());outputs=[]
    previous={(r["scenario"],int(r["horizon"]),float(r["alpha"]),float(r["beta"]),int(r["tau"]),r["query"],float(r["theta"]),float(r["reserve_ratio"]) if r["query"]=="Q2" else .5):r for r in h.read_csv(h.OLD/"health/health_phase.csv")}
    originals={}
    for identity in ["R0","R1","R2","R3","R4"]:
        location=old if identity in ["R0","R1"] else run
        bp=location/("primary_histories.npy" if identity in ["R0","R1"] else identity+"_histories.npy");lp=location/("primary_lawful.npy" if identity in ["R0","R1"] else identity+"_lawful.npy")
        if identity not in ["R0","R1"]:
            for path,name in [(bp,"histories"),(lp,"lawful")]:
                if h.sha(path)!=manifest["identities"][identity][name]["sha256"]:raise ValueError("Corrected bank changed")
        bank=np.load(bp,mmap_mode="r").reshape(-1,4,64,21,6);law=np.load(lp,mmap_mode="r").reshape(-1,4,64)
        structure=bank[:,:,:,:,1]/model["baseline_ba"];reg=bank[:,:,:,:,2]/model["baseline_juveniles"]
        reserve=np.divide(bank[:,:,:,:,2],bank[:,:,:,0:1,2],out=np.zeros_like(reg),where=bank[:,:,:,0:1,2]>0)
        for alpha,beta in itertools.product(design["semantics"]["alpha_grid"],design["semantics"]["beta_grid"]):
            realizes=(structure>=alpha)&(reg>=beta);present=realizes[:,:,0,0];longest=np.zeros(law.shape,int);current=longest.copy()
            for year in range(21):
                current=np.where(realizes[:,:,:,year],0,current+1);longest=np.maximum(longest,current)
                if year not in design["horizons_years"]:continue
                for tau,query,rho in itertools.product(design["semantics"]["recovery_years_grid"],["Q1","Q2"],design["semantics"]["reserve_ratio_grid"]):
                    if query=="Q1" and rho!=.5:continue
                    viable=law&realizes[:,:,:,year]&(longest<=tau);nv=viable.sum(axis=2);nr=(viable&(reserve[:,:,:,year]>=rho)).sum(axis=2)
                    if np.any(nr>nv):raise ValueError("Reserve subset failure")
                    counts=nv if query=="Q1" else nr
                    if identity=="R0" and query=="Q2":
                        lo,hi=h.original.prior.wilson(counts,64,.975)
                    else:lo,hi=wilson_array(counts,64)
                    for si,scenario in enumerate(design["scenarios"]):
                        for theta in design["semantics"]["theta_grid"]:
                            point=present[:,si]&(counts[:,si]/64>=theta);status=np.where(~present[:,si],0,np.where(theta==0,1,np.where(hi[:,si]<theta,0,np.where(lo[:,si]>theta,1,2))))
                            key=(scenario["id"],year,alpha,beta,tau,query,theta,rho)
                            if identity=="R0":
                                if float(point.mean())!=float(previous[key]["health_ensemble_fraction"]):raise ValueError("Original primary Health surface changed")
                                originals[key]=dict(point=point.copy(),nv=nv[:,si].copy(),nr=nr[:,si].copy())
                            original=originals[key]
                            if identity=="R1" and not np.array_equal(point,original["point"]):raise ValueError("Q2 mathematical simplification changed finite labels")
                            outputs.append(dict(result=identity,scenario=scenario["id"],horizon=year,alpha=alpha,beta=beta,tau=tau,query=query,rho=rho,theta=theta,outer_units=len(point),K=64,
                                viable_mass=float(nv[:,si].mean()/64),reserve_mass=float(nr[:,si].mean()/64),finite_health_fraction=float(point.mean()),present_fraction=float(present[:,si].mean()),
                                kernel_TRUE_units=int((status==1).sum()),kernel_FALSE_units=int((status==0).sum()),kernel_MC_UNRESOLVED_units=int((status==2).sum()),
                                finite_health_disagreement_vs_original=float(np.mean(point!=original["point"])),mean_viable_delta_vs_original=float(np.mean(nv[:,si]-original["nv"])/64),mean_reserve_delta_vs_original=float(np.mean(nr[:,si]-original["nr"])/64),
                                local_max_abs_viable_delta=float(np.max(abs(nv[:,si]-original["nv"]))/64),local_q95_abs_viable_delta=float(np.quantile(abs(nv[:,si]-original["nv"]),.95)/64),
                                local_max_abs_reserve_delta=float(np.max(abs(nr[:,si]-original["nr"]))/64),local_q95_abs_reserve_delta=float(np.quantile(abs(nr[:,si]-original["nr"]),.95)/64),
                                status_provenance="Retrospective procedure on immutable counts; not a rewritten historical emitted status" if identity=="R0" else "Corrected single-constraint inference",candidate_measure="Unconditional equal candidate-path weights"))
            print("Primary surfaces",identity,alpha,beta,flush=True)
    h.write_csv("corrected_pipeline/primary_surface.csv",outputs)
    h.write_json("corrected_pipeline/primary_summary.json",dict(query_cells_per_variant=6480,variants=5,rows=len(outputs),original_finite_surface_reproduced=True,
        R1_finite_exactly_preserved=True,all_original_corrected_references_addressable=True,finite_design_weights="Equal declared state/parameter/model units; no ecological posterior interpretation"))

if __name__=="__main__":
    import sys
    {"generate":generate,"evaluate":evaluate}[sys.argv[1]]()
