"""Original-roster oracle attribution, paired fixed prefixes, and source-separated replication."""
import itertools,json,math,hashlib
from collections import defaultdict
import numpy as np
import h05_common as h
import h05_oracle_engine as engine
from h05_inference import classify

def setup():
    _,_,_,core,model=h.original.load_data();worlds=json.loads((h.FOUR/"config/world_design.json").read_text())["worlds"]
    return core,model,worlds,dict(n=core["stats"]["n0"].copy(),d=model["diameter0"].copy())
def observation(state,world,replicate):
    if replicate==0:
        import h04_boundary
        return h04_boundary.make_packet(state,world)
    c=h.original.config()["estimator"];g=h.rng("oracle_observation",world,replicate)
    n=g.binomial(state["n"],c["observation_thinning"]);d=np.clip(state["d"]+g.normal(0,c["dbh_noise_sd"],64),1,300)
    return h.original.ObservationPacket(tuple(map(int,n)),tuple(float(v) if count else None for count,v in zip(n,d)),(c["observation_thinning"],)*64)

def generate():
    h.guard("oracle");c=h.config("oracle_protocol.json");core,model,worlds,anchor=setup();rows=[];states=[]
    import h04_boundary as boundary
    old={(int(r["world"]),int(r["draw"])):r for r in h.read_csv(h.FOUR/"context_ablation/estimator_draws.csv") if r["cohort"]=="boundary" and r["policy"]=="B0_full"}
    for world,cell in enumerate(worlds):
        true,params,_=boundary.world(core,model,cell,world);truth_state=engine.pack_state(true,engine.OracleState);truth_params=engine.pack_parameters(params,engine.OracleParameters)
        replicas=[0]+(c["replicate_panel"]["replicate_ids"] if world in c["replicate_panel"]["worlds"] else [])
        for rep in replicas:
            packet_rep=rep if rep in [1,2] else 0;inference_rep=rep-2 if rep in [3,4] else 0;future_rep=rep-4 if rep in [5,6] else 0
            packet=observation(true,world,packet_rep);marker=boundary.context("B0_full",cell["pressure"],cell["recovery"])
            draws=32 if rep==0 else c["replicate_panel"]["draws"]
            for draw in range(draws):
                sg=h.original.rng("estimate_state",world,draw) if inference_rep==0 else h.rng("oracle_state",world,inference_rep,draw)
                pg=h.original.rng("estimate_params",world,draw) if inference_rep==0 else h.rng("oracle_params",world,inference_rep,draw)
                estimated_state,estimated_params=engine.estimate(packet,marker,anchor,model,sg,pg)
                states.append(dict(world=world,replicate=rep,draw=draw,
                    N_error=sum(estimated_state.counts)-sum(truth_state.counts),BA_error=sum(n*d*d for n,d in zip(estimated_state.counts,estimated_state.diameters))*math.pi/40000-sum(n*d*d for n,d in zip(truth_state.counts,truth_state.diameters))*math.pi/40000,
                    J_error=sum(n for n,d in zip(estimated_state.counts,estimated_state.diameters) if d<10)-sum(n for n,d in zip(truth_state.counts,truth_state.diameters) if d<10)))
                for condition,state,pp in [("E00",estimated_state,estimated_params),("E10",truth_state,estimated_params),("E01",estimated_state,truth_params),("E11",truth_state,truth_params)]:
                    maximum=4096 if rep==0 and draw<8 else (1024 if rep else 256);flags=[];reserve=[];law=[];bank_hashes=[]
                    for chunk in range(maximum//256):
                        if future_rep==0:
                            g=h.original.rng("estimate_future",world,draw) if chunk==0 else h.rng("oracle_future",world,draw,chunk)
                        else:g=h.rng("oracle_future",world,future_rep,draw,chunk,1)
                        bank,lawful=engine.predict(condition,state,pp,256,g);rr=h.original.response(bank,lawful,model)
                        flags.extend(rr["viable"]);reserve.extend(rr["reserve_viable"]);law.extend(lawful);bank_hashes.append(hashlib.sha256(bank.tobytes()+lawful.tobytes()).hexdigest())
                        K=256*(chunk+1)
                        prefixes=[256,1024,4096] if rep==0 else [1024]
                        if K not in prefixes:continue
                        nv,nr=int(np.sum(flags)),int(np.sum(reserve))
                        if nr>nv or np.any(np.array(flags)&~np.array(law)):raise ValueError("Unconditional lawful subset violated")
                        if rep==0 and condition=="E00" and K==256:
                            reference=old[(world,draw)]
                            if (nv,nr,rr["present"])!=(int(reference["nv"]),int(reference["nr"]),reference["present"]=="True"):raise ValueError("Original E00 baseline changed")
                        for query in ["Q1","Q2"]:
                            d=classify(rr["present"],nv,nr,K,.75,query)
                            rows.append(dict(world=world,family=cell["family"],replicate=rep,draw=draw,condition=condition,query=query,K=K,nv=nv,nr=nr,present=rr["present"],finite_health=d["FINITE_BANK_HEALTH"],MC_status=d["KERNEL_MC_STATUS"],
                                adequacy=(nv if query=="Q1" else nr)/K>=.75,source_identity="HFD04_ORIGINAL_BOUNDARY",role_state=type(state).__name__,role_parameters=type(pp).__name__,
                                prefix_hash=hashlib.sha256("".join(bank_hashes).encode()).hexdigest(),finite_candidate_mass=1,success_renormalized=False,
                                variation_axis="BASELINE" if rep==0 else ("OBSERVATION_ONLY" if rep in [1,2] else ("INFERENCE_ONLY" if rep in [3,4] else "FUTURE_ONLY"))))
        if world%3==0:print("Original oracle world",world,"of",len(worlds),flush=True)
    h.write_csv("oracle_attribution/draw_census.csv",rows);h.write_csv("oracle_attribution/state_errors.csv",states)
    h.write_json("oracle_attribution/estimator_information.json",dict(status="STRICT_ROLE_INTERFACES",E00_oracle_inputs=False,
        engine_source_sha256=h.sha(engine.__file__),world_imports_in_engine=False,original_E00_all_2592_counts_reproduced=True,
        misspecification_removed_by_oracles=False,full_roster_retained=True,truth_reference="Original fixed 65536-path bank; finite and kernel-status identities separate"))

def analyze():
    h.guard("oracle");from h04_analysis import score,signed_region
    raw=h.read_csv(h.ROOT/"oracle_attribution/draw_census.csv");refs={(r["world"],r["suite"],r["query"]):r for r in h.read_csv(h.FOUR/"boundary_validation/truth_prefixes.csv") if r["K"]=="65536"}
    groups=defaultdict(list)
    for r in raw:
        for panel in ["BASELINE32","PREFIX8"] if r["replicate"]=="0" and r["K"]=="256" and int(r["draw"])<8 else (["BASELINE32"] if r["replicate"]=="0" and r["K"]=="256" else (["PREFIX8"] if r["replicate"]=="0" else [r["variation_axis"]+"8"])):
            groups[(r["world"],r["family"],r["replicate"],r["condition"],r["query"],r["K"],panel)].append(r)
    rows=[]
    for (world,family,rep,condition,query,K,panel),records in groups.items():
        v=np.array([int(r["nv"])/int(K) for r in records]);reserve=np.array([int(r["nr"])/int(K) for r in records]);mass=v if query=="Q1" else reserve
        health=np.mean([r["finite_health"]=="True" for r in records]);hl=np.mean([r["MC_status"]=="TRUE" for r in records]);hu=np.mean([r["MC_status"]!="FALSE" for r in records])
        vl,vu=np.quantile(v,[.05,.95]);rl,ru=np.quantile(reserve,[.05,.95]);pred_present=np.mean([r["present"]=="True" for r in records])>=.5;pred_adequate=np.mean([r["adequacy"]=="True" for r in records])>=.5
        for suite in h.original.config()["boundary"]["misspecifications"]:
            ref=refs[(world,suite,query)];nv,nr=int(ref["nv"]),int(ref["nr"]);present=ref["present"]=="True";target=(nv if query=="Q1" else nr)/65536;d=classify(present,nv,nr,65536,.75,query);margin=target-.75
            rows.append(dict(world=world,family=family,replicate=rep,condition=condition,query=query,K=int(K),draws=len(records),panel=panel,suite=suite,
                reference_present=present,reference_health=d["FINITE_BANK_HEALTH"],reference_status=d["KERNEL_MC_STATUS"],reference_viable=nv/65536,reference_reserve=nr/65536,
                margin=margin,signed_region=signed_region(margin),estimated_viable=float(v.mean()),estimated_reserve=float(reserve.mean()),viable_lower=float(vl),viable_upper=float(vu),reserve_lower=float(rl),reserve_upper=float(ru),
                health_fraction=float(health),predicted_health=health>=.5,estimate_status="TRUE" if hl>.5 else ("FALSE" if hu<.5 else "MC_UNRESOLVED"),health_MC_lower=float(hl),health_MC_upper=float(hu),
                realization_error=bool(pred_present!=present),adequacy_error=bool(pred_adequate!=(target>=.75)),mass_absolute_error=float(abs(mass.mean()-target)),mass_bias=float(mass.mean()-target),
                interval_interpretation={"E00":"State, parameter and future MC", "E10":"Parameter and future MC", "E01":"State and future MC", "E11":"Future MC only"}[condition]))
    h.write_csv("oracle_attribution/paired_world_losses.csv",rows)
    summaries=[];margins=[];panels=sorted({(r["panel"],r["replicate"],r["K"]) for r in rows})
    for (panel,rep,K),suite,query,condition in itertools.product(panels,h.original.config()["boundary"]["misspecifications"],["Q1","Q2"],["E00","E10","E01","E11"]):
        selected=[r for r in rows if (r["panel"],r["replicate"],r["K"],r["suite"],r["query"],r["condition"])==(panel,rep,K,suite,query,condition)]
        result=score(selected);summaries.append(dict(panel=panel,replicate=rep,K=K,suite=suite,query=query,condition=condition,**result,
            realization_error_fraction=float(np.mean([r["realization_error"] for r in selected])),adequacy_error_fraction=float(np.mean([r["adequacy_error"] for r in selected]))))
        for region in sorted({r["signed_region"] for r in rows}):
            subset=[r for r in selected if r["signed_region"]==region]
            margins.append(dict(panel=panel,replicate=rep,K=K,suite=suite,query=query,condition=condition,region=region,**score(subset)))
    h.write_csv("oracle_attribution/score_summary.csv",summaries);h.write_csv("oracle_attribution/signed_margin_scores.csv",margins)
    pairs=defaultdict(dict)
    for r in rows:pairs[(r["world"],r["panel"],r["replicate"],r["K"],r["suite"],r["query"])][r["condition"]]=r
    contrasts=[]
    for key,values in pairs.items():
        if set(values)!={"E00","E10","E01","E11"}:raise ValueError("Incomplete oracle factorial")
        L={k:v["mass_absolute_error"] for k,v in values.items()}
        contrasts.append(dict(world=key[0],panel=key[1],replicate=key[2],K=key[3],suite=key[4],query=key[5],state_gain=L["E00"]-L["E10"],parameter_gain=L["E00"]-L["E01"],
            interaction=L["E11"]-L["E10"]-L["E01"]+L["E00"],joint_gain=L["E00"]-L["E11"],no_additive_causal_partition=True))
    h.write_csv("oracle_attribution/interactions.csv",contrasts)
    convergence=[]
    for world,query,condition,suite in itertools.product(range(81),["Q1","Q2"],["E00","E10","E01","E11"],h.original.config()["boundary"]["misspecifications"]):
        panel=[r for r in rows if r["world"]==str(world) and r["panel"]=="PREFIX8" and r["query"]==query and r["condition"]==condition and r["suite"]==suite]
        if sorted(r["K"] for r in panel)!=[256,1024,4096]:raise ValueError("Fixed oracle prefixes absent")
        panel=sorted(panel,key=lambda r:r["K"]);a,b,c=panel
        convergence.append(dict(world=world,query=query,condition=condition,suite=suite,mass_delta_256_4096=c["estimated_viable" if query=="Q1" else "estimated_reserve"]-a["estimated_viable" if query=="Q1" else "estimated_reserve"],
            decision_changed_256_4096=a["predicted_health"]!=c["predicted_health"],MC_status_256=a["estimate_status"],MC_status_4096=c["estimate_status"],reference_unresolved=c["reference_status"]=="MC_UNRESOLVED",fixed_draw_units=8))
    h.write_csv("oracle_attribution/numerical_convergence.csv",convergence)
    h.write_json("oracle_attribution/analysis_summary.json",dict(status="EXECUTED_ORIGINAL_ROSTER",worlds=81,paired_world_rows=len(rows),complete_conditions=4,
        loss_contrasts="Signed nonadditive contrasts; model misspecification and finite future MC remain in E11",reference_alteration=False))
    print("Original oracle scored",len(rows),"paired world/query/suite cells",flush=True)

if __name__=="__main__":
    import sys
    {"generate":generate,"analyze":analyze}[sys.argv[1]]()
