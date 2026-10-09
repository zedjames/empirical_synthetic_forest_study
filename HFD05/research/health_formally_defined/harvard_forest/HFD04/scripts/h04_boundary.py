"""Frozen whole-grid truth census, disjoint observation-only estimation."""
import hashlib, itertools, json, math
from statistics import NormalDist
import numpy as np
import h04_common as h
from h04_pilot_v2 import size_boundary_state
from h04_estimator import Context,evaluate
def moments(state,params,years=5):
    n,d=state["n"],state["d"];g=params["growth"];sd=params["sd"];haz=params["hazard"]
    total=float((n*np.exp(-years*haz)*(d*d+2*d*years*g+(years*g)**2+years*sd**2)).sum())
    entry=params["recruit"].reshape(16,4).sum(axis=1);idx=np.arange(16)*4
    for age in range(years):
        total+=float((entry*np.exp(-age*haz[idx])*((1.7+age*g[idx])**2+age*sd[idx]**2)).sum())
    return total*math.pi/40000
def expected_J(state,params,jd,juvenile_mask,years=5):
    j=juvenile_mask;n=state["n"][j];g=params["growth"][j];sd=params["sd"][j];haz=params["hazard"][j]
    cdf=np.array([NormalDist().cdf(float(z)) for z in (10-jd-years*g)/np.maximum(np.sqrt(years)*sd,1e-8)])
    result=float((n*np.exp(-years*haz)*cdf).sum())
    entry=params["recruit"].reshape(16,4).sum(axis=1);idx=np.arange(16)*4
    for age in range(years):
        result+=float((entry*np.exp(-age*params["hazard"][idx])).sum())
    return result
def world(core,model,cell,identity):
    params=h.parameters(model,cell["pressure"],cell["recovery"],h.rng("world_parameters",identity))
    state=size_boundary_state(core,model,cell["structure"],cell["regeneration"],cell.get("juvenile_diameter"))
    meanparams=h.parameters(model,cell["pressure"],cell["recovery"],h.rng("world_state",identity),means=True)
    normalized=False
    if cell["family"]=="Q1":
        target=moments(state,meanparams)
        j=model["diameter0"]<10;d=state["d"].copy()
        # Root of deterministic moment map only; not a simulation/truth score.
        lo,hi=.5,2
        for step in range(48):
            ratio=(lo+hi)/2;state["d"][~j]=d[~j]*ratio
            if moments(state,params)<target:lo=ratio
            else:hi=ratio
        normalized=True
    elif cell["family"]=="Q2":
        juvenile_mask=model["diameter0"]<10
        target=expected_J(state,meanparams,cell["juvenile_diameter"],juvenile_mask)
        lo,hi=1,9.99
        for step in range(48):
            jd=(lo+hi)/2
            if expected_J(state,params,jd,juvenile_mask)>target:lo=jd
            else:hi=jd
        state=size_boundary_state(core,model,cell["structure"],cell["regeneration"],jd);normalized=True
    return state,params,normalized
def make_packet(state,identity):
    c=h.config()["estimator"];g=h.rng("observation",identity)
    counts=g.binomial(state["n"],c["observation_thinning"])
    d=np.clip(state["d"]+g.normal(0,c["dbh_noise_sd"],64),1,300)
    return h.ObservationPacket(tuple(map(int,counts)),tuple(float(v) if n else None for n,v in zip(counts,d)),(c["observation_thinning"],)*64)
def context(policy,p,r):
    if policy=="B2_unknown":return Context(policy)
    if policy=="B1_coarsened":return Context(policy,0 if p in [0,1] else 1,0 if r==0 else 1)
    return Context(policy,p,r)
def chunk_truth(state,params,model,identity,suite):
    c=h.config()["boundary"];folder=h.run_dir()/"truth";folder.mkdir(parents=True,exist_ok=True)
    path=folder/(str(identity)+"_"+suite+".npz")
    if path.exists():
        b=np.load(path);return b["viable"],b["reserve"],b["lawful"],bool(b["present"]),list(b["bank_hashes"])
    viable=[];reserve=[];lawful_flags=[];hashes=[]
    for chunk in range(max(c["truth_K"])//c["truth_chunk"]):
        g=h.rng("truth",identity,chunk)
        if suite=="entry_lognormal_annual_sd_1.2":g=h.EntryRNG(g)
        bank,lawful=h.demography.simulate(state,params,h.scenario(suite),5,c["truth_chunk"],g)
        rr=h.response(bank,lawful,model)
        viable.append(rr["viable"]);reserve.append(rr["reserve_viable"]);lawful_flags.append(lawful)
        hashes.append(hashlib.sha256(bank.tobytes()+lawful.tobytes()).hexdigest())
        present=rr["present"]
    v=np.concatenate(viable);r=np.concatenate(reserve);l=np.concatenate(lawful_flags)
    np.savez_compressed(path,viable=v,reserve=r,lawful=l,present=present,bank_hashes=np.array(hashes))
    return v,r,l,present,hashes
def run_truth():
    h.confirm_guard()
    _,_,_,core,model=h.load_data();c=h.config()["boundary"]
    cells=json.loads((h.ROOT/c["design"]).read_text())["worlds"]
    rows=[];manifest=[]
    for identity,cell in enumerate(cells):
        state,params,normalized=world(core,model,cell,identity)
        for suite in c["misspecifications"]:
            v,r,l,present,hashes=chunk_truth(state,params,model,identity,suite)
            if np.any(r&~v) or np.any(v&~l):raise ValueError("Capacity support violation")
            manifest.append(dict(world=identity,suite=suite,raw_response_sha256=hashlib.sha256(v.tobytes()+r.tobytes()+l.tobytes()).hexdigest(),bank_hashes=hashes))
            for K in c["truth_K"]:
                nv=int(v[:K].sum());nr=int(r[:K].sum())
                for query in ["Q1","Q2"]:
                    counts=[nv] if query=="Q1" else [nv,nr]
                    d=h.classify(present,counts,K,c["theta"])
                    rows.append(dict(world=identity,family=cell["family"],target=cell.get("target"),pressure=cell["pressure"],recovery=cell["recovery"],suite=suite,query=query,
                      K=K,nv=nv,nr=nr,present=present,FINITE_BANK_HEALTH=d["FINITE_BANK_HEALTH"],KERNEL_MC_STATUS=d["KERNEL_MC_STATUS"],
                      margin=min(counts)/K-c["theta"],margin_lower=min(d["lower"])-c["theta"],margin_upper=min(d["upper"])-c["theta"],
                      lower_V=d["lower"][0],upper_V=d["upper"][0],lower_R=d["lower"][-1] if query=="Q2" else None,upper_R=d["upper"][-1] if query=="Q2" else None,
                      failure_class="REALIZATION" if not present else ("Q1_CAPACITY" if nv/K<c["theta"] else ("Q2_RESERVE" if query=="Q2" and nr/K<c["theta"] else "HEALTHY")),
                      normalized_from_private_parameters=normalized))
        if identity%3==0:print("Boundary truth world",identity,"of",len(cells),flush=True)
    h.write_csv("boundary_validation/truth_prefixes.csv",rows)
    h.write_json("provenance/truth_bank_manifest.json",manifest)
    h.write_json("boundary_validation/pre_estimator_census.json",dict(worlds=len(cells),all_worlds_retained=True,truth_census_sha256=h.sha(h.ROOT/"boundary_validation/truth_prefixes.csv"),
        estimator_started=False,max_truth_paths=max(c["truth_K"]),pilot_streams_reused=False))
def run_estimates():
    h.confirm_guard()
    if not (h.ROOT/"boundary_validation/pre_estimator_census.json").exists():raise ValueError("Truth census must precede estimation")
    _,_,_,core,model=h.load_data();c=h.config()["estimator"]
    worlds=json.loads((h.ROOT/"config/world_design.json").read_text())["worlds"]
    anchor=dict(n=core["stats"]["n0"].copy(),d=model["diameter0"].copy())
    draws=[];packets=[]
    for identity,cell in enumerate(worlds):
        state,_,_=world(core,model,cell,identity)
        packets.append(("boundary",identity,identity,make_packet(state,identity),cell["pressure"],cell["recovery"]))
    # Historical worlds and measurements exactly reproduced from frozen seeds.
    import h03_synthetic as previous
    old=h.prior.PROTOCOL["synthetic"]
    grid=list(itertools.product(range(3),range(3),range(3),range(3),range(old["replicates_per_cell"])))
    for identity,(s,r,p,recovery,replicate) in enumerate(grid):
        state=previous.make_state(core,model,s,r,identity);g=h.prior.rng("observation",identity)
        counts=g.binomial(state["n"],old["observation_thinning"])
        d=np.clip(state["d"]+g.normal(0,old["dbh_noise_sd_cm"],64),1,300)
        packet=h.ObservationPacket(tuple(map(int,counts)),tuple(float(v) if n else None for n,v in zip(counts,d)),(old["observation_thinning"],)*64)
        packets.append(("HFD03",identity,10000+identity,packet,p,recovery))
    # All estimators receive only strict packets + policy-specific marker.
    for cohort,identity,seed,packet,p,r in packets:
        for policy in c["contexts"]:
            observation=packet if policy!="B3_label_only" else h.ObservationPacket((0,)*64,(None,)*64,(0.,)*64)
            result=evaluate(observation,context(policy,p,r),anchor,model,seed,c["draws"],c["paths"])
            for row in result:draws.append(dict(cohort=cohort,world=identity,policy=policy,**row))
        if identity%12==0:print("Context estimator",cohort,identity,flush=True)
    h.write_csv("context_ablation/estimator_draws.csv",draws)
    h.write_json("provenance/estimator_information.json",dict(packet_fields=["counts","diameters_cm","observation_probabilities"],
       exact_levels_in_unknown=False,exact_levels_in_coarse=False,true_params_transmitted=False,
       identity_used_only_as_random_nonce=True,truth_module_not_imported=True,label_only_has_no_state=True,
       comparison_precision="HFD0432 draws; historical HFD0364. Paired current B0 controls resolution change."))
if __name__=="__main__":
    import sys
    run_truth() if sys.argv[1]=="truth" else run_estimates()
