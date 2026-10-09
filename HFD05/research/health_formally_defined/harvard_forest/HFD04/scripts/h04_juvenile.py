"""Observed E1 size boundary plus explicitly imputed size-aware dynamics."""
import math
import numpy as np
import h04_common as h
import h03_calibration as cal
def visible_inputs(family):
    registration=h.prior.PROTOCOL
    held={r["stem_id"] for r in h.read_csv(h.THREE/"reconstruction_validation/mask_assignments.csv")
          if r["mask"]==family and r["replicate"]=="0"}
    path=h.prior.RUN/("visible_"+family+"_0.csv")
    return held,path
def profiles(e0,path):
    groups=[[] for _ in range(64)]
    for row in h.read_csv(path):
        old=e0.get(row["stem.id"])
        if cal.secure(old,row) and h.frozen.fate(row)==1 and h.empirical.number(row["dbh"]) is not None:
            groups[h.empirical.cell(old)].append(float(row["dbh"]))
    rows=[];values=[]
    for c,data in enumerate(groups):
        d=np.array(data);young=d[(d>=1)&(d<10)];adult=d[d>=10]
        fraction=len(young)/max(len(d),1)
        # Empty subgroup support is imputed, never silently called observed.
        a=float(young.mean()) if len(young) else None
        b=float(adult.mean()) if len(adult) else None
        values.append((fraction,a,b))
        rows.append(dict(cell=c,n_visible=len(d),young_visible=len(young),fraction=fraction,
                         below10_mean=a,atleast10_mean=b,profile_provenance="OBSERVED_VISIBLE_E1" if len(d) else "IMPUTED_FALLBACK"))
    return values,rows
def enrich(state,profile,generator):
    n=np.zeros(128,int);d=np.ones(128)
    for c,(fraction,a,b) in enumerate(profile):
        dc=state["d"][c];nc=state["n"][c]
        if a is None and b is None: fraction,a,b=.5,.75*dc,1.25*dc
        elif a is None:a=b
        elif b is None:b=a
        young=int(generator.binomial(nc,fraction))
        n[2*c:2*c+2]=[young,nc-young];d[2*c:2*c+2]=[a,b]
        second=float((n[2*c:2*c+2]*d[2*c:2*c+2]**2).sum())
        if nc and second:d[2*c:2*c+2]*=math.sqrt(nc*dc**2/second)
    # Clipping could break BA preservation: no silent claim of commutation.
    before=(n*d*d).sum()
    d=np.clip(d,1,300)
    return dict(n=n,d=d),float((n*d*d).sum()-before)
def simulate_enriched(state,params,scenario,years,count,generator,variant):
    initial=len(state["n"]);width=initial+16*years
    n=np.zeros((count,width),int);d=np.ones((count,width));n[:,:initial]=state["n"];d[:,:initial]=state["d"]
    cells=np.concatenate([np.repeat(np.arange(64),2),np.tile(np.arange(16)*4,years)])
    bank=np.zeros((count,years+1,6));deaths=np.zeros(count);entries=np.zeros(count)
    bank[:,0]=h.demography.metrics(n[:,:initial],d[:,:initial],cells[:initial],deaths,entries)
    lawful=np.ones(count,bool)
    for year in range(1,years+1):
        length=initial+16*(year-1);idx=cells[:length]
        hazard=params["hazard"][idx]*scenario["hazard_factor"]+(idx//16==0)*scenario["hemlock_extra_hazard"]
        survivors=generator.binomial(n[:,:length],np.exp(-hazard))
        deaths+=(n[:,:length]-survivors).sum(axis=1);n[:,:length]=survivors
        mean=params["growth"][idx]*scenario["growth_factor"];sd=params["sd"][idx]*scenario["growth_factor"]
        if variant=="gamma_growth":
            shape=np.maximum((mean/np.maximum(sd,1e-9))**2,1e-6)
            growth=generator.gamma(shape,np.where(mean>0,sd**2/np.maximum(mean,1e-9),0),size=(count,length))
        elif variant=="normal_growth":growth=np.maximum(0,generator.normal(mean,sd,size=(count,length)))
        else:raise ValueError("Unregistered growth kernel")
        d[:,:length]=np.clip(d[:,:length]+growth,1,300)
        born=generator.poisson(params["recruit"].reshape(16,4).sum(axis=1)*scenario["recruitment_factor"],size=(count,16))
        n[:,length:length+16]=born;d[:,length:length+16]=1.7;entries+=born.sum(axis=1)
        ids=[("initial",c,point) for c in range(64) for point in range(2)]+[("birth",t,g) for t in range(1,year+1) for g in range(16)]
        lawful &= h.demography.check_cohorts(n[:,:length+16],d[:,:length+16],ids)
        bank[:,year]=h.demography.metrics(n[:,:length+16],d[:,:length+16],cells[:length+16],deaths,entries)
    return bank,lawful & h.demography.check_history(bank)
def run():
    h.confirm_guard()
    c=h.config()["juvenile"];sources,e0,e1,core,model=h.load_data()
    pool=cal.covariates(e0,e1)
    live=[x for x in pool if h.frozen.fate(e1[x["stem_id"]])==1 and h.empirical.number(e1[x["stem_id"]]["dbh"]) is not None]
    diam=np.array([float(e1[x["stem_id"]]["dbh"]) for x in live]);cells=np.array([x["cell"] for x in live])
    true=(diam>=1)&(diam<10);groups=[]
    for cell in range(64):
        take=cells==cell;values=diam[take];rms=float(np.sqrt(np.mean(values**2))) if take.any() else 0
        classify=np.full(int(take.sum()),rms<10);truth=true[take]
        groups.append(dict(cell=cell,n=int(take.sum()),rms_cm=rms,stem_juveniles=int(truth.sum()),
                           rms_juveniles=int(classify.sum()),FP=int((classify&~truth).sum()),FN=int((~classify&truth).sum()),
                           first_date=min([e1[x["stem_id"]]["exact.date"] for x in live if x["cell"]==cell],default=""),
                           last_date=max([e1[x["stem_id"]]["exact.date"] for x in live if x["cell"]==cell],default=""),
                           support="SECURE_OBSERVED_E1_SUBSET_NOT_2020_ORIGIN"))
    h.write_csv("juvenile_boundary/observed_groups.csv",groups)
    edges=np.array(c["distance_edges_cm"]);dist=np.abs(diam-10);rows=[]
    for si,sigma in enumerate(c["sigma_D_cm"]):
        for rep in range(c["perturbation_replicates"]):
            perturbed=np.clip(diam+h.rng("measurement_error",rep).normal(0,1,len(diam))*sigma,1,300)
            observed=(perturbed<10)
            ba=float((perturbed**2).sum()*math.pi/40000)
            truthba=float((diam**2).sum()*math.pi/40000)
            for lo,hi in zip(edges[:-1],edges[1:]):
                take=(dist>=lo)&(dist<hi)
                rows.append(dict(sigma=sigma,replicate=rep,distance_low=float(lo),distance_high=float(hi),n=int(take.sum()),
                  FP=int((observed&~true&take).sum()),FN=int((~observed&true&take).sum()),deltaJ=int((observed[take].sum()-true[take].sum())),
                  delta_J_over_J0=float((observed[take].sum()-true[take].sum())/model["baseline_juveniles"]),
                  full_subset_realization_before=bool(truthba/model["baseline_ba"]>=.75 and true.sum()/model["baseline_juveniles"]>=.4),
                  full_subset_realization_after=bool(ba/model["baseline_ba"]>=.75 and observed.sum()/model["baseline_juveniles"]>=.4),
                  assumed_error_distribution="REGISTERED_METHODOLOGICAL_GAUSSIAN_NOT_EMPIRICALLY_ESTIMATED"))
    h.write_csv("juvenile_boundary/measurement_error.csv",rows)
    masks=[]
    for family in c["mask_families"]:
        held,path=visible_inputs(family)
        visible_core,visible_model=cal.fit_visible(sources,path,held)
        cov=[x for x in pool if x["stem_id"] in held]
        p,size,alive,_=cal.predict(visible_core,visible_model,cov,int(family[1]),0)
        truth=np.array([h.frozen.fate(e1[x["stem_id"]])==1 and h.empirical.number(e1[x["stem_id"]]["dbh"]) is not None and 1<=float(e1[x["stem_id"]]["dbh"])<10 for x in cov])
        eligible=np.array([h.frozen.fate(e1[x["stem_id"]])==0 or h.empirical.number(e1[x["stem_id"]]["dbh"]) is not None for x in cov])
        pred=alive*(size<10)
        estimate=pred.mean(axis=0)
        masks.append(dict(mask=family,n=int(eligible.sum()),truth_juveniles=int(truth[eligible].sum()),
                          mean_predicted_juveniles=float(pred[:,eligible].sum(axis=1).mean()),
                          delta_J_over_J0=float((pred[:,eligible].sum(axis=1).mean()-truth[eligible].sum())/model["baseline_juveniles"]),
                          membership_brier=float(np.mean((estimate[eligible]-truth[eligible])**2)),
                          masked_DBH_in_profile=False,visible_file_sha256=h.sha(path)))
        print("Juvenile mask",family,"complete",flush=True)
    h.write_csv("juvenile_boundary/masking.csv",masks)
    held,path=visible_inputs("M0");visible_core,visible_model=cal.fit_visible(sources,path,held)
    profile,profile_rows=profiles(e0,path);h.write_csv("juvenile_boundary/visible_size_profiles.csv",profile_rows)
    cfg=h.frozen.load("config/design.json");scenarios=[s for s in cfg["scenarios"] if s["id"] in c["scenarios"]]
    rows=[];initial=[]
    for i in range(c["states"]):
        state=h.demography.reconstruct(visible_core,visible_model,h.rng("representation_state",i,0))
        fine,clip=enrich(state,profile,h.rng("representation_state",i,1))
        coarseba=float((state["n"]*state["d"]**2).sum()*math.pi/40000);fineba=float((fine["n"]*fine["d"]**2).sum()*math.pi/40000)
        initial.append(dict(state=i,coarse_count=int(state["n"].sum()),fine_count=int(fine["n"].sum()),coarse_ba=coarseba,fine_ba=fineba,
                            BA_projection_error=fineba-coarseba,clip_second_moment_error=clip,
                            coarse_J=int((state["n"]*(state["d"]<10)).sum()),fine_J=int((fine["n"]*(fine["d"]<10)).sum()),
                            provenance="IMPUTED_ORIGIN_DISTRIBUTION_WITH_VISIBLE_E1_PROFILE_ANCESTRY"))
        params=h.demography.parameters(visible_model,h.rng("representation_state",i,2))
        for si,scenario in enumerate(scenarios):
            for vi,variant in enumerate(c["variants"]):
                banks=[h.demography.simulate(state,params,scenario,max(c["horizons"]),c["paths"],h.rng("representation_future",i,si,vi),variant),
                       simulate_enriched(fine,params,scenario,max(c["horizons"]),c["paths"],h.rng("representation_future",i,si,vi),variant)]
                for horizon in c["horizons"]:
                    for beta in c["semantic_grid"]["beta"]:
                        for rho in c["semantic_grid"]["rho"]:
                            rr=[h.demography.response(bank,lawful,visible_model,horizon,.75,beta,2,rho) for bank,lawful in banks]
                            for theta in c["semantic_grid"]["theta"]:
                                a,b=rr;nv=[int(x["viable"].sum()) for x in rr];nr=[int(x["reserve_viable"].sum()) for x in rr]
                                decision=[h.classify(x["present"],[v,r],c["paths"],theta) for x,v,r in zip(rr,nv,nr)]
                                flags=[a["present"]!=b["present"],nv[0]/c["paths"]>=theta and nv[1]/c["paths"]<theta or nv[1]/c["paths"]>=theta and nv[0]/c["paths"]<theta,
                                       nr[0]/c["paths"]>=theta and nr[1]/c["paths"]<theta or nr[1]/c["paths"]>=theta and nr[0]/c["paths"]<theta]
                                changed=decision[0]["FINITE_BANK_HEALTH"]!=decision[1]["FINITE_BANK_HEALTH"]
                                decomposition="NONE" if not changed else (["INITIAL_REALIZATION","CONTINUATION","RESERVE"][flags.index(True)] if sum(flags)==1 else "COMBINED")
                                rows.append(dict(state=i,scenario=scenario["id"],variant=variant,horizon=horizon,beta=beta,rho=rho,theta=theta,K=c["paths"],
                                  coarse_present=a["present"],fine_present=b["present"],coarse_nv=nv[0],fine_nv=nv[1],coarse_nr=nr[0],fine_nr=nr[1],
                                  coarse_Health=decision[0]["FINITE_BANK_HEALTH"],fine_Health=decision[1]["FINITE_BANK_HEALTH"],coarse_status=decision[0]["KERNEL_MC_STATUS"],fine_status=decision[1]["KERNEL_MC_STATUS"],change_class=decomposition))
        if i%4==0:print("Representation state",i,flush=True)
    h.write_csv("juvenile_boundary/initial_projection.csv",initial);h.write_csv("juvenile_boundary/propagation.csv",rows)
    h.write_json("juvenile_boundary/empirical_support.json",dict(n=len(live),eligible_pool=h.sha(h.THREE/"reconstruction_validation/eligible_truth_pool.csv"),
       observed_total_J=int(true.sum()),J0=model["baseline_juveniles"],baseline_BA=model["baseline_ba"],never_observed_origin=True,
       conditional_profile="Estimated visible distributions, rescaled to latent second moments; no exact transition commutation claimed"))
if __name__=="__main__":run()
