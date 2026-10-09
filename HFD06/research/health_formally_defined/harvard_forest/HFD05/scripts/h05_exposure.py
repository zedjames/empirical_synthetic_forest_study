"""Per-record exposure inventory and separately addressable reconstruction law."""
import datetime as dt,math,json
from collections import defaultdict
import numpy as np
import h05_common as h

def records(e0,e1,core):
    emp=h.original.empirical;fate=h.original.frozen.fate
    origin=dt.date.fromisoformat(h.original.frozen.load("config/design.json")["origin"]);rows=[]
    def add(identity,group,c,row,epoch,clamp):
        raw=row.get("exact.date")
        try:date=dt.date.fromisoformat(raw)
        except (ValueError,TypeError):raise ValueError("Unidentified exposure: "+str(identity)+" "+str(raw))
        exposure=(origin-date).days/365.25
        if clamp:exposure=max(0,exposure)
        if exposure<0:raise ValueError("Future-dated unresolved observation")
        rows.append(dict(stem_id=identity,component=group,cell=c,exposure_years=exposure,observed_date=raw,date_source=epoch+".exact.date",date_fallback="NONE",raw_status=row["df.status"],origin=origin.isoformat()))
    for identity,old in e0.items():
        if old["df.status"]!="alive":continue
        new=e1.get(identity);same=new is not None and new["tree.id"]==old["tree.id"]
        if (new is None or same) and (fate(new) if same else None) is None:add(identity,"unresolved_fate",emp.cell(old),old,"E0",False)
    for identity,new in e1.items():
        old=e0.get(identity)
        if old is not None and new["tree.id"]!=old["tree.id"]:continue
        if fate(new)==1:
            group="measured_alive" if emp.number(new["dbh"]) is not None else "alive_missing_size"
            add(identity,group,emp.cell(new,old),new,"E1",True)
    rows.sort(key=lambda r:(r["component"],r["cell"],r["stem_id"]))
    groups={name:[[] for _ in range(64)] for name in ["measured_alive","unresolved_fate","alive_missing_size"]}
    for r in rows:groups[r["component"]][r["cell"]].append(r["exposure_years"])
    groups={k:[np.array(v,float) for v in cells] for k,cells in groups.items()}
    for name,count,time in [("measured_alive","n1","tail_years"),("unresolved_fate","unknown","unknown_dt"),("alive_missing_size","alive_missing_size","alive_missing_size_tail")]:
        if not np.array_equal([len(v) for v in groups[name]],core["stats"][count]) or not np.allclose([v.sum() for v in groups[name]],core["stats"][time],atol=1e-8,rtol=1e-12):raise ValueError("Per-record exposure differs from original pooled sufficient statistics")
    return rows,groups

def poisson_binomial(p):
    pmf=np.array([1.])
    for value in p:pmf=np.convolve(pmf,[1-value,value])
    return pmf

def individual(core,model,groups,generator,family="B0"):
    cfg=h.original.frozen.load("config/design.json");rc=cfg["reconstruction"]
    if family not in h.config("exposure_protocol.json")["state_comparison"]["missingness_families"]:raise ValueError("Unregistered missingness family")
    p=h.original.demography.parameters(model,generator);wet=(np.arange(64)//4)%4==3
    mnar=generator.normal(0,rc["mnar_log_hazard_shift_sd"]);wshift=generator.normal(0,rc["wet_log_hazard_shift_sd"])
    fc=h.original.prior.PROTOCOL["missingness"][family]
    if family!="B0":
        mnar=mnar/.45*fc["mnar_sd"];wshift=wshift/.7*fc["wet_sd"]+fc.get("wet_log_shift_mean",0)
    unknown_h=p["hazard"]*np.exp(mnar+wet*wshift);component_h=unknown_h.copy()
    if family=="B1":
        c=np.arange(64);shift=np.array(fc["taxon_log_shift"])[c//16]+np.array(fc["size_log_shift"])[c%4]+np.array(fc["stratum_log_shift"])[(c//4)%4]
        component_h*=np.exp(shift)
    counts={};squares={}
    def component(name,hazard,base,noise_sd):
        noise=generator.normal(0,noise_sd,64);N=np.zeros(64,int);S=np.zeros(64)
        for cell,times in enumerate(groups[name]):
            alive=generator.binomial(1,np.exp(-hazard[cell]*times)).astype(bool)
            sizes=np.clip(base[cell]+p["growth"][cell]*times+noise[cell],1,300)
            N[cell]=int(alive.sum());S[cell]=float((sizes[alive]**2).sum())
        counts[name]=N;squares[name]=S
    component("measured_alive",p["hazard"],model["diameter1"],cfg["fit"]["measurement_dbh_sd_cm"])
    component("unresolved_fate",component_h,model["diameter0"],rc["synthetic_dbh_prior_sd_cm"])
    unseen=generator.poisson(p["recruit"]*cfg["fit"]["recruitment_exposure_years"]*wet*np.exp(generator.normal(0,rc["unseen_recruit_factor_log_sd"])))
    counts["unseen_entry"]=unseen;squares["unseen_entry"]=unseen*1.7**2
    component("alive_missing_size",p["hazard"],model["diameter0"],rc["synthetic_dbh_prior_sd_cm"])
    n=sum(counts.values());d2=sum(squares.values());branches=[0,0];branch_ba=0.
    for conflict in core["conflicts"]:
        branch=int(generator.integers(2));branches[branch]+=1;c=conflict["old_cell"] if branch==0 else conflict["new_cell"]
        if branch==1 and conflict["new_fate"] is not None:alive=bool(conflict["new_fate"]) and bool(generator.random()<math.exp(-p["hazard"][c]*conflict["tail"]))
        else:alive=bool(generator.random()<math.exp(-unknown_h[c]*conflict["years"]))
        diameter=conflict["new_diameter"] if branch==1 else None
        diameter=diameter if diameter is not None else conflict["old_diameter"]+p["growth"][c]*conflict["years"]
        if alive:
            n[c]+=1;square=np.clip(diameter,1,300)**2;d2[c]+=square;branch_ba+=square*math.pi/40000
    d=np.sqrt(np.divide(d2,n,out=model["diameter0"]**2,where=n>0))
    components=dict(measured_projected_ba=float(squares["measured_alive"].sum()*math.pi/40000),
        imputed_ba=float(sum(v.sum() for k,v in squares.items() if k!="measured_alive")*math.pi/40000+branch_ba),
        wet_imputed_stems=int(sum(v[wet].sum() for k,v in counts.items() if k!="measured_alive")),association_branches=branches,mnar_shift=float(mnar),wet_shift=float(wshift))
    return dict(n=n,d=d,components=components,ancestry=core["identity_ledger_sha256"])

def inventory():
    _,e0,e1,core,model=h.original.load_data();rows,groups=records(e0,e1,core);c=h.config("exposure_protocol.json")["inventory"]
    h.write_csv("exposure_model/individual_exposure_records.csv",rows);table=[];dp=[];moments=[]
    for gi,(name,cells) in enumerate(groups.items()):
        first=next((i for i,v in enumerate(cells) if len(v)),None)
        for cell,times in enumerate(cells):
            if not len(times):continue
            post=model["hazard_posterior"][cell];grid=model["hazard_grid"];cdf=post.cumsum()
            quantiles=np.quantile(times,c["exposure_quantiles"])
            for qi,q in enumerate(c["hazard_quantiles"]):
                base=float(grid[min(int(np.searchsorted(cdf,q)),255)])
                for shift in c["shared_log_shift_values"] if name=="unresolved_fate" else [0]:
                    hazard=base*math.exp(shift);p=np.exp(-hazard*times);pooled=math.exp(-hazard*times.mean());mean=float(p.sum());var=float((p*(1-p)).sum());gap=mean-len(times)*pooled
                    if gap < -1e-9:raise ValueError("Conditional Jensen inequality violated")
                    table.append(dict(component=name,cell=cell,n=len(times),exposure_min=float(times.min()),exposure_max=float(times.max()),exposure_mean=float(times.mean()),exposure_variance=float(times.var()),
                        exposure_p05=float(quantiles[1]),exposure_median=float(quantiles[2]),exposure_p95=float(quantiles[3]),source="E0.exact.date" if name=="unresolved_fate" else "E1.exact.date",fallback="NONE",
                        hazard_quantile=q,base_hazard=base,shared_log_shift=shift,conditional_hazard=hazard,PB_mean=mean,PB_variance=var,pooled_mean=len(times)*pooled,pooled_variance=len(times)*pooled*(1-pooled),
                        Jensen_gap=gap,relative_Jensen_gap=gap/max(mean,1e-300),variance_gap=var-len(times)*pooled*(1-pooled),conditional_independence=True))
            if cell==first:
                t=times[:20];p=np.exp(-.12*t);pmf=poisson_binomial(p);n=len(t);pooled=math.exp(-.12*t.mean())
                for count,prob in enumerate(pmf):dp.append(dict(component=name,cell=cell,n=n,count=count,PB_probability=float(prob),pooled_binomial_probability=math.comb(n,count)*pooled**count*(1-pooled)**(n-count),hazard=.12,selection="First source-order records; not selected by outcome"))
                if abs(pmf.sum()-1)>1e-12 or abs(np.dot(np.arange(n+1),pmf)-p.sum())>1e-12:raise ValueError("Independent PB DP identity failed")
                t=times[:c["large_moment"]["max_individuals"]];p=np.exp(-.12*t);N=c["large_moment"]["draws"];K=c["large_moment"]["block"];g=h.rng("exposure_moment",gi)
                draws=np.concatenate([g.binomial(1,p,size=(K,len(p))).sum(axis=1) for _ in range(N//K)]);mean=float(p.sum());var=float((p*(1-p)).sum());se=math.sqrt(var/N)
                fourth=3*var**2+float(np.sum(p*(1-p)*(1-6*p*(1-p))));vse=math.sqrt(max(fourth-var**2,0)/N)
                if abs(draws.mean()-mean)>8*se+1e-12 or abs(np.mean((draws-mean)**2)-var)>8*vse+1e-12:raise ValueError("Individual Bernoulli moment failure")
                moments.append(dict(component=name,cell=cell,n=len(p),draws=N,analytic_mean=mean,empirical_mean=float(draws.mean()),analytic_variance=var,empirical_centered_second_moment=float(np.mean((draws-mean)**2)),mean_SE=se,variance_SE=vse,status="PASS"))
    h.write_csv("exposure_model/cohort_exposure_inventory.csv",table);h.write_csv("exposure_model/poisson_binomial_dp.csv",dp);h.write_csv("exposure_model/bernoulli_moment_tests.csv",moments)
    durations=np.array([t for values in core["durations"] for t,z in values]);survival=core["stats"]["alive"].sum()/(core["stats"]["alive"]+core["stats"]["dead"]).sum()
    h.write_json("exposure_model/pooled_prior_inventory.json",dict(n=len(durations),min=float(durations.min()),max=float(durations.max()),mean=float(durations.mean()),variance=float(durations.var()),
        pooled_hazard=-math.log(survival)/durations.mean(),operation="Approximate pooled prior center only; original individual-duration likelihood retained",corrected_fit_same=True))
    print("Individual exposures audited",len(rows),"records; no fallback dates",flush=True)

def state_comparison():
    _,e0,e1,core,model=h.original.load_data();_,groups=records(e0,e1,core);c=h.config("exposure_protocol.json")["state_comparison"];rows=[];strata=[]
    designs=[("HFD02S_PRIMARY","B0",c["primary_draws"])] + [("HFD03_MISSINGNESS",family,c["missingness_draws"]) for family in c["missingness_families"]]
    for design,family,total in designs:
        for i in range(total):
            generator=lambda:h.original.frozen.Seeds().rng("reconstruct",i) if design=="HFD02S_PRIMARY" else h.original.prior.rng("state",i)
            old=h.original.demography.reconstruct(core,model,generator()) if design=="HFD02S_PRIMARY" else h.original.prior.reconstruct(core,model,i,family)
            new=individual(core,model,groups,generator(),family)
            def values(state,mask=np.ones(64,bool)):
                n,d=state["n"][mask],state["d"][mask]
                return dict(N=int(n.sum()),BA=float((n*d*d).sum()*math.pi/40000),J=int((n*(d<10)).sum()))
            a,b=values(old),values(new)
            row=dict(design=design,family=family,state=i,**{"original_"+k:v for k,v in a.items()},**{"individual_"+k:v for k,v in b.items()},**{k+"_delta":b[k]-a[k] for k in a},
                original_measured_projected_ba=old["components"]["measured_projected_ba"],individual_measured_projected_ba=new["components"]["measured_projected_ba"],
                original_imputed_ba=old["components"]["imputed_ba"],individual_imputed_ba=new["components"]["imputed_ba"],
                original_mnar=old["components"]["mnar_shift"],individual_mnar=new["components"]["mnar_shift"],original_wet=old["components"]["wet_shift"],individual_wet=new["components"]["wet_shift"])
            if row["original_mnar"]!=row["individual_mnar"] or row["original_wet"]!=row["individual_wet"]:raise ValueError("Shared missingness shift changed")
            rows.append(row)
            for axis,assignments in [("taxon",np.arange(64)//16),("E0_sector_proxy",(np.arange(64)//4)%4)]:
                for label in range(4):
                    a,b=values(old,assignments==label),values(new,assignments==label)
                    strata.append(dict(design=design,family=family,state=i,axis=axis,label=label,**{"original_"+k:v for k,v in a.items()},**{"individual_"+k:v for k,v in b.items()},**{k+"_delta":b[k]-a[k] for k in a}))
            if i%32==0:print("Exposure states",design,family,i,flush=True)
    h.write_csv("exposure_model/origin_state_comparison.csv",rows);h.write_csv("exposure_model/origin_stratum_comparison.csv",strata)
    summaries=[]
    for design,family,_ in designs:
        selected=[r for r in rows if r["design"]==design and r["family"]==family]
        for metric in ["N","BA","J"]:
            delta=np.array([r[metric+"_delta"] for r in selected])
            summaries.append(dict(design=design,family=family,metric=metric,n=len(delta),mean_delta=float(delta.mean()),q95_abs_delta=float(np.quantile(abs(delta),.95)),max_abs_delta=float(abs(delta).max()),scope="Coupled-index reconstruction distributions; no exact common-uniform individual-path claim"))
    h.write_csv("exposure_model/origin_comparison_summary.csv",summaries)

if __name__=="__main__":
    import sys
    h.guard("exposure");{"inventory":inventory,"states":state_comparison}[sys.argv[1]]()
