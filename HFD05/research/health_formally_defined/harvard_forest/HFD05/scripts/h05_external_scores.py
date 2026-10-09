"""Frozen adult component predictions and distinct descriptive seedling processes."""
import datetime as dt,itertools,json,math
from collections import defaultdict,Counter
import numpy as np
import h05_common as h
from h05_external_io_v2 import load
from h05_external_align import tag

def bootstrap(rows,cluster,identity):
    keys=sorted({r[cluster] for r in rows});by={k:i for i,k in enumerate(keys)};totals=np.zeros((len(keys),5))
    for r in rows:
        p,y=r["predicted_survival"],r["observed_survival"];totals[by[r[cluster]]]+=[1,1-y,1-p,(p-y)**2,-(y*math.log(max(p,1e-15))+(1-y)*math.log(max(1-p,1e-15)))]
    if not rows:return dict(cluster_units=0,mortality_gap_lower=None,mortality_gap_upper=None,brier_lower=None,brier_upper=None)
    g=h.rng("external_bootstrap",*identity);weights=g.multinomial(len(keys),np.full(len(keys),1/len(keys)),size=1000);sample=weights@totals;rates=(sample[:,1]-sample[:,2])/sample[:,0];brier=sample[:,3]/sample[:,0]
    lo,hi=np.quantile(rates,[.025,.975]);bl,bh=np.quantile(brier,[.025,.975])
    return dict(cluster_units=len(keys),mortality_gap_lower=float(lo),mortality_gap_upper=float(hi),brier_lower=float(bl),brier_upper=float(bh))

def adult():
    h.guard("external_scoring");raw,source=load("HF453","hf453-01");_,_,_,core,model=h.original.load_data();cross=h.read_csv(h.ROOT/"external_sources/HF453/tag_crosswalk.csv")
    histories=defaultdict(dict)
    for r in raw:histories[tag(r["StemTag"])][int(r["year"])]=r
    eligible=[r for r in cross if r["primary_eligible"]=="True"];rows=[];denominators=[]
    windows=[("CUMULATIVE",2021,end) for end in [2022,2023,2024]]+[("ANNUAL",start,start+1) for start in [2021,2022,2023]]
    for wi,(kind,start,end) in enumerate(windows):
        for scheme in ["STRICT","DN_DEATH","DN_ALIVE"]:
            risk=usable=dn=missing=0
            for entry in eligible:
                identity=entry["external_tag"];history=histories[identity];baseline=history.get(start);endpoint=history.get(end)
                if baseline is None or baseline["status"] not in ["A","AU"]:continue
                risk+=1
                if endpoint is None or endpoint["status"] not in ["A","AU","DC","DS","DN"]:missing+=1;continue
                status=endpoint["status"]
                if status=="DN":
                    dn+=1
                    if scheme=="STRICT":continue
                y=1 if status in ["A","AU"] or status=="DN" and scheme=="DN_ALIVE" else 0;cell=int(entry["cell"]);duration=end-start
                grid=model["hazard_grid"];posterior=model["hazard_posterior"][cell];p=float(np.dot(posterior,np.exp(-grid*duration)))
                pmin=float(np.dot(posterior,np.exp(-grid*(duration+30/365.25))));pmax=float(np.dot(posterior,np.exp(-grid*(duration-30/365.25))))
                damage=histories[identity][2021];fields=[damage[k] for k in ["FAD","wound","canker","rot","DWR"]];numeric=[]
                for v in fields:
                    try:numeric.append(float(v))
                    except ValueError:pass
                damage_flag="RECORDED_DAMAGE" if any(v>0 for v in numeric) else ("NO_RECORDED_DAMAGE" if len(numeric)==len(fields) else "UNKNOWN")
                q=entry["E0_quadrat"];cluster="QUADRAT:"+q if q not in ["","NA"] else "UNRESOLVED_QUADRAT_SECTOR:"+entry["sector_proxy"]
                rows.append(dict(target=kind,start_year=start,end_year=end,scheme=scheme,external_tag=identity,linked_stem=entry["linked_stem"],cell=cell,taxon=int(entry["taxon"]),sector_proxy=int(entry["sector_proxy"]),
                    size_group="10_TO_30" if float(entry["E1_dbh"])<30 else "30_PLUS",hemlock="HEMLOCK" if int(entry["taxon"])==0 else "OTHER",damage=damage_flag,spatial_cluster=cluster,
                    raw_start_status=baseline["status"],raw_endpoint_status=status,exposure_years=duration,time_precision="July year only; anniversary approximation",observed_survival=y,predicted_survival=p,predicted_survival_July_lower=pmin,predicted_survival_July_upper=pmax,
                    residual=y-p,brier=(p-y)**2,negative_log_score=-(y*math.log(max(p,1e-15))+(1-y)*math.log(max(1-p,1e-15))),left_truncated2021=True,
                    prior_secure_death_before_endpoint=any(r["status"] in ["DC","DS"] for year,r in history.items() if start<year<end),model_refitted=False))
                usable+=1
            denominators.append(dict(target=kind,start_year=start,end_year=end,scheme=scheme,risk_stems=risk,scored_stems=usable,DN_endpoint_stems=dn,missing_or_blank_endpoint=missing,DN_primary_secure_death=False,
                partition_reconciled=risk==usable+missing+(dn if scheme=="STRICT" else 0)))
    h.write_csv("external_sources/HF453/mortality_predictions.csv",rows);h.write_csv("external_sources/HF453/risk_denominators.csv",denominators)
    summaries=[];calibration=[];intervals=[]
    axes=[("ALL",["ALL"]),("taxon",list(range(4))),("hemlock",["HEMLOCK","OTHER"]),("sector_proxy",list(range(4))),("size_group",["10_TO_30","30_PLUS"]),("damage",["RECORDED_DAMAGE","NO_RECORDED_DAMAGE","UNKNOWN"])]
    for wi,(kind,start,end) in enumerate(windows):
        for si,scheme in enumerate(["STRICT","DN_DEATH","DN_ALIVE"]):
            selected=[r for r in rows if (r["target"],r["start_year"],r["end_year"],r["scheme"])==(kind,start,end,scheme)]
            for ai,(axis,labels) in enumerate(axes):
                for li,label in enumerate(labels):
                    chosen=selected if axis=="ALL" else [r for r in selected if r[axis]==label];n=len(chosen);observed=sum(1-r["observed_survival"] for r in chosen);predicted=sum(1-r["predicted_survival"] for r in chosen)
                    value=dict(target=kind,start_year=start,end_year=end,scheme=scheme,axis=axis,label=label,n=n,observed_deaths=observed,predicted_deaths=predicted,
                        observed_mortality=observed/n if n else None,predicted_mortality=predicted/n if n else None,mortality_gap=(observed-predicted)/n if n else None,
                        brier=sum(r["brier"] for r in chosen)/n if n else None,negative_log_score=sum(r["negative_log_score"] for r in chosen)/n if n else None,
                        source_identity="HFD05_EXTERNAL_COMPONENT_VALIDATION",status="PARTIALLY_ALIGNED_AND_SCORED",target_scope="Conditional observed2021 adult survivor cohort; not full Health or publicly issued prospective forecast")
                    summaries.append(value)
                    for ui,unit in enumerate(["spatial_cluster","external_tag"]):
                        intervals.append(dict(target=kind,start_year=start,end_year=end,scheme=scheme,axis=axis,label=label,bootstrap_unit=unit,replicates=1000,
                            **bootstrap(chosen,unit,(0,wi,si,ai,li,ui)),scope="Conditional source-support cluster bootstrap; descriptive, not unbiased ecological population interval"))
            for lo,hi in zip([0,.2,.4,.6,.8],[.2,.4,.6,.8,1]):
                chosen=[r for r in selected if lo<=r["predicted_survival"] and (r["predicted_survival"]<hi or hi==1)]
                calibration.append(dict(target=kind,start_year=start,end_year=end,scheme=scheme,lower=lo,upper=hi,n=len(chosen),mean_predicted_survival=float(np.mean([r["predicted_survival"] for r in chosen])) if chosen else None,
                    observed_survival=float(np.mean([r["observed_survival"] for r in chosen])) if chosen else None))
        print("Adult independent component scored",kind,start,end,flush=True)
    h.write_csv("external_sources/HF453/mortality_scores.csv",summaries);h.write_csv("external_sources/HF453/hemlock_scores.csv",[r for r in summaries if r["axis"]=="hemlock"])
    h.write_csv("external_sources/HF453/calibration_bins.csv",calibration);h.write_csv("uncertainty/HF453_cluster_intervals.csv",intervals)
    h.write_json("external_sources/HF453/validation_summary.json",dict(status="PARTIALLY_ALIGNED_AND_SCORED",raw_source_sha256=source["sha256"],frozen_hazard_model_only=True,external_outcome_fit=False,
        record_predictions=len(rows),baseline_stems=len(eligible),all_registered_windows_and_schemes=True,DN_primary_death=False,left_truncated2021=True,full_Health_validation=False,publicly_issued_prospective=False))

def seedlings():
    h.guard("external_scoring");raw,source=load("HF355","hf355-02");cross=h.read_csv(h.ROOT/"external_sources/HF355/identity_crosswalk.csv");eligible={r["seedlingID"] for r in cross if r["identity_eligible"]=="True"};groups=defaultdict(list)
    for index,r in enumerate(raw):groups[r["seedlingID"]].append((index,r))
    transitions=[];cohorts=[];graduation=[];excluded=[]
    for identity,visits in sorted(groups.items()):
        if identity not in eligible:continue
        visits=sorted(visits,key=lambda item:int(item[1]["yearOfOb"]));first=visits[0][1];dead=False;reversed_history=False
        for index,r in visits:
            if r["sampled"]=="1" and r["alive"]=="N":dead=True
            elif dead and r["sampled"]=="1" and r["alive"]=="Y":reversed_history=True
        cohort={r["yrFirstObs"] for index,r in visits};recruit={r["trueRecrt"] for index,r in visits}
        cohorts.append(dict(seedlingID=identity,coordinate=first["coordinate"],taxonCode=first["taxonCode"],yrFirstObs=first["yrFirstObs"],trueRecrt=first["trueRecrt"],stable_cohort_fields=len(cohort)==len(recruit)==1,
            baseline2017=first["yrFirstObs"]=="2017",postbaseline_field_recruit=first["yrFirstObs"]!="2017" and first["trueRecrt"]=="1",whole_plot_biological_entry=False,germination_date_identified=False))
        vg=[r for index,r in visits if r["status"]=="VG"]
        if vg:
            graduation.append(dict(seedlingID=identity,coordinate=first["coordinate"],taxonCode=first["taxonCode"],VG_visit_rows=len(vg),first_recorded_VG_year=vg[0]["yearOfOb"],first_recorded_VG_date=vg[0]["date"],
                recorded_at_initial_visit=vg[0]["yearOfOb"]==first["yearOfOb"],postbaseline_graduation_interval_identified=False,whole_plot_entry_identified=False,interpretation="Repeated graduation attribute; first observed already graduated, transition date left-censored"))
        if reversed_history or vg:
            for index,r in visits:excluded.append(dict(seedlingID=identity,source_row=index,year=r["yearOfOb"],reason="SECURE_N_TO_Y_REVERSAL" if reversed_history else "GRADUATED_STATUS_DISTINCT_SIZE_SUPPORT"))
            continue
        usable=[]
        for index,r in visits:
            if r["sampled"]!="1" or r["alive"] not in ["Y","N"]:continue
            try:date=dt.date.fromisoformat(r["date"])
            except (ValueError,TypeError):continue
            if date.year!=int(r["yearOfOb"]):continue
            usable.append((index,r,date))
        for (ia,a,da),(ib,b,db) in zip(usable,usable[1:]):
            duration=(db-da).days/365.25
            if a["alive"]!="Y" or duration<=0:continue
            transitions.append(dict(seedlingID=identity,coordinate=a["coordinate"],taxonCode=a["taxonCode"],start_year=a["yearOfOb"],end_year=b["yearOfOb"],start_date=a["date"],end_date=b["date"],exposure_years=duration,
                observed_survival=1 if b["alive"]=="Y" else 0,flagged=a["flag"]=="1" or b["flag"]=="1",campaign_gap=int(b["yearOfOb"])-int(a["yearOfOb"]),
                model_prediction=None,support="Below1cm subplot process, no fitted matching hazard",whole_plot_expansion=False))
    h.write_csv("external_sources/HF355/transition_records.csv",transitions);h.write_csv("external_sources/HF355/first_observation_cohorts.csv",cohorts)
    h.write_csv("external_sources/HF355/graduation_records.csv",graduation,fields=["seedlingID","coordinate","taxonCode","VG_visit_rows","first_recorded_VG_year","first_recorded_VG_date","recorded_at_initial_visit","postbaseline_graduation_interval_identified","whole_plot_entry_identified","interpretation"])
    h.write_csv("external_sources/HF355/component_exclusions.csv",excluded,fields=["seedlingID","source_row","year","reason"])
    summaries=[];intervals=[]
    panels=[("ALL","ALL",transitions)]+[("taxon",code,[r for r in transitions if r["taxonCode"]==code]) for code in sorted({r["taxonCode"] for r in raw})]+[("end_year",year,[r for r in transitions if r["end_year"]==year]) for year in ["2018","2019","2020","2021"]]+[("UNFLAGGED_ONLY","0",[r for r in transitions if not r["flagged"]])]
    for pi,(axis,label,selected) in enumerate(panels):
        n=len(selected);alive=sum(r["observed_survival"] for r in selected);summaries.append(dict(axis=axis,label=label,n=n,survivors=alive,deaths=n-alive,survival_fraction=alive/n if n else None,
            total_exposure_years=sum(r["exposure_years"] for r in selected),distinct_seedlings=len({r["seedlingID"] for r in selected}),distinct_subplots=len({r["coordinate"] for r in selected}),status="DESCRIPTIVE_ONLY",whole_plot_bridge_identified=False))
        plot=sorted({r["coordinate"] for r in selected});totals=np.zeros((len(plot),2));by={k:i for i,k in enumerate(plot)}
        for r in selected:totals[by[r["coordinate"]]]+=[1,r["observed_survival"]]
        if plot:
            g=h.rng("external_bootstrap",1,pi);weights=g.multinomial(len(plot),np.full(len(plot),1/len(plot)),1000);samples=weights@totals;lo,hi=np.quantile(samples[:,1]/samples[:,0],[.025,.975]);lo,hi=float(lo),float(hi)
        else:lo=hi=None
        intervals.append(dict(axis=axis,label=label,subplot_units=len(plot),replicates=1000,survival_lower=lo,survival_upper=hi,scope="Whole-subplot repeat-history bootstrap, descriptive sampled support only"))
    entry=[]
    for year in ["2017","2018","2019","2020","2021"]:
        chosen=[r for r in cohorts if r["yrFirstObs"]==year and r["stable_cohort_fields"]];n=len(chosen);recruits=sum(r["trueRecrt"]=="1" for r in chosen)
        entry.append(dict(first_observation_year=year,identities=n,field_trueRecrt=recruits,field_recruit_fraction=recruits/n if n and year!="2017" else None,baseline_left_truncation=year=="2017",biological_germination_date=False,whole_plot1to10cm_entry_bridge=False))
    h.write_csv("external_sources/HF355/transition_summary.csv",summaries);h.write_csv("external_sources/HF355/entry_cohort_summary.csv",entry);h.write_csv("uncertainty/HF355_subplot_intervals.csv",intervals)
    h.write_json("external_sources/HF355/validation_summary.json",dict(status="DESCRIPTIVE_ONLY",source_sha256=source["sha256"],dated_survival_transitions=len(transitions),graduated_identities=len(graduation),graduation_time_identified=False,
        whole_plot_J_bridge="IDENTIFICATION_OBSTRUCTION",seedling_specific_hazard_model="ABSENT",source_reversal_exclusions=len({r["seedlingID"] for r in excluded if r["reason"]=="SECURE_N_TO_Y_REVERSAL"}),external_model_fit=False,first_observation_is_germination=False,full_Health_validation=False))
    print("Descriptive independent seedling transitions",len(transitions),flush=True)

if __name__=="__main__":
    import sys
    {"adult":adult,"seedlings":seedlings}[sys.argv[1]]()
