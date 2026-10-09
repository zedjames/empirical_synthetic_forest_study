"""Actual-call audit, analytic moments, fixed moment grid and matched impact."""
import itertools,json,math,hashlib
from collections import defaultdict
import numpy as np
import h05_common as h
from h05_growth import laws,sample,simulate
from h05_inference import classify

def unit_grid():
    c=h.config("gamma_protocol.json")["unit_grid"];rows=[];clip=[];calls=[]
    grid=list(itertools.product(c["means"],c["sds"]))
    # Instrument actual historical simulator calls, not just a copied formula.
    mu=np.zeros(64);sd=np.zeros(64)
    for i,(m,s) in enumerate(grid):mu[i]=m;sd[i]=s
    params=dict(growth=mu,sd=sd,hazard=np.zeros(64),recruit=np.zeros(64))
    scenario=dict(growth_factor=1,hazard_factor=1,recruitment_factor=1,hemlock_extra_hazard=0)
    def capture(m,s,a,b,size):
        for i,(mm,ss,aa,bb) in enumerate(zip(m,s,a,b)):
            calls.append(dict(cell=i,mu=float(mm),sd=float(ss),actual_shape=float(aa),actual_scale=float(bb),actual_mean=float(aa*bb),actual_variance=float(aa*bb*bb)))
    simulate(dict(n=np.ones(64,int),d=np.ones(64)),params,scenario,1,4,h.rng("gamma_moment",0),corrected=False,recorder=capture)
    for i,(mu,sd) in enumerate(grid):
        for corrected in [False,True]:
            law=laws(mu,sd,corrected);a,b,mean,var=[float(law[k]) for k in ["shape","scale","mean","variance"]]
            x=sample(mu,sd,c["draws"],h.rng("gamma_moment",1,i,int(corrected)),corrected)
            se=math.sqrt(var/c["draws"]);vse=math.sqrt((2*a*a+6*a)*b**4/c["draws"]) if var else 0
            observed=float(x.mean());mse=float(np.mean((x-mean)**2));tol=1e-13
            mean_ok=abs(observed-mean)<=8*se+tol;variance_ok=abs(mse-var)<=8*vse+tol
            if not mean_ok or not variance_ok:raise ValueError("Predeclared moment diagnostic failed")
            unresolved=mean>0 and se/mean>.1
            rows.append(dict(case=i,law="CORRECTED" if corrected else "ORIGINAL",mu=mu,sd=sd,shape=a,scale=b,analytic_mean=mean,analytic_variance=var,
                empirical_mean=observed,empirical_centered_second_moment=mse,analytic_mean_se=se,analytic_second_moment_se=vse,
                mean_bound_pass=mean_ok,variance_bound_pass=variance_ok,simulation_precision="MC_UNRESOLVED" if unresolved else "DIAGNOSTIC_RESOLVED",
                mean_error_vs_requested=mean-mu,variance_error_vs_requested=var-sd*sd,mu_zero_positive_variance_infeasible=mu==0 and sd>0))
            for diameter in c["clip_diameters"]:
                y=np.clip(diameter+x,1,300)-diameter
                clip.append(dict(case=i,law=rows[-1]["law"],diameter=diameter,raw_increment_mean=observed,clipped_increment_mean=float(y.mean()),
                    raw_increment_variance=float(x.var()),clipped_increment_variance=float(y.var()),upper_clip_fraction=float((diameter+x>300).mean()),
                    clipping_distinct_from_gamma_parameter_defect=True))
    h.write_csv("gamma_growth/actual_call_unit_inventory.csv",calls)
    h.write_csv("gamma_growth/moment_unit_tests.csv",rows);h.write_csv("gamma_growth/clipping_impact.csv",clip)

def parameter_sets(model,core):
    seed=h.original.frozen.Seeds()
    for i,p in itertools.product(range(256),range(2)):
        yield "HFD02S_PRIMARY",i*2+p,h.original.demography.parameters(model,seed.rng("parameters",i,p)),[1,.65,.85,.5],20
    import h04_boundary as boundary
    cells=json.loads((h.FOUR/"config/world_design.json").read_text())["worlds"]
    for i,cell in enumerate(cells):
        yield "HFD04_BOUNDARY_TRUTH",i,h.original.parameters(model,cell["pressure"],cell["recovery"],h.original.rng("world_parameters",i)),[1],5
    import h03_synthetic as historical
    grid=list(itertools.product(range(3),range(3),range(3),range(3),range(4)))
    for i,(_,_,p,r,_) in enumerate(grid):
        yield "HFD03_SYNTHETIC_TRUTH",i,historical.est.public_parameters(model,historical.est.PublicContext(p,r),h.original.prior.rng("world",i,1)),[1],5
    for i in range(h.original.prior.PROTOCOL["missingness"]["state_draws"]):
        yield "HFD03_MISSINGNESS_B0_B4_SHARED_FUTURE_PARAMETERS",i,h.original.demography.parameters(model,h.original.prior.rng("future",i,0)),[1,.65,.85,.5],20
    for i,(_,_,p,r,_) in enumerate(grid):
        for draw in range(h.original.prior.PROTOCOL["synthetic"]["estimator_draws"]):
            yield "HFD03_ORIGINAL_OBSERVATION_ESTIMATOR",i*64+draw,historical.est.public_parameters(model,historical.est.PublicContext(p,r),h.original.prior.rng("estimate_params",i,draw)),[1],5
    worlds=[("boundary",i,i,c["pressure"],c["recovery"]) for i,c in enumerate(cells)]
    worlds += [("HFD03",i,10000+i,p,r) for i,(_,_,p,r,_) in enumerate(grid)]
    for cohort,i,nonce,p,r in worlds:
        for pi,policy in enumerate(h.original.config()["estimator"]["contexts"]):
            context=boundary.context(policy,p,r);levels=context.levels()
            for draw in range(32):
                g=h.original.rng("estimate_params",nonce,draw);pp,rr=levels[int(g.integers(len(levels)))]
                yield "HFD04_ESTIMATOR_"+cohort+"_"+policy,i*32+draw,h.original.parameters(model,pp,rr,g),[1],5
    import h04_context_prior_estimator as supplemental
    weights=h.original.config() # separate frozen prior, never a chosen outcome
    weights=json.loads((h.FOUR/"config/context_prior_audit.json").read_text())["weights"]
    for i in range(81):
        for draw in range(32):
            g=supplemental.generator("cw_params",i,draw);level=int(g.choice(9,p=weights))
            yield "HFD04_PRIOR_RECONCILED_ESTIMATOR",i*32+draw,h.original.parameters(model,level//3,level%3,g),[1],5
    for cap in h.original.config()["effective_count"]["caps"]:
        empirical,demographic=h.original.isolated_models(cap);cap_model=empirical.fit(core)
        for i in range(16):
            yield "HFD04_CAP_"+str(cap),i,demographic.parameters(cap_model,h.original.rng("cap_state",i,1)),[1,.65,.85,.5],10
    import h04_juvenile as representation
    import h03_calibration as cal
    sources,_,_,_,_=h.original.load_data();held,path=representation.visible_inputs("M0")
    _,visible_model=cal.fit_visible(sources,path,held)
    for i in range(16):
        yield "HFD04_REPRESENTATION_COARSE_AND_ENRICHED",i,h.original.demography.parameters(visible_model,h.original.rng("representation_state",i,2)),[1,.5],10

def census():
    _,_,_,core,model=h.original.load_data();rows=[];totals=defaultdict(lambda:defaultdict(int))
    for family,identity,params,factors,years in parameter_sets(model,core):
        for si,factor in enumerate(factors):
            mu=params["growth"]*factor;sd=params["sd"]*factor;old=laws(mu,sd,False);new=laws(mu,sd,True)
            flags=dict(mu_zero=mu==0,sd_floor=sd<1e-9,mean_denominator_floor=(mu>0)&(mu<1e-9),
                positive_shape_floor=(mu>0)&((mu/np.maximum(sd,1e-9))**2<1e-6))
            affected=np.flatnonzero((old["mean"]!=new["mean"])|(old["variance"]!=new["variance"]))
            # Exact floating-law differences versus material moment defect are
            # distinguished; algebraically equal interior formulas can round.
            material=np.abs(old["mean"]-mu)>.01*np.maximum(mu,1e-300)
            total=totals[family];total["parameter_scenarios"]+=1;total["ancestral_cell_laws"]+=64
            for name,flag in flags.items():total[name]+=int(flag.sum())
            total["material_mean_cells"]+=int(material.sum())
            # Each original initial cell grows every year; four lowest-size
            # ancestral cells recur as birth cohorts in subsequent years.
            mult=np.full(64,years,int);mult[np.arange(16)*4]+=years*(years-1)//2
            total["cohort_year_laws"]+=int(mult.sum());total["material_cohort_year_laws"]+=int(mult[material].sum())
            for cell in np.flatnonzero(material|flags["mu_zero"]|flags["sd_floor"]|flags["mean_denominator_floor"]|flags["positive_shape_floor"]):
                rows.append(dict(family=family,parameter_identity=identity,scenario_index=si,growth_factor=factor,cell=int(cell),years=years,cohort_year_multiplicity=int(mult[cell]),
                    mu=float(mu[cell]),sd=float(sd[cell]),original_shape=float(old["shape"][cell]),original_scale=float(old["scale"][cell]),
                    original_mean=float(old["mean"][cell]),original_variance=float(old["variance"][cell]),corrected_mean=float(new["mean"][cell]),corrected_variance=float(new["variance"][cell]),
                    material_mean=bool(material[cell]),**{k:bool(v[cell]) for k,v in flags.items()}))
        if identity==0:print("Gamma branch census",family,flush=True)
    h.write_csv("gamma_growth/branch_census.csv",rows)
    h.write_csv("gamma_growth/branch_summary.csv",[dict(family=k,**v) for k,v in sorted(totals.items())])
    h.write_json("gamma_growth/branch_scope.json",dict(complete_for=list(totals),
        not_reenumerated=["Original pilots (excluded from confirmatory outcomes)","State-reconstruction recruitment Gamma is a distinct entry parameter law, not the audited growth law"],
        exclusions_status="Listed growth-parameter families complete; reconstruction/entry Gamma deliberately distinct",
        count_law_scope="Cohort-year parameter laws; multiply by recorded path count only for actual simulation calls",
        globals_mutated=False))

def impact():
    c=h.config("gamma_protocol.json")["matched_impact"];seed=h.original.frozen.Seeds();_,_,_,core,model=h.original.load_data()
    scenarios=h.original.frozen.load("config/design.json")["scenarios"];rows=[]
    for i,rep in itertools.product(c["primary_state_ids"],range(c["parameter_replicates"])):
        state=h.original.demography.reconstruct(core,model,seed.rng("reconstruct",i));params=h.original.demography.parameters(model,seed.rng("parameters",i,rep))
        for si,scenario in enumerate(scenarios):
            outputs=[]
            for corrected in [False,True]:
                bank,lawful=simulate(state,params,scenario,20,c["paths"],h.rng("gamma_impact",i,rep,si),corrected=corrected)
                outputs.append((bank,lawful));
            for horizon in c["horizons"]:
                comparisons=[]
                for corrected,(bank,lawful) in enumerate(outputs):
                    rr=h.original.demography.response(bank,lawful,model,horizon,.75,.4,2,1)
                    nv=int(rr["viable"].sum());nr=int(rr["reserve_viable"].sum());d=classify(rr["present"],nv,nr,c["paths"],.75,"Q2")
                    comparisons.append(dict(nv=nv,nr=nr,present=rr["present"],health=d["FINITE_BANK_HEALTH"],status=d["KERNEL_MC_STATUS"],bank_sha256=hashlib.sha256(bank.tobytes()+lawful.tobytes()).hexdigest()))
                a,b=comparisons
                rows.append(dict(state=i,replicate=rep,scenario=scenario["id"],horizon=horizon,K=c["paths"],
                    **{"original_"+k:v for k,v in a.items()},**{"corrected_"+k:v for k,v in b.items()},
                    viable_delta=(b["nv"]-a["nv"])/c["paths"],reserve_delta=(b["nr"]-a["nr"])/c["paths"],health_changed=a["health"]!=b["health"],
                    original_reference="ORIGINAL",corrected_reference="CORRECTED",states_parameters_context_matched=True))
        print("Matched Gamma impact",i,rep,flush=True)
    h.write_csv("gamma_growth/paired_generator_impact.csv",rows)
    h.write_json("gamma_growth/impact_summary.json",dict(cells=len(rows),max_abs_viable_delta=max(abs(r["viable_delta"]) for r in rows),
        max_abs_reserve_delta=max(abs(r["reserve_delta"]) for r in rows),health_changes=sum(r["health_changed"] for r in rows),
        scope="Fixed matched subset, not full original primary surface; all 256 registered horizon cells retained"))

if __name__=="__main__":
    import sys
    h.guard("gamma");{"units":unit_grid,"census":census,"impact":impact}[sys.argv[1]]()
