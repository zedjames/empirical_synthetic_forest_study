"""Frozen E2 comparison plus predeclared perturbation surfaces, no fitting."""
import copy
import itertools
import math
import numpy as np
from h03_common import PROTOCOL, ROOT, frozen, demographic, reconstruct, rng, write_csv, write_json


def two_point_bank(state,params,scenario,years,K,generator):
    # Split each initial size cohort into two deterministic subcohorts, run
    # each independent original kernel with half recruitment, then add counts,
    # BA/deaths/entries. Juvenile support can differ; no size-dependent mortality
    # exists in the old kernel to discover through splitting.
    first=copy.deepcopy(state)
    second=copy.deepcopy(state)
    first["n"]=state["n"]//2
    second["n"]=state["n"]-first["n"]
    n=state["n"]
    a=first["n"]
    b=second["n"]
    denominator=np.divide(.75**2*a+1.25**2*b,n,out=np.ones(64),where=n>0)
    first["d"]=np.clip(state["d"]*.75/np.sqrt(denominator),1,300)
    second["d"]=np.clip(state["d"]*1.25/np.sqrt(denominator),1,300)
    # Clamp can disrupt exact BA at extreme sizes: expose measured difference.
    p=copy.deepcopy(params)
    p["recruit"]=params["recruit"]/2
    x,l1=demographic.simulate(first,p,scenario,years,K,generator)
    y,l2=demographic.simulate(second,p,scenario,years,K,generator)
    return x+y,l1&l2


def run(core,model,audit):
    cfg=PROTOCOL["diagnostics"]
    rows=[]
    reference_rate=audit["annual_candidate_rates"]["high"]
    # Scale the observed-candidate-based rate AND its inherited wet borrowing.
    # This extrapolation is declared sensitivity, not identified wet entry.
    factors={k:v/reference_rate for k,v in audit["annual_candidate_rates"].items()}
    default=dict(entry="high",hemlock_extra_hazard=0,general_hazard_factor=1,
                 wet_count_factor=1,duration_years=5,resolution="cohort_rms")
    designs={}
    for axis in ("entry","hemlock_extra_hazard","general_hazard_factor","wet_count_factor","duration_years","resolution"):
        for value in cfg[axis]:
            design=dict(default);design[axis]=value
            designs[tuple(design.items())]=design
    for e,h,g in itertools.product(cfg["entry"],cfg["hemlock_extra_hazard"],cfg["general_hazard_factor"]):
        design=dict(default,entry=e,hemlock_extra_hazard=h,general_hazard_factor=g)
        designs[tuple(design.items())]=design
    for index,design in enumerate(designs.values()):
        all_end,hem_deaths,juveniles,initial_ba_errors=[],[],[],[]
        for i in range(cfg["states"]):
            state=reconstruct(core,model,i)
            wet=(np.arange(64)//4)%4==3
            state["n"][wet]=np.rint(state["n"][wet]*design["wet_count_factor"]).astype(int)
            params=demographic.parameters(model,rng("diagnostic",i,0))
            params["recruit"]*=factors[design["entry"]]
            scenario=dict(id="diagnostic",hazard_factor=design["general_hazard_factor"],growth_factor=1,
                          recruitment_factor=1,hemlock_extra_hazard=design["hemlock_extra_hazard"])
            years=design["duration_years"]
            K=cfg["futures"]
            simulate=two_point_bank if design["resolution"]=="two_point_size" else demographic.simulate
            bank,lawful=simulate(state,params,scenario,years,K,rng("diagnostic",i,1))
            if not lawful.all():
                raise ValueError("Diagnostic history inadmissible")
            all_end.extend(bank[:,-1,[4,5]].tolist())
            juveniles.extend(bank[:,-1,2].tolist())
            ba0=float(np.sum(state["n"]*state["d"]**2)*math.pi/40000)
            initial_ba_errors.extend((bank[:,0,1]-ba0).tolist())
            # Separate ancestral hemlock survivor simulation, not BA as a
            # mortality proxy. Birth deaths are excluded from this anchor.
            counts=rng("diagnostic",i,2).binomial(
                state["n"][:16],np.exp(-(params["hazard"][:16]*scenario["hazard_factor"]
                                           +scenario["hemlock_extra_hazard"])*years),size=(K,16))
            hem_deaths.extend((state["n"][:16].sum()-counts.sum(axis=1)).tolist())
        x=np.array(all_end)
        rows.append(dict(design_id=index,**design,entry_factor=factors[design["entry"]],
                         predicted_entries=float(np.median(x[:,1])),
                         entries_p05=float(np.quantile(x[:,1],.05)),entries_p95=float(np.quantile(x[:,1],.95)),
                         predicted_all_deaths=float(np.median(x[:,0])),
                         predicted_ancestral_hemlock_deaths=float(np.median(hem_deaths)),
                         hemlock_p05=float(np.quantile(hem_deaths,.05)),hemlock_p95=float(np.quantile(hem_deaths,.95)),
                         predicted_juveniles=float(np.median(juveniles)),
                         initial_ba_error_max=float(np.max(np.abs(initial_ba_errors))),
                         public_entries_approx=5000,public_hemlock_approx=5000,public_all_deaths_lower_bound=11800,
                         comparison_status="DEFINITIONALLY_INCOMPARABLE_PUBLIC_AGGREGATES; directional diagnosis only"))
        if index%8==0:
            print("E2 diagnostic design",index,flush=True)
    write_csv("e2_failure_localization/response_surfaces.csv",rows)
    # Frozen actual comparison file is copied semantically into new provenance,
    # never edited and never substituted by any newly favorable factor setting.
    external=ROOT.parent/"HFD02S/models/e2_comparison.json"
    candidates=list((ROOT.parent/"HFD02S").rglob("*e2*.json"))
    frozen_comparison=[dict(path=str(p.relative_to(ROOT.parent/"HFD02S")),sha256=frozen.sha(p),
                            data=__import__("json").loads(p.read_text())) for p in candidates if ".runs" not in str(p)]
    write_json("e2_failure_localization/summary.json",dict(
        frozen_comparison=frozen_comparison,
        frozen_d0=dict(entries_median=9073,hemlock_deaths_median=2421,all_deaths_median=15296),
        aligned_failure_vector=None,
        public_scope="Approximate all-stem death lower bound, newly recorded stems, approximate hemlock dead since2020; date, frame, size and identity scope not harmonized.",
        residual_policy="No exact signed residual is licensed. Anchors locate diagnostic magnitude, not validation targets.",
        entries_causes=["ENTRY_IDENTIFICATION","SAMPLING_FRAME","OBSERVATION_DEFINITION","MODEL_STRUCTURE","TEMPORAL_ALIGNMENT","UNRESOLVED"],
        hemlock_causes=["SCENARIO","PARAMETER","MODEL_STRUCTURE","STATE_RECONSTRUCTION","TEMPORAL_ALIGNMENT","UNRESOLVED"],
        resolution_limit="Splitting sizes cannot create species hazard absent from size-independent kernel; juvenile threshold changes, not a new individual ecological model.",
        original_replaced=False,parameter_selection=False,diagnostic_designs=len(rows)))
    return rows
