"""Pure numerical-boundary correction from existing exact sufficient counts.

No changed fit, newly sampled state, predictive history, outcome selection or seed.
All point masses/Health, real calibration, entry and E2 results stay unchanged.
"""
import ast
import hashlib
import json
import itertools
import numpy as np
from h03_common import ROOT,RUN,PROTOCOL,PROTOCOL_HASH,frozen,empirical,demographic,reconstruct,wilson,classify,read_csv,write_csv,write_json


def run(expected_functions):
    source=(ROOT/"scripts/h03_common.py").read_text()
    functions={node.name:hashlib.sha256(ast.dump(node,include_attributes=False).encode()).hexdigest()
               for node in ast.parse(source).body if isinstance(node,(ast.FunctionDef,ast.ClassDef)) and node.name!="wilson"}
    if functions!=expected_functions:
        raise ValueError("Non-Wilson model/observation/reconstruction/seed function changed")
    import h03_missingness
    core=empirical.read_core(frozen.verify_frozen_inputs())
    model=empirical.fit(core)
    old=frozen.load("config/design.json")
    cells=h03_missingness.query_grid(old)
    phases=read_csv(ROOT/"missingness/phase_surface.csv")
    index={(row["family"],row["scenario"],int(row["horizon"]),float(row["alpha"]),float(row["beta"]),int(row["tau"]),
            row["query"],float(row["rho"]),float(row["theta"])):row for row in phases}
    K=PROTOCOL["missingness"]["future_draws"]
    invariant_hashes={}
    for family in ["B0","B1","B2","B3","B4"]:
        path=RUN/("missingness_"+family+".npz")
        with np.load(path,allow_pickle=False) as bank:
            masses,reserves,points=bank["mass"],bank["reserve"],bank["point"]
            original_status=bank["status"]
        invariant_hashes[family]={
            key:hashlib.sha256(value.tobytes()).hexdigest()
            for key,value in [("mass",masses),("reserve",reserves),("point",points)]}
        p=np.empty((masses.shape[1],len(cells)),bool)
        for i in range(PROTOCOL["missingness"]["state_draws"]):
            state=reconstruct(core,model,i,family)
            initial=demographic.metrics(state["n"][None,:],state["d"][None,:],np.arange(64),
                                        np.zeros(1),np.zeros(1))[0]
            for j,(_,a,b,_,_,_,_) in enumerate(cells):
                present=initial[1]/model["baseline_ba"]>=a and initial[2]/model["baseline_juveniles"]>=b
                p[i*2:(i+1)*2,j]=present
        theta=np.array([cell[-1] for cell in cells])
        q2=np.array([cell[4]=="Q2" for cell in cells])
        statuses=np.empty(original_status.shape,np.uint8)
        for s,scenario in enumerate(old["scenarios"]):
            v=np.rint(masses[s]*K).astype(int);r=np.rint(reserves[s]*K).astype(int)
            vl,vu=wilson(v,K,.95);ql,qu=wilson(v,K,.975);rl,ru=wilson(r,K,.975)
            lower=np.where(q2,np.minimum(ql,rl),vl);upper=np.where(q2,np.minimum(qu,ru),vu)
            z=np.where(~p|(upper<theta),0,np.where(lower>theta,1,2))
            point=p&(v/K>=theta)&(~q2|(r/K>=theta))
            if not np.array_equal(point,points[s]):
                raise ValueError("Point Health changed during interval-only recount")
            statuses[s]=z
            for j,cell in enumerate(cells):
                row=index[(family,scenario["id"],*cell)]
                row.update(true_units=int((z[:,j]==1).sum()),false_units=int((z[:,j]==0).sum()),
                           unresolved_units=int((z[:,j]==2).sum()),
                           resolved_health_lower=float((z[:,j]==1).mean()),
                           resolved_health_upper=float((z[:,j]!=0).mean()))
        np.savez_compressed(path,mass=masses,reserve=reserves,point=points,status=statuses)
    write_csv("missingness/phase_surface.csv",phases)
    inner=read_csv(ROOT/"monte_carlo/inner_resolution.csv")
    for row in inner:
        count=[int(row["viable_count"])]
        if row["query"]=="Q2":
            count.append(int(row["reserve_count"]))
        status,lower,upper=classify(row["present"]=="True",count,int(row["K"]),float(row["theta"]))
        row.update(status=status,lower=float(min(lower)),upper=float(min(upper)),
                   max_half_width=float(np.max((upper-lower)/2)))
    write_csv("monte_carlo/inner_resolution.csv",inner)
    truth=read_csv(ROOT/"synthetic_validation/truth_census.csv")
    performance=read_csv(ROOT/"synthetic_validation/performance.csv")
    indexed={(int(row["world"]),row["suite"]):row for row in performance}
    Ktruth=PROTOCOL["synthetic"]["truth_paths"]
    for row in truth:
        counts=[int(round(float(row[key])*Ktruth)) for key in ["truth_mass","truth_reserve"]]
        status,lower,upper=classify(row["truth_present"]=="True",counts,Ktruth,PROTOCOL["synthetic"]["theta"])
        row.update(truth_mc_status=status,truth_lower=float(min(lower)),truth_upper=float(min(upper)))
        indexed[(int(row["world"]),row["suite"])].update(
            truth_mc_status=status,truth_lower=float(min(lower)),truth_upper=float(min(upper)))
    write_csv("synthetic_validation/truth_census.csv",truth)
    write_csv("synthetic_validation/performance.csv",performance)
    census=json.loads((ROOT/"synthetic_validation/pre_estimator_census.json").read_text())
    census["truth_census_sha256"]=frozen.sha(ROOT/"synthetic_validation/truth_census.csv")
    write_json("synthetic_validation/pre_estimator_census.json",census)
    write_json("provenance/interval_boundary_amendment.json",dict(
        protocol_sha256=PROTOCOL_HASH,parameters_and_confidence_rules_unchanged=True,
        model_replay_receipt="reports/full_model_replay_receipt.json",
        unchanged_non_wilson_function_hashes=functions,
        sufficient_mass_point_hashes=invariant_hashes,
        reason="IEEE roundoff yielded U(K,K)=0.9999999999999999; standard Wilson identity U(K,K)=1 restored, L(0,K)=0 enforced",
        impact="Only interval representation and boundary classifications; no mass, point Health, model, seed, calibration, entry or E2 change",
        selection=False,preferred_boolean_stopping=False))
    import h03_reports
    h03_reports.run(core,model)
    print("PASS interval-only recount: all underlying forecasts/masses/point Health unchanged")


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--function-hashes",required=True)
    args=parser.parse_args()
    run(json.loads(args.function_hashes))
