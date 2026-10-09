"""Independent analytic audit of completed mathematical phases; no parameter fitting."""
import math,json
from collections import defaultdict
import numpy as np
import h05_common as h

def run():
    for phase in ["q2","gamma","gamma_reached","exposure","exposure_coupling"]:h.guard(phase)
    rows=h.read_csv(h.ROOT/"gamma_growth/moment_unit_tests.csv")
    for r in rows:
        mu,s=map(float,[r["mu"],r["sd"]]);corrected=r["law"]=="CORRECTED"
        if corrected:
            mean=mu;var=0 if mu==0 or s==0 else min(s*s,mu*mu/1e-6)
        else:
            shape=max((mu/max(s,1e-9))**2,1e-6);scale=s*s/max(mu,1e-9) if mu>0 else 0;mean=shape*scale;var=shape*scale*scale
        for expected,name in [(mean,"analytic_mean"),(var,"analytic_variance")]:
            assert math.isclose(expected,float(r[name]),rel_tol=1e-12,abs_tol=1e-25)
        assert r["mean_bound_pass"]==r["variance_bound_pass"]=="True"
    calls=h.read_csv(h.ROOT/"gamma_growth/actual_call_unit_inventory.csv")
    for r in calls:
        mu,s,a,b=map(float,[r["mu"],r["sd"],r["actual_shape"],r["actual_scale"]])
        assert a==max((mu/max(s,1e-9))**2,1e-6) and b==(s*s/max(mu,1e-9) if mu>0 else 0)
    exposures=h.read_csv(h.ROOT/"exposure_model/individual_exposure_records.csv");groups=defaultdict(list)
    for r in exposures:
        assert r["date_fallback"]=="NONE" and r["observed_date"] not in ["","NA"]
        groups[(r["component"],r["cell"])].append(float(r["exposure_years"]))
    inventory=h.read_csv(h.ROOT/"exposure_model/cohort_exposure_inventory.csv")
    for r in inventory:
        times=np.array(groups[(r["component"],r["cell"])]);hazard=float(r["conditional_hazard"]);p=np.exp(-hazard*times)
        mean=float(p.sum());var=float(np.dot(p,1-p));pooled=len(times)*math.exp(-hazard*float(times.mean()))
        assert len(times)==int(r["n"])
        assert abs(mean-float(r["PB_mean"]))<1e-8 and abs(var-float(r["PB_variance"]))<1e-8
        assert abs(mean-pooled-float(r["Jensen_gap"]))<1e-8 and mean>=pooled-1e-8
    dp=defaultdict(list)
    for r in h.read_csv(h.ROOT/"exposure_model/poisson_binomial_dp.csv"):dp[(r["component"],r["cell"])].append(r)
    for key,values in dp.items():
        probabilities=np.array([float(r["PB_probability"]) for r in values]);times=np.array(groups[key])[:int(values[0]["n"])];p=np.exp(-.12*times)
        counts=np.array([int(r["count"]) for r in values]);mean=float(np.dot(counts,probabilities));variance=float(np.dot((counts-mean)**2,probabilities))
        assert abs(probabilities.sum()-1)<1e-12 and abs(mean-p.sum())<1e-12 and abs(variance-np.dot(p,1-p))<1e-12
    controlled=h.read_csv(h.ROOT/"exposure_model/controlled_origin_comparison.csv")
    assert len(controlled)==256+5*128
    for r in controlled:
        assert int(r["N_delta"])==int(r["individual_dated_survivors"])-int(r["original_dated_survivors"])
        assert float(r["conditional_individual_survivor_mean"])>=float(r["conditional_original_survivor_mean"])-1e-8
    h.write_json("mathematical_audit/abc_independent_audit.json",dict(status="PASS",gamma_moment_cells=len(rows),actual_historical_call_cells=len(calls),individual_dated_records=len(exposures),conditional_exposure_cells=len(inventory),exact_DP_distributions=len(dp),controlled_state_draws=len(controlled),
        source_sha256=h.sha(__file__),scope="Mathematics and specified finite audits; downstream factorial, oracle, external scoring and final closure remain separate requirements"))
    print("PASS independent Gamma/exposure analytic audit",flush=True)
if __name__=="__main__":run()
