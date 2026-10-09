"""Resolve rare predictive tails without FFT subtraction or Monte Carlo zeros."""
import math
import numpy as np
import h06_common as h
from h06_mortality import mixture_pmf, discrete_summary

def run():
    h.guard()
    model = h.model()
    q = -np.expm1(-3*model["hazard_grid"])
    cells = h.read_csv(h.ROOT/"mortality/cell_census.csv")
    direct = h.read_csv(h.ROOT/"checks/predictive_pmf_direct.csv")
    checks = h.read_csv(h.ROOT/"checks/predictive_independent.csv")
    original = h.read_csv(h.ROOT/"mortality/predictive_summary.csv")
    final,comparison = [],[]
    for r in original:
        name = r["group"]
        selected = [c for c in cells if int(c["n"]) and (name=="ALL" or c["group"]==name)]
        pmf = np.ones(1)
        for c in sorted(selected,key=lambda c:int(c["n"])):
            c,n = int(c["cell"]),int(c["n"])
            pmf = np.convolve(pmf,mixture_pmf(n,model["hazard_posterior"][c],q))
        other = np.array([float(v["probability"]) for v in direct if v["group"]==name])
        observed = int(r["observed_deaths"])
        tail = math.fsum(float(v) for v in pmf[observed:])
        tail_checked = next(float(v["recurrence_direct_tail"]) for v in checks if v["group"]==name)
        relative_error = abs(tail/tail_checked-1)
        if relative_error > 1e-9 or np.max(abs(pmf-other))>2e-12:
            raise ValueError("Independent positive-convolution rare tail disagreement")
        numeric = {k:float(v) for k,v in r.items() if k in ["predictive_mean","conditional_event_variance","shared_hazard_variance","total_predictive_variance","shared_hazard_variance_fraction","MC_mean","MC_tail","MC_tail_Wilson_lower","MC_tail_Wilson_upper","MC_variance"]}
        integer = {k:int(r[k]) for k in ["adult_stems","observed_deaths","MC_replicates","MC_median","MC_lower95","MC_upper95","MC_tail_hits"]}
        final.append(dict(group=name,**integer,**numeric,**discrete_summary(other,observed),
                          tail_numerical_method="Positive recurrence/direct convolution; checked against loggamma/direct convolution", tail_check_relative_error=relative_error))
        comparison.append(dict(group=name, raw_FFT_tail=float(r["upper_tail"]), positive_recurrence_tail=tail_checked,
                               positive_loggamma_tail=tail, relative_positive_tail_error=relative_error,
                               raw_FFT_negative_roundoff_mass=float(r["fft_negative_roundoff_mass"]),
                               raw_FFT_tail_not_authoritative=name=="HEMLOCK", sampler_or_model_changed=False))
    h.write_csv("mortality/final_predictive_summary.csv",final)
    h.write_csv("checks/rare_tail_precision.csv",comparison)
    h.write_json("checks/rare_tail_precision.json",dict(status="PASS", both_independent_positive_routes=True,
                 tolerance_relative=1e-9, new_hazards_or_events_sampled=False, original_FFT_and_sampler_outputs_preserved=True))
    print("Positive rare-tail precision PASS",flush=True)

if __name__ == "__main__":
    run()
