"""Shared-cell hazards, explicit individual deaths and finite-grid predictive law."""
import hashlib
import math
from collections import Counter
import numpy as np
import h06_common as h

def risk_set():
    cfg = h.config()["mortality"]
    rows = [r for r in h.read_csv(h.FIVE / "external_sources/HF453/mortality_predictions.csv")
            if r["target"] == cfg["target"] and int(r["start_year"]) == cfg["start_year"]
            and int(r["end_year"]) == cfg["end_year"] and r["scheme"] == cfg["scheme"]]
    if len(rows) != 21576 or len({r["external_tag"] for r in rows}) != len(rows):
        raise ValueError("Frozen adult matched denominator mismatch")
    for r in rows:
        if float(r["exposure_years"]) != 3 or r["raw_endpoint_status"] not in ["A","AU","DC","DS"] or r["model_refitted"] != "False":
            raise ValueError("Frozen adult support mismatch")
        if (int(r["cell"]) < 16) != (r["hemlock"] == "HEMLOCK"):
            raise ValueError("Taxon/cell mismatch")
    return rows

def moments(n, posterior, probability):
    eq = np.einsum("cj,j->c", posterior, probability)
    eq2 = np.einsum("cj,j->c", posterior, probability**2)
    # E Var(D|H) + Var E(D|H). Cells are independent; stems within a cell are not.
    return dict(mean=float(np.sum(n*eq)), event_variance=float(np.sum(n*(eq-eq2))),
                hazard_variance=float(np.sum(n*n*(eq2-eq*eq))))

def mixture_pmf(n, weights, probability):
    k = np.arange(n+1)
    coefficient = np.array([math.lgamma(n+1)-math.lgamma(int(v)+1)-math.lgamma(n-int(v)+1) for v in k])
    logs = coefficient[None,:] + k[None,:]*np.log(probability)[:,None] + (n-k)[None,:]*np.log1p(-probability)[:,None]
    pmf = np.einsum("j,jk->k", weights, np.exp(logs))
    if abs(float(pmf.sum())-1) > 1e-10:
        raise ValueError("Cell predictive PMF not normalized")
    return pmf / pmf.sum()

def convolve_fft(distributions):
    length = 1 + sum(len(p)-1 for p in distributions)
    size = 1 << (length-1).bit_length()
    transform = np.ones(size//2+1, dtype=complex)
    for pmf in distributions:
        transform *= np.fft.rfft(pmf, n=size)
    raw = np.fft.irfft(transform, n=size)[:length]
    negative_mass = float(-raw[raw<0].sum())
    if negative_mass > 1e-10: raise ValueError("Predictive convolution unstable")
    pmf = np.maximum(raw,0)
    return pmf/pmf.sum(), negative_mass

def discrete_summary(pmf, observed):
    cumulative = np.cumsum(pmf)
    median, lower, upper = [int(np.searchsorted(cumulative, p)) for p in [.5,.025,.975]]
    return dict(predictive_median=median, predictive_lower95=lower, predictive_upper95=upper,
                upper_tail=float(math.fsum(float(v) for v in pmf[observed:])))

def run():
    h.guard()
    cfg = h.config()["mortality"]
    risk = risk_set()
    model = h.model()
    grid, posterior = model["hazard_grid"], model["hazard_posterior"]
    if posterior.shape != (64,256) or not np.allclose(posterior.sum(axis=1),1,rtol=0,atol=2e-15):
        raise ValueError("Invalid frozen hazard distribution")
    counts = np.bincount([int(r["cell"]) for r in risk], minlength=64)
    deaths = np.bincount([int(r["cell"]) for r in risk], weights=[1-int(r["observed_survival"]) for r in risk], minlength=64).astype(int)
    probability = -np.expm1(-grid*cfg["exposure_years"])
    cell_pmfs = {c:mixture_pmf(int(n),posterior[c],probability) for c,n in enumerate(counts) if n}
    cell_rows, risk_rows = [], []
    for r in risk:
        risk_rows.append({k:r[k] for k in ["external_tag","linked_stem","cell","hemlock","exposure_years","observed_survival","predicted_survival","raw_start_status","raw_endpoint_status"]})
    for c,n in enumerate(counts):
        m = moments(np.array([n]),posterior[c:c+1],probability)
        cell_rows.append(dict(cell=c, group="HEMLOCK" if c<16 else "NONHEMLOCK", n=int(n), observed_deaths=int(deaths[c]),
                              predicted_deaths=m["mean"], conditional_event_variance=m["event_variance"], shared_hazard_variance=m["hazard_variance"]))
    h.write_csv("mortality/matched_risk_set.csv", risk_rows)
    h.write_csv("mortality/cell_census.csv", cell_rows)
    # Store hazard identities and cell counts; actually generate individual events.
    R, B = cfg["replicates"], cfg["chunk_size"]
    sampled_indices = np.zeros((R,64),dtype=np.uint16)
    sampled_deaths = np.zeros((R,64),dtype=np.uint16)
    draws, digest_rows = [], []
    for start in range(0,R,B):
        chunk = start//B
        hazards, events = h.rng("hazards",chunk), h.rng("individual_events",chunk)
        digest = hashlib.sha256()
        for c,n in enumerate(counts):
            if not n: continue
            indices = hazards.choice(256, p=posterior[c], size=B)
            sampled_indices[start:start+B,c] = indices
            # One uniform per actual matched stem, conditional on the shared draw.
            died = events.random((B,int(n))) < probability[indices,None]
            digest.update(np.packbits(died,axis=None).tobytes())
            sampled_deaths[start:start+B,c] = died.sum(axis=1)
        digest_rows.append(dict(chunk=chunk, replicates=B, ordered_individual_event_bits_sha256=digest.hexdigest()))
        if chunk % 64 == 0: print("Individual mortality replicates",start+B,"of",R,flush=True)
    np.savez_compressed(h.output("mortality/replicate_archive.npz"), hazard_indices=sampled_indices, cell_deaths=sampled_deaths, cell_n=counts)
    h.write_csv("mortality/individual_event_digests.csv", digest_rows)
    group_masks = {"HEMLOCK":np.arange(64)<16, "NONHEMLOCK":np.arange(64)>=16, "ALL":np.ones(64,dtype=bool)}
    group_draws = {name:sampled_deaths[:,mask].sum(axis=1).astype(int) for name,mask in group_masks.items()}
    for i in range(R):
        draws.append(dict(replicate=i, hemlock=int(group_draws["HEMLOCK"][i]), nonhemlock=int(group_draws["NONHEMLOCK"][i]), all_adults=int(group_draws["ALL"][i])))
    h.write_csv("mortality/predictive_replicates.csv", draws)
    summaries, pmf_rows, histograms, variance_rows = [], [], [], []
    for name,mask in group_masks.items():
        observed = int(deaths[mask].sum())
        values = group_draws[name]
        pmf, negative_mass = convolve_fft([cell_pmfs[c] for c in np.flatnonzero(mask & (counts>0))])
        m = moments(counts[mask],posterior[mask],probability)
        exact = discrete_summary(pmf,observed)
        total_variance = m["event_variance"] + m["hazard_variance"]
        v = np.arange(len(pmf))
        if abs(float(np.sum(v*pmf))-m["mean"]) > 1e-7 or abs(float(np.sum((v-m["mean"])**2*pmf))-total_variance) > 1e-5:
            raise ValueError("Predictive PMF/moment mismatch")
        tails = int((values>=observed).sum())
        lo,hi = h.prior.original.prior.wilson(np.array([tails]),R,.95)
        quantiles = np.quantile(values,[.5,.025,.975],method="inverted_cdf")
        summaries.append(dict(group=name, adult_stems=int(counts[mask].sum()), observed_deaths=observed,
                              predictive_mean=m["mean"], **exact, conditional_event_variance=m["event_variance"],
                              shared_hazard_variance=m["hazard_variance"], total_predictive_variance=total_variance,
                              shared_hazard_variance_fraction=m["hazard_variance"]/total_variance,
                              MC_replicates=R, MC_mean=float(values.mean()), MC_median=int(quantiles[0]),
                              MC_lower95=int(quantiles[1]), MC_upper95=int(quantiles[2]), MC_tail_hits=tails,
                              MC_tail=tails/R, MC_tail_Wilson_lower=float(lo[0]), MC_tail_Wilson_upper=float(hi[0]),
                              MC_variance=float(values.var()), fft_negative_roundoff_mass=negative_mass,
                              observed_inside_central95=exact["predictive_lower95"]<=observed<=exact["predictive_upper95"]))
        for k,p in enumerate(pmf): pmf_rows.append(dict(group=name,deaths=k,probability=float(p)))
        for k,n in sorted(Counter(values.tolist()).items()): histograms.append(dict(group=name,deaths=k,replicates=n))
        conditional_means = np.sum(counts[mask]*probability[sampled_indices[:,mask]],axis=1)
        conditional_variance = np.sum(counts[mask]*probability[sampled_indices[:,mask]]*(1-probability[sampled_indices[:,mask]]),axis=1)
        variance_rows.append(dict(group=name, analytic_event_variance=m["event_variance"], sampled_conditional_event_variance=float(conditional_variance.mean()),
                                  analytic_shared_hazard_variance=m["hazard_variance"], sampled_shared_conditional_mean_variance=float(conditional_means.var()),
                                  sampled_event_residual_mean=float(np.mean(values-conditional_means)), sampled_event_residual_second_moment=float(np.mean((values-conditional_means)**2))))
        print(name, summaries[-1], flush=True)
    h.write_csv("mortality/predictive_summary.csv", summaries)
    h.write_csv("mortality/predictive_pmf.csv", pmf_rows)
    h.write_csv("mortality/predictive_histogram.csv", histograms)
    h.write_csv("mortality/variance_decomposition.csv", variance_rows)
    h.write_json("mortality/completion.json",dict(status="COMPLETE", scope="STRICT matched2021 adult survivors, cumulative2024, fixed3year July anniversary", replicates=R,
                 individual_events_generated=R*int(counts.sum()), fitted_distribution_sha256=h.sha(h.TWO/"models/fitted_parameters.json"),
                 no_external_refit=True, marginal_independence_not_assumed=True, full_Health_validation=False, public_forecast=False))

if __name__ == "__main__":
    run()
