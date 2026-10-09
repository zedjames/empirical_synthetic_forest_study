"""Independent scalar attribution and recurrence-based predictive distribution."""
import hashlib
import itertools
import math
from collections import Counter
import numpy as np
import h06_common as h

def require(condition, message):
    if not condition: raise ValueError(message)

def binomial_recurrence(n, q):
    """Normalize relative masses anchored at the mode, never factorial/loggamma."""
    if n==0: return np.ones(1)
    if q==0: return np.r_[1.,np.zeros(n)]
    if q==1: return np.r_[np.zeros(n),1.]
    mode = min(n,int(math.floor((n+1)*q)))
    p = np.zeros(n+1)
    p[mode] = 1
    for k in range(mode,n): p[k+1] = p[k]*(n-k)/(k+1)*q/(1-q)
    for k in range(mode,0,-1): p[k-1] = p[k]*k/(n-k+1)*(1-q)/q
    return p / math.fsum(p)

def event_enumeration(counts, hazard_laws):
    """Small fixture enumerates hazards AND every individual outcome bit."""
    result = np.zeros(sum(counts)+1)
    for hazards in itertools.product(*hazard_laws):
        weight = math.prod(item[1] for item in hazards)
        probabilities = [item[0] for n,item in zip(counts,hazards) for _ in range(n)]
        for bits in itertools.product([0,1], repeat=sum(counts)):
            result[sum(bits)] += weight*math.prod(p if bit else 1-p for p,bit in zip(probabilities,bits))
    return result

def fixtures():
    counts = [3,4]
    laws = [[(.1,.6),(.6,.4)],[(.2,.3),(.5,.7)]]
    brute = event_enumeration(counts,laws)
    cell = [sum(w*binomial_recurrence(n,q) for q,w in law) for n,law in zip(counts,laws)]
    recurrence = np.convolve(*cell)
    require(np.max(abs(brute-recurrence))<3e-15,"Individual shared-hazard fixture disagrees")
    marginal = np.convolve(*[binomial_recurrence(n,sum(q*w for q,w in law)) for n,law in zip(counts,laws)])
    require(np.max(abs(brute-marginal))>.01,"Fixture failed to detect incorrect independent marginal hazards")
    # Law of total variance checked against enumeration of the complete fixture.
    mean = sum(n*sum(q*w for q,w in law) for n,law in zip(counts,laws))
    event = sum(n*sum(w*q*(1-q) for q,w in law) for n,law in zip(counts,laws))
    hazard = sum(n*n*(sum(w*q*q for q,w in law)-sum(w*q for q,w in law)**2) for n,law in zip(counts,laws))
    require(abs(sum((k-mean)**2*p for k,p in enumerate(brute))-event-hazard)<1e-13,"Variance fixture fails")
    return dict(individual_patterns=128, joint_hazard_cases=4, PMF_max_error=float(np.max(abs(brute-recurrence))),
                detects_per_stem_hazard_redraw=True, detects_omitted_shared_variance=True)

def attribution():
    model = h.model()
    cfg = h.config()["terminal"]
    archive = np.load(h.ROOT/"attribution/path_archive.npz",allow_pickle=False)
    rows = h.read_csv(h.ROOT/"attribution/paths.csv")
    partition = Counter()
    overlaps = Counter()
    require(len(rows)==8192,"Attribution candidate census incomplete")
    for identity in ["POOLED","INDIVIDUAL"]:
        selected = [r for r in rows if r["reconstruction"]==identity]
        require([int(r["path"]) for r in selected]==list(range(4096)),"Lost or duplicated attribution candidate")
        for i,r in enumerate(selected):
            bank = archive[identity+"_bank"][i]
            ns = archive[identity+"_terminal_n"][i]
            ds = archive[identity+"_terminal_d"][i]
            require(int(ns.sum())==bank[-1,0],"Terminal captured stem census disagrees")
            require(sum(int(n) for n,d in zip(ns,ds) if d<10)==bank[-1,2],"Terminal captured juvenile census disagrees")
            require(abs(math.fsum(float(n)*float(d)**2 for n,d in zip(ns,ds))*math.pi/40000-bank[-1,1])<1e-9,"Terminal captured basal area disagrees")
            lawful = all(math.isfinite(float(v)) and v>=0 for v in bank.flat)
            lawful = lawful and all(math.isfinite(float(d)) and 1<=d<=300 for d in ds)
            lawful = lawful and all(n>=0 and n==int(n) for n in ns)
            lawful = lawful and all(b[2]<=b[0] for b in bank)
            lawful = lawful and bank[1,0]-bank[0,0]==-(bank[1,4]-bank[0,4])+(bank[1,5]-bank[0,5])
            lawful = lawful and bank[1,4]>=bank[0,4] and bank[1,5]>=bank[0,5]
            require(lawful==bool(archive[identity+"_lawful"][i]),"Independent model lawfulness differs")
            ba = bank[-1,1]/model["baseline_ba"] < cfg["alpha"]
            ju = bank[-1,2]/model["baseline_juveniles"] < cfg["beta"]
            require(bank[0,1]/model["baseline_ba"]>=cfg["alpha"] and bank[0,2]/model["baseline_juveniles"]>=cfg["beta"],"Nonterminal excursion unaccounted")
            outcome = "MODEL_UNLAWFULNESS" if not lawful else ("JOINT_JUVENILE_BASAL_AREA" if ba and ju else ("JUVENILE_ONLY" if ju else ("BASAL_AREA_ONLY" if ba else "SUCCESS")))
            require(outcome==r["outcome"] and str(not ba and not ju and lawful)==r["viable"],"Independent attribution differs")
            require(str(ju)==r["terminal_juvenile_failure"] and str(ba)==r["terminal_basal_area_failure"],"Overlapping causes lost")
            partition[(identity,outcome)]+=1
            for name,condition in [("ANY_UNLAWFUL",not lawful),("ANY_JUVENILE_FAILURE",ju),("ANY_BASAL_AREA_FAILURE",ba),("JUVENILE_AND_BASAL",ju and ba),
                                   ("UNLAWFUL_AND_JUVENILE",not lawful and ju),("UNLAWFUL_AND_BASAL",not lawful and ba),("ALL_THREE",not lawful and ju and ba)]:
                overlaps[(identity,name)]+=int(condition)
        require(partition[(identity,"SUCCESS")]==cfg["expected_viable"][identity],"Frozen success count not reproduced")
    for r in h.read_csv(h.ROOT/"attribution/partition.csv"):
        require(partition[(r["reconstruction"],r["outcome"])]==int(r["paths"]),"Partition total corrupted")
    for r in h.read_csv(h.ROOT/"attribution/overlapping_causes.csv"):
        require(overlaps[(r["reconstruction"],r["cause"])]==int(r["paths"]),"Overlap total corrupted")
    h.write_json("checks/attribution.json",dict(status="PASS", candidates_checked=8192, scalar_independent=True,
                 archived_cohort_sizes_checked=True, exact_reconciliation=True, overlapping_causes_checked=True))

def predictive():
    cfg = h.config()["mortality"]
    risk = h.read_csv(h.ROOT/"mortality/matched_risk_set.csv")
    frozen = [r for r in h.read_csv(h.FIVE/"external_sources/HF453/mortality_predictions.csv")
              if r["target"]=="CUMULATIVE" and r["start_year"]=="2021" and r["end_year"]=="2024" and r["scheme"]=="STRICT"]
    require(risk==[{k:r[k] for k in risk[0]} for r in frozen],"Matched risk set selection changed")
    counts = Counter(int(r["cell"]) for r in risk)
    observed = Counter()
    for r in risk: observed[int(r["cell"])]+=1-int(r["observed_survival"])
    model = h.model()
    grid, weights = model["hazard_grid"],model["hazard_posterior"]
    qs = np.array([1-math.exp(-float(v)*3) for v in grid])
    cell_pmfs, scalar = {}, {}
    for c,n in sorted(counts.items()):
        pmf = np.zeros(n+1)
        for q,w in zip(qs,weights[c]): pmf += w*binomial_recurrence(n,float(q))
        require(abs(math.fsum(pmf)-1)<3e-14,"Independent cell PMF normalization failed")
        cell_pmfs[c] = pmf
        mean_q = math.fsum(float(w)*float(q) for w,q in zip(weights[c],qs))
        event = n*math.fsum(float(w)*float(q)*(1-float(q)) for w,q in zip(weights[c],qs))
        shared = n*n*math.fsum(float(w)*(float(q)-mean_q)**2 for w,q in zip(weights[c],qs))
        scalar[c] = (n*mean_q,event,shared)
    archive = np.load(h.ROOT/"mortality/replicate_archive.npz",allow_pickle=False)
    indices,cell_deaths = archive["hazard_indices"],archive["cell_deaths"]
    require(indices.shape==cell_deaths.shape==(cfg["replicates"],64),"Predictive replicate census incomplete")
    require(np.all(indices<256) and np.all(cell_deaths<=archive["cell_n"]),"Impossible shared hazard or event count")
    pmf_table = h.read_csv(h.ROOT/"mortality/predictive_pmf.csv")
    summaries = h.read_csv(h.ROOT/"mortality/predictive_summary.csv")
    check_rows, checked_pmfs = [], []
    for r in summaries:
        name = r["group"]
        cells = [c for c in counts if name=="ALL" or (c<16)==(name=="HEMLOCK")]
        full = np.ones(1)
        # Positive real direct convolution is independent of loggamma/FFT route.
        for c in sorted(cells,key=lambda c:counts[c]): full = np.convolve(full,cell_pmfs[c])
        original = np.array([float(v["probability"]) for v in pmf_table if v["group"]==name])
        require(len(original)==len(full),"Independent PMF support differs")
        error = float(np.max(abs(original-full)))
        require(error<2e-12 and abs(math.fsum(full)-1)<1e-12,"Independent recurrence/direct PMF disagrees")
        for k,p in enumerate(full): checked_pmfs.append(dict(group=name,deaths=k,probability=float(p)))
        mu,event,hazard = [math.fsum(scalar[c][i] for c in cells) for i in range(3)]
        obs = sum(observed[c] for c in cells)
        for field,actual in [("predictive_mean",mu),("conditional_event_variance",event),("shared_hazard_variance",hazard)]:
            require(abs(float(r[field])-actual)<1e-8,"Scalar variance/mean decomposition disagrees")
        cdf = np.cumsum(full)
        quantiles = [int(np.searchsorted(cdf,p)) for p in [.5,.025,.975]]
        require(quantiles==[int(r[f]) for f in ["predictive_median","predictive_lower95","predictive_upper95"]],"Predictive quantiles disagree")
        tail = math.fsum(float(v) for v in full[obs:])
        require(abs(tail-float(r["upper_tail"]))<1e-11,"Independent upper-tail probability disagrees")
        sample = cell_deaths[:,cells].sum(axis=1).astype(int)
        empirical = np.bincount(sample,minlength=len(full))/len(sample)
        cdf_error = float(np.max(abs(np.cumsum(empirical)-cdf)))
        require(cdf_error<.006 and abs(float(sample.mean())-mu)<6*math.sqrt((event+hazard)/len(sample)),"Individual sampler disagrees with exact law")
        # Independent RNG and binomial aggregation sample the SAME hierarchical law.
        generator = h.rng("binomial_check", {"HEMLOCK":0,"NONHEMLOCK":1,"ALL":2}[name])
        alternate = np.zeros(cfg["replicates"],dtype=np.int64)
        for c in sorted(cells):
            draws = generator.choice(256,size=cfg["replicates"],p=weights[c])
            alternate += generator.binomial(counts[c],qs[draws])
        alt_hist = np.bincount(alternate,minlength=len(full))/len(alternate)
        alt_error = float(np.max(abs(np.cumsum(alt_hist)-cdf)))
        require(alt_error<.006 and abs(float(alternate.mean())-mu)<6*math.sqrt((event+hazard)/len(alternate)),"Independent binomial sampler disagrees")
        require(abs(float(sample.var())/(event+hazard)-1)<.03 and abs(float(alternate.var())/(event+hazard)-1)<.03,"Predictive sampler variance disagrees")
        check_rows.append(dict(group=name, scalar_mean=mu, scalar_event_variance=event, scalar_shared_hazard_variance=hazard,
                               recurrence_direct_tail=tail, FFT_vs_direct_max_probability_error=error,
                               primary_CDF_max_error=cdf_error, independent_sampler_CDF_max_error=alt_error,
                               independent_sampler_mean=float(alternate.mean()), independent_sampler_variance=float(alternate.var()),
                               independent_sampler_tail=float(np.mean(alternate>=obs)), status="PASS"))
        print("Independent predictive law",name,"PASS",flush=True)
    # Replay four registered chunks using actual independent per-stem outcomes.
    digests = {int(r["chunk"]):r["ordered_individual_event_bits_sha256"] for r in h.read_csv(h.ROOT/"mortality/individual_event_digests.csv")}
    for chunk in [0,1,255,511]:
        hazards,events = h.rng("hazards",chunk),h.rng("individual_events",chunk)
        digest = hashlib.sha256()
        start = chunk*cfg["chunk_size"]
        for c,n in sorted(counts.items()):
            draws = hazards.choice(256,p=weights[c],size=cfg["chunk_size"])
            require(np.array_equal(draws,indices[start:start+cfg["chunk_size"],c]),"Shared hazard index replay differs")
            # expm1 evaluates small q accurately; tested independently above against exp.
            q = -np.expm1(-grid[draws]*3)
            outcomes = events.random((cfg["chunk_size"],n)) < q[:,None]
            digest.update(np.packbits(outcomes,axis=None).tobytes())
            require(np.array_equal(outcomes.sum(axis=1),cell_deaths[start:start+cfg["chunk_size"],c]),"Individual event count replay differs")
        require(digest.hexdigest()==digests[chunk],"Individual event bitstream replay differs")
    h.write_csv("checks/predictive_independent.csv",check_rows)
    h.write_csv("checks/predictive_pmf_direct.csv",checked_pmfs)
    h.write_json("checks/predictive.json",dict(status="PASS", fixture=fixtures(), numerical_route="mode-recursion plus positive direct convolution and scalar fsum",
                 complete_predictive_support_checked=True, groups=3, independent_binomial_replicates_each=cfg["replicates"],
                 individual_bitstream_replayed_chunks=[0,1,255,511], mean_and_variance_checks=True, no_refit=True))

if __name__ == "__main__":
    import sys
    h.guard()
    {"attribution":attribution,"predictive":predictive}[sys.argv[1]]()
