"""Reconstruct compact evidence and actively reject semantic corruption."""
import copy
import math
from collections import Counter
import numpy as np
import h06_common as h

def require(condition, message):
    if not condition: raise ValueError(message)

def packet():
    direct = h.read_csv(h.ROOT/"checks/predictive_pmf_direct.csv")
    result = dict(partition=h.read_csv(h.ROOT/"attribution/partition.csv"),
                  overlaps=h.read_csv(h.ROOT/"attribution/overlapping_causes.csv"),
                  cells=h.read_csv(h.ROOT/"mortality/cell_census.csv"),
                  summary=h.read_csv(h.ROOT/"mortality/final_predictive_summary.csv"), direct={})
    for name in ["HEMLOCK","NONHEMLOCK","ALL"]:
        pmf = np.array([float(r["probability"]) for r in direct if r["group"]==name])
        cumulative = np.cumsum(pmf)
        observed = next(int(r["observed_deaths"]) for r in result["summary"] if r["group"]==name)
        result["direct"][name] = dict(tail=math.fsum(pmf[observed:]), quantiles=[int(np.searchsorted(cumulative,q)) for q in [.5,.025,.975]],
                                      normalization=math.fsum(pmf))
    return result

def validate(data):
    expected_categories = {"MODEL_UNLAWFULNESS","JUVENILE_ONLY","BASAL_AREA_ONLY","JOINT_JUVENILE_BASAL_AREA","SUCCESS"}
    for name,expected in h.config()["terminal"]["expected_viable"].items():
        chosen = [r for r in data["partition"] if r["reconstruction"]==name]
        require(len(chosen)==5 and {r["outcome"] for r in chosen}==expected_categories,"Incomplete or duplicated outcome partition")
        require(all(int(r["candidate_paths"])==4096 and int(r["paths"])>=0 for r in chosen),"Wrong unconditional candidate denominator")
        require(sum(int(r["paths"]) for r in chosen)==4096,"Failure reconciliation incorrect")
        counts = {r["outcome"]:int(r["paths"]) for r in chosen}
        require(counts["SUCCESS"]==expected,"Frozen viable count changed")
        overlap = {r["cause"]:int(r["paths"]) for r in data["overlaps"] if r["reconstruction"]==name}
        require(overlap["ANY_JUVENILE_FAILURE"]==counts["JUVENILE_ONLY"]+counts["JOINT_JUVENILE_BASAL_AREA"],"Juvenile overlap lost")
        require(overlap["ANY_BASAL_AREA_FAILURE"]==counts["BASAL_AREA_ONLY"]+counts["JOINT_JUVENILE_BASAL_AREA"],"Basal overlap lost")
        require(overlap["JUVENILE_AND_BASAL"]==counts["JOINT_JUVENILE_BASAL_AREA"],"Joint overlap lost")
    model = h.model()
    probability = -np.expm1(-3*model["hazard_grid"])
    require(len(data["summary"])==3 and {r["group"] for r in data["summary"]}=={"HEMLOCK","NONHEMLOCK","ALL"},"Missing predictive group")
    for r in data["summary"]:
        name = r["group"]
        cells = [c for c in data["cells"] if name=="ALL" or c["group"]==name]
        require(int(r["adult_stems"])==sum(int(c["n"]) for c in cells),"Adult support altered")
        require(int(r["observed_deaths"])==sum(int(c["observed_deaths"]) for c in cells),"Observed mortality altered")
        mean=event=hazard=0.
        for c in cells:
            n,cell = int(c["n"]),int(c["cell"])
            p = model["hazard_posterior"][cell]
            eq = math.fsum(float(w)*float(q) for w,q in zip(p,probability))
            mean += n*eq
            event += n*math.fsum(float(w)*float(q)*(1-float(q)) for w,q in zip(p,probability))
            hazard += n*n*math.fsum(float(w)*(float(q)-eq)**2 for w,q in zip(p,probability))
        for key,val in [("predictive_mean",mean),("conditional_event_variance",event),("shared_hazard_variance",hazard),("total_predictive_variance",event+hazard)]:
            require(abs(float(r[key])-val)<1e-8,"Predictive moment/variance law changed: "+key)
        require(abs(float(r["shared_hazard_variance_fraction"])-hazard/(event+hazard))<1e-12,"Shared uncertainty fraction changed")
        exact=data["direct"][name]
        require(abs(exact["normalization"]-1)<1e-12,"Positive predictive PMF not normalized")
        require(abs(float(r["upper_tail"])/exact["tail"]-1)<1e-12,"Reported rare tail not from positive convolution")
        require([int(r[k]) for k in ["predictive_median","predictive_lower95","predictive_upper95"]]==exact["quantiles"],"Wrong central predictive quantiles")
        require(int(r["MC_replicates"])==262144 and int(r["MC_tail_hits"])>=0,"Predictive replicate denominator changed")
        require(float(r["MC_tail"])==int(r["MC_tail_hits"])/262144,"Wrong sampled tail frequency")
        require(float(r["MC_tail_Wilson_upper"])>0 if int(r["MC_tail_hits"])==0 else True,"Zero sampled hits misrepresented as certainty")
    sums = {r["group"]:r for r in data["summary"]}
    require(int(sums["ALL"]["observed_deaths"])==int(sums["HEMLOCK"]["observed_deaths"])+int(sums["NONHEMLOCK"]["observed_deaths"]),"Group deaths do not reconcile")

def audit():
    h.guard()
    data=packet()
    validate(data)
    mutants=[]
    tests=[("missing_joint_category",lambda d:d["partition"].pop(3)),
           ("conditional_candidate_denominator",lambda d:d["partition"][0].update(candidate_paths="43")),
           ("lost_failure",lambda d:d["partition"][1].update(paths="4052")),
           ("lost_overlap",lambda d:d["overlaps"][1].update(paths="0")),
           ("independent_marginal_hazard_variance",lambda d:d["summary"][0].update(shared_hazard_variance="0")),
           ("FFT_roundoff_tail",lambda d:d["summary"][0].update(upper_tail="1.1192966670372426e-15")),
           ("zero_hit_zero_probability",lambda d:d["summary"][0].update(upper_tail="0")),
           ("omitted_observed_death",lambda d:d["summary"][0].update(observed_deaths="684")),
           ("wrong_predictive_interval",lambda d:d["summary"][0].update(predictive_upper95="685")),
           ("wrong_MC_denominator",lambda d:d["summary"][0].update(MC_replicates="65536")),
           ("zero_hit_no_uncertainty",lambda d:d["summary"][0].update(MC_tail_Wilson_upper="0")),
           ("missing_nonhemlock",lambda d:d["summary"].pop(1))]
    for name,mutate in tests:
        corrupted=copy.deepcopy(data)
        mutate(corrupted)
        try: validate(corrupted)
        except ValueError as error: mutants.append(dict(control=name,status="REJECTED",reason=str(error)))
        else: raise ValueError("Hostile control accepted: "+name)
    archive=np.load(h.ROOT/"mortality/replicate_archive.npz",allow_pickle=False)
    deaths=archive["cell_deaths"]
    replicates=h.read_csv(h.ROOT/"mortality/predictive_replicates.csv")
    require(len(replicates)==262144,"Lost individual replicates")
    for i,r in enumerate(replicates):
        require(int(r["replicate"])==i and int(r["hemlock"])==int(deaths[i,:16].sum()) and int(r["nonhemlock"])==int(deaths[i,16:].sum()),"Individual replicate summary altered")
        require(int(r["all_adults"])==int(r["hemlock"])+int(r["nonhemlock"]),"Individual group sums changed")
    h.write_json("checks/hostile_controls.json",dict(status="PASS", active_controls=mutants, accepted_corruption=0, all262144_group_sums_checked=True))
    print("Compact reconstruction and12 active hostile controls PASS",flush=True)

if __name__ == "__main__":
    audit()
