"""Replay the exact first HFD05 stream and keep overlapping predicate failures."""
import math
import numpy as np
import h06_common as h
from h05_exposure import records
from h05_exposure_coupled import reconstruct
from h05_growth import simulate

CATEGORIES = ["MODEL_UNLAWFULNESS", "JUVENILE_ONLY", "BASAL_AREA_ONLY", "JOINT_JUVENILE_BASAL_AREA", "SUCCESS"]

class Capture:
    def __init__(self, generator):
        self.generator = generator
        self.calls = {}
    def __getattr__(self, name):
        method = getattr(self.generator, name)
        if name not in ["binomial", "gamma", "poisson"]:
            return method
        def run(*args, **kwargs):
            result = method(*args, **kwargs)
            self.calls[name] = result.copy()
            return result
        return run

def categories(lawful, basal, juvenile):
    return np.where(~lawful, CATEGORIES[0], np.where(~basal & ~juvenile, CATEGORIES[3],
           np.where(~juvenile, CATEGORIES[1], np.where(~basal, CATEGORIES[2], CATEGORIES[4]))))

def run():
    h.guard()
    cfg = h.config()["terminal"]
    r = h.read_csv(h.FIVE / "corrected_pipeline/local_extreme_selection.csv")[cfg["selection_index"]]
    for k in ["state", "parameter_rep", "variant_index", "horizon", "tau"]:
        if int(r[k]) != cfg[k]: raise ValueError("Selected identity mismatch")
    for k in ["scenario", "query"]:
        if r[k] != cfg[k]: raise ValueError("Selected identity mismatch")
    for k in ["alpha", "beta", "rho", "theta"]:
        if float(r[k]) != cfg[k]: raise ValueError("Selected semantic mismatch")
    _, e0, e1, core, regenerated = h.prior.original.load_data()
    model = h.model()
    for k in model:
        if not np.array_equal(model[k], regenerated[k]): raise ValueError("Frozen fit reconstruction differs")
    _, groups = records(e0, e1, core)
    seeds = h.prior.original.frozen.Seeds()
    a, b, ledger = reconstruct(core, model, groups, seeds.rng("reconstruct", cfg["state"]),
                               h.prior.rng("individual_state", 0, cfg["state"]), "B0")
    params = h.prior.original.demography.parameters(model, seeds.rng("parameters", cfg["state"], cfg["parameter_rep"]))
    design = h.prior.original.frozen.load("config/design.json")
    scenario = next(s for s in design["scenarios"] if s["id"] == cfg["scenario"])
    rows, summaries, overlaps, origins, cohorts = [], [], [], [], []
    archives = {}
    for identity, state in [("POOLED", a), ("INDIVIDUAL", b)]:
        banks, laws, terminal_n, terminal_d = [], [], [], []
        for chunk in range(cfg["chunks"]):
            generator = Capture(h.prior.rng("corrected_future", 100, cfg["selection_index"], chunk))
            bank, lawful = simulate(state, params, scenario, 1, 256, generator, "gamma_growth", False)
            banks.append(bank); laws.append(lawful)
            terminal_n.append(np.concatenate([generator.calls["binomial"], generator.calls["poisson"]], axis=1))
            terminal_d.append(np.concatenate([np.clip(state["d"] + generator.calls["gamma"], 1, 300), np.full((256, 16), 1.7)], axis=1))
        bank, lawful = np.concatenate(banks), np.concatenate(laws)
        n, d = np.concatenate(terminal_n), np.concatenate(terminal_d)
        response = h.prior.original.demography.response(bank, lawful, model, 1, cfg["alpha"], cfg["beta"], cfg["tau"], cfg["rho"])
        if not response["present"]: raise ValueError("Origin no longer realizes organization")
        basal = response["structure"] >= cfg["alpha"]
        juvenile = response["regeneration"] >= cfg["beta"]
        label = categories(lawful, basal, juvenile)
        success = label == "SUCCESS"
        if not np.array_equal(success, response["viable"]): raise ValueError("Partition misses excursion failure")
        if int(success.sum()) != cfg["expected_viable"][identity]: raise ValueError("Exact predecessor viable reconciliation failed")
        if not np.allclose(np.sum(n*d*d, axis=1)*math.pi/40000, bank[:,1,1], rtol=0, atol=1e-9): raise ValueError("Captured cohort BA mismatch")
        if not np.array_equal(np.sum(n*(d<10), axis=1), bank[:,1,2]): raise ValueError("Captured cohort juvenile mismatch")
        for i in range(cfg["K"]):
            rows.append(dict(reconstruction=identity, path=i, chunk=i//256, path_in_chunk=i%256,
                             model_lawful=bool(lawful[i]), terminal_juvenile_failure=bool(~juvenile[i]),
                             terminal_basal_area_failure=bool(~basal[i]), longest_excursion=int(response["longest_excursion"][i]),
                             viable=bool(success[i]), outcome=label[i], terminal_N=int(bank[i,1,0]),
                             terminal_BA=float(bank[i,1,1]), terminal_J=int(bank[i,1,2]),
                             basal_area_ratio=float(response["structure"][i]), juvenile_ratio=float(response["regeneration"][i])))
        for name in CATEGORIES:
            summaries.append(dict(reconstruction=identity, outcome=name, paths=int((label==name).sum()), candidate_paths=cfg["K"]))
        for name, mask in [("ANY_UNLAWFUL", ~lawful), ("ANY_JUVENILE_FAILURE", ~juvenile), ("ANY_BASAL_AREA_FAILURE", ~basal),
                           ("JUVENILE_AND_BASAL", ~juvenile & ~basal), ("UNLAWFUL_AND_JUVENILE", ~lawful & ~juvenile),
                           ("UNLAWFUL_AND_BASAL", ~lawful & ~basal), ("ALL_THREE", ~lawful & ~juvenile & ~basal)]:
            overlaps.append(dict(reconstruction=identity, cause=name, paths=int(mask.sum()), mutually_exclusive=False))
        origins.append(dict(reconstruction=identity, N=int(state["n"].sum()), BA=float(bank[0,0,1]), J=int(bank[0,0,2]),
                            BA_threshold=cfg["alpha"]*model["baseline_ba"], J_threshold=cfg["beta"]*model["baseline_juveniles"], present=True))
        for c in range(64):
            crossings = np.sum(n[:,c]*(state["d"][c]<10)*(d[:,c]>=10), dtype=np.int64)
            cohorts.append(dict(reconstruction=identity, cell=c, initial_n=int(state["n"][c]), initial_d=float(state["d"][c]),
                                initial_juvenile=bool(state["d"][c]<10), terminal_mean_n=float(n[:,c].mean()),
                                terminal_mean_J=float((n[:,c]*(d[:,c]<10)).mean()), mean_crossing10_stems=float(crossings/4096)))
        archives.update({identity+"_"+key: value for key, value in [("bank",bank),("lawful",lawful),("terminal_n",n),("terminal_d",d),("initial_n",state["n"]),("initial_d",state["d"])]})
        print(identity, "attributed", int(success.sum()), "of", cfg["K"], flush=True)
    np.savez_compressed(h.output("attribution/path_archive.npz"), **archives)
    h.write_csv("attribution/paths.csv", rows)
    h.write_csv("attribution/partition.csv", summaries)
    h.write_csv("attribution/overlapping_causes.csv", overlaps)
    h.write_csv("attribution/origins.csv", origins)
    h.write_csv("attribution/cohort_threshold_crossings.csv", cohorts)
    h.write_json("attribution/replay.json", dict(status="PASS", selected=r, reconstruction_ledger=ledger,
                 candidate_paths_per_reconstruction=4096, exact_viable=cfg["expected_viable"], all_candidates_retained=True,
                 historical_stream="HFD05 corrected_future(100,0,chunk)", growth_law="ORIGINAL_GAMMA", terminal_only_partition_earned=True,
                 baseline_hazards=params["hazard"].tolist(), growth=params["growth"].tolist(), recruitment=params["recruit"].tolist()))

if __name__ == "__main__":
    run()
