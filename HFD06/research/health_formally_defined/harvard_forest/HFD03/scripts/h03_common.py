"""HFD03 isolation, provenance, intervals and frozen source interfaces."""
import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from statistics import NormalDist
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[3]
OLD = ROOT.parent / "HFD02S"
for folder in ("scripts", "models", "validation"):
    sys.path.append(str(OLD / folder))
import contracts as frozen
import empirical
import demography as demographic

PROTOCOL = json.loads((ROOT / "config/protocol.json").read_text())
PROTOCOL_HASH = hashlib.sha256((ROOT / "config/protocol.json").read_bytes()).hexdigest()
RUN = ROOT / ".runs" / PROTOCOL_HASH[:16]
PROTECTED = ["research/health_formally_defined/harvard_forest/HFD01/",
             "research/health_formally_defined/harvard_forest/HFD02S/",
             "scripts/verify_hfd01_harvard_forest.py", "scripts/verify_hfd02s_harvard_forest.py"]


def rng(namespace, *indices):
    registered = {"mask", "calibration", "state", "future", "diagnostic", "world",
                  "observation", "truth_future", "estimate_state", "estimate_params",
                  "estimate_future", "inner_future", "outer_future"}
    if namespace not in registered or any(not isinstance(i, int) or i < 0 for i in indices):
        raise ValueError("Unregistered random namespace")
    payload = frozen.canonical([PROTOCOL["master_seed"], namespace, list(indices)])
    seed = int.from_bytes(hashlib.sha256(payload).digest()[:8], "big")
    return np.random.Generator(np.random.PCG64(seed))


def write_json(name, data):
    path = ROOT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+"\n")


def write_csv(name, rows, fields=None):
    rows = list(rows)
    path = ROOT / name
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def frozen_guard():
    baseline = PROTOCOL["planning_scientific_baseline"]
    entries = subprocess.check_output(["git", "ls-tree", "-r", baseline, "--", *PROTECTED],
                                     cwd=REPO, text=True).splitlines()
    identities = {}
    for entry in entries:
        meta, name = entry.split("\t", 1)
        expected = meta.split()[2]
        actual = subprocess.check_output(["git", "hash-object", "--", name],
                                         cwd=REPO, text=True).strip()
        if actual != expected:
            raise ValueError("Frozen tracked source changed: " + name)
        identities[name] = expected
    return identities


def inputs():
    sources = frozen.verify_frozen_inputs()
    tables = []
    for key in ("HF253-E0", "HF253-E1"):
        tables.append({r["stem.id"]: r for r in read_csv(frozen.HFD01 / sources[key]["local_cache_path"])})
    identity = json.loads((ROOT / "provenance/identity_inputs.json").read_text())
    for item in identity["sources"]:
        if frozen.sha(ROOT / ".inputs" / item["name"]) != item["sha256"]:
            raise ValueError("Identity source changed")
    trees = {r["tree.id"]: r for r in read_csv(ROOT / ".inputs/hf253-04-trees-2014.csv")}
    return sources, tables[0], tables[1], trees


def wilson(success, total, confidence=0.95):
    if total <= 0:
        raise ValueError("Empty MC bank")
    x = np.asarray(success, dtype=float)
    if np.any((x < 0) | (x > total)):
        raise ValueError("Invalid success census")
    z = NormalDist().inv_cdf((1+confidence)/2)
    p = x/total
    center = (p+z*z/(2*total))/(1+z*z/total)
    half = z*np.sqrt(p*(1-p)/total+z*z/(4*total*total))/(1+z*z/total)
    lower=np.maximum(0,center-half)
    upper=np.minimum(1,center+half)
    # Wilson's exact boundary identities: floating roundoff must not put
    # U(K,K) below1 (or L(0,K) above0) and force a spurious threshold Boolean.
    return np.where(x==0,0,lower),np.where(x==total,1,upper)


def classify(present, counts, total, theta):
    counts = np.asarray(counts)
    confidence = 0.95 if len(counts) == 1 else 0.975
    lower, upper = wilson(counts, total, confidence)
    if not present or np.any(upper < theta):
        return "FALSE", lower, upper
    if np.all(lower > theta):
        return "TRUE", lower, upper
    return "MC_UNRESOLVED", lower, upper


def phase(bank, lawful, model, horizon, alpha, beta, tau, rho, theta, query):
    response = demographic.response(bank, lawful, model, horizon, alpha, beta, tau, rho)
    counts = [int(response["viable"].sum())]
    if query == "Q2":
        counts.append(int(response["reserve_viable"].sum()))
    status, lo, hi = classify(response["present"], counts, len(bank), theta)
    point = response["present"] and all(c/len(bank) >= theta for c in counts)
    return dict(present=response["present"], viable_mass=counts[0]/len(bank),
                reserve_mass=float(response["reserve_viable"].mean()),
                point_health=bool(point), status=status,
                lower=float(min(lo)), upper=float(min(hi)))


class TraceRNG:
    """Captures frozen draws without giving the reconstruction any truth."""
    def __init__(self, generator):
        self.generator = generator
        self.trace = {}
    def __getattr__(self, method):
        target = getattr(self.generator, method)
        def call(*args, **kwargs):
            result = target(*args, **kwargs)
            self.trace.setdefault(method, []).append(np.array(result, copy=True))
            return result
        return call


class MechanismRNG:
    """Explicit alternative kernel: replace unknown-cohort hazard shifts.

    The old reconstruct function and every other distribution stay unchanged.
    A base variate is always consumed, including deterministic B1, for coupling.
    B1's cell shifts enter the second binomial; 31 association-branch priors
    remain the original separate unresolved-identity model.
    """
    def __init__(self, generator, family):
        self.generator = generator
        self.family = family
        self.normal_index = 0
        self.binomial_index = 0
    def __getattr__(self, method):
        if method == "binomial":
            def draw(n,p,*args,**kwargs):
                index=self.binomial_index
                self.binomial_index+=1
                if self.family=="B1" and index==1:
                    family=PROTOCOL["missingness"]["B1"]
                    c=np.arange(64)
                    shift=(np.array(family["taxon_log_shift"])[c//16]
                           +np.array(family["size_log_shift"])[c%4]
                           +np.array(family["stratum_log_shift"])[(c//4)%4])
                    p=np.exp(np.log(np.clip(p,1e-300,1))*np.exp(shift))
                return self.generator.binomial(n,p,*args,**kwargs)
            return draw
        if method != "normal":
            return getattr(self.generator, method)
        def call(*args, **kwargs):
            result = self.generator.normal(*args, **kwargs)
            index = self.normal_index
            self.normal_index += 1
            # normal #0 is fitted growth; #1 and #2 are MNAR and wet shifts.
            family = PROTOCOL["missingness"][self.family]
            if index == 1:
                if self.family == "B1":
                    return 0.0
                return result / 0.45 * family["mnar_sd"]
            if index == 2:
                return result / 0.7 * family["wet_sd"] + family.get("wet_log_shift_mean", 0)
            return result
        return call


def reconstruct(core, model, state_index, family="B0"):
    generator = rng("state", state_index)
    if family != "B0":
        generator = MechanismRNG(generator, family)
    trace=TraceRNG(generator)
    state = demographic.reconstruct(core, model, trace)
    # All origin fields are inferred; observed support is ANCESTRY, not direct 2020 truth.
    components = state["components"]
    state["provenance"] = dict(observed_anchored_projected_ba=components["measured_projected_ba"],
                               imputed_ba=components["imputed_ba"],
                               current_origin_class="IMPUTED",
                               prior_coordinate_class="RECOVERED_FROM_PRIOR_OBSERVATION")
    measured_count=trace.trace["binomial"][0]
    juvenile=state["d"]<10
    state["provenance"].update(observed_anchored_projected_count=int(measured_count.sum()),
        imputed_count=int(state["n"].sum()-measured_count.sum()),
        observed_anchored_projected_juveniles=int((measured_count*juvenile).sum()),
        imputed_juveniles=int(((state["n"]-measured_count)*juvenile).sum()),
        juvenile_allocation="Source counts apportioned within the final cohort-RMS size classification",
        recovered_current_origin_ba=0)
    return state
