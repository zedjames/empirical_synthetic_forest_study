"""Additive HFD06 outputs and frozen input provenance; no predecessor writes."""
import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[3]
FIVE = ROOT.parent / "HFD05"
TWO = ROOT.parent / "HFD02S"
sys.path.insert(0, str(FIVE / "scripts"))
import h05_common as prior

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def config():
    return json.loads((ROOT / "config/protocol.json").read_text())

def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))

def output(name):
    path = (ROOT / name).resolve()
    if ROOT.resolve() not in path.parents:
        raise ValueError("Output outside HFD06")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path

def write_json(name, value):
    output(name).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")

def write_csv(name, rows):
    rows = list(rows)
    if not rows:
        raise ValueError("Undeclared empty table")
    with output(name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

def rng(domain, *indices):
    cfg = config()
    if domain not in cfg["random_domains"] or any(type(i) is not int or i < 0 for i in indices):
        raise ValueError("Undeclared RNG identity")
    data = json.dumps([cfg["master_seed"], "HFD06", domain, list(indices)], separators=(",", ":")).encode()
    return np.random.Generator(np.random.PCG64(int.from_bytes(hashlib.sha256(data).digest()[:8], "big")))

def protected():
    prefixes = ["research/health_formally_defined/harvard_forest/" + name + "/"
                for name in ["HFD01", "HFD02S", "HFD03", "HFD04", "HFD05"]]
    prefixes += ["scripts/verify_hfd" + name + "_harvard_forest.py"
                 for name in ["01", "02s", "03", "04", "05"]]
    records = subprocess.check_output(["git", "ls-tree", "-r", config()["baseline"], "--", *prefixes], cwd=REPO, text=True)
    # Hash actual bytes in one batch, including every predecessor artifact.
    names, expected = [], []
    for line in records.splitlines():
        meta, name = line.split("\t", 1)
        names.append(name)
        expected.append(meta.split()[2])
    actual = subprocess.check_output(["git", "hash-object", "--stdin-paths"], input="\n".join(names)+"\n", cwd=REPO, text=True).splitlines()
    if actual != expected:
        raise ValueError("Frozen predecessor bytes changed")
    return dict(zip(names, expected))

def guard():
    record = json.loads((ROOT / "provenance/registration.json").read_text())
    for name, digest in record["sources"].items():
        if sha(ROOT / name) != digest:
            raise ValueError("Frozen HFD06 science changed: " + name)
    for name, digest in record["inputs"].items():
        if sha(ROOT.parent / name) != digest:
            raise ValueError("Frozen input changed: " + name)
    if np.__version__ != "2.0.2" or sys.version_info[:3] != (3, 9, 6):
        raise ValueError("Runtime mismatch")

def model():
    data = json.loads((TWO / "models/fitted_parameters.json").read_text())
    return {k: np.asarray(v) if isinstance(v, list) else v for k, v in data.items()}

def register():
    if (ROOT / "provenance/registration.json").exists() or (ROOT / "mortality").exists() or (ROOT / "attribution").exists():
        raise ValueError("Registration must precede outcomes")
    names = ["config/protocol.json"] + sorted(str(p.relative_to(ROOT)) for p in (ROOT / "scripts").glob("*.py"))
    inputs = ["HFD02S/models/fitted_parameters.json", "HFD05/corrected_pipeline/local_extreme_selection.csv",
              "HFD05/corrected_pipeline/local_precision.csv", "HFD05/external_sources/HF453/mortality_predictions.csv",
              "HFD05/external_sources/HF453/tag_crosswalk.csv", "HFD05/external_sources/HF453/mortality_scores.csv"]
    write_json("provenance/registration.json", dict(status="REGISTERED_BEFORE_HFD06_OUTCOMES", sources={n: sha(ROOT/n) for n in names},
               inputs={n: sha(ROOT.parent/n) for n in inputs}, protected_git_blobs=protected(), runtime=config()["runtime"],
               historical_results_known=True, selected_targets_fixed_by_user=True, public_release=False))

if __name__ == "__main__":
    register()
