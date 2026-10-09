"""Frozen design, typed provenance and explicitly registered random streams."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HFD01 = ROOT.parent / "HFD01"
PROVENANCE = {"OBSERVED", "DERIVED", "RECOVERED_FROM_PRIOR_OBSERVATION", "IMPUTED",
              "SIMULATED", "LITERATURE_LICENSED", "EXTERNALLY_CONSTRAINED"}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(name):
    return json.loads((ROOT / name).read_text())


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def quantity(value, provenance_class, source_id, method_id, namespace="latent_origin"):
    if provenance_class not in PROVENANCE or not source_id or not method_id:
        raise ValueError("Missing/unknown provenance")
    if namespace.startswith("raw_") and provenance_class != "OBSERVED":
        raise ValueError("Synthetic/imputed datum cannot enter observed namespace")
    return dict(value=value, provenance_class=provenance_class, source_id=source_id,
                method_id=method_id, namespace=namespace)


class Seeds:
    def __init__(self):
        self.registry = load("provenance/random_seed_registry.json")
        self.used = {}

    def seed(self, namespace, *indices):
        shape = self.registry["namespaces"].get(namespace)
        if shape is None or len(shape) != len(indices) or any(not 0 <= i < n for i, n in zip(indices, shape)):
            raise ValueError("Unregistered random stream")
        payload = [self.registry["registry_version"], self.registry["master_seed"], namespace, list(indices)]
        seed = int.from_bytes(hashlib.sha256(canonical(payload)).digest()[:8], "big")
        self.used[namespace + ":" + ",".join(map(str, indices))] = seed
        return seed

    def rng(self, namespace, *indices):
        import numpy as np
        return np.random.Generator(np.random.PCG64(self.seed(namespace, *indices)))


def recover_coordinate(old, new, name):
    if new.get(name) not in {None, "", "NA"}:
        return quantity(new[name], "OBSERVED", "HF253-E1", "recorded_field", "raw_E1")
    stable = old is not None and new.get("stem.id") == old.get("stem.id") and new.get("tree.id") == old.get("tree.id")
    if not stable or old.get(name) in {None, "", "NA"}:
        return None
    return quantity(old[name], "RECOVERED_FROM_PRIOR_OBSERVATION", "HF253-E0",
                    "stable_tagged_stationary_root_location", "recovered_E1")


def fate(row):
    if row is None or row.get("exact.date") in {None, "", "NA"}:
        return None
    return {"alive": 1, "stem dead": 0}.get(row.get("df.status"))


def verify_frozen_inputs():
    config = load("config/design.json")
    if sha(HFD01 / "manifests/input_manifest.json") != config["hfd01_input_sha256"]:
        raise ValueError("HFD01 input manifest changed")
    manifest = json.loads((HFD01 / "manifests/artifact_manifest.json").read_text())
    sources = {a["source_id"]: a for a in manifest["artifacts"]}
    if config["training_artifacts"] != ["HF253-E0", "HF253-E1"]:
        raise ValueError("E2/context source entered training")
    for key in config["training_artifacts"]:
        if sha(HFD01 / sources[key]["local_cache_path"]) != sources[key]["sha256"]:
            raise ValueError("Frozen census bytes changed")
    return sources
