"""Exact package-only runner; input path adaptation, no scientific rewrites."""
import argparse
import importlib
import json
import subprocess
import sys
from pathlib import Path
import hashlib

ROOT=Path(__file__).resolve().parent
STUDY=ROOT/"research/health_formally_defined/harvard_forest"
SIX=STUDY/"HFD06"

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def audit():
    manifest=json.loads((ROOT/"HFD06_Artifact_Manifest.json").read_text())
    for name,record in manifest["files"].items():
        path=ROOT/name
        if path.is_symlink() or not path.is_file() or sha(path)!=record["sha256"]:
            raise ValueError("Package bytes changed: "+name)
    allowed=set(manifest["files"])|{"HFD06_Artifact_Manifest.json"}
    allowed.update(r["target"] for r in json.loads((ROOT/"Public_Inputs.json").read_text())["sources"])
    for path in ROOT.rglob("*"):
        if path.is_file() and "__pycache__" not in path.parts and str(path.relative_to(ROOT)) not in allowed:
            raise ValueError("Unlisted package file")
    return manifest

def retrieve():
    audit()
    for record in json.loads((ROOT/"Public_Inputs.json").read_text())["sources"]:
        path=ROOT/record["target"]
        if not path.resolve().is_relative_to(ROOT.resolve()): raise ValueError("Input traversal")
        if not path.exists():
            path.parent.mkdir(parents=True,exist_ok=True)
            subprocess.run(["curl","--fail","--location","--silent","--show-error",record["url"],"--output",str(path)],check=True)
        if path.is_symlink() or sha(path)!=record["sha256"]: raise ValueError("Publisher input bytes changed")
    print(json.dumps(dict(status="PASS", narrow_public_inputs=4, public_release=False)))

def modules():
    audit()
    sys.path.insert(0,str(SIX/"scripts"))
    import h06_common as h
    inputs=json.loads((ROOT/"Public_Inputs.json").read_text())["sources"]
    def sources():
        records={}
        for r in inputs:
            path=ROOT/r["target"]
            if sha(path)!=r["sha256"]: raise ValueError("Package input hash mismatch")
            if r["id"] in ["HF253-E0","HF253-E1"]:
                records[r["id"]]=dict(r,local_cache_path=str(path))
        return records
    # Only canonical local input locations are adapted. No private Git guard
    # is called or pretended to be reconstructed by this package.
    h.prior.original.frozen.verify_frozen_inputs=sources
    required=["h06_common","h05_common","h04_common","h03_common","contracts","empirical","demography","estimator"]
    for name in required:
        module=importlib.import_module(name)
        if not Path(module.__file__).resolve().is_relative_to(ROOT.resolve()): raise ValueError("Private module import")
    return h

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("stage",choices=["audit","retrieve","checks","attribution","predictive"])
    stage=parser.parse_args().stage
    if stage=="retrieve": retrieve(); return
    if stage=="audit":
        print(json.dumps(dict(status="PASS",included_files=len(audit()["files"]),private_history_verified=False))); return
    h=modules()
    if stage=="attribution":
        import h06_attribution as a
        import h06_independent as i
        a.run();i.attribution()
    elif stage=="predictive":
        import h06_mortality as m
        import h06_independent as i
        import h06_precision as p
        m.run();i.predictive();p.run()
    else:
        import h06_checks as c
        c.audit()
    for name,module in list(sys.modules.items()):
        if name.startswith(("h06_","h05_","h04_","h03_")) and getattr(module,"__file__",None):
            if not Path(module.__file__).resolve().is_relative_to(ROOT.resolve()): raise ValueError("Research import outside package")
    audit()
    print(json.dumps(dict(status="PASS",stage=stage,scientific_modules_package_only=True,private_Git_or_raw_paths_used=False)))

if __name__=="__main__":
    main()
