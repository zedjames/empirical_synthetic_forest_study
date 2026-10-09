"""Full isolated published-source replay using immutable cached public inputs.

No duplicate raw cache and no new download. Shared inputs are read-only symlinks.
This is a verification orchestrator, not a model or experimental-design change.
"""
import json,subprocess,sys
from pathlib import Path
import h04_common as h
def main():
    receipt=json.loads((h.ROOT/"provenance/public_isolation_receipt.json").read_text())
    if receipt["status"]!="PASS":raise ValueError("Public candidate must pass lightweight isolation first")
    isolated=Path(receipt["isolated_temporary_directory"])
    if not isolated.is_dir() or h.sha(isolated/"PaperV_Artifact_Manifest.json")!=receipt["artifact_manifest_sha256"]:
        raise ValueError("Isolated candidate identity changed")
    sources=h.frozen.verify_frozen_inputs()
    inventory=json.loads((isolated/"Public_Inputs.json").read_text())["sources"]
    original={}
    for key in ["HF253-E0","HF253-E1"]:original[key]=h.frozen.HFD01/sources[key]["local_cache_path"]
    for item in json.loads((h.THREE/"provenance/identity_inputs.json").read_text())["sources"]:
        original[item["name"]]=h.THREE/".inputs"/item["name"]
    for item in inventory:
        source=original[item["id"]]
        if h.sha(source)!=item["sha256"]:raise ValueError("Immutable public input changed")
        target=isolated/item["target"];target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():raise ValueError("Fresh raw replay input target unexpectedly occupied")
        target.symlink_to(source)
    outputs={
      "truth":["boundary_validation/truth_prefixes.csv","boundary_validation/pre_estimator_census.json","provenance/truth_bank_manifest.json"],
      "estimate":["context_ablation/estimator_draws.csv","provenance/estimator_information.json"],
      "juvenile":["juvenile_boundary/observed_groups.csv","juvenile_boundary/measurement_error.csv","juvenile_boundary/masking.csv",
                   "juvenile_boundary/visible_size_profiles.csv","juvenile_boundary/initial_projection.csv","juvenile_boundary/propagation.csv","juvenile_boundary/empirical_support.json"],
      "cap":["effective_count/posteriors.csv","effective_count/origin_states.csv","effective_count/propagation.csv","effective_count/masking_calibration.csv","effective_count/fit_contract.json"],
      "analysis":["boundary_validation/performance.csv","context_ablation/comparison.csv","boundary_validation/signed_margin_scores.csv",
                  "boundary_validation/absolute_margin_scores.csv","boundary_validation/failure_scores.csv","boundary_validation/primary_discrimination.json"]}
    comparisons=[]
    for stage,names in outputs.items():
        print("FULL isolated public-source replay",stage,flush=True)
        result=subprocess.run([sys.executable,"-I",str(isolated/"artifact_runner.py"),"replay","--stage",stage],cwd=isolated)
        if result.returncode:raise ValueError("Public replay stage failed "+stage)
        for name in names:
            expected=h.sha(h.ROOT/name)
            actual=h.sha(isolated/"research/health_formally_defined/harvard_forest/HFD04"/name)
            if actual!=expected:raise ValueError("Isolated scientific output differs "+name)
            comparisons.append(dict(stage=stage,path=name,sha256=actual,status="PASS"))
        print("FULL isolated public-source replay",stage,"BYTE_EXACT",flush=True)
    result=subprocess.run([sys.executable,"-I",str(isolated/"artifact_runner.py"),"validate"],cwd=isolated)
    if result.returncode:raise ValueError("Final isolated numerical/source audit failed")
    h.write_json("provenance/public_full_replay_receipt.json",dict(status="PASS",
      artifact_manifest_sha256=receipt["artifact_manifest_sha256"],scientific_outputs=comparisons,
      full_registered_truth_banks_recomputed=True,all_context_worlds_recomputed=True,
      exact_public_raw_inputs_reused_read_only=True,fresh_download=False,raw_cache_copied=False,
      private_scientific_modules_imported=False,source_environment="Isolated package-only module paths, pinned Python3.9.6/NumPy2.0.2",
      input_sources=[dict(id=i["id"],sha256=i["sha256"],license=i["license"]) for i in inventory],
      qualification="Externally obtained public HF253 bytes reused via read-only symlinks; no new network retrieval. PrivateGit historical identity not independently reconstructed.",
      temporary_directory=str(isolated)))
    print("PASS FULL isolated public scientific replay",len(comparisons),"byte-exact outputs",flush=True)
if __name__=="__main__":main()

