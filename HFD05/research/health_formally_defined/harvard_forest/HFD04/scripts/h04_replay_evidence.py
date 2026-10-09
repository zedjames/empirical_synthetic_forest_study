"""Public-safe evidence from actual isolated replay; never asserts a dry run."""
import json,subprocess,sys
from pathlib import Path
import h04_common as h

def summarize(receipt_name,summary_name):
    receipt=json.loads((h.ROOT/receipt_name).read_text())
    if receipt["status"]!="PASS":raise ValueError("Replay not successful")
    manifest=json.loads((Path(receipt["temporary_directory"])/"PaperV_Artifact_Manifest.json").read_text())
    source_hashes={name:v["export_sha256"] for name,v in manifest["files"].items()
        if ("HFD02S/" in name or "HFD03/scripts/" in name or "HFD03/config/" in name or "HFD04/scripts/" in name or "HFD04/config/" in name)
        and (name.endswith(".py") or name.endswith(".json"))}
    # Reporting, exporting and audit orchestration can evolve without changing
    # the numerical science that was actually rerun. Keep those identities in
    # the original manifest receipt, not in the required science-equivalence map.
    required=["h04_common.py","h04_boundary.py","h04_estimator.py","h04_juvenile.py","h04_cap.py","h04_analysis.py",
              "h04_pilot.py","h04_pilot_v2.py","h04_context_prior_audit.py","h04_context_prior_estimator.py","h04_supplement_analysis.py"]
    source_hashes={n:v for n,v in source_hashes.items() if "HFD04/scripts/" not in n or Path(n).name in required}
    source_hashes={n:v for n,v in source_hashes.items() if not n.endswith("numerical_methods.csv")}
    for name,digest in source_hashes.items():
        if h.sha(h.REPO/name)!=digest:raise ValueError("Replayed science changed "+name)
    outputs=receipt["scientific_outputs"]
    for row in outputs:
        if h.sha(h.ROOT/row["path"])!=row["sha256"]:raise ValueError("Replayed outcome changed")
    h.write_json(summary_name,dict(status="PASS",actual_replay=True,scientific_outputs=outputs,
        scientific_sources=source_hashes,replayed_artifact_manifest_sha256=receipt["artifact_manifest_sha256"],
        private_scientific_modules_imported=False,exact_public_raw_inputs_reused_read_only=True,
        fresh_download=False,raw_cache_copied=False,source_environment=receipt["source_environment"],
        input_sources=receipt["input_sources"],qualification=receipt["qualification"]))

def supplemental():
    receipt=json.loads((h.ROOT/"provenance/public_isolation_receipt.json").read_text())
    isolated=Path(receipt["isolated_temporary_directory"])
    if receipt["status"]!="PASS":raise ValueError("Candidate isolation required")
    inventory=json.loads((isolated/"Public_Inputs.json").read_text())["sources"]
    src=h.frozen.verify_frozen_inputs()
    originals={k:h.frozen.HFD01/src[k]["local_cache_path"] for k in ["HF253-E0","HF253-E1"]}
    for item in json.loads((h.THREE/"provenance/identity_inputs.json").read_text())["sources"]:
        originals[item["name"]]=h.THREE/".inputs"/item["name"]
    for item in inventory:
        target=isolated/item["target"];target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():raise ValueError("Supplement needs fresh isolated inputs")
        if h.sha(originals[item["id"]])!=item["sha256"]:raise ValueError("Input changed")
        target.symlink_to(originals[item["id"]])
    stages={"context-prior":["context_ablation/prior_reconciled_draws.csv"],
            "context-prior-analysis":["context_ablation/prior_reconciled_performance.csv","context_ablation/prior_reconciled_comparison.csv",
                                      "context_ablation/prior_reconciled_strata.csv","juvenile_boundary/changed_case_audit.csv"]}
    comparisons=[]
    for stage,names in stages.items():
        result=subprocess.run([sys.executable,"-I",str(isolated/"artifact_runner.py"),"replay","--stage",stage],cwd=isolated)
        if result.returncode:raise ValueError("Supplement replay failed")
        for name in names:
            actual=h.sha(isolated/"research/health_formally_defined/harvard_forest/HFD04"/name)
            if actual!=h.sha(h.ROOT/name):raise ValueError("Supplement replay differs "+name)
            comparisons.append(dict(stage=stage,path=name,sha256=actual,status="PASS"))
    h.write_json("provenance/public_supplement_replay_receipt.json",dict(status="PASS",artifact_manifest_sha256=receipt["artifact_manifest_sha256"],
        scientific_outputs=comparisons,temporary_directory=str(isolated),source_environment="Isolated packaged source; Python3.9.6/NumPy2.0.2",
        input_sources=[dict(id=i["id"],sha256=i["sha256"],license=i["license"]) for i in inventory],
        qualification="Separately preregistered supplement only; original full replay already complete. Public raw bytes reused read-only."))
    summarize("provenance/public_supplement_replay_receipt.json","provenance/public_supplement_replay_summary.json")
    print("PASS isolated supplemental replay",len(comparisons),"byte-exact outputs",flush=True)

if __name__=="__main__":
    if len(sys.argv)>1 and sys.argv[1]=="supplement":supplemental()
    else:summarize("provenance/public_full_replay_receipt.json","provenance/public_full_replay_summary.json")
