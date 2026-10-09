"""Freeze accepted HFD04 source/output identity without receipt self-hash loops."""
import json,platform,subprocess,sys
import numpy as np
import h04_common as h
def run():
    h.confirm_guard()
    from h04_certificate import gates
    if not all(gates().values()):raise ValueError("Cannot freeze an incomplete evidence certificate")
    selected={}
    for directory in ["scripts","config","boundary_validation","context_ablation","juvenile_boundary","effective_count","decision_rules","figures","reports"]:
        for path in (h.ROOT/directory).glob("*"):
            if not path.is_file() or path.suffix not in [".py",".json",".csv",".svg",".md"] or path.name=="verification_receipt.json":continue
            selected[str(path.relative_to(h.REPO))]=h.sha(path)
    for name in ["registration.json","context_prior_registration.json","truth_bank_manifest.json","estimator_information.json","lightweight_fixture.json","public_full_replay_summary.json","public_supplement_replay_summary.json"]:
        path=h.ROOT/"provenance"/name;selected[str(path.relative_to(h.REPO))]=h.sha(path)
    verifier=h.REPO/"scripts/verify_hfd04_harvard_forest.py";selected[str(verifier.relative_to(h.REPO))]=h.sha(verifier)
    h.write_json("provenance/execution_manifest.json",dict(schema="HFD04-accepted-source-output-v1",files=dict(sorted(selected.items())),
        execution_start=h.config()["actual_execution_start"],
        registration_commit="8586a917a258f8d748f9bb8c276964a2fd677a1d",supplement_registration_commit="251bd60a26",
        runtime=dict(python=platform.python_version(),numpy=np.__version__,system=platform.system(),release=platform.release(),machine=platform.machine(),compiler=platform.python_compiler()),
        exclusions="This manifest, verification/isolation/private replay/publication receipts, export mirrors, ignored raw banks/cache. No scientific outcome is excluded.",
        source_preservation=h.source_guard(),all_registered_worlds_retained=True))
    print("Frozen",len(selected),"accepted source/output identities",flush=True)
if __name__=="__main__":run()
