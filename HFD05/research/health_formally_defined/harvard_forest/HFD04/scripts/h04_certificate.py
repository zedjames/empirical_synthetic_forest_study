"""Evidence-backed completion gates; unfavorable scientific statuses permitted."""
import json
import h04_common as h

def load(name):return json.loads((h.ROOT/name).read_text())
def count(name):return len(h.read_csv(h.ROOT/name))
def gates():
    truth=h.read_csv(h.ROOT/"boundary_validation/truth_prefixes.csv")
    manifests=[load("figures/H04-"+str(i)+".json") for i in range(1,9)]
    isolation=load("provenance/public_isolation_receipt.json") if (h.ROOT/"provenance/public_isolation_receipt.json").exists() else {}
    receipt=load("reports/verification_receipt.json") if (h.ROOT/"reports/verification_receipt.json").exists() else {}
    full=load("provenance/public_full_replay_summary.json") if (h.ROOT/"provenance/public_full_replay_summary.json").exists() else {}
    extra=load("provenance/public_supplement_replay_summary.json") if (h.ROOT/"provenance/public_supplement_replay_summary.json").exists() else {}
    candidate=h.ROOT/"release/candidate/PaperV_Artifact_Manifest.json"
    return dict(
        A=len(truth)==1458 and len(load("config/world_design.json")["worlds"])==81,
        B=count("boundary_validation/signed_margin_scores.csv")==384 and count("boundary_validation/performance.csv")==9720,
        C=count("context_ablation/estimator_draws.csv")==51840 and count("context_ablation/comparison.csv")==48,
        D=count("context_ablation/prior_reconciled_draws.csv")==2592 and count("context_ablation/prior_reconciled_performance.csv")==486 and len(load("config/context_prior_audit.json")["weights"])==9,
        E=count("juvenile_boundary/measurement_error.csv")==240 and count("juvenile_boundary/masking.csv")==4,
        F=count("juvenile_boundary/propagation.csv")==1536 and load("juvenile_boundary/projection_counterexample.json")["actual_response"][0]["Q2"] is False,
        G=count("effective_count/propagation.csv")==12288 and not load("effective_count/fit_contract.json")["cap_selected"],
        H=load("decision_rules/interpretation.json")["historical_unresolved_total"]==963976 and count("decision_rules/hfd03_sufficient_census.csv")>0,
        I=all(h.sha(h.ROOT/m["data"])==m["data_sha256"] and h.sha(h.ROOT/m["figure"])==m["figure_sha256"] and h.sha(h.ROOT/m["generator"])==m["generator_sha256"] for m in manifests)
          and count("config/numerical_methods.csv")>=175 and (h.ROOT/"reports/Numerical_Methods_Contract.md").is_file(),
        J=isolation.get("status")=="PASS" and candidate.is_file() and isolation.get("artifact_manifest_sha256")==h.sha(candidate)
          and full.get("status")=="PASS" and len(full.get("scientific_outputs",[]))==23 and extra.get("status")=="PASS" and len(extra.get("scientific_outputs",[]))==5,
        K=receipt.get("status")=="PASS" and receipt.get("verifier_sha256")==h.sha(h.REPO/"scripts/verify_hfd04_harvard_forest.py") and len(receipt.get("hostile_controls",[]))>=20 and receipt.get("cache",{}).get("archived_banks_checked")==243,
        L=count("reports/referee_map.csv")==12 and (h.ROOT/"reports/PaperV_V3_Revision_Map.md").is_file())

def scientific_status():
    primary=load("boundary_validation/primary_discrimination.json")
    states=[r["status"] for r in primary.values()]
    return "ESTABLISHED" if all(s=="ESTABLISHED" for s in states) else ("PARTIAL" if any(s!="NOT_ESTABLISHED" for s in states) else "NOT_ESTABLISHED")
