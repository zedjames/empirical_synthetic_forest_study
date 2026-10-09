"""Explicit final implementation/output freeze, not automatic verifier repair."""
from pathlib import Path
from h03_common import ROOT, REPO, PROTOCOL_HASH, frozen, write_json


def main():
    code={str(p.relative_to(ROOT)):frozen.sha(p) for p in sorted(ROOT.rglob("*.py"))
          if not any(part.startswith(".") or part=="__pycache__" for part in p.relative_to(ROOT).parts)}
    top=REPO/"scripts/verify_hfd03_harvard_forest.py"
    outputs={}
    excluded={"provenance/execution_manifest.json","provenance/hostile_controls.json",
              "reports/hfd03_contract.json","reports/verification_receipt.json",
              "provenance/publication_receipt.json","provenance/column_provenance.json"}
    columns={}
    source_map={
        "reconstruction_validation":(["HF253-E0","HF253-E1","HFD03-mask-protocol"],"leave-observed-out"),
        "recruitment_audit":(["HF253-E0","HF253-E1","HF253-trees2014","HF253-species-codes","HFD03-entry-assumptions"],"candidate-identity-audit"),
        "missingness":(["HFD02S-frozen-fit","HFD03-B0-B4"],"unconditional-full-pipeline"),
        "synthetic_validation":(["HF-SYNTH-DISC-001","HFD03-observation-operator"],"full-grid-truth-blind-estimation"),
        "monte_carlo":(["HFD03-paired-maximum-banks"],"Wilson-and-state-block-prefix"),
        "e2_failure_localization":(["HFD02S-frozen-fit","HF-third-census-news","HFD03-diagnostic-grid"],"untuned-factor-surface"),
        "figures":(["HFD03-committed-outcomes"],"derived-reader-visualization"),
        "provenance":(["HF253-E0","HF253-E1","HFD03-protocol"],"ancestry-and-namespace-ledger"),
        "specification":(["HFD02S-frozen-design","HFD03-protocol"],"exact-config-exposure"),
        "uncertainty":(["HFD02S-frozen-decomposition"],"noncommensurate-finite-contrast")}
    for p in sorted(ROOT.rglob("*")):
        rel=p.relative_to(ROOT)
        if (not p.is_file() or any(part.startswith(".") or part=="__pycache__" for part in rel.parts)
                or p.suffix==".py" or str(rel) in excluded):
            continue
        outputs[str(rel)]={"sha256":frozen.sha(p),"bytes":p.stat().st_size}
        if p.suffix==".csv":
            import csv
            with p.open(newline="") as handle:
                fields=next(csv.reader(handle))
            sources,method=source_map.get(rel.parts[0],(["HFD03-protocol"],"derived-study-output"))
            columns[str(rel)]={}
            for field in fields:
                cls="DERIVED"
                if rel.parts[0]=="recruitment_audit" and field in {
                    "stem_id","tree_id","tag","stem_tag","taxon","gx","gy","quadrat","dbh","pom","hom","date",
                    "publisher_latin","publisher_synonym","e0_parent_taxon","e0_parent_gx","e0_parent_gy",
                    "e0_parent_quadrat","e0_parent_dbh","e0_parent_pom","e0_parent_hom","e0_parent_date","e0_parent_latin"}:
                    cls="OBSERVED"
                if field.startswith("entry_") or field in {"theta","alpha","beta","tau","rho","reference_weight"}:
                    cls="EXTERNALLY_CONSTRAINED"
                if field.startswith("truth_") and rel.parts[0]=="synthetic_validation":
                    cls="SIMULATED"
                columns[str(rel)][field]=dict(provenance_class=cls,source_ids=sources,
                    method_id=method,namespace="HFD03-derived-study",
                    limitation="Provenance of output column; underlying observation/imputation/simulation graph is retained, not replaced by this class")
    write_json("provenance/column_provenance.json",columns)
    outputs["provenance/column_provenance.json"]={
        "sha256":frozen.sha(ROOT/"provenance/column_provenance.json"),
        "bytes":(ROOT/"provenance/column_provenance.json").stat().st_size}
    docs={str(p.relative_to(REPO)):frozen.sha(p) for p in sorted((REPO/"docs").glob("Health_Harvard_Forest_HFD03_*.md"))}
    write_json("provenance/execution_manifest.json",dict(schema="HFD03-execution-v1",
        protocol_sha256=PROTOCOL_HASH,numpy_version="2.0.2",sources=code,
        verifier_sha256=frozen.sha(top),outputs=outputs,docs=docs,
        registration="Protocol predeclared before comparisons; final implementation/output identities frozen after documented adapter/reporting corrections",
        dependency_policy="Read-only HFD02S modules, HFD01 sources and existing private NumPy runtime"))
    print("Final HFD03 implementation/output freeze:",len(code),"sources",len(outputs),"outputs")


if __name__=="__main__":
    main()
