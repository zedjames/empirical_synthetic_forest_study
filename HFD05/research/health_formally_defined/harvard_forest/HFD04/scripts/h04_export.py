"""Deny-by-default, explicitly enumerated local publication-review candidate."""
import json,re,shutil,subprocess,sys,tempfile
from pathlib import Path
import h04_common as h
from h04_governance import export
def build():
    out=h.ROOT/"release/candidate"
    # Enumerate publication scope before copying. No repository-wide traversal.
    selected={}
    for study,directories in {
      "HFD02S":["scripts","models","validation","config"],
      "HFD03":["scripts","config","synthetic_validation"],
      "HFD04":["scripts","config","boundary_validation","context_ablation","juvenile_boundary","effective_count","decision_rules","figures","reports"]}.items():
        base=h.ROOT.parent/study
        for directory in directories:
            for path in sorted((base/directory).glob("*")):
                if path.is_file() and path.suffix in [".py",".json",".csv",".svg",".md"]:
                    name="research/health_formally_defined/harvard_forest/"+study+"/"+str(path.relative_to(base))
                    selected[name]=path
    explicit=[
      h.THREE/"reconstruction_validation/mask_assignments.csv",
      h.THREE/"provenance/identity_inputs.json",
      h.ROOT/"provenance/registration.json",
      h.ROOT/"provenance/context_prior_registration.json",
      h.ROOT/"provenance/truth_bank_manifest.json",
      h.ROOT/"provenance/estimator_information.json",
      h.ROOT/"provenance/lightweight_fixture.json"]
    explicit += [path for path in [h.ROOT/"provenance/public_full_replay_summary.json",h.ROOT/"provenance/public_supplement_replay_summary.json"] if path.exists()]
    freeze=h.ROOT/"provenance/execution_manifest.json"
    if freeze.exists():explicit.append(freeze)
    # Include every frozen HFD03 manuscript table/figure/report surface, not
    # its raw simulation banks or private publication/worktree receipt.
    previous=json.loads((h.THREE/"provenance/execution_manifest.json").read_text())
    explicit += [h.THREE/name for name in previous["outputs"]]
    for path in explicit:
        selected[str(path.relative_to(h.REPO))]=path
    selected["scripts/verify_hfd04_harvard_forest.py"]=h.REPO/"scripts/verify_hfd04_harvard_forest.py"
    # Original source guard ancestry is included as hashes, not an assertion
    # that every old private record has been independently replayed publicly.
    selected["artifact_runner.py"]=h.ROOT/"scripts/artifact_runner.py"
    allowed=set(selected);manifest={}
    for name,path in selected.items():
        export(name,allowed)
        content=path.read_bytes()
        text=content.decode("utf-8") if path.suffix!=".npz" else ""
        if re.search(r"/Users/[A-Za-z][^\s\"']+",text):
            raise ValueError("Internal absolute path requires explicit sanitization "+name)
        if re.search(r"AKIA[A-Z0-9]{16}|-----BEGIN (?:RSA |OPENSSH )?PRIVATE KEY-----",text):
            raise ValueError("Secret-looking export content "+name)
        target=out/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(content)
        manifest[name]=dict(original_sha256=h.sha(path),export_sha256=h.sha(target),bytes=len(content),
           provenance="PROJECT_AUTHORED_PUBLICATION_REVIEW_SOURCE" if path.suffix==".py" else "PUBLIC_DATA_DERIVED_OR_REGISTERED_PUBLICATION_RECORD",
           redistribution="Owner disclosure/license review pending; no public publication by this tranche")
    sources=h.frozen.verify_frozen_inputs()
    inputs=[]
    for key in ["HF253-E0","HF253-E1"]:
        item=sources[key]
        inputs.append(dict(id=key,url=item["retrieval_url"],sha256=item["sha256"],
                           target="research/health_formally_defined/harvard_forest/HFD01/"+item["local_cache_path"],license="CC0-1.0"))
    for item in json.loads((h.THREE/"provenance/identity_inputs.json").read_text())["sources"]:
        inputs.append(dict(id=item["name"],url=item["url"],sha256=item["sha256"],
                           target="research/health_formally_defined/harvard_forest/HFD03/.inputs/"+item["name"],license="CC0-1.0"))
    payload={"sources":inputs,"citation":"Orwig D, Foster D, Ellison A2023. HF253v6. DOI10.6073/pasta/818789a882a318c1d7f3fc43a2289e12",
      "license_source":"https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF253"}
    target=out/"Public_Inputs.json";target.write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
    manifest["Public_Inputs.json"]=dict(export_sha256=h.sha(target),bytes=target.stat().st_size,provenance="PUBLIC_CANONICAL_INPUT_RETRIEVAL",redistribution="CC0 metadata")
    runtime=out/"requirements.txt";runtime.write_text("numpy==2.0.2\n")
    manifest["requirements.txt"]=dict(export_sha256=h.sha(runtime),bytes=runtime.stat().st_size,provenance="PINNED_RUNTIME",redistribution="Dependency instruction only")
    readme=out/"README.md";readme.write_text("""# Paper V source candidate — REVIEW ONLY

Python3.9.6 and NumPy2.0.2 were used. Install requirements in your own isolated environment.

Lightweight verification: `python -I artifact_runner.py validate` (fixture plus complete published-table audit).
This reads only packaged source/model/table fixtures and verifies independent source-module location, archived kernel hash, exact/kernel endpoints, query-sufficiency counterexample, all published experimental rosters/decisions/confusion scores and export identities. It does NOT establish full empirical raw-data replay or independently verify unavailable privateGit history.

Public-data replay: `python -I artifact_runner.py retrieve`.
Then replay stages truth, estimate, juvenile, cap, analysis, context-prior, context-prior-analysis with `python -I artifact_runner.py replay --stage STAGE`. Truth precedes estimation; context-prior precedes its analysis. The retrieval list contains only four publicly hosted HF253v6 CC0 tables, exact SHA checked. Full bank sizes and all world cells are in config/protocol.json. Raw banks/caches/Git are not distributed. A separate HFD03 decision-census file exposes historical exact counts/statuses; that presentation can be checked without private banks. Public-safe replay summaries document the actually completed full packaged-source replay: 23 original scientific outputs and five separately registered supplemental outputs. Raw bytes were obtained externally earlier and reused read-only, not freshly downloaded during replay. Source/config and output hashes allow comparison against the supplied golden files. Full reproduction regenerates outcomes; validate checks published evidence and does not itself run all banks.

Public adapter replaces only private Git ancestry checking and canonical raw-data location discovery. Exported scientific files are byte-exact; original private-history preservation is an attested hash ledger, not independently reconstructed privateGit history. Actual public source identities are verified against this inclusion manifest.

No unrelated research, Lean formal-world code, credentials or private raw datasets are included. The original repository stays private. HF253 dataset is CC0; cite Orwig/Foster/Ellison2023 HF253v6. Project-authored source disclosure and distribution license require separate owner review. This is a locally validated review candidate, NOT a public release and not aDOI.
""")
    manifest["README.md"]=dict(export_sha256=h.sha(readme),bytes=readme.stat().st_size,provenance="REPRODUCTION_AND_SCOPE",redistribution="Review only")
    record=dict(schema="PaperV-HFD04-public-candidate-v1",files=manifest,explicit_inclusion_list=sorted(manifest),
      deny_by_default=True,raw_data_redistributed=False,unrelated_IP_exported=False,public_publication=False,DOI_created=False,
      source_license="OWNER_REVIEW_REQUIRED_BEFORE_PUBLIC_PUBLICATION",runtime=dict(python="3.9.6",numpy="2.0.2"),
      scientific_registration=h.config()["actual_execution_start"],fixture_scope="Lightweight check; separate full scientific replay summaries retain executed evidence",
      adapter_scope="Private Git/source ancestry audit and canonical raw-data locations only; no scientific model function transformed")
    target=out/"PaperV_Artifact_Manifest.json";target.write_text(json.dumps(record,indent=2,sort_keys=True)+"\n")
    h.write_json("release/PaperV_Artifact_Manifest.json",record)
    # Copy only the explicit small candidate to a fresh isolated temporary dir.
    isolated=Path(tempfile.mkdtemp(prefix="hfd04-artifact-",dir="/private/tmp"))
    for name in [*manifest,"PaperV_Artifact_Manifest.json"]:
        target=isolated/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(out/name,target)
    result=subprocess.run([sys.executable,"-I",str(isolated/"artifact_runner.py"),"validate"],cwd=isolated,text=True,capture_output=True)
    receipt=dict(status="PASS" if result.returncode==0 else "FAIL",returncode=result.returncode,stdout=result.stdout,stderr=result.stderr,
       fixture_only=True,published_table_audit=True,full_public_raw_replay_executed=False,all_included_files_hash_checked=True,private_scientific_imports=False,
       artifact_manifest_sha256=h.sha(out/"PaperV_Artifact_Manifest.json"),files=len(manifest),
       source_bytes=sum(v["bytes"] for v in manifest.values()),isolated_temporary_directory=str(isolated),
       privacy_audit="Explicit scope, no secrets/internal user paths, rawdata not redistributed; source license separately reviewed")
    h.write_json("provenance/public_isolation_receipt.json",receipt)
    if result.returncode:raise ValueError(result.stdout+"\n"+result.stderr)
    print("Isolated public candidate fixture PASS",len(manifest),"files",receipt["source_bytes"],"bytes",flush=True)
if __name__=="__main__":build()
