"""Isolated publication-candidate entry point; no private source imports."""
import argparse,csv,hashlib,importlib,json,os,sys,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parent
STUDY=ROOT/"research/health_formally_defined/harvard_forest"
FOUR=STUDY/"HFD04"
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def audit():
    manifest=json.loads((ROOT/"PaperV_Artifact_Manifest.json").read_text())
    for name,record in manifest["files"].items():
        path=ROOT/name
        if path.is_symlink() or not path.is_file() or sha(path)!=record["export_sha256"]:raise ValueError("Export identity changed "+name)
    return manifest
def module():
    sys.path.insert(0,str(FOUR/"scripts"))
    import h04_common as h
    # Public guard verifies actual published source; private-only historical
    # ancestry remains an attested hash ledger, not a locally verified git tree.
    def public_guard():
        audit()
        return json.loads((FOUR/"provenance/registration.json").read_text())["protected_sources"]
    h.source_guard=public_guard
    h.PUBLIC_MODE=True
    def sources():
        inventory=json.loads((ROOT/"Public_Inputs.json").read_text())["sources"]
        result={}
        for item in inventory:
            path=ROOT/item["target"]
            if sha(path)!=item["sha256"]:raise ValueError("Public input mismatch "+item["id"])
            if item["id"] in ["HF253-E0","HF253-E1"]:
                result[item["id"]]=dict(item,local_cache_path=str(path))
        return result
    h.frozen.verify_frozen_inputs=sources
    # Research modules must resolve in this package, not the private checkout.
    for name in ["h04_common","h03_common","contracts","empirical","demography","estimator"]:
        imported=importlib.import_module(name)
        if not Path(imported.__file__).resolve().is_relative_to(ROOT.resolve()):
            raise ValueError("Private module import "+name)
    return h
def fixture():
    audit();h=module()
    import numpy as np
    data=json.loads((FOUR/"provenance/lightweight_fixture.json").read_text())
    model={k:np.array(v) if isinstance(v,list) else v for k,v in data["model"].items()}
    state={k:np.array(v) for k,v in data["state"].items()}
    params=h.demography.parameters(model,np.random.default_rng(data["parameters_seed"]))
    bank,lawful=h.demography.simulate(state,params,h.scenario(),5,data["K"],np.random.default_rng(data["futures_seed"]))
    actual=hashlib.sha256(bank.tobytes()+lawful.tobytes()).hexdigest()
    if actual!=data["bank_sha256"]:raise ValueError("Kernel fixture replay differs")
    for query,counts in [("Q1",[8]),("Q2",[8,8])]:
        assert h.classify(True,counts,8,0)["KERNEL_MC_STATUS"]=="TRUE"
        assert h.classify(True,counts,8,1)["FINITE_BANK_HEALTH"]
        assert h.classify(True,counts,8,1)["KERNEL_MC_STATUS"]=="MC_UNRESOLVED"
        assert not h.classify(False,counts,8,0)["FINITE_BANK_HEALTH"]
    witness=json.loads((FOUR/"juvenile_boundary/projection_counterexample.json").read_text())
    assert witness["actual_response"][0]["Q2"] is False and witness["actual_response"][1]["Q2"] is True
    print(json.dumps(dict(status="PASS",isolation="Package-only scientific modules; no private raw inputs/Git",kernel_sha256=actual,numpy=np.__version__,
      scope="Lightweight kernel/query fixture, NOT full raw-data scientific replay")))
def retrieve():
    inventory=json.loads((ROOT/"Public_Inputs.json").read_text())["sources"]
    for item in inventory:
        path=ROOT/item["target"]
        if path.exists():
            if sha(path)!=item["sha256"]:raise ValueError("Existing input mismatch; do not overwrite "+item["id"])
            continue
        path.parent.mkdir(parents=True,exist_ok=True)
        payload=urllib.request.urlopen(item["url"],timeout=120).read()
        if hashlib.sha256(payload).hexdigest()!=item["sha256"]:raise ValueError("Publisher changed source "+item["id"])
        path.write_bytes(payload)
    print("Retrieved only exact public CC0 inputs; no unrestricted downloads")
def validate():
    # Fresh actual fixture and published-table audit; no private raw-data replay.
    fixture();h=module()
    h.write_json("provenance/public_isolation_receipt.json",dict(status="PASS",fixture_only=True,full_public_raw_replay_executed=False,
      artifact_manifest_sha256=sha(ROOT/"PaperV_Artifact_Manifest.json"),origin="Fresh public package-only fixture, not privateGit audit"))
    spec=importlib.util.spec_from_file_location("public_hfd04_verifier",ROOT/"scripts/verify_hfd04_harvard_forest.py")
    audit_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit_module)
    result=audit_module.portable()
    print(json.dumps(dict(status="PASS",published_table_audit=result,private_raw_inputs_used=False,private_git_history_verified=False)))
def masks(h):
    THREE=STUDY/"HFD03"
    config=json.loads((THREE/"config/protocol.json").read_text())
    run=THREE/".runs"/sha(THREE/"config/protocol.json")[:16];run.mkdir(parents=True,exist_ok=True)
    inventory=json.loads((ROOT/"Public_Inputs.json").read_text())["sources"]
    e1=next(ROOT/r["target"] for r in inventory if r["id"]=="HF253-E1")
    with e1.open(newline="") as handle:rows=list(csv.DictReader(handle))
    with (THREE/"reconstruction_validation/mask_assignments.csv").open(newline="") as handle:assignments=list(csv.DictReader(handle))
    identity={"stem.id","tree.id","tag","stem.tag","sp"}
    for family in ["M0","M1","M2","M3"]:
        held={r["stem_id"] for r in assignments if r["mask"]==family and r["replicate"]=="0"}
        clean=[{k:v if k in identity or row["stem.id"] not in held else "NA" for k,v in row.items()} for row in rows]
        with (run/("visible_"+family+"_0.csv")).open("w",newline="") as handle:
            out=csv.DictWriter(handle,fieldnames=list(clean[0]),lineterminator="\n");out.writeheader();out.writerows(clean)
def full(stage):
    audit();h=module();masks(h)
    if stage=="truth":
        import h04_boundary as b;b.run_truth()
    elif stage=="estimate":
        import h04_boundary as b;b.run_estimates()
    elif stage=="juvenile":
        import h04_juvenile as j;j.run()
    elif stage=="cap":
        import h04_cap as c;c.run()
    elif stage=="analysis":
        import h04_analysis as a;a.analyze()
    elif stage=="context-prior":
        import h04_context_prior_audit as c;c.run()
    elif stage=="context-prior-analysis":
        import h04_supplement_analysis as a;a.run()
    else:raise ValueError("Select an explicit registered stage")
if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("mode",choices=["fixture","validate","retrieve","replay"]);parser.add_argument("--stage")
    arg=parser.parse_args()
    if arg.mode=="fixture":fixture()
    elif arg.mode=="validate":validate()
    elif arg.mode=="retrieve":retrieve()
    else:full(arg.stage)
