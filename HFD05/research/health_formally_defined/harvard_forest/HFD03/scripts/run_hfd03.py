"""Executable frozen study driver. Workstreams can be replayed independently."""
import argparse
import json
from h03_common import ROOT, RUN, PROTOCOL_HASH, PROTOCOL, inputs, empirical, frozen_guard, write_json


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--stage",choices=["all","calibration","entry","missingness","diagnostics","synthetic","mc","reports"],default="all")
    args=parser.parse_args()
    registration=json.loads((ROOT/"provenance/registration.json").read_text())
    if registration["protocol_sha256"]!=PROTOCOL_HASH or registration["frozen_sources"]!=frozen_guard():
        raise ValueError("Registration or frozen baseline changed")
    manifest_path=ROOT/"provenance/execution_manifest.json"
    if manifest_path.exists():
        from h03_common import frozen
        manifest=json.loads(manifest_path.read_text())
        if manifest["protocol_sha256"]!=PROTOCOL_HASH or any(
                frozen.sha(ROOT/name)!=expected for name,expected in manifest["sources"].items()):
            raise ValueError("Frozen HFD03 implementation changed; cannot overwrite its scientific run")
    RUN.mkdir(parents=True,exist_ok=True)
    sources,e0,e1,trees=inputs()
    core=empirical.read_core(sources)
    model=empirical.fit(core)
    if args.stage in {"entry","all"}:
        import h03_entry
        h03_entry.run(e0,e1,trees)
    if args.stage in {"calibration","all"}:
        import h03_calibration
        h03_calibration.run(sources,e0,e1)
    if args.stage in {"missingness","all"}:
        import h03_missingness
        h03_missingness.run(core,model)
    if args.stage in {"diagnostics","all"}:
        import h03_diagnostics
        h03_diagnostics.run(core,model,json.loads((ROOT/"recruitment_audit/summary.json").read_text()))
    if args.stage in {"synthetic","all"}:
        import h03_synthetic
        h03_synthetic.run(core,model)
    if args.stage in {"mc","all"}:
        import h03_mc
        h03_mc.run(core,model)
    if args.stage in {"reports","all"}:
        import h03_reports
        h03_reports.run(core,model)
    print("HFD03 stage",args.stage,"complete; protocol",PROTOCOL_HASH)


if __name__=="__main__":
    main()
