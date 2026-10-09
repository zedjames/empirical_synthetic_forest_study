"""Explicit pre-confirmatory registration; refuses existing outcome artifacts."""
import json,sys
import h04_common as h
def main():
    if (h.ROOT/"boundary_validation/truth_prefixes.csv").exists():
        raise ValueError("Cannot preregister after confirmatory outcomes")
    pilot=json.loads((h.ROOT/"boundary_validation/pilot_v2_design.json").read_text())
    worlds=[dict(cell,replicate=0) for cell in pilot["cells"]]
    for p in range(3):
        for r in range(3):worlds.append(dict(family="REALIZATION_CONTROL",pressure=p,recovery=r,structure=.55,regeneration=.2,juvenile_diameter=None,target=None,replicate=0))
    h.write_json("config/world_design.json",dict(id="HF-SYNTH-BOUNDARY-001",derived_from=pilot["id"],selection="All72 registered pilot cells plus9 exact realization controls, no confirmatory pruning",worlds=worlds))
    files={str(p.relative_to(h.ROOT)):h.sha(p) for p in sorted((h.ROOT/"scripts").glob("*.py"))}
    configs={str(p.relative_to(h.ROOT)):h.sha(p) for p in sorted((h.ROOT/"config").glob("*.json"))}
    h.write_json("provenance/registration.json",dict(execution_start=h.config()["actual_execution_start"],
       protected_sources=h.source_guard(),scientific_sources=files,configs=configs,
       pilot_outputs={str(p.relative_to(h.ROOT)):h.sha(p) for p in sorted((h.ROOT/"boundary_validation").glob("pilot*"))},
       phase="PILOTS_COMPLETE_CONFIRMATORY_NOT_STARTED",all_worlds_retained=True,
       environment=dict(python=sys.version,numpy=__import__("numpy").__version__),runtime_reused=True))
    print("Registered",len(worlds),"worlds; frozen sources and configs",flush=True)
if __name__=="__main__":main()
