"""Versioned correction of boundary-context mixture, retaining original record."""
import json
from collections import Counter
import h04_common as h
from h04_boundary import world,make_packet
from h04_estimator import Context
from h04_context_prior_estimator import evaluate
def register():
    if (h.ROOT/"context_ablation/prior_reconciled_draws.csv").exists():raise ValueError("Cannot register after supplement outcomes")
    config=dict(id="HFD04-CONTEXT-PRIOR-001",reason="Boundary roster oversamples pressure2; original uniform9 remains a declared sensitivity, not its marginalworld prior.",
      master_seed=202610070404,namespaces=["cw_state","cw_params","cw_future"],worlds=81,draws=32,paths=256,
      prior_source="Entire frozen world roster with equal1/81 design weights; no truth/accuracy enters choice",
      all_worlds_retained=True,prior_interpretation="Experimental-design mixture, NOT ecological frequencies",
      original_results_preserved=True,observation_policy="Same observation packets, exact labels unavailable, new independent estimator streams",
      score_policy="Same Q1/Q2 reference statuses, point classifier>=.5 and metrics; all three truth suites separate")
    worlds=json.loads((h.ROOT/"config/world_design.json").read_text())["worlds"]
    counts=Counter((r["pressure"],r["recovery"]) for r in worlds)
    config["context_counts"]=[counts[(p,r)] for p in range(3) for r in range(3)]
    config["weights"]=[n/81 for n in config["context_counts"]]
    h.write_json("config/context_prior_audit.json",config)
    names=["scripts/h04_context_prior_estimator.py","scripts/h04_context_prior_audit.py","config/context_prior_audit.json","config/world_design.json"]
    h.write_json("provenance/context_prior_registration.json",dict(id=config["id"],phase="REGISTERED_BEFORE_SUPPLEMENT_OUTCOMES",
      sources={name:h.sha(h.ROOT/name) for name in names},original_scientific_registration=json.loads((h.ROOT/"provenance/registration.json").read_text())["configs"],
      outcome_selection=False,prior_counts=counts.total() if hasattr(counts,"total") else sum(counts.values())))
def run():
    record=json.loads((h.ROOT/"provenance/context_prior_registration.json").read_text())
    for name,expected in record["sources"].items():
        if h.sha(h.ROOT/name)!=expected:raise ValueError("Registered supplement changed")
    h.confirm_guard();cfg=json.loads((h.ROOT/"config/context_prior_audit.json").read_text())
    _,_,_,core,model=h.load_data();worlds=json.loads((h.ROOT/"config/world_design.json").read_text())["worlds"]
    anchor=dict(n=core["stats"]["n0"].copy(),d=model["diameter0"].copy());rows=[]
    for identity,cell in enumerate(worlds):
        state,_,_=world(core,model,cell,identity);packet=make_packet(state,identity)
        for row in evaluate(packet,Context("B2_unknown"),cfg["weights"],anchor,model,identity):
            rows.append(dict(world=identity,policy="B2_roster_prior_v2",**row))
        if identity%12==0:print("Prior-reconciled unknown context",identity,flush=True)
    h.write_csv("context_ablation/prior_reconciled_draws.csv",rows)
if __name__=="__main__":
    import sys
    register() if sys.argv[1]=="register" else run()
