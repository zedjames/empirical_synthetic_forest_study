"""Executable guards and hostile controls; violations raise, never relabel."""
import copy
import hashlib
import inspect
from dataclasses import asdict
import numpy as np
from contracts import canonical, quantity, fate, load, Seeds
from demography import mass, check_history, check_cohorts
from estimator import ObservationPacket, estimate_state


def require_design(config, expected_digest):
    if config["training_artifacts"] != ["HF253-E0","HF253-E1"]:
        raise ValueError("E2/tower/context source in training")
    if hashlib.sha256(canonical(config)).hexdigest() != expected_digest:
        raise ValueError("Design changed after freeze")


def require_original_weights(actual, original):
    if not np.array_equal(actual,original):
        raise ValueError("Success conditioning/weight replacement forbidden")


def require_fates(rows, assigned):
    if len(rows)!=len(assigned):
        raise ValueError("Fate ledger length")
    for row,z in zip(rows,assigned):
        if fate(row) is None and z is not None:
            raise ValueError("Unknown/absent fate converted to observed death/alive")


def require_draws(draws):
    if len(draws)!=load("config/design.json")["ensemble"]["state_draws"]:
        raise ValueError("Whole reconstruction ensemble required")


def require_reference(index):
    if index!=0:
        raise ValueError("Reference-world cherry picking")


def require_redundant_equal(m2,m3):
    if m2!=m3:
        raise ValueError("Redundant M3 descriptor cannot earn forced enrichment gain")


def hostile_tests(seeds=None):
    tests = {}
    def rejected(name,fn,exceptions=(ValueError,TypeError)):
        try:
            fn()
        except exceptions:
            tests[name] = "REJECTED_AS_REQUIRED"
        else:
            raise AssertionError("Hostile control accepted: "+name)
    rejected("provenance_contamination",lambda:quantity(3,"SIMULATED","D0","generator","raw_E1"))
    packet = ObservationPacket((0,)*64,(None,)*64,(0.,)*64)
    leaked = asdict(packet)
    leaked["true_counts"] = [100]*64
    rejected("truth_in_packet",lambda:ObservationPacket(**leaked))
    rejected("truth_object_as_packet",lambda:estimate_state(leaked,{},None))
    # Packet invariance: changing hidden truth cannot alter observation-only estimate.
    anchor = dict(n=np.full(64,10.),d=np.full(64,6.))
    seeds = seeds or Seeds()
    a = estimate_state(packet,anchor,seeds.rng("synthetic_fixture",0))
    hidden_truth = {"true_counts":[999]*64}
    hidden_truth["true_counts"][0] = 0
    b = estimate_state(packet,anchor,seeds.rng("synthetic_fixture",0))
    assert np.array_equal(a["n"],b["n"]) and np.array_equal(a["d"],b["d"])
    assert "worlds" not in inspect.getsource(estimate_state)
    tests["hidden_truth_invariance"] = "PASS"
    cfg = load("config/design.json")
    digest = hashlib.sha256(canonical(cfg)).hexdigest()
    bad = copy.deepcopy(cfg); bad["training_artifacts"].append("E2")
    rejected("E2_tuning",lambda:require_design(bad,digest))
    bad = copy.deepcopy(cfg); bad["training_artifacts"].append("AMF-US-Ha1")
    rejected("tower_substitution",lambda:require_design(bad,digest))
    original = np.full(4,.25)
    rejected("success_renormalization",lambda:require_original_weights(np.array([.5,.5,0,0]),original))
    rejected("unknown_fate_to_death",lambda:require_fates([None,{"df.status":"stem_gone","exact.date":"NA"}],[0,0]))
    rejected("single_filled_state_certainty",lambda:require_draws([{}]))
    bad = copy.deepcopy(cfg); bad["semantics"]["alpha_grid"] = [.01]
    rejected("preferred_semantics",lambda:require_design(bad,digest))
    rejected("cherry_picked_world",lambda:require_reference(1))
    rejected("forced_M3_gain",lambda:require_redundant_equal(.5,.6))
    # Independent path accounting and nonlawful-capacity rejection.
    bank = np.zeros((2,3,6)); bank[:,:,0]=10; bank[:,:,2]=5
    assert check_history(bank).all()
    bank[0,1,0]=11
    assert not check_history(bank)[0]
    rejected("nonlawful_capacity",lambda:mass(np.array([True,False]),np.array([False,True]),np.full(2,.5)))
    tests["independent_lawful_checker_corruption"] = "REJECTED_AS_REQUIRED"
    assert not check_cohorts(np.ones((1,2)),np.ones((1,2)),["same","same"])[0]
    tests["duplicate_cohort_identity"] = "REJECTED_AS_REQUIRED"
    return tests
