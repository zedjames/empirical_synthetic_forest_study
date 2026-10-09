"""Production audit predicates, also exercised with adversarial evidence."""
import itertools
import numpy as np
from h03_common import PROTOCOL, classify


def protected(expected,actual):
    if expected!=actual:
        raise ValueError("Frozen source identity differs")


def diagnosis(evidence):
    if evidence.get("parameter_selection") or evidence.get("original_replaced"):
        raise ValueError("E2 optimization is not a diagnostic correction")
    if evidence.get("aligned_failure_vector") is not None:
        raise ValueError("Public E2 anchors are not aligned exact validation targets")


def complete_world_grid(rows):
    cfg=PROTOCOL["synthetic"]
    expected=set(itertools.product(range(3),range(3),range(3),range(3),range(cfg["replicates_per_cell"])))
    correct=[r for r in rows if r["suite"]=="correct"]
    coordinates=[tuple(int(r[k]) for k in ("structure_level","regeneration_level","pressure_level",
                                          "recovery_level","replicate")) for r in correct]
    if len(coordinates)!=len(expected) or set(coordinates)!=expected:
        raise ValueError("World pruning, duplication or outcome selection")
    if set(r["suite"] for r in rows)!=set(cfg["misspecifications"]):
        raise ValueError("Correct and controlled-mismatch suites must remain separate")
    worlds={int(r["world"]) for r in correct}
    if worlds!=set(range(len(expected))) or len(rows)!=len(expected)*3:
        raise ValueError("Suite/world census incomplete")


def claims(claim):
    if claim["scenario_prior"] is not None or claim["semantic_prior"] is not None:
        raise ValueError("No licensed scenario/semantic priors")
    if claim["weighting"]!="finite design-mixture":
        raise ValueError("Mixture weights not ecological posterior")
    if claim["lawful"]!="model-admissibility realization of Lawful":
        raise ValueError("Software admissibility inflated to ecological validation")
    if claim["wet"]!="wet-candidate/protocol-proxy":
        raise ValueError("No certified hydrological mask")
    if claim["M3"]!="redundant supporting query-sufficiency illustration; no information gain":
        raise ValueError("Redundant representation cannot manufacture information")


def capacity(count,total,reported):
    if not 0<=count<=total or total<=0 or abs(count/total-reported)>1e-12:
        raise ValueError("Failure denominator dropped or successful paths renormalized")


def mc_cell(row):
    q2=row["query"]=="Q2"
    counts=[int(row["viable_count"])]
    if q2:
        counts.append(int(row["reserve_count"]))
    present=row["present"] is True or row["present"]=="True"
    theta=float(row["theta"]);K=int(row["K"])
    status,lower,upper=classify(present,counts,K,theta)
    if row["status"]!=status or abs(float(row["lower"])-min(lower))>1e-12 or abs(float(row["upper"])-min(upper))>1e-12:
        raise ValueError("Threshold-overlapping interval forced or altered")
    capacity(counts[0],K,float(row["viable_mass"]))
    return status


def stopping(config):
    if config["optional_extension"] is not None or config["K"]!=[64,256,1024]:
        raise ValueError("No preferred-Boolean adaptive stopping authorized")


def entry(row):
    if row["category"] not in PROTOCOL["entry_audit"]["reference_class_weights"]:
        raise ValueError("New ID is not automatically true ingrowth")
    if row.get("biological_recruitment_confirmed"):
        raise ValueError("Candidate identity insufficient for biological confirmation")
