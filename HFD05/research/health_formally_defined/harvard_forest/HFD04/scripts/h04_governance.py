"""Active scientific-claim validators used by the standalone audit."""
import numpy as np
def exact(expected,actual,message):
    if expected!=actual:raise ValueError(message)
def streams(pilot,final):
    if set(pilot)&set(final):raise ValueError("Pilot/confirmatory stream reuse")
def reference(status,certified):
    if status=="MC_UNRESOLVED" and certified:raise ValueError("Uncertain truth forced certified")
def mechanism(rows):
    if len({r["suite"] for r in rows})!=1:raise ValueError("Mechanisms pooled")
def context(policy,interpretation):
    if policy in {"B1_coarsened","B2_unknown"} and interpretation!="PREDICTION_UNDER_HIDDEN_REALIZED_CONTEXT":
        raise ValueError("Hidden context described as known conditional Health")
def provenance(representation,claim):
    if representation in {"COHORT_RMS","IMPUTED_ORIGIN"} and claim=="OBSERVED_INDIVIDUAL":
        raise ValueError("Derived/imputed distribution mislabeled observation")
def cap_selection(selected,E2_used):
    if selected or E2_used:raise ValueError("Cap selected using exposed E2")
def unconditional_weights(weights,K):
    if len(weights)!=K or not np.allclose(weights,np.full(K,1/K),atol=1e-15,rtol=0):
        raise ValueError("Candidate measure renormalized on successes")
def endpoint(theta,all_success,status):
    if theta==1 and all_success and status=="TRUE":raise ValueError("Finite all success claimed probabilityone")
def mc_status(expected,reported):
    if reported not in {"TRUE","FALSE","MC_UNRESOLVED"} or expected!=reported:
        raise ValueError("Unresolved kernel status forced Boolean")
def claims(value):
    if value.get("scenario_prior") is not None or value.get("semantic_prior") is not None:
        raise ValueError("Ecological prior assigned to design grid")
    if value.get("operational")!="NOT_AUTHORIZED":raise ValueError("Operational assessment claim")
def export(path,allowed):
    from pathlib import PurePosixPath
    p=PurePosixPath(path)
    if path not in allowed or p.is_absolute() or ".." in p.parts or any(x in p.parts for x in [".git",".env",".runs",".runtime",".lake","lean"]):
        raise ValueError("Private/unlicensed/out-of-scope export")
def strata(expected,actual):
    if set(expected)!=set(actual):raise ValueError("Unfavorable registered stratum omitted")

