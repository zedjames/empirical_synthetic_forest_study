"""Observation-only synthetic estimator; no truth generator import."""
import math
from dataclasses import dataclass
import numpy as np
from h03_common import PROTOCOL, demographic, rng
from estimator import ObservationPacket, estimate_state


@dataclass(frozen=True)
class PublicContext:
    pressure_level: int
    recovery_level: int
    def __post_init__(self):
        if self.pressure_level not in range(3) or self.recovery_level not in range(3):
            raise ValueError("Undeclared observed experimental design covariate")


def public_parameters(model,context,generator):
    cfg=PROTOCOL["synthetic"]
    params=demographic.parameters(model,generator)
    base=cfg["pressure_annual_hazard"][context.pressure_level]
    params["hazard"]=base*np.exp(generator.normal(-cfg["hazard_prior_log_sd"]**2/2,
                                                 cfg["hazard_prior_log_sd"],64))
    entry=cfg["recovery_entry_fraction"][context.recovery_level]*model["baseline_juveniles"]
    shares=model["recruitment_mean"]+1/64
    shares=shares/shares.sum()
    params["recruit"]=entry*shares*np.exp(generator.normal(-cfg["entry_prior_log_sd"]**2/2,
                                                          cfg["entry_prior_log_sd"],64))
    return params


def estimator_rng(namespace,*indices):
    if namespace not in {"estimate_state","estimate_params","estimate_future"}:
        raise ValueError("Truth stream is unavailable to estimator")
    return rng(namespace,*indices)


def evaluate(packet,context,anchor,model,observation_id):
    if type(packet) is not ObservationPacket or type(context) is not PublicContext:
        raise ValueError("Estimator rejects latent states, true parameters and unknown packet types")
    cfg=PROTOCOL["synthetic"]
    probabilities,health,reserves=[],[],[]
    scenario=dict(hazard_factor=1,growth_factor=1,recruitment_factor=1,hemlock_extra_hazard=0)
    for draw in range(cfg["estimator_draws"]):
        state=estimate_state(packet,anchor,estimator_rng("estimate_state",observation_id,draw))
        params=public_parameters(model,context,estimator_rng("estimate_params",observation_id,draw))
        bank,lawful=demographic.simulate(state,params,scenario,cfg["horizon"],cfg["estimator_paths"],
                                        estimator_rng("estimate_future",observation_id,draw))
        response=demographic.response(bank,lawful,model,cfg["horizon"],cfg["alpha"],cfg["beta"],
                                      cfg["tau"],cfg["reserve_ratio"])
        p=float(response["viable"].mean())
        r=float(response["reserve_viable"].mean())
        probabilities.append(p)
        reserves.append(r)
        health.append(bool(response["present"] and p>=cfg["theta"] and r>=cfg["theta"]))
    return dict(mass=float(np.mean(probabilities)),
                mass_lower=float(np.quantile(probabilities,.05)),
                mass_upper=float(np.quantile(probabilities,.95)),
                reserve_mass=float(np.mean(reserves)),
                health_probability=float(np.mean(health)),
                interpretation="Fraction of observation-conditioned design draws; not real ecological posterior")
