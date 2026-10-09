"""Observation-only estimator. This module never imports the reference world."""
from dataclasses import dataclass
import numpy as np
from contracts import load
from demography import parameters, simulate, response


@dataclass(frozen=True)
class ObservationPacket:
    counts: tuple
    diameters_cm: tuple
    observation_probabilities: tuple

    def __post_init__(self):
        if not (len(self.counts)==len(self.diameters_cm)==len(self.observation_probabilities)==64):
            raise ValueError("Observation packet shape")
        for n,d,p in zip(self.counts,self.diameters_cm,self.observation_probabilities):
            if not isinstance(n,int) or n<0 or not 0<=p<=1 or (p==0 and n):
                raise ValueError("Invalid observation operator")
            if (n==0 and d is not None) or (n>0 and (d is None or not 1<=d<=300)):
                raise ValueError("Diameter observation/missingness mismatch")


def estimate_state(packet, anchor, rng):
    if type(packet) is not ObservationPacket:
        raise ValueError("Estimator accepts only the strict observation packet")
    cfg = load("config/design.json")["reconstruction"]
    y = np.array(packet.counts)
    p = np.array(packet.observation_probabilities)
    shape = cfg["positive_count_prior_shape"]
    rate = shape/np.maximum(anchor["n"],1)
    intensity = rng.gamma(shape+y,1/(rate+p))
    n = y+rng.poisson((1-p)*intensity)
    prior_var = cfg["synthetic_dbh_prior_sd_cm"]**2
    obs_var = cfg["synthetic_measurement_sd_cm"]**2
    observed = np.array([v is not None for v in packet.diameters_cm])
    obs = np.array([v if v is not None else 0 for v in packet.diameters_cm])
    precision = 1/prior_var+observed/obs_var
    mean = (anchor["d"]/prior_var+observed*obs/obs_var)/precision
    d = np.clip(rng.normal(mean,np.sqrt(1/precision)),1,300)
    return dict(n=n,d=d)


def estimate(packet, anchor, model, seeds, world_index):
    cfg = load("config/design.json")
    v = cfg["validation"]
    draws = cfg["ensemble"]["validation_state_draws"]
    scenarios = [s for s in cfg["scenarios"] if s["id"] in v["scenario_ids"]]
    result = dict(basal_area=[],juveniles=[],probability=[[],[]],health=[[],[]],
                  reserve_probability=[[],[]])
    for i in range(draws):
        # All ablation levels are projections of this same draw, index 2.
        rng = seeds.rng("estimate_state",world_index,2,i)
        state = estimate_state(packet,anchor,rng)
        params = parameters(model,rng)
        variant = cfg["ensemble"]["model_variants"][i%2]
        for j,scenario in enumerate(scenarios):
            bank,lawful = simulate(state,params,scenario,v["horizon_years"],
                                   cfg["ensemble"]["validation_future_draws"],
                                   seeds.rng("estimate_future",world_index,2,i,j),variant)
            rr = response(bank,lawful,model,v["horizon_years"],v["alpha"],v["beta"],v["tau"],v["reserve_ratio"])
            p = float(rr["viable"].mean())
            rp = float(rr["reserve_viable"].mean())
            result["probability"][j].append(p)
            result["reserve_probability"][j].append(rp)
            result["health"][j].append(bool(rr["present"] and p>=v["theta"] and rp>=v["theta"]))
            if j==0:
                result["basal_area"].append(float(bank[0,0,1]))
                result["juveniles"].append(float(bank[0,0,2]))
    return result
