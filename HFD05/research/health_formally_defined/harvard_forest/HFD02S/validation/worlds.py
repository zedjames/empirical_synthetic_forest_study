"""Synthetic truth harness; never imported by estimator.py."""
from dataclasses import dataclass
import numpy as np
from contracts import load
from demography import parameters, simulate, response
from estimator import ObservationPacket, estimate


@dataclass
class SyntheticWorld:
    state: dict
    parameters: dict
    variant: str
    label: str


def public_anchor(core, model):
    s = core["stats"]
    dt = np.divide(s["unknown_dt"],s["unknown"],out=np.zeros(64),where=s["unknown"]>0)
    h = np.sum(model["hazard_posterior"]*model["hazard_grid"],axis=1)
    n = np.maximum(1,s["n1"]+s["alive_missing_size"]+s["unknown"]*np.exp(-h*dt))
    d = np.clip(model["diameter1"],1,300)
    p = np.divide(s["n1"],s["n1"]+s["unknown"]+s["alive_missing_size"],
                  out=np.zeros(64),where=(s["n1"]+s["unknown"]+s["alive_missing_size"])>0)
    p[(np.arange(64)//4)%4==3] = load("config/design.json")["reconstruction"]["synthetic_wet_observation_probability"]
    return dict(n=n,d=d,p=p)


def generate_world(anchor, model, seeds, index):
    cfg = load("config/design.json")
    rng = seeds.rng("world",index)
    shape = cfg["reconstruction"]["positive_count_prior_shape"]
    n = rng.poisson(rng.gamma(shape,anchor["n"]/shape))
    d = np.clip(rng.normal(anchor["d"],cfg["reconstruction"]["synthetic_dbh_prior_sd_cm"]),1,300)
    return SyntheticWorld(dict(n=n,d=d),parameters(model,rng),
                          cfg["ensemble"]["model_variants"][index%2],
                          "HF-SYNTH-REF-001" if index==0 else "HF-SYNTH-WORLD-%03d"%(index+1))


def observe(world, anchor, rng):
    n = rng.binomial(world.state["n"],anchor["p"])
    sd = load("config/design.json")["reconstruction"]["synthetic_measurement_sd_cm"]
    d = np.clip(rng.normal(world.state["d"],sd),1,300)
    return ObservationPacket(tuple(int(x) for x in n),
                             tuple(float(x) if count else None for x,count in zip(d,n)),
                             tuple(float(x) for x in anchor["p"]))


def validate_worlds(core, model, seeds, progress=lambda x:None):
    cfg = load("config/design.json")
    v = cfg["validation"]
    anchor = public_anchor(core,model)
    scenarios = [s for s in cfg["scenarios"] if s["id"] in v["scenario_ids"]]
    rows, reference, packets = [], None, []
    for index in range(cfg["ensemble"]["worlds"]):
        world = generate_world(anchor,model,seeds,index)
        packet = observe(world,anchor,seeds.rng("world_observation",index))
        predictions = estimate(packet,anchor,model,seeds,index)
        packets.append(packet)
        truths = []
        reference_banks = []
        for j,scenario in enumerate(scenarios):
            paths = cfg["ensemble"]["reference_future_draws" if index==0 else "world_future_draws"]
            bank,lawful = simulate(world.state,world.parameters,scenario,v["horizon_years"],paths,
                                   seeds.rng("world_truth_future",index,j),world.variant)
            rr = response(bank,lawful,model,v["horizon_years"],v["alpha"],v["beta"],v["tau"],v["reserve_ratio"])
            true_p = float(rr["viable"].mean())
            true_rp = float(rr["reserve_viable"].mean())
            true_health = bool(rr["present"] and true_p>=v["theta"] and true_rp>=v["theta"])
            true_ba,true_juv = float(bank[0,0,1]),float(bank[0,0,2])
            truths.append(dict(scenario=scenario["id"],lawful_mass=float(lawful.mean()),viable_mass=true_p,
                               reserve_mass=true_rp,present=rr["present"],health=true_health,finite_bank_size=paths))
            if index==0:
                reference_banks.append(dict(scenario=scenario["id"],bank=bank,lawful=lawful))
            predictions_p = np.array(predictions["probability"][j])
            pred_health = float(np.mean(predictions["health"][j]))
            low,high = np.quantile(predictions_p,[.05,.95])
            ba_low,ba_high = np.quantile(predictions["basal_area"],[.05,.95])
            j_low,j_high = np.quantile(predictions["juveniles"],[.05,.95])
            rows.append(dict(world=world.label,index=index,scenario=scenario["id"],true_ba=true_ba,
                             estimated_ba=float(np.mean(predictions["basal_area"])),
                             ba_low=float(ba_low),ba_high=float(ba_high),ba_covered=bool(ba_low<=true_ba<=ba_high),
                             true_juveniles=true_juv,estimated_juveniles=float(np.mean(predictions["juveniles"])),
                             juvenile_covered=bool(j_low<=true_juv<=j_high),
                             true_viable_mass=true_p,estimated_viable_mass=float(predictions_p.mean()),
                             viable_low=float(low),viable_high=float(high),viable_covered=bool(low<=true_p<=high),
                             true_health=true_health,estimated_health_fraction=pred_health,
                             brier=(pred_health-int(true_health))**2,
                             failure_provenance="SIMULATED-world-and-observation / IMPUTED-estimator"))
        if index==0:
            reference = dict(label=world.label,selection="FIXED_INDEX_ZERO_NOT_SELECTED_BY_HEALTH",
                             latent_counts=world.state["n"].tolist(),latent_diameters_cm=world.state["d"].tolist(),
                             parameters={k:x.tolist() for k,x in world.parameters.items()},variant=world.variant,
                             synthetic_observations=dict(counts=packet.counts,diameters_cm=packet.diameters_cm,
                                                         observation_probabilities=packet.observation_probabilities),
                             truth=truths,estimation=predictions,history_banks=reference_banks)
        if (index+1)%4==0:
            progress("synthetic worlds %d/%d"%(index+1,cfg["ensemble"]["worlds"]))
    return rows, reference, anchor, packets
