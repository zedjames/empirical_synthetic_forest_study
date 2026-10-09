"""Explicit information interfaces. No world/truth generator import in this module."""
from dataclasses import dataclass
import numpy as np
from h04_common import ObservationPacket,estimate_state,parameters,demography
from h04_estimator import Context

@dataclass(frozen=True)
class EstimatedState:
    counts:tuple
    diameters:tuple
@dataclass(frozen=True)
class OracleState:
    counts:tuple
    diameters:tuple
@dataclass(frozen=True)
class EstimatedParameters:
    hazard:tuple
    growth:tuple
    sd:tuple
    recruit:tuple
@dataclass(frozen=True)
class OracleParameters:
    hazard:tuple
    growth:tuple
    sd:tuple
    recruit:tuple

def pack_state(value,role):return role(tuple(map(int,value["n"])),tuple(map(float,value["d"])))
def pack_parameters(value,role):return role(**{k:tuple(map(float,value[k])) for k in ["hazard","growth","sd","recruit"]})
def estimate(packet,marker,anchor,model,state_rng,param_rng):
    if type(packet)is not ObservationPacket or type(marker)is not Context or marker.policy!="B0_full":raise ValueError("E00 only accepts observed full-context interface")
    state=estimate_state(packet,anchor,state_rng)
    levels=marker.levels();p,r=levels[int(param_rng.integers(len(levels)))]
    params=parameters(model,p,r,param_rng)
    return pack_state(state,EstimatedState),pack_parameters(params,EstimatedParameters)
def predict(condition,state,params,paths,generator):
    roles={"E00":(EstimatedState,EstimatedParameters),"E10":(OracleState,EstimatedParameters),"E01":(EstimatedState,OracleParameters),"E11":(OracleState,OracleParameters)}
    if condition not in roles or (type(state),type(params))!=roles[condition]:raise ValueError("Oracle leakage or mislabeled condition")
    n=np.asarray(state.counts);d=np.asarray(state.diameters);p={k:np.asarray(getattr(params,k)) for k in ["hazard","growth","sd","recruit"]}
    if len(n)!=64 or len(d)!=64 or any(len(v)!=64 for v in p.values()) or np.any(n<0) or np.any(n!=np.floor(n)) or np.any((d<1)|(d>300)) or any(np.any(v<0)|np.any(~np.isfinite(v)) for v in p.values()):raise ValueError("Invalid role payload")
    return demography.simulate(dict(n=n,d=d),p,dict(hazard_factor=1,growth_factor=1,recruitment_factor=1,hemlock_extra_hazard=0),5,paths,generator)
