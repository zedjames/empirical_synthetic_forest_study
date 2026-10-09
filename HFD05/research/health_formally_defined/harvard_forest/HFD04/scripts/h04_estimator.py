"""Only strict observation packets and actually disclosed context labels enter.

No world design, truth implementation or original labels are imported.
"""
from dataclasses import dataclass
import itertools
import numpy as np
from h04_common import ObservationPacket, estimate_state, parameters, demography, rng, response, classify

@dataclass(frozen=True)
class Context:
    policy: str
    pressure: object=None
    recovery: object=None
    def __post_init__(self):
        if self.policy in {"B0_full","B3_label_only"}:
            if type(self.pressure) is not int or type(self.recovery) is not int or self.pressure not in range(3) or self.recovery not in range(3):
                raise ValueError("Full context requires disclosed labels")
        elif self.policy=="B1_coarsened":
            if type(self.pressure) is not int or type(self.recovery) is not int or self.pressure not in range(2) or self.recovery not in range(2):
                raise ValueError("Coarse context may not carry exact labels")
        elif self.policy=="B2_unknown":
            if self.pressure is not None or self.recovery is not None:raise ValueError("Hidden exact context in unknown marker")
        else:raise ValueError("Unknown information condition")
    def levels(self):
        if self.policy=="B1_coarsened":
            return list(itertools.product([[0,1],[2]][self.pressure],[[0],[1,2]][self.recovery]))
        if self.policy=="B2_unknown":return list(itertools.product(range(3),range(3)))
        return [(self.pressure,self.recovery)]

def evaluate(packet,context,anchor,model,identity,draws,paths):
    if type(packet) is not ObservationPacket or type(context) is not Context:
        raise ValueError("Estimator rejects truth and latent packets")
    if context.policy=="B3_label_only":
        if any(packet.counts) or any(p!=0 for p in packet.observation_probabilities):
            raise ValueError("Label-only packet contains state observations")
    levels=context.levels();rows=[]
    for draw in range(draws):
        state=estimate_state(packet,anchor,rng("estimate_state",identity,draw))
        generator=rng("estimate_params",identity,draw)
        # Draw mixing under disclosed design weights, never true level lookup.
        p,r=levels[int(generator.integers(len(levels)))]
        params=parameters(model,p,r,generator)
        bank,lawful=demography.simulate(state,params,dict(hazard_factor=1,growth_factor=1,recruitment_factor=1,hemlock_extra_hazard=0),5,paths,rng("estimate_future",identity,draw))
        rr=response(bank,lawful,model)
        nv=int(rr["viable"].sum());nr=int(rr["reserve_viable"].sum())
        q1=classify(rr["present"],[nv],paths,.75);q2=classify(rr["present"],[nv,nr],paths,.75)
        rows.append(dict(draw=draw,present=rr["present"],nv=nv,nr=nr,K=paths,Q1=q1["FINITE_BANK_HEALTH"],Q2=q2["FINITE_BANK_HEALTH"],Q1_status=q1["KERNEL_MC_STATUS"],Q2_status=q2["KERNEL_MC_STATUS"]))
    return rows
