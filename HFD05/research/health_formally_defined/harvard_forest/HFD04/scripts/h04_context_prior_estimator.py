"""Unknown-context estimator variant; no truth generator imports."""
import hashlib,json
import numpy as np
from h04_common import ObservationPacket,estimate_state,parameters,demography,response,classify
from h04_estimator import Context
def generator(domain,identity,draw):
    if domain not in {"cw_state","cw_params","cw_future"}:raise ValueError("Unregistered context-prior stream")
    payload=json.dumps([202610070404,"HFD04-CONTEXT-PRIOR-001",domain,identity,draw],separators=(",",":")).encode()
    return np.random.Generator(np.random.PCG64(int.from_bytes(hashlib.sha256(payload).digest()[:8],"big")))
def evaluate(packet,marker,weights,anchor,model,identity):
    if type(packet) is not ObservationPacket or type(marker)is not Context or marker.policy!="B2_unknown" or marker.pressure is not None or marker.recovery is not None:
        raise ValueError("Only observed packet and unknown marker allowed")
    if len(weights)!=9 or any(w<0 for w in weights) or not np.isclose(sum(weights),1,atol=1e-15):
        raise ValueError("Invalid disclosed experimental-design prior")
    rows=[]
    for draw in range(32):
        state=estimate_state(packet,anchor,generator("cw_state",identity,draw))
        g=generator("cw_params",identity,draw);level=int(g.choice(9,p=weights))
        params=parameters(model,level//3,level%3,g)
        bank,lawful=demography.simulate(state,params,dict(hazard_factor=1,growth_factor=1,recruitment_factor=1,hemlock_extra_hazard=0),5,256,generator("cw_future",identity,draw))
        rr=response(bank,lawful,model);nv=int(rr["viable"].sum());nr=int(rr["reserve_viable"].sum())
        q1=classify(rr["present"],[nv],256,.75);q2=classify(rr["present"],[nv,nr],256,.75)
        rows.append(dict(draw=draw,nv=nv,nr=nr,K=256,present=rr["present"],Q1=q1["FINITE_BANK_HEALTH"],Q2=q2["FINITE_BANK_HEALTH"],
             Q1_status=q1["KERNEL_MC_STATUS"],Q2_status=q2["KERNEL_MC_STATUS"],sampled_context_level=level,
             sampled_context_provenance="DESIGN_PRIOR_DRAW_NOT_HIDDEN_REALIZED_LABEL"))
    return rows
