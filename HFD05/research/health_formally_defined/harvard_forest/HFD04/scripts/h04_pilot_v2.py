"""Independent second pilot; juvenile diameter controls reserve margin."""
import json, math
import numpy as np
import h04_common as h
def size_boundary_state(core,model,structure,regeneration,juvenile_diameter):
    state=h.make_state(core,model,structure,regeneration)
    if juvenile_diameter is not None:
        j=model["diameter0"]<10
        state["d"][j]=juvenile_diameter
        jb=float((state["n"][j]*state["d"][j]**2).sum()*math.pi/40000)
        ab=float((state["n"][~j]*state["d"][~j]**2).sum()*math.pi/40000)
        state["d"][~j]*=math.sqrt(max(structure*model["baseline_ba"]-jb,1)/max(ab,1))
    return state
def run():
    _,_,_,core,model=h.load_data();c=json.loads((h.ROOT/"config/pilot_v2.json").read_text())
    rows=[];cells=[]
    for family in ["Q1","Q2"]:
        for p in (range(3) if family=="Q1" else [2]):
            for r in range(3):
                params=h.parameters(model,p,r,h.rng("pilot_v2_"+family,p,r),means=True)
                for ti,target in enumerate(h.PILOT["targets"]):
                    low,high=h.PILOT["Q1_structure_bracket"] if family=="Q1" else c["Q2_diameter_bracket_cm"]
                    for step in range(c["bisection_steps"]):
                        value=(low+high)/2
                        s=value if family=="Q1" else c["Q2_structure"]
                        jd=None if family=="Q1" else value
                        state=size_boundary_state(core,model,s,.9,jd)
                        bank,lawful=h.demography.simulate(state,params,h.scenario(),5,c["paths"],h.rng("pilot_v2_"+family,p,r,ti,step))
                        rr=h.response(bank,lawful,model)
                        mass=float(rr["viable" if family=="Q1" else "reserve_viable"].mean())
                        rows.append(dict(family=family,pressure=p,recovery=r,target=target,step=step,structure=s,juvenile_diameter=jd,mass=mass,present=rr["present"]))
                        if (mass<target)==(family=="Q1"):low=value
                        else:high=value
                    cells.append(dict(family=family,pressure=p,recovery=r,target=target,structure=s,regeneration=.9,juvenile_diameter=jd,pilot_mass=mass))
        print("PILOT002",family,"complete",flush=True)
    h.write_csv("boundary_validation/pilot_v2_probes.csv",rows)
    h.write_json("boundary_validation/pilot_v2_design.json",dict(id=c["id"],cells=cells,pilot_only=True,all_probes_retained=True,confirmatory_streams_used=False))
if __name__=="__main__":run()
