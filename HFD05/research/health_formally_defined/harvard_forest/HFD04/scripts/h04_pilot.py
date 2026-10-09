"""Separately registered pilot; no estimator and no confirmatory streams."""
import h04_common as h
import numpy as np
def run():
    _,_,_,core,model=h.load_data();c=h.PILOT
    rows=[];design=[]
    for family in ["Q1","Q2"]:
        pressures=range(3) if family=="Q1" else [c["Q2_pressure"]]
        for p in pressures:
            for r in range(3):
                params=h.parameters(model,p,r,h.rng("pilot_"+family,p,r),means=True)
                for ti,target in enumerate(c["targets"]):
                    low,high=c[family+"_structure_bracket"] if family=="Q1" else c["Q2_regeneration_bracket"]
                    # Q1 mass rises with structure; Q2 reserve fraction decreases as
                    # own starting juvenile support increases. All probes retained.
                    for step in range(c["bisection_steps"]):
                        middle=(low+high)/2
                        structure=middle if family=="Q1" else c["Q2_structure"]
                        regeneration=c["Q1_regeneration"] if family=="Q1" else middle
                        state=h.make_state(core,model,structure,regeneration)
                        bank,lawful=h.demography.simulate(state,params,h.scenario(),c["horizon"],c["paths"],h.rng("pilot_"+family,p,r,ti,step))
                        rr=h.response(bank,lawful,model)
                        mass=float(rr["viable" if family=="Q1" else "reserve_viable"].mean())
                        rows.append(dict(family=family,pressure=p,recovery=r,target=target,step=step,structure=structure,regeneration=regeneration,mass=mass,present=rr["present"]))
                        if (mass<target)==(family=="Q1"):low=middle
                        else:high=middle
                    design.append(dict(family=family,pressure=p,recovery=r,target=target,
                                       structure=middle if family=="Q1" else c["Q2_structure"],
                                       regeneration=c["Q1_regeneration"] if family=="Q1" else middle,
                                       pilot_mass=mass,bracket_low=low,bracket_high=high))
        print("PILOT family",family,flush=True)
    h.write_csv("boundary_validation/pilot_probes.csv",rows)
    h.write_json("boundary_validation/pilot_design.json",dict(id=c["id"],cells=design,pilot_only=True,all_probes_retained=True,confirmatory_streams_used=False))
    print("PILOT complete",len(rows),"probes",len(design),"candidate cells",flush=True)
if __name__=="__main__":run()
