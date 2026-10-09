"""Small source-contained reproducibility fixture, not confirmatory data."""
import hashlib,json,math
import numpy as np
import h04_common as h
def build():
    _,_,_,core,model=h.load_data()
    serial={k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in model.items()}
    state=dict(n=core["stats"]["n0"].tolist(),d=model["diameter0"].tolist())
    params=h.demography.parameters(model,np.random.default_rng(2704))
    bank,lawful=h.demography.simulate({k:np.array(v) for k,v in state.items()},params,h.scenario(),5,16,np.random.default_rng(2705))
    h.write_json("provenance/lightweight_fixture.json",dict(id="HFD04-LIGHTWEIGHT-001",parameters_seed=2704,futures_seed=2705,K=16,model=serial,state=state,
       bank_sha256=hashlib.sha256(bank.tobytes()+lawful.tobytes()).hexdigest(),purpose="Explicit source-contained derived-model fixture, not extra benchmark worlds"))
    # Actual query witness: equal count and BA/RMS, different juvenile support.
    ba=100*100*math.pi/40000
    models=dict(baseline_ba=ba,baseline_juveniles=100)
    histories=[]
    for J in [0,50]:
        b=np.array([[[100,ba,J,0,0,0],[100,ba,J,0,0,0]]],float)
        r=h.demography.response(b,np.array([True]),models,1,.75,.4,0,.5)
        histories.append(dict(J=J,present=r["present"],Q2=bool(r["present"] and r["reserve_viable"].all())))
    h.write_json("juvenile_boundary/projection_counterexample.json",dict(coarse=dict(n=100,RMS_d=10,BA=ba),
      fine_A=dict(n=[100],d=[10],J=0),fine_B=dict(n=[50,50],d=[9,math.sqrt(119)],J=50),
      exact_second_moment=100,actual_response=histories,
      conclusion="Same RMS/count/BA does not determine Realizes or Q2; explicit numerical witness, not a conservative dynamics theorem"))
    print("Lightweight kernel fixture and query counterexample generated",flush=True)
if __name__=="__main__":build()

