"""HFD04 private outputs, immutable predecessor interfaces and seed domains."""
import csv, hashlib, importlib.util, json, math, subprocess, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[3]
OLD=ROOT.parent/"HFD02S"
THREE=ROOT.parent/"HFD03"
for directory in [THREE/"scripts",OLD/"scripts",OLD/"models",OLD/"validation"]:
    sys.path.append(str(directory))
import contracts as frozen
import empirical
import demography
import h03_common as prior
from estimator import ObservationPacket, estimate_state
PILOT=json.loads((ROOT/"config/pilot.json").read_text())
DOMAINS={"pilot_Q1","pilot_Q2","pilot_v2_Q1","pilot_v2_Q2","world_parameters","world_state","observation","truth","estimate_state","estimate_params","estimate_future","measurement_error","representation_state","representation_future","cap_state","cap_future","cap_calibration"}
def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def rng(domain,*indices):
    if domain not in DOMAINS or any(type(i) is not int or i<0 for i in indices):
        raise ValueError("Unregistered random identity")
    payload=json.dumps([PILOT["master_seed"],"HFD04",domain,list(indices)],separators=(",",":")).encode()
    return np.random.Generator(np.random.PCG64(int.from_bytes(hashlib.sha256(payload).digest()[:8],"big")))
def config(): return json.loads((ROOT/"config/protocol.json").read_text())
def run_dir():
    path=ROOT/".runs"/sha(ROOT/"config/protocol.json")[:16]
    path.mkdir(parents=True,exist_ok=True)
    return path
def write_json(name,value):
    path=ROOT/name;path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n")
def write_csv(name,rows):
    rows=list(rows)
    if not rows: raise ValueError("Empty table "+name)
    path=ROOT/name;path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as handle:
        out=csv.DictWriter(handle,fieldnames=list(rows[0]),lineterminator="\n")
        out.writeheader();out.writerows(rows)
def read_csv(path):
    with Path(path).open(newline="") as handle: return list(csv.DictReader(handle))
def source_guard():
    paths=["research/health_formally_defined/harvard_forest/"+p+"/" for p in ["HFD01","HFD02S","HFD03"]]
    paths+=["scripts/verify_hfd"+p+"_harvard_forest.py" for p in ["01","02s","03"]]
    entries=subprocess.check_output(["git","ls-tree","-r","0f28602d910d9a8841bf397e9ec842778a13e3f1","--",*paths],cwd=REPO,text=True).splitlines()
    result={}
    for line in entries:
        meta,name=line.split("\t",1);expected=meta.split()[2]
        actual=subprocess.check_output(["git","hash-object","--",name],cwd=REPO,text=True).strip()
        if actual!=expected: raise ValueError("Frozen scientific mutation "+name)
        result[name]=expected
    return result
def load_data():
    sources,e0,e1,trees=prior.inputs()
    core=empirical.read_core(sources);model=empirical.fit(core)
    return sources,e0,e1,core,model
def confirm_guard():
    record=json.loads((ROOT/"provenance/registration.json").read_text())
    for group in ["scientific_sources","configs"]:
        for name,expected in record[group].items():
            if sha(ROOT/name)!=expected:raise ValueError("Registered implementation/config changed "+name)
    if record["protected_sources"]!=source_guard():raise ValueError("Old scientific identity changed")
    if __import__("numpy").__version__!="2.0.2":raise ValueError("Unregistered NumPy runtime")
def isolated_models(cap):
    """Independent module globals; never mutate original modules or configs."""
    cfg=frozen.load("config/design.json")
    cfg["fit"]["parameter_effective_sample_cap"]=cap
    modules=[]
    for label,path in [("empirical",OLD/"scripts/empirical.py"),("demography",OLD/"models/demography.py")]:
        name="h04_cap_"+str(cap)+"_"+label
        spec=importlib.util.spec_from_file_location(name,path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        def load(name,original=frozen.load,design=cfg):
            return json.loads(json.dumps(design)) if name=="config/design.json" else original(name)
        module.load=load;modules.append(module)
    return modules
def classify(present,counts,K,theta):
    counts=np.asarray(counts,dtype=int)
    if K<=0 or np.any(counts<0) or np.any(counts>K): raise ValueError("Invalid unconditional census")
    lower,upper=prior.wilson(counts,K,.95 if len(counts)==1 else .975)
    finite=bool(present and np.all(counts/K>=theta))
    if not present: status="FALSE"
    elif theta==0: status="TRUE"
    elif np.any(upper<theta): status="FALSE"
    elif np.all(lower>theta): status="TRUE"
    else: status="MC_UNRESOLVED"
    return dict(FINITE_BANK_HEALTH=finite,KERNEL_MC_STATUS=status,
                lower=list(map(float,lower)),upper=list(map(float,upper)))
def parameters(model,p,r,generator,means=False):
    params=demography.parameters(model,generator)
    h=[.005,.03,.12][p]
    entry=[.025,.008,.001][r]*model["baseline_juveniles"]
    shares=model["recruitment_mean"]+1/64;shares=shares/shares.sum()
    params["hazard"]=np.full(64,h) if means else h*np.exp(generator.normal(-.15**2/2,.15,64))
    params["recruit"]=entry*shares if means else entry*shares*np.exp(generator.normal(-.2**2/2,.2,64))
    if means:params["growth"]=model["growth_mean"].copy()
    return params
def make_state(core,model,structure,regeneration):
    n=core["stats"]["n0"].copy();d=model["diameter0"].copy();juvenile=d<10
    n[juvenile]=np.rint(n[juvenile]*regeneration*model["baseline_juveniles"]/max(n[juvenile].sum(),1)).astype(int)
    jba=float((n[juvenile]*d[juvenile]**2).sum()*math.pi/40000)
    aba=float((n[~juvenile]*d[~juvenile]**2).sum()*math.pi/40000)
    d[~juvenile]*=math.sqrt(max(structure*model["baseline_ba"]-jba,1)/max(aba,1))
    return dict(n=n,d=np.clip(d,1,300))
def response(bank,lawful,model):
    c=PILOT
    return demography.response(bank,lawful,model,c["horizon"],c["alpha"],c["beta"],c["tau"],c["reserve_ratio"])
def scenario(suite="correct"):
    return dict(hazard_factor=1,growth_factor=1,recruitment_factor=1,hemlock_extra_hazard=.12 if suite=="hemlock_extra_hazard_0.12" else 0)
class EntryRNG:
    def __init__(self,generator):self.generator=generator
    def __getattr__(self,key):
        if key!="poisson":return getattr(self.generator,key)
        def call(lam,size=None):
            scale=self.generator.lognormal(-1.2**2/2,1.2,size=(size[0],1))
            return self.generator.poisson(np.broadcast_to(lam,size)*scale)
        return call
