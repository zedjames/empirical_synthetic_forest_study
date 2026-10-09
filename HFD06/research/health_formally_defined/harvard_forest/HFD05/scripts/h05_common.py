"""Additive HFD05 paths, exact predecessor guards and independent streams."""
import csv,hashlib,importlib.util,json,subprocess,sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[3]
FOUR=ROOT.parent/"HFD04"
THREE=ROOT.parent/"HFD03"
OLD=ROOT.parent/"HFD02S"
sys.path.insert(0,str(FOUR/"scripts"))
import h04_common as original
for module,path in [(original,FOUR/"scripts/h04_common.py"),(original.demography,OLD/"models/demography.py"),(original.empirical,OLD/"scripts/empirical.py")]:
    if Path(module.__file__).resolve()!=path.resolve():raise ValueError("Historical import resolution changed")
BASELINE="2b217b36347dab531f3c5e3bd5f85a3e65d0fe67"
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read_csv(path):
    with Path(path).open(newline="") as handle:return list(csv.DictReader(handle))
def output_path(name):
    path=(ROOT/name).resolve()
    if ROOT.resolve() not in path.parents:raise ValueError("Output outside additive study")
    return path
def write_csv(name,rows,fields=None):
    rows=list(rows);path=output_path(name);path.parent.mkdir(parents=True,exist_ok=True)
    if not rows and fields is None:raise ValueError("Empty table without declared schema")
    with path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields or list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
def write_json(name,value):
    path=output_path(name);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+"\n")
def config(name="protocol.json"):return json.loads((ROOT/"config"/name).read_text())
def rng(domain,*indices):
    c=config()
    if domain not in c["random_namespaces"] or any(type(i)is not int or i<0 for i in indices):raise ValueError("Unregistered random identity")
    data=json.dumps([c["master_seed"],"HFD05",domain,list(indices)],separators=(",",":")).encode()
    return np.random.Generator(np.random.PCG64(int.from_bytes(hashlib.sha256(data).digest()[:8],"big")))
def predecessor_guard():
    prefixes=["research/health_formally_defined/harvard_forest/"+s+"/" for s in ["HFD01","HFD02S","HFD03","HFD04"]]
    prefixes += ["scripts/verify_hfd"+s+"_harvard_forest.py" for s in ["01","02s","03","04"]]
    entries=subprocess.check_output(["git","ls-tree","-r",BASELINE,"--",*prefixes],cwd=REPO,text=True).splitlines();result={}
    for row in entries:
        meta,name=row.split("\t",1);expected=meta.split()[2]
        actual=subprocess.check_output(["git","hash-object","--",name],cwd=REPO,text=True).strip()
        if actual!=expected:raise ValueError("Frozen predecessor mutation "+name)
        result[name]=expected
    return result
def guard(phase):
    record=json.loads((ROOT/"provenance"/(phase+"_registration.json")).read_text())
    for name,digest in record["sources"].items():
        if sha(ROOT/name)!=digest:raise ValueError("Registered phase changed "+name)
    predecessor_guard()
    if np.__version__!="2.0.2":raise ValueError("Unregistered runtime")
def register(phase,names,outcomes):
    if (ROOT/"provenance"/(phase+"_registration.json")).exists():raise ValueError("Phase already registered")
    if any((ROOT/name).exists() for name in outcomes):raise ValueError("Cannot preregister after phase outcomes")
    record=dict(phase=phase,status="REGISTERED_BEFORE_PHASE_OUTCOMES",sources={name:sha(ROOT/name) for name in names},
        predecessor_baseline=BASELINE,protected_sources=predecessor_guard(),outcome_selection=False,
        runtime=dict(python="3.9.6",numpy=np.__version__),execution_start=config()["actual_execution_start"])
    write_json("provenance/"+phase+"_registration.json",record)
