"""Package-only review runner. Operational adapters are declared, never scientific rewrites."""
import argparse,csv,hashlib,importlib.util,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
STUDY=ROOT/'research/health_formally_defined/harvard_forest'
FIVE=STUDY/'HFD05'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def audit():
    manifest=json.loads((ROOT/'HFD05_Artifact_Manifest.json').read_text())
    expected=set(manifest['files'])|{'HFD05_Artifact_Manifest.json'}
    runtime_outputs=set(manifest.get('runtime_output_allowlist',[]))
    permitted_runtime={'research/health_formally_defined/harvard_forest/HFD05/corrected_pipeline/'+name for name in ['reference_counts_partial.csv','context_draws_partial.csv','cap_surface_partial.csv','local_precision_partial.csv']}
    if not runtime_outputs.issubset(permitted_runtime):raise ValueError('Out-of-scope runtime output permission')
    for name,record in manifest['files'].items():
        path=ROOT/name
        if path.is_symlink() or not path.is_file() or sha(path)!=record['sha256']:raise ValueError('Artifact mutation '+name)
    for path in ROOT.rglob('*'):
        if not path.is_file():continue
        name=str(path.relative_to(ROOT))
        if name in expected or name in runtime_outputs:continue
        if any(part in ['.inputs','.cache','.runs','__pycache__'] for part in path.relative_to(ROOT).parts):continue
        raise ValueError('Unlisted export '+name)
    return manifest
def module():
    manifest=audit()
    spec=importlib.util.spec_from_file_location('original_public_runner',ROOT/'artifact_runner.py');legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy);legacy.module()
    sys.path.insert(0,str(FIVE/'scripts'));import h05_common as h
    def own_predecessors():
        # Local package integrity is verified; absent private ancestry is an
        # attested ledger, NOT a reproduced private Git-history check.
        audit();return json.loads((FIVE/'provenance/q2_registration.json').read_text())['protected_sources']
    def phase_guard(phase):
        audit();record=json.loads((FIVE/'provenance'/(phase+'_registration.json')).read_text())
        for name,digest in record['sources'].items():
            if sha(FIVE/name)!=digest:raise ValueError('Preregistered public phase bytes differ '+name)
        if __import__('numpy').__version__!='2.0.2':raise ValueError('Unregistered NumPy')
    h.predecessor_guard=own_predecessors;h.guard=phase_guard
    for name in ['h05_common','h04_common','h03_common','contracts','empirical','demography','estimator']:
        imported=__import__(name)
        if not Path(imported.__file__).resolve().is_relative_to(ROOT.resolve()):raise ValueError('Private import '+name)
    return h
def tier1():
    h=module();import numpy as np
    from h05_inference import classify
    from h05_growth import laws,sample
    from fractions import Fraction
    evaluated=0
    for K in range(1,25):
        for nv in range(K+1):
            for nr in range(nv+1):
                for theta in [0,.5,.75,.9,1]:
                    for present in [False,True]:
                        d=classify(present,nv,nr,K,theta,'Q2');expected=present and Fraction(nr,K)>=Fraction(str(theta));assert d['FINITE_BANK_HEALTH']==expected;evaluated+=1
    mu=np.array([0,1e-10,.02,.7]);sd=np.array([.1,.1,.1,0]);law=laws(mu,sd,True);assert np.allclose(law['mean'],mu,rtol=1e-14,atol=0)
    x=sample(mu,sd,(16,4),h.rng('fixture',0),True);assert np.all(x[:,0]==0) and np.all(x[:,3]==.7)
    # An actual original/revised generator fixture, using the exact included model.
    data=json.loads((STUDY/'HFD04/provenance/lightweight_fixture.json').read_text());model={k:np.array(v) if isinstance(v,list) else v for k,v in data['model'].items()};state={k:np.array(v) for k,v in data['state'].items()};p=h.original.demography.parameters(model,np.random.default_rng(data['parameters_seed']))
    from h05_growth import simulate
    outputs=[]
    for revised in [False,True]:
        bank,lawful=simulate(state,p,h.original.scenario(),5,data['K'],np.random.default_rng(data['futures_seed']),'gamma_growth',revised);rr=h.original.response(bank,lawful,model);assert np.all(~rr['reserve_viable']|rr['viable']);outputs.append(dict(corrected=revised,nv=int(rr['viable'].sum()),nr=int(rr['reserve_viable'].sum()),bank_sha256=hashlib.sha256(bank.tobytes()+lawful.tobytes()).hexdigest()))
    print(json.dumps(dict(status='PASS',tier=1,exhaustive_queries=evaluated,generator_fixture=outputs,private_inputs=False)))
def tier2():
    h=module();from collections import Counter,defaultdict
    def rows(name):
        with (FIVE/name).open(newline='') as handle:return list(csv.DictReader(handle))
    primary=rows('corrected_pipeline/primary_surface.csv');missing=rows('corrected_pipeline/missingness_surface.csv')
    assert Counter(r['result'] for r in primary)==dict.fromkeys(['R0','R1','R2','R3','R4'],6480)
    assert Counter(r['result'] for r in missing)==dict.fromkeys(['R0','R1','R2','R3','R4'],24300)
    for table in [primary,missing]:
        for r in table:
            assert float(r['reserve_mass'])<=float(r['viable_mass'])
            assert sum(int(r[k]) for k in ['kernel_TRUE_units','kernel_FALSE_units','kernel_MC_UNRESOLVED_units'])==int(r['outer_units'])
        for result in ['R1','R2']:assert all(float(r['finite_health_disagreement_vs_original'])==0 for r in table if r['result']==result)
    cross=rows('corrected_pipeline/reference_prediction_crossmatrix.csv');groups=defaultdict(set)
    for r in cross:groups[(r['cohort'],r['world'],r['suite'],r['reference_K'],r['query'],r['policy'])].add(r['crossmatrix'])
    assert all(len(v)==4 for v in groups.values())
    adult=rows('external_sources/HF453/mortality_scores.csv');target=next(r for r in adult if (r['target'],r['end_year'],r['scheme'],r['axis'])==('CUMULATIVE','2024','STRICT','ALL'));assert int(target['n'])==21576 and int(target['observed_deaths'])==1477
    assert len(rows('external_sources/HF355/transition_records.csv'))==8840
    for i in range(1,11):
        metadata=json.loads((FIVE/'figures'/('H05-'+str(i)+'.json')).read_text());assert sha(FIVE/metadata['data'])==metadata['data_sha256']
        for name,digest in metadata['sources'].items():assert sha(FIVE/name)==digest
    print(json.dumps(dict(status='PASS',tier=2,primary_cells=len(primary),missingness_cells=len(missing),crossmatrix_cells=len(cross),ten_figures=True,private_inputs=False)))
def retrieve():
    audit();inventory=json.loads((ROOT/'Public_Inputs.json').read_text())['sources'];external=json.loads((FIVE/'external_sources/source_inventory.json').read_text())['sources']
    entries=inventory+[dict(id=r['name'],target=str(Path('research/health_formally_defined/harvard_forest/HFD05')/r['local_input']),url=r['url'],sha256=r['sha256']) for r in external]
    for entry in entries:
        path=(ROOT/entry['target']).resolve()
        if not path.is_relative_to(ROOT.resolve()):raise ValueError('Unisolated public input')
        if path.exists():
            if sha(path)!=entry['sha256']:raise ValueError('Existing input mismatch, refusing overwrite')
            continue
        path.parent.mkdir(parents=True,exist_ok=True)
        subprocess.run(['curl','--fail','--location','--silent','--show-error',entry['url'],'--output',str(path)],check=True)
        if sha(path)!=entry['sha256']:raise ValueError('Publisher bytes changed, not accepted')
    print(json.dumps(dict(status='PASS',public_inputs=len(entries),only_exact_public_CC0_sources=True)))
def tier3(stage):
    h=module()
    if stage=='external-audit':
        import h05_external_audit as a;a.audit()
    elif stage=='external-scoring':
        import h05_external_scores as a;a.adult();a.seedlings()
    elif stage=='oracle':
        import h05_oracle as a;a.generate();a.analyze()
    elif stage=='oracle-analysis':
        # Recovery of a completely generated, byte-verified frozen census after
        # the original runtime-only partial-file allowlist omission. audit()
        # already verifies every included census byte against its frozen hash.
        import h05_oracle as a;a.analyze()
    elif stage=='references':
        import h05_boundary_replay as a;a.references()
    elif stage=='context':
        import h05_boundary_replay as a;a.predictions();a.analyze()
    elif stage=='figures':
        import h05_figures as a;a.run()
    else:raise ValueError('Select an explicit supported stage')
    print(json.dumps(dict(status='PASS',tier=3,stage=stage,private_raw_or_git_used=False)))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['audit','tier1','tier2','retrieve','tier3']);p.add_argument('--stage');a=p.parse_args();start=time.monotonic()
    if a.mode=='audit':print(json.dumps(dict(status='PASS',files=len(audit()['files']))))
    elif a.mode=='tier1':tier1()
    elif a.mode=='tier2':tier2()
    elif a.mode=='retrieve':retrieve()
    else:tier3(a.stage)
    print(json.dumps(dict(observed_runtime_seconds=time.monotonic()-start)))
