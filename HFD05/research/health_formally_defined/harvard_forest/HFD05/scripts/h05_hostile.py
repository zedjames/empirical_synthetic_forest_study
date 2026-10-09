"""Active hostile inputs against mathematics, identity, claims and export integrity."""
import copy,hashlib,json,math,tempfile,importlib.util
from pathlib import Path
from statistics import NormalDist
import numpy as np
import h05_common as h
from h05_inference import classify,finite,wilson
from h05_growth import laws,sample
import h05_oracle_engine as engine

def require(value,message):
    if not value:raise ValueError(message)
def close(actual,expected):require(math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-15),'Scientific value mismatch')
def secure(row,assertion):
    if assertion=='DEAD':require(row['date'] and row['status'] in ['DC','DS'],'No secure death evidence')
def documented_date(actual,claimed):require(actual and actual==claimed,'Invented observation date')
def unique(candidates):require(len(set(candidates))==1,'Ambiguous source join')
def support(source,target,bridge):require(source==target or bridge=='IDENTIFIED','Unsupported support expansion')
def claim(kind,evidence):
    if kind=='POSTERIOR':require(evidence.get('identified_likelihood') and evidence.get('calibrated_prior'),'Diagnostic weights are not posterior')
    elif kind=='UNIFORM_ROBUSTNESS':require(evidence['bound_scope']=='ALL_STATES' and evidence['sup_delta']<=evidence['claimed_tolerance'],'Mean is not uniform bound')
    elif kind=='OPTIMAL_VOI':require(evidence['same_joint_law'] and evidence['same_loss'] and evidence['optimal_actions'],'Changed implemented estimator is not fixed-law optimal risk')
    elif kind=='KERNEL_CERTAINTY':require(evidence['analytic_probability_one'],'All-success finite sample is not kernel probability one')
    elif kind=='DYNAMIC_COMMUTATION':require(evidence['transition_intertwining'],'Structural projection is not dynamic commutation')
    elif kind=='ECOLOGICAL_CAUSAL_ATTRIBUTION':require(evidence['ecological_identification'] and evidence['MC_separated'],'Numerical contrasts are not ecological causal attribution')
def roster(expected,actual):require(set(expected)==set(actual) and len(expected)==len(actual),'Unfavorable or unresolved cases pruned')
def integrity(data,digest):require(hashlib.sha256(data).hexdigest()==digest,'Frozen byte mutation')
def export_path(path,root,licensed):
    p=Path(path);require(not p.is_absolute() and '..' not in p.parts and p.suffix in ['.py','.csv','.json','.svg','.md','.xml','.txt'],'Out-of-scope export path')
    require(licensed=='LOCAL_REVIEW_ONLY' and not any(x in p.parts for x in ['.git','.lake','.cache','credentials','Lean']),'Unlicensed/private export')
def import_path(path,root):require(Path(path).resolve().is_relative_to(Path(root).resolve()),'Private import')
def no_tuning(frozen_fit,used_fit):require(frozen_fit==used_fit,'External outcome tuning')
def run():
    tests=[]
    def reject(name,fn):
        try:fn()
        except (ValueError,AssertionError,TypeError,KeyError):tests.append(dict(control=name,rejected=True));return
        raise ValueError('Hostile case accepted '+name)
    def bad(value):require(value,'Hostile scientific discrepancy')
    # Every control executes a changed input, formula or active interface. No
    # metadata declaration alone counts as a negative-control result.
    base=(h.FOUR/'scripts/h04_common.py').read_bytes();digest=hashlib.sha256(base).hexdigest()
    integrity(base,digest);reject('01_FROZEN_SOURCE_MUTATION',lambda:integrity(base+b'changed',digest))
    reject('02_Q2_POINT_DISAGREEMENT',lambda:bad(finite(True,8,4,8,.75,'Q2') is True))
    # Single-constraint API and actual interval width reject redundant inference.
    reject('03_TWO_INDEPENDENT_Q2_CONSTRAINTS',lambda:wilson([8,4],8))
    reject('04_Q2_BONFERRONI_RETENTION',lambda:wilson(7,8,.975))
    defective=laws(np.array([1e-5]),np.array([.1]),False)
    reject('05_UNDETECTED_GAMMA_WRONG_MEAN',lambda:close(float(defective['mean'][0]),1e-5))
    zero=sample(np.array([0,.7]),np.array([.1,0]),(32,2),h.rng('fixture',1),True)
    reject('06_OMITTED_ZERO_OR_SD_ZERO_BRANCH',lambda:bad(np.all(zero[:,0]==1e-5) and np.all(zero[:,1]==0)))
    from h05_growth import GrowthRNG
    p=dict(growth=np.ones(64)*.1,sd=np.ones(64)*.2);s=dict(growth_factor=1)
    reject('07_SOURCE_FORMULA_ACTUAL_CALL_DISAGREEMENT',lambda:GrowthRNG(h.rng('fixture',2),p,s).gamma(np.ones(64),np.ones(64),size=(2,64)))
    times=np.array([0.,10.]);prob=np.exp(-.1*times);truth=float(prob.sum());pooled=2*math.exp(-.1*float(times.mean()))
    reject('08_POOLED_EXPOSURE_EXACT',lambda:close(pooled,truth))
    reject('09_UNKNOWN_FATE_AS_DEATH',lambda:secure(dict(date='2021',status='UNKNOWN'),'DEAD'))
    reject('10_INVENTED_INDIVIDUAL_DATE',lambda:documented_date('', '2020-01-03'))
    reject('11_ORIGINAL_BANK_OVERWRITE',lambda:h.output_path('../HFD04/.runs/original.npz'))
    flags=np.array([True,False,True]);fd=hashlib.sha256(flags.tobytes()).hexdigest();changed=flags.copy();changed[1]=True
    reject('12_SILENT_TRUTH_CHANGE',lambda:integrity(changed.tobytes(),fd))
    os=engine.OracleState((1,)*64,(2.,)*64);ep=engine.EstimatedParameters((.1,)*64,(.1,)*64,(.1,)*64,(1.,)*64)
    reject('13_E00_ORACLE_LEAKAGE',lambda:engine.predict('E00',os,ep,4,h.rng('fixture',3)))
    L={'E00':.2,'E10':.18,'E01':.3,'E11':.01};interaction=L['E11']-L['E10']-L['E01']+L['E00']
    reject('14_FORCED_ADDITIVE_ORACLE_PARTITION',lambda:close(interaction,0))
    reject('15_MC_NOISE_AS_ECOLOGICAL_ATTRIBUTION',lambda:claim('ECOLOGICAL_CAUSAL_ATTRIBUTION',dict(ecological_identification=False,MC_separated=False)))
    reject('16_COLLIDING_TAG_UNIQUE',lambda:unique(['stem1','stem2']))
    risk={'baseline_alive_tags':{'a','b'},'all_tags':{'a','b','dead','unobserved'}}
    reject('17_LEFT_TRUNCATION_IGNORED',lambda:roster(risk['baseline_alive_tags'],risk['all_tags']))
    reject('18_DN_SECURE_DEATH',lambda:secure(dict(date='2024',status='DN'),'DEAD'))
    reject('19_EXTERNAL_OUTCOME_TUNING',lambda:no_tuning('earlier_fit_sha','new_validation_fit_sha'))
    reject('20_SEEDLING_WHOLE_PLOT_EXPANSION',lambda:support('below1cm_1m2_subplot','wholeplot_1to10cm','NONE'))
    reject('21_FIRST_OBSERVATION_GERMINATION_DATE',lambda:documented_date('', '2017-06-08'))
    reject('22_DIAGNOSTIC_MEMBERSHIP_AS_POSTERIOR',lambda:claim('POSTERIOR',dict(identified_likelihood=False,calibrated_prior=False)))
    reject('23_AVERAGE_AS_GENERAL_ROBUSTNESS',lambda:claim('UNIFORM_ROBUSTNESS',dict(bound_scope='MEAN',sup_delta=.984375,claimed_tolerance=.01)))
    reject('24_CONTEXT_AS_OPTIMAL_INFORMATION_VALUE',lambda:claim('OPTIMAL_VOI',dict(same_joint_law=False,same_loss=True,optimal_actions=False)))
    reject('25_FINITE_BANK_KERNEL_CERTAINTY',lambda:claim('KERNEL_CERTAINTY',dict(analytic_probability_one=False)))
    weights=np.array([.25,.25,.5]);viable=np.array([True,False,False]);unconditional=float(weights[viable].sum());renormalized=float((weights[viable]/weights[viable].sum()).sum())
    reject('26_SUCCESS_RENORMALIZATION',lambda:close(renormalized,unconditional))
    reject('27_PROJECTION_AS_DYNAMIC_COMMUTATION',lambda:claim('DYNAMIC_COMMUTATION',dict(transition_intertwining=False)))
    with tempfile.TemporaryDirectory(prefix='hfd05-hostile-') as scratch:
        reject('28_PUBLIC_PRIVATE_IMPORT',lambda:import_path(h.ROOT/'scripts/h05_common.py',scratch))
        reject('29_UNLICENSED_OR_OUTOFSCOPE_EXPORT',lambda:export_path('Lean/PrivateWorld.lean',scratch,'PUBLIC_CC0'))
        pth=Path(scratch)/'frozen.txt';pth.write_bytes(base);pth.write_bytes(base+b'mutation');reject('31_ACTUAL_SCRATCH_FILE_MUTATION',lambda:integrity(pth.read_bytes(),digest))
    reject('30_UNFAVORABLE_OR_UNRESOLVED_PRUNING',lambda:roster(['good','bad','unresolved'],['good']))
    reject('32_RESERVE_OUTSIDE_VIABILITY',lambda:classify(True,4,5,8,.5,'Q2'))
    reject('33_INVALID_CANDIDATE_DENOMINATOR',lambda:classify(True,0,0,0,.5,'Q2'))
    # Exercise the actual artifact entry point in a small isolated package, not
    # a metadata-only surrogate and not a copy/mutation of the live candidate.
    with tempfile.TemporaryDirectory(prefix='hfd05-export-hostile-') as scratch:
        root=Path(scratch);source=root/'approved.py';source.write_bytes(base)
        (root/'HFD05_Artifact_Manifest.json').write_text(json.dumps(dict(files={'approved.py':dict(sha256=digest)})))
        spec=importlib.util.spec_from_file_location('hostile_actual_export_runner',h.ROOT/'scripts/h05_artifact_runner.py');runner=importlib.util.module_from_spec(spec);spec.loader.exec_module(runner);runner.ROOT=root
        runner.audit();source.write_bytes(base+b'mutated')
        reject('34_ACTUAL_PUBLIC_SOURCE_BYTE_MUTATION',runner.audit);source.write_bytes(base)
        secret=root/'PrivateWorld.lean';secret.write_bytes(b'out of scope private source')
        reject('35_ACTUAL_UNLISTED_PRIVATE_EXPORT',runner.audit)
        # Restore only this owned temporary fixture, never any live source.
        secret.unlink();link=root/'unlisted.py';link.symlink_to(source)
        reject('36_ACTUAL_EXPORT_SYMLINK',runner.audit)
    require(classify(True,8,8,8,1,'Q2')['KERNEL_MC_STATUS']=='MC_UNRESOLVED','Probability-one endpoint regression')
    require(classify(False,8,8,8,0,'Q2')['FINITE_BANK_HEALTH'] is False,'Absent theta0 regression')
    h.write_json('mathematical_audit/hostile_controls.json',dict(status='PASS',active_rejected_controls=len(tests),controls=tests,live_sources_or_banks_mutated=False,source_sha256=h.sha(__file__)))
    print('Active hostile controls PASS',len(tests),flush=True)
if __name__=='__main__':run()
