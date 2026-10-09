"""Independent full sufficient-count and table audit, without rerunning old writers."""
import json,math
from collections import defaultdict
from statistics import NormalDist
import numpy as np
import h05_common as h

def independent(count,K,theta,present):
    z=NormalDist().inv_cdf(.975);den=K+z*z
    center=(count+z*z/2)/den;radius=z*np.sqrt(count*(1-count/K)+z*z/4)/den
    lower=np.where(count==0,0,np.maximum(0,center-radius));upper=np.where(count==K,1,np.minimum(1,center+radius))
    finite=present&(count/K>=theta);status=np.where(~present,0,np.where(theta==0,1,np.where(upper<theta,0,np.where(lower>theta,1,2))))
    return finite,status
def audit():
    h.predecessor_guard();import h03_missingness as prior
    cells=prior.query_grid(h.original.frozen.load('config/design.json'));theta=np.array([c[6] for c in cells]);q2=np.array([c[4]=='Q2' for c in cells]);manifest=json.loads((h.ROOT/'corrected_pipeline/missingness_bank_manifest.json').read_text());rows=h.read_csv(h.ROOT/'corrected_pipeline/missingness_surface.csv');by={(r['family'],r['result'],r['scenario'],int(r['horizon']),float(r['alpha']),float(r['beta']),int(r['tau']),r['query'],float(r['rho']),float(r['theta'])):r for r in rows};checked=0
    for record in manifest['banks']:
        path=h.ROOT/record['private_path']
        if h.sha(path)!=record['sha256']:raise ValueError('Private corrected bank mutation')
        bank=np.load(path);nv,nr,p=bank['nv'].astype(int),bank['nr'].astype(int),bank['present'];n=np.where(q2,nr,nv)
        if np.any(nr>nv) or np.any(nv>256):raise ValueError('Unconditional subset/count failure')
        finite,status=independent(n,256,theta,p)
        if not np.array_equal(finite,bank['finite']) or not np.array_equal(status,bank['status']):raise ValueError('Independent Wilson/count disagreement')
        for si,scenario in enumerate(h.original.frozen.load('config/design.json')['scenarios']):
            for j,cell in enumerate(cells):
                r=by[(record['family'],record['result'],scenario['id'],*cell)]
                expected={'finite_health_fraction':float(finite[si,:,j].mean()),'present_fraction':float(p[si,:,j].mean()),'viable_mass':float(nv[si,:,j].mean()/256),'reserve_mass':float(nr[si,:,j].mean()/256),'kernel_TRUE_units':int((status[si,:,j]==1).sum()),'kernel_FALSE_units':int((status[si,:,j]==0).sum()),'kernel_MC_UNRESOLVED_units':int((status[si,:,j]==2).sum())}
                if any(float(r[k])!=v for k,v in expected.items()):raise ValueError('Missingness published aggregate mismatch')
        checked+=n.size
    refs=h.read_csv(h.ROOT/'corrected_pipeline/reference_counts.csv');draws=h.read_csv(h.ROOT/'corrected_pipeline/context_draws.csv')
    for r in refs+draws:
        nv,nr,K=int(r['nv']),int(r['nr']),int(r['K']);p=r['present']=='True';query=r['query'];n=nr if query=='Q2' else nv
        if not(0<=nr<=nv<=K):raise ValueError('Reference subset failure')
        f,s=independent(np.array([n]),K,.75,np.array([p]));name=['FALSE','TRUE','MC_UNRESOLVED'][int(s[0])]
        if bool(f[0])!=(r['finite_health']=='True') or name!=r['MC_status']:raise ValueError('Reference/context independent inference mismatch')
    original={(r['cohort'],r['world'],r['policy'],r['draw']):r for r in h.read_csv(h.FOUR/'context_ablation/estimator_draws.csv')};changed=0
    paired=defaultdict(dict)
    for r in draws:
        key=(r['cohort'],r['world'],r['policy'],r['draw'],r['query']);paired[key][r['result']]=r
        if r['result']=='ORIGINAL_PRED':
            old=original[key[:4]]
            if any(r[k]!=old[k] for k in ['nv','nr','present']):raise ValueError('Protected context counts replaced')
    for key,pair in paired.items():
        if set(pair)!={'ORIGINAL_PRED','CORRECTED_PRED'}:raise ValueError('Incomplete original/corrected prediction pair')
        changed+=any(pair['ORIGINAL_PRED'][k]!=pair['CORRECTED_PRED'][k] for k in ['nv','nr','present'])
    refpair=defaultdict(dict)
    for r in refs:
        if r['K']=='65536':continue
        refpair[(r['cohort'],r['world'],r['suite'],r['K'],r['query'])][r['result']]=r
    refchanged=0
    for key,pair in refpair.items():
        if set(pair)!={'ORIGINAL_REF','CORRECTED_REF'}:raise ValueError('Incomplete reference pair')
        refchanged+=any(pair['ORIGINAL_REF'][k]!=pair['CORRECTED_REF'][k] for k in ['nv','nr','present'])
    cross=h.read_csv(h.ROOT/'corrected_pipeline/reference_prediction_crossmatrix.csv');groups=defaultdict(set)
    for r in cross:groups[(r['cohort'],r['world'],r['suite'],r['reference_K'],r['query'],r['policy'])].add(r['crossmatrix'])
    expected={'ORIGINAL_PRED_ORIGINAL_REF','CORRECTED_PRED_ORIGINAL_REF','ORIGINAL_PRED_CORRECTED_REF','CORRECTED_PRED_CORRECTED_REF'}
    if any(v!=expected for v in groups.values()):raise ValueError('Incomplete crossmatrix')
    for name in ['cap_surface.csv','representation_surface.csv']:
        data=h.read_csv(h.ROOT/'corrected_pipeline'/name)
        if len(data)!=(61440 if name=='cap_surface.csv' else 4608):raise ValueError('Selected grid pruned')
    h.write_json('mathematical_audit/pipeline_independent_audit.json',dict(status='PASS',missingness_conditional_counts=checked,missingness_rows=len(rows),reference_rows=len(refs),context_rows=len(draws),context_count_pairs_changed=changed,reference_count_pairs_changed=refchanged,crossmatrix_rows=len(cross),all_crossmatrix_pairs_complete=True,all_selected_grids_retained=True,source_sha256=h.sha(__file__)))
    print('Independent pipeline audit PASS',checked,changed,refchanged,flush=True)
if __name__=='__main__':audit()
