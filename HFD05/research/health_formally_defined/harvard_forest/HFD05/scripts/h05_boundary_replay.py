"""Complete corrected synthetic references and context predictions, no empirical-date invention."""
import itertools,json,hashlib
from collections import defaultdict
import numpy as np
import h05_common as h
from h05_growth import simulate
from h05_inference import classify
import h04_boundary as boundary
import h03_synthetic as prior

def setup():
    _,_,_,core,model=h.original.load_data()
    anchor=dict(n=core['stats']['n0'].copy(),d=model['diameter0'].copy())
    worlds=json.loads((h.FOUR/'config/world_design.json').read_text())['worlds']
    return core,model,anchor,worlds
def historical(core,model):
    cfg=h.original.prior.PROTOCOL['synthetic']
    grid=list(itertools.product(range(3),range(3),range(3),range(3),range(cfg['replicates_per_cell'])))
    for world,(s,r,p,recovery,rep) in enumerate(grid):
        state=prior.make_state(core,model,s,r,world)
        params=prior.est.public_parameters(model,prior.est.PublicContext(p,recovery),h.original.prior.rng('world',world,1))
        g=h.original.prior.rng('observation',world);n=g.binomial(state['n'],cfg['observation_thinning']);d=np.clip(state['d']+g.normal(0,cfg['dbh_noise_sd_cm'],64),1,300)
        packet=h.original.ObservationPacket(tuple(map(int,n)),tuple(float(v) if nn else None for nn,v in zip(n,d)),(cfg['observation_thinning'],)*64)
        yield world,state,params,packet,p,recovery
def row(cohort,world,suite,K,present,nv,nr,identity,query,**extra):
    d=classify(bool(present),int(nv),int(nr),int(K),.75,query)
    return dict(cohort=cohort,world=world,suite=suite,K=K,query=query,result=identity,present=bool(present),nv=int(nv),nr=int(nr),finite_health=d['FINITE_BANK_HEALTH'],MC_status=d['KERNEL_MC_STATUS'],mass=(nv if query=='Q1' else nr)/K,margin=(nv if query=='Q1' else nr)/K-.75,**extra)
def references():
    h.guard('boundary_replay');core,model,anchor,worlds=setup();outputs=[];hashes=[]
    original=h.read_csv(h.FOUR/'boundary_validation/truth_prefixes.csv')
    for r in original:
        if int(r['K']) not in [4096,16384,65536]:continue
        outputs.append(row('boundary',int(r['world']),r['suite'],int(r['K']),r['present']=='True',int(r['nv']),int(r['nr']),'ORIGINAL_REF',r['query'],family=r['family'],bank_source='Protected original reference'))
    for world,cell in enumerate(worlds):
        state,params,_=boundary.world(core,model,cell,world)
        for suite in h.original.config()['boundary']['misspecifications']:
            nv=nr=0;digest=hashlib.sha256()
            for chunk in range(16):
                g=h.original.rng('truth',world,chunk)
                if suite=='entry_lognormal_annual_sd_1.2':g=h.original.EntryRNG(g)
                bank,lawful=simulate(state,params,h.original.scenario(suite),5,1024,g,'gamma_growth',True)
                rr=h.original.response(bank,lawful,model);nv+=int(rr['viable'].sum());nr+=int(rr['reserve_viable'].sum());digest.update(bank.tobytes());digest.update(lawful.tobytes())
                K=(chunk+1)*1024
                if K in [4096,16384]:
                    for query in ['Q1','Q2']:outputs.append(row('boundary',world,suite,K,rr['present'],nv,nr,'CORRECTED_REF',query,family=cell['family'],bank_source=digest.hexdigest()))
            hashes.append(dict(cohort='boundary',world=world,suite=suite,K=16384,bank_sha256=digest.hexdigest()))
        if world%3==0:print('Corrected boundary references',world,'of81',flush=True)
        h.write_csv('corrected_pipeline/reference_counts_partial.csv',outputs)
    old={(int(r['world']),r['suite']):r for r in h.read_csv(h.THREE/'synthetic_validation/truth_census.csv')}
    suites=h.original.prior.PROTOCOL['synthetic']['misspecifications']
    for world,state,params,packet,p,recovery in historical(core,model):
        for si,suite in enumerate(suites):
            reference=old[(world,suite)]
            for query in ['Q1','Q2']:
                outputs.append(row('HFD03',world,suite,4096,reference['truth_present']=='True',round(float(reference['truth_mass'])*4096),round(float(reference['truth_reserve'])*4096),'ORIGINAL_REF',query,family='HISTORICAL_FACTORIAL',bank_source='Protected original reference'))
            g=h.original.prior.rng('truth_future',world,si)
            if si==2:g=h.original.EntryRNG(g)
            bank,lawful=simulate(state,params,h.original.scenario(suite),5,4096,g,'gamma_growth',True)
            rr=h.original.response(bank,lawful,model);nv=int(rr['viable'].sum());nr=int(rr['reserve_viable'].sum());digest=hashlib.sha256(bank.tobytes()+lawful.tobytes()).hexdigest()
            for query in ['Q1','Q2']:outputs.append(row('HFD03',world,suite,4096,rr['present'],nv,nr,'CORRECTED_REF',query,family='HISTORICAL_FACTORIAL',bank_source=digest))
            hashes.append(dict(cohort='HFD03',world=world,suite=suite,K=4096,bank_sha256=digest))
        if world%12==0:print('Corrected historical references',world,'of324',flush=True)
        h.write_csv('corrected_pipeline/reference_counts_partial.csv',outputs)
    h.write_csv('corrected_pipeline/reference_counts.csv',outputs);h.write_csv('corrected_pipeline/reference_bank_hashes.csv',hashes)
    h.write_json('corrected_pipeline/reference_completion.json',dict(status='COMPLETE',boundary_worlds=81,historical_worlds=324,suites=3,corrected_boundary_K=[4096,16384],corrected_historical_K=4096,corrected65536='NOT_REEVALUATED',source_sha256=h.sha(__file__)))
def predictions():
    h.guard('boundary_replay');core,model,anchor,worlds=setup();outputs=[]
    original=h.read_csv(h.FOUR/'context_ablation/estimator_draws.csv')
    by={(r['cohort'],int(r['world']),r['policy'],int(r['draw'])):r for r in original}
    packets=[]
    for world,cell in enumerate(worlds):
        state,params,_=boundary.world(core,model,cell,world);packets.append(('boundary',world,world,boundary.make_packet(state,world),cell['pressure'],cell['recovery']))
    for world,state,params,packet,p,recovery in historical(core,model):packets.append(('HFD03',world,10000+world,packet,p,recovery))
    for cohort,world,identity,packet,p,r in packets:
        for policy in ['B0_full','B1_coarsened','B2_unknown','B3_label_only']:
            observation=packet if policy!='B3_label_only' else h.original.ObservationPacket((0,)*64,(None,)*64,(0.,)*64)
            marker=boundary.context(policy,p,r);levels=marker.levels()
            for draw in range(32):
                ref=by[(cohort,world,policy,draw)]
                state=h.original.estimate_state(observation,anchor,h.original.rng('estimate_state',identity,draw))
                g=h.original.rng('estimate_params',identity,draw);pp,rrr=levels[int(g.integers(len(levels)))];params=h.original.parameters(model,pp,rrr,g)
                bank,lawful=simulate(state,params,h.original.scenario(),5,256,h.original.rng('estimate_future',identity,draw),'gamma_growth',True)
                rr=h.original.response(bank,lawful,model);nv=int(rr['viable'].sum());nr=int(rr['reserve_viable'].sum())
                for query in ['Q1','Q2']:
                    for identity2,present,nvv,nrr in [('ORIGINAL_PRED',ref['present']=='True',int(ref['nv']),int(ref['nr'])),('CORRECTED_PRED',rr['present'],nv,nr)]:
                        outputs.append(row(cohort,world,'nominal_prediction',256,present,nvv,nrr,identity2,query,policy=policy,draw=draw,exposure='INAPPLICABLE_NO_INDIVIDUAL_DATES',oracle_input=False))
        if world%12==0:print('Corrected context predictions',cohort,world,flush=True)
        h.write_csv('corrected_pipeline/context_draws_partial.csv',outputs)
    h.write_csv('corrected_pipeline/context_draws.csv',outputs)
    h.write_json('corrected_pipeline/context_completion.json',dict(status='COMPLETE',worlds=405,policies=4,draws=32,paths=256,rows=len(outputs),oracle_inputs=False,exposure_factor='INAPPLICABLE'))
def analyze():
    h.guard('boundary_replay');from h04_analysis import score
    if any(json.loads((h.ROOT/'corrected_pipeline'/name).read_text())['status']!='COMPLETE' for name in ['reference_completion.json','context_completion.json']):raise ValueError('Incomplete crossmatrix')
    refs=h.read_csv(h.ROOT/'corrected_pipeline/reference_counts.csv');raw=h.read_csv(h.ROOT/'corrected_pipeline/context_draws.csv');groups=defaultdict(list)
    for r in raw:groups[(r['cohort'],r['world'],r['policy'],r['query'],r['result'])].append(r)
    outputs=[]
    for ref in refs:
        if ref['K']=='65536':continue
        for policy,pred in itertools.product(['B0_full','B1_coarsened','B2_unknown','B3_label_only'],['ORIGINAL_PRED','CORRECTED_PRED']):
            records=groups[(ref['cohort'],ref['world'],policy,ref['query'],pred)];health=sum(r['finite_health']=='True' for r in records)/32;mass=np.mean([float(r['mass']) for r in records]);pl=sum(r['MC_status']=='TRUE' for r in records)/32;pu=sum(r['MC_status']!='FALSE' for r in records)/32
            v=np.array([int(r['nv'])/256 for r in records]);reserve=np.array([int(r['nr'])/256 for r in records]);vl,vu=np.quantile(v,[.05,.95]);rl,ru=np.quantile(reserve,[.05,.95])
            outputs.append(dict(cohort=ref['cohort'],world=ref['world'],family=ref['family'],suite=ref['suite'],reference_K=int(ref['K']),query=ref['query'],policy=policy,crossmatrix=pred+'_'+ref['result'],reference_present=ref['present']=='True',reference_health=ref['finite_health']=='True',reference_status=ref['MC_status'],reference_mass=float(ref['mass']),margin=float(ref['margin']),estimated_mass=float(mass),mass_absolute_error=abs(float(mass)-float(ref['mass'])),reference_viable=int(ref['nv'])/int(ref['K']),reference_reserve=int(ref['nr'])/int(ref['K']),estimated_viable=float(v.mean()),estimated_reserve=float(reserve.mean()),viable_lower=float(vl),viable_upper=float(vu),reserve_lower=float(rl),reserve_upper=float(ru),health_fraction=health,predicted_health=health>=.5,estimate_status='TRUE' if pl>.5 else ('FALSE' if pu<.5 else 'MC_UNRESOLVED'),health_MC_lower=pl,health_MC_upper=pu,uncertainty_scope='Constructed-world paired current32-draw design'))
    h.write_csv('corrected_pipeline/reference_prediction_crossmatrix.csv',outputs)
    grouped=defaultdict(list)
    for r in outputs:
        key=(r['cohort'],r['suite'],r['reference_K'],r['query'],r['policy'],r['crossmatrix'])
        grouped[key+('ALL',)].append(r)
        if r['reference_present'] and abs(r['margin'])<=.1:grouped[key+('PRESENT_NEAR_0.1',)].append(r)
    summaries=[]
    for key,records in sorted(grouped.items()):summaries.append(dict(cohort=key[0],suite=key[1],reference_K=key[2],query=key[3],policy=key[4],crossmatrix=key[5],panel=key[6],**score(records)))
    h.write_csv('corrected_pipeline/reference_crossmatrix_scores.csv',summaries)
    h.write_json('corrected_pipeline/reference_crossmatrix_summary.json',dict(status='COMPLETE',paired_cells=len(outputs),all_worlds_retained=True,reference65536_corrected='NOT_REEVALUATED',exposure_factor='INAPPLICABLE_NO_DATES',old_truths_overwritten=False))
if __name__=='__main__':
    import sys
    {'references':references,'predictions':predictions,'analyze':analyze}[sys.argv[1]]()
