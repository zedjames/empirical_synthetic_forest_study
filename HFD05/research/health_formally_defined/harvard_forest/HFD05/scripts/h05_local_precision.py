"""Adverse primary capacity extremes at fixed higher Monte Carlo prefixes."""
import json,math
import numpy as np
import h05_common as h
from h05_growth import simulate
from h05_exposure import records
from h05_exposure_coupled import reconstruct
from h05_bank_counts import count
from h05_inference import classify
from h05_pipeline_primary import paths

def select():
    old,run=paths();surface=h.read_csv(h.ROOT/'corrected_pipeline/primary_surface.csv');selected=[]
    design=h.original.frozen.load('config/design.json');scenario_ids=[s['id'] for s in design['scenarios']]
    _,_,_,core,model=h.original.load_data();original=np.load(old/'primary_histories.npy',mmap_mode='r');ol=np.load(old/'primary_lawful.npy',mmap_mode='r')
    for result in ['R3','R4']:
        rows=[r for r in surface if r['result']==result]
        field=lambda r:'local_max_abs_viable_delta' if r['query']=='Q1' else 'local_max_abs_reserve_delta'
        maximum=max(float(r[field(r)]) for r in rows);attaining=[r for r in rows if float(r[field(r)])==maximum]
        attaining.sort(key=lambda r:(r['scenario'],int(r['horizon']),float(r['alpha']),float(r['beta']),int(r['tau']),r['query'],float(r['rho']),float(r['theta'])))
        bank=np.load(run/(result+'_histories.npy'),mmap_mode='r');law=np.load(run/(result+'_lawful.npy'),mmap_mode='r');seen=set()
        for r in attaining:
            si=scenario_ids.index(r['scenario']);cell=(int(r['horizon']),float(r['alpha']),float(r['beta']),int(r['tau']),r['query'],float(r['rho']),float(r['theta']))
            ob=original[:,:,:,si].reshape(-1,64,21,6);cb=bank[:,:,:,si].reshape(-1,64,21,6);olf=ol[:,:,:,si].reshape(-1,64);clf=law[:,:,:,si].reshape(-1,64)
            for index in range(len(ob)):
                nv0,nr0,p0=count(ob[index],olf[index],model,[cell]);nv,nr,p=count(cb[index],clf[index],model,[cell]);delta=int((nv-nv0)[0] if r['query']=='Q1' else (nr-nr0)[0])/64
                if abs(delta)!=maximum:continue
                state,rep,vi=np.unravel_index(index,(256,2,2));physical=(int(state),int(rep),int(vi),si)
                if physical in seen:continue
                seen.add(physical);selected.append(dict(selection_index=len(selected),selected_result=result,state=physical[0],parameter_rep=physical[1],variant_index=physical[2],scenario_index=si,scenario=r['scenario'],horizon=cell[0],alpha=cell[1],beta=cell[2],tau=cell[3],query=cell[4],rho=cell[5],theta=cell[6],original_K=64,original_target_delta=delta,global_abs_maximum=maximum))
                if len(seen)==8:break
            if len(seen)==8:break
        if len(seen)!=8:raise ValueError('Fewer than8 unique extremes; requires explicit selection amendment')
    h.write_csv('corrected_pipeline/local_extreme_selection.csv',selected);return selected,core,model,design
def run():
    h.guard('local_precision');selected,core,model,design=select();_,e0,e1,_,_=h.original.load_data();_,groups=records(e0,e1,core);seed=h.original.frozen.Seeds();rows=[]
    for record in selected:
        i,rep,vi,si=[record[k] for k in ['state','parameter_rep','variant_index','scenario_index']];a,b,ledger=reconstruct(core,model,groups,seed.rng('reconstruct',i),h.rng('individual_state',0,i),'B0');params=h.original.demography.parameters(model,seed.rng('parameters',i,rep));scenario=design['scenarios'][si];variant=['gamma_growth','normal_growth'][vi]
        cell=tuple(record[k] for k in ['horizon','alpha','beta','tau','query','rho','theta'])
        for identity,state,corrected in [('ORIGINAL',a,False),('INDIVIDUAL_ORIGINAL_GAMMA',b,False),('INDIVIDUAL_CORRECTED_GAMMA',b,True)]:
            nv=nr=0
            for chunk in range(16):
                bank,lawful=simulate(state,params,scenario,record['horizon'],256,h.rng('corrected_future',100,record['selection_index'],chunk),variant,corrected)
                vv,rr,p=count(bank,lawful,model,[cell]);nv+=int(vv[0]);nr+=int(rr[0]);K=(chunk+1)*256
                if K not in [256,1024,4096]:continue
                d=classify(bool(p[0]),nv,nr,K,record['theta'],record['query']);n=nv if record['query']=='Q1' else nr
                lo,hi=h.original.prior.wilson(np.array([n]),K,.975)
                rows.append(dict(**record,comparison=identity,K=K,nv=nv,nr=nr,present=bool(p[0]),mass=n/K,margin=n/K-record['theta'],finite_health=d['FINITE_BANK_HEALTH'],MC_status=d['KERNEL_MC_STATUS'],distinct_law_difference_lower=float(lo[0]),distinct_law_difference_upper=float(hi[0]),N=int(state['n'].sum()),BA=float(np.sum(state['n']*state['d']**2)*math.pi/40000),J=int(np.sum(state['n']*(state['d']<10))),targeted_not_representative=True))
        print('Local extreme higher precision',record['selection_index'],'of16',flush=True)
        h.write_csv('corrected_pipeline/local_precision_partial.csv',rows)
    h.write_csv('corrected_pipeline/local_precision.csv',rows)
    comparisons=[]
    for record in selected:
        for K in [256,1024,4096]:
            values={r['comparison']:r for r in rows if r['selection_index']==record['selection_index'] and r['K']==K};old=values['ORIGINAL']
            for name in ['INDIVIDUAL_ORIGINAL_GAMMA','INDIVIDUAL_CORRECTED_GAMMA']:
                new=values[name];comparisons.append(dict(selection_index=record['selection_index'],K=K,comparison=name,mass_delta=new['mass']-old['mass'],difference_lower=new['distinct_law_difference_lower']-old['distinct_law_difference_upper'],difference_upper=new['distinct_law_difference_upper']-old['distinct_law_difference_lower'],finite_health_changed=new['finite_health']!=old['finite_health'],original_status=old['MC_status'],new_status=new['MC_status'],targeted_not_representative=True,interval_scope='Conservative two-distinct-law marginal MC bracket, not general ecological robustness'))
    h.write_csv('corrected_pipeline/local_precision_contrasts.csv',comparisons);h.write_json('corrected_pipeline/local_precision_completion.json',dict(status='COMPLETE',selected_units=16,prefixes=[256,1024,4096],no_uniform_robustness_claim=True))
if __name__=='__main__':run()
