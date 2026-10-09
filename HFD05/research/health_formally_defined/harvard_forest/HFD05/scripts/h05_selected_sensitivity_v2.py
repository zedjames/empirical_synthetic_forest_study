"""All original selected cap and representation future grids under explicit corrections."""
import itertools,json,types
import numpy as np
import h05_common as h
from h05_growth import simulate,laws,sample
from h05_exposure import records
from h05_exposure_coupled import reconstruct
from h05_inference import classify
import h04_juvenile as juvenile
import h03_calibration as cal

class FineGrowth:
    def __init__(self,generator,params,scenario):self.g=generator;self.params=params;self.scenario=scenario
    def __getattr__(self,k):return getattr(self.g,k)
    def gamma(self,shape,scale,size=None):
        cells=np.concatenate([np.repeat(np.arange(64),2),np.tile(np.arange(16)*4,(size[1]-128)//16)]);mu=self.params['growth'][cells]*self.scenario['growth_factor'];sd=self.params['sd'][cells]*self.scenario['growth_factor'];old=laws(mu,sd,False)
        if not np.array_equal(shape,old['shape']) or not np.array_equal(scale,old['scale']):raise ValueError('Fine growth actual call mismatch')
        return sample(mu,sd,size,self.g,True)
def cap():
    h.guard('selected_sensitivity_v2');sources,e0,e1,core,base=h.original.load_data();_,groups=records(e0,e1,core);c=h.original.config()['effective_count'];design=h.original.frozen.load('config/design.json');outputs=[]
    original=h.read_csv(h.FOUR/'effective_count/propagation.csv')
    for r in original:outputs.append(dict(result='R0',**r,exposure='POOLED',historical_status=True))
    for ci,cap in enumerate(c['caps']):
        empirical,demographic=h.original.isolated_models(cap);model=empirical.fit(core)
        # Clone the unchanged composing function with an explicitly local demographic
        # dependency. Original globals/modules remain untouched. Reconstruction DOES
        # sample recruitment parameters at this cap, so cap200 must not be substituted.
        local_globals=dict(reconstruct.__globals__)
        local_globals['h']=types.SimpleNamespace(original=types.SimpleNamespace(prior=h.original.prior,demography=demographic))
        local_reconstruct=types.FunctionType(reconstruct.__code__,local_globals,'cap_specific_controlled_reconstruct',reconstruct.__defaults__)
        for i in range(16):
            old,new,ledger=local_reconstruct(core,model,groups,h.original.rng('cap_state',i,0),h.rng('individual_state',10+ci,i),'B0')
            params=demographic.parameters(model,h.original.rng('cap_state',i,1))
            for si,scenario in enumerate(design['scenarios']):
                for vi,variant in enumerate(c['variants']):
                    for identity,state,corrected in [('R1',old,False),('R2',old,True),('R3',new,False),('R4',new,True)]:
                        bank,lawful=simulate(state,params,scenario,10,256,h.original.rng('cap_future',i,si,vi),variant,corrected)
                        for horizon,rho,query,theta in itertools.product(c['horizons'],[.5,1],['Q1','Q2'],[.5,.75,.9]):
                            rr=demographic.response(bank,lawful,model,horizon,.75,.4,2,rho);nv=int(rr['viable'].sum());nr=int(rr['reserve_viable'].sum());d=classify(rr['present'],nv,nr,256,theta,query)
                            outputs.append(dict(result=identity,cap=str(cap),state=str(i),scenario=scenario['id'],variant=variant,horizon=str(horizon),rho=str(rho),query=query,theta=str(theta),present=str(rr['present']),nv=str(nv),nr=str(nr),K='256',FINITE_BANK_HEALTH=str(d['FINITE_BANK_HEALTH']),KERNEL_MC_STATUS=d['KERNEL_MC_STATUS'],exposure='CONTROLLED_INDIVIDUAL' if identity in ['R3','R4'] else 'POOLED',historical_status=False))
            if i%4==0:print('Corrected cap',cap,i,flush=True)
        h.write_csv('corrected_pipeline/cap_surface_partial.csv',outputs)
    oldmap={tuple(r[k] for k in ['cap','state','scenario','variant','horizon','rho','query','theta']):r for r in original}
    for r in outputs:
        if r['result']!='R1':continue
        key=tuple(r[k] for k in ['cap','state','scenario','variant','horizon','rho','query','theta']);ref=oldmap[key]
        if any(r[k]!=ref[k] for k in ['nv','nr','present','FINITE_BANK_HEALTH']):raise ValueError('Original cap replay differs')
    h.write_csv('corrected_pipeline/cap_surface.csv',outputs);h.write_json('corrected_pipeline/cap_completion.json',dict(status='COMPLETE',caps=4,states=16,variants=5,rows=len(outputs),old_finite_reproduced=True))
def representation():
    h.guard('selected_sensitivity_v2');sources,e0,e1,core,base=h.original.load_data();held,path=juvenile.visible_inputs('M0');visible,model=cal.fit_visible(sources,path,held);profile,_=juvenile.profiles(e0,path);c=h.original.config()['juvenile'];design=h.original.frozen.load('config/design.json');scenarios=[s for s in design['scenarios'] if s['id'] in c['scenarios']];outputs=[]
    originals=h.read_csv(h.FOUR/'juvenile_boundary/propagation.csv')
    for r in originals:outputs.append(dict(result='R0',**r,exposure='NOT_REEVALUATED_FINE_PROFILE_OBSTRUCTION'))
    for i in range(16):
        coarse=h.original.demography.reconstruct(visible,model,h.original.rng('representation_state',i,0));fine,clip=juvenile.enrich(coarse,profile,h.original.rng('representation_state',i,1));params=h.original.demography.parameters(model,h.original.rng('representation_state',i,2))
        for si,scenario in enumerate(scenarios):
            for vi,variant in enumerate(c['variants']):
                for identity,corrected in [('R1',False),('R2',True)]:
                    g=h.original.rng('representation_future',i,si,vi)
                    cb,cl=simulate(coarse,params,scenario,10,256,g,variant,corrected)
                    g=h.original.rng('representation_future',i,si,vi)
                    if corrected and variant=='gamma_growth':g=FineGrowth(g,params,scenario)
                    fb,fl=juvenile.simulate_enriched(fine,params,scenario,10,256,g,variant)
                    for horizon,beta,rho,theta in itertools.product(c['horizons'],c['semantic_grid']['beta'],c['semantic_grid']['rho'],c['semantic_grid']['theta']):
                        a,b=[h.original.demography.response(bank,law,model,horizon,.75,beta,2,rho) for bank,law in [(cb,cl),(fb,fl)]];nv=[int(x['viable'].sum()) for x in [a,b]];nr=[int(x['reserve_viable'].sum()) for x in [a,b]];d=[classify(x['present'],v,r,256,theta,'Q2') for x,v,r in zip([a,b],nv,nr)]
                        flags=[a['present']!=b['present'],(nv[0]/256>=theta)!=(nv[1]/256>=theta),(nr[0]/256>=theta)!=(nr[1]/256>=theta)];changed=d[0]['FINITE_BANK_HEALTH']!=d[1]['FINITE_BANK_HEALTH'];reason='NONE' if not changed else (['INITIAL_REALIZATION','CONTINUATION','RESERVE'][flags.index(True)] if sum(flags)==1 else 'COMBINED')
                        outputs.append(dict(result=identity,state=str(i),scenario=scenario['id'],variant=variant,horizon=str(horizon),beta=str(beta),rho=str(rho),theta=str(theta),K='256',coarse_present=str(a['present']),fine_present=str(b['present']),coarse_nv=str(nv[0]),fine_nv=str(nv[1]),coarse_nr=str(nr[0]),fine_nr=str(nr[1]),coarse_Health=str(d[0]['FINITE_BANK_HEALTH']),fine_Health=str(d[1]['FINITE_BANK_HEALTH']),coarse_status=d[0]['KERNEL_MC_STATUS'],fine_status=d[1]['KERNEL_MC_STATUS'],change_class=reason,exposure='NOT_REEVALUATED_FINE_PROFILE_OBSTRUCTION'))
        if i%4==0:print('Corrected representation',i,flush=True)
    oldmap={tuple(r[k] for k in ['state','scenario','variant','horizon','beta','rho','theta']):r for r in originals}
    for r in outputs:
        if r['result']!='R1':continue
        old=oldmap[tuple(r[k] for k in ['state','scenario','variant','horizon','beta','rho','theta'])]
        if any(r[k]!=old[k] for k in ['coarse_nv','fine_nv','coarse_nr','fine_nr','coarse_Health','fine_Health']):raise ValueError('Original representation differs')
    h.write_csv('corrected_pipeline/representation_surface.csv',outputs);h.write_json('corrected_pipeline/representation_completion.json',dict(status='COMPLETE',old_finite_reproduced=True,rows=len(outputs),exposure='NOT_REEVALUATED_FINE_PROFILE_OBSTRUCTION',dynamic_commutation=False))
if __name__=='__main__':
    import sys
    {'cap':cap,'representation':representation}[sys.argv[1]]()
