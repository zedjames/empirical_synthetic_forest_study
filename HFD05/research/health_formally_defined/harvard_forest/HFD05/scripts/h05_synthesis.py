"""Derived full-census sensitivity, exact claims and figure-source tables."""
import itertools,json,math
from collections import defaultdict
import numpy as np
import h05_common as h

def entry():
    original=h.read_csv(h.THREE/'recruitment_audit/candidate_audit.csv')
    if len(original)!=6992:raise ValueError('Entry candidate census changed')
    requirements=[('THRESHOLD_CROSSING','Secure below1cm predecessor matched to above1cm dated successor','Missing E0 stem record cannot distinguish below-threshold from missed eligible stem'),('PREVIOUSLY_UNRECORDED_ELIGIBLE_STEM','Earlier exhaustive frame and linked earlier size or verified omission','A large E1 diameter does not supply missing E0 size'),('NEW_STEM_EXISTING_PLANT','Stable plant association plus dated stem emergence','E1 plant table absent; new stem ID not emergence date'),('IDENTITY_TAG_CHANGE','Resolved old/new tag and plant relation','Candidate collisions and association conflicts retained'),('SAMPLING_FRAME_CHANGE','Comparable frame footprints and documented coverage','No invented locations or completion of unobserved area'),('UNRESOLVED_ENTRY','Additional record resolving competing mechanisms','Observed later live candidate compatible with multiple latent histories')]
    h.write_csv('entry_missingness/entry_identification_matrix.csv',[dict(latent_class=a,necessary_observables=b,identification_obstruction=c,candidate_census=6992,diagnostic_weights_are_posterior=False,full_biological_entry_identified=False) for a,b,c in requirements])
    h.write_csv('entry_missingness/candidate_provenance.csv',[dict(**r,source_identity='Protected original6992 candidate audit',biological_class_certified=False,weight_interpretation='Diagnostic membership bound, not posterior') for r in original])
    # A concrete observation-map witness, not fabricated field evidence.
    witness=[dict(assignment='A_THRESHOLD_CROSSING',prior_stem_state='Below detection threshold, unrecorded',later_record='Same fixed live tagged size/date/identity record',earlier_record='No recorded stem',observation_equal=True,biological_entry=True),dict(assignment='B_EARLIER_UNRECORDED_STEM',prior_stem_state='Already eligible but missed in earlier frame',later_record='Same fixed live tagged size/date/identity record',earlier_record='No recorded stem',observation_equal=True,biological_entry=False)]
    h.write_csv('entry_missingness/observational_nonidentification_witness.csv',witness)
def missingness():
    import h03_missingness as old
    cells=old.query_grid(h.original.frozen.load('config/design.json'));manifest=json.loads((h.ROOT/'corrected_pipeline/missingness_bank_manifest.json').read_text())
    if manifest['status']!='COMPLETE':raise ValueError('Incomplete missingness census')
    by={(r['family'],r['result']):r for r in manifest['banks']};rows=[];defaults=[]
    for family in ['B0','B1','B2','B3','B4']:
        original=np.load(h.original.prior.RUN/('missingness_'+family+'.npz'));ov=np.rint(original['mass']*256).astype(int);orr=np.rint(original['reserve']*256).astype(int);op=np.repeat(original['point'][:,:,::5],5,axis=2)
        for result in ['R2','R3','R4']:
            record=by[(family,result)];path=h.ROOT/record['private_path']
            if h.sha(path)!=record['sha256']:raise ValueError('Corrected missingness bank mismatch')
            bank=np.load(path)
            for si,scenario in enumerate(h.original.frozen.load('config/design.json')['scenarios']):
                for j,(horizon,alpha,beta,tau,query,rho,theta) in enumerate(cells):
                    if theta!=.75:continue
                    n=bank['nv'][si,:,j].astype(int) if query=='Q1' else bank['nr'][si,:,j].astype(int);on=ov[si,:,j] if query=='Q1' else orr[si,:,j];present=bank['present'][si,:,j];margin=n/256-theta;original_margin=on/256-theta;near=op[si,:,j]&(abs(original_margin)<=.1);delta=(n-on)/256
                    rows.append(dict(family=family,result=result,scenario=scenario['id'],horizon=horizon,alpha=alpha,beta=beta,tau=tau,query=query,rho=rho,theta=theta,units=len(n),K=256,present_units=int(present.sum()),realization_flip_units=int((present!=op[si,:,j]).sum()),original_present_near_units=int(near.sum()),near_finite_health_flip_units=int(((bank['finite'][si,:,j]!=original['point'][si,:,j])&near).sum()),mean_margin=float(margin.mean()),minimum_margin=float(margin.min()),maximum_margin=float(margin.max()),q05_margin=float(np.quantile(margin,.05)),q95_margin=float(np.quantile(margin,.95)),local_max_abs_capacity_delta=float(abs(delta).max()),local_q95_abs_capacity_delta=float(np.quantile(abs(delta),.95)),near_max_abs_delta=float(abs(delta[near]).max()) if near.any() else None,near_q95_abs_delta=float(np.quantile(abs(delta[near]),.95)) if near.any() else None,design_scope='All original semantic cells, latent families kept separate; no ecological posterior'))
                    if (alpha,beta,tau,rho)!=(.75,.4,2,.5):continue
                    for unit in range(len(n)):
                        defaults.append(dict(family=family,result=result,state=unit//2,variant=['gamma_growth','normal_growth'][unit%2],scenario=scenario['id'],horizon=horizon,query=query,present=bool(present[unit]),original_present=bool(op[si,unit,j]),capacity=n[unit]/256,original_capacity=on[unit]/256,tipping_theta=n[unit]/256,margin=margin[unit],capacity_delta=delta[unit],MC_status=int(bank['status'][si,unit,j]),provenance_link='exposure_model/controlled_origin_comparison.csv',stratum_link='exposure_model/controlled_stratum_comparison.csv'))
    h.write_csv('entry_missingness/boundary_conditioned_missingness.csv',rows);h.write_csv('entry_missingness/default_conditional_tipping.csv',defaults)
    h.write_json('entry_missingness/full_missingness_synthesis.json',dict(status='COMPLETE',families=5,all_corrected_results=3,all_semantic_cells_retained=True,conditional_default_rows=len(defaults),no_uniform_robustness_claim=True,origin_provenance='IMPUTED, anchored to dated measured/imputed components; no observed2020 origin',taxon_and_proxy_strata='exposure_model/controlled_stratum_comparison.csv',posterior_interpretation=False))
def headlines():
    rows=[]
    for source in ['primary','missingness']:
        records=h.read_csv(h.ROOT/'corrected_pipeline'/(source+'_surface.csv'))
        for result in ['R1','R2','R3','R4']:
            selected=[r for r in records if r['result']==result];mx=max(float(r['finite_health_disagreement_vs_original']) for r in selected);local=max(max(float(r['local_max_abs_viable_delta']),float(r['local_max_abs_reserve_delta'])) for r in selected)
            status='EXACTLY_PRESERVED' if mx==0 and local==0 else ('NUMERICALLY_CHANGED_CONCLUSION_STABLE' if mx==0 else 'MODIFIED')
            rows.append(dict(headline=source+'_finite_health_and_capacity',result=result,status=status,cells=len(selected),maximum_finite_health_disagreement=mx,maximum_local_capacity_delta=local,scope='Executed full original grid only, not uniform robustness',evidence='corrected_pipeline/'+source+'_surface.csv'))
    for name,description,status,evidence in [('Q2_MATHEMATICAL_EQUIVALENCE','General nonnegative unrenormalized viable measure; reserve subset','EXACTLY_PRESERVED','mathematical_audit/q2_equivalence.md'),('Q2_KERNEL_INFERENCE','Single95 reserve constraint replaces two97.5 constraints','SUPERSEDED','q2_correction/resolution_comparison.csv'),('GAMMA_GENERATOR_LAW','Reached shape floor changes actual moments','MODIFIED','gamma_growth/branch_summary.csv'),('UNKNOWN_ENTRY','6992 candidates do not identify biological mechanisms','EXACTLY_PRESERVED','entry_missingness/entry_identification_matrix.csv'),('CORRECTED_MAX_REFERENCE65536','Original max remains, corrected max not executed','NOT_REEVALUATED','corrected_pipeline/reference_completion.json'),('FINE_INDIVIDUAL_EXPOSURE','No empirical individual origin-size profile adapter','NOT_REEVALUATED','corrected_pipeline/representation_completion.json'),('EXTERNAL_FULL_HEALTH','Component validation does not identify full Health','NOT_REEVALUATED','external_sources/alignment_matrix.csv')]:rows.append(dict(headline=name,result=description,status=status,cells='',maximum_finite_health_disagreement='',maximum_local_capacity_delta='',scope='Source-qualified explicit claim',evidence=evidence))
    h.write_csv('corrected_pipeline/headline_decisions.csv',rows)
if __name__=='__main__':entry();missingness();headlines()
