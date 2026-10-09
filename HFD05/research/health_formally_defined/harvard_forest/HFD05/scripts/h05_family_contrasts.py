"""All five latent families paired on corrected conditional state/model units."""
import json
import numpy as np
import h05_common as h
def run():
    import h03_missingness as old
    cells=old.query_grid(h.original.frozen.load('config/design.json'));manifest=json.loads((h.ROOT/'corrected_pipeline/missingness_bank_manifest.json').read_text());lookup={(r['family'],r['result']):r for r in manifest['banks']};outputs=[]
    for result in ['R2','R3','R4']:
        base=np.load(h.ROOT/lookup[('B0',result)]['private_path'])
        for family in ['B1','B2','B3','B4']:
            bank=np.load(h.ROOT/lookup[(family,result)]['private_path'])
            for si,scenario in enumerate(h.original.frozen.load('config/design.json')['scenarios']):
                for j,(horizon,alpha,beta,tau,query,rho,theta) in enumerate(cells):
                    if theta!=.75:continue
                    n=bank['nv'][si,:,j] if query=='Q1' else bank['nr'][si,:,j];b=base['nv'][si,:,j] if query=='Q1' else base['nr'][si,:,j];margin=b.astype(float)/256-theta;near=base['present'][si,:,j]&(abs(margin)<=.1);d=n.astype(float)/256-b.astype(float)/256
                    flips=(bank['finite'][si,:,j]!=base['finite'][si,:,j]);outputs.append(dict(family=family,result=result,reference_family='B0',scenario=scenario['id'],horizon=horizon,alpha=alpha,beta=beta,tau=tau,query=query,rho=rho,theta=theta,units=256,reference_present_near_units=int(near.sum()),near_health_flip_units=int((flips&near).sum()),all_health_flip_units=int(flips.sum()),realization_flip_units=int((bank['present'][si,:,j]!=base['present'][si,:,j]).sum()),mean_capacity_delta=float(d.mean()),maximum_local_abs_delta=float(abs(d).max()),q95_local_abs_delta=float(np.quantile(abs(d),.95)),near_max_abs_delta=float(abs(d[near]).max()) if near.any() else None,conditional_pairing='Original shared state/future identities; missingness and survivor RNG consumption can differ',weights='Finite design, not posterior'))
    h.write_csv('entry_missingness/family_boundary_contrasts.csv',outputs)
    summary=[]
    for result in ['R2','R3','R4']:
        for family in ['B1','B2','B3','B4']:
            chosen=[r for r in outputs if r['result']==result and r['family']==family];summary.append(dict(result=result,family=family,semantic_cells=len(chosen),cells_with_present_near_reversals=sum(r['near_health_flip_units']>0 for r in chosen),maximum_near_reversal_units=max(r['near_health_flip_units'] for r in chosen),maximum_all_reversal_units=max(r['all_health_flip_units'] for r in chosen),maximum_local_capacity_delta=max(r['maximum_local_abs_delta'] for r in chosen),scope='Conditional finite design; no uniform state-space or causal ecology claim'))
    h.write_csv('entry_missingness/family_boundary_summary.csv',summary)
if __name__=='__main__':run()
