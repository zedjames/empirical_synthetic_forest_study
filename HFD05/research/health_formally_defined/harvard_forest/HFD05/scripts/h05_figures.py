"""Small source-qualified vector figures with exact machine-readable plotting rows."""
import json,html
from collections import defaultdict
import h05_common as h

def table(path):return h.read_csv(h.ROOT/path)
def plot(number,title,units,records,sources,selection,caption):
    name='H05-'+str(number);h.write_csv('figures/data/'+name+'.csv',records)
    width=1050;height=150+45*len(records);values=[abs(float(r['value'])) for r in records];scale=max(values or [1]) or 1
    lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">','<rect width="100%" height="100%" fill="white"/>',f'<text x="24" y="32" font-family="sans-serif" font-size="21">{html.escape(title)}</text>',f'<text x="24" y="58" font-family="sans-serif" font-size="13">{html.escape(units)}</text>']
    for i,r in enumerate(records):
        y=90+45*i;value=float(r['value']);color='#a44130' if value<0 else '#35708c'
        lines.extend([f'<text x="24" y="{y+17}" font-family="sans-serif" font-size="13">{html.escape(r["label"])}</text>',f'<rect x="430" y="{y}" width="{abs(value)/scale*420:.3f}" height="22" fill="{color}"/>',f'<text x="865" y="{y+17}" font-family="sans-serif" font-size="13">{value:.6g}</text>'])
    lines.append(f'<text x="24" y="{height-15}" font-family="sans-serif" font-size="11">{html.escape(name+": source-qualified descriptive plot; see caption and data")}</text></svg>')
    path=h.output_path('figures/'+name+'.svg');path.write_text('\n'.join(lines)+'\n')
    h.write_json('figures/'+name+'.json',dict(id=name,title=title,units=units,data='figures/data/'+name+'.csv',data_sha256=h.sha(h.ROOT/'figures/data'/(name+'.csv')),svg_sha256=h.sha(path),sources={p:h.sha(h.ROOT/p) for p in sources},selection=selection,caption=caption,protocol_sha256=h.sha(h.ROOT/'config/protocol.json'),result_identity='HFD05_SOURCE_QUALIFIED_CORRECTION',ecological_population_inference=False))
def run():
    src='q2_correction/resolution_comparison.csv';records=[]
    for r in table(src):
        records.append(dict(label=r['family']+' unresolved decrease',value=int(r['old_unresolved'])-int(r['corrected_unresolved'])))
    plot(1,'Single reserve constraint: numerical resolution','Count of newly resolved count/query evaluations',records,[src],'All archived count families','Finite-bank Health is exactly preserved; this plot compares pointwise approximate confidence procedures, not ecological performance.')
    src='gamma_growth/branch_summary.csv';records=[]
    for r in table(src):records.append(dict(label=r['family']+' material mean laws',value=r['material_mean_cells']))
    plot(2,'Reached Gamma shape-floor moment defect','Number of ancestral-cell/scenario laws',records,[src],'All confirmatory parameter families','Actual-call branch census; rare reached branches have material moment errors even where finite Health is unchanged. Zero and denominator-floor counts remain in the source table.')
    src='exposure_model/controlled_origin_comparison.csv';groups=defaultdict(list)
    for r in table(src):groups[(r['design'],r['family'])].append(float(r['conditional_individual_survivor_mean'])-float(r['conditional_original_survivor_mean']))
    plot(3,'Unequal exposure: conditional Jensen gap','Mean conditional expected extra dated survivors',[dict(label=' / '.join(k),value=sum(v)/len(v)) for k,v in sorted(groups.items())],[src],'Every primary/missingness reconstruction draw','Conditional on shared hazards, not a marginal independence claim. Finite survivor draws and other provenance remain separate.')
    src='corrected_pipeline/headline_decisions.csv';records=[dict(label=r['headline']+' '+r['result'],value=r['maximum_finite_health_disagreement']) for r in table(src) if r['maximum_finite_health_disagreement']!='']
    plot(4,'Complete correction grids: finite Health changes','Maximum conditional ensemble disagreement per semantic cell',records,[src,'corrected_pipeline/primary_surface.csv','corrected_pipeline/missingness_surface.csv'],'All semantic cells; each correction factor separately','These maxima summarize finite-design changes, not uniform ecological robustness. Continuous local capacity changes are retained in the full tables.')
    src='oracle_attribution/near_boundary_scores.csv';records=[]
    for r in table(src):
        if (r['suite'],r['panel'],r['K'],r['subset'])!=('correct','BASELINE32','256','PRESENT_NEAR_0.1'):continue
        records.append(dict(label=r['query']+' '+r['condition']+' target mass MAE',value=r['viable_MAE'] if r['query']=='Q1' else r['reserve_MAE']))
    plot(5,'Four oracle interfaces: nonadditive error attribution','Constructed present near-boundary target-mass MAE',records,[src,'oracle_attribution/interactions.csv'],'Nominal suite, original max references; both queries and all four roles','Supplying parameters alone can worsen loss. Both-oracle residual includes future MC and, in other suites, predictive misspecification; no additive causal partition.')
    src='external_sources/HF453/mortality_scores.csv';records=[]
    for r in table(src):
        if r['target']=='CUMULATIVE' and r['scheme']=='STRICT' and r['axis']=='hemlock':records.append(dict(label=r['end_year']+' '+r['label']+' mortality gap',value=r['mortality_gap']))
    plot(6,'Independent adult component: retained mortality underprediction','Observed minus predicted mortality, eligible2021 survivor support',records,[src,'external_sources/HF453/risk_denominators.csv'],'All cumulative windows, primary secure endpoints, hemlock/other','Retrospective left-truncated cohort, no external refit or publicly issued forecast. DN bounds, annual windows, cluster intervals and exclusions remain available.')
    src='external_sources/HF355/transition_summary.csv';records=[dict(label=r['axis']+' '+r['label']+' survival',value=r['survival_fraction']) for r in table(src) if r['axis'] in ['ALL','end_year'] and r['survival_fraction']!='']
    plot(7,'Seedling transitions on their actual observational support','Descriptive sampled below1cm-subplot survival fraction',records,[src,'external_sources/alignment_matrix.csv'],'All and end-year strata; no whole-plot expansion','Repeated actual-date transitions, not fitted adult-model forecasts. No identified bridge to whole-plot1–10cm juveniles, germination date or postbaseline graduation date.')
    src='entry_missingness/diagnostic_tipping.csv';groups=defaultdict(list)
    for r in table(src):
        if r['theta']=='0.75':groups[(r['query'],r['scenario'],r['unknown_log_hazard_shift'],r['entry_membership_multiplier'])].append(r['finite_health']=='True')
    plot(8,'Entry and unknown-fate diagnostic tipping','Finite-design Health fraction at theta0.75',[dict(label=' / '.join(k),value=sum(v)/len(v)) for k,v in sorted(groups.items())],[src,'entry_missingness/boundary_conditioned_missingness.csv'],'All preregistered16 states and3x3 diagnostic ranges, D0/D3 and both queries','Diagnostic membership and latent-fate ranges are not posteriors. Complete five-family surfaces retain margin/Realizes/provenance and local upper tails.')
    src='corrected_pipeline/reference_crossmatrix_scores.csv';records=[]
    for r in table(src):
        if r['cohort']=='boundary' and r['suite']=='correct' and r['policy']=='B0_full' and r['panel']=='PRESENT_NEAR_0.1':records.append(dict(label=r['query']+' K'+r['reference_K']+' '+r['crossmatrix'],value=r['reference_unresolved_fraction']))
    plot(9,'Reference precision and decision-boundary uncertainty','Fraction of near-boundary reference labels unresolved',records,[src,'corrected_pipeline/reference_counts.csv'],'All four prediction/reference combinations at corrected available prefixes','Original65,536 references are separately retained; corrected65,536 is not evaluated. Unresolved labels remain in the census, not silently assigned truth.')
    src='corrected_pipeline/local_precision_contrasts.csv';records=[dict(label='case'+r['selection_index']+' K'+r['K']+' '+r['comparison'],value=r['mass_delta']) for r in table(src)]
    plot(10,'Local failure diagnosis: fixed-prefix exposure extremes','Signed target-capacity change in adverse selected physical cases',records,[src,'corrected_pipeline/local_extreme_selection.csv'],'Every available global-max attaining case, all prefixes and both individual-time laws','Large local effects persist at higher precision. Adverse selection is not representative. Distinct-law MC brackets, threshold decisions and original states remain in source data; no ecological causal error partition.')
    h.write_json('figures/index.json',dict(status='COMPLETE',figures=['H05-'+str(i) for i in range(1,11)],source_hashes_required=True,all_svg_and_data_retained=True))
if __name__=='__main__':run()
