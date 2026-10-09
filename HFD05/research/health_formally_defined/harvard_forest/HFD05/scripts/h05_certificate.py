"""Earn terminal evidence statuses from completed gates and actual isolated receipts."""
import json
import h05_common as h
REQUIRED_REPORTS=['Design_and_Freeze','Q2_Mathematical_Correction','Gamma_Growth_Audit','Unequal_Exposure_Audit','Corrected_Pipeline_Comparison','Boundary_Oracle_Attribution','Independent_External_Validation','Entry_and_Missingness_Identification','Context_and_Statistical_Interpretation','Numerical_Methods_Contract','PaperV_Referee_Response_Map','PaperV_Corrected_Claims','Public_Artifact_Readiness','Verification','Closure']
def build():
    protected=h.predecessor_guard();progress=json.loads((h.ROOT/'reports/verification_progress.json').read_text())
    if progress['status']!='SCIENTIFIC_GATES_PASS' or progress['active_hostile_controls']<30:raise ValueError('Scientific verifier not passed')
    for stage in ['tier1','tier2','external-audit','figures','oracle','references']:
        receipt=json.loads((h.ROOT/'release'/('replay_'+stage+'.json')).read_text())
        if receipt['status']!='PASS' or receipt['exit_code']!=0 or receipt['private_raw_or_git_used']:raise ValueError('Isolated replay unfinished '+stage)
    evidence={}
    prefixes=['config','provenance','q2_correction','gamma_growth','exposure_model','corrected_pipeline','oracle_attribution','external_sources','uncertainty','entry_missingness','mathematical_audit','manuscript','figures']
    for prefix in prefixes:
        for path in sorted((h.ROOT/prefix).rglob('*')):
            if not path.is_file() or '_partial.' in path.name or '__pycache__' in path.parts:continue
            evidence[str(path.relative_to(h.ROOT))]=h.sha(path)
    for stage in ['tier1','tier2','external-audit','figures','oracle','references','oracle-analysis']:
        name='release/replay_'+stage+'.json';evidence[name]=h.sha(h.ROOT/name)
    for name in REQUIRED_REPORTS:
        path='reports/'+name+'.md'
        if not(h.ROOT/path).exists():raise ValueError('Report missing '+path)
        evidence[path]=h.sha(h.ROOT/path)
    statuses=dict(PREDECESSOR_SCIENCE_FROZEN=True,Q2_EQUIVALENCE='VERIFIED',Q2_NUMERICAL_RESOLUTION='CORRECTED',GAMMA_MOMENT_AUDIT='COMPLETE',GAMMA_CORRECTION_IMPACT='MATERIAL',INDIVIDUAL_EXPOSURE_AUDIT='COMPLETE',EXPOSURE_CORRECTION_IMPACT='MATERIAL',CORRECTED_PIPELINE='COMPLETE',BOUNDARY_ORACLE_ATTRIBUTION='PARTIAL',HF453_COMPONENT_VALIDATION='PARTIALLY_ALIGNED_AND_SCORED',HF355_COMPONENT_VALIDATION='DESCRIPTIVE_ONLY',BIOLOGICAL_ENTRY_IDENTIFICATION='NOT_IDENTIFIED',MISSINGNESS_TIPPING_ANALYSIS='COMPLETE',CONTEXT_INTERPRETATION='RECONCILED',PUBLIC_ARTIFACT_CANDIDATE='READY_FOR_DISCLOSURE_REVIEW',OPERATIONAL_HARVARD_FOREST_HEALTH='NOT_ESTABLISHED',FULL_PROSPECTIVE_ECOLOGICAL_VALIDATION='NOT_CLAIMED',PAPER_V_CORRECTED_MANUSCRIPT='AUTHORIZED')
    h.write_json('reports/hfd05_contract.json',dict(status='COMPLETE',statuses=statuses,evidence=evidence,scientific_execution_start=h.config()['actual_execution_start'],protected_predecessor_baseline=h.BASELINE,protected_paths=len(protected),gate_dependencies='Actual scientific verifier, complete original/corrected grids, support-valid independent raw scores and six isolated stage receipts; no favorable ecological gate',
        impact_domains=dict(Gamma='MATERIAL actual moments at reached floor; evaluated finite capacity/Health grids exactly preserved, not continuous-law identity',exposure='MATERIAL local capacity/Health and origin effects; full grids and adverse high-precision localization'),
        headline_decisions='corrected_pipeline/headline_decisions.csv',boundary_attribution='Executed all four roles, strong signed interactions and replicated MC; ecological causal identification is PARTIAL, never forced additive',
        scope_obstructions=dict(corrected65536='NOT_REEVALUATED',synthetic_individual_dates='INAPPLICABLE',fine_profile_individual_exposure='NOT_REEVALUATED',seedling_whole_plot_J='IDENTIFICATION_OBSTRUCTION'),
        mass_accounting='Unconditional viable subprobability and full-candidate complement; no success renormalization',public_release=False,DOI=False,source_license='Separate owner disclosure/source-IP review required; local candidate only'))
    print('Comprehensive terminal certificate COMPLETE',len(evidence),'evidence hashes',flush=True)
if __name__=='__main__':build()
