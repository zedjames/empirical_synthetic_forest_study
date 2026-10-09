"""Second implementation: raw-source crosswalk, endpoint and score reconstruction."""
import csv,datetime,io,json,math,re
from collections import defaultdict
import numpy as np
import h05_common as h

def raw(name,encoding):
    sources=json.loads((h.ROOT/'external_sources/source_inventory.json').read_text())
    if isinstance(sources,dict):sources=sources.get('sources',sources)
    record=next(r for r in sources if r['name']==name)
    path=h.ROOT/record['local_input']
    if h.sha(path)!=record['sha256']:raise ValueError('External raw hash mismatch')
    return list(csv.DictReader(io.StringIO(path.read_bytes().decode(encoding))))
def canonical(value):return str(int(value)) if re.fullmatch('[0-9]+',value) else value
def check(condition,message):
    if not condition:raise ValueError(message)
def audit():
    h.predecessor_guard();_,e0,e1,core,model=h.original.load_data()
    adult=raw('hf453-01-tree-census-2021-2024.csv','latin1');seed=raw('hf355-02-seedlings.csv','ascii')
    # Independent relational join: ambiguous candidates never collapsed.
    aliases=defaultdict(set)
    for epoch in [e0,e1]:
        for stem,r in epoch.items():
            for value in [r['stem.tag'],r['tag']+'.'+r['stem.tag']]:
                if value not in ['', 'NA']:aliases[canonical(value)].add(stem)
    history=defaultdict(dict)
    for r in adult:
        key=canonical(r['StemTag']);year=int(r['year']);check(year not in history[key],'Duplicate adult visit');history[key][year]=r
    cross={r['external_tag']:r for r in h.read_csv(h.ROOT/'external_sources/HF453/tag_crosswalk.csv')}
    for key,r in cross.items():
        check(set(r['candidate_stems'].split('|'))- {''}==aliases[key],'Adult candidate census mismatch')
        if r['primary_eligible']!='True':continue
        stem=r['linked_stem'];check(aliases[key]=={stem},'Ambiguous adult tag accepted')
        check(stem in e0 and stem in e1 and e0[stem]['df.status']=='alive','Left truncation/two-epoch baseline mismatch')
        check(e0[stem]['tree.id']==e1[stem]['tree.id'] and e0[stem]['sp']==e1[stem]['sp'],'Association mismatch')
        check(float(e1[stem]['dbh'])>=10 and history[key][2021]['status'] in ['A','AU'],'Adult support mismatch')
        dead=False
        for year,item in sorted(history[key].items()):
            check(not(dead and item['status'] in ['A','AU']),'Secure reversal accepted')
            dead=dead or item['status'] in ['DC','DS']
    predictions=h.read_csv(h.ROOT/'external_sources/HF453/mortality_predictions.csv');cache={};by=defaultdict(list)
    keys=set()
    for r in predictions:
        key=(r['target'],r['start_year'],r['end_year'],r['scheme'],r['external_tag']);check(key not in keys,'Duplicate scored adult');keys.add(key)
        entry=cross[r['external_tag']];check(entry['primary_eligible']=='True','Excluded adult scored');a,b=int(r['start_year']),int(r['end_year']);start,end=history[r['external_tag']][a],history[r['external_tag']][b]
        check(start['status'] in ['A','AU'] and r['left_truncated2021']=='True','Incorrect adult risk baseline')
        status=end['status'];check(status in ['A','AU','DC','DS','DN'],'Missing status treated secure')
        check(not(r['scheme']=='STRICT' and status=='DN'),'DN treated as secure death')
        observed=int(status in ['A','AU'] or status=='DN' and r['scheme']=='DN_ALIVE')
        check(observed==int(r['observed_survival']),'Adult endpoint mismatch')
        cell=int(entry['cell']);duration=b-a
        if (cell,duration) not in cache:
            cache[(cell,duration)]=math.fsum(float(weight)*math.exp(-float(hazard)*duration) for hazard,weight in zip(model['hazard_grid'],model['hazard_posterior'][cell]))
        p=cache[(cell,duration)];y=observed
        check(abs(float(r['predicted_survival'])-p)<2e-15,'Frozen hazard prediction mismatch')
        check(abs(float(r['brier'])-(p-y)**2)<2e-15,'Brier mismatch')
        check(abs(float(r['negative_log_score'])+y*math.log(max(p,1e-15))+(1-y)*math.log(max(1-p,1e-15)))<1e-12,'Log score mismatch')
        check(r['model_refitted']=='False','External-outcome refit')
        by[key[:4]].append(r)
    for summary in h.read_csv(h.ROOT/'external_sources/HF453/mortality_scores.csv'):
        selected=by[(summary['target'],summary['start_year'],summary['end_year'],summary['scheme'])];axis,label=summary['axis'],summary['label']
        if axis!='ALL':selected=[r for r in selected if r[axis]==label]
        n=len(selected);check(n==int(summary['n']),'Adult score denominator mismatch')
        if n:
            check(sum(1-int(r['observed_survival']) for r in selected)==int(summary['observed_deaths']),'Adult deaths mismatch')
            for field in ['brier','negative_log_score']:
                check(abs(math.fsum(float(r[field]) for r in selected)/n-float(summary[field]))<1e-12,'Adult aggregate mismatch')
    for denominator in h.read_csv(h.ROOT/'external_sources/HF453/risk_denominators.csv'):
        a,b=int(denominator['start_year']),int(denominator['end_year']);eligible=[key for key,r in cross.items() if r['primary_eligible']=='True' and history[key][a]['status'] in ['A','AU']]
        check(len(eligible)==int(denominator['risk_stems']),'Risk pool mismatch')
        check(len(by[(denominator['target'],str(a),str(b),denominator['scheme'])])==int(denominator['scored_stems']),'Scored pool mismatch')
    # Derive every seedling transition from actual sampled Y/N dates.
    histories=defaultdict(list)
    for r in seed:histories[r['seedlingID']].append(r)
    seed_cross=h.read_csv(h.ROOT/'external_sources/HF355/identity_crosswalk.csv');expected={};cohorts=0;reversals=graduations=0
    for identity in seed_cross:
        if identity['identity_eligible']!='True':continue
        visits=sorted(histories[identity['seedlingID']],key=lambda r:int(r['yearOfOb']));cohorts+=1;dead=False;bad=False
        for r in visits:
            if r['sampled']=='1' and r['alive']=='N':dead=True
            if dead and r['sampled']=='1' and r['alive']=='Y':bad=True
        vg=any(r['status']=='VG' for r in visits);reversals+=bad;graduations+=vg
        if bad or vg:continue
        usable=[]
        for r in visits:
            if r['sampled']!='1' or r['alive'] not in ['Y','N']:continue
            try:date=datetime.date.fromisoformat(r['date'])
            except ValueError:continue
            if date.year==int(r['yearOfOb']):usable.append((r,date))
        for (a,da),(b,db) in zip(usable,usable[1:]):
            if a['alive']!='Y' or db<=da:continue
            expected[(identity['seedlingID'],a['date'],b['date'])]=(int(b['alive']=='Y'),(db-da).days/365.25)
    transitions=h.read_csv(h.ROOT/'external_sources/HF355/transition_records.csv');actual={}
    for r in transitions:
        key=(r['seedlingID'],r['start_date'],r['end_date']);check(key not in actual,'Duplicate seedling transition');actual[key]=(int(r['observed_survival']),float(r['exposure_years']))
        check(r['model_prediction']=='' and r['whole_plot_expansion']=='False','Invented seedling model/support')
    check(actual==expected,'Seedling actual-date transition census mismatch')
    for r in h.read_csv(h.ROOT/'external_sources/HF355/transition_summary.csv'):
        chosen=transitions
        if r['axis']=='taxon':chosen=[v for v in chosen if v['taxonCode']==r['label']]
        elif r['axis']=='end_year':chosen=[v for v in chosen if v['end_year']==r['label']]
        elif r['axis']=='UNFLAGGED_ONLY':chosen=[v for v in chosen if v['flagged']=='False']
        check(len(chosen)==int(r['n']) and sum(int(v['observed_survival']) for v in chosen)==int(r['survivors']),'Seedling summary mismatch')
    check(len(h.read_csv(h.ROOT/'external_sources/HF355/first_observation_cohorts.csv'))==cohorts,'Cohort denominator mismatch')
    h.write_json('mathematical_audit/external_independent_audit.json',dict(status='PASS',raw_adult_rows=len(adult),raw_seedling_rows=len(seed),adult_predictions=len(predictions),all_endpoint_and_score_rows_recomputed=True,seedling_transitions=len(actual),cohorts=cohorts,seedling_reversal_identities=reversals,VG_identities=graduations,source_sha256=h.sha(__file__),raw_join_and_scores_independent=True,bootstrap_separately_checked=True))
    print('Independent external audit PASS',len(predictions),len(actual),flush=True)
if __name__=='__main__':audit()
