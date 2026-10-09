"""Save actual isolated command receipts, not aspirational replay declarations."""
import json,subprocess,sys,time
import h05_common as h
def run(stage):
    candidate=h.ROOT/'release/candidate';runner=candidate/'hfd05_runner.py';args=[sys.executable,'-I',str(runner)]
    if stage in ['tier1','tier2']:args.append(stage)
    else:args+=['tier3','--stage',stage]
    previous=h.ROOT/'release'/('replay_'+stage+'.json')
    if previous.exists() and json.loads(previous.read_text())['status']=='FAIL':
        failed=json.loads(previous.read_text());h.write_json('release/replay_'+stage+'_preserved_failure.json',failed)
    started=time.monotonic();result=subprocess.run(args,cwd=candidate,capture_output=True,text=True);elapsed=time.monotonic()-started
    receipt=dict(status='PASS' if result.returncode==0 else 'FAIL',stage=stage,exit_code=result.returncode,stdout=result.stdout,stderr=result.stderr,observed_runtime_seconds=elapsed,artifact_manifest_sha256=h.sha(candidate/'HFD05_Artifact_Manifest.json'),public_inputs_sha256=h.sha(candidate/'Public_Inputs.json'),scientific_sources_package_only=True,private_raw_or_git_used=False,source_sha256=h.sha(__file__))
    h.write_json('release/replay_'+stage+'.json',receipt)
    if stage=='oracle-analysis' and result.returncode==0:
        first=json.loads((h.ROOT/'release/replay_oracle.json').read_text())
        if first['status']!='FAIL' or 'Original oracle world 78 of 81' not in first['stdout'] or 'reference_counts_partial.csv' not in first['stderr']:raise ValueError('Recovery must identify the actual completed-generation bookkeeping failure')
        h.write_json('release/replay_oracle_preserved_failure.json',first)
        # Confirm the retained generated census and state errors are exactly the
        # frozen private scientific outputs, independently of the new adapter.
        digests={}
        for name in ['draw_census.csv','state_errors.csv']:
            old=h.ROOT/'oracle_attribution'/name;new=candidate/'research/health_formally_defined/harvard_forest/HFD05/oracle_attribution'/name
            if h.sha(old)!=h.sha(new):raise ValueError('Completed generated oracle census differs')
            digests[name]=h.sha(new)
        h.write_json('release/replay_oracle.json',dict(receipt,stage='oracle',complete_generation_then_analysis_recovery=True,generation_census_sha256=digests,generation_receipt='release/replay_oracle_preserved_failure.json',analysis_receipt='release/replay_oracle-analysis.json',total_observed_runtime_seconds=first['observed_runtime_seconds']+receipt['observed_runtime_seconds']))
    print(result.stdout,end='');print(result.stderr,end='',file=sys.stderr)
    if result.returncode:raise SystemExit(result.returncode)
if __name__=='__main__':run(sys.argv[1])
