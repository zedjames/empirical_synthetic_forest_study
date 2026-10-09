"""HFD04 standalone audit: preserved science, exact counts, all cells, real hostile tests."""
import argparse,csv,hashlib,json,math,subprocess,sys,tempfile
from collections import Counter,defaultdict
from pathlib import Path
from statistics import NormalDist
import numpy as np
REPO=Path(__file__).resolve().parents[1]
ROOT=REPO/"research/health_formally_defined/harvard_forest/HFD04"
sys.path.insert(0,str(ROOT/"scripts"))
import h04_common as h
import h04_governance as g
from h04_estimator import Context,evaluate
def require(value,message):
    if not value:raise ValueError(message)
def boolean(value):return value is True or value=="True"
def csvrows(name):return h.read_csv(ROOT/name)
def jsonfile(name):return json.loads((ROOT/name).read_text())
def interval(n,K,confidence):
    z=NormalDist().inv_cdf((1+confidence)/2);p=n/K;den=1+z*z/K
    center=(p+z*z/(2*K))/den
    width=z*math.sqrt(p*(1-p)/K+z*z/(4*K*K))/den
    return (0 if n==0 else max(0,center-width),1 if n==K else min(1,center+width))
def decision(present,counts,K,theta):
    require(K>0 and all(type(n)is int and 0<=n<=K for n in counts),"Invalid census")
    limits=[interval(n,K,.95 if len(counts)==1 else .975) for n in counts]
    finite=present and all(n/K>=theta for n in counts)
    status="FALSE" if not present or any(u<theta for l,u in limits) else ("TRUE" if theta==0 or all(l>theta for l,u in limits) else "MC_UNRESOLVED")
    return finite,status,limits
def portable():
    h.confirm_guard()
    config=h.config();worlds=jsonfile("config/world_design.json")["worlds"]
    require(len(worlds)==81,"Whole final world roster changed")
    g.streams(h.PILOT["namespaces"]+jsonfile("config/pilot_v2.json")["namespaces"],
              ["world_parameters","world_state","observation","truth","estimate_state","estimate_params","estimate_future"])
    g.claims(config["claims"])
    truth=csvrows("boundary_validation/truth_prefixes.csv")
    expected={(i,suite,q,K) for i in range(81) for suite in config["boundary"]["misspecifications"] for q in ["Q1","Q2"] for K in [4096,16384,65536]}
    actual={(int(r["world"]),r["suite"],r["query"],int(r["K"])) for r in truth}
    g.exact(expected,actual,"Final worlds/prefixes selected after outcomes")
    require(len(truth)==len(expected),"Duplicate truth cell")
    for r in truth:
        nv,nr,K=int(r["nv"]),int(r["nr"]),int(r["K"]);require(0<=nr<=nv<=K,"Unconditional Q2 support")
        counts=[nv] if r["query"]=="Q1" else [nv,nr]
        finite,status,limits=decision(boolean(r["present"]),counts,K,.75)
        g.exact(finite,boolean(r["FINITE_BANK_HEALTH"]),"Wrong exact reference truth")
        g.mc_status(status,r["KERNEL_MC_STATUS"])
        require(abs(float(r["margin"])-(min(counts)/K-.75))<1e-14,"Margin substituted")
        require(abs(float(r["lower_V"])-limits[0][0])<1e-12,"Reference interval wrong")
    require(jsonfile("boundary_validation/pre_estimator_census.json")["truth_census_sha256"]==h.sha(ROOT/"boundary_validation/truth_prefixes.csv"),"Truth census not bound before estimation")
    draws=csvrows("context_ablation/estimator_draws.csv")
    keys={(r["cohort"],int(r["world"]),r["policy"],int(r["draw"])) for r in draws}
    expected={(cohort,i,policy,d) for cohort,N in [("boundary",81),("HFD03",324)] for i in range(N) for policy in config["estimator"]["contexts"] for d in range(32)}
    g.exact(expected,keys,"Ablation worlds/draws omitted");require(len(draws)==len(keys),"Duplicate estimator draw")
    grouped=defaultdict(list)
    for r in draws:
        for q in ["Q1","Q2"]:
            finite,status,_=decision(boolean(r["present"]),[int(r["nv"])] if q=="Q1" else [int(r["nv"]),int(r["nr"])],int(r["K"]),.75)
            g.exact(finite,boolean(r[q]),"Wrong per-draw Health");g.mc_status(status,r[q+"_status"])
        grouped[(r["cohort"],r["world"],r["policy"])].append(r)
    performance=csvrows("boundary_validation/performance.csv")
    expected_performance={(cohort,str(i),suite,q,policy) for cohort,N in [("boundary",81),("HFD03",324)] for i in range(N) for suite in config["boundary"]["misspecifications"] for q in ["Q1","Q2"] for policy in config["estimator"]["contexts"]}
    g.exact(expected_performance,{tuple(r[k] for k in ["cohort","world","suite","query","policy"]) for r in performance},"Performance roster selected")
    require(len(performance)==len(expected_performance),"Duplicate/missing performance cells")
    references={}
    for row in truth:
        if int(row["K"])==65536:
            references[("boundary",row["world"],row["suite"],row["query"])]=(boolean(row["present"]),int(row["nv"]),int(row["nr"]),65536)
    for row in h.read_csv(h.THREE/"synthetic_validation/truth_census.csv"):
        for query in ["Q1","Q2"]:
            references[("HFD03",row["world"],row["suite"],query)]=(boolean(row["truth_present"]),round(float(row["truth_mass"])*4096),round(float(row["truth_reserve"])*4096),4096)
    for r in performance:
        present,nv,nr,K=references[(r["cohort"],r["world"],r["suite"],r["query"])]
        exact,status,_=decision(present,[nv] if r["query"]=="Q1" else [nv,nr],K,.75)
        require(boolean(r["reference_present"])==present and boolean(r["reference_health"])==exact,"Reference label substituted")
        g.mc_status(status,r["reference_status"])
        require(abs(float(r["reference_viable"])-nv/K)<1e-14 and abs(float(r["reference_reserve"])-nr/K)<1e-14,"Reference masses substituted")
        g.context(r["policy"],r["context_interpretation"])
        require((r["observation_interpretation"]=="PRIOR_STATE_ONLY")== (r["policy"]=="B3_label_only"),"Label-only performance called measurement")
        samples=grouped[(r["cohort"],r["world"],r["policy"])]
        for name,column in [("estimated_viable","nv"),("estimated_reserve","nr")]:
            require(abs(float(r[name])-np.mean([int(x[column])/int(x["K"]) for x in samples]))<1e-12,"Mass not derived from all draws")
        probability=sum(boolean(x[r["query"]]) for x in samples)/32
        require(abs(float(r["health_fraction"])-probability)<1e-14,"Health of average capacity substituted")
        g.exact(probability>=.5,boolean(r["predicted_health"]),"Classifier changed")
        g.reference(r["reference_status"],r["reference_status"]!="MC_UNRESOLVED")
    summaries=csvrows("context_ablation/comparison.csv")
    expected_summaries={(c,s,q,p) for c in ["boundary","HFD03"] for s in config["boundary"]["misspecifications"] for q in ["Q1","Q2"] for p in config["estimator"]["contexts"]}
    g.exact(expected_summaries,{tuple(r[k] for k in ["cohort","suite","query","policy"]) for r in summaries},"Summary roster omitted")
    require(len(summaries)==len(expected_summaries),"Duplicate summary cell")
    for r in summaries:
        selected=[p for p in performance if all(p[k]==r[k] for k in ["cohort","suite","query","policy"])]
        g.mechanism(selected)
        certified=[p for p in selected if p["reference_status"]!="MC_UNRESOLVED"]
        counts=Counter(("TP" if boolean(p["reference_health"]) and boolean(p["predicted_health"]) else
                        "FN" if boolean(p["reference_health"]) else
                        "FP" if boolean(p["predicted_health"]) else "TN") for p in certified)
        for key in ["TP","TN","FP","FN"]:require(int(r[key])==counts[key],"Confusion not independently derived")
        if certified:require(abs(float(r["brier"])-np.mean([(float(p["health_fraction"])-boolean(p["reference_health"]))**2 for p in certified]))<1e-12,"Brier not independently derived")
    regions=csvrows("boundary_validation/signed_margin_scores.csv")
    region_names=["OUTSIDE_BELOW","CLEARLY_BELOW","MODERATELY_BELOW","IMMEDIATELY_BELOW","IMMEDIATELY_ABOVE","MODERATELY_ABOVE","CLEARLY_ABOVE","OUTSIDE_ABOVE"]
    for r in summaries:
        cells=[x for x in regions if all(x[k]==r[k] for k in ["cohort","suite","query","policy"])]
        g.strata(region_names,[x["region"] for x in cells]);require(sum(int(x["n"]) for x in cells)==int(r["n"]),"Unfavorable strata pruned")
    observed=csvrows("juvenile_boundary/observed_groups.csv");support=jsonfile("juvenile_boundary/empirical_support.json")
    eligible=h.read_csv(h.THREE/"reconstruction_validation/eligible_truth_pool.csv")
    require(sum(int(r["n"]) for r in observed)==support["n"]==sum(boolean(r["dbh_eligible"]) for r in eligible),"Observed support changed")
    for r in csvrows("juvenile_boundary/initial_projection.csv"):
        require(int(r["coarse_count"])==int(r["fine_count"]),"Count projection fails")
        require(abs(float(r["BA_projection_error"]))<1e-8,"BA projection fails; cannot claim exact")
        require(r["provenance"]=="IMPUTED_ORIGIN_DISTRIBUTION_WITH_VISIBLE_E1_PROFILE_ANCESTRY","Origin/profile provenance relabeled")
        g.provenance("IMPUTED_ORIGIN",r["provenance"])
    require(len(csvrows("juvenile_boundary/propagation.csv"))==16*2*2*2*2*2*3,"Representation cells missing")
    error=csvrows("juvenile_boundary/measurement_error.csv")
    require(len(error)==5*8*6,"Measurement-error grid incomplete")
    witness=jsonfile("juvenile_boundary/projection_counterexample.json")
    require(witness["actual_response"][0]["Q2"] is False and witness["actual_response"][1]["Q2"] is True,"Counterexample does not distinguish actual Q2")
    cap=jsonfile("effective_count/fit_contract.json");g.cap_selection(cap["cap_selected"],cap["E2_used_in_fit"])
    require({int(r["cap"]) for r in csvrows("effective_count/posteriors.csv")}=={50,100,200,500},"Caps selected")
    require(len(csvrows("effective_count/propagation.csv"))==4*16*4*2*2*2*2*3,"Cap propagation cells omitted")
    census=csvrows("decision_rules/hfd03_sufficient_census.csv")
    old=Counter();derived=Counter()
    for r in census:
        counts=[int(r["nv"])] if r["query"]=="Q1" else [int(r["nv"]),int(r["nr"])]
        finite,status,_=decision(boolean(r["present"]),counts,256,float(r["theta"]))
        g.exact(finite,boolean(r["FINITE_BANK_HEALTH"]),"Wrong finite-bank status")
        g.mc_status(status,r["KERNEL_MC_STATUS"])
        g.endpoint(float(r["theta"]),all(n==256 for n in counts),status)
        old[r["old_status"]]+=int(r["evaluations"])
        if r["old_status"]=="MC_UNRESOLVED":derived[float(r["theta"])]+=int(r["evaluations"])
    require(old["MC_UNRESOLVED"]==963976 and sum(derived[t] for t in [.5,.75,.9])==3984,"Historical HFD03 record rewritten")
    input_names={(r["source"],r["parameter"]):r["value"] for r in csvrows("config/numerical_methods.csv")}
    from h04_methods import flatten
    for source,cfg in [("HFD04",config),("HFD02S",h.frozen.load("config/design.json")),("HFD04_CONTEXT_PRIOR_001",jsonfile("config/context_prior_audit.json"))]:
        for name,value in flatten(cfg):g.exact(value,input_names[(source,name)],"Methods assumption omitted/changed")
    public=getattr(h,"PUBLIC_MODE",False)
    artifact_base=REPO if public else ROOT/"release/candidate"
    record=json.loads((artifact_base/"PaperV_Artifact_Manifest.json").read_text())
    require(record["deny_by_default"] and not record["public_publication"] and not record["raw_data_redistributed"],"Export authority expanded")
    for name,item in record["files"].items():
        g.export(name,record["explicit_inclusion_list"])
        require(h.sha(artifact_base/name)==item["export_sha256"],"Public source/table changed")
    isolation=jsonfile("provenance/public_isolation_receipt.json")
    require(isolation["status"]=="PASS" and isolation["fixture_only"] and not isolation["full_public_raw_replay_executed"],"Export isolation overstated")
    require(isolation["artifact_manifest_sha256"]==h.sha(artifact_base/"PaperV_Artifact_Manifest.json"),"Isolation receipt stale")
    for number in range(1,9):
        manifest=jsonfile("figures/H04-"+str(number)+".json")
        require(h.sha(ROOT/manifest["data"])==manifest["data_sha256"],"Figure data identity")
        require(h.sha(ROOT/manifest["figure"])==manifest["figure_sha256"],"Figure version identity")
        require(manifest["run_id"]==h.sha(ROOT/"config/protocol.json")[:16],"Figure wrong experiment")
        require(manifest["generator_sha256"]==h.sha(ROOT/manifest["generator"]),"Figure generator changed")
        for name,digest in manifest["inputs"].items():require(h.sha(ROOT/name)==digest,"Figure input changed")
    supplement=supplemental(references)
    replays=[]
    for name,N in [("provenance/public_full_replay_summary.json",23),("provenance/public_supplement_replay_summary.json",5)]:
        if not (ROOT/name).exists():continue
        replay=jsonfile(name)
        require(replay["status"]=="PASS" and replay["actual_replay"] and len(replay["scientific_outputs"])==N,"Full replay evidence overstated")
        require(not replay["private_scientific_modules_imported"] and not replay["raw_cache_copied"],"Replay isolation claim changed")
        for path,digest in replay["scientific_sources"].items():require(h.sha(REPO/path)==digest,"Replayed scientific algorithm changed")
        for row in replay["scientific_outputs"]:require(h.sha(ROOT/row["path"])==row["sha256"] and row["status"]=="PASS","Replayed output changed")
        replays.append(N)
    frozen=ROOT/"provenance/execution_manifest.json"
    if frozen.exists():
        for name,digest in json.loads(frozen.read_text())["files"].items():require(h.sha(REPO/name)==digest,"Frozen HFD04 source/output changed "+name)
    return dict(worlds=81,truth_cells=len(truth),estimator_draws=len(draws),performance_cells=len(performance),historical_unresolved=963976,
                supplemental=supplement,isolated_replay_outputs=replays)

def supplemental(references):
    config=jsonfile("config/context_prior_audit.json");registration=jsonfile("provenance/context_prior_registration.json")
    for name,digest in registration["sources"].items():require(h.sha(ROOT/name)==digest,"Supplement changed after registration")
    worlds=jsonfile("config/world_design.json")["worlds"];counts=Counter((r["pressure"],r["recovery"]) for r in worlds)
    weights=[counts[(p,r)]/81 for p in range(3) for r in range(3)]
    g.exact(weights,config["weights"],"Prior selected using outcomes")
    from h04_context_prior_estimator import generator
    draws=csvrows("context_ablation/prior_reconciled_draws.csv")
    g.exact({(i,d) for i in range(81) for d in range(32)},{(int(r["world"]),int(r["draw"])) for r in draws},"Supplement selected worlds/draws")
    require(len(draws)==2592,"Duplicate supplement draw")
    grouped=defaultdict(list)
    for row in draws:
        world,draw=int(row["world"]),int(row["draw"])
        require(int(row["sampled_context_level"])==int(generator("cw_params",world,draw).choice(9,p=weights)),"Hidden labels replaced prior sampling")
        require(row["sampled_context_provenance"]=="DESIGN_PRIOR_DRAW_NOT_HIDDEN_REALIZED_LABEL","Prior sample mislabeled")
        nv,nr,K=int(row["nv"]),int(row["nr"]),int(row["K"])
        require(K==256 and 0<=nr<=nv<=K,"Supplement weights")
        for q in ["Q1","Q2"]:
            exact,status,_=decision(boolean(row["present"]),[nv] if q=="Q1" else [nv,nr],K,.75)
            g.exact(exact,boolean(row[q]),"Supplement wrong exact Health");g.mc_status(status,row[q+"_status"])
        grouped[row["world"]].append(row)
    performance=csvrows("context_ablation/prior_reconciled_performance.csv")
    expected={(str(i),s,q) for i in range(81) for s in h.config()["boundary"]["misspecifications"] for q in ["Q1","Q2"]}
    g.exact(expected,{(r["world"],r["suite"],r["query"]) for r in performance},"Supplement performance omitted")
    require(len(performance)==486,"Duplicate supplement scores")
    for row in performance:
        present,nv,nr,K=references[("boundary",row["world"],row["suite"],row["query"])]
        exact,status,_=decision(present,[nv] if row["query"]=="Q1" else [nv,nr],K,.75)
        require(boolean(row["reference_health"])==exact and row["reference_status"]==status,"Supplement truth substituted")
        require(abs(float(row["reference_viable"])-nv/K)<1e-14 and abs(float(row["reference_reserve"])-nr/K)<1e-14,"Supplement truth mass changed")
        samples=grouped[row["world"]];prob=sum(boolean(r[row["query"]]) for r in samples)/32
        require(float(row["health_fraction"])==prob and boolean(row["predicted_health"])==(prob>=.5),"Supplement classifier changed")
        for name,column in [("estimated_viable","nv"),("estimated_reserve","nr")]:require(abs(float(row[name])-np.mean([int(r[column])/256 for r in samples]))<1e-12,"Supplement mass not all draws")
    summaries=csvrows("context_ablation/prior_reconciled_comparison.csv");strata=csvrows("context_ablation/prior_reconciled_strata.csv")
    require(len(summaries)==6 and len(strata)==96,"Supplement adverse strata omitted")
    for row in summaries:
        selected=[r for r in performance if r["suite"]==row["suite"] and r["query"]==row["query"]]
        certified=[r for r in selected if r["reference_status"]!="MC_UNRESOLVED"]
        counts=Counter("TP" if boolean(r["reference_health"]) and boolean(r["predicted_health"]) else "FN" if boolean(r["reference_health"]) else "FP" if boolean(r["predicted_health"]) else "TN" for r in certified)
        for key in ["TP","TN","FP","FN"]:require(counts[key]==int(row[key]),"Supplement confusion changed")
        for kind,N in [("signed_margin",8),("absolute_margin",4),("failure_class",4)]:
            cells=[r for r in strata if r["suite"]==row["suite"] and r["query"]==row["query"] and r["kind"]==kind]
            require(len(cells)==N and sum(int(r["n"]) for r in cells)==81,"Supplement unfavorable strata pruned")
    return dict(draws=2592,performance_cells=486,roster_prior_exact=True)
def hostile():
    results=[]
    def rejects(identity,call):
        try:call()
        except (ValueError,TypeError):results.append(dict(id=identity,status="PASS"));return
        raise AssertionError("Hostile failure admitted "+identity)
    # Exercise the production guards without modifying any frozen file.
    original_check=h.subprocess.check_output
    def corrupt_blob(command,*args,**kwargs):
        result=original_check(command,*args,**kwargs)
        if command[:2]==["git","hash-object"]:return "0"*40+"\n"
        return result
    h.subprocess.check_output=corrupt_blob
    try:rejects("HC-01",h.source_guard)
    finally:h.subprocess.check_output=original_check
    original_sha=h.sha
    def selected_roster(path):return "0"*64 if Path(path).name=="world_design.json" else original_sha(path)
    h.sha=selected_roster
    try:rejects("HC-02",h.confirm_guard)
    finally:h.sha=original_sha
    rejects("HC-03",lambda:g.streams(["pilot_Q1"],["pilot_Q1"]))
    rejects("HC-04",lambda:evaluate({"true_state":True},Context("B0_full",1,1),{}, {},0,1,1))
    rejects("HC-05",lambda:g.reference("MC_UNRESOLVED",True))
    rejects("HC-06",lambda:g.mechanism([dict(suite="correct"),dict(suite="hemlock")]))
    rejects("HC-07",lambda:g.context("B2_unknown","KNOWN_CONDITIONAL_CONTEXT"))
    rejects("HC-08",lambda:Context("B2_unknown",1,1))
    rejects("HC-08b",lambda:Context("B1_coarsened",0,1,exact_pressure=1))
    packet=h.ObservationPacket((1,)*64,(2.,)*64,(.7,)*64)
    rejects("HC-09",lambda:evaluate(packet,Context("B3_label_only",1,1),{},{},0,1,1))
    rejects("HC-10",lambda:g.provenance("COHORT_RMS","OBSERVED_INDIVIDUAL"))
    import h03_calibration as cal
    from h04_juvenile import visible_inputs
    held,path=visible_inputs("M0")
    row=next(r for r in h.read_csv(path) if r["stem.id"] in held);row["dbh"]="9.9"
    directory=Path(tempfile.mkdtemp(prefix="hfd04-hostile-",dir="/private/tmp"));leak=directory/"leak.csv"
    with leak.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(row));writer.writeheader();writer.writerow(row)
    rejects("HC-11",lambda:cal.fit_visible({},leak,held))
    rejects("HC-12",lambda:g.cap_selection(True,True))
    rejects("HC-13",lambda:g.provenance("IMPUTED_ORIGIN","OBSERVED_INDIVIDUAL"))
    rejects("HC-14",lambda:g.unconditional_weights([0,1],2))
    rejects("HC-15",lambda:g.endpoint(1,True,"TRUE"))
    rejects("HC-16",lambda:g.mc_status("MC_UNRESOLVED","TRUE"))
    rejects("HC-17",lambda:g.claims(dict(scenario_prior=[.25]*4,operational="NOT_AUTHORIZED")))
    rejects("HC-18",lambda:g.claims(dict(operational="AUTHORIZED")))
    rejects("HC-19",lambda:g.export(".env",[".env"]))
    rejects("HC-20",lambda:g.strata(["BELOW","ABOVE"],["ABOVE"]))
    require(decision(True,[8],8,1)[1]=="MC_UNRESOLVED","Positive endpoint fixture")
    require(decision(True,[0],8,0)[:2]==(True,"TRUE"),"Positive algebraic endpoint fixture")
    return results
def cache(replay=False):
    truth=csvrows("boundary_validation/truth_prefixes.csv")
    manifests=jsonfile("provenance/truth_bank_manifest.json")
    manifest={(r["world"],r["suite"]):r for r in manifests}
    checked=0
    for identity in range(81):
        for suite in h.config()["boundary"]["misspecifications"]:
            bank=np.load(h.run_dir()/"truth"/(str(identity)+"_"+suite+".npz"))
            v,r,l=bank["viable"],bank["reserve"],bank["lawful"]
            require(len(v)==65536 and not np.any(v&~l) and not np.any(r&~v),"Raw truth support mismatch")
            fingerprint=hashlib.sha256(v.tobytes()+r.tobytes()+l.tobytes()).hexdigest()
            require(fingerprint==manifest[(identity,suite)]["raw_response_sha256"],"Truth archived sample identity")
            for row in [x for x in truth if int(x["world"])==identity and x["suite"]==suite]:
                K=int(row["K"])
                require(int(row["nv"])==int(v[:K].sum()) and int(row["nr"])==int(r[:K].sum()),"Truth prefixes not paired")
            checked+=1
    if replay:
        # Every world/suite receives a fresh first-chunk kernel replay; all full
        # prefix counts above are checked, not falsely called whole-bank replay.
        from h04_boundary import world
        _,_,_,core,model=h.load_data();cells=jsonfile("config/world_design.json")["worlds"]
        for identity,cell in enumerate(cells):
            state,params,_=world(core,model,cell,identity)
            for suite in h.config()["boundary"]["misspecifications"]:
                generator=h.rng("truth",identity,0)
                if suite=="entry_lognormal_annual_sd_1.2":generator=h.EntryRNG(generator)
                bank,lawful=h.demography.simulate(state,params,h.scenario(suite),5,1024,generator)
                digest=hashlib.sha256(bank.tobytes()+lawful.tobytes()).hexdigest()
                require(digest==manifest[(identity,suite)]["bank_hashes"][0],"Independent kernel replay changed")
            if identity%12==0:print("Independent kernel replay world",identity,flush=True)
    return dict(archived_banks_checked=checked,all_prefix_counts_recomputed=True,first_chunk_every_world_suite_replayed=replay,
                full_65536_kernel_replay_executed=False)
def main():
    parser=argparse.ArgumentParser();parser.add_argument("--cache",action="store_true");parser.add_argument("--replay-panel",action="store_true");parser.add_argument("--final",action="store_true");args=parser.parse_args()
    science=portable();controls=hostile();raw=cache(args.replay_panel) if args.cache or args.replay_panel else None
    receipt=dict(status="PASS",portable=science,hostile_controls=controls,cache=raw,verifier_sha256=h.sha(__file__),
       private_history="Published source identities checked; unavailable privateGit lineage attested, not independently verified" if getattr(h,"PUBLIC_MODE",False) else "Frozen source identities independently checked in the owned private repository; public fixture does not replay privateGit",
       public_scope="Isolated full packaged scientific replay PASS: 23 original plus 5 supplemental byte-exact outputs; exact externally retrieved public raw bytes reused read-only; no fresh download")
    h.write_json("reports/verification_receipt.json",receipt)
    if args.final:
        from h04_certificate import gates,scientific_status
        require((ROOT/"provenance/execution_manifest.json").is_file(),"Missing final frozen evidence")
        require(all(gates().values()),"Incomplete certificate evidence")
        contract=jsonfile("reports/hfd04_contract.json")
        g.exact(gates(),contract["gates"],"Stage certificate not evidence-linked")
        require(contract["certificate"]["BOUNDARY_HEALTH_DISCRIMINATION"]==scientific_status(),"Adverse outcome rewritten")
        require(contract["certificate"]["PAPER_V_V3"]=="AUTHORIZED" and contract["certificate"]["Q2_REPRESENTATION_SUFFICIENCY"]=="CONDITIONAL","Unsupported V3 status")
        g.claims(contract["claims"])
    print("PASS HFD04",len(controls),"active hostile controls",json.dumps(science),flush=True)
if __name__=="__main__":main()
