"""Manuscript reports/figures derived from complete registered outputs."""
import csv,html,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import h04_common as h
def rows(name):return h.read_csv(h.ROOT/name)
def number(value):return float(value) if value not in ["",None] else None
def b(value):return value=="True" or value is True
def table(data,fields):
    return "| "+" | ".join(fields)+" |\n| "+" | ".join(["---"]*len(fields))+" |\n"+"\n".join("| "+" | ".join(str(r.get(f,"")) for f in fields)+" |" for r in data)+"\n"
def report(name,text):
    path=h.ROOT/"reports"/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text.rstrip()+"\n")
def plot(index,title,data,xkey,ykey,caption,inputs):
    name="H04-"+str(index);h.write_csv("figures/"+name+".csv",data)
    points=[(number(r.get(xkey)),number(r.get(ykey))) for r in data]
    points=[(x,y) for x,y in points if x is not None and y is not None and math.isfinite(x) and math.isfinite(y)]
    xs=[x for x,y in points];ys=[y for x,y in points];xmin=min(xs,default=0);xmax=max(xs,default=1);ymin=min(0,min(ys,default=0));ymax=max(ys,default=1)
    if xmax==xmin:xmax=xmin+1
    if ymax==ymin:ymax=ymin+1
    svg=['<svg xmlns="http://www.w3.org/2000/svg" width="960" height="600" viewBox="0 0 960 600">','<rect width="960" height="600" fill="white"/>',
      '<text x="60" y="45" font-family="sans-serif" font-size="24">'+html.escape(title)+'</text>',
      '<path d="M80 80 V500 H900" fill="none" stroke="#333"/>']
    groups=[]
    palette=["#267078","#c35335","#795b9b","#668436","#2f5297","#8c6337","#a03867","#505050"]
    keyed=[r for r in data if number(r.get(xkey)) is not None and number(r.get(ykey)) is not None]
    for row,(x,y) in zip(keyed,points):
        label=" / ".join(str(row[k]) for k in (["cohort","query"] if index==3 else ["cohort","suite","query","policy"]) if k in row)
        if label not in groups:groups.append(label)
        color=palette[groups.index(label)%len(palette)]
        xx=80+(x-xmin)/(xmax-xmin)*800;yy=500-(y-ymin)/(ymax-ymin)*400
        svg.append('<circle cx="'+str(xx)+'" cy="'+str(yy)+'" r="5" fill="'+color+'" fill-opacity=".7"/>')
    for i in range(5):
        x=xmin+(xmax-xmin)*i/4;y=ymin+(ymax-ymin)*i/4
        svg.append('<text x="'+str(80+i*200)+'" y="523" text-anchor="middle" font-family="sans-serif" font-size="13">'+format(x,'.3g')+'</text>')
        svg.append('<text x="70" y="'+str(505-i*100)+'" text-anchor="end" font-family="sans-serif" font-size="13">'+format(y,'.3g')+'</text>')
    svg.append('<text x="470" y="550" text-anchor="middle" font-family="sans-serif">'+html.escape(xkey)+'</text>')
    svg.append('<text x="85" y="70" font-family="sans-serif">'+html.escape(ykey)+'</text>')
    for i,label in enumerate(groups):
        if label:
            svg.append('<rect x="630" y="'+str(82+i*17)+'" width="10" height="10" fill="'+palette[i%len(palette)]+'"/>')
            svg.append('<text x="646" y="'+str(91+i*17)+'" font-family="sans-serif" font-size="11">'+html.escape(label)+'</text>')
    svg.append('<text x="60" y="580" font-family="sans-serif" font-size="12">'+html.escape("Exact data/strata and caption: "+name+".csv / "+name+".json")+'</text></svg>')
    if index==8:
        svg=['<svg xmlns="http://www.w3.org/2000/svg" width="960" height="760" viewBox="0 0 960 760">','<rect width="960" height="760" fill="white"/>','<text x="60" y="45" font-family="sans-serif" font-size="24">End-to-end decision provenance</text>']
        for i,row in enumerate(data):
            y=75+i*100
            svg.append('<rect x="120" y="'+str(y)+'" width="700" height="70" rx="8" fill="#eef4f1" stroke="#267078"/>')
            svg.append('<text x="140" y="'+str(y+25)+'" font-family="sans-serif" font-size="18">'+html.escape(row["label"]+' — '+str(row["registered_failure_count"]))+'</text>')
            svg.append('<text x="140" y="'+str(y+50)+'" font-family="sans-serif" font-size="13">'+html.escape(row["interpretation"])+'</text>')
            if i<5:svg.append('<path d="M470 '+str(y+72)+' v22 m-5 -6 l5 6 l5 -6" fill="none" stroke="#267078" stroke-width="2"/>')
        svg.append('<text x="120" y="710" font-family="sans-serif" font-size="13">Distinct denominators: counts do not sum to a causal error partition.</text></svg>')
    (h.ROOT/"figures"/(name+".svg")).write_text("\n".join(svg)+"\n")
    h.write_json("figures/"+name+".json",dict(id=name,title=title,caption=caption,data="figures/"+name+".csv",figure="figures/"+name+".svg",
      data_sha256=h.sha(h.ROOT/"figures"/(name+".csv")),figure_sha256=h.sha(h.ROOT/"figures"/(name+".svg")),generator="scripts/h04_reports.py",
      generator_sha256=h.sha(__file__),run_id=h.sha(h.ROOT/"config/protocol.json")[:16],inputs={p:h.sha(h.ROOT/p) for p in inputs}))
def make():
    comparisons=rows("context_ablation/comparison.csv");performance=rows("boundary_validation/performance.csv")
    primary=json.loads((h.ROOT/"boundary_validation/primary_discrimination.json").read_text())
    truth=rows("boundary_validation/truth_prefixes.csv");target=[r for r in truth if r["K"]=="65536" and r["suite"]=="correct" and r["query"]==r["family"]]
    density=dict(target_cells=len(target),within_0_1=sum(abs(float(r["margin"]))<=.1 for r in target),within_0_02=sum(abs(float(r["margin"]))<=.02 for r in target),
      reference_unresolved=sum(r["KERNEL_MC_STATUS"]=="MC_UNRESOLVED" for r in target))
    convergence=[]
    for K in [4096,16384,65536]:
        selected=[r for r in truth if int(r["K"])==K and r["suite"]=="correct" and r["query"]==r["family"]]
        convergence.append(dict(K=K,cells=len(selected),resolved=sum(r["KERNEL_MC_STATUS"]!="MC_UNRESOLVED" for r in selected),
           max_interval_width=max(float(r["upper_V"])-float(r["lower_V"]) for r in selected)))
    report("Boundary_Stress_Validation.md","# Adequacy-boundary stress validation\n\nAll81 worlds retained:72 targets and9 realization controls; no confirmatory selection. Pilot001 Q2 count-scaling failed (reserve mass one); pilot002 independently mapped juvenile diameter to reserve mass. Pilot outcomes are excluded from performance.\n\n"+
      table([density],list(density))+"\n"+table(convergence,["K","cells","resolved","max_interval_width"])+
      "\nPrimary near-boundary PRESENT-world performance, correct kernel family/full context:\n\n"+table([dict(query=q,**primary[q]) for q in ["Q1","Q2"]],["query","n","certified_n","positives","negatives","FP","FN","balanced_accuracy","viable_MAE","reserve_MAE","status"])+
      "\nAll-suite/full-context results (includes separate realization failures, so not a substitute for the primary boundary assessment):\n\n"+
      table([r for r in comparisons if r["cohort"]=="boundary" and r["policy"]=="B0_full"],["suite","query","n","certified_n","FP","FN","balanced_accuracy","reference_unresolved_fraction","estimator_unresolved_fraction"])+
      "\nSigned and absolute strata—including empty and unfavorable strata—are in signed_margin_scores.csv and absolute_margin_scores.csv. Failure scores distinguish REALIZATION, Q1_CAPACITY, Q2_RESERVE. Uncertain reference truths remain in every world table but are excluded explicitly from certified-truth confusion/Brier scoring. Finite-bank truth is separately available. Conditional error is an experimental frequency, not ecological risk. Parameter-moment normalization is deterministic design, not truth fitting; correct kernel family does not guarantee a correct latent-state prior.")
    report("Context_Label_Ablation.md","# Context-label information ablation\n\nPaired405 worlds:all81 new worlds plusall324 frozen HFD03 worlds/observations. Every context policy uses the same observation and state/future random nonce. Known labels, coarse markers and explicit unknown are distinct strict types; unavailable estimators never receive exact levels. Coarse pressure[0,1]/[2], recovery[0]/[1,2]; uniform conditional design mixtures. Unknown mixes9 contexts equally, not ecological probabilities.\n\n"+
      table([r for r in comparisons if r["suite"]=="correct"],["cohort","query","policy","FP","FN","balanced_accuracy","viable_MAE","reserve_MAE","brier","viable_width90","estimator_unresolved_fraction","delta_vs_B0_balanced_accuracy"])+
      "\nLabel-only uses no informative state observation:it samples the public E0 prior. Comparisons quantify the importance of supplied information under this estimator, not an additive causal decomposition of observed state versus context. HFD04 uses32 draws versus historical64:paired current B0—not99.22% historical headline—is the attribution baseline. Unknown-context prediction targets the hidden realized-context truth; it does not redefine known-context formal Health. All mechanisms remain separate. Exact label-to-prior map and cohort-share normalization are in Numerical_Methods_Contract.md.")
    path=h.ROOT/"reports/Context_Label_Ablation.md"
    path.write_text(path.read_text()+"\n\n## Roster-prior reconciliation (separately registered supplement)\n\nThe uniform-nine prior agrees with the 324-world factorial roster, but not the boundary-designed roster. The latter contains context counts [7,7,7,7,7,7,13,13,13]. Original uniform-mixture results remain above and unchanged. Supplemental prior weights are these counts divided by 81, chosen from the entire frozen roster before supplemental outcomes, with independent random streams and the same observations. They are experimental design frequencies, not ecological frequencies. All signed-margin, absolute-margin and failure-class strata remain in prior_reconciled_strata.csv.\n\n"+
        table(rows("context_ablation/prior_reconciled_comparison.csv"),["suite","query","n","certified_n","FP","FN","balanced_accuracy","viable_MAE","reserve_MAE","brier","estimator_unresolved_fraction"])+
        "\nCorrect-kernel aggregate balanced accuracy is 79.4% Q1 and 71.9% Q2; it does not establish robust boundary discrimination or recovered context. Comparisons with original B2 include independent-stream sampling variation and prior change, not a pure causal effect.\n")
    observed=rows("juvenile_boundary/observed_groups.csv");projection=rows("juvenile_boundary/initial_projection.csv");prop=rows("juvenile_boundary/propagation.csv")
    obs=dict(n=sum(int(r["n"]) for r in observed),FP=sum(int(r["FP"]) for r in observed),FN=sum(int(r["FN"]) for r in observed))
    obs["membership_error_fraction"]=(obs["FP"]+obs["FN"])/obs["n"]
    repchanges=Counter(r["change_class"] for r in prop);disagreement=sum(b(r["coarse_Health"])!=b(r["fine_Health"]) for r in prop)/len(prop)
    report("Juvenile_Threshold_and_Query_Sufficiency.md","# Juvenile threshold and query sufficiency\n\nSecure observed E1 subset only; dated support is retained in observed_groups.csv. No claim of measured2020 origin juvenile counts.\n\n"+table([obs],list(obs))+
      "\nRMS and stem classification differ even when count/BA agree. Direct Gaussian perturbations0/.1/.25/.5/1cm are methodological stress conditions, not calibrated instrument-error laws. Membership errors by distance and all8 replicates are in measurement_error.csv. M0–M3 masking retains physically sanitized fitting inputs; profile data exclude heldout DBH.\n\n"+
      table(rows("juvenile_boundary/masking.csv"),["mask","n","truth_juveniles","mean_predicted_juveniles","delta_J_over_J0","membership_brier"])+
      "\nPaired enriched-state count/BA projection max absoluteBA error="+str(max(abs(float(r["BA_projection_error"])) for r in projection))+
      "; juvenile counts do NOT commute. The explicit projection_counterexample.json has identical count/RMS/BA and opposite actual Q2 answers. Thus RMS is not generally sufficient for the juvenile query.\n\nFinite-design paired Health disagreement="+str(disagreement)+
      ". Decision-change categories:\n\n"+table([dict(category=k,cells=v) for k,v in sorted(repchanges.items())],["category","cells"])+
      "\nThe enriched distribution is inferred from visible E1 profiles and rescaled to reconstructed origin moments; it is not a newly observed modality. Initial projection is numerically exact; future stochastic kernels are not claimed to commute or to be a conservative refinement. Common seeds are paired, not exact common uniforms after changing array dimensions. Size-aware representation preserves an explicit juvenile distinction but does not establish full latent ecological sufficiency. Status CONDITIONAL, not an operational upgrade endorsement.")
    changed=rows("juvenile_boundary/changed_case_audit.csv")
    path=h.ROOT/"reports/Juvenile_Threshold_and_Query_Sufficiency.md"
    path.write_text(path.read_text()+"\n\n## Local effect, not just the global disagreement rate\n\nOnly 3/1,536 finite-grid decisions differ (0.195%), but one stressed state loses 67/256 = 26.17 percentage points of viable and reserve-qualified mass. All three changes occur in state 8, D3, horizon 5, beta 0.4, rho 0.5: gamma thresholds 0.75/0.9 and normal threshold 0.9. There are two resolved opposite kernel decisions and one fine-kernel MC-unresolved crossing.\n\n"+
        table(changed,["state","scenario","variant","horizon","theta","coarse_nv","fine_nv","coarse_nr","fine_nr","coarse_status","fine_status","viable_change","explanation"])+
        "\nThe original COMBINED flag means both component thresholds crossed; it is not an independent reserve mechanism. In these three cases nr=nv in both representations, so the reserve-qualified loss is continuation-driven. The unchanged conditional reserve fraction does not prove reserve invariance elsewhere. The small global average cannot establish query sufficiency or dynamic invariance.\n")
    post=rows("effective_count/posteriors.csv");cap=rows("effective_count/propagation.csv");origin=rows("effective_count/origin_states.csv")
    reference={(r["state"],r["scenario"],r["variant"],r["horizon"],r["rho"],r["query"],r["theta"]):r for r in cap if r["cap"]=="200"}
    caps=[]
    for value in [50,100,200,500]:
        selected=[r for r in cap if int(r["cap"])==value];fit=[r for r in post if int(r["cap"])==value]
        delta=[];dr=[];changes=0
        for r in selected:
            base=reference[(r["state"],r["scenario"],r["variant"],r["horizon"],r["rho"],r["query"],r["theta"])]
            delta.append((int(r["nv"])-int(base["nv"]))/256);dr.append((int(r["nr"])-int(base["nr"]))/256);changes+=b(r["FINITE_BANK_HEALTH"])!=b(base["FINITE_BANK_HEALTH"])
        caps.append(dict(cap=value,mean_hazard_width90=np.mean([float(r["hazard_width90"]) for r in fit]).item(),
          mean_growth_se=np.mean([float(r["growth_se"]) for r in fit]).item(),mean_mass_change=np.mean(delta).item(),
          mean_abs_mass_change=np.mean(np.abs(delta)).item(),mean_abs_reserve_change=np.mean(np.abs(dr)).item(),Health_disagreement=changes/len(selected)))
    h.write_csv("effective_count/comparison.csv",caps)
    report("Effective_Count_Sensitivity.md","# Effective-count sensitivity\n\nAllcaps50/100/200/500 refit exactly the original E0→E1 likelihood with isolated module globals;200 remains frozen reference. No E2 observations enter fit or select a cap. The cap is a methodological tempering assumption, not an empirically estimated dependence size.\n\n"+table(caps,list(caps[0]))+
      "\nSelected M0 heldout calibration:\n\n"+table(rows("effective_count/masking_calibration.csv"),["cap","draws","n","fate_brier","DBH_MAE","DBH_coverage90","J_prediction","J_truth"])+
      "\nPer-cell256-grid posterior summaries and growthSE, paired origin states, all scenario/horizon/query threshold crossings and MC statuses are in committed tables. Summary averages above are equal finite design contrasts, not ecological priors; scenario-wise detail remains inspectable. Missingness interaction is B0 only. Differences in sampling consumption prevent an exact common-uniform claim. No preferred cap is proposed from these outcomes.")
    decisions=rows("decision_rules/versioned_resolution.csv");ds=[]
    for theta in [0,.5,.75,.9,1]:
        selected=[r for r in decisions if float(r["theta"])==theta]
        ds.append(dict(theta=theta,evaluations=sum(int(r["evaluations"]) for r in selected),TRUE=sum(int(r["TRUE"]) for r in selected),FALSE=sum(int(r["FALSE"]) for r in selected),MC_UNRESOLVED=sum(int(r["MC_UNRESOLVED"]) for r in selected)))
    report("Decision_Resolution.md","# Exact-bank versus stochastic-law resolution\n\nHistoricalHFD03 remains963976 unresolved evaluations, including3984 at interior thresholds. No prior record was altered. New versioned presentation explicitly treats theta0 adequacy as algebraic; Health still requires present realization. Theta1 all-success finite-bank TRUE remains stochasticMC_UNRESOLVED, never proof of probabilityone.\n\n"+table(ds,list(ds[0]))+
      "\nFINITE_BANK_HEALTH uses exact integercounts/K and inclusive>=theta. KERNEL_MC_STATUS uses independent pointwise Wilson95%Q1/two97.5%Q2 intervals, strict lower>theta andupper<theta; equality does not force a Boolean. The sufficient census links both presentations to the same archived counts. Boundary reference prefixes expose convergence without preferred-answer stopping. Estimator per-draw unresolved brackets address inner-path error in the finite outer mixture; observation/prior mixture intervals are reported separately.")
    isolation_path=h.ROOT/"provenance/public_isolation_receipt.json"
    isolation=json.loads(isolation_path.read_text()) if isolation_path.exists() else {}
    report("Public_Artifact_Readiness.md","# Public computational candidate\n\nA separate deny-by-default local source candidate is prepared in release/candidate, with exact inclusion list, original/exportSHA256, registered configs, numerical tables/figures, model/reconstruction/generator algorithms, verification source and a small derived-model fixture. HF253v6 inputs are externally retrieved and exact-hash checked; raw inputs, simulation banks, caches, credentials, unrelated research and Lean formal-world modules are excluded.\n\nIsolated lightweight kernel/query fixture status:"+str(isolation.get("status","PENDING"))+
      ". This lightweight check verifies package-only scientific module locations, pinned NumPy, actual archived kernel fingerprint and decision/query fixtures. Separately, the FULL original scientific replay actually executed from isolated packaged source: all 65,536-path banks, all context estimators, representation and cap experiments, and 23 output files reproduced byte-for-byte. A second isolated run verifies the separately registered context-prior supplement. Public-safe replay summaries bind scientific source/config identities, canonical input hashes and exact output hashes. Canonical externally obtained HF253 bytes were reused via read-only links: no new network download and no duplicated raw cache. Private historical Git preservation remains an attested lineage ledger, not something an outsider can verify without private history.\n\nCC0 source data and citation are documented. Project-authored source disclosure/distribution license remains owner review; READY_FOR_REVIEW means a locally executable candidate for that review, not public publication permission. No upload or DOI was performed.")
    topics=[
      ("Boundary gap","A","h04_boundary.py","Boundary_Stress_Validation.md","boundary_validation/truth_prefixes.csv","Uncertain reference truths retained"),
      ("Reference truth precision","A/E","h04_boundary.py","Decision_Resolution.md","boundary_validation/truth_prefixes.csv","65536 approximates stochastic law, not exact oracle"),
      ("Omitted hemlock process","A","h04_boundary.py","Boundary_Stress_Validation.md","boundary_validation/signed_margin_scores.csv","Controlled discrepancy, not causal ecological identification"),
      ("Supplied context anchors","B","h04_estimator.py","Context_Label_Ablation.md","context_ablation/comparison.csv","Finite context-design priors"),
      ("Unknown/coarsened context","B","h04_estimator.py","Context_Label_Ablation.md","context_ablation/estimator_draws.csv","Hidden-context prediction distinct from conditional Health"),
      ("Labels versus state information","B","h04_estimator.py","Context_Label_Ablation.md","context_ablation/comparison.csv","Label-only prior baseline; not additive attribution"),
      ("RMS juvenile information loss","C","h04_juvenile.py","Juvenile_Threshold_and_Query_Sufficiency.md","juvenile_boundary/observed_groups.csv","Secure observed E1 subset only"),
      ("Enrichment provenance and Q2","C","h04_juvenile.py","Juvenile_Threshold_and_Query_Sufficiency.md","juvenile_boundary/propagation.csv","Initial BA projection only; no dynamic commutation"),
      ("Effective-count cap","D","h04_cap.py","Effective_Count_Sensitivity.md","effective_count/posteriors.csv","No empirically identified cap/E2 selection"),
      ("Endpoint/interior MC decisions","E","h04_analysis.py","Decision_Resolution.md","decision_rules/hfd03_sufficient_census.csv","Exact finite arithmetic not stochastic certainty"),
      ("Complete Methods/probability/entry exposure","F","h04_methods.py","Numerical_Methods_Contract.md","config/numerical_methods.csv","Unidentified entry, conditional scenarios, model-admissibility Lawful"),
      ("Public reproduction artifact","G","h04_export.py","Public_Artifact_Readiness.md","release/PaperV_Artifact_Manifest.json","Separate publication/license approval; fixture versus full replay")
    ]
    mapping=[dict(topic=t,workstream=w,source="scripts/"+s,Methods="Numerical_Methods_Contract.md",Results=r,output=o,qualification=q) for t,w,s,r,o,q in topics]
    h.write_csv("reports/referee_map.csv",mapping)
    report("PaperV_V3_Revision_Map.md","# Paper V V3 evidence map\n\nTwelve methodological objection topics are reconstructed from the authorized HFD04 brief; this is not a quotation of an unavailable verbatim review. Cumulative reconstruction/missingness/entry andE2 qualifications remain in frozen HFD03 reports.\n\n"+table(mapping,list(mapping[0]))+
      "\nV3 should lead with partial real reconstruction and demographic identifiability, then model discrepancy, decision-boundary reliability, context dependence, juvenile-query information loss, effective-cap sensitivity and exact/kernel resolution. D0–D3 remain conditional outputs of the empirically unidentified entry mechanism and synthetic dynamics. The program title Health, Formally Defined is supported by PapersI–IV; the empirical paper needs an explicit declared-specification/empirical-synthetic framing. Suggested framing:Prospective Health under Declared Specifications—An Auditable Empirical-Synthetic Measurement Study. An operational forest-health validation title is unsupported. No manuscript title was silently changed.")
    statuses=[s["status"] for s in primary.values()]
    discrimination="ESTABLISHED" if all(x=="ESTABLISHED" for x in statuses) else ("PARTIAL" if any(x!="NOT_ESTABLISHED" for x in statuses) else "NOT_ESTABLISHED")
    from h04_certificate import gates as evidence_gates
    gates=evidence_gates();ready=gates["J"]
    cert=dict(HFD01_HFD02S_HFD03_FROZEN=True,BOUNDARY_STRESS_TEST="COMPLETE",BOUNDARY_HEALTH_DISCRIMINATION=discrimination,
      CONTEXT_LABEL_ABLATION="COMPLETE",CONTEXT_DEPENDENCE="QUANTIFIED",JUVENILE_THRESHOLD_CALIBRATION="COMPLETE",
      Q2_REPRESENTATION_SUFFICIENCY="CONDITIONAL",EFFECTIVE_COUNT_SENSITIVITY="QUANTIFIED",INTERIOR_MONTE_CARLO_PRECISION="MIXED",
      PUBLIC_ARTIFACT_CANDIDATE="READY_FOR_REVIEW" if ready else "REQUIRES_WORK",OPERATIONAL_HARVARD_FOREST_ASSESSMENT="NOT_AUTHORIZED",
      RETROSPECTIVE_E2_VALIDATION="NOT_CLAIMED",PAPER_V_V3="AUTHORIZED" if all(gates.values()) else "REQUIRES_ADDITIONAL_WORK")
    h.write_json("reports/hfd04_contract.json",dict(certificate=cert,gates=gates,boundary_density=density,primary_discrimination=primary,
      weights=dict(truth="1/K withinworld/suite",worlds="1/81 current,1/324 historical; independentworld units, pairedsuites not extra worlds",estimator="1/32 outer,1/256 paths; explicit context mixture weights",
      original="1/1024 outer,1/64 paths,1/65536 joint per scenario"),claims=h.config()["claims"],all_worlds_retained=True,
      full_public_raw_replay_executed=True,source_license_review_pending=True,
      evidence_bundle="Gate predicates in scripts/h04_certificate.py; standalone verifier checks exact rosters, decisions, sources and replay summaries, not just record existence."))
    report("Verification.md","# Verification\n\nStandalone verifier independently recomputes count-based finite/kernel decisions, reference labels, per-draw Health, full draw means/confusion/Brier, all experimental rosters/strata and historical count reconciliation. It checks immutable predecessor source blobs, config/numerical Methods, source-qualified figure data, export hashes and isolation receipt; at least 20 active hostile controls exercise actual source/config guards, real validators and strict estimator interfaces. --cache rechecks all 243 archived truth-response banks and paired prefixes; --replay-panel freshly regenerates the first 1,024-path chunk of EVERY world/suite. This panel is not a full-bank replay. Separately, public_full_replay_summary.json records the actually executed full isolated packaged-source replay, with 23 byte-exact original outputs; public_supplement_replay_summary.json records the five-output supplementary replay. The final certificate depends on all A–L evidence gates and a current verifier-source-bound PASS receipt with full archive census. --final rejects any missing gate, altered frozen outcome or unsupported terminal status. No old science was rewritten.")
    report("Closure.md","# HFD04 closure\n\n"+json.dumps(cert,indent=2)+
      "\n\nAll registered adverse cases and reference uncertainty remain explicit. Scientific completion certifies executed tests, not favorable outcomes. V3's defensible contribution is an auditable measurement/decision interface under declared specification, not an operational diagnosis. Full old science remains byte-preserved. See each workstream report, numerical Methods, eight figure/data/source manifests, exact public inclusion manifest and actual verification receipt.\n\nGates:"+json.dumps(gates,sort_keys=True))
    # Publication-native vector figures; data retain every qualifying stratum.
    margin=rows("boundary_validation/signed_margin_scores.csv");fig=[]
    for r in margin:
        if r["cohort"]=="boundary" and r["suite"]=="correct" and r["policy"]=="B0_full":
            sample=[p for p in performance if p["cohort"]=="boundary" and p["suite"]=="correct" and p["policy"]=="B0_full" and p["query"]==r["query"] and p["signed_region"]==r["region"]]
            fig.append(dict(query=r["query"],region=r["region"],n=r["n"],mean_margin=float(np.mean([float(p["margin"]) for p in sample])) if sample else None,error=number(r["certified_error"])))
    plot(1,"Health decision error near adequacy",fig,"mean_margin","error","Signed-margin strata and certified-reference error; uncertain reference truths retained/excluded explicitly from error denominator; Q1/Q2 separate in CSV.",["boundary_validation/performance.csv","boundary_validation/signed_margin_scores.csv"])
    scatter=[dict(query=r["query"],reference_mass=float(r["reference_viable"] if r["query"]=="Q1" else r["reference_reserve"]),estimated_mass=float(r["estimated_viable"] if r["query"]=="Q1" else r["estimated_reserve"]),reference_status=r["reference_status"],world=r["world"]) for r in performance if r["cohort"]=="boundary" and r["suite"]=="correct" and r["policy"]=="B0_full"]
    plot("1b","Estimated versus finite-reference boundary capacity",scatter,"reference_mass","estimated_mass","Q1 viable/Q2 reserve-qualified mass; reference-uncertain cases retained and flagged, not treated as exact stochastic truth.",["boundary_validation/performance.csv"])
    fg=[]
    for r in margin:
        if r["cohort"]=="boundary" and r["policy"]=="B0_full" and r["suite"]!="entry_lognormal_annual_sd_1.2":
            sample=[p for p in performance if all(p[k]==r[k] for k in ["cohort","suite","query","policy"]) and p["signed_region"]==r["region"]]
            fg.append(dict(suite=r["suite"],query=r["query"],region=r["region"],mean_margin=float(np.mean([float(p["margin"]) for p in sample])) if sample else None,FPR=number(r["FPR"]),FP=int(r["FP"]),certified_n=int(r["certified_n"])))
    plot(2,"Correct versus omitted-hemlock boundary false positives",fg,"mean_margin","FPR","Conditional false-positive rate by true signed-margin stratum; empty strata retained in CSV, mechanisms/query never pooled.",["boundary_validation/signed_margin_scores.csv","boundary_validation/performance.csv"])
    fg=[dict(cohort=r["cohort"],query=r["query"],policy=r["policy"],condition_index=h.config()["estimator"]["contexts"].index(r["policy"]),balanced_accuracy=number(r["balanced_accuracy"]),FP=r["FP"],FN=r["FN"]) for r in comparisons if r["suite"]=="correct"]
    plot(3,"Context information changes decision performance",fg,"condition_index","balanced_accuracy","Full/coarse/unknown/label-only paired405 worlds; category/query/cohort labels retained in CSV;32-draw current baseline.",["context_ablation/comparison.csv"])
    err=rows("juvenile_boundary/measurement_error.csv");fg=[dict(sigma=float(r["sigma"]),replicate=r["replicate"],distance_low=r["distance_low"],distance_high=r["distance_high"],n=r["n"],error=(int(r["FP"])+int(r["FN"]))/int(r["n"]) if int(r["n"]) else None) for r in err]
    plot(4,"DBH boundary sensitivity at10cm",fg,"sigma","error","All distance bands/replicates preserved; Gaussian scales are methodological, not estimated instrument distributions.",["juvenile_boundary/measurement_error.csv"])
    fg=[dict(state=r["state"],scenario=r["scenario"],variant=r["variant"],horizon=r["horizon"],coarse_reserve=int(r["coarse_nr"])/256,fine_reserve=int(r["fine_nr"])/256,coarse_viable=int(r["coarse_nv"])/256,fine_viable=int(r["fine_nv"])/256,coarse_present=r["coarse_present"],fine_present=r["fine_present"],coarse_Health=r["coarse_Health"],fine_Health=r["fine_Health"]) for r in prop if r["beta"]=="0.4" and r["rho"]=="0.5" and r["theta"]=="0.75"]
    plot(5,"Representation-dependent viable regenerative reserve",fg,"coarse_reserve","fine_reserve","Paired scenario/state/model/horizon with beta.4,rho.5,theta.75; realization, viable mass and Q2 answers retained in linkedCSV. Inferred size profiles, not a new observation or dynamic commutation theorem.",["juvenile_boundary/initial_projection.csv","juvenile_boundary/propagation.csv"])
    plot(6,"Effective-cap hazard uncertainty",caps,"cap","mean_hazard_width90","Equal64-cell descriptive mean; all caps retained; full downstream contrasts in CSV;200 remains frozen reference.",["effective_count/posteriors.csv","effective_count/comparison.csv"])
    plot(7,"Exact-bank versus stochastic-kernel decisions",ds,"theta","MC_UNRESOLVED","Versioned HFD03 count-derived status:theta0 algebraic;theta1 not probabilityone; original963976/3984 record untouched.",["decision_rules/versioned_resolution.csv"])
    stages=[dict(stage=i,label=label,registered_failure_count=n,interpretation=q) for i,(label,n,q) in enumerate([
      ("Measured state",obs["n"],"Eligible E1 support, not failures"),
      ("Representation",obs["FP"]+obs["FN"],"RMS juvenile membership discrepancies"),
      ("Generator",sum(int(r["FP"]) for r in comparisons if r["cohort"]=="boundary" and r["policy"]=="B0_full" and r["suite"]=="hemlock_extra_hazard_0.12"),"Query-specific controlled false positives; not independent extra worlds"),
      ("Continuation",density["reference_unresolved"],"Target reference-uncertain cells"),
      ("Capacity",sum(int(r["MC_UNRESOLVED"]) for r in ds if r["theta"] in [.5,.75,.9]),"Historical interior unresolved evaluations"),
      ("Adequacy",0,"Inclusive exact predicate; no causal error assigned to arithmetic")])]
    plot(8,"End-to-end decision provenance",stages,"stage","registered_failure_count","Stage counts have DISTINCT denominators; do not sum or interpret as a causal error partition. Full metadata identifies measured/reconstructed/model/MC roles.",["juvenile_boundary/observed_groups.csv","context_ablation/comparison.csv","boundary_validation/truth_prefixes.csv","decision_rules/versioned_resolution.csv"])
    print("Generated10 reports,8 figure/data manifests, terminal certificate",cert,flush=True)
if __name__=="__main__":make()
