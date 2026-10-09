"""Data-derived reports, manuscript specification, provenance and nine SVG plots."""
import csv
import html
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
from h03_common import ROOT, REPO, RUN, PROTOCOL, PROTOCOL_HASH, frozen, reconstruct, read_csv, write_csv, write_json


def load(name):
    return json.loads((ROOT/name).read_text())


def markdown(name,text):
    path=ROOT/name
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(text.rstrip()+"\n")


def table(rows,fields):
    body=["| "+" | ".join(fields)+" |","| "+" | ".join("---" for _ in fields)+" |"]
    for row in rows:
        body.append("| "+" | ".join(str(row[k]).replace("|","/") for k in fields)+" |")
    return "\n".join(body)


def bars(name,title,labels,values,unit,second=None):
    width=920;height=120+len(labels)*34
    maximum=max(max(values,default=0),max(second,default=0) if second is not None else 0,1e-12)
    items=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="#faf8f3"/>',
           f'<text x="24" y="30" font-family="sans-serif" font-size="18">{html.escape(title)}</text>',
           f'<text x="24" y="54" font-family="sans-serif" font-size="12">{html.escape(unit)}</text>']
    for i,(label,value) in enumerate(zip(labels,values)):
        y=76+i*34
        items.append(f'<text x="24" y="{y+14}" font-family="sans-serif" font-size="11">{html.escape(str(label))}</text>')
        if second is None:
            bar_width=560*value/maximum
            items.append(f'<rect x="300" y="{y}" width="{max(bar_width,0):.4f}" height="18" fill="#587c59"/>')
        else:
            bar_width=560*value/maximum
            other_width=560*second[i]/maximum
            items.append(f'<rect x="300" y="{y}" width="{max(bar_width,0):.4f}" height="8" fill="#587c59"/>')
            items.append(f'<rect x="300" y="{y+10}" width="{max(other_width,0):.4f}" height="8" fill="#b56848"/>')
        items.append(f'<text x="866" y="{y+14}" font-family="sans-serif" font-size="10">{value:.3f}</text>')
    items.append("</svg>")
    result="\n".join(items)
    if name:
        markdown("figures/"+name+".svg",result)
    return result


def scatter(name,title,rows,x,y,anchors=None,diagonal=False):
    width=920;height=620
    xmax=max(max(float(r[x]) for r in rows),max([a[0] for a in anchors],default=0) if anchors else 0,1)
    ymax=max(max(float(r[y]) for r in rows),max([a[1] for a in anchors],default=0) if anchors else 0,1)
    if diagonal:
        xmax=ymax=max(xmax,ymax)
    items=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
           '<rect width="100%" height="100%" fill="#faf8f3"/>',
           f'<text x="24" y="30" font-family="sans-serif" font-size="17">{html.escape(title)}</text>',
           f'<text x="24" y="55" font-family="sans-serif" font-size="12">x: {html.escape(x)} [0,{xmax:.2f}] — y: {html.escape(y)} [0,{ymax:.2f}]</text>',
           '<path d="M80 90 V540 H850" fill="none" stroke="#555"/>']
    for row in rows:
        px=80+float(row[x])/xmax*740
        py=540-float(row[y])/ymax*430
        color={"hemlock":"#587c59","other_canopy":"#326b9c","shrub":"#b56848","other":"#896ba4"}.get(
            row.get("taxon_group"),"#587c59" if str(row.get("truth_health","False"))=="True" else "#b56848")
        tip=" / ".join(str(row.get(k,"")) for k in ("taxon_group","E0_size_class","stratum","mask","replicate","world"))
        items.append(f'<circle cx="{px:.3f}" cy="{py:.3f}" r="3" fill="{color}" opacity="0.65"><title>{html.escape(tip)}</title></circle>')
    if diagonal:
        items.append('<path d="M80 540 L820 110" stroke="#555" stroke-dasharray="4,4" fill="none"/>')
    for a,b,label in anchors or []:
        px=80+a/xmax*740;py=540-b/ymax*430
        items.append(f'<path d="M{px:.3f} 90 V540 M80 {py:.3f} H850" stroke="#326b9c" stroke-dasharray="5,5"/>')
        items.append(f'<text x="90" y="580" font-family="sans-serif" font-size="12">{html.escape(label)}</text>')
    items.append("</svg>")
    markdown("figures/"+name+".svg","\n".join(items))


def run(core,model):
    a=load("reconstruction_validation/summary.json")
    b=load("missingness/summary.json")
    c=load("recruitment_audit/summary.json")
    d=load("e2_failure_localization/summary.json")
    d["frozen_original_comparison_csv"]=dict(
        path="HFD02S/validation/e2_model_comparison.csv",
        sha256=frozen.sha(ROOT.parent/"HFD02S/validation/e2_model_comparison.csv"),
        rows=read_csv(ROOT.parent/"HFD02S/validation/e2_model_comparison.csv"))
    write_json("e2_failure_localization/summary.json",d)
    e=load("synthetic_validation/summary.json")
    f=load("monte_carlo/summary.json")
    calibration=read_csv(ROOT/"reconstruction_validation/calibration.csv")
    all_cal=[r for r in calibration if r["subgroup"]=="ALL"]
    phase=read_csv(ROOT/"missingness/phase_surface.csv")
    inner=read_csv(ROOT/"monte_carlo/inner_resolution.csv")
    outer=read_csv(ROOT/"monte_carlo/outer_resolution.csv")
    performance=read_csv(ROOT/"synthetic_validation/performance.csv")
    diagnostics=read_csv(ROOT/"e2_failure_localization/response_surfaces.csv")
    # Produce complete lineage for all 1024 MC states even when an earlier
    # running numerical worker predated the supplementary provenance writer.
    provenance=[]
    for i in range(PROTOCOL["monte_carlo"]["outer_N"][-1]):
        state=reconstruct(core,model,i)
        total_ba=float((state["n"]*state["d"]**2).sum()*math.pi/40000)
        provenance.append(dict(family="MC_B0",state=i,count=int(state["n"].sum()),
            juveniles=int((state["n"]*(state["d"]<10)).sum()),basal_area=total_ba,
            observed_anchored_projected_ba=state["components"]["measured_projected_ba"],
            imputed_ba=state["components"]["imputed_ba"],
            observed_anchor_fraction=state["components"]["measured_projected_ba"]/max(total_ba,1e-9),
            directly_observed_origin_ba=0,origin_provenance="IMPUTED",
            **{k:state["provenance"][k] for k in ["observed_anchored_projected_count","imputed_count",
                "observed_anchored_projected_juveniles","imputed_juveniles","recovered_current_origin_ba"]}))
    write_csv("provenance/mc_state_contributions.csv",provenance)
    support=read_csv(ROOT/"provenance/state_contributions.csv")
    for row in support:
        state=reconstruct(core,model,int(row["state"]),row["family"])
        row.update({k:state["provenance"][k] for k in ["observed_anchored_projected_count","imputed_count",
            "observed_anchored_projected_juveniles","imputed_juveniles","recovered_current_origin_ba"]})
        row["observed_anchor_count_fraction"]=state["provenance"]["observed_anchored_projected_count"]/max(int(state["n"].sum()),1)
        row["observed_anchor_juvenile_fraction"]=state["provenance"]["observed_anchored_projected_juveniles"]/max(
            int((state["n"]*(state["d"]<10)).sum()),1)
    write_csv("provenance/state_contributions.csv",support)
    old=frozen.load("config/design.json")
    write_json("specification/frozen_hfd02s_numeric_design.json",old)
    spec=[]
    for scope,data in [("HFD02S",old),("HFD03",PROTOCOL)]:
        def flatten(prefix,value):
            if isinstance(value,dict):
                for k,v in value.items():
                    flatten(prefix+"."+k,v)
            else:
                spec.append(dict(scope=scope,parameter=prefix,value=json.dumps(value,ensure_ascii=False)))
        flatten("design",data)
    write_csv("specification/numerical_assumptions.csv",spec)
    write_csv("specification/scenarios.csv",old["scenarios"],fields=list(old["scenarios"][0]))
    comparison=read_csv(ROOT.parent/"HFD02S/uncertainty/decomposition.csv")
    write_csv("uncertainty/frozen_finite_design_contrasts.csv",comparison)
    finite={r["source"]:float(r["conditional_variance"]) for r in comparison}
    contrast=(min(finite["scenario"],finite["semantic"])>max(finite["parameter"],finite["model"]))
    weights=[
        dict(study="HFD02S frozen",states=256,parameters=2,models=2,K=64,outer_units=1024,history_weight="1/(1024*64)"),
        dict(study="HFD03 missingness, per family/scenario",states=128,parameters=1,models=2,K=256,outer_units=256,history_weight="1/(256*256)"),
        dict(study="HFD03 MC inner panel, per scenario",states=16,parameters=2,models=2,K=1024,outer_units=64,history_weight="1/(64*1024)"),
        dict(study="HFD03 MC outer maximum, per scenario",states=1024,parameters=2,models=2,K=64,outer_units=4096,history_weight="1/(4096*64)")]
    write_csv("specification/finite_mixture_weights.csv",weights)
    raw_fields=read_csv(ROOT.parent/"HFD02S/reconstruction/field_missingness.csv") if (ROOT.parent/"HFD02S/reconstruction/field_missingness.csv").exists() else core["fields"]
    write_csv("provenance/raw_field_support.csv",raw_fields)
    write_json("provenance/namespace_rules.json",dict(core_classes=sorted(frozen.PROVENANCE),
        display_alias={"DIRECTLY_OBSERVED":"OBSERVED"},raw_current="Only E1 recorded fields OBSERVED",
        prior_location="Stable stem/tree coordinates RECOVERED_FROM_PRIOR_OBSERVATION",
        origin2020="Every date-projected/fate-imputed/size-imputed origin field IMPUTED",
        anchored_contribution="Observed-source projected BA is ancestry support, not directly observed 2020 state",
        recovered_ba="Zero current BA is licensed merely by recovered root coordinates; location recovery is a separate field-support measure",
        count_support="Measured-source projected counts versus imputed counts; juvenile source allocation follows final cohort RMS classification",
        future="SIMULATED",thresholds="EXTERNALLY_CONSTRAINED",
        headline_links=["provenance/state_contributions.csv","provenance/mc_state_contributions.csv"]))
    write_json("provenance/source_catalog.json",dict(
        primary_observed=["HF253-E0","HF253-E1"],
        archived_metadata="HFD01/manifests/artifact_manifest.json",
        new_identity_sources={"HF253-trees2014":"identity_inputs.json:sources[0]",
                              "HF253-species-codes":"identity_inputs.json:sources[1]"},
        scientific_prior="HFD02S config/design.json and immutable empirical fit",
        new_method_sources={"HFD03-mask-protocol":"config/protocol.json:calibration",
                            "HFD03-entry-assumptions":"config/protocol.json:entry_audit",
                            "HFD03-B0-B4":"config/protocol.json:missingness",
                            "HF-SYNTH-DISC-001":"config/protocol.json:synthetic",
                            "HFD03-paired-maximum-banks":"config/protocol.json:monte_carlo",
                            "HFD03-diagnostic-grid":"config/protocol.json:diagnostics"},
        exposed_external_source="HF-third-census-news from immutable HFD01 manifest; never training"))
    headline=[]
    for row in outer:
        if int(row["N_X"])==1024:
            headline.append(dict(row,provenance_link="provenance/mc_state_contributions.csv",
                                 state_prefix="0..1023",interpretation="finite-design Health fraction bracket"))
    write_csv("provenance/headline_health_links.csv",headline)
    editorial=[
        ("Health in empirical section","Health under declared demonstration specification / H_S"),
        ("Lawful in empirical section","model-admissibility realization of Lawful; bookkeeping is not ecological validation"),
        ("wet-area uncertainty","wet-candidate/protocol-proxy uncertainty; no certified hydrological geography"),
        ("ingrowth/recruitment for new IDs","candidate threshold entrants; membership bounds are assumptions"),
        ("scenario uncertainty dominates","Across declared finite ranges, scenario/semantic contrasts numerically larger than fitted parameter/model contrasts"),
        ("M3 new empirical information","Supporting worked query-relative sufficiency illustration; redundant M3 adds no information"),
        ("ensemble probability","finite design-mixture fraction; not ecological frequency or posterior over scenarios/semantics")]
    write_csv("specification/editorial_mapping.csv",[dict(old=x,revised=y) for x,y in editorial])
    calibration_table=[{k:r[k] for k in ["mask","replicate","n","brier","dbh_mae","coverage_50","coverage_80","coverage_90"]} for r in all_cal]
    poorest=sorted([r for r in calibration if r["subgroup"]!="ALL" and int(r["dbh_n"])>=20],
                   key=lambda r:float(r["coverage_90"]))[:12]
    masks=Counter(r["mask"] for r in all_cal)
    all_unresolved=sum(int(r["unresolved_units"]) for r in phase)
    full_total=sum(int(r["true_units"])+int(r["false_units"])+int(r["unresolved_units"]) for r in phase)
    unresolved_by_theta=[]
    for theta in [0,.5,.75,.9,1]:
        selected=[row for row in phase if float(row["theta"])==theta]
        unresolved_by_theta.append(dict(theta=theta,conditional_units=len(selected)*256,
                                       MC_UNRESOLVED=sum(int(row["unresolved_units"]) for row in selected)))
    write_csv("monte_carlo/full_grid_threshold_resolution.csv",unresolved_by_theta)
    certificate=dict(HFD02S_FROZEN=True,REAL_DATA_RECONSTRUCTION_CALIBRATION="COMPLETE",
        MISSINGNESS_SENSITIVITY="COMPLETE",THRESHOLD_ENTRY_IDENTIFICATION="BOUNDED",
        E2_GENERATOR_FAILURE_LOCALIZED="PARTIAL",SYNTHETIC_HEALTH_DISCRIMINATION=e["discrimination"],
        CONTROLLED_MISSPECIFICATION_TEST="COMPLETE",
        PUBLICATION_MONTE_CARLO_PRECISION="MIXED" if all_unresolved else f["precision"],
        OPERATIONAL_HARVARD_FOREST_ASSESSMENT="NOT_AUTHORIZED",PAPER_V_V2="AUTHORIZED")
    contract=dict(schema="HFD03-contract-v1",protocol_sha256=PROTOCOL_HASH,
        planning_scientific_baseline=PROTOCOL["planning_scientific_baseline"],
        actual_execution_start=PROTOCOL["actual_execution_start"],
        final_hfd03_commit=None,
        commit_identity_semantics="Closure fills final_hfd03_commit with the final scientific-content commit. Its publication-metadata successor is identified from git log on provenance/publication_receipt.json, avoiding a self-referential hash.",
        terminal_certificate=certificate,selected_mc_precision=f["precision"],
        full_grid_unresolved_conditional_cells=all_unresolved,full_grid_conditional_cells=full_total,
        gates=dict(A=len(all_cal)==12,B=len(b["results"])==5,C=c["measured_living_candidates"]==6992,
                   D=d["diagnostic_designs"]>0 and not d["parameter_selection"],
                   E=e["results"]["correct"]["positives"]>0 and e["results"]["correct"]["negatives"]>0,
                   F=len(e["results"])==3,G=len(inner)==10368 and len(outer)==216,
                   H=len(spec)>100,I=contrast,J="PENDING_VERIFIER"),
        claims=PROTOCOL["claims"],scientific_limits=[
            "Coarse RMS DBH intervals under-cover actual individual heldout measurements.",
            "Secure death pool does not turn all G statuses or whole-plant disappearance into death.",
            "Entry membership weights are sensitivity assumptions, not identified probabilities.",
            "Public E2 definitions/exposures are incomparable: no exact residual or validated correction.",
            "Correct-model discrimination is controlled self-consistency; hemlock discrepancy causes false positives.",
            "No primary synthetic truth masses within 0.1 of theta: mass-boundary performance not established.",
            "Selected MC panel meets targets; broader missingness/semantic surface retains unresolved cells.",
            "Mixture fractions are conditional finite-design summaries, not operational ecological Health."])
    publication=ROOT/"provenance/publication_receipt.json"
    if publication.exists():
        contract["final_hfd03_commit"]=json.loads(publication.read_text())["final_hfd03_commit"]
    write_json("reports/hfd03_contract.json",contract)
    if (ROOT/"provenance/interval_boundary_amendment.json").exists():
        contract["full_model_replay_receipt"]="reports/full_model_replay_receipt.json"
        contract["interval_recount_receipt"]="provenance/interval_boundary_amendment.json"
        write_json("reports/hfd03_contract.json",contract)
    design=f"""# HFD03 Design

Scientific baseline `{PROTOCOL['planning_scientific_baseline']}`; actual execution start `{PROTOCOL['actual_execution_start']}`.
Registration hash `{PROTOCOL_HASH}`. All original HFD01/HFD02S tracked bytes and verifiers remain frozen.

Five workstreams test reconstruction, missingness, entry identifiability, generator discrepancy and discriminative recovery; a separate numerical study tests MC resolution. Config was written before outcome comparisons. Targeted public E0 tree/species inputs have publisher MD5, length and SHA256 receipts. No published E1 tree table exists in HF253 v6. Derived E1 grouping is not an observed plant census.

Read [protocol](../config/protocol.json), [registration](../provenance/registration.json), [source receipts](../provenance/identity_inputs.json). Three masks per family, 256 reconstruction draws; B0–B4 use 128 states × 1 parameter × 2 growth models × 4 scenarios × 256 histories through year20. Full old semantic grid, horizons5/10/20. Seed coupling is not asserted to be identical common uniforms where generator consumption varies.

No new Lean theorem or Swift change is needed for this empirical tranche: the existing formal interface is instantiated by conditional realization, model-admissible continuation and unconditioned capacity. Scientific validation is separate from formal implication.
"""
    aggregate_coverage=[]
    for mask in ["M0","M1","M2","M3"]:
        subgroup=[r for r in calibration if r["mask"]==mask and r["subgroup"]!="ALL"]
        for nominal in [50,80,90]:
            aggregate_coverage.append(dict(mask=mask,nominal=nominal,cohort_intervals=len(subgroup),
                living_coverage=sum(r["living_covered_"+str(nominal)]=="True" for r in subgroup)/len(subgroup),
                ba_coverage=sum(r["ba_covered_"+str(nominal)]=="True" for r in subgroup)/len(subgroup),
                juvenile_coverage=sum(r["juvenile_covered_"+str(nominal)]=="True" for r in subgroup)/len(subgroup)))
    write_csv("reconstruction_validation/aggregate_interval_coverage.csv",aggregate_coverage)
    reconstruction=f"""# Real-data reconstruction calibration

Eligible pool: {a['eligible']} stable tagged transitions. Alive/A and dated stem-dead records corroborated by explicit D codes qualify; Miss/MT/S and ambiguous associations, absent/gone/broken/missing fates do not. Secure death is stem-level, never whole-plant death.

For each mask, all nonidentity measurement fields (including otherwise-unused AGB/measurement-ID fields) are physically removed from a new visible CSV before unchanged HFD02S read_core, fit and reconstruct run. Only stem/tree/tag/taxon identity metadata remains. All fitted priors are rebuilt from visible transitions, not a fit on full E1. Assignment uses only E0 taxon/size/stratum, known campaign month and observed missing-field pattern; the truth evaluator is downstream.

The first reconstruction remains the 2020-origin kernel. Its origin medians are recorded but NOT scored against 2018–2019 observations as though simultaneous. Calibration scores use a tagged-endpoint projection of the exact sampled unknown-fate kernel at actual exposure. This retains cohort RMS DBH, MNAR shifts, 1cm prior spread and original fit; it does not silently replace the estimator with individual-growth inference.

{table(calibration_table,list(calibration_table[0]))}

50/80/90% coverage has no invented pass gate. Individual DBH and annualized growth bias/MAE/RMSE and interval widths, living count/BA/juvenile errors AND predictive intervals, taxon-share errors, fate calibration/Brier/log scores and coordinate recovery are committed. Growth coverage is the same affine interval event as DBH coverage, with annualized widths. Missing alive DBH cannot provide observed BA truth: those few records are excluded at scoring ONLY, with an explicit census, rather than comparing a partial observed total against a full prediction.

Aggregate heldout-cohort interval coverage:

{table(aggregate_coverage,list(aggregate_coverage[0]))}

These are descriptive pooled cohort/replicate intervals with shared fitted draws; they are not independent binomial coverage trials. The frozen state kernel is cohort RMS, so poor individual-size coverage does not by itself establish poor aggregate BA calibration. Aggregate intervals test the latter separately.

Worst declared cohort slices (minimum20 observed diameters) by 90% coverage:

{table(poorest,['mask','replicate','taxon_group','E0_size_class','stratum','dbh_n','dbh_mae','coverage_90'])}

Cell code: (taxon×4+stratum)×4+E0 size, as frozen in HFD02S. Wider size bins show poor individual coverage: cohort RMS is a structural limitation, not evidence that unseen diameters are accurately reconstructed. Leave-observed-out is conditional on securely observed records and cannot validate nonignorable fate in the unobserved proxy region.

Data: [calibration](../reconstruction_validation/calibration.csv), [fate curve](../reconstruction_validation/fate_curve.csv), [coordinate recovery](../reconstruction_validation/coordinate_recovery.csv), [taxon composition](../reconstruction_validation/taxon_composition.csv), [eligible IDs](../reconstruction_validation/eligible_truth_pool.csv), [assignments](../reconstruction_validation/mask_assignments.csv).
"""
    missing=f"""# Missingness sensitivity

All five frozen families propagated state → future histories → unconditional capacity → conditional Health across {b['results']['B0']['phase_cells']} phase specifications per family. B0 uses the original reconstruction unchanged. B1 has deterministic observed-cell hazard shifts, B2 narrower latent shifts, B3 broader shifts, B4 stronger wet-proxy-specific shifts. All preserve original measured projections, entry borrowing and association branches. B1 cell shifts apply ordinary unknown-fate cohorts; the 31 unresolved association-branch priors remain frozen, an explicit residual limitation.

These are declared conditional latent-fate reconstruction families, not identified likelihood models of the data-acquisition mechanism or fitted probabilities that a stem goes unobserved. B1 depends on recorded taxon/size/stratum and exposure/protocol covariates; it does not learn unseen fate from the missingness pattern.

{table([dict(family=k,**v) for k,v in b['results'].items()],['family','phase_cells','point_disagreement','mean_mass_change'])}

Disagreement averages paired conditional unit/specification labels, not a probability over mechanisms. Continuous changes precede thresholding. Phase-cell means, p05/p95 mass distributions, scenarios and semantic cells are in [full surface](../missingness/phase_surface.csv); [state distributions](../missingness/state_distribution.csv) and [linked support](../provenance/state_contributions.csv) expose the reconstructed burden.

TRUE/FALSE/MC_UNRESOLVED classifications produce a lower/upper bracket for the finite-design Health fraction. Point labels are descriptive. Across all five families, {all_unresolved}/{full_total} conditional cells remain MC-unresolved (including strict-rule endpoints theta0/1). Do not equate disagreement with established biological reversal.
"""
    recruitment=f"""# Candidate threshold-entry audit

All6992 measured living E1-only stems audited, not declared biological ingrowth. E0 tree table has83813 rows; E1 tree table is unavailable. Tree/stem IDs, both tags, taxon, positions/quadrat, size/POM/date and protocol proxy evidence are recorded individually.

{table(read_csv(ROOT/'recruitment_audit/class_census.csv'),['category','n','reference_weight'])}

Annual candidate exposure bounds (six-year convention): {json.dumps(c['annual_candidate_rates'],sort_keys=True)}. No E0 prior-coded root licenses certain threshold entry, so lower bound0 is honest. Most candidates are small, compatible with entry but not proof. 1337 attach to previously alive roots,24 to previously dead roots,5631 have no E0 parent-tree record. A new stem on an existing root is not new whole-plant recruitment.

Reference class memberships are externally stipulated sensitivity weights, not calibrated probabilities. Upper bound assumes every measured candidate contributes exposure, not that all are biologically identified. These bounds propagate to the diagnostic generator INCLUDING inherited wet entry borrowing; that extrapolation is not identified wet-region entry. See [individual audit](../recruitment_audit/candidate_audit.csv).
"""
    failure=f"""# Frozen E2 generator failure localization

Original D0 median entries9073, hemlock deaths2421 and all-stem deaths15296 remain immutable. Public E2 reports approximately5000 newly recorded threshold stems, almost5000 hemlock deaths since2020, and OVER11800 dead stems. These differ in exposure dates, frame, identity and status interpretation. No exact aligned failure vector is licensed; particularly11800 is a lower bound, not an exact all-death target.

{d['diagnostic_designs']} predeclared designs perturb candidate entry interpretation, hemlock hazard, general hazard, wet-proxy counts, duration4/5/6years and two-point initial-size resolution; joint entry×hemlock×general-hazard grid also executed. [Response surface](../e2_failure_localization/response_surfaces.csv) overlays approximate public anchors without optimizing parameters or calling a close point a validated correction.

Entry causes: {', '.join(d['entries_causes'])}. Hemlock causes: {', '.join(d['hemlock_causes'])}. Attribution remains PARTIAL: factor sensitivity localizes vulnerabilities, not unique ecological causes. Extra hemlock hazard represents a process absent from D0. Count-size splitting exposes juvenile threshold sensitivity but cannot introduce missing size-dependent/species hazards; original mortality is size-independent. Two subcohort histories retain unconditional weights and count balance.

Diagnostic scale, not a chosen correction: the predeclared reference candidate membership gives median entries4533.5 versus9127.5 under the all-candidate rate in this separate eight-state diagnostic sample. Thus candidate interpretation can account for a large portion of the entry magnitude, but the membership weight is uncalibrated and comparison scope remains unresolved. Do not replace frozen9073 with4533.5. Hemlock excess annual hazard0.04 changes ancestral hemlock death median2250 to5889 in this panel; raising general hazard to1.5 alone gives3253. Four-to-six-year exposure gives1826–2660 without excess hazard. This points to a missing hemlock-specific process, not a uniquely identified fitted cause. These hemlock diagnostics count initially present hemlocks only; original frozen cumulative deaths can also include new cohorts. The two-point-size run has a maximum initial BA discrepancy0.2914m² from diameter clamping, explicitly archived.
"""
    correct=e["results"]["correct"]
    primary_worlds=[row for row in performance if row["suite"]=="correct"]
    present_count=sum(row["truth_present"]=="True" for row in primary_worlds)
    present_inadequate=sum(row["truth_present"]=="True" and row["truth_health"]=="False" for row in primary_worlds)
    synthetic=f"""# Balanced synthetic Health discrimination

New benchmark `HF-SYNTH-DISC-001`; old all-negative HFD02S test untouched. Full3×3×3×3×4 grid gives324 independent world states. Truth census persisted before estimator performance: {correct['positives']} positive and {correct['negatives']} negative primary reference cases. No selection by truth Health; misspecification suites paired on the same states, not additional independent state replicates.

Nonvacuity: {present_count} primary worlds presently realize the organization; {present_inadequate} of them fail prospective adequacy. The test therefore exercises both halves of Health, not only a current-state realization label.

Correct-model metrics:

```json
{json.dumps(correct,indent=2)}
```

The estimator receives only strict noisy/thinned observation packets and public imposed pressure/recovery context labels; no actual generating rates, latent state, structural/regenerative level or truth stream. Its state anchor is fixed public E0. Truth and estimator rates share a declared family in the correct suite, not hidden point parameters.

“Correct model” refers to the demographic response kernel, rate-family and observation operator. The fixed E0 state prior is not claimed to be the generating distribution of the Cartesian synthetic state grid.

Truth is finite4096-path reference with Wilson status, not exact infinite-law Health. Estimator uses64 observation-conditioned draws ×256 futures. Nominal90% mass intervals mix state/parameter uncertainty with finite predictive sampling and are not pure MC intervals.

Under the predeclared discrimination gate the status is {e['discrimination']}. This is controlled recovery under a known model class, not real-forest validation. No primary viable masses lie within0.1 of theta: mass-boundary behavior remains untested despite current-state boundary examples. All misses/false positives remain in [performance](../synthetic_validation/performance.csv).
"""
    misspec=f"""# Controlled misspecification

Do not pool correct-model recovery with robustness. Same initial observations and estimator are used across three paired future mechanisms; estimator never receives the mechanism label.

{table([dict(suite=k,**v) for k,v in e['results'].items()],['suite','positives','negatives','balanced_accuracy','brier','viable_mass_mae','nominal90_mass_coverage'])}

The absent hemlock-specific process raises false positives from4 to{e['results']['hemlock_extra_hazard_0.12']['confusion']['FP']}; it degrades Brier score and inflates capacity. Broad annual path-level entry dispersion (log SD1.2 at unchanged mean) also reduces interval coverage. Neither test supports invariance to arbitrary model discrepancy. No benchmark was retuned to improve these outcomes.
"""
    mc=f"""# Monte Carlo resolution

Fixed maximum banks first, then paired K64/256/1024 prefixes on16 B0 states,2 parameter draws,2 model variants,D0/D3,h5/10/20. Separate outer N128/256/512/1024 prefixes at K64 give4096 maximum outer units. Inner and outer studies are NOT a universal N1024×K1024 run.

Selected study: {f['precision']}; maximum-K resolved fraction{f['inner_resolved_fraction']}, maximum interval half-width{f['inner_max_half_width']:.6f}, maximum outer state-block Health SE{f['outer_max_state_block_health_se']:.6f}. Targets were predeclared; no4096 extension or preferred-Boolean stopping.

Wilson score95% Q1 and two-component Bonferroni97.5% Q2 intervals have approximate pointwise binomial coverage, not a simultaneous guarantee across the entire grid. Initial realization is conditional on each sampled state. TRUE requires all lower bounds strictly above theta; FALSE requires absence or an upper bound strictly below; otherwise MC_UNRESOLVED. Strict inequalities mean finite samples cannot confidently equal theta0/1 endpoints.

The broader missingness/semantic surface retains{all_unresolved} unresolved conditional cells, so overall publication precision is {certificate['PUBLICATION_MONTE_CARLO_PRECISION']}, even though the selected maximum-K panel meets its targets. Read [inner](../monte_carlo/inner_resolution.csv), [outer](../monte_carlo/outer_resolution.csv), and portable sufficient counts. Outer state blocks preserve dependence among parameter/model draws. MC intervals do not cover reconstruction-model or ecological uncertainty.

Full-grid resolution by threshold:

{table(unresolved_by_theta,list(unresolved_by_theta[0]))}

Most unresolved labels occur at theta0/1 under the predeclared strict interval rule. At theta0, nonnegative capacity already establishes adequacy algebraically: its conservative MC label is NOT evidence of uncertainty about that algebraic fact. Interior thresholds retain3984 unresolved conditional cases at K256, so the MIXED certificate is not based solely on endpoint artifacts. No post-outcome rule change is used to make the surface appear more resolved.

Method source: [NIST Wilson interval discussion](https://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm).
"""
    uncertainty=f"""# Finite-design uncertainty interpretation

Across the declared finite ranges, between-scenario and between-semantic contrasts were numerically larger than fitted parameter/model contrasts: comparison verified={contrast}. Frozen contrast values are [re-exposed unchanged](../uncertainty/frozen_finite_design_contrasts.csv).

These are noncommensurate, differently conditioned finite-design variances, not independent Sobol components, ecology-wide uncertainty rankings or a causal decomposition. Measurement/date projection and imputation remain coupled. No distributions over scenarios or semantic thresholds are scientifically licensed here; none are invented. No posterior probability over forests, scenario mixtures or semantic specifications is claimed.

M2/M3 remains a worked illustration of query-relative sufficiency inherited from the formal papers, not a new principal empirical discovery. Redundant M3 cannot add information.
"""
    numerical=f"""# Reader-inspectable numerical specification

{table(old['scenarios'],['id','hazard_factor','growth_factor','recruitment_factor','hemlock_extra_hazard'])}

Frozen MNAR log-hazard SD.45, wet-proxy SD.7, unseen entry log SD.35, measurement SD.1cm, synthetic prior DBH SD1cm and observation SD.25cm. Association branches: equiprobable recorded E0/E1 identity alternatives, not confirmed reassociations. Fit:256-point annual hazard grid, survival-prior strength32, effective cap200, six-year candidate-entry exposure, growth inclusion[-.1,1]cm/y, growth means[0,.7], POM tolerance.05m. Origin2020-01-03. R1 uses BA/E0BA and juvenile1–10cm count/E0juveniles; C1 uses annual path excursions plus endpoint realization, not endpoint alone.

Alpha[.5,.75,1], beta[.1,.4,.8], tau[0,2,5]years, theta[0,.5,.75,.9,1], reserve rho[.5,1]. Model-admissibility checks identity/count balance/nonnegative bounded size, not scientific ecological Lawful.

{table(weights,list(weights[0]))}

Each outer unit has equal design weight1/outer_units; each conditional history1/K; their product gives the listed joint weight. Successful histories retain original weights: failures are never dropped or renormalized. Scenario and semantic cells remain separately indexed, not given posterior weights. For each current state Health requires present realization AND capacity adequacy; averaging conditional Health yields a finite-design fraction, not Health of an averaged organization.

Counts/K is the exact finite-atom capacity of the archived bank, with a definite finite-bank Health value. Wilson intervals address sampling resolution as an approximation to the separately declared stochastic demographic kernel; they are not uncertainty in finite-bank arithmetic or in the formal implication itself. Both meanings are kept distinct, and neither is an ecologically established probability law.

Every configuration leaf is exposed in [numerical assumptions](../specification/numerical_assumptions.csv), with both frozen HFD02S and HFD03 protocols. [Editorial map](../specification/editorial_mapping.csv) fixes V2 terminology. Sources: [HF253 publisher archive](https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF253) and the read-only frozen EML/news receipts.
"""
    issues=[
        ("Reconstruction insufficiently validated","Reconstruction_Calibration.md","real E1 masks; honest undercoverage"),
        ("Missingness model carries too much burden","Missingness_Sensitivity.md","B0–B4 full pipeline"),
        ("E2 discrepancy indicates model failure","Generator_Failure_Localization.md","frozen mismatch and diagnostic surfaces, no correction"),
        ("Scenario parameters absent","Numerical_Specification.md","complete numeric tables"),
        ("Noncommensurate uncertainty","Uncertainty_Interpretation.md","finite contrasts only; no invented priors"),
        ("Demonstration semantics mistaken for ecological Health","Numerical_Specification.md","H_S claim boundaries"),
        ("Class-degenerate synthetic Health test","Balanced_Synthetic_Validation.md","324 full-grid worlds with both classes"),
        ("MC resolution inadequate","Monte_Carlo_Resolution.md","paired inner/outer; unresolved status retained"),
        ("Wet terminology overstated","Numerical_Specification.md","protocol proxy, not certified swamp"),
        ("Lawful mostly software invariants","Numerical_Specification.md","model-admissibility realization"),
        ("M3 ablation tautological","Uncertainty_Interpretation.md","supporting query-sufficiency illustration, no gain"),
        ("Ensemble probability unclear","Numerical_Specification.md","explicit finite-design weights"),
        ("Entry identifiability unresolved","Recruitment_Audit.md","complete identity audit and assumption bounds")]
    revision="# Paper V V2 referee-response map\n\n"+table([dict(concern=x,evidence=y,response=z) for x,y,z in issues],["concern","evidence","response"])
    revision+="\n\nV2 should center the auditable observation → reconstruction → model-admissible prospective histories → unconditional capacity → Health-under-specification architecture. Lead empirical limitations with measured calibration and missingness, preserve exposed generator failure, then report synthetic recovery and controlled discrepancy. Formal query sufficiency is supporting. No operational ecological diagnosis or retrospective E2 validation. Accurate narrowed claims are authorized; not every diagnostic needs a favorable outcome.\n"
    closure=f"""# HFD03 Closure

```json
{json.dumps(certificate,indent=2)}
```

All scientific gates A–I derive from executed evidence; gate J must be earned by the standalone verifier before publication. Closure is conditional on that final receipt. The new positive/negative benchmark improves controlled discrimination, but coarse DBH calibration, entry ambiguity and hemlock process discrepancy remain substantive limitations. Selected MC precision is good; full semantic/missingness cells can remain unresolved. Diagnose rather than retrofit.

Read [machine contract](hfd03_contract.json), [referee map](PaperV_Revision_Map.md), nine figure/data pairs and [verification](Verification.md). Original source banks are untouched; no Swift/Lean source changes. Scientific baseline, execution start and publication receipt are distinct.
"""
    reports={"Design":design,"Reconstruction_Calibration":reconstruction,"Missingness_Sensitivity":missing,
        "Recruitment_Audit":recruitment,"Generator_Failure_Localization":failure,
        "Balanced_Synthetic_Validation":synthetic,"Misspecification_Validation":misspec,
        "Monte_Carlo_Resolution":mc,"Uncertainty_Interpretation":uncertainty,
        "Numerical_Specification":numerical,"PaperV_Revision_Map":revision,
        "Verification":"# HFD03 Verification\n\nStandalone verifier recomputes interval labels, prefixes, candidate census, synthetic scores and immutable source identities, then executes16 adversarial controls. Portable audit uses committed sufficient statistics; --cache adds source, sanitized-mask and raw missingness checks; --replay reruns each workstream under frozen protocol, not a new design. Only a PASS receipt earns gate J.\n",
        "Closure":closure}
    for name,text in reports.items():
        markdown("reports/"+name+".md",text)
    summaries={"Design":"Design","Reconstruction_Calibration":"Reconstruction_Calibration",
        "Missingness_Sensitivity":"Missingness_Sensitivity","Recruitment_Audit":"Recruitment_Audit",
        "Generator_Failure":"Generator_Failure_Localization","Synthetic_Validation":"Balanced_Synthetic_Validation",
        "Monte_Carlo":"Monte_Carlo_Resolution","Specification":"Numerical_Specification",
        "Verification":"Verification","Closure":"Closure"}
    for label,target in summaries.items():
        path=REPO/"docs"/("Health_Harvard_Forest_HFD03_"+label+".md")
        path.parent.mkdir(exist_ok=True)
        path.write_text("# HFD03 — "+label.replace("_"," ")+"\n\n"+
            "[Full report](../research/health_formally_defined/harvard_forest/HFD03/reports/"+target+".md)\n\n"+
            "Additive revision against frozen4053fa; real masking calibration, five missingness families, bounded entry audit, untuned E2 diagnosis, balanced synthetic discrimination and explicit MC uncertainty. Operational assessment NOT_AUTHORIZED.\n")
    # Figure1 data keep all strata, not just convenient high-coverage slices.
    write_csv("figures/H03-1.csv",calibration)
    scatter("H03-1","Real E1 masking: observed versus predicted cohort BA",
            [row for row in calibration if row["subgroup"]!="ALL"],
            "ba_truth","ba_prediction",diagonal=True)
    figure2=[dict(family=k,**v) for k,v in b["results"].items()]
    write_csv("figures/H03-2.csv",figure2)
    bars("H03-2","Missingness changes downstream Health labels",[r["family"] for r in figure2],
         [r["point_disagreement"] for r in figure2],"Paired point-label disagreement versus B0; full MC brackets in phase surface")
    figure3=read_csv(ROOT/"recruitment_audit/class_census.csv")
    write_csv("figures/H03-3.csv",figure3)
    bars("H03-3","E1-only candidates: evidence classes",[r["category"].replace("POSSIBLE_","P_") for r in figure3],
         [int(r["n"]) for r in figure3],"Candidate counts; class membership is not confirmed biological recruitment")
    write_csv("figures/H03-4.csv",diagnostics)
    scatter("H03-4","E2 diagnostic grid — approximate public anchors, not validation",
            diagnostics,"predicted_entries","predicted_ancestral_hemlock_deaths",
            [(5000,5000,"Approximate public entries/hemlock anchors; definitions not harmonized")])
    figure5=[r for r in performance if r["suite"]=="correct"]
    write_csv("figures/H03-5.csv",figure5)
    scatter("H03-5","Correct-model reference truth versus estimated Health",
            figure5,"world","estimate_health_probability")
    figure6=[dict(suite=k,**{x:v[x] for x in ("brier","balanced_accuracy","viable_mass_mae","nominal90_mass_coverage")})
             for k,v in e["results"].items()]
    write_csv("figures/H03-6.csv",figure6)
    bars("H03-6","Controlled generator mismatch: Brier score",[r["suite"] for r in figure6],
         [r["brier"] for r in figure6],"Lower is better. Correct model and misspecification never pooled.")
    figure7=[]
    for K in PROTOCOL["monte_carlo"]["K"]:
        selected=[r for r in inner if int(r["K"])==K]
        figure7.append(dict(study="inner",resolution=K,
                           unresolved_fraction=sum(r["status"]=="MC_UNRESOLVED" for r in selected)/len(selected),
                           max_half_width=max(float(r["max_half_width"]) for r in selected)))
    for n in PROTOCOL["monte_carlo"]["outer_N"]:
        selected=[r for r in outer if int(r["N_X"])==n]
        figure7.append(dict(study="outer",resolution=n,
                           unresolved_fraction=float(np.mean([float(r["unresolved_fraction"]) for r in selected])),
                           max_half_width=max(float(r["state_block_health_se"]) for r in selected)))
    write_csv("figures/H03-7.csv",figure7)
    import xml.etree.ElementTree as ET
    top=bars(None,"Inner paired K: maximum conditional interval half-width",
             [str(row["resolution"]) for row in figure7[:3]],[row["max_half_width"] for row in figure7[:3]],
             "Wilson95% Q1 / Bonferroni97.5% Q2; target<=0.04 at maximumK")
    bottom=bars(None,"Outer paired N: maximum state-block Health fraction SE",
                [str(row["resolution"]) for row in figure7[3:]],[row["max_half_width"] for row in figure7[3:]],
                "K64 throughout; state-block SE, not inner interval width. Target<=0.03 at N1024.")
    aplot,bplot=ET.fromstring(top),ET.fromstring(bottom)
    first_height=int(aplot.attrib["height"])
    second_height=int(bplot.attrib["height"])
    parent=ET.Element("{http://www.w3.org/2000/svg}svg",width="920",height=str(first_height+second_height))
    aplot.attrib.update(x="0",y="0")
    bplot.attrib.update(x="0",y=str(first_height))
    parent.extend([aplot,bplot])
    markdown("figures/H03-7.svg",ET.tostring(parent,encoding="unicode"))
    grouped=defaultdict(Counter)
    for row in inner:
        key=(row["K"],row["scenario"],row["horizon"],row["query"],row["rho"],row["theta"])
        grouped[key][row["status"]]+=1
    figure8=[dict(K=k[0],scenario=k[1],horizon=k[2],query=k[3],rho=k[4],theta=k[5],
                  TRUE=v["TRUE"],FALSE=v["FALSE"],MC_UNRESOLVED=v["MC_UNRESOLVED"])
             for k,v in sorted(grouped.items())]
    full_map=[]
    for family in ["B0","B1","B2","B3","B4"]:
        for theta in [0,.5,.75,.9,1]:
            selected=[row for row in phase if row["family"]==family and float(row["theta"])==theta]
            full_map.append(dict(family=family,theta=theta,K=256,
                TRUE=sum(int(row["true_units"]) for row in selected),
                FALSE=sum(int(row["false_units"]) for row in selected),
                MC_UNRESOLVED=sum(int(row["unresolved_units"]) for row in selected)))
    write_csv("figures/H03-8.csv",full_map)
    items=['<svg xmlns="http://www.w3.org/2000/svg" width="920" height="960">',
           '<rect width="100%" height="100%" fill="#faf8f3"/>',
           '<text x="24" y="28" font-family="sans-serif" font-size="17">Full-grid MC resolution: family / threshold</text>',
           '<text x="24" y="53" font-family="sans-serif" font-size="12">Green TRUE, gray FALSE, rust MC_UNRESOLVED. Equal finite query-design counts, not ecological frequencies.</text>']
    for index,row in enumerate(full_map):
        y=80+index*34;x=270
        total=row["TRUE"]+row["FALSE"]+row["MC_UNRESOLVED"]
        items.append(f'<text x="24" y="{y+15}" font-family="sans-serif" font-size="12">{row["family"]} / theta{row["theta"]}</text>')
        for field,color in [("TRUE","#587c59"),("FALSE","#aaa"),("MC_UNRESOLVED","#b56848")]:
            width=560*row[field]/total
            items.append(f'<rect x="{x:.4f}" y="{y}" width="{width:.4f}" height="20" fill="{color}"><title>{field}: {row[field]} / {total}</title></rect>')
            x+=width
    items.append("</svg>")
    markdown("figures/H03-8.svg","\n".join(items))
    source=read_csv(ROOT/"provenance/state_contributions.csv")
    figure9=[]
    for family in b["results"]:
        selected=[r for r in source if r["family"]==family]
        figure9.append(dict(family=family,observed_anchored_ba_fraction=float(np.mean([float(r["observed_anchor_fraction"]) for r in selected])),
                           directly_observed_origin_fraction=0,
                           observed_anchored_count_fraction=float(np.mean([float(r["observed_anchor_count_fraction"]) for r in selected])),
                           observed_anchored_juvenile_fraction=float(np.mean([float(r["observed_anchor_juvenile_fraction"]) for r in selected])),
                           imputed_ancestry_ba_fraction=1-float(np.mean([float(r["observed_anchor_fraction"]) for r in selected]))))
    write_csv("figures/H03-9.csv",figure9)
    bars("H03-9","Origin state: measured ancestry versus imputed ancestry",
         [r["family"] for r in figure9],[r["observed_anchored_ba_fraction"] for r in figure9],
         "Green measured-ancestry/date-projected; rust imputed-ancestry. ALL2020 origin fields remain IMPUTED.",
         [r["imputed_ancestry_ba_fraction"] for r in figure9])
    markdown("figures/README.md","# HFD03 figures\n\nNine reproducible SVG/CSV pairs, generated from committed outcomes. Plots expose selected readable summaries; companion CSVs and reports preserve full surfaces. No raster generator, outcome pruning or favorable axes selection. H03-4 approximate anchors are definitionally incomparable diagnostics. H03-7 compares separate numerical studies. H03-9 distinguishes ancestry from directly observed current state.\n")
    markdown("README.md",f"""# HFD03 — Reconstruction calibration, failure localization and discriminative validation

Additive revision; HFD01/HFD02S frozen at4053fa. [Closure](reports/Closure.md), [contract](reports/hfd03_contract.json), [protocol](config/protocol.json), [referee map](reports/PaperV_Revision_Map.md), [figures](figures/README.md).

Use Python3.9+ with NumPy2.0.2, preferably the existing private HFD02S runtime. No Lean rebuild, new worktree or copied cache required.

```sh
python scripts/verify_hfd03_harvard_forest.py
python scripts/verify_hfd03_harvard_forest.py --cache
python scripts/verify_hfd03_harvard_forest.py --replay
```

Portable audit reads committed evidence. Cache audit reuses read-only HFD01 inputs and HFD03 isolated identity/sanitized caches. Full replay regenerates only HFD03 outputs under the frozen design; no original HFD02S bank is written.

Private generated masked CSVs and compressed missingness sufficient banks are in .runs/{PROTOCOL_HASH[:16]}; public extra identity inputs in .inputs. These are not copied historical scientific banks. Source IDs, namespace/provenance rules, seed algorithm, all candidate IDs, all mask assignments and sufficient MC counts are committed. The manifests bind implementation and output hashes.

Scientific limits are first-class results: coarse DBH undercoverage; entry ambiguity; incomparable exposed E2 anchors; model-dependent discrimination; no mass-near-threshold synthetic cases; unresolved broad-grid cells. Operational Harvard Forest Health remains NOT_AUTHORIZED.
""")
    write_json("provenance/implementation_notes.json",dict(
        protocol_unchanged=True,changes=[
            "B1 vector hazard shifts enter the ordinary unknown-cohort binomial because old scalar logging cannot serialize vectors; 31 identity branches unchanged.",
            "Added annualized growth-error reporting by affine transformation of existing heldout DBH predictions; no model, seed, mask or outcome altered.",
            "Supplemental deterministic MC state-provenance replay covers all1024 outer states; no numerical study rerun or change.",
            "Ancillary E1 AGB/measurement identifiers removed from masked fitting files even though the unchanged fitter never reads them; predictions unchanged.",
            "Added aggregate predictive intervals, readable cohort labels and matched observed-BA scoring population; no latent missing DBH promoted to truth.",
            "Read-only RNG trace adds measured-ancestry count/juvenile decomposition without consuming or changing random draws; coordinate recovery never licenses current BA.",
            "Endpoint projection now applies the frozen parameter sampler's [0,.7] clipping to traced raw growth draws; fixes adapter fidelity, not a new parameter or favorable tuning.",
            "After full model replay, Wilson exact boundary identities L(0,K)=0 and U(K,K)=1 enforced against IEEE roundoff; saved-count interval recount changes no forecast, mass, point Health or parameter.",
            "No statistical parameter changes after observing outcomes."]))
    print("Reports, specification and nine figure/data pairs generated",flush=True)
