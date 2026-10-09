"""Record-level independent-source support, before frozen-model scoring."""
import csv,datetime as dt,json,re
from collections import defaultdict,Counter
import h05_common as h
from h05_external_io_v2 import load

def tag(value):
    if re.fullmatch(r"[0-9]+",value):return str(int(value))
    return value
def reconcile(rows,key,dataset):
    grouped=defaultdict(list);excluded=[];accepted=[];bad=set()
    for index,row in enumerate(rows):grouped[key(row)].append((index,row))
    for identity,items in grouped.items():
        distinct={json.dumps(row,sort_keys=True) for index,row in items}
        if len(distinct)>1:
            bad.add(identity[0]);reason="NONIDENTICAL_DUPLICATE_VISIT"
            for index,row in items:excluded.append(dict(dataset=dataset,source_row=index,identity=identity[0],year=identity[1],reason=reason))
        else:
            accepted.append(items[0])
            for index,row in items[1:]:excluded.append(dict(dataset=dataset,source_row=index,identity=identity[0],year=identity[1],reason="EXACT_DUPLICATE_COLLAPSED"))
    return accepted,excluded,bad

def run():
    h.guard("external_alignment");mort,mort_source=load("HF453","hf453-01");seed,seed_source=load("HF355","hf355-02");taxa,taxa_source=load("HF355","hf355-01");canopy,canopy_source=load("HF355","hf355-05")
    _,e0,e1,core,model=h.original.load_data();aliases=defaultdict(set)
    for epoch in [e0,e1]:
        for identity,row in epoch.items():
            for candidate in [row["stem.tag"],row["tag"]+"."+row["stem.tag"]]:
                if candidate not in ["","NA"]:aliases[tag(candidate)].add(identity)
    visits,excluded,bad=reconcile(mort,lambda r:(tag(r["StemTag"]),r["year"]),"HF453")
    bytag=defaultdict(list)
    for index,row in visits:bytag[tag(row["StemTag"])].append((index,row))
    for identity in bad:bytag.setdefault(identity,[])
    cross=[];reversal=set()
    for identity,items in sorted(bytag.items()):
        sequence=sorted(items,key=lambda x:int(x[1]["year"]));dead_seen=False
        for index,row in sequence:
            if row["status"] in ["DC","DS"]:dead_seen=True
            elif dead_seen and row["status"] in ["A","AU"]:reversal.add(identity)
        candidates=sorted(aliases[identity]);stem=candidates[0] if len(candidates)==1 else None;reason="ALIGNED";base=next((r for index,r in items if r["year"]=="2021"),None)
        if identity in bad:reason="NONIDENTICAL_DUPLICATE_VISIT"
        elif identity in reversal:reason="SECURE_DEAD_TO_ALIVE_REVERSAL"
        elif not candidates:reason="NO_RECORDED_STEM_TAG_LINK"
        elif len(candidates)!=1:reason="AMBIGUOUS_STEM_TAG_COLLISION"
        elif stem not in e0 or stem not in e1:reason="NO_TWO_EPOCH_STABLE_IDENTITY"
        elif e0[stem]["tree.id"]!=e1[stem]["tree.id"] or e0[stem]["sp"]!=e1[stem]["sp"]:reason="TREE_OR_TAXON_ASSOCIATION_CONFLICT"
        elif e0[stem]["df.status"]!="alive":reason="E0_NOT_RECORDED_ALIVE"
        elif h.original.empirical.number(e1[stem]["dbh"]) is None or float(e1[stem]["dbh"])<10:reason="E1_ADULT_SIZE_SUPPORT_ABSENT"
        elif base is None or base["status"] not in ["A","AU"]:reason="NOT_OBSERVED_ALIVE_BASELINE2021"
        elif h.original.frozen.fate(e1[stem])==0:reason="PRIOR_SECURE_DEATH_LATER_ALIVE_CONTRADICTION"
        old=e0.get(stem,{});new=e1.get(stem,{});cell=h.original.empirical.cell(old) if old else None
        cross.append(dict(dataset="HF453",external_tag=identity,linked_stem=stem or "",candidate_stems="|".join(candidates),candidate_count=len(candidates),primary_eligible=reason=="ALIGNED",reason=reason,
            raw_E0_stem_tag=old.get("stem.tag",""),raw_E1_stem_tag=new.get("stem.tag",""),raw_species=old.get("sp",""),cell=cell if cell is not None else "",taxon=cell//16 if cell is not None else "",sector_proxy=(cell//4)%4 if cell is not None else "",
            E1_dbh=new.get("dbh",""),E0_quadrat=old.get("quadrat",""),baseline_year=2021,baseline_status=base["status"] if base else "",left_truncated=True))
        if reason!="ALIGNED":
            for index,row in items:excluded.append(dict(dataset="HF453",source_row=index,identity=identity,year=row["year"],reason=reason))
    identities=defaultdict(list);seed_visits,se,bad_seed=reconcile(seed,lambda r:(r["seedlingID"],r["yearOfOb"]),"HF355");excluded.extend(se)
    for index,row in seed_visits:identities[row["seedlingID"]].append((index,row))
    for identity in bad_seed:identities.setdefault(identity,[])
    seed_cross=[]
    for identity,items in sorted(identities.items()):
        reason="ALIGNED_SUBPLOT_PROCESS"
        if identity in bad_seed:reason="NONIDENTICAL_DUPLICATE_VISIT"
        elif len({r["coordinate"] for index,r in items})>1:reason="INCONSISTENT_SUBPLOT"
        elif len({r["taxonCode"] for index,r in items})>1:reason="INCONSISTENT_TAXON"
        elif any(r["status"] in ["ID","IH"] for index,r in items):reason="INVALID_IDENTITY_STATUS"
        first=items[0][1] if items else next(r for r in seed if r["seedlingID"]==identity);eligible=reason=="ALIGNED_SUBPLOT_PROCESS"
        seed_cross.append(dict(dataset="HF355",seedlingID=identity,coordinate=first["coordinate"],taxonCode=first["taxonCode"],visits=len(items),identity_eligible=eligible,reason=reason,
            below1cm_support=True,whole_plot1to10cm_bridge_identified=False,flagged_visits=sum(r["flag"]=="1" for index,r in items)))
        for index,row in items:
            if not eligible:excluded.append(dict(dataset="HF355",source_row=index,identity=identity,year=row["yearOfOb"],reason=reason));continue
            try:date=dt.date.fromisoformat(row["date"])
            except (ValueError,TypeError):excluded.append(dict(dataset="HF355",source_row=index,identity=identity,year=row["yearOfOb"],reason="MISSING_OR_INVALID_OBSERVATION_DATE"));continue
            if date.year!=int(row["yearOfOb"]):excluded.append(dict(dataset="HF355",source_row=index,identity=identity,year=row["yearOfOb"],reason="DATE_YEAR_MISMATCH"))
            if row["sampled"]!="1" or row["alive"] not in ["Y","N"]:excluded.append(dict(dataset="HF355",source_row=index,identity=identity,year=row["yearOfOb"],reason="NOT_SAMPLED_SECURE_ALIVE_DEAD_VISIT"))
    protected=json.loads((h.THREE/"provenance/identity_inputs.json").read_text())["sources"];species_record=next(r for r in protected if r["name"]=="hf253-02-species-codes.csv");species_path=h.THREE/".inputs"/species_record["name"]
    if h.sha(species_path)!=species_record["sha256"]:raise ValueError("Frozen botanical code source changed")
    species={r["sp"]:r for r in h.read_csv(species_path)};botanical=[]
    for row in taxa:
        code=row["taxonCode"].lower();candidate=species.get(code)
        match=candidate is not None and candidate["genus"].strip().lower()==row["genus"].strip().lower() and candidate["species"].strip().lower()==row["species"].strip().lower()
        botanical.append(dict(HF355_taxonCode=row["taxonCode"],HF355_latin=row["latin"],HF253_code=code if match else "",HF253_latin=candidate["latin"] if candidate else "",documented_botanical_match=match,reason="DOCUMENTED_CODE_AND_BOTANICAL_MATCH" if match else "NO_EXACT_DOCUMENTED_MATCH",original_taxon_group=h.original.empirical.taxon(code) if match else "",model_below1cm_hazard_identified=False))
    plot_keys={tag(r["plot"]) for r in canopy};matched=sum(tag(r["coordinate"]) in plot_keys for r in seed_cross)
    matrix=[dict(dataset="HF453",source_unit="Tagged adult>=10cm stems",model_component="Earlier-fit adult survival hazard",alignment="PARTIALLY_ALIGNED_PENDING_SCORING",source_rows=len(mort),source_identities=len(bytag),eligible_identities=sum(r["primary_eligible"] for r in cross),time="July2021 observed survivors to July2022/23/24; year-precision anniversary approximation",obstruction="Excluded collisions/nonadult/unlinked/dead-at-baseline/identity reversals; no2020-to2021 validation",full_Health_identified=False),
        dict(dataset="HF355",source_unit="Below1cm seedlings in1m2 subplots",model_component="Descriptive survival/graduation/first-recorded cohort; not original1-10cm whole-plot J",alignment="DESCRIPTIVE_ONLY",source_rows=len(seed),source_identities=len(identities),eligible_identities=sum(r["identity_eligible"] for r in seed_cross),time="Actual dated repeated seedling visits",obstruction="No identified support/scale bridge or fitted seedling hazard; no whole-plot expansion or germination date",full_Health_identified=False),
        dict(dataset="HF355_CANOPY",source_unit="Recorded subplot canopy photos",model_component="Recorded observational context only",alignment="DESCRIPTIVE_ONLY",source_rows=len(canopy),source_identities=len(plot_keys),eligible_identities=matched,time="Only recorded acquisition year/date; no missing-year extrapolation",obstruction="No learned light-response law; unmatched subplot identifiers preserved",full_Health_identified=False)]
    h.write_csv("external_sources/HF453/tag_crosswalk.csv",cross);h.write_csv("external_sources/HF355/identity_crosswalk.csv",seed_cross);h.write_csv("external_sources/HF355/botanical_crosswalk.csv",botanical)
    h.write_csv("external_sources/exclusion_ledger.csv",excluded,fields=["dataset","source_row","identity","year","reason"]);h.write_csv("external_sources/alignment_matrix.csv",matrix)
    reasons=defaultdict(set)
    for item in excluded:reasons[(item["dataset"],item["source_row"])].add(item["reason"])
    gate=[]
    for dataset,raw,key,yearkey in [("HF453",mort,"StemTag","year"),("HF355",seed,"seedlingID","yearOfOb")]:
        for index,row in enumerate(raw):
            why=reasons[(dataset,index)]
            gate.append(dict(dataset=dataset,source_row=index,raw_identity=row[key],year=row[yearkey],alignment_eligible=not why,reasons="|".join(sorted(why)),scope="Source-support gate, not per-target outcome scoring eligibility"))
    h.write_csv("external_sources/record_gate_ledger.csv",gate)
    h.write_json("external_sources/alignment_certificate.json",dict(status="SUPPORT_COMMITTED_BEFORE_SCORES",sources={r["name"]:r["sha256"] for r in [mort_source,seed_source,taxa_source,canopy_source]},
        HF453=dict(raw_rows=len(mort),unique_tags=len(bytag),eligible_baseline_stems=sum(r["primary_eligible"] for r in cross),reasons=dict(Counter(r["reason"] for r in cross)),reversal_tags=len(reversal)),
        HF355=dict(raw_rows=len(seed),unique_seedling_IDs=len(identities),eligible_identities=sum(r["identity_eligible"] for r in seed_cross),reasons=dict(Counter(r["reason"] for r in seed_cross)),whole_plot_bridge=False),
        record_partition={dataset:dict(raw_n=sum(r["dataset"]==dataset for r in gate),eligible_n=sum(r["dataset"]==dataset and r["alignment_eligible"] for r in gate),excluded_n=sum(r["dataset"]==dataset and not r["alignment_eligible"] for r in gate)) for dataset in ["HF453","HF355"]},
        exclusion_reasons_may_overlap=True,flag1="Metadata inspection request, not automatic death; preserve strata/sensitivity",lineage="Frozen earlier-fit model, later independent public component observations; no external scores or fit yet",scores_computed=False))
    print("External support alignment",[(r["dataset"],r["eligible_identities"]) for r in matrix],flush=True)
if __name__=="__main__":run()
