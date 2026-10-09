"""All candidate stems retained; identity evidence does not certify recruitment."""
from collections import Counter, defaultdict
from h03_common import ROOT, read_csv, empirical, frozen, PROTOCOL, write_csv, write_json


def run(e0,e1,trees):
    cfg = PROTOCOL["entry_audit"]
    species={r["sp"]:r for r in read_csv(ROOT/".inputs/hf253-02-species-codes.csv")}
    categories = list(cfg["reference_class_weights"])
    by_tag = defaultdict(list)
    for old in e0.values():
        for field in ("tag","stem.tag"):
            if old[field] not in {"NA",""}:
                by_tag[(field,old[field])].append(old)
    rows = []
    for key in sorted(e1,key=int):
        new = e1[key]
        if key in e0 or frozen.fate(new)!=1 or empirical.number(new["dbh"]) is None:
            continue
        parent = trees.get(new["tree.id"])
        diameter = float(new["dbh"])
        c = empirical.cell(new)
        matches = {r["stem.id"] for f in ("tag","stem.tag") for r in by_tag.get((f,new[f]),[])}
        conflicting = any(e0[k]["tree.id"]!=new["tree.id"] for k in matches)
        coordinate_missing = any(empirical.number(new[f]) is None for f in ("gx","gy","quadrat"))
        if conflicting:
            category = "IDENTITY_OR_TAGGING_AMBIGUITY"
        elif coordinate_missing or empirical.sector(new)==3:
            category = "SAMPLING_FRAME_DIFFERENCE"
        elif parent and parent["df.status"]=="prior" and diameter<cfg["small_dbh_cm"]:
            category = "CERTAIN_OR_HIGH_CONFIDENCE_THRESHOLD_ENTRY"
        elif diameter>=cfg["large_dbh_cm"]:
            category = "POSSIBLE_E0_MISS"
        elif diameter<cfg["small_dbh_cm"]:
            category = "POSSIBLE_THRESHOLD_ENTRY"
        else:
            category = "UNRESOLVED"
        low = float(category=="CERTAIN_OR_HIGH_CONFIDENCE_THRESHOLD_ENTRY")
        ref = cfg["reference_class_weights"][category]
        rows.append(dict(stem_id=key,tree_id=new["tree.id"],tag=new["tag"],stem_tag=new["stem.tag"],
                         taxon=new["sp"],gx=new["gx"],gy=new["gy"],quadrat=new["quadrat"],dbh=diameter,
                         publisher_latin=species.get(new["sp"],{}).get("latin","UNLISTED"),
                         publisher_synonym=species.get(new["sp"],{}).get("syn","UNLISTED"),
                         pom=new["pom"],hom=new["hom"],date=new["exact.date"],cell=c,
                         e0_parent_status=parent["df.status"] if parent else "ABSENT",
                         e0_parent_taxon=parent["sp"] if parent else "ABSENT",
                         e0_parent_gx=parent["gx"] if parent else "ABSENT",
                         e0_parent_gy=parent["gy"] if parent else "ABSENT",
                         e0_parent_quadrat=parent["quadrat"] if parent else "ABSENT",
                         e0_parent_dbh=parent["dbh"] if parent else "ABSENT",
                         e0_parent_pom=parent["pom"] if parent else "ABSENT",
                         e0_parent_hom=parent["hom"] if parent else "ABSENT",
                         e0_parent_date=parent["exact.date"] if parent else "ABSENT",
                         e0_parent_latin=species.get(parent["sp"],{}).get("latin","UNLISTED") if parent else "ABSENT",
                         e0_tag_matches=";".join(sorted(matches,key=int)),category=category,
                         entry_low=low,entry_reference=ref,entry_high=1,
                         provenance="OBSERVED fields; DERIVED classification; EXTERNALLY_CONSTRAINED memberships"))
    if len(rows)!=6992:
        raise ValueError("Candidate population changed")
    census = Counter(r["category"] for r in rows)
    write_csv("recruitment_audit/candidate_audit.csv",rows)
    write_csv("recruitment_audit/class_census.csv",[
        dict(category=c,n=census[c],reference_weight=cfg["reference_class_weights"][c]) for c in categories])
    rates = {}
    for bound in ("low","reference","high"):
        rates[bound] = sum(r["entry_"+bound] for r in rows)/cfg["exposure_years"]
    summary = dict(measured_living_candidates=len(rows),class_census=dict(census),
                   annual_candidate_rates=rates,reference_weights=cfg["reference_class_weights"],
                   interpretation=cfg["weight_interpretation"],
                   certain_entry_count=census["CERTAIN_OR_HIGH_CONFIDENCE_THRESHOLD_ENTRY"],
                   caution="Low bound can be zero: no E0 prior-coded root establishes entry. Bounds do not cover unmeasured-region entrants.",
                   e1_tree_table="Unavailable in published HF253 v6; E1 plant-level grouping is DERIVED",
                   e0_tree_table_rows=len(trees),e0_parent_status=dict(Counter(r["e0_parent_status"] for r in rows)))
    write_json("recruitment_audit/summary.json",summary)
    return summary
