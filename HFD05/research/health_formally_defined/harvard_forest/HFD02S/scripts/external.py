"""Exposure-aware check only. Called after primary result artifact is frozen."""
import html
import re
from contracts import HFD01, sha


def external_check(sources, primary_path):
    if not primary_path.is_file() or primary_path.name!="health_phase.csv":
        raise ValueError("Primary results must already be frozen")
    primary_digest = sha(primary_path)
    source = sources["HF-third-census-news"]
    path = HFD01/source["local_cache_path"]
    if sha(path)!=source["sha256"]:
        raise ValueError("Exposed public aggregate source changed")
    text = html.unescape(re.sub("<[^>]+>"," ",path.read_text()))
    text = re.sub(r"\s+"," ",text)
    # All figures below are already-exposed published aggregates, not fit targets.
    for fragment in ["11,800","5,000","8000"]:
        if fragment not in text:
            raise ValueError("Published aggregate wording requires renewed audit: "+fragment)
    return dict(status="EXPOSURE_AWARE_EXTERNAL_CONSISTENCY_CHECK",source_id=source["source_id"],
                source_sha256=source["sha256"],primary_result_sha256=primary_digest,
                exposed_before_design=True,training_used=False,adjustments_from_check=False,
                measurements_over=70000,new_threshold_stems_approximately=5000,
                dead_stems_since_2020_over=11800,hemlock_dead_since_2020_almost=5000,
                hemlock_dead_since_initial_2014_over=8000,
                comparison="Different coverage, raw-status interpretation and elapsed-time baseline prevent direct score or validation",
                retrospective_validation="NOT_CLAIMED")
