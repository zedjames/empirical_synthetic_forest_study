# Paper V evidence map

The table below is an **index of scientific claims to expected reviewed release files**, not a claim that the underlying modules are already publicly available in this repository. Relative paths name the computational identities preserved in the original isolated candidates.

| Published scientific question | Artifact evidence required for public inspection |
| --- | --- |
| Census-to-state reconstruction, fit and unknown fate | `HFD02S/scripts/empirical.py`, `HFD02S/models/demography.py`, `HFD03/scripts/h03_common.py` and input hashes |
| Corrected individual-exposure model | `HFD05/exposure_model/` and `HFD05/corrected_pipeline/` reports and scripts |
| Gamma parameter-law consistency | `HFD05/gamma_growth/` audit, declared moment tests |
| Reserve query, failure mass and Q1/Q2 inference | `HFD05/q2_correction/`, `HFD05/mathematical_audit/` |
| Main corrected prospective capacity surfaces | `HFD05/corrected_pipeline/primary_surface.csv` and configurations |
| Deliberately near-boundary synthetic test and oracle roles | `HFD04/`, `HFD05/oracle_attribution/` |
| Candidate entry / seedling support nonidentification | `HFD05/entry_missingness/` and source support crosswalks |
| Independent adult mortality matched risk set | `HFD05/external_sources/HF453/` |
| One-year terminal juvenile failure attribution | `HFD06/attribution/partition.csv`, `HFD06/attribution/overlapping_causes.csv` |
| Shared-hazard mortality predictive probability | `HFD06/mortality/final_predictive_summary.csv`, `HFD06/mortality/variance_decomposition.csv` and numerical checks |
| Source disclosure and replay | Audited allowlist, per-file hashes, README, runners and isolated replay receipts |

Paper findings are identified in the manuscript and its DOI. This map becomes an executable source index only when every required referenced path has been published and checked in a public release.

Formal Lean verification of Paper IV is a **distinct** scientific artifact. Private Tier-series Lean code must not be copied into this Paper V repository.
