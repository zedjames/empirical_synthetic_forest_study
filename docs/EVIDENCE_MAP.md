# Paper V evidence map

The scientific source and full frozen evidence packages are now published in [v1.0.0](https://github.com/zedjames/empirical_synthetic_forest_study/releases/tag/v1.0.0). This table indexes the paper's claims to corresponding scientific files within the **extracted HFD-05/HFD-06 release ZIP roots**; the public Git source snapshot omits 35 larger evidence files already present in the complete ZIP assets. Relative names below retain the original study path suffixes. See [reproduction instructions](REPRODUCIBILITY.md) and [source inventory](../SOURCE_INVENTORY.json).

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

Paper findings are identified in the manuscript and its DOI. The complete executable release is tagged at public commit [`a365fbffb6138a0c798276fe7b6f4577e9132305`](https://github.com/zedjames/empirical_synthetic_forest_study/commit/a365fbffb6138a0c798276fe7b6f4577e9132305); both scientific package SHA-256 identities are recorded in [`SHA256SUMS`](../SHA256SUMS). The certified 24-stage scientific replay predates distribution, so the release is not represented as a new third-party public replay.

Formal Lean verification of Paper IV is a **distinct** scientific artifact. Private Tier-series Lean code must not be copied into this Paper V repository.
