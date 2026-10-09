# Paper V source candidate — REVIEW ONLY

Python3.9.6 and NumPy2.0.2 were used. Install requirements in your own isolated environment.

Lightweight verification: `python -I artifact_runner.py validate` (fixture plus complete published-table audit).
This reads only packaged source/model/table fixtures and verifies independent source-module location, archived kernel hash, exact/kernel endpoints, query-sufficiency counterexample, all published experimental rosters/decisions/confusion scores and export identities. It does NOT establish full empirical raw-data replay or independently verify unavailable privateGit history.

Public-data replay: `python -I artifact_runner.py retrieve`.
Then replay stages truth, estimate, juvenile, cap, analysis, context-prior, context-prior-analysis with `python -I artifact_runner.py replay --stage STAGE`. Truth precedes estimation; context-prior precedes its analysis. The retrieval list contains only four publicly hosted HF253v6 CC0 tables, exact SHA checked. Full bank sizes and all world cells are in config/protocol.json. Raw banks/caches/Git are not distributed. A separate HFD03 decision-census file exposes historical exact counts/statuses; that presentation can be checked without private banks. Public-safe replay summaries document the actually completed full packaged-source replay: 23 original scientific outputs and five separately registered supplemental outputs. Raw bytes were obtained externally earlier and reused read-only, not freshly downloaded during replay. Source/config and output hashes allow comparison against the supplied golden files. Full reproduction regenerates outcomes; validate checks published evidence and does not itself run all banks.

Public adapter replaces only private Git ancestry checking and canonical raw-data location discovery. Exported scientific files are byte-exact; original private-history preservation is an attested hash ledger, not independently reconstructed privateGit history. Actual public source identities are verified against this inclusion manifest.

No unrelated research, Lean formal-world code, credentials or private raw datasets are included. The original repository stays private. HF253 dataset is CC0; cite Orwig/Foster/Ellison2023 HF253v6. Project-authored source disclosure and distribution license require separate owner review. This is a locally validated review candidate, NOT a public release and not aDOI.
