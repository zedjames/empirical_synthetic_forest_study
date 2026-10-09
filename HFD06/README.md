# Paper V terminal computational review package

This local package reproduces exactly two completed scientific calculations:
selected-state terminal failure attribution and frozen matched adult mortality
calibration with shared hazards. Its current manuscript is
`research/health_formally_defined/harvard_forest/HFD06/manuscript/Manuscript.md`.
It is not a public release or a grant of source redistribution rights.

Use Python3.9.6 and NumPy2.0.2. Commands from the extracted package root:

```sh
python hfd06_runner.py audit
python hfd06_runner.py checks
python hfd06_runner.py retrieve
python hfd06_runner.py attribution
python hfd06_runner.py predictive
```

`audit` verifies exact included bytes. `checks` reconstructs summaries and runs
12 active hostile controls; no raw census retrieval is needed. `retrieve`
fetches only four exact public census/identity files listed in
`Public_Inputs.json`, refuses mismatched bytes, and never downloads whole banks.
`attribution` replays both4,096-path banks and independently classifies every
path. `predictive` replays262,144 explicit individual-event replicates, checks
the independent hierarchical sampler and two positive rare-tail calculations.
All scientific code is byte-preserved. Only canonical public input locations
are adapted; no private Git history or raw location is accessed or claimed.

The earlier full Paper V review artifact and its broader evidence remain
separate. This minimal terminal package contains the dependencies required by
these two obligations, not duplicate caches, private Lean sources, or all older
scientific ensembles. Archive and file hashes are in the inclusion manifest.

The source license and disclosure scope require owner approval. The publisher
CC0 data terms do not confer a license over scientific source IP. No dataset
raw bytes are shipped; exact public locators, versions and hashes are retained.
No public upload, DOI or visibility change is part of preparation.
