# Public data, provenance and support

All external data should be cited to its publisher. Published data access is separate from software source licensing. This repository does not redistribute raw census bytes.

## HF253 — ForestGEO plot census and identity inputs

- Publisher: Harvard Forest Data Archive
- Dataset: HF253, version 6
- DOI: [10.6073/pasta/818789a882a318c1d7f3fc43a2289e12](https://doi.org/10.6073/pasta/818789a882a318c1d7f3fc43a2289e12)
- Archive: https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF253
- Source license recorded in the audited release candidate: CC0 1.0
- Four exact raw input URLs and immutable hashes: [public input manifest](../data/public_inputs.json)

## HF453 — independent adult-tree health/mortality carrier

- Annual adult-tree status records, 2021–2024, within the ForestGEO setting
- Archive: https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF453
- DOI: https://doi.org/10.6073/pasta/0bc6f43d7406c9c01e84a3dd4ed65193

Its strict scored 2021-survivor to 2024 risk set is a conditional adult-stem component test, not an unconditional forecast issued in 2020.

## HF355 — seedling survey

- Separate small-plot seedling observation and graduation carrier, initially below 1 cm diameter
- Archive: https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF355
- DOI: https://doi.org/10.6073/pasta/6c4f69748959fff364ecf89fc293922f

The seedling sampling frame does not identify whole-plot biological threshold entry into the paper's 1–10 cm stem-support proxy.

## Rules

Never silently substitute a different dataset version, unmatched tag, source date, or endpoint crosswalk. Missing/ambiguous status classes must retain their stated exclusion and sensitivity treatment. Raw public inputs are retrieved separately and checked against the cited hashes; the absence of raw Git copies is deliberate.
