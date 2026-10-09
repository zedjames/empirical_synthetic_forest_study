# Paper V computational verification

Paper: *Prospective Health under Declared Specifications: An Auditable Empirical-Synthetic Forest Measurement Study* (Health, Formally Defined, Paper V). Paper DOI: https://doi.org/10.5281/zenodo.23268987

## Independently verified scientific reproduction

The frozen studies completed **24 of 24 prescribed replay stages** in the original Apple Silicon ARM64/Darwin 27.0.0 reference environment (Python 3.9.6, NumPy 2.0.2). A separate private inspection confirmed **69 of 69 evidence-file checksums**, the original HFD-05 kernel identity, the original HFD-06 fitted hazard posterior, and **nine of nine** formerly divergent HFD-06 predictive output identities. No scientific source, original inputs, frozen expected results, or acceptance thresholds were changed.

The original independent Intel run exhibited one-ULP floating-point differences. Reference-runtime exact reproduction is established; cross-platform byte identity is not claimed. These computations verify the declared model's implementation and leave empirical/ecological licensing and prospective external validation as separate scientific obligations.

## Package identities

- HFD-05 evidence ZIP SHA-256: `45a10000799de1b7974e650d943adc7859428d830ab8e21f123506d92e715738`
- HFD-06 evidence ZIP SHA-256: `2a36eeccd582a25db9b08e0fb258ddaaaa346e8dfea7809588cf751d3dcbc272`

## Distribution status

The reviewed source and byte-preserved evidence packages are distributed in [v1.0.0](https://github.com/zedjames/empirical_synthetic_forest_study/releases/tag/v1.0.0). [The public verification record](../verification/scientific-verification.json) records the scientific verdict, exact reference identities and environment; [nine-output verification](../verification/nine-output-verification.json) and [runtime comparison](../verification/runtime-comparison.json) preserve exact output identities and portability limits. Private local paths, full runtime fingerprints and private Git ancestry are not part of these public records.

The scientific replay was completed before publication and is not being repeated or represented as a new public reader's experiment. Downloaded release bytes are checked against the frozen originals. Reproduction uses the public ZIPs and hash-pinned public publisher inputs, without the private research repository. See [full instructions](REPRODUCIBILITY.md).
