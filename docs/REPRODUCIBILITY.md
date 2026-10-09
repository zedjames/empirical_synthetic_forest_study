# Reproducibility contract

## Scientific object and claim boundary

Paper V studies a declared prospective Health computation on a 35-ha ForestGEO forest substrate. The published model is a specific empirical-synthetic measurement specification. Exact numerical replay does not establish that its state reconstruction, recruitment semantics, demographic law or Health criteria are ecologically valid.

## Runtime recorded in the research artifacts

- Python 3.9.6
- NumPy 2.0.2
- Deterministic, namespace-indexed PCG64 streams
- Source bytes, configuration objects and output hashes recorded in the scientific inclusion manifests

The published packages include their scientific dependencies and do not import from the private research checkout. The exact-byte reference is native ARM64, Darwin 27.0.0, Python 3.9.6 / NumPy 2.0.2, with Clang 21.0.0. Numerical differences on Intel are retained and documented; matching Python/NumPy version strings alone does not guarantee cross-platform byte identity.

## Download and verify the complete packages

Use the two frozen ZIPs from [v1.0.0](https://github.com/zedjames/empirical_synthetic_forest_study/releases/tag/v1.0.0), not the partial browseable Git source roots. From a public checkout:

```bash
gh release download v1.0.0 --repo zedjames/empirical_synthetic_forest_study \
  --dir downloaded --pattern 'PaperV_*.zip' --pattern SHA256SUMS
python3 scripts/verify_release_assets.py downloaded
unzip downloaded/PaperV_HFD05_v1.0.0.zip -d hfd05-run
unzip downloaded/PaperV_HFD06_v1.0.0.zip -d hfd06-run
```

The verifier checks both archive hashes and all 492 HFD05 / 67 HFD06 manifest entries. Keep separate working copies: scientific stages regenerate their outputs. Use a Python 3.9.6 environment with `numpy==2.0.2`, and run each command from the corresponding extracted root. No other Python package is required. `curl` is needed for the explicit publisher retrieval stages. Full stages take substantial computation.

## Full reproduction commands

The independently isolated HFD-05 candidate defines:

```bash
# In hfd05-run, with the reference Python environment activated:
python -I hfd05_runner.py audit
python -I artifact_runner.py validate
python -I hfd05_runner.py tier1
python -I hfd05_runner.py tier2
python -I artifact_runner.py retrieve
python -I artifact_runner.py replay --stage truth
python -I artifact_runner.py replay --stage estimate
python -I artifact_runner.py replay --stage juvenile
python -I artifact_runner.py replay --stage cap
python -I artifact_runner.py replay --stage analysis
python -I artifact_runner.py replay --stage context-prior
python -I artifact_runner.py replay --stage context-prior-analysis
python -I hfd05_runner.py retrieve
python -I hfd05_runner.py tier3 --stage external-audit
python -I hfd05_runner.py tier3 --stage external-scoring
python -I hfd05_runner.py tier3 --stage oracle
python -I hfd05_runner.py tier3 --stage references
python -I hfd05_runner.py tier3 --stage context
python -I hfd05_runner.py tier3 --stage figures
```

The HFD-06 terminal package defines, from its package root:

```bash
# In hfd06-run:
python -I hfd06_runner.py audit
python -I hfd06_runner.py retrieve
python -I hfd06_runner.py checks
python -I hfd06_runner.py attribution
python -I hfd06_runner.py predictive
```

These are the same 24 prescribed stages already certified in the isolated reference environment. Publication transfers the identical scientific bytes and verifies the public distribution; it does not assert a new independent reader's scientific replay or rerun the studies.

The lightweight validation verifies archived numeric identities; the broader stages regenerate simulation results. The HFD-06 `attribution` step checks 8,192 selected candidate trajectories, and `predictive` evaluates a shared-hazard adult mortality distribution using the strict matched risk set.

## Input acquisition

The authoritative input allowlist is [data/public_inputs.json](../data/public_inputs.json). Retrieve its four publisher-hosted HF253 sources directly, verify exact SHA-256 before parsing, and keep raw files out of Git. HF453 and HF355 supply separately published adult and seedling observation carriers with their own sampling frames and alignment requirements. See [data sources](DATA_SOURCES.md).

## Expected scientific verification boundaries

- The one-year selected path audit classifies 4,053 pooled-exposure failures and 137 individual-exposure failures as terminal juvenile-support failures, without basal-area or model-unlawful failures, within two 4,096-history conditional banks.
- In the 2021-survivor to 2024 adult risk-set test, hemlock deaths: observed 685; fitted mean approximately 280.348; shared-hazard predictive central interval 220–347; upper-tail probability approximately 4.6308e-21.
- These are reference findings from the certified frozen implementation, not ecological validation or a new experiment performed during publication.

The paper and the computational artifact should identify the same strict risk set, input data versions, code modules, path weights and result identities. Golden files must be marked as reference outputs, while independent regeneration must use the source and public inputs supplied.
