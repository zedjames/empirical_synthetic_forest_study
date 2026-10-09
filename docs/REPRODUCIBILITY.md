# Reproducibility contract

## Scientific object and claim boundary

Paper V studies a declared prospective Health computation on a 35-ha ForestGEO forest substrate. The published model is a specific empirical-synthetic measurement specification. Exact numerical replay does not establish that its state reconstruction, recruitment semantics, demographic law or Health criteria are ecologically valid.

## Runtime recorded in the research artifacts

- Python 3.9.6
- NumPy 2.0.2
- Deterministic, namespace-indexed PCG64 streams
- Source bytes, configuration objects and output hashes recorded in the scientific inclusion manifests

The published release must include runnable dependencies, not import from the private research checkout.

## Command contract to be enabled after source disclosure

The independently isolated HFD-05 candidate defines:

```bash
python -I artifact_runner.py validate
python -I artifact_runner.py retrieve
python -I artifact_runner.py replay --stage truth
python -I artifact_runner.py replay --stage estimate
python -I artifact_runner.py replay --stage juvenile
python -I artifact_runner.py replay --stage cap
python -I artifact_runner.py replay --stage analysis
python -I artifact_runner.py replay --stage context-prior
python -I artifact_runner.py replay --stage context-prior-analysis
```

The HFD-06 terminal package defines, from its package root:

```bash
python hfd06_runner.py audit
python hfd06_runner.py checks
python hfd06_runner.py retrieve
python hfd06_runner.py attribution
python hfd06_runner.py predictive
```

**These commands are documentation of the intended release interface; the runners are not yet in the public repository. Do not report these steps as externally executed until the files are published and actually run in a clean checkout.**

The lightweight validation verifies archived numeric identities; the broader stages regenerate simulation results. The HFD-06 `attribution` step checks 8,192 selected candidate trajectories, and `predictive` evaluates a shared-hazard adult mortality distribution using the strict matched risk set.

## Input acquisition

The authoritative input allowlist is [data/public_inputs.json](../data/public_inputs.json). Retrieve its four publisher-hosted HF253 sources directly, verify exact SHA-256 before parsing, and keep raw files out of Git. HF453 and HF355 supply separately published adult and seedling observation carriers with their own sampling frames and alignment requirements. See [data sources](DATA_SOURCES.md).

## Expected scientific verification boundaries

- The one-year selected path audit classifies 4,053 pooled-exposure failures and 137 individual-exposure failures as terminal juvenile-support failures, without basal-area or model-unlawful failures, within two 4,096-history conditional banks.
- In the 2021-survivor to 2024 adult risk-set test, hemlock deaths: observed 685; fitted mean approximately 280.348; shared-hazard predictive central interval 220–347; upper-tail probability approximately 4.6308e-21.
- These numbers are **reported paper findings**, not newly re-executed by the present metadata repository.

The paper and the computational artifact should identify the same strict risk set, input data versions, code modules, path weights and result identities. Golden files must be marked as reference outputs, while independent regeneration must use the source and public inputs supplied.
