# Prospective Health under Declared Specifications

**An Auditable Empirical-Synthetic Forest Measurement Study**  
Zed James · *Health, Formally Defined*, Paper V (2026)  
**Paper DOI:** [10.5281/zenodo.23268987](https://doi.org/10.5281/zenodo.23268987)

This repository is the public companion for a measurement-science study using Harvard Forest ForestGEO records, explicit stochastic reconstruction, prospective continuation capacities, synthetic decision tests, and an independent adult-mortality component evaluation.

## Publication status

**The reviewed computational artifact is published as [v1.0.0](https://github.com/zedjames/empirical_synthetic_forest_study/releases/tag/v1.0.0).** Original scientific verification confirmed 24/24 ARM stages, 69/69 evidence checksums, and all nine formerly divergent predictive outputs. The public source contains every scientific Python module and the bounded configurations, documentation and reference results; the two frozen release ZIPs contain the complete evidence needed for reproduction. See [Release status](docs/RELEASE_STATUS.md) and [verification records](verification/scientific-verification.json).

The underlying research repository remains private. This repository contains only the reviewed scientific export and public documentation. No private Git history, unrelated formal code, confidential files, credentials, or cached raw census datasets are included.

## Reproduce the artifact

Download the two complete ZIPs and their checksums from [v1.0.0](https://github.com/zedjames/empirical_synthetic_forest_study/releases/tag/v1.0.0):

```bash
gh release download v1.0.0 --repo zedjames/empirical_synthetic_forest_study \
  --dir downloaded --pattern 'PaperV_*.zip' --pattern SHA256SUMS
python3 scripts/verify_release_assets.py downloaded
unzip downloaded/PaperV_HFD05_v1.0.0.zip -d hfd05-run
unzip downloaded/PaperV_HFD06_v1.0.0.zip -d hfd06-run
```

Run the [full reproduction commands](docs/REPRODUCIBILITY.md) from the separate extracted package roots. The exact-byte reference environment is native ARM64, Darwin 27.0.0, Python 3.9.6 and NumPy 2.0.2. Other platforms can differ numerically; Intel exact-byte identity is not claimed. Public publisher inputs are retrieved by the package runners and checked against fixed hashes. No private repository is required.

The Git source directories `HFD05/` and `HFD06/` retain files up to 1 MiB. They are a browseable source and reference snapshot, not complete execution directories: 35 larger evidence files are intentionally provided only inside the two ZIPs. [SOURCE_INVENTORY.json](SOURCE_INVENTORY.json) records the full partition; [source review](verification/source-review.json) records all 536 candidate-file dispositions and resolution of the 98 review flags.

## The scientific claim

The study evaluates a *declared* prospective Health specification and investigates which pieces of that computation are supported by observations, mathematical invariants, model calibration and finite numerical evidence. It **does not** establish an operational ecological Health diagnosis for Harvard Forest.

Selected results include identification limits for census-derived recruitment, partial reconstruction calibration, near-boundary synthetic decision failures, exact juvenile-support failure attribution in a selected one-year model experiment, and independent evidence of underpredicted adult hemlock mortality.

## Navigate

- [Citation and manuscript identity](CITATION.cff)
- [Licensing by material type](LICENSE.md)
- [Reproducibility contract](docs/REPRODUCIBILITY.md)
- [Verified scientific reproducibility and reference runtime](docs/SCIENTIFIC_VERIFICATION.md)
- [Public-source provenance and licenses](docs/DATA_SOURCES.md)
- [Scientific evidence index](docs/EVIDENCE_MAP.md)
- [Release status and disclosure gates](docs/RELEASE_STATUS.md)
- [GitHub Releases + Zenodo complete-evidence procedure](docs/GITHUB_ZENODO_RELEASE.md)
- [Manuscript / source correspondence](docs/MANUSCRIPT_MAP.md)

## Licensing

As in the public RMMO research series, [LICENSE.md](LICENSE.md) specifies **CC BY-NC-ND 4.0** for the manuscript, documentation, figures, and original derived result tables, and **PolyForm Noncommercial 1.0.0** for released paper-specific executable scientific code and repository build/audit helpers. These terms apply to the reviewed source and the two unchanged release ZIPs; third-party materials retain their own terms. [Third-party notices](THIRD_PARTY_NOTICES.md) preserve the publisher attribution and CC0 exceptions.

Cite the **paper DOI** above for the scientific findings, and cite the computational artifact's `v1.0.0` tag and public commit separately for the released implementation.
