# Public artifact status and scope

## Current public artifact

The reviewed scientific source and complete evidence assets are published as [v1.0.0](https://github.com/zedjames/empirical_synthetic_forest_study/releases/tag/v1.0.0). The existing landing page, manuscript identifier, source references, evidence map and publication documentation remain available. Reproduction uses the complete release ZIPs plus separately retrieved public inputs; the browseable Git source snapshot intentionally omits 35 evidence files larger than 1 MiB.

The underlying research repository remains private. Scientific verification certified 24/24 exact ARM stages, 69/69 evidence checksums and closure of all three original reproducibility blockers. The release uses those same immutable packages, not newly generated results. It does not assert that a third-party public reader has executed the science. See [scientific verification](SCIENTIFIC_VERIFICATION.md).

## Scientific export requirements

The source export is **deny-by-default**. It may include only explicitly authorized reproducibility files after owner review of each filename, its content, license, upstream dependencies, and provenance. Git history of the private repository is out of scope.

The authorized release should provide: exact source modules and configuration, environment pin, seed/config manifests, independent numerical checkers, source input URLs and SHA-256 hashes, published numeric outputs and their checksums, and an independently runnable package for principal results.

The repository's [licensing policy](../LICENSE.md) follows the RMMO split: **CC BY-NC-ND 4.0** for manuscripts, documentation, figures and original derived tables; **PolyForm Noncommercial 1.0.0** for published project-authored executable code and workflow helpers. The owner explicitly authorized publication following final security and licensing inspection. All 536 prepared candidate files were inspected and all 98 recorded attention flags resolved; [the final dispositions](../verification/source-review.json) distinguish unchanged scientific files, preserved public wrappers and the excluded obsolete preparation-status file. CC0 permissions for external factual fields and publisher metadata remain separate.

## Release safeguards

1. Exact candidate and scientific package inventories constrain every published file.
2. Security and licensing inspection excludes private paths, credentials, raw census copies and unrelated source. Historical guard code and source-hash ledgers are retained as provenance, not a requirement for private Git access.
3. The already certified public-package-only replay establishes reference-runtime execution without private imports.
4. Downloaded release ZIPs must match both original archive identities and all scientific member hashes.
5. The public-only commit and `v1.0.0` release identify the implementation; private ancestry is never pushed.

The scientific manuscript and existing paper DOI are unchanged. No artifact DOI or archival deposit is asserted by this GitHub publication.

The scientific claims remain conditional on the declared model and observation support; a downloadable computational artifact is not an ecological validation certificate.
