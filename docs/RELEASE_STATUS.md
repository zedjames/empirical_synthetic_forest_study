# Public artifact status and scope

## Current state: private scientific replay certified; public source release not yet uploaded

This repository currently provides a public and citable landing page, manuscript identifier, source-data retrieval references, methodological evidence map, and the review contract for the future executable source package. **A full independent scientific reproduction is not yet possible from this repository alone.**

The private research repository remains private. The HFD-05 and HFD-06 isolated scientific packages exist in a separately controlled research environment. The private scientific reproducibility record has been independently certified: 24/24 exact ARM scientific stages, 69/69 evidence checksums, and closure of all three historical mismatch investigations. See [scientific verification](SCIENTIFIC_VERIFICATION.md). This does not yet establish reproduction by a public reader from the source and evidence packages.

## Scientific export requirements

The source export is **deny-by-default**. It may include only explicitly authorized reproducibility files after owner review of each filename, its content, license, upstream dependencies, and provenance. Git history of the private repository is out of scope.

The authorized release should provide: exact source modules and configuration, environment pin, seed/config manifests, independent numerical checkers, source input URLs and SHA-256 hashes, published numeric outputs and their checksums, and an independently runnable package for principal results.

The repository's [licensing policy](../LICENSE.md) now follows the RMMO split: **CC BY-NC-ND 4.0** for manuscripts, documentation, figures and original derived tables; **PolyForm Noncommercial 1.0.0** for published project-authored executable code and workflow helpers. The intended source candidates were previously internally audited. Their inclusion inventories mark private authored scientific code **not authorized for source redistribution pending file-level owner disclosure decisions**. Selecting the public license does not override that status. The explicit file-review gate remains controlling. CC0 permissions for externally published field data do not license private research software.

## Release gates

1. Owner approves the public file allowlist, use of the documented PolyForm Noncommercial 1.0.0 terms for the approved scientific source, and exclusions.
2. Static checks reject credentials, private research paths, unreviewed private dependencies and disallowed files. This is an additional safeguard, not a replacement for human review.
3. A clean-room replay operates using only the intended public artifact and documented external public inputs.
4. Public release contents and hashes are independently checked against the audited manifest.
5. A publicly accessible immutable commit/tag is verified; only then update the article to identify that **scientific** source commit.
6. If desired, archive that exact public release and cite its separately assigned DOI (none is asserted here).

An empty or metadata-only repository commit is **not** an executable reproduction certificate.
