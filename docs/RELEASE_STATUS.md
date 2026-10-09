# Public artifact status and scope

## Current state: metadata and provenance landing published; scientific source under disclosure review

This repository currently provides a public and citable landing page, manuscript identifier, source-data retrieval references, methodological evidence map, and the review contract for the future executable source package. **A full independent scientific reproduction is not yet possible from this repository alone.**

The private research repository remains private. The HFD-05 and HFD-06 isolated scientific packages exist in a separately controlled research environment. Their recorded internal test passes are evidence of internal verification, not proof that a public reader has reproduced them.

## Scientific export requirements

The source export is **deny-by-default**. It may include only explicitly authorized reproducibility files after owner review of each filename, its content, license, upstream dependencies, and provenance. Git history of the private repository is out of scope.

The authorized release should provide: exact source modules and configuration, environment pin, seed/config manifests, independent numerical checkers, source input URLs and SHA-256 hashes, published numeric outputs and their checksums, and an independently runnable package for principal results.

The intended source candidates were previously internally audited. Their inclusion inventories mark authored scientific code **not authorized for source redistribution pending owner decisions on license and disclosure**. That status remains controlling until the explicit source-review gate is completed. CC0 permissions for externally published field data do not license private research software.

## Release gates

1. Owner approves the public file allowlist, software license and exclusions.
2. Static checks reject credentials, private research paths, unreviewed private dependencies and disallowed files. This is an additional safeguard, not a replacement for human review.
3. A clean-room replay operates using only the intended public artifact and documented external public inputs.
4. Public release contents and hashes are independently checked against the audited manifest.
5. A publicly accessible immutable commit/tag is verified; only then update the article to identify that **scientific** source commit.
6. If desired, archive that exact public release and cite its separately assigned DOI (none is asserted here).

An empty or metadata-only repository commit is **not** an executable reproduction certificate.
