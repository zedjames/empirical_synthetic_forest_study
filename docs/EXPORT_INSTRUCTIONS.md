# Building a disclosure-safe computational release

The prepublication review packages contain research sources and derivative outputs. This public repository must **never** import the private repository's Git history or indiscriminately mirror a research directory.

## Review and export procedure

1. Work in an isolated local folder, outside the private checkout and outside this public repository's worktree.
2. Extract the independently audited HFD-05 and HFD-06 candidates there. Their original inventories record exact per-file SHA-256 hashes. Retain the original archives unchanged for internal evidence.
3. Review each required source, model fixture, numerical result and dependency separately. Reconcile overlapping HFD-05/HFD-06 modules; preserve the original scientific file bytes and explicit result identities. Avoid replacing frozen published outputs by similarly named corrected outputs.
4. Approve a **new** explicitly enumerated release allowlist and source license. Do not carry forward inventory entries marked `NOT_AUTHORIZED_PENDING_DISCLOSURE_AND_LICENSE_APPROVAL` as though they granted redistribution.
5. Exclude any raw HF253 census bytes, privately cached data, credentials, private research Git history, unrelated Lean/Tier modules, unreviewed third-party licensed material and irrelevant research.
6. Run `python scripts/preflight_export.py PATH_TO_EXPORT --inventory PATH_TO_APPROVED_INVENTORY --report audit.json`. Independently review every listed file and the report. A static scan cannot prove source confidentiality or scientific correctness.
7. Recreate an environment using only the proposed public export plus the exact publisher-hosted data retrieved under the [public input manifest](../data/public_inputs.json). Run the full isolated replay commands and compare output SHA-256 to the approved golden results. Verify paths do not resolve into any private checkout.
8. Only after those gates pass, commit the approved export to this repository, run the public replay action, tag the exact scientific commit and update the manuscript's source-citation sentence.
9. Optionally deposit the frozen public release in an archival repository for a distinct computational artifact DOI.

## Important constraints

A reviewer needs executable source, the declared input provenance and outputs sufficient for the paper's **principal** quantitative claims, not every private development artifact. The HFD-05 broad experiment and the HFD-06 terminal checks must remain separately reproducible and their claim boundaries intact.

The original HFD-06 review inventory contains 67 files, including large derived prediction archives. Its current redistribution status is explicitly `NOT_AUTHORIZED_PENDING_DISCLOSURE_AND_LICENSE_APPROVAL`. The source license and publication decision are the owner's. Publishing a public README does not turn this internal package into a licensed archive.

The article's CC BY-NC-ND 4.0 license covers the article as stated in the manuscript. A separate scientific code license must be selected before code disclosure; publisher CC0 data licenses do not cover project-authored code.
