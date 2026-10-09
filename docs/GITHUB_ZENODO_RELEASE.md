# Paper V — GitHub Release and Zenodo archival procedure

This is a release procedure, not a claim that the scientific packages have already been made public. The private development repository and its Git history must remain private.

## Reviewed scientific packages

- Full HFD-05 candidate: **492** manifest-listed files; the separate 229-file Paper V subset does not replace this broader candidate.
- HFD-06 candidate: **67** manifest-listed files, frozen review ZIP SHA-256 \`2a36eeccd582a25db9b08e0fb258ddaaaa346e8dfea7809588cf751d3dcbc272\` (31,157,413 bytes).
- A large HFD-06 derived mortality table is ~104 MB uncompressed. It belongs inside the compressed evidence ZIP, not a normal tracked Git file.

The \`PaperV_Release_Execution_Kit.zip\` in the accompanying author's review package contains \`scripts/package_release.py\`, its unit tests, and exact instructions to build a deterministic HFD-05 ZIP and verify the original HFD-06 ZIP. The audited scientific files themselves must come from the existing *frozen isolated packages*, not a whole-repository export.

## Gates

1. Confirm file-by-file disclosure scope and licensing. \`LICENSE.md\` follows the RMMO split: **PolyForm Noncommercial 1.0.0** on released project-authored executable code; **CC BY-NC-ND 4.0** on manuscript, documentation, figures and original derived tables. Publisher datasets retain original terms.
2. Verify candidate manifests against every source, CSV, NPZ and configuration byte. Deny unexpected files, links, private Git history, credentials, raw publisher input caches and unrelated code.
3. Run isolated scientific replays from *each candidate's own root*, using its documented runner and Python 3.9.6/NumPy 2.0.2. Retrieve only exact publisher HF253 public inputs with hash checks. Both corrected HFD-05 results and HFD-06 attribution/predictive checks must reproduce their frozen controls.
4. Track only the reviewed public source/configurations in Git. Attach the two approved compressed evidence ZIP files plus \`SHA256SUMS\` and the reviewed replay receipt as GitHub Release assets.
5. Download the public assets without private access. Check \`python scripts/verify_release_assets.py DOWNLOAD_DIRECTORY\`; the separate GitHub Actions release-integrity workflow performs the same byte checks. **Neither check substitutes for scientific replay.**
6. Freeze a public source commit and tag. Identify the actually released scientific code commit in the article, not the repository's earlier documentation-only commits.
7. Archive the **complete two-ZIP evidence package manually on Zenodo**, verifying both assets are present. Zenodo's GitHub integration normally archives the repository snapshot and should **not** be presumed to copy attached binary release assets. Source-only and complete-evidence DOIs, if both created, must be described separately. Do not invent an artifact DOI or change the article DOI \`10.5281/zenodo.23268987\`.

### GitHub release command (only after all gates pass)

\`\`\`sh
gh release create v1.0.0 \
  PaperV_HFD05_v1.0.0.zip PaperV_HFD06_v1.0.0.zip \
  SHA256SUMS release-review.json \
  -R zedjames/empirical_synthetic_forest_study \
  --verify-tag --title "Paper V computational evidence v1.0.0" \
  --notes-file release-notes.md
\`\`\`

### Zenodo owner action

Connect the GitHub account under Zenodo **Linked accounts** if necessary, but use **New upload** for the complete evidence packages and hashes; enter software/data provenance, article DOI as a related identifier, exact source tag, and applicable license information. Publishing an archive and assigning its DOI require an authorized account action. A separate optional GitHub integration can archive the source-code snapshot.

**Release incomplete until:** source/asset public retrieval, independently reproduced principal outputs, SHA-256 verification, and complete Zenodo evidence deposit are all confirmed.
