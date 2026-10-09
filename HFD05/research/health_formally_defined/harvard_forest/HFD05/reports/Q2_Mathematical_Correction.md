# Exact reserve-query correction

For a nonnegative, unrenormalized history measure of total mass at most one,
the reserve-qualified continuation set is the intersection of viable histories
with the reserve requirement. Set inclusion and measure monotonicity therefore
give `0 <= reserve mass <= viable mass <= 1`. At every inclusive threshold,
requiring both masses is exactly equivalent to requiring reserve-qualified
mass alone. This does not change Health's realization requirement or its
adequacy semantics. The general argument is in
`mathematical_audit/q2_equivalence.md`; finite weighted fixtures supplement,
rather than replace, that argument.

The corrected margin is reserve-qualified mass minus threshold. Its kernel
resolution uses one pointwise Wilson95 interval for the reserve count, not two
Bonferroni components. Certification is strict: lower bound above threshold
for TRUE, upper bound below threshold for FALSE. Present realization and
threshold zero give algebraic adequacy; all-success finite banks at threshold
one do not certify kernel probability one. Confidence is pointwise in a query,
not simultaneous across the entire semantic grid.

## Executed identity and resolution audit

The independent audit checked75,449 archived count cells representing
4,212,549 evaluations. Every finite-bank point label was identical. All4,320
original primary reserve-query surface cells were recomputed directly from
the hash-checked accepted trajectory archive and reproduced the reported
finite Health fractions. All243 original boundary response banks were checked
directly for reserve-support inclusion and archived prefix count agreement.

| Archived family | Evaluations | Old-procedure unresolved | Corrected unresolved |
| --- | ---: | ---: | ---: |
| Full historical missingness census | 4,147,200 | 160,626 | 160,292 |
| Historical synthetic truth | 972 | 0 | 0 |
| Boundary references, all prefixes | 729 | 51 | 48 |
| Current context estimator draws | 51,840 | 253 | 222 |
| Design-prior supplemental draws | 2,592 | 26 | 23 |
| Cap propagation | 6,144 | 2 | 2 |
| Coarse size representation | 1,536 | 0 | 0 |
| Enriched size representation | 1,536 | 2 | 2 |

These old-procedure comparisons use the original two-component procedure with
the subsequently declared algebraic zero-threshold handling. Historical emitted
statuses are retained in a separate column: earlier endpoint conventions are
not overwritten or quietly relabeled. The narrower interval changes numerical
certification, not point Health or the demographic generator.

For the original correct-model PRESENT near-boundary reserve query, six
maximum-bank references remain uncertified. All64 binary truth assignments
are retained; the original point classifier's descriptive balanced-accuracy
identification bracket is0.5 to0.5961538462. This is reference-truth ambiguity,
not an ecological confidence interval or estimator sampling variance.

Evidence: `finite_bank_identity.csv`, `resolution_comparison.csv`,
`primary_surface_identity.csv`, `boundary_rescoring.csv`,
`boundary_margin_strata.csv`, `design_conditional_intervals.csv` and
`reference_ambiguity.json` under `q2_correction/`; independently checked in
`mathematical_audit/q2_independent_audit.json`.
