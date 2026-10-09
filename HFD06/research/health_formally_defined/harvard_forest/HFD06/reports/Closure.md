# Terminal scientific findings for Paper V

Both identified computational obligations are complete. The selected state
contrast comes exclusively from the juvenile continuation predicate, and
hemlock mortality remains exceptional after integrating the frozen model's
shared-hazard uncertainty. These numerical findings are entered directly in
the [Paper V revision](../manuscript/Manuscript.md) and its
[Methods](../manuscript/Methods.md). No new ecological model or world was added.

## Terminal failure attribution

State205, parameter replicate1, original Gamma growth, D0, one year,
α0.5/β0.8/τ0 uses the exact first frozen local-precision stream. All paths are
retained; the initial realization and the excursion condition are checked.

| Outcome | Pooled exposure | Individual exposure |
| --- | ---: | ---: |
| Candidate histories | 4,096 | 4,096 |
| Model unlawfulness | 0 | 0 |
| Juvenile support failure only | 4,053 | 137 |
| Basal area failure only | 0 | 0 |
| Joint juvenile and basal area failure | 0 | 0 |
| Successful continuation | 43 | 3,959 |

The five exclusive outcomes reconcile exactly to each unconditional candidate
census. All overlapping raw causes are retained in
[path classifications](../attribution/paths.csv) and
[overlap counts](../attribution/overlapping_causes.csv), including zero joint
causes. The juvenile requirement is65,661.6; origin juvenile counts are66,600
and66,900. Changing reconstruction changes distance to a sensitive juvenile
threshold, not basal-area adequacy or lawfulness. This is a targeted local
result, not a uniform capacity claim or a pathwise causal matching.

## Shared hazard mortality calibration

The frozen strict2021–2024 matched adult cohort contains7,557 hemlock and
14,019 nonhemlock stems. Each predictive replicate samples one hazard per
occupied cell, shared by all its stems, then generates actual individual
mortality outcomes conditional on those hazards. No external outcome refits
the model, adjusts the cap, changes support, or tunes the experiment.

| Adult group | Observed deaths | Predictive mean | Median | Central95% interval | P of at least observed deaths |
| --- | ---: | ---: | ---: | --- | ---: |
| Hemlock | 685 | 280.348 | 279 | 220–347 | 4.6308×10⁻²¹ |
| Nonhemlock | 792 | 769.938 | 768 | 654–896 | 0.353712 |
| All | 1,477 | 1,050.285 | 1,049 | 918–1,192 | 3.2045×10⁻⁸ |

| Adult group | Conditional event variance | Shared hazard variance | Total variance | Shared fraction |
| --- | ---: | ---: | ---: | ---: |
| Hemlock | 266.318 | 782.132 | 1,048.450 | 74.60% |
| Nonhemlock | 693.360 | 3,121.328 | 3,814.688 | 81.82% |
| All | 959.678 | 3,903.460 | 4,863.138 | 80.27% |

Variances have units deaths². Shared hazards widen the model's predictive
distribution materially, yet cannot account for the observed hemlock count.
Nonhemlock does not show the same predictive extremeness. These conclusions
are conditional on the existing cell-independent hazard model and matched
adult support, not population-wide ecological truth or full Health validation.

## Independent numerical evidence

There are262,144 explicit individual-event replicates, with cell hazard
indices, cell death counts and event bitstream digests retained. Exact
hemlock/nonhemlock sums form every all-adult replicate. The independent
Binomial-count sampler adds262,144 numerical-check draws per group under
the same hierarchical law; it does not add worlds or fit a second model.

Mode-centered recurrence with positive direct convolution checks the complete
predictive support independently of loggamma masses and FFT convolution.
Another positive loggamma/direct route checks rare tails to relative error
below2.3×10⁻¹³, much tighter than the registered10⁻⁹ tolerance. The original
FFT hemlock tail is roundoff-dominated and retained only as a diagnostic.
Zero simulation exceedances are retained with Wilson95 upper bound
1.4654×10⁻⁵, not misrepresented as zero probability. The reported small
tail comes from two agreeing positive deterministic calculations.

Scalar checks reconstruct all8,192 terminal histories, while an exhaustive
two-cell individual-event fixture detects marginal per-stem hazard redraw
and omission of shared variance. Four independently replayed chunks reproduce
the ordered individual-event bitstreams. All numerical mean, variance and
distribution checks pass. The pre-execution metadata syntax failure and
outcome-aware rare-tail numerical correction are both preserved in provenance.

## Completion and approval boundary

The scientific stop rule is satisfied by these two results and their numerical
checks, irrespective of their unfavorable findings. The whole-plot juvenile
and entry identification limitations remain scientifically explicit but are
not invitations to expand this tranche. Operational forest Health remains
not established; full prospective ecological validation is not claimed.

The computational artifact is prepared locally for disclosure/source-license
review. No public release, DOI, visibility change or publication is authorized
by this tranche. Source IP and publisher-data terms require separate review.
The existing retained worktree and runtime will be reused; no Lean rebuild,
new cache or new worktree is required.
