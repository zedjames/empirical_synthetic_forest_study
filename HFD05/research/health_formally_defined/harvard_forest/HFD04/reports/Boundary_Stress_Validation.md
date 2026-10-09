# Adequacy-boundary stress validation

All81 worlds retained:72 targets and9 realization controls; no confirmatory selection. Pilot001 Q2 count-scaling failed (reserve mass one); pilot002 independently mapped juvenile diameter to reserve mass. Pilot outcomes are excluded from performance.

| target_cells | within_0_1 | within_0_02 | reference_unresolved |
| --- | --- | --- | --- |
| 72 | 67 | 23 | 6 |

| K | cells | resolved | max_interval_width |
| --- | --- | --- | --- |
| 4096 | 72 | 59 | 0.030437935527116045 |
| 16384 | 72 | 62 | 0.015207585964144954 |
| 65536 | 72 | 66 | 0.007602584706981608 |

Primary near-boundary PRESENT-world performance, correct kernel family/full context:

| query | n | certified_n | positives | negatives | FP | FN | balanced_accuracy | viable_MAE | reserve_MAE | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Q1 | 37 | 31 | 14 | 17 | 10 | 4 | 0.5630252100840336 | 0.25089325775971283 | 0.2512965846706081 | PARTIAL |
| Q2 | 52 | 46 | 22 | 24 | 17 | 4 | 0.5549242424242424 | 0.1961951622596154 | 0.22428483229417068 | PARTIAL |

All-suite/full-context results (includes separate realization failures, so not a substitute for the primary boundary assessment):

| suite | query | n | certified_n | FP | FN | balanced_accuracy | reference_unresolved_fraction | estimator_unresolved_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| correct | Q1 | 81 | 75 | 20 | 4 | 0.7082111436950147 | 0.07407407407407407 | 0.024691358024691357 |
| correct | Q2 | 81 | 75 | 27 | 4 | 0.6519607843137255 | 0.07407407407407407 | 0.04938271604938271 |
| hemlock_extra_hazard_0.12 | Q1 | 81 | 81 | 34 | 0 | 0.734375 | 0.0 | 0.024691358024691357 |
| hemlock_extra_hazard_0.12 | Q2 | 81 | 81 | 51 | 0 |  | 0.0 | 0.04938271604938271 |
| entry_lognormal_annual_sd_1.2 | Q1 | 81 | 77 | 21 | 4 | 0.7072230014025245 | 0.04938271604938271 | 0.024691358024691357 |
| entry_lognormal_annual_sd_1.2 | Q2 | 81 | 76 | 34 | 4 | 0.5942173479561316 | 0.06172839506172839 | 0.04938271604938271 |

Signed and absolute strata—including empty and unfavorable strata—are in signed_margin_scores.csv and absolute_margin_scores.csv. Failure scores distinguish REALIZATION, Q1_CAPACITY, Q2_RESERVE. Uncertain reference truths remain in every world table but are excluded explicitly from certified-truth confusion/Brier scoring. Finite-bank truth is separately available. Conditional error is an experimental frequency, not ecological risk. Parameter-moment normalization is deterministic design, not truth fitting; correct kernel family does not guarantee a correct latent-state prior.
