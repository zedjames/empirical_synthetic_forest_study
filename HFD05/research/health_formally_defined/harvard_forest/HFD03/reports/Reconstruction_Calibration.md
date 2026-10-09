# Real-data reconstruction calibration

Eligible pool: 66284 stable tagged transitions. Alive/A and dated stem-dead records corroborated by explicit D codes qualify; Miss/MT/S and ambiguous associations, absent/gone/broken/missing fates do not. Secure death is stem-level, never whole-plant death.

For each mask, all nonidentity measurement fields (including otherwise-unused AGB/measurement-ID fields) are physically removed from a new visible CSV before unchanged HFD02S read_core, fit and reconstruct run. Only stem/tree/tag/taxon identity metadata remains. All fitted priors are rebuilt from visible transitions, not a fit on full E1. Assignment uses only E0 taxon/size/stratum, known campaign month and observed missing-field pattern; the truth evaluator is downstream.

The first reconstruction remains the 2020-origin kernel. Its origin medians are recorded but NOT scored against 2018–2019 observations as though simultaneous. Calibration scores use a tagged-endpoint projection of the exact sampled unknown-fate kernel at actual exposure. This retains cohort RMS DBH, MNAR shifts, 1cm prior spread and original fit; it does not silently replace the estimator with individual-growth inference.

| mask | replicate | n | brier | dbh_mae | coverage_50 | coverage_80 | coverage_90 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M0 | 0 | 9989 | 0.14243739386668147 | 2.9626133669676693 | 0.27429486389949503 | 0.4488237467668432 | 0.49993841606109124 |
| M0 | 1 | 10030 | 0.1416050383076546 | 2.9516727752283036 | 0.26179995108828563 | 0.4486426999266324 | 0.4986549278552213 |
| M0 | 2 | 9839 | 0.14496155186928075 | 3.0421237762461013 | 0.2626021370207417 | 0.4358265241986172 | 0.5027027027027027 |
| M1 | 0 | 13136 | 0.14666626120333848 | 3.280611780203699 | 0.23964975730465404 | 0.40934615018559056 | 0.4678785571523746 |
| M1 | 1 | 13112 | 0.14929190856499963 | 3.262023508614585 | 0.2461435278336687 | 0.41716968477531857 | 0.47398677780971543 |
| M1 | 2 | 13115 | 0.14782836357267612 | 3.235484299117112 | 0.2581537287812321 | 0.4191302689300019 | 0.4734884608048827 |
| M2 | 0 | 11526 | 0.14563983384048781 | 3.05874276694112 | 0.2718113612004287 | 0.4352625937834941 | 0.49110396570203646 |
| M2 | 1 | 11753 | 0.14516287437114989 | 2.9248281233100952 | 0.2741766310048248 | 0.4451436962450178 | 0.5027270820222363 |
| M2 | 2 | 11715 | 0.14771171976620004 | 3.053555066275846 | 0.27878402711577166 | 0.4428556297002436 | 0.4990996716449529 |
| M3 | 0 | 8131 | 0.1560193414293665 | 2.6873606004268513 | 0.3085122817858798 | 0.49019002008342344 | 0.5510582419280087 |
| M3 | 1 | 8151 | 0.15543032480131566 | 2.7021441150801886 | 0.32609030667283095 | 0.5074741870858376 | 0.5604869779627061 |
| M3 | 2 | 8118 | 0.15109401129244743 | 2.758453145999054 | 0.3009230769230769 | 0.5013846153846154 | 0.5595384615384615 |

50/80/90% coverage has no invented pass gate. Individual DBH and annualized growth bias/MAE/RMSE and interval widths, living count/BA/juvenile errors AND predictive intervals, taxon-share errors, fate calibration/Brier/log scores and coordinate recovery are committed. Growth coverage is the same affine interval event as DBH coverage, with annualized widths. Missing alive DBH cannot provide observed BA truth: those few records are excluded at scoring ONLY, with an explicit census, rather than comparing a partial observed total against a full prediction.

Aggregate heldout-cohort interval coverage:

| mask | nominal | cohort_intervals | living_coverage | ba_coverage | juvenile_coverage |
| --- | --- | --- | --- | --- | --- |
| M0 | 50 | 160 | 0.75625 | 0.7 | 0.7875 |
| M0 | 80 | 160 | 0.975 | 0.90625 | 0.93125 |
| M0 | 90 | 160 | 0.99375 | 0.9375 | 0.95625 |
| M1 | 50 | 164 | 0.7804878048780488 | 0.7987804878048781 | 0.7560975609756098 |
| M1 | 80 | 164 | 0.9817073170731707 | 0.9573170731707317 | 0.8963414634146342 |
| M1 | 90 | 164 | 1.0 | 0.9695121951219512 | 0.9207317073170732 |
| M2 | 50 | 168 | 0.7916666666666666 | 0.7142857142857143 | 0.7916666666666666 |
| M2 | 80 | 168 | 0.9642857142857143 | 0.8809523809523809 | 0.9464285714285714 |
| M2 | 90 | 168 | 1.0 | 0.9226190476190477 | 0.9702380952380952 |
| M3 | 50 | 169 | 0.7869822485207101 | 0.7455621301775148 | 0.7514792899408284 |
| M3 | 80 | 169 | 0.9704142011834319 | 0.9171597633136095 | 0.9112426035502958 |
| M3 | 90 | 169 | 0.9881656804733728 | 0.9467455621301775 | 0.9349112426035503 |

These are descriptive pooled cohort/replicate intervals with shared fitted draws; they are not independent binomial coverage trials. The frozen state kernel is cohort RMS, so poor individual-size coverage does not by itself establish poor aggregate BA calibration. Aggregate intervals test the latter separately.

Worst declared cohort slices (minimum20 observed diameters) by 90% coverage:

| mask | replicate | taxon_group | E0_size_class | stratum | dbh_n | dbh_mae | coverage_90 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| M0 | 1 | other | 10-30cm | central | 32 | 5.783724337786948 | 0.03125 |
| M3 | 0 | hemlock | 30+cm | east | 41 | 6.4912514376166985 | 0.04878048780487805 |
| M3 | 0 | hemlock | 30+cm | central | 35 | 6.063733435935564 | 0.05714285714285714 |
| M1 | 1 | other | 10-30cm | west | 114 | 5.985882338253981 | 0.06140350877192982 |
| M2 | 0 | other | 10-30cm | central | 32 | 5.883017565394283 | 0.0625 |
| M3 | 2 | other | 10-30cm | east | 32 | 5.383259506154337 | 0.0625 |
| M3 | 1 | other_canopy | 30+cm | east | 135 | 7.608068263229282 | 0.06666666666666667 |
| M0 | 2 | other_canopy | 30+cm | west | 140 | 8.976210554126704 | 0.07142857142857142 |
| M3 | 1 | other_canopy | 30+cm | west | 82 | 9.215362189423177 | 0.07317073170731707 |
| M0 | 0 | other | 10-30cm | west | 48 | 5.616719709848602 | 0.08333333333333333 |
| M0 | 0 | other | 10-30cm | east | 36 | 5.237615970041017 | 0.08333333333333333 |
| M0 | 2 | other | 10-30cm | west | 36 | 5.780429086804322 | 0.08333333333333333 |

Cell code: (taxon×4+stratum)×4+E0 size, as frozen in HFD02S. Wider size bins show poor individual coverage: cohort RMS is a structural limitation, not evidence that unseen diameters are accurately reconstructed. Leave-observed-out is conditional on securely observed records and cannot validate nonignorable fate in the unobserved proxy region.

Data: [calibration](../reconstruction_validation/calibration.csv), [fate curve](../reconstruction_validation/fate_curve.csv), [coordinate recovery](../reconstruction_validation/coordinate_recovery.csv), [taxon composition](../reconstruction_validation/taxon_composition.csv), [eligible IDs](../reconstruction_validation/eligible_truth_pool.csv), [assignments](../reconstruction_validation/mask_assignments.csv).
