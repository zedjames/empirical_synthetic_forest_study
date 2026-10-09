# Effective-count sensitivity

Allcaps50/100/200/500 refit exactly the original E0→E1 likelihood with isolated module globals;200 remains frozen reference. No E2 observations enter fit or select a cap. The cap is a methodological tempering assumption, not an empirically estimated dependence size.

| cap | mean_hazard_width90 | mean_growth_se | mean_mass_change | mean_abs_mass_change | mean_abs_reserve_change | Health_disagreement |
| --- | --- | --- | --- | --- | --- | --- |
| 50 | 0.045667171062356064 | 0.03115298558572168 | -0.0005035400390625 | 0.0005035400390625 | 0.00254058837890625 | 0.001953125 |
| 100 | 0.03951313778722444 | 0.028606177963223216 | 0.000701904296875 | 0.00140380859375 | 0.001312255859375 | 0.0013020833333333333 |
| 200 | 0.034802119211089194 | 0.026920153372541583 | 0.0 | 0.0 | 0.0 | 0.0 |
| 500 | 0.03115214085480423 | 0.025747453052116592 | 0.0026702880859375 | 0.0032196044921875 | 0.00627899169921875 | 0.001953125 |

Selected M0 heldout calibration:

| cap | draws | n | fate_brier | DBH_MAE | DBH_coverage90 | J_prediction | J_truth |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 50 | 64 | 9989 | 0.14360382297714525 | 2.96344695982264 | 0.49316418278113067 | 4411.193660358287 | 4636 |
| 100 | 64 | 9989 | 0.1439516474559527 | 2.966879556160568 | 0.4969823869934721 | 4391.126155177166 | 4636 |
| 200 | 64 | 9989 | 0.14416045170945735 | 2.9695068377969074 | 0.4799852198546619 | 4385.829604005807 | 4636 |
| 500 | 64 | 9989 | 0.14429746459568155 | 2.9701976783244466 | 0.49180933612513855 | 4378.710556661665 | 4636 |

Per-cell256-grid posterior summaries and growthSE, paired origin states, all scenario/horizon/query threshold crossings and MC statuses are in committed tables. Summary averages above are equal finite design contrasts, not ecological priors; scenario-wise detail remains inspectable. Missingness interaction is B0 only. Differences in sampling consumption prevent an exact common-uniform claim. No preferred cap is proposed from these outcomes.
