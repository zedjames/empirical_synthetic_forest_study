# Controlled misspecification

Do not pool correct-model recovery with robustness. Same initial observations and estimator are used across three paired future mechanisms; estimator never receives the mechanism label.

| suite | positives | negatives | balanced_accuracy | brier | viable_mass_mae | nominal90_mass_coverage |
| --- | --- | --- | --- | --- | --- | --- |
| correct | 69 | 255 | 0.9921568627450981 | 0.006293402777777778 | 0.012590949918016976 | 0.9876543209876543 |
| entry_lognormal_annual_sd_1.2 | 68 | 256 | 0.990234375 | 0.007740162037037037 | 0.014314627941743827 | 0.8827160493827161 |
| hemlock_extra_hazard_0.12 | 39 | 285 | 0.9403508771929825 | 0.09724633487654322 | 0.10705491054205248 | 0.8827160493827161 |

The absent hemlock-specific process raises false positives from4 to34; it degrades Brier score and inflates capacity. Broad annual path-level entry dispersion (log SD1.2 at unchanged mean) also reduces interval coverage. Neither test supports invariance to arbitrary model discrepancy. No benchmark was retuned to improve these outcomes.
