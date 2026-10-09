# Balanced synthetic Health discrimination

New benchmark `HF-SYNTH-DISC-001`; old all-negative HFD02S test untouched. Full3×3×3×3×4 grid gives324 independent world states. Truth census persisted before estimator performance: 69 positive and 255 negative primary reference cases. No selection by truth Health; misspecification suites paired on the same states, not additional independent state replicates.

Nonvacuity: 139 primary worlds presently realize the organization; 70 of them fail prospective adequacy. The test therefore exercises both halves of Health, not only a current-state realization label.

Correct-model metrics:

```json
{
  "balanced_accuracy": 0.9921568627450981,
  "brier": 0.006293402777777778,
  "calibration": [
    {
      "lower": 0,
      "n": 248,
      "predicted": 0.00031502016129032257,
      "truth": 0.0
    },
    {
      "lower": 0.2,
      "n": 3,
      "predicted": 0.234375,
      "truth": 0.0
    },
    {
      "lower": 0.6,
      "n": 5,
      "predicted": 0.68125,
      "truth": 0.2
    },
    {
      "lower": 0.8,
      "n": 68,
      "predicted": 1.0,
      "truth": 1.0
    }
  ],
  "clustering": "324 independent world states; suites paired on same state, not 972 independent worlds",
  "confusion": {
    "FN": 0,
    "FP": 4,
    "TN": 251,
    "TP": 69
  },
  "negatives": 255,
  "nominal90_mass_coverage": 0.9876543209876543,
  "positives": 69,
  "sensitivity": 1.0,
  "specificity": 0.984313725490196,
  "threshold_adjacent_mae": null,
  "threshold_adjacent_n": 0,
  "unresolved_truth_cells": 0,
  "viable_mass_bias": 0.011253827883873456,
  "viable_mass_mae": 0.012590949918016976,
  "worlds": 324
}
```

The estimator receives only strict noisy/thinned observation packets and public imposed pressure/recovery context labels; no actual generating rates, latent state, structural/regenerative level or truth stream. Its state anchor is fixed public E0. Truth and estimator rates share a declared family in the correct suite, not hidden point parameters.

“Correct model” refers to the demographic response kernel, rate-family and observation operator. The fixed E0 state prior is not claimed to be the generating distribution of the Cartesian synthetic state grid.

Truth is finite4096-path reference with Wilson status, not exact infinite-law Health. Estimator uses64 observation-conditioned draws ×256 futures. Nominal90% mass intervals mix state/parameter uncertainty with finite predictive sampling and are not pure MC intervals.

Under the predeclared discrimination gate the status is ESTABLISHED. This is controlled recovery under a known model class, not real-forest validation. No primary viable masses lie within0.1 of theta: mass-boundary behavior remains untested despite current-state boundary examples. All misses/false positives remain in [performance](../synthetic_validation/performance.csv).
