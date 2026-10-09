# Generator diagnostic interface

The diagnostic generator is implemented in
[h03_diagnostics.py](../scripts/h03_diagnostics.py). Its factor grid is frozen in
[protocol.json](../config/protocol.json), under `diagnostics`.

The shared, untuned outputs are
[response_surfaces.csv](../e2_failure_localization/response_surfaces.csv) and
[the failure-localization report](../reports/Generator_Failure_Localization.md).
This directory does not duplicate trajectories or introduce a corrected
primary generator.

The factors are candidate-entry interpretation, hemlock excess hazard, general
hazard, wet-candidate/protocol-proxy state counts, exposure duration, and initial
size resolution. One-factor perturbations and the declared joint
entry/hemlock/general-hazard grid are retained in full.

Two-point size splitting is a diagnostic of aggregation sensitivity, not a
newly validated individual-tree model. It cannot create a size-dependent
mortality process absent from the frozen kernel. Ancestral hemlock mortality
excludes deaths of subsequently added cohorts; that scope differs from the
original cumulative hemlock comparison.

All public E2 anchors remain exposed, approximate, and definitionally
incomparable. No parameter set is selected to minimize their discrepancy, and
the original D0 result is never replaced.
