# Implemented Gamma growth and correction

The original Gamma branch has three numerical floors: the SD in the shape
denominator, shape itself, and mean in the scale denominator. Its actual
conditional mean is shape times scale and variance is shape times scale squared.
They are not generally the nominal growth and variance. The instrumented actual
simulator calls agree exactly with the audited formulas; the60-case moment
grid reports120 original/corrected law cells, including zero means, zero SD,
each floor boundary and conventional growth. Clipping to the diameter range
is audited separately from pre-clipping Gamma moments.

The corrected positive law sets shape to the larger of `(mean/SD)^2` and
`10^-6`, then scale to mean divided by shape. It preserves mean. Its variance
is mean squared divided by shape: it equals requested SD squared only outside
the shape-floor branch. Zero mean gives deterministic zero; zero SD gives
deterministic mean. Positive variance at zero mean is impossible for a
nonnegative random variable. The exact argument and clipping distinction are
in `mathematical_audit/gamma_moments.md`.

## Reached branches

The source-based census reproduces all named confirmatory parameter identities
for the primary, historical missingness, historical synthetic truth and
observation estimator, current boundary/factorial context estimators, cap,
representation and design-prior supplemental families. It counts64 ancestral
cell laws and separately records cohort-year multiplicities. It does not
mislabel a mathematical possibility as an executed branch or interpret an
empty-cohort draw as a live biological stem.

In the primary ensemble, twelve positive shape-floor cell/scenario laws were
reached among131,072 ancestral cell laws. They arise from three parameter
identities across all four scenarios. All twelve exceed the declared1% relative
conditional mean discrepancy. The SD and mean denominator floors were not
reached in these confirmatory parameter families; zero-mean branches were
reached and remain explicitly degenerate. Other estimator/missingness families
also reach positive shape-floor branches; `branch_summary.csv` gives each
family's full denominator and material-mean count.

The defect is **REACHED_MATERIAL_MOMENTS**. For example, one primary nominal
mean0.0000423206 has original actual mean0.0003186235 and variance0.1015209;
the corrected mean is0.0000423206 and variance0.0017910. A rare tiny-shape law
can materially distort its moments without producing frequent finite-bank
events. The moment simulation explicitly marks insufficient relative precision
MC_UNRESOLVED instead of treating many zero draws as proof of correctness.

## Prospective impact

The initial fixed panel has256 horizon cells. A separately registered,
source-detected panel covers every reached primary parameter identity at
16,384 paths,48 horizon cells, with matched state/parameters/context and verified
identical terminal RNG states across the scale-only correction. Neither panel
changes viable/reserve counts or finite Health labels. The complete primary
factorial subsequently preserves all6,480 finite query surfaces and their
local counts for the growth-only correction. Those are executed finite-bank
facts, not proof of identical continuous laws or universal threshold robustness.

The normal-growth alternative remains the original rectified normal law;
its actual moments are not silently equated with the Gamma target moments.
Corrected missingness/context propagation and broader headline qualification
are reported separately in the final pipeline comparison.

Evidence: all six law/branch/impact tables and both panel manifests in
`gamma_growth/`, corrected primary surfaces, and
`mathematical_audit/abc_independent_audit.json`.
