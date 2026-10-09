# Actual-time reconstruction

The original reconstruction pools times in three survivor components:
living measured endpoint stems, unresolved baseline fates and living endpoint
stems without measured size. The audit recovered98,018 dated component records
and reproduced each original cell's count and exposure sum without date
fallbacks. Existing unresolved-identity branches already use individual dates
and retain their separate law. The original fitting likelihood also already
uses individual survival durations; its pooled hazard prior center is a
distinct approximate operation, inventoried without an external refit.

Conditional on the cohort hazard and shared missingness shifts, individual
survival probabilities are exponential in actual exposure. Independent
Bernoulli survivors form a Poisson-binomial count. Its mean is the sum of the
probabilities and variance the sum of probability times its complement.
Convexity makes the individual mean at least the pooled-time mean. This is a
conditional expectation inequality, not a guaranteed ordering of sampled
counts. Mixing shared hazard uncertainty introduces dependence. Neither a
global Poisson-binomial law nor a universal variance ordering is claimed.

`cohort_exposure_inventory.csv` gives the durations' range/variance/quantiles,
source, hazard quantile, shared shift, both survivor means/variances, absolute
and relative Jensen gaps. Independent dynamic-programming distributions and
larger Bernoulli moment fixtures verify the formula and implementation. The
general proof is `mathematical_audit/unequal_exposure.md`.

## State comparisons and controlled interpretation

All256 original primary state identities and128 identities in each of the five
original missingness families are compared. The individual-time carrier keeps
the64 ancestral cohorts, shared component-size noise and baseline RMS sizes;
each record's actual exposure enters survival and projected diameter, then
survivor counts and squared sizes are reaggregated. Observed DBHs, dates or
positions are not invented. Unknown fate does not become death.

The first variable-consumption RNG comparison remains reported. A subsequent
controlled comparison holds the original parameters, missingness shifts,
size-noise, unseen entry and association draws fixed. Only dated survivor laws
and their growth exposures change; sampled Bernoulli survivors still have MC
variation. The controlled primary mean count change is+26.34375, with absolute
95th percentile286.25 and maximum471; mean basal-area change0.02128165 has
absolute95th percentile3.3601542 and maximum5.4908329. Juvenile mean change
is+24.375,95th percentile260.75 and maximum465. These local quantities are not
replaced by averages or declared deterministic bias estimates.

All taxon/proxy strata, component provenance and missingness families remain
in the state tables. Reconstructed origin fields are inferred from earlier
observations, not directly observed2020 state. The complete primary factorial
shows that small average state changes can accompany large conditional capacity
differences and changed Health answers. Higher-precision localization and
the full missingness propagation remain distinct downstream analyses.

Evidence: record/cohort inventories, dynamic-programming and Bernoulli tests,
initial and controlled origin/stratum comparisons under `exposure_model/`,
plus the independent analytic audit. No aggregate average establishes uniform
prospective robustness.
