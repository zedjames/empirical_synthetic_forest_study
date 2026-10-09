# Actual observation times and conditional survivor laws

Condition on a cohort hazard h≥0 and on all shared latent shifts. For individual
dated exposures Δᵢ≥0, pᵢ=exp(-hΔᵢ). Conditional independence gives the
Poisson-binomial survivor count S=Σ Bᵢ, with

    E[S | h,Δ]=Σ pᵢ; Var[S | h,Δ]=Σ pᵢ(1-pᵢ).

The original pooled approximation instead uses Binomial(n,exp(-h mean(Δ))).
Since the second derivative of exp(-hΔ) is h² exp(-hΔ)≥0, Jensen's inequality
gives mean(pᵢ)≥exp(-h mean(Δ)). Hence the individual conditional survivor mean
is at least the pooled mean. Equality holds for equal exposures, h=0, or an
empty group. This is a conditional law comparison, not a claim that any finite
sample of the individual law must have more survivors. The two variances are
reported separately; no universal variance ordering is asserted. Mixing over
shared hazards/shifts introduces dependence and is not an unconditional
independent Poisson-binomial law.

The empirical hazard likelihood already uses each training duration: living
records contribute -hΣΔ, and dead records contribute Σlog(1-exp(-hΔ)). Its
pooled approximate prior center is a separate operation. This tranche retains
the fitted model and inventories that approximation instead of implying the
entire fit has become exact.

Three reconstruction components require actual times: measured living
endpoint records project from their dated endpoint to origin; unresolved
endpoint fates retain the actual dated baseline exposure; dated living records
without DBH project survival while retaining latent size. Unknown or absent
fate never becomes recorded death. Existing unresolved association branches
already use individual dates and retain their original law. All record counts
and exposure sums must reproduce the original sufficient statistics.

The individual-time implementation retains the 64-cell count/RMS carrier.
Component baseline RMS size and shared measurement/prior-size noise remain
unchanged, but growth uses each record's exposure; squared projected size is
summed only for sampled survivors, then reaggregated. These are modeled
component sizes, not invented observed individual DBHs or wet-area positions.
Changing exposure can change the final RMS juvenile threshold classification;
this is explicitly reported and is not dynamic commutation.

Two reconstruction comparisons are retained. The initial variable-consumption
comparison changes unrelated future RNG consumption. The controlled comparison
captures original parameters, shifts, size-noise, unseen entry and association
draws, and substitutes only the three dated component laws. Its individual
Bernoulli draws still have MC error. Analytic Jensen gaps, finite differences,
local maxima and 95th percentiles remain separate evidence.
