# Nonnegative growth law and the implemented floors

For a Gamma increment with shape a and scale b, E[G]=ab and Var[G]=ab².
Let μ≥0 be the requested annual growth and s≥0 its requested SD. The original
implementation sets

    a₀=max((μ/max(s,10⁻⁹))²,10⁻⁶)
    b₀=1{μ>0} s²/max(μ,10⁻⁹).

Thus its actual mean and variance are a₀b₀ and a₀b₀², not in general μ and s².
For μ≥10⁻⁹, s≥10⁻⁹ and an inactive shape floor, these equal the requested
moments algebraically. A reached positive shape floor increases the mean
above μ in that region. For 0<μ<10⁻⁹ and fixed positive s, the scale denominator
floor bounds b₀, but the positive-μ limit can still differ discontinuously from
the value zero at μ=0. The SD denominator floor is a third, separate branch.
At s=0 the original law has scale zero even when μ>0 and returns zero growth.
At μ=0 it returns zero, notwithstanding the positive shape floor.

The corrected positive law is

    a₁=max((μ/s)²,10⁻⁶), b₁=μ/a₁  for μ>0 and s>0.

Consequently E[G]=μ and Var[G]=μ²/a₁. The requested variance s² is achieved
only when the shape floor is inactive. When the floor is active, variance is
reduced to μ²/10⁻⁶; this is a declared regularization, not a claim that both
requested moments remain exact. At μ=0 the corrected law is deterministic
zero. At s=0 it is deterministic μ. A nonnegative random variable with mean
zero is almost surely zero, so positive variance in the zero-mean branch is
mathematically infeasible. These deterministic cases are not invalid
zero-shape Gamma calls.

The simulator subsequently clips diameter to [1,300]. Therefore the increment
in the state is min(G,300-d) for d in that interval, with generally lower mean
and changed variance. Pre-clipping Gamma moments and post-clipping diameter
moments are different targets. The registered clipping grid reports both.

The source-based census reproduces original parameter stream identities.
It distinguishes zero means, each denominator floor, positive shape-floor
activation, and a >1% conditional mean discrepancy. It counts parameter laws
and cohort-year multiplicities, not all live biological stems. Birth cohorts
inherit the lowest-size ancestral cell for each of the 16 taxon/sector groups.
Pilots are excluded; recruitment-parameter Gamma laws are a distinct construct.

Tiny shapes have very large relative sampling error and many numerical zero
draws. The registered simulation diagnostics use analytic Chebyshev bounds,
and mark inadequate relative precision MC_UNRESOLVED. Successful diagnostics
do not replace the analytic moment calculation or prove that rare tails have
been resolved. Both fixed and all-reached-primary-identity impact panels remain
reported; neither alone proves full-grid or ecological robustness.
