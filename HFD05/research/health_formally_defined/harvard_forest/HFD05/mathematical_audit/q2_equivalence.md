# Exact reserve-qualified continuation requirement

Let `(Hist, Σ, μ)` be the declared unconditional nonnegative finite history measure, with `μ(Hist) ≤ 1`. Let measurable `V` denote lawful continuing histories and measurable `Rρ` denote the declared final reserve requirement. Reserve-qualified support is **defined** as `V ∩ Rρ`, not a separately normalized successful-history law.

For every history, membership in `V ∩ Rρ` implies membership in `V`. The disjoint partition

\[
V=(V\cap R_\rho)\;\dot\cup\;(V\setminus R_\rho)
\]

and finite additivity give `μ(V)=μ(V∩Rρ)+μ(V\Rρ)`. Nonnegativity therefore establishes

\[
0\le m_R(\rho)\le m_V\le \mu(\mathrm{Hist})\le1.
\]

For any inclusive real threshold `θ`, `(mV≥θ ∧ mR≥θ) → mR≥θ` by conjunction elimination. Conversely, `mR≥θ` and the derived inequality `mR≤mV` imply `mV≥θ`. Hence

\[
Q_2(\theta,\rho)\iff m_R(\rho)\ge\theta,
\qquad \delta_2=m_R(\rho)-\theta.
\]

The Health definition is unchanged: `Realizes(x) ∧ Adequate(C(x;d,t))`. Q2 is one reserve-qualified continuation constraint; it is not two independent adequacy constraints.

For a finite candidate bank, define `μ(A)=Σ_i w_i 1_A(γ_i)` with `w_i≥0` and total candidate weight at most one. The same pointwise subset relation and termwise order establish the inequality and equivalence for arbitrary nonnegative weights; no iid assumption is needed. Equal-weight integer census gives `0≤nr≤nv≤K`; comparison against rational inclusive thresholds can be performed exactly using integer cross-products. Zero total weight is permitted for the measure theorem; an MC interval nevertheless needs a positive sample size.

The interval correction is statistically distinct from query equivalence. A single pointwise Wilson 95% interval for reserve-qualified successes replaces the old pair of 97.5% intervals. Interior kernel certification uses strict `L>θ` / `U<θ`; remaining cases are unresolved. At threshold zero, adequacy follows algebraically while Health still depends on current realization. Finite all-success at threshold one does not establish stochastic probability one.

Verification combines this general measure proof with exact weighted finite truth tables and exhaustive/rescored archived candidate counts. The raw HFD04 bank support and count/prefix identities are checked directly. Historical source and outputs remain unchanged.
