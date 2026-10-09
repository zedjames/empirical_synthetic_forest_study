# Methods

## Health, specification and the measured organization

We study Health as current organizational realization together with adequate
prospective capacity under a declared query. A specification fixes the
organization, realization relation, context, horizon, lawful response histories,
continuation condition, candidate-history measure and capacity requirement.
Health is not equivalent to a current vital sign or a successful sampled future.
The general object is query-relative: different questions may inspect different
aspects of the same viable response measure, including the support of response
costs. A scalar continuation probability is one query of that object, not the
entire lawful definition.

Writing Ξ for the specification, x for state, d for context and h for horizon,
the general interface is

$$\operatorname{Health}_{\Xi,q}(x,d,h)
=\operatorname{Realizes}_{\Xi}(x,d)\;\land\;
\operatorname{Adequate}_{q}(\Gamma_{\Xi}(x,d,h)).$$

The viable-cost measure is the pushforward of the original candidate measure
restricted to lawful continuing histories:

$$\Gamma_{\Xi}=(\operatorname{cost}_{\Xi})_*
\bigl(\mu_{\Xi}|_{\operatorname{Lawful}_{\Xi}\cap
\operatorname{Continues}_{\Xi}}\bigr).$$

The present prototype uses terminal reserve as its response annotation for this
pushforward. It does not identify a comprehensive energetic/ecological cost
model. Its scalar queries inspect total mass and reserve-qualified mass of that
restricted measure.

The forest organization in this measurement study is a finite cohort state:
nonnegative integer stem counts, diameters and ancestral taxon/stratum/size
identities. Basal area and juvenile abundance represent two explicitly selected
organizational functions. Lawful means admissible under this model: finite
counts and sizes, unique cohort identities, diameters within 1–300 cm, stem
balance, nondecreasing cumulative deaths and births, and juvenile count no
greater than total count. Model admissibility is distinct from empirical
validation of ecological lawfulness. We do not claim that these selected
variables exhaust the forest's organization or all possible Health queries.

Let BA₀ and J₀ denote the earlier observed basal-area and juvenile anchors.
Realization requires BA/BA₀ ≥ α and J/J₀ ≥ β. Juveniles are stems with
1 ≤ diameter < 10 cm. Continuation requires realization at the terminal horizon
and no run of annual nonrealization longer than τ. A viable history is lawful
and continuing. Restricting the unconditional candidate measure to viable
histories produces a subprobability measure: failed or unlawful histories retain
their probability in its complement. Collection never filters for success and
never renormalizes the survivors. Costs in this prototype are the declared
terminal reserve functional, not a validated comprehensive ecological budget.

The continuation query asks whether viable mass meets an inclusive threshold θ.
The reserve query also requires terminal juvenile abundance divided by that
history's initial juvenile abundance to meet ρ; a zero initial denominator has
reserve zero. Its capacity is the unconditional mass of viable histories
meeting that reserve condition. This event is a subset of viability. Therefore
requiring both its mass and viable mass to meet the same θ is exactly equivalent
to requiring the reserve-qualified mass alone. This equivalence holds for every
nonnegative unrenormalized finite measure, not just equal-weight samples.
Health additionally requires present realization, independently of adequacy.

With V the viable event and R the reserve event,

$$\Gamma(C)=\sum_i w_i\,1_{V_i}\,1_{\operatorname{cost}_i\in C},\qquad
1=\mu(\neg\operatorname{Lawful})+
\mu(\operatorname{Lawful}\cap\neg\operatorname{Continues})+\mu(V),$$

$$\mu(V\cap R)\leq\mu(V),\qquad
Q_2=[\mu(V)\geq\theta\land\mu(V\cap R)\geq\theta]
\iff[\mu(V\cap R)\geq\theta].$$

## Observational support and identity

The training observations are the version-six Harvard Forest large-plot tree
census. File labels designate campaigns, not simultaneous visits. The earlier
dry census occurred in 2010–2011 with a central winter component in 2012–2014;
later dated observations occurred in 2018–2019 without the same central-swamp
coverage. The modeled origin is 2020-01-03. It is reconstructed, never directly
observed. Exact publisher files, versions, access dates, licenses and hashes are
in the input inventory. [Publisher dataset](https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF253).

Transition fitting requires stable stem identity, stable plant association and
unchanged recorded taxon. Unknown, absent, gone or ambiguously associated stems
are not certified deaths. Dated living records with missing diameter retain
alive evidence while size remains imputed. Thirty-one unresolved associations
select one branch, never duplicate a stem. Coordinates are carried or recovered
from actual recorded predecessors only. No missing individual wet locations are
generated. A central-third winter-census indicator is a protocol/coverage proxy,
not an independently certified ecological wet mask.

The fitted cells combine four taxon groups, four spatial/protocol strata and
four earlier-diameter classes, with edges 3, 10 and 30 cm assigned to their
upper bins. Taxon groups are hemlock, the registered canopy-code list, the
registered shrub-code list and other. Growth fitting requires measurement
height agreement within 0.05 m and annual diameter change between −0.1 and
1 cm/year. Original measurements are retained even when ineligible for this
specific fit.

## Actual likelihood and parameter laws

Each cell's annual hazard is drawn from a discrete 256-point grid between
−log(0.9995) and −log(0.5). The earlier-fit prior is a Beta law for annual
survival with strength 32, centered using pooled observed survival and pooled
actual exposure. Under the hazard transformation its log density is
−a h + (b−1) log(1−exp(−h)); the survival-to-hazard Jacobian is included.
The cell likelihood is −h times the sum of living exposures plus the sum over
secure deaths of log(1−exp(−h Δᵢ)). Thus this likelihood already uses individual
durations. Its pooled prior center is a separately identified approximation.
The likelihood is tempered by min(1, cap/n_fit), then normalized over the grid.
Caps 50, 100, 200 and 500 are reported; 200 is a declared reference, not an
estimated ecological effective sample size.

The growth mean is (sum of eligible changes + 32 times pooled mean)/(n_growth+32),
clipped to 0–0.7 cm/year. Growth SD is at least 0.02; sparse cells borrow pooled
SD. Parameter SE is SD/sqrt(min(max(n_growth,1),cap)). A parameter growth mean
is drawn from the corresponding Normal distribution and clipped to 0–0.7.
Recruitment intensity uses later-only living measured candidates divided by a
six-year declared exposure. Its draw is Gamma with shape min(6 mean+1,cap) and
scale (mean+1/6)/shape. Zero observed candidates therefore do not imply a
structural zero. Candidate appearance is not identified biological recruitment.

The actual parameter laws, including their supports, can be summarized as

$$H_c\sim\operatorname{Categorical}(h_j,p_{cj}),\quad
G_c=\operatorname{clip}_{[0,0.7]}\mathcal N(m_c,\operatorname{SE}_c^2),\quad
\Lambda_c\sim\operatorname{Gamma}\left(a_c,
\frac{\bar\lambda_c+1/6}{a_c}\right),\quad
a_c=\min(6\bar\lambda_c+1,\operatorname{cap}).$$

## Reconstruction with individual exposure

The original reconstruction projects each component by a Binomial law using
its cell-average duration and a component RMS diameter. The revised
reconstruction retains every actual recorded duration. Conditional on the same
sampled cell hazard and shared missingness effects, each recorded stem survives
independently with probability exp(−h Δᵢ). The sum is Poisson-binomial, with mean
Σpᵢ and conditional variance Σpᵢ(1−pᵢ). Convexity gives
Σexp(−h Δᵢ) ≥ n exp(−h mean Δ). It does not imply a universal variance ordering
or independence after shared hazards are marginalized.

Measured living stems project from their later dates; unresolved-fate stems
project from their recorded earlier dates; living missing-size stems retain
their own later-date exposure. No missing dates are invented. Unknown-fate
hazards include shared Normal log shifts with SD 0.45 and a wet-proxy shift with
SD 0.7. Other declared missingness families modify these shared mechanisms,
not the raw observations. Sizes use component RMS anchors plus sampled growth
times individual duration and the original shared component size-noise draw,
clipped to 1–300 cm. Missing sizes remain model-imputed.

Unobserved wet entry uses Poisson intensity sampled recruitment × six years ×
the wet proxy × exp(Normal(0,0.35)), initially at 1.7 cm. Association branches
remain the original branch draws. The controlled comparison captures the
original nuisance draws, unobserved-entry draw and association residual, and
changes only dated-component survivor/time calculations. Component count and
squared-diameter sums are recombined as an RMS state. This preserves the
declared aggregate basal area but not individual juvenile membership. A first
variable-consumption comparison is retained separately; its differences are not
attributed entirely to exposure because unrelated RNG consumption also changed.

## Annual growth, entry and balance

Annual survival is Binomial with probability exp(−scenario hazard factor ×
ancestral-cell hazard − hemlock excess). Surviving cohorts grow and retain their
ancestral hazard cell. Births are Poisson from the registered cell-group entry
intensity and begin at 1.7 cm. Every history records live count, basal area,
juveniles, hemlock basal area, cumulative deaths and cumulative births; annually
N_next = N − deaths + births.

For positive growth mean μ and SD s, the original Gamma shape is
max((μ/max(s,10⁻⁹))²,10⁻⁶), with scale s²/max(μ,10⁻⁹).
Consequently activating the shape floor can change its mean away from μ.
The revised Gamma uses shape a=max((μ/s)²,10⁻⁶) and scale μ/a. Its mean is μ;
its variance is μ²/a, equal to s² only when the shape floor is inactive. At
μ=0 it is deterministically zero; at s=0 it is deterministically μ. The adapter
still consumes its declared Gamma draw before deterministic replacement,
making the coupling explicit. All final diameters remain clipped to 1–300 cm.
This is a mean-preserving nonnegative law, not a claim to preserve arbitrary
unattainable target moments at the floor. [Gamma parameterization](https://numpy.org/doc/2.0/reference/random/generated/numpy.random.Generator.gamma.html).

The alternative growth law is **rectified**, not conditionally truncated,
Normal: max(0, μ+sZ). For s>0 its actual mean is
μ Φ(μ/s)+s φ(μ/s), and its second moment is
(μ²+s²) Φ(μ/s)+μs φ(μ/s). For s=0 it is μ. These are not generally the same
moments as the Gamma law. Diameter clipping changes either law's moments again.
We test pre-clipping law moments separately from post-clipping diagnostics.
Rare, extremely skewed Gamma branches may be analytically identified while
their finite-sample empirical moments remain sampling-unresolved.

For the revised positive Gamma branch,

$$a=\max((\mu/s)^2,10^{-6}),\qquad b=\mu/a,\qquad
\mathbb E X=\mu,\qquad\operatorname{Var}(X)=\mu^2/a.$$

The four demonstration scenarios have hazard/growth/entry multipliers
(1,1,1), (1.6,0.65,0.65), (1,0.85,0.9), and (1.8,0.5,0.5).
The latter two add annual hemlock hazard 0.08. These are declared scenarios,
not inferred ecological probabilities or evidence-selected forecasts.

## Numerical inference and paired correction designs

Conditional individual reconstruction and the annual count generator are

$$A_i\mid H,M,W\sim\operatorname{Bernoulli}(e^{-H_i\Delta_i}),\quad
N_c=\sum_{i\in c}A_i+U_c+B_c^{\mathrm{association}},\quad
d_c^2=\frac{\sum_{i\in c}A_i d_i^2+
U_c(1.7)^2+S_c^{\mathrm{association}}}{N_c},$$

with the earlier RMS anchor used when N_c=0, and

$$N_{t+1}=N_t-D_t+B_t,\qquad
N_{t,c}-D_{t,c}\sim\operatorname{Binomial}
(N_{t,c},e^{-H_c f_H-1_{c\in\mathrm{hemlock}}H_{\rm excess}}),\quad
B_{t,g}\sim\operatorname{Poisson}(f_B\Lambda_g).$$

For a fixed bank, exact integer counts and inclusive thresholds determine
finite-bank Health. Kernel inference instead uses one approximate pointwise
95% Wilson interval for viable count or reserve-qualified viable count.
Present Health is certified only if its lower bound is strictly above θ;
absence or an upper bound strictly below θ gives false, otherwise it remains
Monte Carlo unresolved. At θ=0 adequacy is algebraic; present realization is
still required. At θ=1 an all-success bank establishes its finite answer but
never certifies stochastic probability one. Earlier emitted confidence
statuses are preserved separately from corrected retrospective count inference.

The factorial correction separates original calculation, single-query
inference, mean-preserving Gamma, individual exposure, and both model
corrections. Complete primary surfaces use 256 reconstructed states × two
parameter replicates × two growth models, with 64 futures per scenario and
the original full semantic grid. Five missingness families use 128 states ×
two models × 256 futures per scenario. Every conditional unit retains equal
declared weight; the reported fraction of healthy units is not Health of their
average capacity and is not an ecological posterior. No prior over scenarios
or semantic thresholds is supplied. Maxima, upper quantiles and conditional
boundary effects accompany ensemble means.

Adverse primary cases are selected from the largest observed local capacity
change before generating higher-precision futures. Every available attaining
physical case, up to the planned quota, is retained. Fixed prefixes of 256,
1,024 and 4,096 new futures diagnose whether the extreme is only coarse
sampling noise. Marginal difference brackets may combine two 97.5% intervals
for **two distinct model laws**; these are not the abolished two redundant
constraints of the reserve query. Targeted extremes do not estimate global
robustness.

## Synthetic decision reliability, oracle interfaces and context

All 81 boundary-designed worlds and all 324 earlier factorial worlds remain
in the analysis. World design uses deterministic moment approximations, never
realized future outcomes, to approach selected decision boundaries. Original
reference banks remain fixed at 65,536 paths with 4,096 and 16,384 prefixes.
Corrected references use 4,096/16,384 for boundary worlds and 4,096 for the
factorial panel; corrected 65,536-path results are not claimed. Every original
and corrected prediction/reference combination is separately scored. Synthetic
strict packets contain no individual observation dates, so the exposure
correction is explicitly inapplicable there.

The estimator receives only thinned counts, noisy observed diameters and
policy-permitted context. It does not receive true state or parameters. Full
pressure labels anchor hazards at 0.005, 0.03 or 0.12; recovery labels anchor
total entry at 0.025, 0.008 or 0.001 times baseline juveniles. Lognormal hazard
and entry variation has log SD 0.15 and 0.2, centered to preserve means. Unknown
context mixes all nine contexts uniformly; coarse context mixes only within
its disclosed groups; label-only discards state observations. These are
specified design priors. The boundary roster's nonuniform context frequencies
and a separately registered roster-prior experiment remain explicitly distinct.
They are not different information levels of one optimal Bayes procedure.

The current comparison uses 32 estimator draws × 256 futures. State count
inference uses the registered Gamma–Poisson anchor model, and diameter
inference combines prior SD 1 cm with observation SD 0.25 cm. The estimated
Health fraction is classified at 0.5. Its 5/95% between-draw mass range mixes
state, parameter and future variability; it is not an inner-Monte-Carlo-only
confidence interval.

Four strictly typed oracle conditions distinguish estimated state/parameters,
true state only, true parameters only, and both. The observation-only interface
rejects oracle payloads. Fixed future prefixes and replicate panels vary
observation, inference or future RNG sources one at a time. Signed state and
parameter gains and their factorial interaction are reported; they are not
forced into an additive causal error budget. Supplying both oracles still
retains future Monte Carlo and predictive-law misspecification. Extra hemlock
hazard 0.12 and annual lognormal entry dispersion 1.2 remain separate suites.

Classification rates use only explicitly certified reference labels; every
unresolved reference remains in the tables and in sensitivity brackets. Paired
world bootstraps respect constructed-world strata and use 1,000 replicates.
Empty strata remain undefined, not zero. These intervals are conditional on
the declared design, not ecological-population confidence statements.

For oracle losses L₀₀,L₁₀,L₀₁,L₁₁, the reported signed interaction is

$$I=L_{11}-L_{10}-L_{01}+L_{00}.$$

Neither I=0 nor nonnegative stage gains are imposed. For each capacity query,
the signed margin and Monte Carlo decision relation are

$$m=\widehat\Gamma(C_q)-\theta,\qquad
\text{TRUE if }L_q>\theta;\quad
\text{FALSE if absent or }U_q<\theta;\quad
\text{otherwise unresolved},$$

with the separate algebraic theta-zero rule already specified. Between-draw
mixture ranges and paired bootstrap contrasts have the distinct scopes above.

## Representation, entry identification and latent tipping

The size-enriched representation splits each cohort into two visible-data
conditional sizes and rescales its second moment. Missing groups use an
explicit imputed fallback. Initial count and basal-area preservation do not
imply juvenile preservation or commuting stochastic transitions. All original
selected states, scenarios and cap values remain in their corrected future
comparisons. Individual-date reconstruction of fine latent size profiles is
not asserted where that bridge has not been implemented.

Query insufficiency is witnessed even by two-stem states: diameters (6,12) and
(√90,√90) have the same count and basal-area second moment, but juvenile counts
one and two. Thus, for structural projection π and juvenile query J,

$$\pi(x)=\pi(y)\quad\text{while}\quad J(x)\ne J(y),$$

so J cannot factor through that projection. Preserving π at one time supplies
no transition-intertwining equation for the stochastic kernels.

The 6,992 later-only living measured candidates retain their original identity,
plant, taxon, size, date and location audit. A new recorded stem may reflect
threshold crossing, an earlier missed stem, a new stem on an existing plant,
tag/identity change or frame change. A later record with no earlier stem record
is compatible with both threshold crossing and an already-eligible missed stem.
The observation map is therefore noninjective for biological entry. Diagnostic
membership weights are not posteriors. A fixed small diagnostic panel varies
unknown-fate log hazard by −0.5/0/+0.5 and future candidate-entry intensity by
0/0.5/1, retaining all combinations. Its tipping ranges complement, but do not
replace, the complete five-family semantic surfaces and measured/imputed,
taxon and coverage-proxy provenance.

## Independent external components

The later adult census uses the fixed version-one 2021–2024 public dataset.
Exact recorded stem tags link to stable earlier identities; ambiguous matches,
association conflicts, nonadult earlier sizes, baseline absence/death and secure
dead-to-alive reversals are excluded with a full ledger. Of 30,808 tags, 21,578
form the eligible 2021-survivor adult cohort. Cumulative and annual windows
retain distinct risk denominators. Alive codes A/AU and secure death codes
DC/DS define primary endpoints. DN (“dead and not found”) is not a secure
biological death in the primary analysis; all-DN-death and all-DN-alive bounds
are retained. No implicit forward filling supplies an endpoint. Year-only July
precision is an anniversary approximation with ±30-day prediction diagnostics,
not an invented exact date. [Adult census](https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF453).

Frozen earlier-fit hazard predictions are integrated over the existing hazard
posterior with no external refit or outcome-selected cap. We report Brier and
log scores, calibration bins, mortality residuals and hemlock, taxon, size,
damage and spatial-proxy strata. Probability clipping at 10⁻¹⁵ is solely for
log arithmetic, not a changed prediction. Whole-quadrat and whole-stem-history
bootstraps provide conditional descriptive intervals. A separate scalar-checked
aggregation replay tests platform matrix-product warnings; original intervals
are retained and compared, never silently replaced.

The fixed version-four seedling dataset observes below-1-cm stems in 1-m²
subplots. Its 7,193 eligible identities support descriptive dated survival
transitions, not the forest model's whole-plot 1–10-cm juvenile target. Repeated
visit fields determine sampled alive/dead status; secure dead-to-alive histories
are excluded. Actual visit dates determine exposure, including genuine campaign
gaps. First observation is not germination. The graduation-coded identity was
already recorded graduated at its initial visit, so its transition date is
left-censored. Recorded canopy photos supply only their recorded subplot/year
context. No unsupported spatial expansion, missing-year extrapolation or
seedling hazard bridge is fitted. [Seedling dataset](https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF355).

These independent components test identifiable parts of the measurement model.
They are retrospective tests under fixed support, not a publicly issued
prospective forecast and not validation of the full Health proposition.

## Reproducibility and interpretation

Source/configuration hashes are frozen before each newly executed phase.
Scientific predecessors and raw inputs remain byte-identical; new banks have
separate identities. Python 3.9.6 and NumPy 2.0.2 are pinned. Independent count,
moment, raw-source and score implementations test the active calculations.
Replay separates lightweight mathematical/unit fixtures, published table and
headline reconstruction, and full ensembles/external retrieval. Disclosure and
source-license review are distinct from numerical reproducibility. Public
publication is not implied by preparing a local review artifact.

The conclusions answer three distinct questions: whether the measurement model
is empirically identified and calibrated; whether specified decisions are
reliable near their boundaries; and whether a representation preserves a
particular query. Successful execution of a declared specification answers none
of those empirical questions by itself. Exact deterministic set refinement
under fixed semantics is separate from changes to stochastic candidate laws.
