# Paper V V3 numerical Methods contract

## Scientific object and measured support

HF253 v6 (Orwig, Foster, Ellison2023, DOI10.6073/pasta/818789a882a318c1d7f3fc43a2289e12) supplies E0 stems2014 and E1 stems2019; file labels do not imply simultaneous dates. E0 dry census2010–2011, central winter component2012–2014; E1 dated2018–2019 excludes the central swamp. CC0 dataset and citation: https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF253 . Frozen inputs are byte-hashed; canonical inputs manifest specifies exact filenames/URLs. No E2 enters training.

Identity uses stem.id + consistent tree.id; raw taxon must agree for transition fitting. Absent, gone, missing and ambiguous associations do not become deaths. Secure calibration death additionally requires D and excludes Miss/MT/S. Stem-level death is not whole-plant death. DBH missing from dated alive records retains survival evidence, not an invented measurement. Stable tagged stationary-root coordinates are recovered only from recorded E0 fields. No missing wet individual positions are generated. Wet-candidate/protocol-proxy is a central-third gx winter2012–2014 flag, not a certified ecological wet mask.

Four taxon groups: tsugca; enumerated canopycodes; enumerated shrubcodes; other. Exact code lists are empirical.taxon. Stratum is floor(gx/(700/3)), clipped0..2, plus fixed E0 winter central proxy3. Size uses E0 DBH edges3,10,30 with right-side insertion (boundary belongs to upper bin); fitted cell=(taxon×4+stratum)×4+size. Source hom field checks same measurement height within0.05m. Only growth between−0.1 and1cm/year enters fit; original recorded values are not overwritten. Coordinate and DBH eligibility remain separate.

## Demographic likelihood and parameters

Annual hazard grid:256 equally spaced points between−log(.9995) and−log(.5). Let pooled h=−log(global observed survival)/mean(actual exposure). Beta annual-survival prior has strength32: a=32exp(−pooled h), b=32−a. Its hazard-space log density is−a h+(b−1)log(1−exp(−h)), including the survival-to-hazard Jacobian. Per-cell log likelihood is−hΣ_alive duration+Σ_dead log(1−exp(−h duration)). Multiply likelihood by min(1,n_cap/n_fit) with caps50,100,200,500;200 is frozen reference, not an estimated effective sample size. Normalize discrete posterior over the256 gridpoints.

Cell growth mean=(sum growth+32 pooledmean)/(n_growth+32), clipped0..0.7. SD=max(.02, sample populationSD), sparse cells borrow pooledSD. SE=SD/sqrt(min(max(n_growth,1),n_cap)). Parameter draws use the grid posterior and normal(mean,SE) clipped0..0.7. Candidate entry mean is E1-only living measured count/6years, not biologically identified recruitment. Wet proxy borrows dry species/size share. Entry draw is Gamma(shape=min(6mean+1,n_cap),scale=(mean+1/6)/shape).

## Origin reconstruction

Origin is2020-01-03, not the E1 observation dates. Observed living measured cohorts are projected by Binomial(n1,exp(−h tail)) and RMS_d1+g tail+Normal(0,.1). Unknown-fate cohorts use Binomial(n_unknown,exp(−h exp(M+wetW) exposure)) with M~Normal(0,.45),W~Normal(0,.7), and RMS_d0+g exposure+Normal(0,1). Dated living missing-size cohorts are separately projected and size-imputed. Unseen wet count is Poisson(entry×6×wet×exp(Normal(0,.35))) at1.7cm. Thirty-one unresolved associations choose one branch, never duplicate a stem. Cohort diameter is sqrt(sum component n d²/sum n), preserving declared component BA, not observed individual juvenile membership. Every origin-state field is inferred; observation ancestry does not make a2020 value observed.

## Annual predictive kernel and continuation

N_next=N−D+B. D is binomial mortality with survival exp(−hazard×scenariofactor−hemlock excess). Growth is either gamma with shape=max((mean/SD)²,1e−6),scale=SD²/mean, or Normal(mean,SD) truncated atzero; both have DBH clipped1..300. Entry is Poisson of cohort-group entry×scenariofactor. Each annual entry cohort starts at1.7cm and retains its ancestral cell identity. Original model has no size-dependent mortality beyond original cell indexing.

| scenario | hazard | growth | entry | hemlock excess |
| --- | ---: | ---: | ---: | ---: |
| D0 |1|1|1|0|
| D1 |1.6|.65|.65|0|
| D2 |1|.85|.9|.08|
| D3 |1.8|.5|.5|.08|

These are synthetic demonstration assumptions, not identified ecological futures. Annual path records retain living count, BA, juvenile count, hemlock BA, cumulative deaths, cumulative entries. All candidate histories retain denominator weight, including failures.

## Formal measurement interface

X is the declared finite cohort state (counts/diameters/ancestral cells); d is the fixed scenario and horizon-relevant context; t is annual modeled horizon. Hist is the unconditional candidate bank. Lawful is a model-admissibility realization of Lawful: finite nonnegative counts/diameters, unique cohort identities, bounds1..300, count balance, monotone deaths/entries and juvenile<=total. This is not empirical validation or proof of ecological lawfulness.

Realizes is BA/BA0>=alpha AND J/J0>=beta, where J counts1<=DBH<10. Continues requires final realization and longest contiguous annual excursion below realization<=tau. o is either tagged field observation/missingness or explicitly simulated thinning+diameter noise. collect is exhaustive candidate-history generation, not success filtering. Capacity is unconditional lawful-continuing mass; Q2 adds unconditional viable-with-final-J/own-initial-J>=rho mass. r is the declared threshold requirement; Adequate uses inclusive>=theta. Health requires current realization AND adequacy; absence is independently FALSE.

## Boundary and context designs

Two independently streamed pilots are excluded from confirmatory performance; pilot001 failed Q2 count-scaling and is retained. Pilot002 maps structure and juvenile diameter to expected boundary locations. The final81-world roster is frozen without final outcome selection. Analytic parameter-moment normalization is a world-design approximation, not a reference oracle; latent parameters remain hidden from estimators. Initial states may be correlated with context/parameters by experimental construction. Correct-model means shared response-kernel class, NOT that the observation-prior joint distribution is the true world-design distribution.

Truth uses fixed65536-path banks, with paired chunk prefixes4096/16384/65536; reference-mass uncertainty is retained. Correct, extra hemlock hazard.12, and annual lognormal entrysd1.2 suites stay separate. Q1 signed margin=mV−.75; Q2=min(mV,mR)−.75. Present-realization failures are separately identified. Reference-uncertain worlds remain in all tables, while certified-truth scores exclude them explicitly.

Pressure labels anchor hazards[.005,.03,.12]. Recovery anchors totalentry=[.025,.008,.001]×baselineJ. Cell share=(fittedentry_c+1/64)/sum(fittedentry+1/64). Multiply each hazard by exp(Normal(−.15²/2,.15)); each cellentry by exp(Normal(−.20²/2,.20)). No true latent rates enter the estimator. Full context carries exact labels; coarsened pressure groups[0,1]/[2], recovery[0]/[1,2]. Unknown context mixes allnine equally. Coarse mixture is uniform inside its disclosed group. These are design priors, not ecological scenario probabilities. Label-only replaces all observations with missing markers and uses the public E0 prior anchor.

Estimator uses32 state/parameter draws×256 futures, not historical64-draw HFD03 precision. Current B0 is the paired attribution reference. Gamma–Poisson count inference uses shape16 and rate16/public E0 count; diameters combine Normal priorSD1 with observationSD.25. Prediction p_H is fraction of finite design draws satisfying Health, NOT Health of averaged capacity. 5/95% mixture intervals describe observation/parameter/design dispersion, not just MC error. Per-draw inner MC statuses yield a separate resolved-fraction bracket; its relation to.5 classifies the finite-draw prediction interface, not a population confidence theorem.

### Separately registered context-prior reconciliation

The original 324-world factorial roster has uniform nine-context frequencies. The boundary-designed 81-world roster instead has counts [7,7,7,7,7,7,13,13,13], ordered pressure-major/recovery-minor. Original uniform-nine boundary predictions are retained as a registered prior sensitivity, not described as that roster's marginal prior. HFD04-CONTEXT-PRIOR-001 was registered after this design mismatch was identified and before any supplemental outcomes. It uses those counts divided by 81, chosen only from the frozen world roster, with independent PCG64 namespaces cw_state/cw_params/cw_future and master seed 202610070404. The actual pressure/recovery labels remain hidden; the sampled context column is a simulated prior draw, not a recovered realized label. All 81 worlds, 32 draws, 256 paths, three separate truth mechanisms, both queries, and every signed/absolute/failure stratum are retained. This is a versioned methodological correction, not a replacement or outcome-tuned prior. Independent streams mean contrasts with the original mixture also contain finite-draw sampling variation.

## Juvenile representation and effective cap

Real E1 audit is securely observed dated subset only. Methodological Gaussian DBH perturbations0/.1/.25/.5/1cm are not estimated instrument errors. Masked DBH never enters fitted size profiles. Enriched reconstructed states use two visible-E1 conditional mean sizes per cohort, below10 andatleast10; fallback is imputed. Integer split and second-moment rescaling establish initial count/BA projection to numerical tolerance. J does not commute, and transition kernels are NOT claimed to commute. Distributional ancestry is distinct from newly measured modality. Paired seeds are shared, but changed RNG array dimensions preclude exact common-uniform coupling.

Effective-cap experiments reuse exact original fit/parameter/reconstruction code in isolated module globals. They report allfour caps with16 states×two growthmodels, four scenarios and256 paths, horizons5/10. M0 heldout calibration uses64 draws percap. No cap is chosen by E2 agreement. Missingness interaction is limited to B0; no broader interaction is claimed.

## Exact finite-bank answer versus stochastic inference

Finite counts/K and inclusive thresholds give exact archived-bank Health. Wilson95% Q1, two Bonferroni97.5% Wilson components Q2 give approximate pointwise stochastic-kernel intervals. TRUE requires all lower bounds strictly above theta; FALSE follows absence or an upper bound belowtheta; otherwise MC_UNRESOLVED. At theta0 adequacy is algebraic, present still required. At theta1 all-success finite-bank TRUE never certifies kernel probabilityone. Historical HFD03 unresolved counts are untouched; derived versioned presentation changes only endpoint interpretation.

## Weights and claim scope

Frozen HFD02S has256states×2parameters×2growthmodels=1024outer units;64paths perunit/scenario: weights1/1024 and1/64, joint1/65536. HFD04 truth weights1/K withinworld/suite; benchmark world weight1/81 (and1/324 historical panel). Estimator outerdraw weight1/32 andpath weight1/256. Cap ensemble16×2=32outer units percap/scenario; representation ensemble16×2=32units perrepresentation/scenario. Context mixtures are conditional design distributions, not priors on realized ecology. No scenario or semantic prior is assigned.

Finite-design Health fraction=(1/N)sum indicator[unit Health]. This is NOT indicator[(1/N)sum unit continuation mass>=theta]. Manuscript wording: finite-design Health fraction under the declared demonstration specification. Operational Harvard Forest assessment NOT_AUTHORIZED; retrospective E2 validation NOT_CLAIMED.

All numerical leaves, executable bindings and pinned runtime are exposed in numerical_methods.csv and the source inclusion manifest. Source package is prepared for separate disclosure/license review; no public release/DOI is authorized here.
