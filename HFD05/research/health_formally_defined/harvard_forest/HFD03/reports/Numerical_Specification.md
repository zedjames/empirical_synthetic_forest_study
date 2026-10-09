# Reader-inspectable numerical specification

| id | hazard_factor | growth_factor | recruitment_factor | hemlock_extra_hazard |
| --- | --- | --- | --- | --- |
| D0 | 1 | 1 | 1 | 0 |
| D1 | 1.6 | 0.65 | 0.65 | 0 |
| D2 | 1 | 0.85 | 0.9 | 0.08 |
| D3 | 1.8 | 0.5 | 0.5 | 0.08 |

Frozen MNAR log-hazard SD.45, wet-proxy SD.7, unseen entry log SD.35, measurement SD.1cm, synthetic prior DBH SD1cm and observation SD.25cm. Association branches: equiprobable recorded E0/E1 identity alternatives, not confirmed reassociations. Fit:256-point annual hazard grid, survival-prior strength32, effective cap200, six-year candidate-entry exposure, growth inclusion[-.1,1]cm/y, growth means[0,.7], POM tolerance.05m. Origin2020-01-03. R1 uses BA/E0BA and juvenile1–10cm count/E0juveniles; C1 uses annual path excursions plus endpoint realization, not endpoint alone.

Alpha[.5,.75,1], beta[.1,.4,.8], tau[0,2,5]years, theta[0,.5,.75,.9,1], reserve rho[.5,1]. Model-admissibility checks identity/count balance/nonnegative bounded size, not scientific ecological Lawful.

| study | states | parameters | models | K | outer_units | history_weight |
| --- | --- | --- | --- | --- | --- | --- |
| HFD02S frozen | 256 | 2 | 2 | 64 | 1024 | 1/(1024*64) |
| HFD03 missingness, per family/scenario | 128 | 1 | 2 | 256 | 256 | 1/(256*256) |
| HFD03 MC inner panel, per scenario | 16 | 2 | 2 | 1024 | 64 | 1/(64*1024) |
| HFD03 MC outer maximum, per scenario | 1024 | 2 | 2 | 64 | 4096 | 1/(4096*64) |

Each outer unit has equal design weight1/outer_units; each conditional history1/K; their product gives the listed joint weight. Successful histories retain original weights: failures are never dropped or renormalized. Scenario and semantic cells remain separately indexed, not given posterior weights. For each current state Health requires present realization AND capacity adequacy; averaging conditional Health yields a finite-design fraction, not Health of an averaged organization.

Counts/K is the exact finite-atom capacity of the archived bank, with a definite finite-bank Health value. Wilson intervals address sampling resolution as an approximation to the separately declared stochastic demographic kernel; they are not uncertainty in finite-bank arithmetic or in the formal implication itself. Both meanings are kept distinct, and neither is an ecologically established probability law.

Every configuration leaf is exposed in [numerical assumptions](../specification/numerical_assumptions.csv), with both frozen HFD02S and HFD03 protocols. [Editorial map](../specification/editorial_mapping.csv) fixes V2 terminology. Sources: [HF253 publisher archive](https://harvardforest1.fas.harvard.edu/exist/apps/datasets/showData.html?id=HF253) and the read-only frozen EML/news receipts.
