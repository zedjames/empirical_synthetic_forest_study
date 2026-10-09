# Monte Carlo resolution

Fixed maximum banks first, then paired K64/256/1024 prefixes on16 B0 states,2 parameter draws,2 model variants,D0/D3,h5/10/20. Separate outer N128/256/512/1024 prefixes at K64 give4096 maximum outer units. Inner and outer studies are NOT a universal N1024×K1024 run.

Selected study: ADEQUATE; maximum-K resolved fraction1.0, maximum interval half-width0.031999, maximum outer state-block Health SE0.004651. Targets were predeclared; no4096 extension or preferred-Boolean stopping.

Wilson score95% Q1 and two-component Bonferroni97.5% Q2 intervals have approximate pointwise binomial coverage, not a simultaneous guarantee across the entire grid. Initial realization is conditional on each sampled state. TRUE requires all lower bounds strictly above theta; FALSE requires absence or an upper bound strictly below; otherwise MC_UNRESOLVED. Strict inequalities mean finite samples cannot confidently equal theta0/1 endpoints.

The broader missingness/semantic surface retains963976 unresolved conditional cells, so overall publication precision is MIXED, even though the selected maximum-K panel meets its targets. Read [inner](../monte_carlo/inner_resolution.csv), [outer](../monte_carlo/outer_resolution.csv), and portable sufficient counts. Outer state blocks preserve dependence among parameter/model draws. MC intervals do not cover reconstruction-model or ecological uncertainty.

Full-grid resolution by threshold:

| theta | conditional_units | MC_UNRESOLVED |
| --- | --- | --- |
| 0 | 1244160 | 631746 |
| 0.5 | 1244160 | 1377 |
| 0.75 | 1244160 | 1017 |
| 0.9 | 1244160 | 1590 |
| 1 | 1244160 | 328246 |

Most unresolved labels occur at theta0/1 under the predeclared strict interval rule. At theta0, nonnegative capacity already establishes adequacy algebraically: its conservative MC label is NOT evidence of uncertainty about that algebraic fact. Interior thresholds retain3984 unresolved conditional cases at K256, so the MIXED certificate is not based solely on endpoint artifacts. No post-outcome rule change is used to make the surface appear more resolved.

Method source: [NIST Wilson interval discussion](https://www.itl.nist.gov/div898/handbook/prc/section2/prc241.htm).
