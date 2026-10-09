# Original-roster inference-stage attribution

The original81-world roster is evaluated before any corrected-generator
comparison. The four strict interfaces distinguish estimated state/parameters,
known state only, known parameters only, and both known. The baseline retains
all32 estimator draws and256 future paths; all2,592 original baseline draws
reproduce archived counts exactly. Eight fixed draws per world also use paired
256/1024/4096 prefixes. A fixed eight-world panel varies observation, inference
and future seeds one source at a time. Three truth suites remain paired on
the same world, not counted as three independent world samples.

The independent audit checks34,176 draw/query/prefix records, every exact
baseline count, information-role identity and finite query rule. Numerically
identical oracle payloads supplied under either wrong role are rejected by
the baseline interface. The forecasting kernel remains nominal in all roles:
the hemlock hazard supplement and annual entry-dispersion mismatch are not
magically removed by knowledge of state and nominal parameters.

## Correct-model PRESENT near-boundary baseline

| Query / information | Certified worlds | FP | FN | Balanced accuracy | Target mass MAE |
| --- | ---: | ---: | ---: | ---: | ---: |
| Viable / estimated both | 31 | 10 | 4 | 0.563025 | 0.250893 |
| Viable / known state | 31 | 7 | 4 | 0.651261 | 0.240601 |
| Viable / known parameters | 31 | 8 | 7 | 0.514706 | 0.367734 |
| Viable / known both | 31 | 0 | 0 | 1.000000 | 0.003966 |
| Reserve / estimated both | 46 | 17 | 4 | 0.554924 | 0.224285 |
| Reserve / known state | 46 | 11 | 8 | 0.589015 | 0.205603 |
| Reserve / known parameters | 46 | 15 | 7 | 0.528409 | 0.317096 |
| Reserve / known both | 46 | 0 | 0 | 1.000000 | 0.003819 |

Uncertified reference worlds are retained but not treated as known truth in
these certified classifier rates. The results strongly interact: supplying
parameters alone worsens target mass error, whereas supplying both nearly
removes it in the correctly specified synthetic experiment. This is not an
additive partition of real ecological failure or a guarantee that every oracle
intervention benefits an approximate estimator.

Using absolute target-mass loss, the reserve-query state-only gain is0.01868145
with design-conditional paired-family bootstrap interval[-0.02054188,0.05758398];
parameter-only gain is-0.09281158 with interval[-0.14573570,-0.04111485]. The
interaction `L11-L10-L01+L00` is-0.29459616 with interval
[-0.34838481,-0.23601683]. Each interval conditions on the constructed world
design, not an ecological population. All signed-margin strata, reference
ambiguity and fixed-prefix numerical changes remain reported.

Even the known-both condition retains finite future MC error and, in the
misspecified suites, generator mismatch. Its correct-model classifier result
does not establish observation-conditioned operational forest Health.
An empty hemlock/reserve PRESENT-near stratum is retained with undefined
contrast estimates; no favorable roster or precision is substituted.

Evidence: `oracle_attribution/` draw census, losses, scores, interactions,
near-boundary tables, one-source contrasts, convergence and independent audit;
`uncertainty/oracle_*` contains qualified intervals. The first empty-stratum
bootstrap implementation and failure receipt are preserved separately from
the successful explicit second implementation.
