# Context-label information ablation

Paired405 worlds:all81 new worlds plusall324 frozen HFD03 worlds/observations. Every context policy uses the same observation and state/future random nonce. Known labels, coarse markers and explicit unknown are distinct strict types; unavailable estimators never receive exact levels. Coarse pressure[0,1]/[2], recovery[0]/[1,2]; uniform conditional design mixtures. Unknown mixes9 contexts equally, not ecological probabilities.

| cohort | query | policy | FP | FN | balanced_accuracy | viable_MAE | reserve_MAE | brier | viable_width90 | estimator_unresolved_fraction | delta_vs_B0_balanced_accuracy |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| boundary | Q1 | B0_full | 20 | 4 | 0.7082111436950147 | 0.1809159267095872 | 0.20118656864872686 | 0.2303515625 | 0.32124083719135804 | 0.024691358024691357 | 0.0 |
| boundary | Q1 | B1_coarsened | 14 | 4 | 0.7763929618768328 | 0.18315727328076775 | 0.20498167438271606 | 0.17450520833333333 | 0.5202907986111112 | 0.012345679012345678 | 0.06818181818181812 |
| boundary | Q1 | B2_unknown | 16 | 2 | 0.7859237536656891 | 0.1499711024908372 | 0.1743348675009645 | 0.14953125 | 0.6669825424382716 | 0.024691358024691357 | 0.07771260997067442 |
| boundary | Q1 | B3_label_only | 32 | 24 | 0.24926686217008798 | 0.5572387318552277 | 0.5198157039689429 | 0.7458723958333333 | 0.0 | 0.0 | -0.4589442815249267 |
| boundary | Q2 | B0_full | 27 | 4 | 0.6519607843137255 | 0.1809159267095872 | 0.20118656864872686 | 0.30212239583333333 | 0.32124083719135804 | 0.04938271604938271 | 0.0 |
| boundary | Q2 | B1_coarsened | 21 | 5 | 0.6899509803921569 | 0.18315727328076775 | 0.20498167438271606 | 0.24696614583333334 | 0.5202907986111112 | 0.037037037037037035 | 0.03799019607843135 |
| boundary | Q2 | B2_unknown | 23 | 2 | 0.732843137254902 | 0.1499711024908372 | 0.1743348675009645 | 0.21912760416666666 | 0.6669825424382716 | 0.024691358024691357 | 0.08088235294117652 |
| boundary | Q2 | B3_label_only | 32 | 17 | 0.33210784313725494 | 0.5572387318552277 | 0.5198157039689429 | 0.6533723958333333 | 0.0 | 0.0 | -0.31985294117647056 |
| HFD03 | Q1 | B0_full | 4 | 0 | 0.9921568627450981 | 0.011463683328510802 | 0.011463683328510802 | 0.006004050925925926 | 0.02628882137345679 | 0.0030864197530864196 | 0.0 |
| HFD03 | Q1 | B1_coarsened | 19 | 6 | 0.9192668371696504 | 0.07340118031442901 | 0.07340118031442901 | 0.043354552469135804 | 0.1575611255787037 | 0.0030864197530864196 | -0.07289002557544766 |
| HFD03 | Q1 | B2_unknown | 33 | 18 | 0.8048593350383632 | 0.19815628616898148 | 0.19815628616898148 | 0.1054175106095679 | 0.4439597800925926 | 0.027777777777777776 | -0.1872975277067349 |
| HFD03 | Q1 | B3_label_only | 147 | 0 | 0.711764705882353 | 0.4482877754870756 | 0.4482877754870756 | 0.4533420138888889 | 0.001177300347222222 | 0.0 | -0.2803921568627451 |
| HFD03 | Q2 | B0_full | 4 | 0 | 0.9921568627450981 | 0.011463683328510802 | 0.011463683328510802 | 0.006004050925925926 | 0.02628882137345679 | 0.0030864197530864196 | 0.0 |
| HFD03 | Q2 | B1_coarsened | 19 | 6 | 0.9192668371696504 | 0.07340118031442901 | 0.07340118031442901 | 0.043354552469135804 | 0.1575611255787037 | 0.0030864197530864196 | -0.07289002557544766 |
| HFD03 | Q2 | B2_unknown | 33 | 18 | 0.8048593350383632 | 0.19815628616898148 | 0.19815628616898148 | 0.1054175106095679 | 0.4439597800925926 | 0.027777777777777776 | -0.1872975277067349 |
| HFD03 | Q2 | B3_label_only | 147 | 0 | 0.711764705882353 | 0.4482877754870756 | 0.4482877754870756 | 0.4533420138888889 | 0.001177300347222222 | 0.0 | -0.2803921568627451 |

Label-only uses no informative state observation:it samples the public E0 prior. Comparisons quantify the importance of supplied information under this estimator, not an additive causal decomposition of observed state versus context. HFD04 uses32 draws versus historical64:paired current B0—not99.22% historical headline—is the attribution baseline. Unknown-context prediction targets the hidden realized-context truth; it does not redefine known-context formal Health. All mechanisms remain separate. Exact label-to-prior map and cohort-share normalization are in Numerical_Methods_Contract.md.


## Roster-prior reconciliation (separately registered supplement)

The uniform-nine prior agrees with the 324-world factorial roster, but not the boundary-designed roster. The latter contains context counts [7,7,7,7,7,7,13,13,13]. Original uniform-mixture results remain above and unchanged. Supplemental prior weights are these counts divided by 81, chosen from the entire frozen roster before supplemental outcomes, with independent random streams and the same observations. They are experimental design frequencies, not ecological frequencies. All signed-margin, absolute-margin and failure-class strata remain in prior_reconciled_strata.csv.

| suite | query | n | certified_n | FP | FN | balanced_accuracy | viable_MAE | reserve_MAE | brier | estimator_unresolved_fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| correct | Q1 | 81 | 75 | 11 | 5 | 0.7943548387096775 | 0.20979949574411652 | 0.23318745177469136 | 0.12037760416666667 | 0.024691358024691357 |
| correct | Q2 | 81 | 75 | 18 | 5 | 0.7193627450980392 | 0.20979949574411652 | 0.23318745177469136 | 0.18376302083333335 | 0.024691358024691357 |
| hemlock_extra_hazard_0.12 | Q1 | 81 | 81 | 23 | 0 | 0.8203125 | 0.33733716423128857 | 0.46217628761574076 | 0.17487702546296297 | 0.024691358024691357 |
| hemlock_extra_hazard_0.12 | Q2 | 81 | 81 | 40 | 0 |  | 0.33733716423128857 | 0.46217628761574076 | 0.34467833719135804 | 0.024691358024691357 |
| entry_lognormal_annual_sd_1.2 | Q1 | 81 | 77 | 12 | 6 | 0.7727910238429172 | 0.2117624165099344 | 0.25074824580439814 | 0.1351461038961039 | 0.024691358024691357 |
| entry_lognormal_annual_sd_1.2 | Q2 | 81 | 76 | 25 | 6 | 0.611665004985045 | 0.2117624165099344 | 0.25074824580439814 | 0.2611019736842105 | 0.024691358024691357 |

Correct-kernel aggregate balanced accuracy is 79.4% Q1 and 71.9% Q2; it does not establish robust boundary discrimination or recovered context. Comparisons with original B2 include independent-stream sampling variation and prior change, not a pure causal effect.
