# Juvenile threshold and query sufficiency

Secure observed E1 subset only; dated support is retained in observed_groups.csv. No claim of measured2020 origin juvenile counts.

| n | FP | FN | membership_error_fraction |
| --- | --- | --- | --- |
| 53611 | 1187 | 18 | 0.022476730521721288 |

RMS and stem classification differ even when count/BA agree. Direct Gaussian perturbations0/.1/.25/.5/1cm are methodological stress conditions, not calibrated instrument-error laws. Membership errors by distance and all8 replicates are in measurement_error.csv. M0–M3 masking retains physically sanitized fitting inputs; profile data exclude heldout DBH.

| mask | n | truth_juveniles | mean_predicted_juveniles | delta_J_over_J0 | membership_brier |
| --- | --- | --- | --- | --- | --- |
| M0 | 9988 | 4636 | 4500.37890625 | -0.0016523641671844732 | 0.11755414633546583 |
| M1 | 13133 | 5329 | 5039.02734375 | -0.00353293439392278 | 0.11024708865139642 |
| M2 | 11522 | 5270 | 5007.375 | -0.003199739269222803 | 0.12478963477477868 |
| M3 | 8131 | 3967 | 3578.7734375 | -0.00473002866210997 | 0.13269558288326505 |

Paired enriched-state count/BA projection max absoluteBA error=6.821210263296962e-13; juvenile counts do NOT commute. The explicit projection_counterexample.json has identical count/RMS/BA and opposite actual Q2 answers. Thus RMS is not generally sufficient for the juvenile query.

Finite-design paired Health disagreement=0.001953125. Decision-change categories:

| category | cells |
| --- | --- |
| COMBINED | 3 |
| NONE | 1533 |

The enriched distribution is inferred from visible E1 profiles and rescaled to reconstructed origin moments; it is not a newly observed modality. Initial projection is numerically exact; future stochastic kernels are not claimed to commute or to be a conservative refinement. Common seeds are paired, not exact common uniforms after changing array dimensions. Size-aware representation preserves an explicit juvenile distinction but does not establish full latent ecological sufficiency. Status CONDITIONAL, not an operational upgrade endorsement.


## Local effect, not just the global disagreement rate

Only 3/1,536 finite-grid decisions differ (0.195%), but one stressed state loses 67/256 = 26.17 percentage points of viable and reserve-qualified mass. All three changes occur in state 8, D3, horizon 5, beta 0.4, rho 0.5: gamma thresholds 0.75/0.9 and normal threshold 0.9. There are two resolved opposite kernel decisions and one fine-kernel MC-unresolved crossing.

| state | scenario | variant | horizon | theta | coarse_nv | fine_nv | coarse_nr | fine_nr | coarse_status | fine_status | viable_change | explanation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 8 | D3 | gamma_growth | 5 | 0.75 | 256 | 189 | 256 | 189 | TRUE | MC_UNRESOLVED | -0.26171875 | CONTINUATION_DRIVEN |
| 8 | D3 | gamma_growth | 5 | 0.9 | 256 | 189 | 256 | 189 | TRUE | FALSE | -0.26171875 | CONTINUATION_DRIVEN |
| 8 | D3 | normal_growth | 5 | 0.9 | 256 | 198 | 256 | 198 | TRUE | FALSE | -0.2265625 | CONTINUATION_DRIVEN |

The original COMBINED flag means both component thresholds crossed; it is not an independent reserve mechanism. In these three cases nr=nv in both representations, so the reserve-qualified loss is continuation-driven. The unchanged conditional reserve fraction does not prove reserve invariance elsewhere. The small global average cannot establish query sufficiency or dynamic invariance.
