# MCS-NMF standalone validation report

Status date: 2026-11-15.

## Completion

- 9 datasets × 20 common masks × 12 included methods = 2160 method-runs.
- All 2160 predictions are finite, match the corresponding saved floor-10% mask, and were independently re-evaluated on unlabeled samples.
- Six metrics use one implementation: ACC, geometric-entropy NMI, Purity, pair-counting F-score, Precision, and Recall.
- PNLPSCONMF is excluded from ranking because its reproduction remains unresolved; its files and diagnostic report are preserved separately.
- No ScaleSem-NMF manuscript file was changed.

## Position of MCS-NMF among included methods

MCS-NMF ranks first on ACC and Purity for 8 of 9 datasets, and first on NMI, F-score, and Recall for all 9 datasets. On AR, LpCNMF has the highest mean ACC (65.44% versus MCS-NMF 65.18%); MCS-NMF remains first on AR NMI, F-score, and Recall. These ranks apply only to the 12 included, successfully evaluated methods and must not be generalized to PNLPSCONMF.

Against HLPNMF, the mean paired ACC difference is positive on all nine datasets. The 95% paired bootstrap interval crosses zero on Optdigits ACC and NMI; the MSRA25 NMI interval also includes values very close to zero. The other reported ACC/NMI intervals are positive under this protocol.

## Palm qualification

Palm MCS-NMF and HLPNMF use the same existing floor-protocol masks and p=10 prediction files. The p=10 setting was transferred from the reference Palm setting after the initial p=5 audit and must be disclosed as a reference-transferred sensitivity setting. It is not evidence of a preregistered Palm parameter choice.

## Interpretation boundary

This package supports a standalone common-budget comparison under the explicitly stated floor-10% protocol. It does not reproduce the reference paper's ceil-rounding sampler, does not require numerical identity with its tables, and does not establish a comparison against the unresolved PNLPSCONMF implementation. Any manuscript claim must preserve these limits.

## Files

- `per_run_six_metrics.csv`: all 2160 independently recomputed rows.
- `summary_long.csv`: mean and standard deviation by dataset and method.
- `main_table_formatted.csv`: publication-style six-metric entries.
- `paired_differences.csv`: MCS-NMF minus every included baseline with paired bootstrap intervals.
- `MCSNMF_vs_HLPNMF.csv`: focused paired comparison.
- `prediction_audit.csv`: source file, mask match, cluster count, and largest-cluster fraction.
- `MCSNMF_vs_HLPNMF_paired.png`: paired ACC/NMI visualization.
