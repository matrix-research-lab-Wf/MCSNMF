# MCSNMF

Official MATLAB implementation and reproducibility material for **Mass-Conserving Semantic Nonnegative Matrix Factorization (MCS-NMF)**.

## Contents

- `code/mcsnmf.m`: core MCS-NMF optimizer.
- `code/build_high_order_laplacian.m`: high-order hypergraph operator.
- `code/run_mcsnmf.m`: dataset-level entry point.
- `experiments/run_post_sinkhorn_comparison.m`: paired Row-only, Post-Sinkhorn, and MCS-NMF experiment.
- `experiments/run_all_proposed_method.m`: nine-dataset reproduction entry point.
- `protocol/`: the 20 paired labeled/unlabeled splits used in the main comparison.
- `results/`: six-metric summaries, paired statistics, runtime, robustness, ablation, and sensitivity results.
- `manuscript/`: LaTeX source, figures, Highlights, and compiled manuscript.

## Requirements

- MATLAB R2019a or later.
- No GPU is required.

## Data format

Place each dataset in `data/` as a MAT file containing either:

- `X` (features by samples) and `y` (samples by 1), or
- `fea` (samples by features) and `gnd` (samples by 1).

The `data/` directory contains the nine datasets used in the experiments: `AR.mat`, `PIE.mat`, `MSRA25.mat`, `YaleB.mat`, `COIL20.mat`, `COIL100.mat`, `MNIST.mat`, `Palm.mat`, and `Optdigits.mat`.

## Quick start

```matlab
addpath('code');
result = run_mcsnmf('data/COIL20.mat', 'protocol/COIL20/COIL20_seed20260627.mat', 2, 3);
disp(result.metrics)
```

The last two arguments are the hyperedge neighborhood size `p` and diffusion order `k`. The main experiments use `alpha = 1000`, at most 100 block iterations, and the saved paired supervision masks.

The mathematical basis update is the standard multiplicative update. Numerical
division uses a `1e-12` denominator floor, and a basis candidate is accepted only
when the complete objective does not increase. The denominator audit is provided
in `results/optimizer_denominator_audit.csv`.

## Main protocol

The main comparison uses class-stratified supervision with `ceil(0.1*n_q)` labeled samples in each class, 20 paired splits, and evaluation on unlabeled samples. The six reported measures are ACC, NMI, Purity, F-score, Precision, and Recall.

## Reproducing the paper tables

The exact paired splits are under `protocol/`. Aggregated values used by the manuscript are under `results/`. Third-party baseline implementations are not republished; their original implementations and method-specific settings are identified in the manuscript.

To reproduce the proposed-method runs on all nine datasets:

```matlab
addpath('code');
addpath('experiments');
run_all_proposed_method;
summarize_post_sinkhorn;
```

For one dataset:

```matlab
run_post_sinkhorn_comparison('COIL20');
```

This experiment uses the same 20 saved masks for three paired variants:

- `Row-only`: simplex rows without the class-mass column marginal;
- `Post-Sinkhorn`: one marginal projection applied after Row-only training;
- `MCS-NMF`: the class-mass marginal enforced throughout optimization.

The completed paired runs and 10,000-resample interval summary are stored in
`results/post_sinkhorn_paired_runs.csv` and
`results/post_sinkhorn_summary.csv`, respectively.

## Repository scope

This repository contains the proposed method, reproducibility protocol, result summaries, figures, and manuscript source. Temporary files, QA renderings, and machine-specific paths are excluded.
