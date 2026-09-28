# AUC diagnostics

The procedure for selecting the number of changepoints uses AUC
confidence-interval lower bounds to assess whether neighboring segments can
be distinguished. Short segments or widely spaced held-out observations can
leave too few observations to compute an interval, even when an AUC point
estimate is available. `auc_diagnostics` records how often intervals are
unavailable, helping users interpret the selection results and review their
settings.

## Why use this selection approach?

`bcmlrNET_model_select()` follows the approach described in
[Section 3.4 of the bcmlr paper](https://arxiv.org/html/2603.25668v1#S3.SS4).
It fits a fixed, conservatively large number of changepoints using Gibbs
sampling and assesses neighboring segments with held-out observations. The
result is a **posterior distribution of the number of fitted changepoints
that appear to correspond to true changepoints**. When its mode is positive,
the function refits on the full series with that number of changepoints.

An alternative is to assign a prior to the number of changepoints and use
[reversible-jump MCMC](https://doi.org/10.1093/biomet/82.4.711) to explore models
with different parameter dimensions. Our approach avoids those
between-dimension moves and their proposal design, providing a computationally
tractable alternative based on fixed-dimensional Gibbs fits. Relative runtime
depends on the models and implementations; this is not a general speedup claim.
The distribution described here is constructed from posterior summaries of
the initial fit.

The trade-off is the need for enough observations for both fitting and
evaluation. Choose `min_size`, `max_num_cp`, and `thinning` together rather
than treating them as independent settings.

## Why can an interval be unavailable?

For nonconstant scores, the implementation uses pROC's DeLong calculation,
which [returns unavailable confidence limits when either group has fewer than two observations](https://raw.githubusercontent.com/xrobin/pROC/master/R/delong.R).
An AUC point estimate may still exist. Here, evaluation observations are
networks at evaluation time points, not nodes or posterior samples.

Other unavailable intervals can arise from calculation errors or non-finite
limits. Constant finite scores with both classes present are treated as
non-informative, with lower bound 0.5; non-finite prediction scores stop the
fit. These diagnostics assess interval availability, not interval coverage or
sampler convergence. An available interval need not have adequate precision.

## How unavailable intervals enter selection

For each posterior sample, a fitted changepoint receives selection indicator
**1** when its available AUC lower bound exceeds `threshold`; otherwise, its
indicator is **0**. The sum of these indicators is the number of fitted
changepoints that appear to correspond to true changepoints in that sample.

When an interval is unavailable, its indicator is set to zero: that fitted
changepoint is not included in this number. The posterior sample itself is
retained. This is a decision rule for an unevaluable comparison, not an AUC
value of zero or evidence that no change occurred. Unavailable intervals can
therefore affect the resulting posterior distribution.

## Read the diagnostic table

```r
fit$auc_diagnostics
```

Each row refers to an initially fitted changepoint in chronological order,
not an observation index or a changepoint from the final refit. Each retained
posterior sample contributes one interval evaluation per initial changepoint.
The summary is computed once after Stage 1 without recomputing intervals.

| Column | Meaning |
| --- | --- |
| `fitted_changepoint` | Position in the initial fit: `1, ..., max_num_cp`. |
| `n_evaluations` | Interval evaluations for this changepoint: `num_iter - num_warmup`. |
| `n_unavailable` | Evaluations that did not produce a usable AUC interval. |
| `pct_unavailable` | `100 * n_unavailable / n_evaluations`, on a 0--100 scale. |
| `n_unavailable_limited_samples` | Unavailable intervals with fewer than two evaluation observations in at least one neighboring segment. |
| `n_unavailable_other` | Unavailable intervals despite at least two evaluation observations in both segments; the cause requires further investigation. |

The last two columns sum to `n_unavailable`. They count interval evaluations,
not observations or segments. "Limited samples" means limited evaluation
observations, not too few posterior samples.

## Choose settings together

| Setting | Practical trade-off |
| --- | --- |
| `min_size` | Larger values exclude shorter segments, but require more observations and can rule out short-lived regimes. |
| `max_num_cp` | Fits this many changepoints initially, creating one more segment. A larger value allows more fitted changepoints but can leave fewer evaluation observations per segment. |
| `thinning` | For values above one, every `thinning`-th observation is held out. Larger values leave more data for fitting but fewer for evaluation. |

When `thinning = 1`, the same observations are used for fitting and evaluation;
this is not an independent holdout assessment. More posterior samples do not
supply additional evaluation observations within a given segmentation.

There is no universal best combination, especially for short series. Consider
the series length and plausible segment durations, inspect `auc_diagnostics`,
and compare nearby feasible settings. Average segment length alone does not
guarantee enough held-out observations in every segment. Preserve the paper's
settings for reproduction; report sensitivity analyses separately.

## Warning control

`warn_small_auc = FALSE` suppresses only the warning associated with fewer
than two evaluation observations in a neighboring segment. The diagnostic
table is still returned, other interval failures still warn, and an entirely
unavailable set stops selection. This changes reporting, not calculations.

[Back to the README](../README.md)
