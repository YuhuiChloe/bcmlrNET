
# bcmlrNET

**Bayesian changepoint analysis for dynamic networks**

Identify when network structure changes, quantify uncertainty, and examine
which network features are associated with the detected changes.

`bcmlrNET` is research software for analyzing time-ordered networks using
adjacency matrices and network summary statistics. Its application workflows
connect changepoint inference with network-level coefficient summaries,
dyad-level coefficient differences, and node-level involvement measured
through row-wise L2 norms.

> **Development status:** The research implementation is being organized into
> an R package for distribution through GitHub. Installation and the
> package-qualified examples below describe the target package workflow;
> they have not yet been validated against a packaged release.

<!-- Add one representative application figure here after its image file
has been committed. Add badges only for checks that actually exist. -->

## What the software provides

- **Changepoint inference:** Analyze a specified number of changepoints or
  assess uncertainty in their number through the model-selection workflow.
- **Interpretation at multiple levels:** Examine coefficients for network
  statistics, coefficient differences for individual dyads, and posterior
  summaries of node-level involvement.
- **Application workflows:** Follow the analyses of Enron email communication
  networks and POLECAT weekly conflict networks.

## Planned GitHub installation

Once the installable package is available in this repository, the installation
command will be:

```r
install.packages("remotes")
remotes::install_github("YuhuiChloe/bcmlrNET")
```

The examples use `bcmlrNET::function_name()` to identify package functions
explicitly. These calls require an installed package with the corresponding
functions exported.

## Main interface

The two intended entry points are:

| Analysis | Entry point |
| --- | --- |
| Fit conditional on a specified number of changepoints | `bcmlrNET::bcmlrNET()` |
| Analyze an unknown number of changepoints | `bcmlrNET::bcmlrNET_model_select()` |

The example below follows the model-selection call in the Enron research
script.

<!-- Add a fixed-number example after verifying the bcmlrNET() signature. -->

## Input data

The application workflow below supplies two aligned inputs:

| Input | Representation |
| --- | --- |
| `A` | A list of adjacency matrices in chronological order, with a common node set and consistent node ordering. |
| `S` | A numeric matrix or data frame with one row per network and one column per network statistic, in the same chronological order as `A`. |

The Enron example uses binary, undirected networks without self-loops and the
statistics `edges`, `triangles`, and `isolates`. Keep node labels in the same
order as the rows and columns of every adjacency matrix. Keep weekly interval
metadata in the same order as the network list.

## Model-selection workflow

This example assumes that `A_list` already contains the prepared Enron
network series. The settings are taken from the Enron application; they are
not generic defaults or a short demonstration run. The data-preparation
workflow must be completed before this example can be run.

```r
S <- bcmlrNET::compute_netstats(
  A_list,
  network_stats = c("edges", "triangles", "isolates")
)

fit <- bcmlrNET::bcmlrNET_model_select(
  A = A_list,
  S = S,
  min_size = 6,
  latent_dim = 10,
  eta_dir = 1.01,
  b_tau = 0.5,
  thinning = 5,
  prior_kappa = "default",
  alpha_f = 0.05,
  threshold = 0.5,
  max_num_cp = 8,
  sd_beta = 3,
  num_iter = 10000,
  num_warmup = 5000,
  print_progress = TRUE,
  print_outputs = FALSE,
  prior_for_U = "Spike_Slab",
  lambda_init = "equal",
  sparse_networks = FALSE,
  center_std_A = TRUE,
  seed = 2026
)

fit$Kappa_mode
fit$probs_num_cp
```

The application scripts use the following result components:

| Component | Use in the application workflow |
| --- | --- |
| `Kappa_mode` | Report estimated changepoint locations and map them to observation weeks. |
| `Kappa` | Supply changepoint samples to the location-uncertainty plots. |
| `num_cp_dist` | Plot the distribution of the number of fitted changepoints that appear to correspond to true changepoints. |
| `probs_num_cp` | Inspect the probabilities associated with different changepoint counts. |
| `Beta` | Summarize coefficient differences for the network statistics. |

With an aligned vector of weekly labels named `time_labels`, the reported
locations can be displayed as:

```r
time_labels[fit$Kappa_mode]
```

## From fitted results to figures

The planned plotting interface preserves the research functions used in the
application scripts:

| Figure or summary | Research function |
| --- | --- |
| Network-statistic trajectories | `plot_netstats_paper()` |
| Distribution of the number of changepoints | `plot_histogram()` |
| Credible intervals for network-statistic coefficient differences | `plot_beta_diff_CIs()` |
| Heatmaps of dyad-level coefficient differences | `plot_B_diff_heatmaps()` |
| Node involvement through row-wise L2 norms | `plot_node_involvement_L2()` |
| Network statistics with changepoint estimates and credible intervals | `plot_netstats_cp_cis()` |

For example, the count-distribution plot will use:

```r
p_count <- bcmlrNET::plot_histogram(
  num_vec = fit$num_cp_dist,
  max_num = 8,
  y_max = 1,
  title = "Posterior distribution of the number of changepoints"
)

print(p_count)
```

Here, the plotted count is the number of fitted changepoints that appear to
correspond to true changepoints under the selection procedure. It should not
be confused with the initial number supplied through `max_num_cp`.

The node-level plots summarize row-wise L2 norms of coefficient differences.
Nonzero posterior median coefficient differences are not, by themselves,
statements of statistical significance or causal effects.

## Applications and reproducibility

**Enron.** The research workflow prepares email communication networks,
aggregates email aliases to person-level nodes, and examines changepoints
alongside network-statistic, dyad-level, and node-level summaries.

**POLECAT.** The research workflow reads a weekly network array, aligns it with
node labels and a week-index map, and produces the corresponding changepoint
and coefficient-summary figures.

The package release is intended to distinguish two reproduction routes:

| Route | Scope |
| --- | --- |
| Full recomputation | Prepare the data, fit the models, and generate the figures. |
| Figure regeneration | Load saved fitted results and recreate the figures without refitting. |

Figure regeneration is not a substitute for rerunning model fitting. The
release documentation will identify the inputs and settings used for each
route, including seeds, software versions, node ordering, weekly interval
conventions, and output files. Availability of application datasets and saved
results will be documented separately from package installation.

## Development priorities

- Package the model and reusable summaries with documented inputs and outputs.
- Add a small runnable example, automated tests, and package checks.
- Separate application preparation, model fitting, competitor comparisons,
  and figure generation into reproducible workflows.

## Related work and attribution

This project forms part of dissertation research on changepoint analysis.
The related [bcmlr repository](https://github.com/YuhuiChloe/bcmlr) contains the
companion research implementation.

The Enron preprocessing script credits the
[CPD_generative Enron workflow](https://github.com/allenkei/CPD_generative/blob/main/R_code/Enron/Enron.R).
The research analyses use competing methods as comparison baselines, not as
components of the `bcmlrNET` method.

<!-- Before release: add the verified network-method manuscript citation,
software citation, approved license, and data-source documentation.
Do not infer the network paper's citation from the earlier bcmlr paper. -->

## Questions and reproducibility reports

When reporting an issue, include the command that failed, the error message,
a minimal example where possible, and the output of `sessionInfo()`. Identify
whether the issue arose during installation, model fitting, or figure
generation.
