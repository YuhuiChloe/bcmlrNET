# bcmlrNET

## When did the network change, and what changed with it?

`bcmlrNET` is R software for detecting and interpreting changepoints in dynamic
networks. It analyzes a sequence of adjacency matrices, optionally together
with network summary statistics, to identify changes in network structure
and quantify uncertainty.

A changepoint is only the beginning of the analysis. Posterior summaries help
examine the network characteristics, connections between pairs of nodes, and
individual nodes associated with each detected change.

- **Detect changes:** Estimate changepoint locations conditional on a specified
  number of changepoints, or use the procedure for selecting the number of
  changepoints when it is unknown.
- **Quantify uncertainty:** Examine uncertainty in the number and locations
  of changepoints, together with posterior summaries of coefficient changes.
- **Interpret changes:** Connect the detected changes to network statistics,
  dyad-level coefficient differences, and node-level involvement.

> **Status:** R package in development. The package-qualified examples below
> describe the interface being prepared for GitHub distribution; installation
> and end-to-end examples have not yet been validated in a packaged release.

<!-- Add a representative figure here when its image file is committed. -->

## Installation

Install the package from GitHub with:

```r
install.packages("remotes")
remotes::install_github("YuhuiChloe/bcmlrNET")
```
To accelerate posterior sampling, `bcmlrNET` implements computationally
intensive matrix operations and node-level parameter updates in C++ through
`Rcpp` and `RcppArmadillo`. The sampler also caches intermediate results to
avoid redundant calculations.

Users fit models and generate visualizations through R without writing C++
code. Building from source requires an R-compatible C++ toolchain.

### Installation troubleshooting

If `remotes::install_github()` fails with **HTTP error 401: Bad credentials**,
GitHub has rejected the authentication credentials sent by the installation
request. A common cause is an invalid, expired, or revoked personal access
token (PAT) saved in the user's R/Git environment. This error occurs before
the package is downloaded and compiled; changing the package's R or C++
functions will not resolve it.

`bcmlrNET` is public. Users can either update their saved credentials or
install from the public ZIP archive without a token.

#### Option 1: Replace the saved token

[Create a new token in your own GitHub account](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens),
then run these commands in an interactive R session:

```r
if (!requireNamespace("gitcreds", quietly = TRUE)) {
  install.packages("gitcreds")
}
gitcreds::gitcreds_set(url = "https://github.com")
```

Choose **Replace these credentials** if prompted, and enter the new token.
This updates the saved GitHub credential, which may also be used by other
GitHub tools on your computer. Enter the token only at the prompt; do not
put it in example scripts, README files, or issue reports. `gitcreds` requires
Git to be installed.

If you previously defined `GITHUB_PAT` or `GITHUB_TOKEN` in `.Renviron` or
another startup file, remove or update the obsolete entry there as well.
Otherwise, R may continue to use that token instead of the replacement in
the credential store. These commands show only whether the variables are
set and then clear them for the current R session:

```r
nzchar(Sys.getenv(c("GITHUB_PAT", "GITHUB_TOKEN")))
Sys.unsetenv(c("GITHUB_PAT", "GITHUB_TOKEN"))
```

`Sys.unsetenv()` does not edit startup files. After correcting any saved
settings, restart R/RStudio and retry installation:

```r
remotes::install_github(
  "YuhuiChloe/bcmlrNET",
  dependencies = NA,
  force = TRUE,
  upgrade = "never"
)
```

See the [gitcreds documentation](https://gitcreds.r-lib.org/reference/gitcreds_get.html)
for details about credential storage and environment variables.

#### Option 2: Install from the public ZIP archive

This method downloads the source archive directly and avoids the GitHub API
authentication used by `install_github()`:

```r
remotes::install_url(
  "https://github.com/YuhuiChloe/bcmlrNET/archive/refs/heads/main.zip",
  dependencies = NA,
  force = TRUE,
  upgrade = "never"
)
```

GitHub generates this ZIP from the `main` branch. R downloads and extracts
it into temporary storage, then builds and installs the package. No manual
download or token is required for this route, and the stored credentials
are not changed.

If R cannot download the archive, [download the ZIP in your browser](https://github.com/YuhuiChloe/bcmlrNET/archive/refs/heads/main.zip),
then run the following and select the downloaded ZIP:

```r
remotes::install_local(
  file.choose(),
  dependencies = NA,
  force = TRUE,
  upgrade = "never"
)
```

Both methods build the package from source and require the C++ toolchain
described above. After installation succeeds, check that the package and
both prepared datasets load:

```r
library(bcmlrNET)
data("POLECAT", package = "bcmlrNET")
data("Enron", package = "bcmlrNET")

length(POLECAT$A_list)  # 184
length(Enron$A_list)    # 100
dim(Enron$A_list[[1]])  # 89 89
```


## Input data

The prepared POLECAT dataset supplies the networks and aligned country and
weekly labels used in the dissertation example:

```r
data("POLECAT", package = "bcmlrNET")
A_list <- POLECAT$A_list
node_labels <- POLECAT$node_labels
time_labels <- POLECAT$time_labels
date_list <- POLECAT$date_list
```

Let `A_list` be a list of `T` adjacency matrices in chronological order. Each
matrix has dimensions `n x n` and represents a network on the same `n` nodes.
The node ordering must be identical across all matrices.

The current adjacency-matrix implementation requires symmetric matrices with
zero diagonal entries and finite values: undirected networks without
self-loops. The network-statistic example below uses binary, zero-one
adjacency matrices.

An optional input, `S`, is a numeric matrix with `T` rows and one column per
network statistic. Row `t` of `S` must correspond to `A_list[[t]]`.

For example, compute the number of edges, triangles, and isolated nodes at
each time point:

```r
S <- bcmlrNET::compute_netstats(
  A_list,
  network_stats = c("edges", "triangles", "isolates")
)
```

To analyze adjacency matrices without additional network statistics, pass
`S = NULL` explicitly to either fitting function.

## Fit the method

The following calls assume that `A_list` contains your network sequence and
that `S` has been constructed as above. The settings illustrate the interface;
they are not universal tuning recommendations. The series must be long enough
to accommodate the requested number of changepoints and minimum segment length.

### A specified number of changepoints

Use `bcmlrNET()` to fit conditional on a chosen number of changepoints. Here,
`num_CP = 2` specifies two changepoints and three segments.

```r
fit_fixed <- bcmlrNET::bcmlrNET(
  A = A_list,
  S = S,
  num_CP = 2,
  min_size = 10,
  latent_dim = 10,
  prior_for_U = "Spike_Slab",
  num_iter = 5000,
  num_warmup = 2500,
  seed = 2026
)

fit_fixed$Kappa_mode
```

`Kappa_mode` contains the point estimates (marginal posterior modes) of the
changepoint locations, expressed as indices in the original network sequence.

### An unknown number of changepoints

Use `bcmlrNET_model_select()` when the number of changepoints is unknown.
`max_num_cp` specifies the number initially fitted for the selection procedure,
not the final number the method must return.

```r
fit <- bcmlrNET::bcmlrNET_model_select(
  A = A_list,
  S = S,
  max_num_cp = 5,
  min_size = 10,
  latent_dim = 10,
  prior_for_U = "Spike_Slab",
  thinning = 5,
  alpha_f = 0.05,
  threshold = 0.5,
  num_iter = 5000,
  num_warmup = 2500,
  seed = 2026
)

fit$num_cp
fit$probs_num_cp
```

The procedure assesses whether networks on opposite sides of each fitted
changepoint can be distinguished using held-out observations. Here,
`thinning = 5` holds out every fifth observation for this assessment; it does
not mean retaining every fifth MCMC sample.

The selected number of changepoints is the mode of the distribution of the
number of fitted changepoints supported by the selection criterion. When this
number is positive, the function automatically refits using the full network
sequence. Changepoint locations and coefficient summaries from that refit are
conditional on the selected number of changepoints, rather than averaged over
different numbers of changepoints.

If zero changepoints are selected, the function returns summaries of the
number of changepoints without a refitted model. Location and coefficient plots
should then be skipped.

### Main outputs

| Component | Contents |
| --- | --- |
| `num_cp` | Selected number of changepoints; returned by model selection. |
| `probs_num_cp` | Probabilities for the number of changepoints, from zero to `max_num_cp`. |
| `num_cp_dist` | Posterior samples of the number of fitted changepoints that appear to correspond to true changepoints. |
| `Kappa_mode` | Point estimates (marginal posterior modes) of changepoint locations. |
| `Kappa` | Posterior samples of changepoint locations. |
| `Beta` | Posterior samples of coefficients for the network statistics. |
| `U`, `Lambda` | Posterior samples used to construct dyad-level coefficient matrices. |

In model-selection output, `Kappa_mode`, `Kappa`, `Beta`, `U`, and `Lambda`
are present only when the selected number of changepoints is positive.

The model-selection output also includes `fit$auc_diagnostics`, which
summarizes the availability of AUC confidence intervals during the initial
fit. See [AUC diagnostics](documentation/auc-diagnostics.md) for column
definitions, the treatment of unavailable intervals, and warning control.

## What changed, and how certain are we?

Start with uncertainty in the number of changepoints:

```r
p_num_cp <- bcmlrNET::plot_histogram(
  num_vec = fit$num_cp_dist,
  max_num = 5,
  title = "Posterior distribution of the number of changepoints"
)
print(p_num_cp)
```

The following figure shows the posterior distribution of the number of fitted
changepoints that appear to correspond to true changepoints in the POLECAT
application.

![POLECAT posterior distribution of the number of changepoints](man/figures/POLECAT_DistNumCP.png)

Next, examine changes in the coefficients of the network statistics. This plot
shows posterior means and 95% credible intervals for differences between
adjacent segments:

```r
if (fit$num_cp > 0) {
  p_beta <- bcmlrNET::plot_beta_diff_CIs(
    Beta = fit$Beta,
    coef_labels = colnames(S),
    mfrow = c(1, fit$num_cp),
    ci_level = 0.95
  )
  print(p_beta)
}
```

The visualization functions cover complementary aspects of the analysis:

| Question | Function |
| --- | --- |
| How do network statistics evolve over time? | `plot_netstats_paper()` |
| How uncertain is the number of changepoints? | `plot_histogram()` |
| Where are the changepoints, and how uncertain are their locations? | `plot_netstats_cp_cis()` |
| Which network-statistic coefficients change? | `plot_beta_diff_CIs()` |
| Which dyad-level coefficients change? | `plot_B_diff_heatmaps()` |
| Which nodes have the largest row-wise coefficient changes? | `plot_node_involvement_L2()` |

The optional competing-method panel in `plot_netstats_cp_cis()` requires
`patchwork` (`install.packages("patchwork")`). Omit `other_cps`, or use
`other_cps = list()`, when no competing methods are available.

Node involvement is summarized using row-wise L2 norms of differences
between coefficient matrices, calculated for each posterior sample. These summaries
characterize changes associated with nodes; they do not establish causation.

## Applications

### Enron: changes in email communication

The [Enron email corpus](https://www.cs.cmu.edu/~enron/) contains email
communications involving Enron employees. Our application uses the network
data supplied by `igraphdata` to construct email communication networks over
non-overlapping weekly intervals. Nodes represent unique employees (89 in total) after
email aliases are combined, and edges represent email communication between
them during the corresponding week.

The analysis examines when communication patterns change and which employees
and employee pairs are most involved in the detected changes.

See [Example_code.R](Example_code.R) for the Enron fitting and plotting workflow.

### POLECAT: changes in international conflict networks

[POLECAT](https://andrewhalterman.com/publication/plover-polecat-new-event-data/)
is a political event dataset constructed from international news reports and
coded according to the PLOVER event ontology. Our application focuses on
material-conflict events, including assaults, coercive actions, mobilizations,
protests, and reductions in relations (including sanctions). Throughout this
application, we refer to these material-conflict events as **conflicts**.

The networks are constructed over non-overlapping weekly intervals. Nodes
represent countries (197 in total), and edges represent reported conflicts between them
during the corresponding week.

The analysis examines when these conflict relations change and uses posterior
summaries to characterize the countries and country pairs associated with
the detected changes.

See [Example_code.R](Example_code.R) for the POLECAT fitting and plotting workflow.

## Reproduce an application

Follow the analysis from its inputs to its figures:

1. **Prepare the networks.** Construct the adjacency-matrix list and align its
   node labels and weekly intervals.
2. **Fit the method.** Compute the chosen network statistics and run the
   fitting function with the application's settings and random seed.
3. **Generate the figures.** Use the returned object to summarize changepoints
   and visualize network-, dyad-, and node-level coefficient changes.

For a reproducible run, use the same data version, preprocessing steps,
software versions, and computational settings. Record the package release or
commit, the random seed, and the R session information. A seed alone does not
guarantee identical numerical results across different computing environments.

## Related work and attribution

The companion [bcmlr project](https://github.com/YuhuiChloe/bcmlr) provides the
related research implementation for multivariate series.

The Enron preprocessing builds on the
[CPD_generative Enron workflow](https://github.com/allenkei/CPD_generative/blob/main/R_code/Enron/Enron.R).
Source datasets and comparison methods should be credited separately from
the `bcmlrNET` software.

<!-- Before release, add the verified method/software citations and approved
license. Add links to application scripts only after they are committed. -->

## Questions and issues

When reporting a problem, include a minimal example, the error message, the
package version or commit, and the output of `sessionInfo()`.
