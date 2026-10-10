# LPA and BCH helper functions for tidySEM

Reusable R functions for running a **latent profile analysis (LPA)** with
[tidySEM](https://cjvanlissa.github.io/tidySEM/) and then comparing profiles on
distal outcomes with the **BCH** method.

Helper functions were created as I needed to run the BCH analysis with over 10 variables, and I found it more practical to use functions to do so.

## Files

| File | What it is |
|---|---|
| `LPA_BCH_HelperFunctions.R` | The functions. Load with `source("LPA_BCH_HelperFunctions.R")`. |
| `Example - LPA and BCH Helper Functions.Rmd` | A worked example using simulated data. Knit it to see every function in use. |

## Installation

```r
install.packages(c("tidySEM", "OpenMx", "dplyr", "ggplot2"))
install.packages("ggpattern")  # optional, for black-and-white bar charts
```

## Functions

**LPA**

| Function | Purpose |
|---|---|
| `run_lpa()` | Fit models across numbers of profiles and variance structures; optionally cache to disk |
| `fit_table()` | Fit indices (AIC, BIC, entropy, ...) with EEI/VVI labels |
| `class_sizes()` | Posterior-weighted size of each profile |
| `avg_posterior()` | Average posterior probability of the most likely class |
| `extract_profile_means()` | Indicator means for each profile |
| `plot_lpa_profiles()` | Line or bar plot of profile means |

**BCH**

| Function | Purpose |
|---|---|
| `run_bch()` | BCH means for one outcome |
| `bch_omnibus()` | Overall test across profiles, with an η² approximation |
| `bch_pairwise()` | Pairwise profile comparisons with significance stars |
| `bch_effect_sizes()` | Pairwise Cohen's *d* |
| `extract_bch_means()` | Means and 95% CIs as a data frame |
| `plot_bch_outcomes()` | Bar chart of several outcomes, in colour or black-and-white patterns |
| `bch_full()` | Runs all BCH steps for one outcome in one call |

## Quick start
## Replace my_data, var1-var3 and outcome with your own data and variable names

```r
source("LPA_BCH_HelperFunctions.R")

profiles <- run_lpa(my_data, indicators = c("var1", "var2", "var3"), classes = 1:6)
fit_table(profiles)

chosen <- profiles[["equal var 4"]]
profile_labels <- c(class1 = "Low", class2 = "Moderate", class3 = "High", class4 = "Mixed")

plot_lpa_profiles(chosen, profile_labels = profile_labels)

res <- bch_full(chosen, my_data$outcome, "Outcome", profile_labels = profile_labels)
```

## Notes

- **Profile labels** are a named vector matching tidySEM's class names, e.g.
  `c(class1 = "Label 1", class2 = "Label 2")`. Check the profile means first,
  since class numbering can change between runs.
- **Outcome order:** the outcome passed to the BCH functions must have the same
  rows, in the same order, as the data used to fit the LPA.
- **Caching:** if `cache_file` already exists, `run_lpa()` loads it instead of
  refitting. Delete the file, or use a new name, after changing the data or settings.
- **Decimal places:** estimates are shown to 3 decimal places by default. Change
  this with the `digits` argument, e.g. `extract_profile_means(chosen, digits = 4)`.
- **Effect sizes:**
  - η² for the omnibus test is approximated as (χ² − df) / (χ² + N).
  - Cohen's *d* uses the pooled SD of the BCH-estimated variances, weighted by the
    posterior class sizes.

## Data

- No actual research data is included in this repository.
- See bolded text in the example document to understand the code and identify where you need to insert your dataset and your variable names.

## Citation

If you use these functions, please cite tidySEM:
Van Lissa, C.J. tidySEM: Tidy Structural Equation Modeling. R package.
https://cjvanlissa.github.io/tidySEM/