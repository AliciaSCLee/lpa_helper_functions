
# =============================================================================
# Helper functions for Latent Profile Analysis (LPA) and BCH auxiliary-variable
# analysis with tidySEM
#
# Author: Alicia S.C. Lee
#
# Usage: source("LPA_BCH_HelperFunctions.R")
#
# Requires: tidySEM, OpenMx, dplyr, ggplot2, ggpattern (only for
#           plot_bch_outcomes(style = "pattern"))
# =============================================================================

library(tidySEM)
library(OpenMx)
library(dplyr)
library(ggplot2)

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

# Replace class names ("class1", ...) with readable labels if supplied.
.apply_labels <- function(x, profile_labels = NULL) {
  if (is.null(profile_labels)) return(x)
  out <- unname(profile_labels[as.character(x)])
  ifelse(is.na(out), as.character(x), out)
}

# Order profiles as given in profile_labels (or class order if NULL).
.profile_factor <- function(x, profile_labels = NULL) {
  levels <- if (is.null(profile_labels)) sort(unique(as.character(x))) else unname(profile_labels)
  factor(x, levels = levels)
}

# Remove significance stars and convert to numeric.
.num <- function(x) as.numeric(gsub("\\*+", "", x))

# Accept either a raw BCH object or the list returned by run_bch().
.get_bch <- function(x) if (is.list(x) && !is.null(x$bch)) x$bch else x

# -----------------------------------------------------------------------------
# 1. LPA
# -----------------------------------------------------------------------------

#' Fit a set of latent profile models (optionally cached to disk)
#'
#' @param data         Data frame.
#' @param indicators   Character vector of indicator column names.
#' @param classes      Numbers of profiles to fit, e.g. 1:6.
#' @param variances    "equal", "varying", or both.
#' @param covariances  "zero" (local independence), "equal", or "varying".
#' @param cache_file   Optional .rds path. If it exists, models are loaded from
#'                     it instead of being refit; otherwise results are saved.
#' @param seed         Random seed for reproducibility.
#' @return A list of fitted mx models (as returned by tidySEM::mx_profiles).

run_lpa <- function(data,
                    indicators,
                    classes     = 1:6, ## Default, can be overriden
                    variances   = c("equal", "varying"),
                    covariances = "zero", ## Covariance of 0 = I
                    cache_file  = NULL,
                    seed        = 12345) {
  
  if (!is.null(cache_file) && file.exists(cache_file)) {
    message("Loading cached models from ", cache_file)
    return(readRDS(cache_file))
  }

  ## Check if indicators exist in dataset
  
  missing_vars <- setdiff(indicators, names(data))
  if (length(missing_vars) > 0) {
    stop("Indicators not found in data: ", paste(missing_vars, collapse = ", "))
  }
  
  start_time <- Sys.time()
  set.seed(seed)
  
  profiles <- do.call(mx_profiles, list(
    data        = data[, indicators, drop = FALSE],
    classes     = classes,
    variances   = variances,
    covariances = covariances,
    expand_grid = TRUE
  ))
  
  if (!is.null(cache_file)) saveRDS(profiles, cache_file)
  
  message("LPA finished in ", format(round(Sys.time() - start_time, 1)))
  profiles
}


#' Model fit table with short model-type codes (EEI / VVI)
#'
#' @param profiles  Output of run_lpa().
#' @param columns   Fit columns to keep. NULL keeps all.

fit_table <- function(profiles,
                      columns = c("Model", "Name", "warning", "LL", "Parameters",
                                  "n", "AIC", "BIC", "saBIC", "Entropy",
                                  "prob_min", "n_min", "n_max")) {
  fit <- as.data.frame(table_fit(profiles)) %>%
    mutate(Model = Name %>%
             gsub("equal var", "EEI", .) %>%
             gsub("free var", "VVI", .))
  
  if (!is.null(columns)) fit <- fit[, intersect(columns, names(fit)), drop = FALSE]
  fit
}

#' Estimated class sizes (posterior-weighted counts)
#'
#' @return Named numeric vector, e.g. c(class1 = 1813.2, class2 = 955.0, ...).

class_sizes <- function(profile) {
  sp <- class_prob(profile)$sum.posterior
  setNames(sp$count, sp$class)
}

#' Average posterior probability for the most likely class (diagonal)

avg_posterior <- function(profile, digits = 3) {
  m <- class_prob(profile, type = "avg.mostlikely")$avg.mostlikely
  round(diag(as.matrix(m)), digits)
}

#' Extract indicator means for each profile
#'
#' @param profile         A single fitted model, e.g. profiles[["equal var 5"]].
#' @param strip_suffix    Optional regex removed from variable names (e.g. "_W1$").
#' @param profile_labels  Optional named vector: c(class1 = "Label 1", ...).

extract_profile_means <- function(profile,
                                  strip_suffix   = NULL,
                                  profile_labels = NULL,
                                  digits = 3) {
  out <- as.data.frame(table_results(profile, columns = c("label", "est", "confint", "class"),
                                     digits = digits)) %>%
    filter(grepl("^Means", label)) %>%
    mutate(
      Variable = gsub("^Means\\.", "", label),
      est      = .num(est)
    )
  
  if (!is.null(strip_suffix)) out$Variable <- gsub(strip_suffix, "", out$Variable)
  
  sizes <- class_sizes(profile)
  props <- sizes / sum(sizes)
  
  out %>%
    mutate(
      Proportion = round(100 * props[class], 1),
      Profile    = .apply_labels(class, profile_labels),
      Label      = paste0(Profile, " (", Proportion, "%)")
    ) %>%
    select(class, Profile, Label, Proportion, Variable, est)
}


#' Plot profile means
#'
#' @param profile         A single fitted model.
#' @param type            "line" (profile plot) or "bar".
#' @param var_labels      Optional named vector to rename indicators,
#'                        e.g. c(SPP_W1 = "SPP"). Also sets their order.
#' @param profile_labels  Optional named vector: c(class1 = "Label 1", ...).
#' @param strip_suffix    Optional regex removed from variable names.
#' @param title,y_label   Plot text.

plot_lpa_profiles <- function(profile,
                              type           = c("line", "bar"),
                              var_labels     = NULL,
                              profile_labels = NULL,
                              strip_suffix   = NULL,
                              title          = NULL,
                              y_label        = "Estimated mean") {
  type <- match.arg(type)
  d <- extract_profile_means(profile, strip_suffix, profile_labels)
  
  if (!is.null(var_labels)) {
    d$Variable <- ifelse(d$Variable %in% names(var_labels),
                         var_labels[d$Variable], d$Variable)
    d$Variable <- factor(d$Variable, levels = unique(c(unname(var_labels), d$Variable)))
  } else {
    d$Variable <- factor(d$Variable, levels = unique(d$Variable))
  }
  
  # Keep profile order consistent with profile_labels / class order
  
  label_order <- d %>% mutate(P = .profile_factor(Profile, profile_labels)) %>%
    arrange(P) %>% pull(Label) %>% unique()
  d$Label <- factor(d$Label, levels = label_order)
  
  p <- if (type == "line") {
    ggplot(d, aes(x = Variable, y = est, group = Label, colour = Label, shape = Label)) +
      geom_line(linewidth = 1.1) +
      geom_point(size = 3) +
      labs(colour = "Profile", shape = "Profile") +
      theme_classic() +
      theme(legend.position = "bottom")
  } else {
    ggplot(d, aes(x = Label, y = est, fill = Variable)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.65) +
      labs(fill = NULL) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  p +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    labs(x = NULL, y = y_label, title = title)
}


# -----------------------------------------------------------------------------
# 2. BCH analysis of distal outcomes
# -----------------------------------------------------------------------------

#' Run BCH for one outcome and return the profile means
#'
#' @param profile         A single fitted model.
#' @param outcome         Numeric vector, same rows/order as the LPA data.
#' @param outcome_name    Name used in printed output and plots.
#' @param profile_labels  Optional named vector: c(class1 = "Label 1", ...).
#' @return list(bch = BCH model, results = data frame of means, name = outcome_name)

run_bch <- function(profile, outcome, outcome_name = "Outcome",
                    profile_labels = NULL, verbose = TRUE) {
  
  bch <- BCH(profile, data = outcome)
  
  results <- as.data.frame(table_results(bch)) %>%
    filter(grepl("^Means", label)) %>%
    mutate(
      Profile = .apply_labels(group, profile_labels),
      est     = .num(est_sig)
    ) %>%
    select(Profile, group, est, est_sig, se, pval, confint)
  
  if (verbose) {
    cat("\nBCH means for", outcome_name, ":\n")
    print(results)
  }
  
  list(bch = bch, results = results, name = outcome_name)
}


#' Omnibus Wald/LR test across profiles, with an eta-squared approximation
#'
#' eta2 = (chi2 - df) / (chi2 + N)
#'
#' @param bch  Output of run_bch() (or a raw BCH object).
#' @param N    Sample size, e.g. nrow(your_data).

bch_omnibus <- function(bch, N, outcome_name = NULL, verbose = TRUE) {
  lr <- lr_test(.get_bch(bch))
  omnibus <- as.data.frame(lr$overall)
  omnibus$eta2 <- round((omnibus$LL_dif - omnibus$df) / (omnibus$LL_dif + N), 3)
  
  if (verbose) {
    cat("\nOmnibus test", if (!is.null(outcome_name)) paste("for", outcome_name), ":\n")
    print(omnibus)
  }
  omnibus
}

#' Pairwise comparisons between profiles
bch_pairwise <- function(bch, profile_labels = NULL, outcome_name = NULL, verbose = TRUE) {
  lr <- lr_test(.get_bch(bch))
  
  pairwise <- as.data.frame(lr$pairwise) %>%
    mutate(
      Profile1 = .apply_labels(Model1, profile_labels),
      Profile2 = .apply_labels(Model2, profile_labels),
      sig = case_when(
        p < .001 ~ "***",
        p < .01  ~ "**",
        p < .05  ~ "*",
        TRUE     ~ ""
      )
    ) %>%
    select(Profile1, Profile2, p, sig)
  
  if (verbose) {
    cat("\nPairwise comparisons", if (!is.null(outcome_name)) paste("for", outcome_name), ":\n")
    print(pairwise)
  }
  pairwise
}

#' Pairwise Cohen's d using BCH means/variances and weighted class sizes
#'
#' @param class_n  Named vector of class sizes; use class_sizes(profile).
bch_effect_sizes <- function(bch, class_n, profile_labels = NULL,
                             outcome_name = NULL, verbose = TRUE) {
  res <- as.data.frame(table_results(.get_bch(bch)))
  
  means <- res %>% filter(grepl("^Means", label)) %>%
    transmute(group, mean = .num(est_sig))
  vars  <- res %>% filter(grepl("^Variances", label)) %>%
    transmute(group, variance = .num(est_sig))
  
  s <- left_join(means, vars, by = "group") %>%
    mutate(sd = sqrt(variance), n = class_n[group])
  
  pairs <- combn(s$group, 2, simplify = FALSE)
  
  d <- do.call(rbind, lapply(pairs, function(p) {
    a <- s[s$group == p[1], ]
    b <- s[s$group == p[2], ]
    sd_pooled <- sqrt(((a$n - 1) * a$sd^2 + (b$n - 1) * b$sd^2) / (a$n + b$n - 2))
    data.frame(
      Profile1 = .apply_labels(p[1], profile_labels),
      Profile2 = .apply_labels(p[2], profile_labels),
      d        = round((a$mean - b$mean) / sd_pooled, 3)
    )
  }))
  
  if (verbose) {
    cat("\nPairwise Cohen's d", if (!is.null(outcome_name)) paste("for", outcome_name), ":\n")
    print(d)
  }
  d
}

#' Extract BCH means with confidence intervals (for plotting / tables)
extract_bch_means <- function(bch, outcome_name = NULL, profile_labels = NULL) {
  if (is.null(outcome_name) && is.list(bch) && !is.null(bch$name)) outcome_name <- bch$name
  
  as.data.frame(table_results(.get_bch(bch))) %>%
    filter(grepl("^Means", label)) %>%
    mutate(
      est      = .num(est_sig),
      CI_lower = as.numeric(gsub("\\[(.+),.*\\]", "\\1", confint)),
      CI_upper = as.numeric(gsub("\\[.+, (.+)\\]", "\\1", confint)),
      Outcome  = outcome_name,
      Profile  = .profile_factor(.apply_labels(group, profile_labels), profile_labels)
    ) %>%
    select(Outcome, Profile, group, est, CI_lower, CI_upper, se, pval)
}

#' Bar chart of BCH means (with 95% CIs) for one or more outcomes
#'
#' @param bch_list        List of run_bch() outputs.
#' @param outcome_names   Optional names; defaults to the names stored by run_bch().
#' @param style           "colour" (default) or "pattern" (black & white, needs ggpattern).
plot_bch_outcomes <- function(bch_list,
                              outcome_names  = NULL,
                              profile_labels = NULL,
                              style          = c("colour", "pattern"),
                              title          = "BCH-estimated means",
                              y_label        = "BCH-estimated mean") {
  style <- match.arg(style)
  
  if (is.null(outcome_names)) {
    outcome_names <- vapply(bch_list, function(x) x$name, character(1))
  }
  
  combined <- bind_rows(Map(extract_bch_means, bch_list, outcome_names,
                            MoreArgs = list(profile_labels = profile_labels))) %>%
    mutate(Outcome = factor(Outcome, levels = outcome_names))
  
  dodge <- position_dodge(width = 0.75)
  
  if (style == "pattern") {
    if (!requireNamespace("ggpattern", quietly = TRUE)) {
      stop("style = 'pattern' needs the ggpattern package.")
    }
    patterns <- rep(c("stripe", "crosshatch", "circle", "wave", "none"),
                    length.out = length(outcome_names))
    p <- ggplot(combined, aes(x = Profile, y = est, pattern = Outcome)) +
      ggpattern::geom_col_pattern(
        position = dodge, width = 0.65, fill = "white", colour = "black",
        pattern_fill = "black", pattern_colour = "black",
        pattern_density = 0.25, pattern_spacing = 0.03
      ) +
      ggpattern::scale_pattern_manual(values = setNames(patterns, outcome_names))
  } else {
    p <- ggplot(combined, aes(x = Profile, y = est, fill = Outcome)) +
      geom_col(position = dodge, width = 0.65)
  }
  
  p +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, group = Outcome),
                  position = dodge, width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(x = NULL, y = y_label, title = title)
}

#' Convenience wrapper: run all BCH steps for one outcome
#'
#' @return list(bch, means, omnibus, pairwise, effect_sizes)
bch_full <- function(profile, outcome, outcome_name, N = length(outcome),
                     profile_labels = NULL, verbose = TRUE) {
  b <- run_bch(profile, outcome, outcome_name, profile_labels, verbose)
  list(
    bch          = b,
    means        = b$results,
    omnibus      = bch_omnibus(b, N, outcome_name, verbose),
    pairwise     = bch_pairwise(b, profile_labels, outcome_name, verbose),
    effect_sizes = bch_effect_sizes(b, class_sizes(profile), profile_labels,
                                    outcome_name, verbose)
  )
}
