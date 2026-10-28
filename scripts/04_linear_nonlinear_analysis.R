# ---- Linear and non-linear peak classification (SNITCH-based) ----
# Classifies each consensus peak's H3K27ac signal trajectory across age as
# linear-increasing/decreasing (LI/LD) or non-linear-increasing/decreasing
# (NLI/NLD), by comparing a linear model and a GAM fit per peak. Run
# separately on a resected-only sample set and a combined resected-and-
# post-mortem sample set; peaks classified in both become the epigenetic
# ageing signature (EpiAgeSig).

library(data.table)
library(dplyr)
library(mgcv)

prepare_data <- function(data, age) {
  scaled_data <- log2(data + 1)
  list(dat_scaled = scaled_data, Age = age, ages_grid = seq(min(age), max(age)))
}

classify_peak <- function(beta_values, age, cpg_name, covariates = NULL) {

  data <- data.frame(beta = as.numeric(beta_values), Age = age)
  if (!is.null(covariates)) data <- cbind(data, as.data.frame(covariates))
  data <- data[stats::complete.cases(data), , drop = FALSE]
  if (nrow(data) < 5) return(list(CpG = cpg_name, lm_pval = NA, lm_coef = NA, gam_pval = NA, gam_slope = NA))

  covar_terms <- if (!is.null(covariates)) colnames(covariates) else character(0)

  # linear model: beta ~ Age + covariates
  lm_formula <- reformulate(c("Age", covar_terms), response = "beta")
  lm_model <- stats::lm(lm_formula, data = data)
  sm <- summary(lm_model)
  lm_pval <- if ("Age" %in% rownames(sm$coefficients)) sm$coefficients["Age", 4] else NA
  lm_coef <- if ("Age" %in% rownames(sm$coefficients)) sm$coefficients["Age", 1] else NA

  # GAM: beta ~ s(Age, k) + covariates
  k_use <- max(3, min(5, length(unique(data$Age)) - 1))
  gam_formula <- stats::as.formula(paste("beta ~", paste(c(sprintf("s(Age, k=%d)", k_use), covar_terms), collapse = " + ")))
  gam_model <- mgcv::gam(gam_formula, data = data, method = "REML")
  gsum <- summary(gam_model)
  gam_pval <- if (!is.null(gsum$s.table)) gsum$s.table[1, "p-value"] else NA

  # direction/magnitude of non-linear change: GAM-predicted signal at the
  # oldest age minus the youngest age, holding covariates at their median/
  # reference level
  newdata_endpoint <- data.frame(Age = c(min(data$Age), max(data$Age)))
  if (!is.null(covariates)) {
    for (covar in colnames(covariates)) {
      v <- data[[covar]]
      newdata_endpoint[[covar]] <- if (is.factor(v) || is.character(v)) {
        factor(rep(levels(as.factor(v))[1], 2), levels = levels(as.factor(v)))
      } else rep(stats::median(v, na.rm = TRUE), 2)
    }
  }
  endpoint_pred <- stats::predict(gam_model, newdata = newdata_endpoint)
  gam_slope <- as.numeric(endpoint_pred[2] - endpoint_pred[1])

  list(CpG = cpg_name, lm_pval = lm_pval, lm_coef = lm_coef, gam_pval = gam_pval, gam_slope = gam_slope)
}

run_classification <- function(dat_scaled, age, covariates = NULL) {

  results <- lapply(seq_len(ncol(dat_scaled)), function(i) {
    classify_peak(dat_scaled[, i], age, colnames(dat_scaled)[i], covariates)
  })

  results_df <- do.call(rbind, lapply(results, as.data.frame))

  # classification is based on nominal (uncorrected) p-values
  results_df$class <- NA
  results_df[results_df$lm_pval <= 0.05 & results_df$lm_coef > 0, "class"] <- "LI"
  results_df[results_df$lm_pval <= 0.05 & results_df$lm_coef < 0, "class"] <- "LD"
  results_df[results_df$gam_pval < 0.05 & results_df$gam_slope > 0, "class"] <- "NLI"
  results_df[results_df$gam_pval < 0.05 & results_df$gam_slope < 0, "class"] <- "NLD"

  results_df[!is.na(results_df$class), ]
}

# ---- Run separately on each sample set ----

sample_sets <- c("resected_only", "resected_and_postmortem")
classified_all <- list()

for (set_name in sample_sets) {

  dgList <- readRDS(paste0("/path/to/", set_name, "_dgList.rds"))

  covariates <- dgList$samples[, c("sex", "library_structure", "brain_bank"), drop = FALSE]
  covariates[] <- lapply(covariates, as.factor)

  scaled <- prepare_data(data = t(dgList$counts), age = dgList$samples$age)

  classified_all[[set_name]] <- run_classification(
    dat_scaled = scaled$dat_scaled,
    age = scaled$Age,
    covariates = covariates
  )
}

# ---- EpiAgeSig: peaks classified (LI/LD/NLI/NLD) in both sample sets ----

common_peaks <- Reduce(intersect, lapply(classified_all, function(df) df$CpG))
