#' spm1d R Package - Core SPM Classes
#'
#' One-Dimensional Statistical Parametric Mapping in R
#' Port of Python spm1d package by Todd Pataky
#'
#' This module implements the core SPM object classes

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Machine epsilon for numerical stability
eps <- .Machine$double.eps

#' Convert degrees of freedom to string
#' @param df degrees of freedom (vector)
dflist2str <- function(df) {
  paste0("(", paste(df, collapse = ", "), ")")
}

#' Convert p-value to formatted string
#' @param p p-value
p2string <- function(p) {
  if (p < 0.0005) {
    return("p < 0.001")
  } else {
    return(sprintf("p = %.3f", p))
  }
}

# ==============================================================================
# BASE SPM CLASSES
# ==============================================================================

#' Base SPM Parent Class
#'
#' @description
#' Parent class for all SPM objects (both 0D and 1D)
SPMParent <- function() {
  structure(
    list(
      STAT = NULL,
      z = NULL,
      df = NULL,
      dim = NULL
    ),
    class = "SPMParent"
  )
}

# ==============================================================================
# 0D SPM CLASSES (scalar statistics)
# ==============================================================================

#' Base 0D SPM Class
#' @param STAT statistic type ('T', 'F', 'X2', 'T2')
#' @param z test statistic value
#' @param df degrees of freedom
#' @param beta regression coefficients (optional)
#' @param residuals residuals (optional)
#' @param sigma2 variance (optional)
SPM0D <- function(STAT, z, df, beta = NULL, residuals = NULL, sigma2 = NULL) {
  obj <- list(
    STAT = STAT,
    z = as.numeric(z),
    df = df,
    dim = 0,
    beta = beta,
    residuals = residuals,
    sigma2 = sigma2,
    effect = paste("SPM{", STAT, "}", sep = ""),
    effect_short = STAT
  )
  class(obj) <- c(paste0("SPM0D_", STAT), "SPM0D", "SPMParent")
  return(obj)
}

#' SPM{t} 0D Class
#' @noRd
SPM0D_T <- function(z, df, beta = NULL, residuals = NULL, sigma2 = NULL) {
  SPM0D("T", z, df, beta, residuals, sigma2)
}

#' SPM{F} 0D Class
#' @noRd
SPM0D_F <- function(z, df, X = NULL, beta = NULL, residuals = NULL, sigma2 = NULL) {
  obj <- SPM0D("F", z, df, beta, residuals, sigma2)
  obj$X <- X
  return(obj)
}

#' SPM{T2} 0D Class (Hotelling's T-squared)
#' @noRd
SPM0D_T2 <- function(z, df, residuals = NULL) {
  SPM0D("T2", z, df, residuals = residuals)
}

#' SPM{X2} 0D Class (Chi-squared)
#' @noRd
SPM0D_X2 <- function(z, df, residuals = NULL) {
  SPM0D("X2", z, df, residuals = residuals)
}

# ==============================================================================
# 1D SPM CLASSES (continuum statistics)
# ==============================================================================

#' Base 1D SPM Class
#' @param STAT statistic type ('T', 'F', 'X2', 'T2')
#' @param z test statistic continuum (vector)
#' @param df degrees of freedom
#' @param fwhm full-width at half-maximum (field smoothness)
#' @param resels resel counts
#' @param X design matrix (optional)
#' @param beta fitted parameters (optional)
#' @param residuals residuals (optional)
#' @param sigma2 variance (optional)
#' @param roi region of interest mask (optional)
SPM <- function(STAT, z, df, fwhm, resels, X = NULL, beta = NULL,
                residuals = NULL, sigma2 = NULL, roi = NULL) {

  # Handle masked arrays (ROI)
  if (!is.null(roi)) {
    z <- ifelse(roi, z, NA)
  }

  obj <- list(
    STAT = STAT,
    z = as.numeric(z),
    df = df,
    dim = 1,
    Q = length(z),  # number of nodes (time points)
    fwhm = fwhm,
    resels = resels,
    X = X,
    beta = beta,
    residuals = residuals,
    sigma2 = sigma2,
    roi = roi,
    effect = paste("SPM{", STAT, "}", sep = ""),
    effect_short = STAT,
    ismasked = !is.null(roi)
  )

  class(obj) <- c(paste0("SPM_", STAT), "SPM", "SPMParent")
  return(obj)
}

#' SPM{t} 1D Class
#' @noRd
SPM_T <- function(z, df, fwhm, resels, X = NULL, beta = NULL,
                  residuals = NULL, sigma2 = NULL, roi = NULL) {
  SPM("T", z, df, fwhm, resels, X, beta, residuals, sigma2, roi)
}

#' SPM{F} 1D Class
#' @noRd
SPM_F <- function(z, df, fwhm, resels, X = NULL, beta = NULL,
                  residuals = NULL, X0 = NULL, sigma2 = NULL, roi = NULL) {
  obj <- SPM("F", z, df, fwhm, resels, X, beta, residuals, sigma2, roi)
  obj$X0 <- X0  # reduced design matrix for F-test
  return(obj)
}

#' SPM{T2} 1D Class (Hotelling's T-squared)
#' @noRd
SPM_T2 <- function(z, df, fwhm, resels, X = NULL, beta = NULL,
                   residuals = NULL, sigma2 = NULL, roi = NULL) {
  SPM("T2", z, df, fwhm, resels, X, beta, residuals, sigma2, roi)
}

#' SPM{X2} 1D Class (Chi-squared)
#' @noRd
SPM_X2 <- function(z, df, fwhm, resels, X = NULL, beta = NULL,
                   residuals = NULL, sigma2 = NULL, roi = NULL) {
  SPM("X2", z, df, fwhm, resels, X, beta, residuals, sigma2, roi)
}

# ==============================================================================
# PRINT METHODS
# ==============================================================================

#' @export
print.SPM0D <- function(x, ...) {
  cat(sprintf("SPM{%s} 0D\n", x$STAT))
  cat(sprintf("  statistic: %s = %.3f\n", x$STAT, x$z))
  cat(sprintf("  df: %s\n", dflist2str(x$df)))
  invisible(x)
}

#' @export
print.SPM <- function(x, ...) {
  cat(sprintf("SPM{%s} 1D\n", x$STAT))
  cat(sprintf("  statistic: %s\n", x$STAT))
  cat(sprintf("  nodes: %d\n", x$Q))
  cat(sprintf("  df: %s\n", dflist2str(x$df)))
  cat(sprintf("  FWHM: %.3f\n", x$fwhm))
  if (x$ismasked) {
    cat("  (ROI masked)\n")
  }
  invisible(x)
}

# ==============================================================================
# SUMMARY METHODS
# ==============================================================================

#' @export
summary.SPM0D <- function(object, ...) {
  cat(sprintf("\n%s\n", object$effect))
  cat(sprintf("%s = %.3f, df = %s\n",
              object$STAT, object$z, dflist2str(object$df)))
}

#' @export
summary.SPM <- function(object, ...) {
  cat(sprintf("\n%s (1D)\n", object$effect))
  cat(sprintf("  Nodes (Q): %d\n", object$Q))
  cat(sprintf("  df: %s\n", dflist2str(object$df)))
  cat(sprintf("  FWHM: %.3f\n", object$fwhm))
  cat(sprintf("  %s range: [%.3f, %.3f]\n",
              object$STAT, min(object$z, na.rm = TRUE),
              max(object$z, na.rm = TRUE)))
}

# ==============================================================================
# INFERENCE CLASSES (for storing inference results)
# ==============================================================================

#' SPM Inference Result (0D)
#' @param spm original SPM object
#' @param alpha significance level
#' @param zstar critical threshold
#' @param p p-value
SPM0Di <- function(spm, alpha, zstar, p) {
  obj <- list(
    spm = spm,
    STAT = spm$STAT,
    alpha = alpha,
    zstar = zstar,
    z = spm$z,
    p = p,
    h = p < alpha,  # null hypothesis rejection
    df = spm$df,
    dim = 0
  )
  class(obj) <- c(paste0("SPM0Di_", spm$STAT), "SPM0Di")
  return(obj)
}

#' SPM Inference Result (1D)
#' @param spm original SPM object
#' @param alpha significance level
#' @param zstar critical threshold
#' @param clusters cluster information
#' @param p_set p-value for set-level inference
#' @param p_cluster p-values for cluster-level inference
SPMi <- function(spm, alpha, zstar, clusters, p_set = NULL, p_cluster = NULL) {
  obj <- list(
    spm = spm,
    STAT = spm$STAT,
    alpha = alpha,
    zstar = zstar,
    z = spm$z,
    clusters = clusters,
    p_set = p_set,
    p_cluster = p_cluster,
    h0reject = !is.null(clusters) && nrow(clusters) > 0,
    df = spm$df,
    fwhm = spm$fwhm,
    dim = 1
  )
  class(obj) <- c(paste0("SPMi_", spm$STAT), "SPMi")
  return(obj)
}

# ==============================================================================
# PRINT METHODS FOR INFERENCE OBJECTS
# ==============================================================================

#' @export
print.SPM0Di <- function(x, ...) {
  cat(sprintf("SPM{%s} 0D Inference\n", x$STAT))
  cat(sprintf("  alpha: %.3f\n", x$alpha))
  cat(sprintf("  %s: %.3f\n", x$STAT, x$z))
  cat(sprintf("  %s*: %.3f\n", x$STAT, x$zstar))
  cat(sprintf("  %s\n", p2string(x$p)))
  cat(sprintf("  H0 rejected: %s\n", ifelse(x$h, "TRUE", "FALSE")))
  invisible(x)
}

#' @export
print.SPMi <- function(x, ...) {
  cat(sprintf("SPM{%s} 1D Inference\n", x$STAT))
  cat(sprintf("  alpha: %.3f\n", x$alpha))
  cat(sprintf("  threshold: %s* = %.3f\n", x$STAT, x$zstar))
  cat(sprintf("  FWHM: %.3f\n", x$fwhm))

  if (x$h0reject) {
    cat(sprintf("  Clusters: %d suprathreshold cluster(s)\n",
                nrow(x$clusters)))
    print(x$clusters)
  } else {
    cat("  No suprathreshold clusters\n")
  }

  invisible(x)
}

# ==============================================================================
# EXPORTS
# ==============================================================================

# Export main constructor functions
# These will be called by stats functions (glm, ttest, anova, etc.)
