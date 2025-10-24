#' spm1d R Package - Inference Module
#'
#' Implements statistical inference for SPM objects
#' Adds inference() methods to SPM classes

# Source required modules (in practice, these would be in the package)
# source("spm_classes.R")
# source("rft.R")
# source("geometry.R")

# ==============================================================================
# INFERENCE METHODS FOR SPM OBJECTS
# ==============================================================================

#' Perform statistical inference on SPM object
#' @param spm SPM object (0D or 1D)
#' @param alpha significance level (default 0.05)
#' @param two_tailed logical, perform two-tailed test (for T statistics)
#' @param withBonf Logical, use Bonferroni correction
#' @param cluster_threshold Logical
#' @export
inference <- function(spm, alpha = 0.05, two_tailed = TRUE, withBonf = TRUE,
                      cluster_threshold = NULL) {
  UseMethod("inference")
}

# ==============================================================================
# 0D INFERENCE
# ==============================================================================

#' Inference for 0D T statistic
#' @export
#' @noRd
inference.SPM0D_T <- function(spm, alpha = 0.05, two_tailed = T, withBonf = T, cluster_threshold = F) {

  df <- spm$df[2]  # degrees of freedom for T

  # Calculate threshold
  if (two_tailed) {
    zstar <- t_isf0d(alpha / 2, df)
    # Two-tailed p-value
    p <- 2 * t_sf0d(abs(spm$z), df)
  } else {
    zstar <- t_isf0d(alpha, df)
    # One-tailed p-value
    p <- t_sf0d(spm$z, df)
  }

  # Create inference object
  result <- SPM0Di(spm, alpha, zstar, p)
  result$two_tailed <- two_tailed

  return(result)
}

#' Inference for 0D F statistic
#' @export
#' @noRd
inference.SPM0D_F <- function(spm, alpha = 0.05, two_tailed = T, withBonf = T, cluster_threshold = F) {

  df <- spm$df  # c(df_num, df_den)

  # Calculate threshold
  zstar <- f_isf0d(alpha, df)

  # Calculate p-value
  p <- f_sf0d(spm$z, df)

  # Create inference object
  result <- SPM0Di(spm, alpha, zstar, p)

  return(result)
}

#' Inference for 0D Chi-squared statistic
#' @export
#' @noRd
inference.SPM0D_X2 <- function(spm, alpha = 0.05, two_tailed = T, withBonf = T, cluster_threshold = F) {

  df <- spm$df

  # Calculate threshold
  zstar <- chi2_isf0d(alpha, df)

  # Calculate p-value
  p <- chi2_sf0d(spm$z, df)

  # Create inference object
  result <- SPM0Di(spm, alpha, zstar, p)

  return(result)
}

#' Inference for 0D T^2 statistic
#' @export
#' @noRd
inference.SPM0D_T2 <- function(spm, alpha = 0.05, two_tailed = T, withBonf = T, cluster_threshold = F) {

  df <- spm$df  # c(p, n)

  # Calculate threshold
  zstar <- t2_isf0d(alpha, df)

  # Calculate p-value
  p <- t2_sf0d(spm$z, df)

  # Create inference object
  result <- SPM0Di(spm, alpha, zstar, p)

  return(result)
}

# ==============================================================================
# 1D INFERENCE
# ==============================================================================

#' Inference for 1D T statistic
#' @export
#' @noRd
inference.SPM_T <- function(spm, alpha = 0.05, two_tailed = T, withBonf = T, cluster_threshold = F) {

  df <- spm$df
  fwhm <- spm$fwhm
  resels <- spm$resels

  # Calculate critical threshold using RFT
  # Python spm1d divides alpha by 2 for two-tailed tests
  # and uses Bonferroni correction by default
  if (two_tailed) {
    a <- alpha / 2
  } else {
    a <- alpha
  }

  zstar <- t_isf(a, df, fwhm, resels, withBonf = withBonf)

  # Find clusters
  if (two_tailed) {
    # For two-tailed, look at absolute values
    clusters <- find_clusters(abs(spm$z), zstar)
  } else {
    clusters <- find_clusters(spm$z, zstar)
  }

  # Calculate set-level p-value
  p_set <- t_sf(spm$z, df, fwhm, resels)
  if (two_tailed) {
    p_set <- min(1, 2 * p_set)
  }

  # Calculate cluster-level p-values
  p_cluster <- NULL
  if (nrow(clusters) > 0) {
    p_cluster <- cluster_pvalues(clusters, fwhm, resels)
    clusters$p <- p_cluster
  }

  # Create inference object
  result <- SPMi(spm, alpha, zstar, clusters, p_set, p_cluster)
  result$two_tailed <- two_tailed

  return(result)
}

#' Inference for 1D F statistic
#' @export
#' @noRd
inference.SPM_F <- function(spm, alpha = 0.05, two_tailed = T, withBonf = T, cluster_threshold = F) {

  df <- spm$df
  fwhm <- spm$fwhm
  resels <- spm$resels

  # Calculate critical threshold using RFT with Bonferroni
  zstar <- f_isf(alpha, df, fwhm, resels, withBonf = withBonf)

  # Find clusters
  clusters <- find_clusters(spm$z, zstar)

  # Calculate set-level p-value
  p_set <- f_sf(spm$z, df, fwhm, resels)

  # Calculate cluster-level p-values
  p_cluster <- NULL
  if (nrow(clusters) > 0) {
    p_cluster <- cluster_pvalues(clusters, fwhm, resels)
    clusters$p <- p_cluster
  }

  # Create inference object
  result <- SPMi(spm, alpha, zstar, clusters, p_set, p_cluster)

  return(result)
}

#' Inference for 1D Chi-squared statistic
#' @export
#' @noRd
inference.SPM_X2 <- function(spm, alpha = 0.05, two_tailed = T, withBonf = T, cluster_threshold = F) {

  df <- spm$df
  fwhm <- spm$fwhm
  resels <- spm$resels

  # Calculate critical threshold using RFT
  zstar <- chi2_isf(alpha, df, fwhm, resels)

  # Find clusters
  clusters <- find_clusters(spm$z, zstar)

  # Calculate set-level p-value
  p_set <- chi2_sf(spm$z, df, fwhm, resels)

  # Calculate cluster-level p-values
  p_cluster <- NULL
  if (nrow(clusters) > 0) {
    p_cluster <- cluster_pvalues(clusters, fwhm, resels)
    clusters$p <- p_cluster
  }

  # Create inference object
  result <- SPMi(spm, alpha, zstar, clusters, p_set, p_cluster)

  return(result)
}

#' Inference for 1D T^2 statistic
#' @export
#' @noRd
inference.SPM_T2 <- function(spm, alpha = 0.05, two_tailed = T, withBonf = T, cluster_threshold =  F) {

  df <- spm$df
  fwhm <- spm$fwhm
  resels <- spm$resels

  # Calculate critical threshold using RFT
  zstar <- t2_isf(alpha, df, fwhm, resels)

  # Find clusters
  clusters <- find_clusters(spm$z, zstar)

  # Calculate set-level p-value
  p_set <- t2_sf(spm$z, df, fwhm, resels)

  # Calculate cluster-level p-values
  p_cluster <- NULL
  if (nrow(clusters) > 0) {
    p_cluster <- cluster_pvalues(clusters, fwhm, resels)
    clusters$p <- p_cluster
  }

  # Create inference object
  result <- SPMi(spm, alpha, zstar, clusters, p_set, p_cluster)

  return(result)
}

# ==============================================================================
# PLOTTING METHODS
# ==============================================================================

#' Plot SPM inference results
#' @param x SPMi object (inference result)
#' @param xlab String x axis label
#' @param ylab String y axis label
#' @param main String Main plot title
#' @param plot_threshold Logical
#' @param plot_clusters Logical
#' @param col String
#' @param col_threshold String
#' @param col_cluster String
#' @param lwd Integer Line width
#' @param ... additional plotting parameters
#' @importFrom graphics abline polygon grid legend
#' @importFrom grDevices adjustcolor
#' @export
plot.SPMi <- function(x,
                      xlab = "Time / Position",
                      ylab = NULL,
                      main = NULL,
                      plot_threshold = TRUE,
                      plot_clusters = TRUE,
                      col = "black",
                      col_threshold = "red",
                      col_cluster = "lightblue",
                      lwd = 2,
                      ...) {

  if (is.null(ylab)) {
    ylab <- x$STAT
  }

  if (is.null(main)) {
    main <- sprintf("SPM{%s} Inference (alpha = %.3f)", x$STAT, x$alpha)
  }

  Q <- length(x$z)
  time_points <- 1:Q

  # Set up plot limits
  ylim <- range(c(x$z, x$zstar, -x$zstar), na.rm = TRUE)
  ylim <- ylim + c(-0.1, 0.1) * diff(ylim)

  # Create base plot
  plot(time_points, x$z, type = "l",
       xlab = xlab, ylab = ylab, main = main,
       ylim = ylim, col = col, lwd = lwd, ...)

  # Add threshold line(s)
  if (plot_threshold) {
    abline(h = x$zstar, col = col_threshold, lty = 2, lwd = 1.5)
    if (!is.null(x$two_tailed) && x$two_tailed) {
      abline(h = -x$zstar, col = col_threshold, lty = 2, lwd = 1.5)
    }
  }

  # Highlight clusters
  if (plot_clusters && x$h0reject && nrow(x$clusters) > 0) {
    for (i in 1:nrow(x$clusters)) {
      cluster_range <- x$clusters$start[i]:x$clusters$end[i]
      polygon(c(cluster_range, rev(cluster_range)),
              c(rep(ylim[1], length(cluster_range)),
                rev(x$z[cluster_range])),
              col = adjustcolor(col_cluster, alpha.f = 0.3),
              border = NA)
    }
  }

  # Add grid
  grid(col = "gray90")

  # Add legend
  legend_text <- c("Statistic")
  legend_col <- c(col)
  legend_lty <- c(1)
  legend_lwd <- c(lwd)

  if (plot_threshold) {
    legend_text <- c(legend_text, sprintf("Threshold (alpha=%.3f)", x$alpha))
    legend_col <- c(legend_col, col_threshold)
    legend_lty <- c(legend_lty, 2)
    legend_lwd <- c(legend_lwd, 1.5)
  }

  legend("topright", legend = legend_text,
         col = legend_col, lty = legend_lty, lwd = legend_lwd,
         bg = "white")

  invisible(x)
}

#' Plot SPM object (before inference)
#' @param x SPM object
#' @param xlab String x axis label
#' @param ylab String y axis label
#' @param main String Main plot title
#' @param col String
#' @param lwd Integer Line width
#' @param ... additional plotting parameters
#' @export
#' @importFrom graphics grid
plot.SPM <- function(x,
                     xlab = "Time / Position",
                     ylab = NULL,
                     main = NULL,
                     col = "black",
                     lwd = 2,
                     ...) {

  if (is.null(ylab)) {
    ylab <- x$STAT
  }

  if (is.null(main)) {
    main <- sprintf("SPM{%s}", x$STAT)
  }

  Q <- length(x$z)
  time_points <- 1:Q

  plot(time_points, x$z, type = "l",
       xlab = xlab, ylab = ylab, main = main,
       col = col, lwd = lwd, ...)

  grid(col = "gray90")

  invisible(x)
}

# ==============================================================================
# SUMMARY METHODS FOR INFERENCE OBJECTS
# ==============================================================================

#' @export
summary.SPM0Di <- function(object, ...) {
  cat(sprintf("\nSPM{%s} 0D Inference Results\n", object$STAT))
  cat(sprintf("================================\n"))
  cat(sprintf("Significance level: alpha = %.3f\n", object$alpha))
  cat(sprintf("Test statistic: %s = %.4f\n", object$STAT, object$z))
  cat(sprintf("Critical value: %s* = %.4f\n", object$STAT, object$zstar))
  cat(sprintf("P-value: %s\n", p2string(object$p)))
  cat(sprintf("Null hypothesis rejected: %s\n",
              ifelse(object$h, "YES", "NO")))
  if (!is.null(object$two_tailed) && object$two_tailed) {
    cat("(Two-tailed test)\n")
  }
  invisible(object)
}

#' @export
summary.SPMi <- function(object, ...) {
  cat(sprintf("\nSPM{%s} 1D Inference Results\n", object$STAT))
  cat(sprintf("================================\n"))
  cat(sprintf("Significance level: alpha = %.3f\n", object$alpha))
  cat(sprintf("Critical threshold: %s* = %.4f\n", object$STAT, object$zstar))
  cat(sprintf("FWHM: %.3f\n", object$fwhm))
  cat(sprintf("Set-level p-value: %s\n", p2string(object$p_set)))

  if (object$h0reject) {
    cat(sprintf("\nSuprathreshold clusters: %d\n", nrow(object$clusters)))
    cat("----------------------------\n")
    print(object$clusters)
  } else {
    cat("\nNo suprathreshold clusters detected.\n")
  }

  if (!is.null(object$two_tailed) && object$two_tailed) {
    cat("\n(Two-tailed test)\n")
  }

  invisible(object)
}
