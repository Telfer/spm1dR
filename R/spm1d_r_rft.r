#' spm1d R Package - Random Field Theory (RFT) Module
#'
#' Implements Random Field Theory calculations for 1D statistical inference
#' Based on Worsley et al. (1996) and Friston et al. (1994)

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

eps <- .Machine$double.eps

#' Euler characteristic for different distributions
#' @param u threshold value
#' @param df degrees of freedom
#' @param stat statistic type
euler_characteristic <- function(u, df, stat = "T") {
  if (stat == "T") {
    return(euler_characteristic_T(u, df))
  } else if (stat == "F") {
    return(euler_characteristic_F(u, df))
  } else if (stat == "X2") {
    return(euler_characteristic_X2(u, df))
  } else if (stat == "T2") {
    return(euler_characteristic_T2(u, df))
  }
}

#' Euler characteristic for T field
#' @noRd
euler_characteristic_T <- function(t, df) {
  # For 1D T field, EC density (per resel)
  # From Worsley et al. (1996) equation for T fields
  # ρ(t) = (4 log 2)^(1/2) / (2π) * (1 + t²/ν)^(-(ν-1)/2)
  #
  # Note: Some formulations omit the (4 log 2)^(1/2) factor
  # Python spm1d uses the simpler form without this factor

  if (t <= 0) return(0)

  # Standard EC density for T field (matches SPM/spm1d)
  rho <- (2 * pi)^(-0.5) * (1 + t^2 / df)^(-(df - 1) / 2)

  return(rho)
}

#' Euler characteristic for F field
#' @noRd
euler_characteristic_F <- function(f, df) {
  # df is a vector: c(df_num, df_den)
  # Match Python spm1d implementation exactly

  if (f < 0) return(c(1, Inf))  # To bypass warnings

  df_num <- df[1]
  df_den <- df[2]

  # Ensure df_num >= 1 to avoid NaN
  df_num <- max(df_num, 1.0)

  # Python formula from prob.py:
  # a = FOUR_LOG2/TWO_PI
  # EC[1] = a**0.5 * exp(gammaln((v+k-1)/2)-b) * 2**0.5 * (k*z/v)**(0.5*(k-1)) * (1+k*z/v)**(-0.5*(v+k-2))
  # where b = gammaln(v/2) + gammaln(k/2)

  FOUR_LOG2 <- 4 * log(2)
  TWO_PI <- 2 * pi

  a <- FOUR_LOG2 / TWO_PI
  b <- lgamma(df_den/2) + lgamma(df_num/2)

  ec0d <- 1 - pf(f, df_num, df_den)

  ec1d <- sqrt(a) * exp(lgamma((df_den + df_num - 1)/2) - b) * sqrt(2) *
          (df_num * f / df_den)^(0.5 * (df_num - 1)) *
          (1 + df_num * f / df_den)^(-0.5 * (df_den + df_num - 2))

  return(ec1d)
}

#' Euler characteristic for Chi-squared field
#' @noRd
euler_characteristic_X2 <- function(x2, df) {
  (2 * pi)^(-0.5) * (x2 / 2)^((df - 1) / 2) * exp(-x2 / 2)
}

#' Euler characteristic for T² (Hotelling's) field
#' @noRd
euler_characteristic_T2 <- function(t2, df) {
  # df is a vector: c(p, n)
  # p = number of variables, n = sample size
  p <- df[1]
  n <- df[2]

  (2 * pi)^(-0.5) * (t2)^((p - 1) / 2) *
    (1 + t2 / (n - p))^(-(n - 1) / 2)
}

# ==============================================================================
# T DISTRIBUTION RFT FUNCTIONS
# ==============================================================================

#' T distribution: survival function (0D)
#' @param t test statistic
#' @param df degrees of freedom
#' @importFrom stats pt
#' @noRd
t_sf0d <- function(t, df) {
  pt(t, df, lower.tail = FALSE)
}

#' T distribution: inverse survival function (0D)
#' @param alpha significance level
#' @param df degrees of freedom
#' @importFrom stats qt
#' @noRd
t_isf0d <- function(alpha, df) {
  qt(alpha, df, lower.tail = FALSE)
}

#' T distribution: survival function (1D)
#' @param t test statistic continuum
#' @param df degrees of freedom
#' @param fwhm field smoothness
#' @param resels resel counts
#' @noRd
t_sf <- function(t, df, fwhm, resels) {
  # 1D Random Field Theory survival function for T field
  # P(max(T) > t) using Euler characteristic heuristic
  # From: Friston et al. (1994), Worsley et al. (1996)

  # Get maximum of the field (absolute value for proper tail probability)
  t_max <- max(abs(t), na.rm = TRUE)

  # Expected Euler characteristic at this threshold
  resel_count <- resels[1]
  df_val <- df[2]

  # EC density
  rho <- euler_characteristic_T(t_max, df_val)

  # Expected EC = resel count * EC density
  # P(max > threshold) ≈ E[EC(threshold)]
  p <- resel_count * rho

  # Ensure p is between 0 and 1
  p <- max(0, min(p, 1))

  return(p)
}

#' T distribution: inverse survival function (1D)
#' @param alpha significance level
#' @param df degrees of freedom
#' @param fwhm field smoothness
#' @param resels resel counts
#' @param withBonf Logical
#' @importFrom stats pt qt
#' @noRd
t_isf <- function(alpha, df, fwhm, resels, withBonf = TRUE) {
  # Find threshold t such that P(max(T) > t) = alpha
  # Python spm1d uses Bonferroni correction by default

  resel_count <- resels[1]
  df_val <- df[2]
  Q <- round((resel_count * fwhm) + 1)  # Approximate number of nodes

  # Objective: find t where corrected p-value = alpha
  objective <- function(t) {
    if (t <= 0) return(1e10)

    # RFT p-value
    rho <- euler_characteristic_T(t, df_val)
    p_rft <- resel_count * rho

    # Bonferroni p-value (if enabled)
    if (withBonf && !is.null(Q)) {
      p_bonf <- Q * pt(t, df_val, lower.tail = FALSE)
      p_final <- min(p_rft, p_bonf)
    } else {
      p_final <- p_rft
    }

    return((p_final - alpha)^2)
  }

  # Initial guess from univariate
  t0 <- qt(alpha, df_val, lower.tail = FALSE)

  # Search bounds
  lower_bound <- t0 * 0.5
  upper_bound <- t0 * 3

  # Optimize
  result <- optimize(
    objective,
    interval = c(lower_bound, upper_bound),
    tol = 1e-9
  )

  return(result$minimum)
}

# ==============================================================================
# F DISTRIBUTION RFT FUNCTIONS
# ==============================================================================

#' F distribution: survival function (0D)
#' @importFrom stats pf
#' @noRd
f_sf0d <- function(f, df) {
  pf(f, df[1], df[2], lower.tail = FALSE)
}

#' F distribution: inverse survival function (0D)
#' @importFrom stats qf
#' @noRd
f_isf0d <- function(alpha, df) {
  qf(alpha, df[1], df[2], lower.tail = FALSE)
}

#' F distribution: survival function (1D)
#' @noRd
f_sf <- function(f, df, fwhm, resels) {
  f_max <- max(f, na.rm = TRUE)
  resel_count <- resels[1]
  ec <- euler_characteristic_F(f_max, df)
  p <- resel_count * ec
  return(min(max(p, 0), 1))
}

#' F distribution: inverse survival function (1D) with Bonferroni correction
#' @param alpha significance level
#' @param df degrees of freedom (vector: c(df_num, df_den))
#' @param fwhm field smoothness
#' @param resels resel counts
#' @param withBonf use Bonferroni correction if less severe (default TRUE)
#' @importFrom stats pf qf optimize
f_isf <- function(alpha, df, fwhm, resels, withBonf = TRUE) {
  # Find threshold f such that P(max(F) > f) = alpha

  resel_count <- resels[1]
  df_num <- df[1]
  df_den <- df[2]
  Q <- round((resel_count * fwhm) + 1)

  # Objective: find f where corrected p-value = alpha
  objective <- function(f) {
    if (f <= 0) return(1e10)

    # RFT p-value
    rho <- euler_characteristic_F(f, df)
    p_rft <- resel_count * rho

    # Bonferroni p-value (if enabled)
    if (withBonf && !is.null(Q)) {
      p_bonf <- Q * pf(f, df_num, df_den, lower.tail = FALSE)
      p_final <- min(p_rft, p_bonf)
    } else {
      p_final <- p_rft
    }

    return((p_final - alpha)^2)
  }

  # Initial guess from univariate
  f0 <- qf(alpha, df_num, df_den, lower.tail = FALSE)

  # Search bounds
  lower_bound <- f0 * 0.5
  upper_bound <- f0 * 5

  # Optimize
  result <- optimize(
    objective,
    interval = c(lower_bound, upper_bound),
    tol = 1e-9
  )

  return(result$minimum)
}

# ==============================================================================
# CHI-SQUARED DISTRIBUTION RFT FUNCTIONS
# ==============================================================================

#' Chi-squared distribution: survival function (0D)
#' @importFrom stats pchisq
#' @noRd
chi2_sf0d <- function(x2, df) {
  pchisq(x2, df, lower.tail = FALSE)
}

#' Chi-squared distribution: inverse survival function (0D)
#' @importFrom stats qchisq
#' @noRd
chi2_isf0d <- function(alpha, df) {
  qchisq(alpha, df, lower.tail = FALSE)
}

#' Chi-squared distribution: survival function (1D)
#' @noRd
chi2_sf <- function(x2, df, fwhm, resels) {
  x2_max <- max(x2, na.rm = TRUE)
  resel_count <- resels[1]
  ec <- euler_characteristic_X2(x2_max, df)
  p <- resel_count * ec
  return(min(max(p, 0), 1))
}

#' Chi-squared distribution: inverse survival function (1D)
#' @importFrom stats optim
#' @noRd
chi2_isf <- function(alpha, df, fwhm, resels) {
  resel_count <- resels[1]

  obj_fn <- function(x2) {
    ec <- euler_characteristic_X2(x2, df)
    abs(resel_count * ec - alpha)
  }

  x2_0 <- chi2_isf0d(alpha, df)
  result <- optim(x2_0, obj_fn, method = "Brent",
                  lower = 0, upper = 10 * x2_0)

  return(result$par)
}

# ==============================================================================
# T² (HOTELLING'S) DISTRIBUTION RFT FUNCTIONS
# ==============================================================================

#' T² distribution: survival function (0D)
#' Hotelling's T² follows an F distribution after transformation
#' importFrom stats pf
#' @noRd
t2_sf0d <- function(t2, df) {
  # df = c(p, n) where p = # variables, n = sample size
  # T² ~ (n-1)p/(n-p) * F(p, n-p)
  p <- df[1]
  n <- df[2]
  f_stat <- t2 * (n - p) / ((n - 1) * p)
  pf(f_stat, p, n - p, lower.tail = FALSE)
}

#' T² distribution: inverse survival function (0D)
#' @importFrom stats qf
#' @noRd
t2_isf0d <- function(alpha, df) {
  p <- df[1]
  n <- df[2]
  f_crit <- qf(alpha, p, n - p, lower.tail = FALSE)
  t2_crit <- f_crit * (n - 1) * p / (n - p)
  return(t2_crit)
}

#' T² distribution: survival function (1D)
#' @noRd
t2_sf <- function(t2, df, fwhm, resels) {
  t2_max <- max(t2, na.rm = TRUE)
  resel_count <- resels[1]
  ec <- euler_characteristic_T2(t2_max, df)
  p <- resel_count * ec
  return(min(max(p, 0), 1))
}

#' T² distribution: inverse survival function (1D)
#' @noRd
t2_isf <- function(alpha, df, fwhm, resels) {
  resel_count <- resels[1]

  obj_fn <- function(t2) {
    ec <- euler_characteristic_T2(t2, df)
    abs(resel_count * ec - alpha)
  }

  t2_0 <- t2_isf0d(alpha, df)
  result <- optim(t2_0, obj_fn, method = "Brent",
                  lower = 0, upper = 10 * t2_0)

  return(result$par)
}

# ==============================================================================
# CLUSTER-LEVEL INFERENCE
# ==============================================================================

#' Find clusters in a 1D continuum that exceed threshold
#' @param z statistic continuum
#' @param zstar threshold
#' @return data frame with cluster information
find_clusters <- function(z, zstar) {
  # Identify suprathreshold regions
  above_threshold <- z > zstar

  # Handle NAs (from ROI masking)
  above_threshold[is.na(above_threshold)] <- FALSE

  if (!any(above_threshold)) {
    return(data.frame(
      cluster = integer(0),
      start = integer(0),
      end = integer(0),
      extent = integer(0),
      peak_height = numeric(0),
      peak_location = integer(0)
    ))
  }

  # Find cluster starts and ends
  diff_thresh <- diff(c(0, as.integer(above_threshold), 0))
  starts <- which(diff_thresh == 1)
  ends <- which(diff_thresh == -1) - 1

  # Create cluster data frame
  n_clusters <- length(starts)
  clusters <- data.frame(
    cluster = 1:n_clusters,
    start = starts,
    end = ends,
    extent = ends - starts + 1,
    peak_height = numeric(n_clusters),
    peak_location = integer(n_clusters)
  )

  # Find peak in each cluster
  for (i in 1:n_clusters) {
    cluster_region <- starts[i]:ends[i]
    cluster_values <- z[cluster_region]
    peak_idx <- which.max(cluster_values)
    clusters$peak_height[i] <- cluster_values[peak_idx]
    clusters$peak_location[i] <- cluster_region[peak_idx]
  }

  return(clusters)
}

#' Calculate cluster p-values
#' @param clusters cluster data frame
#' @param fwhm field smoothness
#' @param resels resel counts
cluster_pvalues <- function(clusters, fwhm, resels) {
  # Simple approximation for cluster p-values
  # More sophisticated methods would use cluster extent probabilities

  if (nrow(clusters) == 0) {
    return(numeric(0))
  }

  # Bonferroni correction as conservative estimate
  n_clusters <- nrow(clusters)
  p_values <- rep(1.0, n_clusters)

  # Could implement more sophisticated cluster extent probabilities here
  # For now, return conservative values

  return(p_values)
}
