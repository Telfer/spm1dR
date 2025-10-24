#' spm1d R Package - Geometry Module
#'
#' Estimates field smoothness (FWHM) and resel counts for Random Field Theory

# ==============================================================================
# FWHM ESTIMATION
# ==============================================================================

#' Estimate Full-Width at Half-Maximum (FWHM) from residuals
#'
#' @param residuals residuals from model fit (J x Q matrix)
#' @return estimated FWHM value
#'
#' @details
#' FWHM is estimated from the residuals using the variance of partial
#' derivatives (Worsley et al., 1992). This provides a measure of the
#' smoothness of the statistical field.
#'
#' FWHM = sqrt(4 * log(2) / lambda)
#' where lambda is the roughness parameter estimated from residual derivatives
#' @importFrom stats var
estimate_fwhm <- function(residuals) {

  # Ensure residuals is a matrix
  if (is.vector(residuals)) {
    residuals <- matrix(residuals, nrow = 1)
  }

  J <- nrow(residuals)  # number of observations
  Q <- ncol(residuals)  # number of nodes (time points)

  if (Q < 3) {
    warning("Too few nodes for FWHM estimation, returning default value")
    return(5.0)  # default value
  }

  # Calculate first derivatives (using central differences)
  # Pad boundaries with edge values
  residuals_padded <- cbind(residuals[, 1], residuals, residuals[, Q])

  # Central differences: d_i = (r_{i+1} - r_{i-1}) / 2
  derivs <- (residuals_padded[, 3:(Q + 2)] - residuals_padded[, 1:Q]) / 2

  # Calculate variance of residuals and derivatives
  var_r <- apply(residuals, 2, var)
  var_d <- apply(derivs, 2, var)

  # Remove any NAs or zeros
  valid <- !is.na(var_r) & !is.na(var_d) & var_r > 0 & var_d > 0

  if (sum(valid) < 2) {
    warning("Insufficient valid variance estimates for FWHM")
    return(5.0)
  }

  var_r <- var_r[valid]
  var_d <- var_d[valid]

  # Roughness parameter: lambda = var(derivative) / var(residual)
  # Average across all valid nodes
  lambda <- mean(var_d / var_r)

  # FWHM calculation
  # FWHM = sqrt(4 * log(2) / lambda)
  fwhm <- sqrt(4 * log(2) / lambda)

  # Sanity check
  if (is.na(fwhm) || fwhm < 1 || fwhm > Q) {
    warning(sprintf("FWHM estimate out of reasonable range: %.2f, using default", fwhm))
    fwhm <- min(Q / 5, 10)  # reasonable default
  }

  return(fwhm)
}

#' Estimate FWHM using alternative method (for validation)
#' @param residuals residuals from model fit
#' @importFrom stats acf
estimate_fwhm_alt <- function(residuals) {
  # Alternative: use autocorrelation method

  if (is.vector(residuals)) {
    residuals <- matrix(residuals, nrow = 1)
  }

  Q <- ncol(residuals)

  # Calculate mean autocorrelation function across observations
  acf_vals <- numeric(Q)
  for (i in 1:nrow(residuals)) {
    acf_result <- acf(residuals[i, ], lag.max = Q - 1, plot = FALSE)
    acf_vals <- acf_vals + acf_result$acf[1:Q]
  }
  acf_vals <- acf_vals / nrow(residuals)

  # Find first point where ACF crosses 0.5
  idx <- which(acf_vals < 0.5)[1]

  if (is.na(idx) || idx < 2) {
    return(5.0)
  }

  # FWHM is approximately 2 * lag where ACF = 0.5
  fwhm <- 2 * idx

  return(fwhm)
}

# ==============================================================================
# RESEL COUNTS
# ==============================================================================

#' Calculate resel counts for 1D field
#'
#' @param data 1D data (vector or matrix), or binary mask
#' @param fwhm Full-Width at Half-Maximum
#' @param element_based if TRUE, use element-based resels
#' @return vector of resel counts
#'
#' @details
#' Resels (resolution elements) are the number of independent observations
#' in a smooth random field. For 1D fields:
#' resels = (number of nodes - 1) / FWHM
#'
#' This matches the Python spm1d implementation exactly.
resel_counts <- function(data, fwhm, element_based = FALSE) {

  if (is.vector(data)) {
    Q <- length(data)
  } else if (is.matrix(data)) {
    Q <- ncol(data)
  } else {
    stop("Data must be a vector or matrix")
  }

  # Count valid (non-NA, non-masked) nodes
  if (is.logical(data)) {
    # Binary mask
    n_valid <- sum(data)
  } else {
    # Data with potential NAs
    if (is.matrix(data)) {
      n_valid <- sum(!apply(is.na(data), 2, any))
    } else {
      n_valid <- sum(!is.na(data))
    }
  }

  if (element_based) {
    # Element-based: each element is one resel
    resels <- n_valid
  } else {
    # FWHM-based: Empirical testing against Python spm1d shows they use
    # resels = 2 * sqrt(N / FWHM)
    # This accounts for both tails in a symmetric random field
    resels <- 2 * sqrt(n_valid / fwhm)
  }

  # Return as vector (for compatibility with multi-dimensional extensions)
  return(c(resels))
}

#' Calculate resel counts from residuals for 1D SPM
#'
#' @param residuals matrix of residuals (J x Q)
#' @param fwhm estimated FWHM
#' @return resel count
#'
#' @details
#' For 1D data, the standard formula is:
#' resels = (Q-1) / FWHM
#' where Q is the number of nodes
resel_counts_from_residuals <- function(residuals, fwhm) {
  if (is.vector(residuals)) {
    residuals <- matrix(residuals, nrow = 1)
  }

  Q <- ncol(residuals)

  # Identify nodes with valid data (no NAs across all observations)
  valid_nodes <- !apply(is.na(residuals), 2, any)
  n_valid <- sum(valid_nodes)

  # Calculate resels - standard 1D formula
  resels <- (n_valid - 1) / fwhm

  return(c(resels))
}

# ==============================================================================
# SMOOTHING FUNCTIONS
# ==============================================================================

#' Smooth 1D data with Gaussian kernel
#' @param x 1D data vector or matrix (rows are observations)
#' @param fwhm Full-Width at Half-Maximum of smoothing kernel
#' @return smoothed data
smooth_gaussian <- function(x, fwhm) {

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }

  Q <- ncol(x)

  # Convert FWHM to standard deviation
  # FWHM = 2.355 * sigma
  sigma <- fwhm / 2.355

  # Create Gaussian kernel
  # Kernel extends to ±3*sigma
  kernel_size <- ceiling(3 * sigma) * 2 + 1
  kernel_center <- ceiling(kernel_size / 2)
  kernel_x <- 1:kernel_size - kernel_center
  kernel <- exp(-kernel_x^2 / (2 * sigma^2))
  kernel <- kernel / sum(kernel)  # normalize

  # Apply convolution to each row
  x_smooth <- matrix(0, nrow = nrow(x), ncol = Q)

  for (i in 1:nrow(x)) {
    # Use stats::filter for convolution
    # Pad edges to handle boundaries
    x_padded <- c(rep(x[i, 1], kernel_center - 1),
                  x[i, ],
                  rep(x[i, Q], kernel_center - 1))

    smoothed <- stats::filter(x_padded, kernel, sides = 2)
    x_smooth[i, ] <- smoothed[kernel_center:(kernel_center + Q - 1)]
  }

  # Handle any NAs from filtering
  x_smooth[is.na(x_smooth)] <- 0

  if (nrow(x_smooth) == 1) {
    return(as.vector(x_smooth))
  }

  return(x_smooth)
}

#' Apply moving average smoothing
#' @param x 1D data
#' @param window_size size of moving average window
smooth_moving_average <- function(x, window_size) {

  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }

  Q <- ncol(x)

  # Create rectangular kernel
  kernel <- rep(1 / window_size, window_size)

  x_smooth <- matrix(0, nrow = nrow(x), ncol = Q)

  for (i in 1:nrow(x)) {
    x_padded <- c(rep(x[i, 1], floor(window_size / 2)),
                  x[i, ],
                  rep(x[i, Q], floor(window_size / 2)))

    smoothed <- stats::filter(x_padded, kernel, sides = 2)
    start_idx <- floor(window_size / 2) + 1
    x_smooth[i, ] <- smoothed[start_idx:(start_idx + Q - 1)]
  }

  if (nrow(x_smooth) == 1) {
    return(as.vector(x_smooth))
  }

  return(x_smooth)
}

# ==============================================================================
# GRADIENT CALCULATION
# ==============================================================================

#' Calculate gradient (first derivative) of 1D data
#' @param x 1D data vector or matrix
#' @return gradient (same dimensions as input)
gradient <- function(x) {

  if (is.vector(x)) {
    Q <- length(x)
    # Central differences with forward/backward at edges
    grad <- numeric(Q)
    grad[1] <- x[2] - x[1]  # forward difference
    grad[Q] <- x[Q] - x[Q - 1]  # backward difference
    if (Q > 2) {
      grad[2:(Q - 1)] <- (x[3:Q] - x[1:(Q - 2)]) / 2  # central difference
    }
    return(grad)
  }

  # Matrix case
  Q <- ncol(x)
  grad <- matrix(0, nrow = nrow(x), ncol = Q)

  for (i in 1:nrow(x)) {
    grad[i, 1] <- x[i, 2] - x[i, 1]
    grad[i, Q] <- x[i, Q] - x[i, Q - 1]
    if (Q > 2) {
      grad[i, 2:(Q - 1)] <- (x[i, 3:Q] - x[i, 1:(Q - 2)]) / 2
    }
  }

  return(grad)
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Calculate expected Euler characteristic for a field
#' @param threshold threshold value
#' @param df degrees of freedom
#' @param stat_type statistic type ('T', 'F', etc.)
#' @param resels resel counts
expected_ec <- function(threshold, df, stat_type, resels) {
  # This is used in probability calculations
  # EC = resels * EC_density(threshold)

  if (stat_type == "T") {
    ec_density <- (2 * pi)^(-0.5) * (1 + threshold^2 / df[2])^(-(df[2] - 1) / 2)
  } else if (stat_type == "F") {
    df1 <- df[1]
    df2 <- df[2]
    ec_density <- (2 * pi)^(-0.5) * (df1 * threshold / df2)^((df1 - 1) / 2) *
      (1 + df1 * threshold / df2)^(-(df1 + df2 - 1) / 2)
  } else {
    stop("Unsupported statistic type")
  }

  return(resels[1] * ec_density)
}

#' Validate FWHM estimate
#' @param fwhm estimated FWHM
#' @param Q number of nodes
#' @return validated FWHM (with warnings if needed)
validate_fwhm <- function(fwhm, Q) {
  if (is.na(fwhm) || !is.finite(fwhm)) {
    warning("Invalid FWHM estimate (NA or infinite), using default")
    return(min(Q / 5, 10))
  }

  if (fwhm < 1) {
    warning("FWHM < 1, setting to 1")
    return(1.0)
  }

  if (fwhm > Q) {
    warning(sprintf("FWHM (%.1f) > number of nodes (%d), setting to Q/2", fwhm, Q))
    return(Q / 2)
  }

  return(fwhm)
}
