#' spm1d R Package - Statistical Functions Module
#'
#' Implements main statistical tests: GLM, t-tests, ANOVA, regression

# Ensure required packages
#if (!require("MASS")) install.packages("MASS")
#if (!require("Matrix")) install.packages("Matrix")

#library(MASS)  # for ginv (generalized inverse)
#library(Matrix)  # for rankMatrix

# Machine epsilon
eps <- .Machine$double.eps

# ==============================================================================
# GENERAL LINEAR MODEL (GLM)
# ==============================================================================

#' General Linear Model for t contrasts
#'
#' @param Y dependent variable (J x Q matrix: J observations, Q nodes)
#' @param X design matrix (J x B: J observations, B parameters)
#' @param c contrast vector (length B)
#' @param Q non-sphericity specifiers (not currently supported)
#' @param roi region of interest mask (logical vector length Q)
#' @importFrom MASS ginv
#' @importFrom Matrix rankMatrix
#' @return SPM_T or SPM0D_T object
#' @export
#'
#' @examples
#' # Example: two-sample t-test via GLM
#' Y <- rbind(matrix(rnorm(50), 5, 10), matrix(rnorm(50, 1), 5, 10))
#' X <- cbind(1, c(rep(0, 5), rep(1, 5)))
#' c <- c(0, 1)  # test second parameter
#' spm <- glm(Y, X, c)
#' spmi <- inference(spm, alpha = 0.05, two_tailed = TRUE, withBonf = TRUE)
glm <- function(Y, X, c, Q = NULL, roi = NULL) {

  ### Assemble data
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  c <- matrix(c, ncol = 1)  # ensure column vector

  ### Solve the GLM
  b <- ginv(X) %*% Y  # parameters
  eij <- Y - X %*% b  # residuals
  df <- nrow(Y) - rankMatrix(X)[1]  # degrees of freedom

  # Calculate variance (fast method using column sums)
  diagR <- colSums(eij^2)  # residual sum of squares
  sigma2 <- diagR / df  # variance

  ### Compute t statistic
  cXXc <- t(c) %*% solve(t(X) %*% X) %*% c
  cXXc <- as.numeric(cXXc)  # convert to scalar

  t_stat <- as.vector(t(c) %*% b) / (sqrt(sigma2 * cXXc + eps))

  ### Estimate df due to non-sphericity (if specified)
  if (!is.null(Q)) {
    warning("Non-sphericity estimation not yet implemented")
    # df <- estimate_df_T(Y, X, eij, Q)
  }

  eij <- as.matrix(eij)

  ### Return appropriate SPM object
  if (ncol(Y) > 1) {
    # 1D case: continuum data

    ### Estimate field smoothness
    fwhm <- estimate_fwhm(eij)

    ### Compute resel counts
    if (is.null(roi)) {
      resels <- resel_counts_from_residuals(eij, fwhm)
    } else {
      # Apply ROI mask
      B <- apply(is.na(eij), 2, any)  # nodes with NaN
      B <- !B & roi  # nodes in ROI and not NaN
      mask <- !B  # TRUE for masked-out regions
      resels <- resel_counts(mask, fwhm, element_based = FALSE)
      t_stat <- ifelse(roi, t_stat, NA)  # mask t values outside ROI
    }

    ### Assemble SPM{t} object
    s <- as.vector(sigma2)
    result <- SPM_T(t_stat, c(1, df), fwhm, resels, X, b, eij,
                    sigma2 = s, roi = roi)

  } else {
    # 0D case: single observation
    b_vec <- as.vector(b)
    r <- as.vector(eij)
    s2 <- as.numeric(sigma2[1])

    result <- SPM0D_T(t_stat, c(1, df), beta = b_vec,
                      residuals = r, sigma2 = s2)
  }

  return(result)
}

# ==============================================================================
# T-TESTS
# ==============================================================================

#' One-sample t-test
#'
#' @param Y data matrix (J x Q: J observations, Q nodes)
#' @param mu null hypothesis mean (default 0)
#' @param roi region of interest mask
#' @return SPM_T or SPM0D_T object
#' @export
#'
#' @examples
#' Y <- matrix(rnorm(50, mean = 0.5), 5, 10)
#' spm <- ttest(Y)
#' spmi <- inference(spm, alpha = 0.05, two_tailed = TRUE)
ttest <- function(Y, mu = 0, roi = NULL) {

  Y <- as.matrix(Y)
  J <- nrow(Y)
  Q <- ncol(Y)

  # Center data by subtracting mu
  Y <- Y - mu

  # Design matrix: intercept only
  X <- matrix(1, nrow = J, ncol = 1)

  # Contrast: test if intercept differs from zero
  c <- 1

  # Use GLM
  result <- glm(Y, X, c, roi = roi)

  return(result)
}

#' Two-sample t-test
#'
#' @param YA data for group A (JA x Q matrix)
#' @param YB data for group B (JB x Q matrix)
#' @param equal_var assume equal variance (default TRUE)
#' @param roi region of interest mask
#' @return SPM_T or SPM0D_T object
#' @export
#'
#' @examples
#' YA <- matrix(rnorm(50), 5, 10)
#' YB <- matrix(rnorm(50, mean = 1), 5, 10)
#' spm <- ttest2(YA, YB)
#' spmi <- inference(spm, alpha = 0.05)
ttest2 <- function(YA, YB, equal_var = TRUE, roi = NULL) {

  YA <- as.matrix(YA)
  YB <- as.matrix(YB)

  JA <- nrow(YA)
  JB <- nrow(YB)

  if (ncol(YA) != ncol(YB)) {
    stop("YA and YB must have the same number of columns (nodes)")
  }

  # Combine data
  Y <- rbind(YA, YB)

  # Design matrix: intercept + group indicator
  X <- cbind(1, c(rep(0, JA), rep(1, JB)))

  # Contrast: test difference between groups
  c <- c(0, 1)

  # Use GLM
  result <- glm(Y, X, c, roi = roi)

  return(result)
}

#' Paired t-test
#'
#' @param YA data for condition A (J x Q matrix)
#' @param YB data for condition B (J x Q matrix, same subjects)
#' @param roi region of interest mask
#' @return SPM_T or SPM0D_T object
#' @export
#'
#' @examples
#' YA <- matrix(rnorm(50), 5, 10)
#' YB <- YA + matrix(rnorm(50, mean = 0.5), 5, 10)
#' spm <- ttest_paired(YA, YB)
#' spmi <- inference(spm, alpha = 0.05, two_tailed = TRUE)
ttest_paired <- function(YA, YB, roi = NULL) {

  YA <- as.matrix(YA)
  YB <- as.matrix(YB)

  if (nrow(YA) != nrow(YB) || ncol(YA) != ncol(YB)) {
    stop("YA and YB must have the same dimensions")
  }

  # Compute differences
  Y_diff <- YA - YB

  # One-sample t-test on differences
  result <- ttest(Y_diff, mu = 0, roi = roi)

  return(result)
}

# ==============================================================================
# REGRESSION
# ==============================================================================

#' Linear regression (correlation analysis)
#'
#' @param Y dependent variable (J x Q matrix)
#' @param x independent variable (length J vector)
#' @param roi region of interest mask
#' @return SPM_T or SPM0D_T object
#' @export
#'
#' @examples
#' x <- 1:10
#' Y <- outer(x, rep(1, 5)) + matrix(rnorm(50), 10, 5)
#' spm <- regress(Y, x)
#' spmi <- inference(spm, alpha = 0.05)
regress <- function(Y, x, roi = NULL) {

  Y <- as.matrix(Y)
  x <- as.vector(x)

  if (nrow(Y) != length(x)) {
    stop("Number of rows in Y must equal length of x")
  }

  # Design matrix: intercept + regressor
  X <- cbind(1, x)

  # Contrast: test slope parameter
  c <- c(0, 1)

  # Use GLM
  result <- glm(Y, X, c, roi = roi)

  return(result)
}

# ==============================================================================
# ONE-WAY ANOVA
# ==============================================================================

#' One-way ANOVA
#'
#' @param Y list of data matrices, one per group
#' @param roi region of interest mask
#' @importFrom MASS ginv
#' @return SPM_F or SPM0D_F object
#' @export
#'
#' @examples
#' Y1 <- matrix(rnorm(30), 3, 10)
#' Y2 <- matrix(rnorm(30, mean = 0.5), 3, 10)
#' Y3 <- matrix(rnorm(30, mean = 1), 3, 10)
#' spm <- anova1(list(Y1, Y2, Y3))
#' spmi <- inference(spm, alpha = 0.05, withBonf = TRUE)
anova1 <- function(Y, roi = NULL) {

  if (!is.list(Y)) {
    stop("Y must be a list of matrices, one per group")
  }

  # Get group sizes and total N
  K <- length(Y)  # number of groups
  J <- sapply(Y, nrow)  # observations per group
  J_total <- sum(J)

  # Check that all groups have same number of columns
  Q <- ncol(Y[[1]])
  if (!all(sapply(Y, ncol) == Q)) {
    stop("All groups must have the same number of columns (nodes)")
  }

  # Convert to matrices
  Y <- lapply(Y, as.matrix)

  # Combine all data
  Y_all <- do.call(rbind, Y)

  # Create design matrix (one column per group)
  X <- matrix(0, nrow = J_total, ncol = K)
  start_idx <- 1
  for (k in 1:K) {
    end_idx <- start_idx + J[k] - 1
    X[start_idx:end_idx, k] <- 1
    start_idx <- end_idx + 1
  }

  # Reduced model: grand mean only
  X0 <- matrix(1, nrow = J_total, ncol = 1)

  # Fit full model
  b_full <- ginv(X) %*% Y_all
  eij_full <- Y_all - X %*% b_full

  # Fit reduced model
  b_reduced <- ginv(X0) %*% Y_all
  eij_reduced <- Y_all - X0 %*% b_reduced

  # Calculate F statistic
  # SS_effect = SS_reduced - SS_full
  SS_full <- colSums(eij_full^2)
  SS_reduced <- colSums(eij_reduced^2)
  SS_effect <- SS_reduced - SS_full

  # Degrees of freedom
  df_effect <- K - 1
  df_error <- J_total - K

  # Mean squares
  MS_effect <- SS_effect / df_effect
  MS_error <- SS_full / df_error

  # F statistic
  F_stat <- MS_effect / (MS_error + eps)

  # Variance
  sigma2 <- MS_error

  ### Return appropriate SPM object
  if (Q > 1) {
    # 1D case

    ### Estimate field smoothness
    fwhm <- estimate_fwhm(eij_full)

    ### Compute resel counts
    if (is.null(roi)) {
      resels <- resel_counts_from_residuals(eij_full, fwhm)
    } else {
      B <- apply(is.na(eij_full), 2, any)
      B <- !B & roi
      mask <- !B
      resels <- resel_counts(mask, fwhm, element_based = FALSE)
      F_stat <- ifelse(roi, F_stat, NA)
    }

    result <- SPM_F(F_stat, c(df_effect, df_error), fwhm, resels,
                    X, b_full, eij_full, X0, sigma2 = as.vector(sigma2),
                    roi = roi)

  } else {
    # 0D case
    result <- SPM0D_F(F_stat, c(df_effect, df_error),
                      X, as.vector(b_full), as.vector(eij_full),
                      sigma2 = as.numeric(sigma2[1]))
  }

  return(result)
}

#' Two-way ANOVA (factorial design)
#'
#' @param Y data matrix (J x Q)
#' @param A factor A (length J vector)
#' @param B factor B (length J vector)
#' @param roi region of interest mask
#' @return list of SPM_F objects for main effects and interaction
#' @importFrom stats model.matrix
#' @export
anova2 <- function(Y, A, B, roi = NULL) {

  Y <- as.matrix(Y)
  A <- as.factor(A)
  B <- as.factor(B)

  if (nrow(Y) != length(A) || nrow(Y) != length(B)) {
    stop("Factor lengths must match number of rows in Y")
  }

  J <- nrow(Y)

  # Create design matrix using model.matrix
  # Include main effects and interaction
  formula_full <- ~ A * B
  X <- model.matrix(formula_full)

  # Get indices for different effects
  # This is simplified - in practice need better parsing
  n_a <- nlevels(A)
  n_b <- nlevels(B)

  # Main effect A
  idx_a <- 2:n_a
  # Main effect B
  idx_b <- (n_a + 1):(n_a + n_b - 1)
  # Interaction
  idx_ab <- (n_a + n_b):ncol(X)

  # Fit full model
  b_full <- MASS::ginv(X) %*% Y
  eij_full <- Y - X %*% b_full
  SS_full <- colSums(eij_full^2)
  df_error <- J - Matrix::rankMatrix(X)[1]
  MS_error <- SS_full / df_error

  # Function to test subset of parameters
  test_effect <- function(idx_effect) {
    # Reduced model without these parameters
    X_reduced <- X[, -idx_effect, drop = FALSE]
    b_reduced <- MASS::ginv(X_reduced) %*% Y
    eij_reduced <- Y - X_reduced %*% b_reduced
    SS_reduced <- colSums(eij_reduced^2)

    SS_effect <- SS_reduced - SS_full
    df_effect <- length(idx_effect)
    MS_effect <- SS_effect / df_effect

    F_stat <- MS_effect / (MS_error + eps)

    return(list(F = F_stat, df = c(df_effect, df_error)))
  }

  # Test each effect
  effect_A <- test_effect(idx_a)
  effect_B <- test_effect(idx_b)
  effect_AB <- test_effect(idx_ab)

  # Create SPM objects for each effect
  Q <- ncol(Y)

  create_spm_f <- function(F_stat, df_vec) {
    if (Q > 1) {
      fwhm <- estimate_fwhm(eij_full)
      if (is.null(roi)) {
        resels <- resel_counts_from_residuals(eij_full, fwhm)
      } else {
        B <- apply(is.na(eij_full), 2, any)
        B <- !B & roi
        mask <- !B
        resels <- resel_counts(mask, fwhm, element_based = FALSE)
        F_stat <- ifelse(roi, F_stat, NA)
      }
      return(SPM_F(F_stat, df_vec, fwhm, resels, X, b_full, eij_full,
                   sigma2 = as.vector(MS_error), roi = roi))
    } else {
      return(SPM0D_F(F_stat, df_vec, X, as.vector(b_full),
                     as.vector(eij_full), sigma2 = as.numeric(MS_error[1])))
    }
  }

  result <- list(
    A = create_spm_f(effect_A$F, effect_A$df),
    B = create_spm_f(effect_B$F, effect_B$df),
    AB = create_spm_f(effect_AB$F, effect_AB$df)
  )

  return(result)
}

# ==============================================================================
# REPEATED MEASURES ANOVA (ONE-WAY)
# ==============================================================================

#' One-way repeated measures ANOVA
#'
#' @param Y list of data matrices, one per condition (same subjects)
#' @param roi region of interest mask
#' @return SPM_F or SPM0D_F object
#' @export
anova1rm <- function(Y, roi = NULL) {

  if (!is.list(Y)) {
    stop("Y must be a list of matrices, one per condition")
  }

  K <- length(Y)  # number of conditions
  J <- nrow(Y[[1]])  # number of subjects
  Q <- ncol(Y[[1]])  # number of nodes

  # Check dimensions
  if (!all(sapply(Y, nrow) == J)) {
    stop("All conditions must have the same number of subjects")
  }
  if (!all(sapply(Y, ncol) == Q)) {
    stop("All conditions must have the same number of nodes")
  }

  # Convert to matrices
  Y <- lapply(Y, as.matrix)

  # Stack data: each subject contributes K observations
  Y_stacked <- do.call(rbind, Y)

  # Design matrix with subject effects
  # Columns: conditions (K) + subjects (J)
  X <- matrix(0, nrow = J * K, ncol = K + J)

  # Condition effects
  for (k in 1:K) {
    idx <- ((k - 1) * J + 1):(k * J)
    X[idx, k] <- 1
  }

  # Subject effects
  for (j in 1:J) {
    idx <- seq(j, J * K, by = J)
    X[idx, K + j] <- 1
  }

  # Reduced model: subjects only (no condition effect)
  X0 <- X[, (K + 1):(K + J), drop = FALSE]

  # Fit models
  b_full <- MASS::ginv(X) %*% Y_stacked
  eij_full <- Y_stacked - X %*% b_full

  b_reduced <- MASS::ginv(X0) %*% Y_stacked
  eij_reduced <- Y_stacked - X0 %*% b_reduced

  # Calculate F statistic
  SS_full <- colSums(eij_full^2)
  SS_reduced <- colSums(eij_reduced^2)
  SS_effect <- SS_reduced - SS_full

  df_effect <- K - 1
  df_error <- (J - 1) * (K - 1)

  MS_effect <- SS_effect / df_effect
  MS_error <- SS_full / df_error

  F_stat <- MS_effect / (MS_error + eps)

  ### Return appropriate SPM object
  if (Q > 1) {
    fwhm <- estimate_fwhm(eij_full)

    if (is.null(roi)) {
      resels <- resel_counts_from_residuals(eij_full, fwhm)
    } else {
      B <- apply(is.na(eij_full), 2, any)
      B <- !B & roi
      mask <- !B
      resels <- resel_counts(mask, fwhm, element_based = FALSE)
      F_stat <- ifelse(roi, F_stat, NA)
    }

    result <- SPM_F(F_stat, c(df_effect, df_error), fwhm, resels,
                    X, b_full, eij_full, X0,
                    sigma2 = as.vector(MS_error), roi = roi)

  } else {
    result <- SPM0D_F(F_stat, c(df_effect, df_error),
                      X, as.vector(b_full), as.vector(eij_full),
                      sigma2 = as.numeric(MS_error[1]))
  }

  return(result)
}

# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================

#' Center data (remove mean)
#' @param Y data matrix
#' @return centered data
center_data <- function(Y) {
  Y <- as.matrix(Y)
  Y_centered <- sweep(Y, 2, colMeans(Y), "-")
  return(Y_centered)
}

#' Standardize data (mean 0, sd 1)
#' @param Y data matrix
#' @return standardized data
#' @importFrom stats sd
standardize_data <- function(Y) {
  Y <- as.matrix(Y)
  Y_std <- sweep(Y, 2, colMeans(Y), "-")
  Y_std <- sweep(Y_std, 2, apply(Y, 2, sd), "/")
  return(Y_std)
}
