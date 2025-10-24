#' spm1d R vs Python Validation Tests
#'
#' Generate test datasets and run analyses in both R and Python
#' to validate equivalence of implementations

# ==============================================================================
# TEST 1: TWO-SAMPLE T-TEST
# ==============================================================================

cat("\n=== TEST 1: Two-Sample T-Test ===\n")

# Set seed for reproducibility
set.seed(42)

# Generate two groups
n_per_group <- 8
n_timepoints <- 101

YA_2samp <- matrix(rnorm(n_per_group * n_timepoints), n_per_group, n_timepoints)
YB_2samp <- matrix(rnorm(n_per_group * n_timepoints, mean = 0.5),
                   n_per_group, n_timepoints)

# Save for Python
write.csv(YA_2samp, "test1_groupA.csv", row.names = FALSE)
write.csv(YB_2samp, "test1_groupB.csv", row.names = FALSE)

# Run R analysis
spm_r <- ttest2(YA_2samp, YB_2samp)
spmi_r <- inference(spm_r, alpha = 0.05, two_tailed = TRUE)

# Save R results
results_r <- data.frame(
  t_values = spm_r$z,
  fwhm = spm_r$fwhm,
  resels = spm_r$resels,
  threshold = spmi_r$zstar,
  p_value = spmi_r$p_set,
  df1 = spm_r$df[1],
  df2 = spm_r$df[2]
)
write.csv(results_r, "test1_results_r.csv", row.names = FALSE)

cat("Saved: test1_groupA.csv, test1_groupB.csv, test1_results_r.csv\n")
cat("R Threshold:", spmi_r$zstar, "\n")
cat("R FWHM:", spm_r$fwhm, "\n")
cat("R Resels:", spm_r$resels, "\n")

# Python comparison code:
cat("\n### Python code to run:\n")
cat("import spm1d\nimport pandas as pd\nimport numpy as np\n\n")
cat("YA = pd.read_csv('test1_groupA.csv').values\n")
cat("YB = pd.read_csv('test1_groupB.csv').values\n")
cat("spm = spm1d.stats.ttest2(YA, YB)\n")
cat("spmi = spm.inference(alpha=0.05, two_tailed=True)\n")
cat("results_py = pd.DataFrame({\n")
cat("    't_values': spm.z,\n")
cat("    'fwhm': [spm.fwhm] * len(spm.z),\n")
cat("    'resels': [float(spm.resels[1])] * len(spm.z),\n")  # Ensure float conversion
cat("    'threshold': [spmi.zstar] * len(spm.z),\n")
cat("    'p_value': [spmi.p_set if hasattr(spmi, 'p_set') else 0] * len(spm.z),\n")
cat("    'df1': [spm.df[0]] * len(spm.z),\n")
cat("    'df2': [float(spm.df[1])] * len(spm.z)\n")  # Ensure float conversion
cat("})\n")
cat("results_py.to_csv('test1_results_py.csv', index=False)\n")
cat("print(f'Threshold: {spmi.zstar}')\n")
cat("print(f'FWHM: {spm.fwhm}')\n")
cat("print(f'Resels: {spm.resels[0]}')\n\n")

# ==============================================================================
# TEST 2: PAIRED T-TEST
# ==============================================================================

cat("\n=== TEST 2: Paired T-Test ===\n")

set.seed(123)

# Generate paired data (pre/post design)
n_subjects <- 10
Y_pre <- matrix(rnorm(n_subjects * n_timepoints), n_subjects, n_timepoints)
Y_post <- Y_pre + matrix(rnorm(n_subjects * n_timepoints, mean = 0.8, sd = 0.5),
                         n_subjects, n_timepoints)

# Save for Python
write.csv(Y_pre, "test2_pre.csv", row.names = FALSE)
write.csv(Y_post, "test2_post.csv", row.names = FALSE)

# Run R analysis
spm_r <- ttest_paired(Y_pre, Y_post)
spmi_r <- inference(spm_r, alpha = 0.05, two_tailed = TRUE)

# Save R results
results_r <- data.frame(
  t_values = spm_r$z,
  fwhm = spm_r$fwhm,
  resels = spm_r$resels,
  threshold = spmi_r$zstar,
  p_value = spmi_r$p_set,
  df = spm_r$df[2]
)
write.csv(results_r, "test2_results_r.csv", row.names = FALSE)

cat("Saved: test2_pre.csv, test2_post.csv, test2_results_r.csv\n")
cat("R Threshold:", spmi_r$zstar, "\n")

cat("\n### Python code to run:\n")
cat("Y_pre = pd.read_csv('test2_pre.csv').values\n")
cat("Y_post = pd.read_csv('test2_post.csv').values\n")
cat("spm = spm1d.stats.ttest_paired(Y_pre, Y_post)\n")
cat("spmi = spm.inference(alpha=0.05, two_tailed=True)\n")
cat("results_py = pd.DataFrame({\n")
cat("    't_values': spm.z,\n")
cat("    'fwhm': [spm.fwhm] * len(spm.z),\n")
cat("    'resels': [spm.resels[0]] * len(spm.z),\n")
cat("    'threshold': [spmi.zstar] * len(spm.z),\n")
cat("    'p_value': [spmi.p] * len(spm.z),\n")
cat("    'df': [spm.df[1]] * len(spm.z)\n")
cat("})\n")
cat("results_py.to_csv('test2_results_py.csv', index=False)\n\n")

# ==============================================================================
# TEST 3: LINEAR REGRESSION
# ==============================================================================

cat("\n=== TEST 3: Linear Regression ===\n")

set.seed(456)

# Generate data with linear relationship to predictor
n_subjects <- 15
predictor <- seq(20, 80, length.out = n_subjects)  # e.g., age

# Create age effect
Y_reg <- matrix(rnorm(n_subjects * n_timepoints), n_subjects, n_timepoints)
age_effect <- outer(scale(predictor)[,1],
                    ifelse(1:n_timepoints > 40 & 1:n_timepoints < 70, 1.5, 0))
Y_reg <- Y_reg + age_effect

# Save for Python
write.csv(Y_reg, "test3_Y.csv", row.names = FALSE)
write.csv(data.frame(predictor = predictor), "test3_predictor.csv", row.names = FALSE)

# Run R analysis
spm_r <- regress(Y_reg, predictor)
spmi_r <- inference(spm_r, alpha = 0.05, two_tailed = FALSE)  # one-tailed

# Save R results
results_r <- data.frame(
  t_values = spm_r$z,
  fwhm = spm_r$fwhm,
  resels = spm_r$resels,
  threshold = spmi_r$zstar,
  p_value = spmi_r$p_set,
  df = spm_r$df[2]
)
write.csv(results_r, "test3_results_r.csv", row.names = FALSE)

cat("Saved: test3_Y.csv, test3_predictor.csv, test3_results_r.csv\n")
cat("R Threshold:", spmi_r$zstar, "\n")

cat("\n### Python code to run:\n")
cat("Y = pd.read_csv('test3_Y.csv').values\n")
cat("x = pd.read_csv('test3_predictor.csv')['predictor'].values\n")
cat("spm = spm1d.stats.regress(Y, x)\n")
cat("spmi = spm.inference(alpha=0.05, two_tailed=False)\n")
cat("results_py = pd.DataFrame({\n")
cat("    't_values': spm.z,\n")
cat("    'fwhm': [spm.fwhm] * len(spm.z),\n")
cat("    'resels': [spm.resels[0]] * len(spm.z),\n")
cat("    'threshold': [spmi.zstar] * len(spm.z),\n")
cat("    'p_value': [spmi.p] * len(spm.z),\n")
cat("    'df': [spm.df[1]] * len(spm.z)\n")
cat("})\n")
cat("results_py.to_csv('test3_results_py.csv', index=False)\n\n")

# ==============================================================================
# TEST 4: ONE-WAY ANOVA
# ==============================================================================

cat("\n=== TEST 4: One-Way ANOVA ===\n")

set.seed(789)

# Generate three groups
n_per_group <- 6

Y1_anova <- matrix(rnorm(n_per_group * n_timepoints), n_per_group, n_timepoints)
Y2_anova <- matrix(rnorm(n_per_group * n_timepoints, mean = 0.5),
                   n_per_group, n_timepoints)
Y3_anova <- matrix(rnorm(n_per_group * n_timepoints, mean = 1.0),
                   n_per_group, n_timepoints)

# Save for Python
write.csv(Y1_anova, "test4_group1.csv", row.names = FALSE)
write.csv(Y2_anova, "test4_group2.csv", row.names = FALSE)
write.csv(Y3_anova, "test4_group3.csv", row.names = FALSE)

# Run R analysis
spm_r <- anova1(list(Y1_anova, Y2_anova, Y3_anova))
spmi_r <- inference(spm_r, alpha = 0.05)

# Save R results
results_r <- data.frame(
  F_values = spm_r$z,
  fwhm = spm_r$fwhm,
  resels = spm_r$resels,
  threshold = spmi_r$zstar,
  p_value = spmi_r$p_set,
  df1 = spm_r$df[1],
  df2 = spm_r$df[2]
)
write.csv(results_r, "test4_results_r.csv", row.names = FALSE)

cat("Saved: test4_group1.csv, test4_group2.csv, test4_group3.csv, test4_results_r.csv\n")
cat("R Threshold:", spmi_r$zstar, "\n")
cat("R FWHM:", spm_r$fwhm, "\n")

cat("\n### Python code to run:\n")
cat("Y1 = pd.read_csv('test4_group1.csv').values\n")
cat("Y2 = pd.read_csv('test4_group2.csv').values\n")
cat("Y3 = pd.read_csv('test4_group3.csv').values\n")
cat("spm = spm1d.stats.anova1([Y1, Y2, Y3], equal_var=True)\n")
cat("spmi = spm.inference(alpha=0.05)\n")
cat("results_py = pd.DataFrame({\n")
cat("    'F_values': spm.z,\n")
cat("    'fwhm': [spm.fwhm] * len(spm.z),\n")
cat("    'resels': [spm.resels[0]] * len(spm.z),\n")
cat("    'threshold': [spmi.zstar] * len(spm.z),\n")
cat("    'p_value': [spmi.p] * len(spm.z),\n")
cat("    'df1': [spm.df[0]] * len(spm.z),\n")
cat("    'df2': [spm.df[1]] * len(spm.z)\n")
cat("})\n")
cat("results_py.to_csv('test4_results_py.csv', index=False)\n\n")

# ==============================================================================
# TEST 5: REPEATED MEASURES ANOVA
# ==============================================================================

cat("\n=== TEST 5: Repeated Measures ANOVA ===\n")

set.seed(101112)

# Generate repeated measures data (3 conditions, same subjects)
n_subjects <- 8

Y_cond1 <- matrix(rnorm(n_subjects * n_timepoints), n_subjects, n_timepoints)
Y_cond2 <- matrix(rnorm(n_subjects * n_timepoints, mean = 0.3),
                  n_subjects, n_timepoints)
Y_cond3 <- matrix(rnorm(n_subjects * n_timepoints, mean = 0.8),
                  n_subjects, n_timepoints)

# Save for Python
write.csv(Y_cond1, "test5_cond1.csv", row.names = FALSE)
write.csv(Y_cond2, "test5_cond2.csv", row.names = FALSE)
write.csv(Y_cond3, "test5_cond3.csv", row.names = FALSE)

# Run R analysis
spm_r <- anova1rm(list(Y_cond1, Y_cond2, Y_cond3))
spmi_r <- inference(spm_r, alpha = 0.05)

# Save R results
results_r <- data.frame(
  F_values = spm_r$z,
  fwhm = spm_r$fwhm,
  resels = spm_r$resels,
  threshold = spmi_r$zstar,
  p_value = spmi_r$p_set,
  df1 = spm_r$df[1],
  df2 = spm_r$df[2]
)
write.csv(results_r, "test5_results_r.csv", row.names = FALSE)

cat("Saved: test5_cond1.csv, test5_cond2.csv, test5_cond3.csv, test5_results_r.csv\n")
cat("R Threshold:", spmi_r$zstar, "\n")

cat("\n### Python code to run:\n")
cat("Y1 = pd.read_csv('test5_cond1.csv').values\n")
cat("Y2 = pd.read_csv('test5_cond2.csv').values\n")
cat("Y3 = pd.read_csv('test5_cond3.csv').values\n")
cat("spm = spm1d.stats.anova1rm([Y1, Y2, Y3], equal_var=True)\n")
cat("spmi = spm.inference(alpha=0.05)\n")
cat("results_py = pd.DataFrame({\n")
cat("    'F_values': spm.z,\n")
cat("    'fwhm': [spm.fwhm] * len(spm.z),\n")
cat("    'resels': [spm.resels[0]] * len(spm.z),\n")
cat("    'threshold': [spmi.zstar] * len(spm.z),\n")
cat("    'p_value': [spmi.p] * len(spm.z),\n")
cat("    'df1': [spm.df[0]] * len(spm.z),\n")
cat("    'df2': [spm.df[1]] * len(spm.z)\n")
cat("})\n")
cat("results_py.to_csv('test5_results_py.csv', index=False)\n\n")

# ==============================================================================
# TEST 6: GLM WITH CUSTOM CONTRAST
# ==============================================================================

cat("\n=== TEST 6: GLM with Custom Contrast ===\n")

set.seed(131415)

# Generate 2x2 factorial design data
n_per_cell <- 5
n_subjects <- n_per_cell * 4

# Factors: A (2 levels), B (2 levels)
factor_A <- rep(c(0, 0, 1, 1), each = n_per_cell)
factor_B <- rep(c(0, 1, 0, 1), each = n_per_cell)

# Design matrix
X <- cbind(1, factor_A, factor_B, factor_A * factor_B)

# Generate data with interaction effect
Y_glm <- matrix(rnorm(n_subjects * n_timepoints), n_subjects, n_timepoints)
interaction_effect <- ifelse(1:n_timepoints > 50 & 1:n_timepoints < 80, 1.2, 0)
for (i in (3*n_per_cell+1):(4*n_per_cell)) {  # A=1, B=1 cell
  Y_glm[i, ] <- Y_glm[i, ] + interaction_effect
}

# Test interaction contrast
c_interaction <- c(0, 0, 0, 1)

# Save for Python
write.csv(Y_glm, "test6_Y.csv", row.names = FALSE)
write.csv(as.data.frame(X), "test6_X.csv", row.names = FALSE)
write.csv(data.frame(contrast = c_interaction), "test6_contrast.csv", row.names = FALSE)

# Run R analysis
spm_r <- glm(Y_glm, X, c_interaction)
spmi_r <- inference(spm_r, alpha = 0.05, two_tailed = TRUE)

# Save R results
results_r <- data.frame(
  t_values = spm_r$z,
  fwhm = spm_r$fwhm,
  resels = spm_r$resels,
  threshold = spmi_r$zstar,
  p_value = spmi_r$p_set,
  df = spm_r$df[2]
)
write.csv(results_r, "test6_results_r.csv", row.names = FALSE)

cat("Saved: test6_Y.csv, test6_X.csv, test6_contrast.csv, test6_results_r.csv\n")
cat("R Threshold:", spmi_r$zstar, "\n")

cat("\n### Python code to run:\n")
cat("Y = pd.read_csv('test6_Y.csv').values\n")
cat("X = pd.read_csv('test6_X.csv').values\n")
cat("c = pd.read_csv('test6_contrast.csv')['contrast'].values\n")
cat("spm = spm1d.stats.glm(Y, X, c)\n")
cat("spmi = spm.inference(alpha=0.05, two_tailed=True)\n")
cat("results_py = pd.DataFrame({\n")
cat("    't_values': spm.z,\n")
cat("    'fwhm': [spm.fwhm] * len(spm.z),\n")
cat("    'resels': [spm.resels[0]] * len(spm.z),\n")
cat("    'threshold': [spmi.zstar] * len(spm.z),\n")
cat("    'p_value': [spmi.p] * len(spm.z),\n")
cat("    'df': [spm.df[1]] * len(spm.z)\n")
cat("})\n")
cat("results_py.to_csv('test6_results_py.csv', index=False)\n\n")

# ==============================================================================
# AUTOMATED COMPARISON FUNCTION
# ==============================================================================

compare_results <- function(test_name) {
  # Read R and Python results
  r_file <- paste0(test_name, '_results_r.csv')
  py_file <- paste0(test_name, '_results_py.csv')

  results_r <- read.csv(r_file)
  results_py <- read.csv(py_file)

  cat('\\n====================\\n')
  cat('Test:', test_name, '\\n')
  cat('====================\\n')

  # Compare statistic values (t or F)
  stat_col <- if('t_values' %in% names(results_r)) 't_values' else 'F_values'

  # Check for sign flip (Python and R might use opposite conventions)
  cor_stat_pos <- cor(results_r[[stat_col]], results_py[[stat_col]])
  cor_stat_neg <- cor(results_r[[stat_col]], -results_py[[stat_col]])

  if (abs(cor_stat_neg) > abs(cor_stat_pos)) {
    cat('NOTE: Sign convention differs between R and Python\\n')
    cor_stat <- cor_stat_neg
    max_diff_stat <- max(abs(results_r[[stat_col]] - (-results_py[[stat_col]])))
  } else {
    cor_stat <- cor_stat_pos
    max_diff_stat <- max(abs(results_r[[stat_col]] - results_py[[stat_col]]))
  }

  cat('Statistic correlation:', cor_stat, '\\n')
  cat('Max statistic difference:', max_diff_stat, '\\n')

  # Compare FWHM
  fwhm_r <- results_r$fwhm[1]
  fwhm_py <- results_py$fwhm[1]
  fwhm_diff_pct <- abs(fwhm_r - fwhm_py) / fwhm_py * 100

  cat('FWHM - R:', fwhm_r, ', Python:', fwhm_py)
  cat(' (', round(fwhm_diff_pct, 2), '% diff)\\n', sep='')

  # Compare Resels
  resels_r <- results_r$resels[1]
  resels_py <- results_py$resels[1]
  resels_diff_pct <- abs(resels_r - resels_py) / resels_py * 100

  cat('Resels - R:', resels_r, ', Python:', resels_py)
  cat(' (', round(resels_diff_pct, 2), '% diff)\\n', sep='')

  # Compare Threshold
  thresh_r <- results_r$threshold[1]
  thresh_py <- results_py$threshold[1]
  thresh_diff_pct <- abs(thresh_r - thresh_py) / thresh_py * 100

  cat('Threshold - R:', thresh_r, ', Python:', thresh_py)
  cat(' (', round(thresh_diff_pct, 2), '% diff)\\n', sep='')

  # Compare degrees of freedom
  if('df' %in% names(results_r)) {
    cat('df - R:', results_r$df[1], ', Python:', results_py$df[1], '\\n')
  } else {
    cat('df1 - R:', results_r$df1[1], ', Python:', results_py$df1[1], '\\n')
    cat('df2 - R:', results_r$df2[1], ', Python:', results_py$df2[1], '\\n')
  }

  # Overall assessment
  cat('\\nAssessment: ')
  if(cor_stat > 0.9999 && thresh_diff_pct < 2) {
    cat('PASS - Excellent agreement\\n')
  } else if(cor_stat > 0.999 && thresh_diff_pct < 5) {
    cat('PASS - Good agreement\\n')
  } else {
    cat('FAIL - Significant differences\\n')
  }

  return(list(
    stat_cor = cor_stat,
    stat_max_diff = max_diff_stat,
    fwhm_diff_pct = fwhm_diff_pct,
    resels_diff_pct = resels_diff_pct,
    thresh_diff_pct = thresh_diff_pct
  ))
}

# Run comparisons for all tests
compare_results('test1')
compare_results('test2')
compare_results('test3')
compare_results('test4')
compare_results('test5')
compare_results('test6')

cat("\n=== All test data generated! ===\n")
cat("Run the Python code snippets above, then use compare_results() to validate.\n")
