###############################################################################
#  STAT 6545 - Section 4.7: Canonical Correlation Analysis (CCA)
#  Credit-card customer data: Financial capacity  vs  Spending / usage behavior
#
#  Follows Johnson & Wichern (6th ed.), Chapter 10:
#    10.2  Result 10.1/10.2  -> canonical correlations & coefficient vectors
#    10.3  Canonical structure (variate-variable correlations, 10-14/10-29)
#    10.4  Sample canonical variates, scatterplots, split-sample validation
#    10.5  Matrices of errors of approximation (10-34)
#          Proportion of explained sample variance (10-37) + redundancy
#    10.6  Bartlett large-sample LR tests (10-39) and sequential tests (10-41)
#

###############################################################################

args     <- commandArgs(trailingOnly = TRUE)
data_csv <- if (length(args) >= 1) args[1] else "synthetic_credit_card_customer_behavior_dataset.csv"
out_dir  <- if (length(args) >= 2) args[2] else "CCA_output"

set.seed(6545)
tab_dir  <- file.path(out_dir, "tables")
plot_dir <- file.path(out_dir, "plots")
dir.create(tab_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

save_tab <- function(x, name, digits = 4, row_names = TRUE) {
  if (is.matrix(x) || is.data.frame(x)) {
    num <- vapply(as.data.frame(x), is.numeric, logical(1))
    x <- as.data.frame(x)
    x[num] <- lapply(x[num], round, digits)
  }
  write.csv(x, file.path(tab_dir, name), row.names = row_names)
}

###############################################################################
# 0. VARIABLE SETS
#    Set 1  X(1): financial capacity / creditworthiness   (p = 3, the smaller set)
#    Set 2  X(2): spending & card-usage behavior          (q = 11)
#
#  Variables deliberately EXCLUDED to avoid exact / near-exact linear
#  dependencies (which make S11 or S22 singular, see EDA checks):
#    Monthly_Spending      = sum of the 7 category spendings (exact identity)
#    Avg_Transaction_Value ~ Monthly_Spending / Monthly_Transactions
#    Reward_Points_Earned  ~ r = .99 with spending (a proxy for it)
#    Balance / payment variables (Outstanding_Balance, Statement_Balance,
#      Payment_Amount, Payment_Ratio, Credit_Utilization) are ratio-defined
#      from each other and from Credit_Limit; they describe repayment, not
#      spending, and are kept out of this analysis.
###############################################################################
set1 <- c("Annual_Income", "Credit_Limit", "Credit_Score")
set2 <- c("Online_Shopping_Spending", "Grocery_Spending", "Fuel_Spending",
          "Dining_Spending", "Travel_Spending", "Entertainment_Spending",
          "Utility_Bill_Spending", "Monthly_Transactions", "Mobile_App_Login",
          "International_Transactions", "EMI_Count")

# right-skewed monetary variables -> Box-Cox, using the same lambda-estimation
# procedure as the EDA section (Dejene's code). Both sets are transformed so
# that the relationship between the sets stays linear (see write-up).
# Not transformed: Credit_Score, Monthly_Transactions, Mobile_App_Login,
# EMI_Count (near-symmetric); International_Transactions (count with zeros).
bc_vars <- c("Annual_Income", "Credit_Limit",
              "Online_Shopping_Spending", "Grocery_Spending", "Fuel_Spending",
              "Dining_Spending", "Travel_Spending", "Entertainment_Spending",
              "Utility_Bill_Spending")

###############################################################################
# 1. PREPROCESSING
###############################################################################
raw <- read.csv(data_csv, stringsAsFactors = FALSE)
vars <- c(set1, set2)
stopifnot(all(vars %in% names(raw)))
dat <- raw[, vars]

# 1a. data-quality checks
n_raw      <- nrow(dat)
n_missing  <- sum(!complete.cases(dat))
n_dup      <- sum(duplicated(raw[, setdiff(names(raw), "Customer_ID")]))
n_negative <- sum(dat < 0, na.rm = TRUE)
dat <- dat[complete.cases(dat), ]

# 1b. Box-Cox transformation (identical to the EDA code)
#     x^(lambda) = (x^lambda - 1) / lambda,  lambda = 0 -> log(x)
#     lambda chosen by maximizing the profile likelihood: minimizing the log
#     variance of the transformed, geometric-mean-scaled data (centered logs).
box_cox <- function(log_values, lambda) {
  if (lambda == 0) return(log_values)
  expm1(lambda * log_values) / lambda
}
estimate_lambda <- function(x) {
  stopifnot(all(x > 0))                       # Box-Cox needs strictly positive data
  centered_logs <- log(x) - mean(log(x))
  loss   <- function(power) log(var(box_cox(centered_logs, power)))
  best   <- optimize(loss, interval = c(-2, 2))$minimum
  powers <- c(-2, 0, best, 2)
  powers[which.min(vapply(powers, loss, numeric(1)))]
}

skew <- function(x) { m <- mean(x); s <- sd(x); mean((x - m)^3) / s^3 }
prep    <- dat
lambdas <- setNames(rep(NA_real_, length(vars)), vars)
for (v in bc_vars) {
  lambdas[v] <- estimate_lambda(dat[[v]])
  prep[[v]]  <- box_cox(log(dat[[v]]), lambdas[v])
}
skew_tab <- data.frame(
  variable    = vars,
  set         = ifelse(vars %in% set1, "X1 (financial)", "X2 (behavior)"),
  transformed = ifelse(vars %in% bc_vars, "Box-Cox", "none"),
  lambda      = lambdas,
  skew_before = sapply(dat[vars],  skew),
  skew_after  = sapply(prep[vars], skew))
save_tab(skew_tab, "01_preprocessing_boxcox_lambda_skewness.csv", digits = 4,
         row_names = FALSE)

# 1c. standardize (CCA on the correlation matrix R; results 10-30)
Z  <- scale(as.matrix(prep[vars]))
Z1 <- Z[, set1, drop = FALSE]
Z2 <- Z[, set2, drop = FALSE]
n  <- nrow(Z); p <- length(set1); q <- length(set2)
stopifnot(p <= q)

# random subsample of 5000 customers, used only to keep scatterplots readable
idx <- sort(sample(n, min(n, 5000)))

###############################################################################
# 2. CCA ENGINE  (Result 10.2 with R in place of S, standardized variables)
###############################################################################
mat_pow <- function(M, pow) {                # symmetric matrix power, e.g. -1/2
  e <- eigen(M, symmetric = TRUE)
  e$vectors %*% diag(e$values^pow, nrow(M)) %*% t(e$vectors)
}

cca_fit <- function(Z1, Z2) {
  R   <- cor(cbind(Z1, Z2))
  p   <- ncol(Z1); q <- ncol(Z2)
  i1  <- 1:p; i2 <- (p + 1):(p + q)
  R11 <- R[i1, i1]; R22 <- R[i2, i2]; R12 <- R[i1, i2]; R21 <- t(R12)
  R11h <- mat_pow(R11, -1/2); R22h <- mat_pow(R22, -1/2)

  # eigen-problem of R11^{-1/2} R12 R22^{-1} R21 R11^{-1/2}   (p x p)
  M  <- R11h %*% R12 %*% solve(R22) %*% R21 %*% R11h
  eM <- eigen(M, symmetric = TRUE)
  rho2 <- pmax(eM$values, 0); rho <- sqrt(rho2)

  # eigen-problem of R22^{-1/2} R21 R11^{-1} R12 R22^{-1/2}   (q x q)
  N  <- R22h %*% R21 %*% solve(R11) %*% R12 %*% R22h
  eN <- eigen(N, symmetric = TRUE)

  A <- t(R11h %*% eM$vectors)          # rows = a_k'   (p x p)
  B <- t(R22h %*% eN$vectors)          # rows = b_k'   (q x q)

  # sign convention: (i) U_k positively correlated with its own set on balance,
  # (ii) Corr(U_k, V_k) = +rho_k
  for (k in 1:p) {
    if (sum(A[k, ] %*% R11) < 0) A[k, ] <- -A[k, ]
    if (drop(A[k, ] %*% R12 %*% B[k, ]) < 0) B[k, ] <- -B[k, ]
  }
  for (k in seq_len(q)[-(1:p)])
    if (sum(B[k, ] %*% R22) < 0) B[k, ] <- -B[k, ]

  dimnames(A) <- list(paste0("U", 1:p), colnames(Z1))
  dimnames(B) <- list(paste0("V", 1:q), colnames(Z2))
  list(R = R, R11 = R11, R22 = R22, R12 = R12, R21 = R21,
       rho = rho, rho2 = rho2, A = A, B = B, p = p, q = q, n = nrow(Z1))
}

fit <- cca_fit(Z1, Z2)
rho <- fit$rho; A <- fit$A; B <- fit$B
R11 <- fit$R11; R22 <- fit$R22; R12 <- fit$R12; R21 <- fit$R21

# sanity check against base R's stats::cancor()
rho_cancor <- cancor(Z1, Z2)$cor
stopifnot(max(abs(rho - rho_cancor)) < 1e-8)

# partitioned correlation matrix R = [R11 R12; R21 R22]
save_tab(fit$R, "02_correlation_matrix_R_partitioned.csv")

# canonical variate scores (unit sample variance, uncorrelated, 10-25/10-26)
U <- Z1 %*% t(A)
V <- Z2 %*% t(B[1:p, , drop = FALSE])
colnames(V) <- paste0("V", 1:p)
score_check <- round(cor(cbind(U, V)), 6)
save_tab(score_check, "03_check_corr_of_canonical_scores.csv", digits = 6)

###############################################################################
# 3. CANONICAL CORRELATIONS & COEFFICIENTS
###############################################################################
cc_tab <- data.frame(pair = paste0("(U", 1:p, ", V", 1:p, ")"),
                     canonical_corr = rho, squared_canonical_corr = rho^2)
save_tab(cc_tab, "04_canonical_correlations.csv", row_names = FALSE)

save_tab(t(A),        "05_standardized_coefficients_financial_capacity_U.csv")   # columns = U_k
save_tab(t(B[1:p, ]), "06_standardized_coefficients_spending_behavior_V.csv")   # columns = V_k

###############################################################################
# 4. CANONICAL STRUCTURE: LOADINGS & CROSS-LOADINGS   (10-15 / 10-29)
###############################################################################
R_U_z1 <- A %*% R11                      # U with own set    (p x p)
R_V_z2 <- B[1:p, ] %*% R22               # V with own set    (p x q)
R_U_z2 <- A %*% R12                      # U with other set  (cross-loadings)
R_V_z1 <- B[1:p, ] %*% R21               # V with other set  (cross-loadings)

save_tab(t(R_U_z1), "07_loadings_financial_capacity_variables_with_U.csv")
save_tab(t(R_V_z2), "08_loadings_spending_behavior_variables_with_V.csv")
save_tab(t(R_V_z1), "09_crossloadings_financial_capacity_variables_with_V.csv")
save_tab(t(R_U_z2), "10_crossloadings_spending_behavior_variables_with_U.csv")

###############################################################################
# 5. PROPORTION OF EXPLAINED VARIANCE (10-37) AND REDUNDANCY
###############################################################################
prop1 <- rowSums(R_U_z1^2) / p           # by U_i, of set 1
prop2 <- rowSums(R_V_z2^2) / q           # by V_i, of set 2
# fin = financial-capacity set, spend = spending-behavior set
# own set  (10-37): fin_by_U, spend_by_V      redundancy: fin_by_V, spend_by_U
# _cum = cumulative over pairs 1..k, as in (10-37)
var_tab <- data.frame(
  pair           = paste0("(U", 1:p, ", V", 1:p, ")"),
  canonical_corr = rho,
  fin_by_U       = prop1,          fin_by_U_cum   = cumsum(prop1),
  spend_by_V     = prop2,          spend_by_V_cum = cumsum(prop2),
  fin_by_V       = prop1 * rho^2,  fin_by_V_cum   = cumsum(prop1 * rho^2),
  spend_by_U     = prop2 * rho^2,  spend_by_U_cum = cumsum(prop2 * rho^2))
save_tab(var_tab, "11_variance_explained_and_redundancy.csv", row_names = FALSE)

###############################################################################
# 6. MATRICES OF ERRORS OF APPROXIMATION (10-34), standardized version
###############################################################################
Ainv <- solve(A)                          # columns a_z^(i)
Binv <- solve(B)                          # columns b_z^(i)
err_summary <- NULL
for (r in 1:p) {
  E12 <- R12; E11 <- R11; E22 <- R22
  for (i in 1:r) {
    E12 <- E12 - rho[i] * Ainv[, i] %o% Binv[, i]
    E11 <- E11 - Ainv[, i] %o% Ainv[, i]
    E22 <- E22 - Binv[, i] %o% Binv[, i]
  }
  dimnames(E12) <- dimnames(R12); dimnames(E11) <- dimnames(R11)
  dimnames(E22) <- dimnames(R22)
  save_tab(E12, sprintf("12_error_R12_first%d_pairs.csv", r))
  save_tab(E11, sprintf("13_error_R11_first%d_pairs.csv", r))
  save_tab(E22, sprintf("14_error_R22_first%d_pairs.csv", r))
  err_summary <- rbind(err_summary, data.frame(
    r = r,
    max_abs_E12 = max(abs(E12)), max_abs_E11 = max(abs(E11)),
    max_abs_E22 = max(abs(E22)),
    trace_E11_over_p = sum(diag(E11)) / p,     # = 1 - R2_{z1 | U1..Ur}
    trace_E22_over_q = sum(diag(E22)) / q))    # = 1 - R2_{z2 | V1..Vr}
}
save_tab(err_summary, "15_error_matrices_summary.csv", row_names = FALSE)

###############################################################################
# 7. LARGE-SAMPLE INFERENCE: BARTLETT LR TESTS (10-39) & SEQUENTIAL (10-41)
###############################################################################
bart <- n - 1 - (p + q + 1) / 2
test_tab <- do.call(rbind, lapply(0:(p - 1), function(k) {
  stat <- -bart * sum(log(1 - rho[(k + 1):p]^2))
  df   <- (p - k) * (q - k)
  data.frame(
    H0 = if (k == 0) "Sigma12 = 0 (all rho* = 0)"
         else sprintf("only first %d rho* nonzero (rho*_%d..rho*_%d = 0)", k, k + 1, p),
    wilks_lambda = prod(1 - rho[(k + 1):p]^2),
    bartlett_chisq = stat, df = df,
    crit_05 = qchisq(.95, df), crit_01 = qchisq(.99, df),
    p_value = pchisq(stat, df, lower.tail = FALSE),
    decision_01 = ifelse(stat > qchisq(.99, df), "Reject H0", "Do not reject H0"))
}))
save_tab(test_tab, "16_bartlett_sequential_tests.csv", digits = 6, row_names = FALSE)

###############################################################################
# 8. STABILITY CHECKS (Section 10.4 advice)
#    (a) split-sample validation  (b) outlier sensitivity
###############################################################################
half  <- sample(n) <= n / 2
fitA  <- cca_fit(Z1[half, ],  Z2[half, ])
fitB  <- cca_fit(Z1[!half, ], Z2[!half, ])
loadA <- t(fitA$A %*% fitA$R11); loadB <- t(fitB$A %*% fitB$R11)
# loading_similarity = Tucker congruence of the U-loadings in half 1 vs half 2
# (1 = identical pattern;
split_tab <- data.frame(
  pair = paste0("(U", 1:p, ", V", 1:p, ")"),
  canonical_corr_full_sample = rho,
  canonical_corr_half1       = fitA$rho,
  canonical_corr_half2       = fitB$rho,
  loading_similarity_half1_vs_half2 = sapply(1:p, function(k)
    abs(sum(loadA[, k] * loadB[, k])) / sqrt(sum(loadA[, k]^2) * sum(loadB[, k]^2))))
save_tab(split_tab, "17_split_half_validation.csv", row_names = FALSE)

# (b) outlier sensitivity: rerun the CCA without customers whose squared
#     generalized distance exceeds the chi-square(p+q) 99th percentile
d2   <- mahalanobis(Z, colMeans(Z), cov(Z))
keep <- d2 <= qchisq(0.99, df = p + q)
fit_out  <- cca_fit(Z1[keep, ], Z2[keep, ])
sens_tab <- data.frame(
  pair = paste0("(U", 1:p, ", V", 1:p, ")"),
  canonical_corr_all_customers     = rho,
  canonical_corr_without_flagged   = fit_out$rho,
  n_customers_removed              = sum(!keep))
save_tab(sens_tab, "18_outlier_sensitivity.csv", row_names = FALSE)

###############################################################################
# 9. PLOTS
###############################################################################
pal <- c("#2B6CB0", "#C05621", "#2F855A", "#6B46C1")

# 9b. scatterplots of canonical variate pairs (random 5000 points)

draw_scatter <- function() {
  par(mfrow = c(1, p), mar = c(4.2, 4.2, 2.2, 0.8), oma = c(0, 0, 0, 0),
      mgp = c(2.4, 0.6, 0), tcl = -0.3, las = 1, family = "sans",
      cex.axis = 0.85, cex.lab = 0.95, pty = "s")
  for (k in 1:p) {
    lim <- extendrange(c(U[idx, k], V[idx, k]), f = 0.03)   # same range on x and y
    plot(U[idx, k], V[idx, k], type = "n", xlim = lim, ylim = lim, axes = FALSE,
         xlab = bquote("Financial-capacity variate " * hat(U)[.(k)]),
         ylab = bquote("Spending-behavior variate " * hat(V)[.(k)]))
    abline(h = 0, v = 0, col = "grey88", lwd = 0.8)
    points(U[idx, k], V[idx, k], pch = 16, cex = 0.3, col = rgb(0.15, 0.30, 0.55, 0.30))
    fit <- coef(lm(V[, k] ~ U[, k]))                         # fitted line, clipped
    xr  <- range(U[idx, k])                                  # to the data range
    segments(xr[1], fit[1] + fit[2] * xr[1], xr[2], fit[1] + fit[2] * xr[2],
             col = "black", lwd = 1.4)
    axis(1, col = "grey30"); axis(2, col = "grey30"); box(bty = "l", col = "grey30")
    mtext(sprintf("(%s) Pair %d", letters[k], k), side = 3, line = 0.6,
          adj = 0, font = 2, cex = 0.9)
    text(lim[1], lim[2], bquote(italic(r) == .(sprintf("%.3f", rho[k]))),
         adj = c(0, 1), cex = 0.9)
  }
}
png(file.path(plot_dir, "01_canonical_variate_scatterplots.png"),
    width = 9, height = 3.3, units = "in", res = 300)
draw_scatter(); invisible(dev.off())
pdf(file.path(plot_dir, "01_canonical_variate_scatterplots.pdf"),
    width = 9, height = 3.3, useDingbats = FALSE)
draw_scatter(); invisible(dev.off())

# 9c. structure loadings (own-set correlations) for each pair

pretty_name <- function(x) {
  map <- c(Annual_Income = "Annual income", Credit_Limit = "Credit limit",
           Credit_Score = "Credit score", Online_Shopping_Spending = "Online shopping",
           Grocery_Spending = "Grocery", Fuel_Spending = "Fuel",
           Dining_Spending = "Dining", Travel_Spending = "Travel",
           Entertainment_Spending = "Entertainment", Utility_Bill_Spending = "Utility bills",
           Monthly_Transactions = "Monthly transactions", Mobile_App_Login = "Mobile app logins",
           International_Transactions = "International transactions", EMI_Count = "EMI count")
  ifelse(x %in% names(map), map[x], gsub("_", " ", x))
}
draw_loading_panel <- function(load, letter, set_label, sym, k, col_pos, col_neg,
                               show_names = TRUE) {
  v   <- rev(load)
  xl  <- c(-1.3, 1.3)                       # room for value labels beyond +/-1
  nms <- if (show_names) rev(pretty_name(names(load))) else rep("", length(v))
  # empty frame first, then guide lines, then bars on top
  y <- barplot(v, horiz = TRUE, xlim = xl, axes = FALSE, col = NA, border = NA,
               names.arg = nms, las = 1, cex.names = 0.8, space = 0.35)
  abline(v = c(-1, -0.5, 0.5, 1), col = "grey90", lwd = 0.7)
  barplot(v, horiz = TRUE, add = TRUE, axes = FALSE, names.arg = rep("", length(v)),
          border = NA, col = ifelse(v >= 0, col_pos, col_neg), space = 0.35)
  abline(v = 0, col = "grey30", lwd = 0.8)
  axis(1, at = seq(-1, 1, 0.5), col = "grey30", cex.axis = 0.75)
  off <- 0.05
  pos <- v >= 0     # text() takes one alignment per call: draw +/- separately
  if (any(pos))  text(v[pos] + off,  y[pos],  sprintf("%.2f", v[pos]),
                      adj = c(0, 0.5), cex = 0.65, col = "grey20", xpd = TRUE)
  if (any(!pos)) text(v[!pos] - off, y[!pos], sprintf("%.2f", v[!pos]),
                      adj = c(1, 0.5), cex = 0.65, col = "grey20", xpd = TRUE)
  title_txt <- sprintf("(%s) %s", letter, set_label)
  mtext(bquote(bold(.(title_txt)) %->% bold(hat(.(as.name(sym)))[.(k)])),
        side = 3, line = 0.5, adj = 0, cex = 0.75)
}
draw_loadings <- function() {
  # variable names only in the first column; the first column is wider so
  # all panels have the same plotting width.
  # top row (3 variables) shorter than bottom row (11 variables)
  name_lines <- 9.5                         # left margin (lines) for the names
  lay <- matrix(1:(2 * p), 2, byrow = TRUE)
  hts <- c(length(set1) + 1.5, length(set2) + 1.5)
  # layout() shrinks the text size, so measure the margin-line height after it
  layout(lay, heights = hts)
  fig_w <- par("din")[1]; csi <- par("csi") # device width (in), inches per margin line
  plot_w <- (fig_w - csi * (name_lines + 1.2) - (p - 1) * csi * (1 + 1.2)) / p
  col_w  <- c(csi * (name_lines + 1.2) + plot_w, rep(csi * 2.2 + plot_w, p - 1))
  layout(lay, widths = col_w, heights = hts)
  par(mgp = c(2, 0.5, 0), tcl = -0.25, family = "sans", oma = c(1.5, 0, 0, 0))
  mar_for <- function(k) c(2.5, if (k == 1) name_lines else 1, 2.2, 1.2)
  for (k in 1:p) {
    par(mar = mar_for(k))
    draw_loading_panel(R_U_z1[k, ], letters[k], "Financial", "U", k,
                       "#2C4A75", "#B5653B", show_names = (k == 1))
  }
  for (k in 1:p) {
    par(mar = mar_for(k))
    draw_loading_panel(R_V_z2[k, ], letters[p + k], "Spending", "V", k,
                       "#2C4A75", "#B5653B", show_names = (k == 1))
  }
  mtext("Canonical loading (correlation between variable and canonical variate)",
        side = 1, outer = TRUE, line = 0.3, cex = 0.75)
}
png(file.path(plot_dir, "02_structure_loadings.png"),
    width = 6.5, height = 4.5, units = "in", res = 300)
draw_loadings(); invisible(dev.off())
pdf(file.path(plot_dir, "02_structure_loadings.pdf"),
    width = 6.5, height = 4.5, useDingbats = FALSE)
draw_loadings(); invisible(dev.off())

###############################################################################
# 10. SAVE SCORES AND A TEXT SUMMARY
###############################################################################
scores <- data.frame(Customer_ID = raw$Customer_ID[complete.cases(raw[, vars])],
                     round(U, 5), round(V, 5))
write.csv(scores, file.path(out_dir, "canonical_variate_scores.csv"), row.names = FALSE)

sink(file.path(out_dir, "CCA_summary.txt"))
cat("CANONICAL CORRELATION ANALYSIS - credit-card customer data\n")
cat("Johnson & Wichern, Chapter 10\n", strrep("=", 70), "\n\n", sep = "")
cat(sprintf("Rows read: %d | missing removed: %d | duplicate rows: %d | negative values: %d\n",
            n_raw, n_missing, n_dup, n_negative))
cat(sprintf("n = %d, p = %d (set 1), q = %d (set 2)\n", n, p, q))
cat("Set 1 X(1):", paste(set1, collapse = ", "), "\n")
cat("Set 2 X(2):", paste(set2, collapse = ", "), "\n")
cat("Box-Cox applied to (lambda):\n")
cat(sprintf("  %-26s %.4f\n", bc_vars, lambdas[bc_vars]), sep = "")
cat("All variables standardized -> analysis on correlation matrix R.\n")
cat("\n")

cat("Canonical correlations (checked against stats::cancor):\n")
print(round(cc_tab[, -1], 4), row.names = FALSE)
cat("\nStandardized coefficients, set 1 (columns = U_k):\n"); print(round(t(A), 3))
cat("\nStandardized coefficients, set 2 (columns = V_k):\n"); print(round(t(B[1:p, ]), 3))
cat("\nLoadings: corr(U_k, X1):\n"); print(round(t(R_U_z1), 3))
cat("\nLoadings: corr(V_k, X2):\n"); print(round(t(R_V_z2), 3))
cat("\nCross-loadings: corr(U_k, X2):\n"); print(round(t(R_U_z2), 3))
cat("\nCross-loadings: corr(V_k, X1):\n"); print(round(t(R_V_z1), 3))
cat("\nProportion of variance explained & redundancy:\n")
print(round(var_tab[, -1], 4), row.names = FALSE)
cat("\nErrors of approximation (max |entry|):\n")
print(round(err_summary, 4), row.names = FALSE)
cat("\nBartlett sequential tests:\n")
print(test_tab, row.names = FALSE, digits = 4)
cat("\nSplit-half validation:\n")
print(cbind(split_tab[1], round(split_tab[-1], 4)), row.names = FALSE)
cat("\nOutlier sensitivity:\n")
print(cbind(sens_tab[1], round(sens_tab[-1], 4)), row.names = FALSE)
cat("\nNOTE: with n =", n, "even trivially small canonical correlations are\n",
    "'significant' (cf. Example 10.8). Judge practical importance by the size\n",
    "of rho*, the variance explained, and redundancy, not the p-values alone.\n")
sink()

cat("Done. Results written to:", normalizePath(out_dir), "\n")