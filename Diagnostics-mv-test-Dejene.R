
###***************************************************
#
#
#          Data diagnostics
#
#
#*******************************************************
#*
#*
#*
# The following code is fully working on my computer 
# In order to get the corresponding result on your computer you need to have your won corresponding folder and 
# save the data into that folder and initialize the code'
#
# The kernel of the code is the lecture note and R code from Dr. Paige 
# Some modification were made to make our plots attractive and handy for dispaying

# Data and flder preparation

DATA_FILE <- "C:/Users/ddcdq7/OneDrive - University of Missouri/Fall2026/MVA/Project/synthetic_credit_card_customer_behavior_dataset.csv"
out <- file.path(dirname(DATA_FILE), "Proj-results-dejene", "diagnostics results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("writexl", quietly = TRUE)) {
  stop("Install the Excel package first: install.packages('writexl')")
}
for (folder in c("PNG", "PDF", "Excel", "Text", "CSV", "LaTeX")) {
  dir.create(file.path(out, folder), recursive = TRUE, showWarnings = FALSE)
}

# Keep the CSV and save the same table as an Excel workbook.
save_result_table <- function(x, file, row.names = TRUE) {
  table <- as.data.frame(x)
  utils::write.csv(table, file, row.names = row.names)
  if (row.names) table <- data.frame(Row = rownames(table), table, check.names = FALSE)
  workbook <- paste0(tools::file_path_sans_ext(basename(file)), ".xlsx")
  writexl::write_xlsx(table, file.path(out, "Excel", workbook))
}

EXTRA_DIAGNOSTICS <- TRUE # Produce both original diagnostic PDFs AND report figures
set.seed(123)
alpha <- .05
B <- 1999L
vars <- c("Age", "Annual_Income", "Credit_Limit", "Monthly_Spending",
          "Credit_Utilization", "Credit_Score")
money <- c("Annual_Income", "Credit_Limit", "Monthly_Spending")
save_csv <- function(x, name) save_result_table(x, file.path(out, "CSV", name), row.names = FALSE)
d <- read.csv(DATA_FILE, stringsAsFactors = FALSE, check.names = TRUE,
              na.strings = c("", "NA", "N/A"))
stopifnot(all(vars %in% names(d)), all(vapply(d[vars], is.numeric, logical(1))))

# 1. Audit; do not equate duplicate IDs with erroneous records.
audit <- data.frame(Variable = names(d),
  Type = vapply(d, function(x) class(x)[1], character(1)),
  Missing = vapply(d, function(x) sum(is.na(x)), numeric(1)),
  Nonfinite_nonmissing = vapply(d, function(x)
    if (is.numeric(x)) sum(!is.finite(x) & !is.na(x)) else 0, numeric(1)))
save_csv(audit, "Variable_Audit.csv")
dup <- duplicated(d) | duplicated(d, fromLast = TRUE)
ids <- intersect(c("Customer_ID", "Customer.ID"), names(d))
dup_id <- rep(FALSE, nrow(d))
for (v in ids) dup_id <- dup_id | (!is.na(d[[v]]) &
  (duplicated(d[[v]]) | duplicated(d[[v]], fromLast = TRUE)))
save_csv(data.frame(Original_Row = seq_len(nrow(d)), Duplicate_record = dup,
                    Duplicate_ID = dup_id), "Duplicate_Review.csv")
num <- names(d)[vapply(d, is.numeric, logical(1))]
save_csv(do.call(rbind, lapply(num, function(v) {
  x <- d[[v]][is.finite(d[[v]])]
  data.frame(Variable = v, Minimum = if(length(x)) min(x) else NA_real_,Maximum = if(length(x)) max(x) else NA_real_, Negative = sum(x < 0))
})), "Range_Review.csv")
# Verify ranges against the data dictionary; negative balances may be valid.
# No undocumented range correction, imputation or duplicate removal is applied.
ok <- apply(as.matrix(d[vars]), 1, function(x) all(is.finite(x)))
save_csv(data.frame(Original_Row = seq_len(nrow(d)), Included = ok), "Case_Audit.csv")
X <- d[ok, vars, drop = FALSE]
rows <- which(ok)
stopifnot(nrow(X) > ncol(X) + 2L)
keep <- vapply(X, sd, numeric(1)) > 0
save_csv(data.frame(Variable = names(X), Retained = keep), "Constant_Columns.csv")
X <- X[, keep, drop = FALSE]
stopifnot(ncol(X) >= 2L)
if (qr(scale(X))$rank < ncol(X)) stop("Resolve linearly dependent columns first.")
original_data <- X
analysis_variables <- names(X)
monetary <- intersect(money, analysis_variables)
if (any(as.matrix(X[monetary]) < 0)) stop("Check the monetary columns for negative values.")
transformed_data <- original_data
transformations <- data.frame(Variable = analysis_variables, Method = "Original", Lambda = NA_real_,
                              stringsAsFactors = FALSE)

# transformation

box_cox <- function(log_values, lambda) { if (lambda == 0) return(log_values)
  expm1(lambda * log_values) / lambda}

for (variable in monetary) {
  values <- original_data[[variable]]
  index <- match(variable, analysis_variables)
  if (any(values == 0) || min(values) == max(values)) {
    transformations$Method[index] <- if (any(values == 0)) "Original (contains zero)" else "Original (constant)"
    next
  }

  log_values <- log(values)
  # Centering the logs incorporates the Jacobian in the profile likelihood.
  centered_logs <- log_values - mean(log_values)
  loss <- function(power) log(var(box_cox(centered_logs, power)))
  best_power <- optimize(loss, interval = c(-2, 2))$minimum
  powers <- c(-2, 0, best_power, 2)
  lambda <- powers[which.min(vapply(powers, loss, numeric(1)))]
  transformed_data[[variable]] <- box_cox(log_values, lambda)
  transformations$Method[index] <- "Box-Cox"
  transformations$Lambda[index] <- lambda
}

save_csv(transformations, "Variable_Transformations.csv")
save_csv(data.frame(Original_Row = rows, transformed_data), "BoxCox_Selected_Data.csv")

# One fixed random subsample for tests on both scales (Shapiro limit: 5000).
ii <- sort(sample.int(nrow(X), min(5000L, nrow(X))))
save_csv(data.frame(Original_Row = rows[ii]), "Normality_Test_Rows.csv")
q <- qnorm(ppoints(length(ii), a = .5))
null_rQ <- replicate(B, cor(sort(rnorm(length(ii))), q))
plot_rows <- sort(sample.int(nrow(X), min(2000L, nrow(X))))

# 2. Screening on the common rows and variable set defined above.
screen <- function(x, label) {
  n <- nrow(x); p <- ncol(x)
  Z <- scale(x)
  if (rcond(cov(Z)) < 1e-12) stop("Ill-conditioned covariance: review redundancy.")
  tests <- do.call(rbind, lapply(names(x), function(v) {
    y <- x[[v]][ii]
    if (sd(y) == 0) return(data.frame(Variable=v, rQ=NA, rQ_p=NA, SW_p=NA))
    r <- cor(sort(y), q)
    data.frame(Variable = v, rQ = r,
      rQ_p = (1 + sum(null_rQ <= r))/(B + 1), SW_p = shapiro.test(y)$p.value)
  }))
  tests$rQ_Bonferroni <- p.adjust(tests$rQ_p, "bonferroni")
  tests$SW_Bonferroni <- p.adjust(tests$SW_p, "bonferroni")
  tests$n_test <- length(ii)
  save_csv(tests, paste0(label, "_Normality_Tests.csv"))
  D2 <- mahalanobis(Z, colMeans(Z), cov(Z))
  iqr_flag <- sapply(x, function(y) {
    b <- quantile(y, c(.25, .75)); w <- b[2] - b[1]
    y < b[1] - 1.5*w | y > b[2] + 1.5*w
  })
  save_csv(data.frame(Original_Row = rows, iqr_flag), paste0(label, "_IQR_Flags.csv"))
  flags <- data.frame(Original_Row = rows, D2 = D2,
    Flag_99pct = D2 > qchisq(.99, p),
    Flag_Bonferroni = D2 > qchisq(1-alpha/n, p))
  save_csv(flags, paste0(label, "_Mahalanobis_Flags.csv"))
  # These chi-square thresholds are approximate screens, not exact outlier tests.
  sets <- c(combn(names(x), 2, simplify = FALSE),
    if(p > 3) list(names(x)[1:3]) else list(), list(names(x)))
  radial <- do.call(rbind, lapply(sets, function(v) {
    A <- Z[, v, drop = FALSE]
    ds <- mahalanobis(A, colMeans(A), cov(A))
    data.frame(Variables = paste(v, collapse = "; "),
      Proportion_inside_50pct = mean(ds <= qchisq(.5, length(v))))}))
  save_csv(radial, paste0(label, "_Central_50pct_Check.csv"))
  if (EXTRA_DIAGNOSTICS) {
    save_figure(paste0(label, "_Diagnostics"), 11, 8.5, function() {
      par(mfrow = c(2, 4), mar = c(4, 4, 3.5, 1),
          cex.main = 1, cex.lab = .85, cex.axis = .8)
      for (v in names(x)) {
        hist(x[[v]], breaks = "FD", main = v, xlab = "", col = "grey80")
        qqnorm(x[[v]], main = v, pch = 16, cex = .2); qqline(x[[v]], col = 2)
      }
      par(mfrow = c(1, 1))
      pairs(Z[plot_rows, ], pch = 16, cex = .25,
            col = ifelse(flags$Flag_99pct[plot_rows], "red", "grey50"))
      par(mfrow = c(2, 3), mar = c(4, 4, 3.5, 1),
          cex.main = 1, cex.lab = .85, cex.axis = .8)
      for (v in sets) {
        A <- Z[, v, drop = FALSE]
        ds <- mahalanobis(A, colMeans(A), cov(A))
        plot(qchisq(ppoints(n, a = .5), length(v)), sort(ds), pch = 16, cex = .2,
          main = paste(v, collapse = " / "), xlab = "Chi-square quantiles", ylab = "Ordered D-squared")
        abline(0, 1, col = 2)
      }
    }, multipage = TRUE)
  }
  invisible(list(tests = tests, flags = flags, iqr = iqr_flag, radial = radial))}
save_figure <- function(name, width, height, draw, multipage = FALSE) {
  for (ext in c("pdf", "png")) {
    suffix <- if (multipage && ext == "png") "_%02d" else ""
    path <- file.path(out, toupper(ext), paste0(name, suffix, ".", ext))
    if (ext == "pdf") pdf(path, width = width, height = height,
                          useDingbats = FALSE, pointsize = 14)
    else png(path, width = width, height = height, units = "in", res = 300, pointsize = 14)
    tryCatch(draw(), finally = dev.off())  }}

original <- screen(X, "Original")
transformed <- screen(transformed_data, "BoxCox")
raw_flags <- original$flags
transformed_flags <- transformed$flags
cat("Results saved in:", out, "\n")


# Plots

blue <- "#195D9C"
red <- "#B83B36"
scales <- list(Original = original_data, Transformed = transformed_data)
variables <- names(original_data)
pages <- split(variables, ceiling(seq_along(variables) / 8))
normal_quantiles <- qnorm(ppoints(nrow(original_data), a = .5))

draw_comparison <- function(page_variables, kind) {
  par(mfrow = c(4, 4), mar = c(3.5, 4, 3.6, 1),
      mgp = c(2.3, .7, 0), cex.main = 1, cex.lab = .85, cex.axis = .8)
  for (scale_name in names(scales)) {
    for (slot in seq_len(8)) {
      if (slot > length(page_variables)) {
        plot.new()
        next
      }
      variable <- page_variables[slot]
      values <- scales[[scale_name]][[variable]]
      title <- paste(scale_name, gsub("_", " ", variable), sep = "\n")
      if (kind == "histograms") {
        hist(values, breaks = "FD", main = title, xlab = "Value",
             ylab = "Frequency", col = blue, border = "white")
      } else {
        plot(normal_quantiles, sort(values), main = title, pch = 16,
             cex = .22, col = blue, xlab = "Normal quantile", ylab = "Sample quantile")
        qqline(values, col = red, lwd = 1.3)
      }
    }
  }
}

for (kind in c("histograms", "normal_qq")) {
  for (page in seq_along(pages)) {
    name <- if (length(pages) == 1L) kind else sprintf("%s_%02d", kind, page)
    save_figure(name, 14, 12, function() draw_comparison(pages[[page]], kind))
  }
}
save_figure("mahalanobis_qq", 10, 4.2, function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1),
      mgp = c(2.8, .8, 0), cex.main = 1, cex.lab = .85, cex.axis = .8)
  theoretical <- qchisq(ppoints(nrow(X), a = .5), ncol(X))
  for (label in c("Original", "Transformed")) {
    distances <- if (label == "Original") raw_flags$D2 else transformed_flags$D2
    plot(theoretical, sort(distances), pch = 16, cex = .25, col = blue,
         main = paste(label, "scale"),
         xlab = paste0("Chi-square quantile (df = ", ncol(X), ")"),
         ylab = "Ordered squared Mahalanobis distance")
    abline(0, 1, col = red, lwd = 1.3)  }})
# save
report_table <- function(tab, name, caption, label, digits = 4) {
  save_result_table(tab, file.path(out, "CSV", paste0(name, ".csv")), row.names = FALSE)
  display <- lapply(tab, function(x) {
    if (is.numeric(x)) {
      vapply(x, function(z) {
        if (is.na(z)) return("NA")
        if (z != 0 && abs(z) < 0.001) return(formatC(z, format = "e", digits = 2))
        format(round(z, digits), trim = TRUE, scientific = FALSE)
      }, character(1))
    } else as.character(x)  })
  display <- as.data.frame(display, check.names = FALSE)
  # All report labels are controlled here; replace underscores for LaTeX.
  display[] <- lapply(display, function(z) gsub("_", " ", z, fixed = TRUE))
  lines <- c("\\begin{table}[htbp]", "\\centering", "\\small",
    paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"),
    "\\setlength{\\tabcolsep}{4pt}",
    paste0("\\begin{tabular}{l", paste(rep("r", ncol(tab)-1), collapse = ""), "}"),
    "\\toprule", paste0(paste(names(tab), collapse = " & "), " \\\\"),
    "\\midrule",
    apply(display, 1, function(z) paste0(paste(z, collapse = " & "), " \\\\")),
    "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file.path(out, "LaTeX", paste0(name, ".tex")))}
a <- original$tests; b <- transformed$tests
normality <- data.frame(Variable = a$Variable,`Original rQ` = a$rQ, `Transformed rQ` = b$rQ,
  `Original SW adj. p` = a$SW_Bonferroni, `Transformed SW adj. p` = b$SW_Bonferroni, check.names = FALSE)
report_table(normality, "Normality_Comparison", "Normality diagnostics before and after Box-Cox transformation.",
             "tab:normality")
# Full test files in the root also contain Monte Carlo rQ p-values.
normality_decisions <- data.frame(Variable = a$Variable,
  Original_rQ_reject = a$rQ_Bonferroni < alpha,
  Transformed_rQ_reject = b$rQ_Bonferroni < alpha,
  Original_SW_reject = a$SW_Bonferroni < alpha,
  Transformed_SW_reject = b$SW_Bonferroni < alpha)
save_result_table(normality_decisions, file.path(out, "CSV", "Normality_Decisions.csv"), row.names = FALSE)
joint <- data.frame(Diagnostic = "Central 50 percent: full vector",
  Original = tail(original$radial$Proportion_inside_50pct, 1),
  Transformed = tail(transformed$radial$Proportion_inside_50pct, 1))
if (ncol(X) > 3) joint <- rbind(joint, data.frame(
  Diagnostic = "Central 50 percent: first three variables",
  Original = original$radial$Proportion_inside_50pct[nrow(original$radial)-1],
  Transformed = transformed$radial$Proportion_inside_50pct[nrow(transformed$radial)-1]))
report_table(joint, "Multivariate_Comparison",
             "Joint diagnostics; central-region entries are proportions.", "tab:joint")
ni <- c(colSums(original$iqr), sum(rowSums(original$iqr) > 0))
nt <- c(colSums(transformed$iqr), sum(rowSums(transformed$iqr) > 0))
iqr_summary <- data.frame(Variable = c(names(X), "At least one variable"),
  `Original n` = ni, `Original percent` = 100*ni/nrow(X),
  `Transformed n` = nt, `Transformed percent` = 100*nt/nrow(X), check.names = FALSE)
report_table(iqr_summary, "IQR_Comparison", "Univariate outlier screening; counts overlap across variables.","tab:iqr", digits = 2)
mo <- colSums(raw_flags[c("Flag_99pct", "Flag_Bonferroni")])
mt <- colSums(transformed_flags[c("Flag_99pct", "Flag_Bonferroni")])
mahal <- data.frame(Rule = c("99 percent", "Bonferroni"),
  Cutoff = qchisq(c(.99, 1-alpha/nrow(X)), ncol(X)),
  `Original n` = mo, `Original percent` = 100*mo/nrow(X),
  `Transformed n` = mt, `Transformed percent` = 100*mt/nrow(X), check.names = FALSE)
report_table(mahal, "Mahalanobis_Comparison", "Approximate multivariate outlier screening.",
             "tab:mahalanobis", digits = 3)
audit_summary <- data.frame(Item = c("Input rows", "Retained rows", "Excluded rows",
    "Missing cells", "Nonfinite nonmissing cells", "Rows in duplicate records",
    "Rows with duplicate IDs", "Selected constant columns"),
  Count = c(nrow(d), sum(ok), sum(!ok), sum(audit$Missing),
    sum(audit$Nonfinite_nonmissing), sum(dup), sum(dup_id), sum(!keep)))
save_result_table(audit_summary, file.path(out, "CSV", "Audit_Summary.csv"), row.names = FALSE)
save_result_table(transformed_flags[order(-transformed_flags$D2), ],
          file.path(out, "CSV", "Transformed_Outliers_Ranked.csv"), row.names = FALSE)
writeLines(c("All complete observations are used for Q-Q plots and outlier screening. The scatterplot matrix uses up to 2000 rows.",
  paste("Normality test sample size:", length(ii)), paste("Normal simulations:", B),
  paste("Minimum Monte Carlo p-value:", 1/(B+1)),
  "Chi-square thresholds are approximate screens.",
  "No outliers or duplicate rows were automatically deleted.",
  "Box-Cox powers and exclusions are recorded in Variable_Transformations.csv.",
  "Repeat checks for each analysis variable set and within relevant groups."),
  file.path(out, "Text", "Analysis_Notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "Text", "R_Session_Info.txt"))
report_sections <- list("Data audit" = audit_summary, "Box-Cox powers" = transformations,
  "Normality tests" = normality, "Adjusted test decisions" = normality_decisions,
  "Multivariate diagnostics" = joint, "IQR outliers" = iqr_summary,
  "Mahalanobis outliers" = mahal)
capture.output({
  cat("Data diagnostics\n")
  for (section in names(report_sections)) {
    cat("\n", section, "\n", sep = "")
    print(report_sections[[section]], row.names = FALSE, digits = 6)
  }
}, file = file.path(out, "Text", "Diagnostics_Summary.text"))
cat("\nDiagnostics results folder:", normalizePath(out), "\n")
