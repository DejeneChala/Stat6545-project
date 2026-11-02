# =============================================================================
# STAT 6545 - Multivariate Statistical Methods (Fall 2026)
# Multivariate Profiling and Segmentation of Credit Card Customers
# Sections 4.6 and 4.7: Discriminant analysis, logistic regression, and
#                       tree-based predictive comparison
# Script author: Tobechukwu Kingsley Abor
#
# Contents
#   A. Setup, data loading, preparation
#   B. Section 4.7  Linear regression vs regression tree vs random forest vs
#                   XGBoost  (outcomes: log Monthly_Spending; Credit_Score)
#        - 5-fold CV and 30% holdout: RMSE, MAE, R2
#        - variable importance and partial-dependence plots
#   C. Section 4.6  Discriminant analysis
#        - normality and covariance checks WITHIN card types
#        - LDA / QDA for Card_Type with CV and holdout
#   D. Section 4.7  Logistic regression for Premium card (Platinum/Signature)
#        - odds ratios, LR test, calibration, ROC; compared with LDA and QDA
#
# Predictors that define or duplicate a response are excluded (see comments).
# All figures and tables are written to OUT_DIR.
# =============================================================================
setwd(r"(C:\Users\tjgfx9\OneDrive - University of Missouri\Masters\Stat6545-project)")
credit_data <- read.csv("synthetic_credit_card_customer_behavior_dataset.csv")
## A. Setup -------------------------------------------------------------------
pkgs <- c("readxl", "MASS", "rpart", "rpart.plot", "ranger", "xgboost",
          "ggplot2", "psych", "car")
need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) {
    if (is.na(getOption("repos")["CRAN"]) || getOption("repos")["CRAN"] == "@CRAN@")
        options(repos = c(CRAN = "https://cloud.r-project.org"))
    install.packages(need)
}
for (p in c("MASS", "rpart", "ranger", "xgboost", "ggplot2"))
    suppressPackageStartupMessages(library(p, character.only = TRUE))
 
DATA_PATH <- "synthetic_credit_card_customer_behavior_dataset.xlsx"  # data file
SHEET     <- 1            # sheet name or number (Excel only)
OUT_DIR   <- "output_4_6_4_7"
SEED      <- 6545
K_FOLDS   <- 5
TRAIN_FRAC <- 0.70
N_SUB     <- 3000         # subsample for O(n^2) / small-sample tests
NT        <- max(1, parallel::detectCores() - 1)   # threads for ranger/xgboost
dir.create(OUT_DIR, showWarnings = FALSE)
set.seed(SEED)
 
fig <- function(name, plot, w = 10, h = 6.5)
    ggsave(file.path(OUT_DIR, name), plot, width = w, height = h, dpi = 150, bg = "white")
save_png <- function(name, expr, w = 1600, h = 1200) {
    png(file.path(OUT_DIR, name), width = w, height = h, res = 150)
    on.exit(dev.off()); expr
}
theme_set(theme_minimal(base_size = 12) +
              theme(plot.title = element_text(face = "bold"), legend.position = "bottom"))
 
model_cols <- c("Linear regression" = "#5B7073", "Regression tree" = "#F0A53A",
                "Random forest" = "#1C7C85", "XGBoost" = "#C8553D")
raw <- credit_data
names(raw) <- gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", trimws(names(raw))))                      
names(raw)[names(raw) == "Mobile_App_Login"] <- "Mobile_App_Logins"   # dataset spelling
 
cat_vars  <- c("Gender", "Occupation", "Card_Type")
money_vars <- c("Annual_Income", "Credit_Limit", "Monthly_Spending", "Avg_Transaction_Value",
                "Online_Shopping_Spending", "Grocery_Spending", "Fuel_Spending",
                "Dining_Spending", "Travel_Spending", "Entertainment_Spending",
                "Utility_Bill_Spending", "Outstanding_Balance", "Statement_Balance",
                "Payment_Amount", "Cash_Advance_Amount", "Reward_Points_Earned",
                "Reward_Points_Redeemed")
need_cols <- c(cat_vars, money_vars, "Age", "Card_Age_Months", "Monthly_Transactions",
               "EMI_Count", "International_Transactions", "Mobile_App_Logins",
               "Payment_Ratio", "Credit_Utilization", "Credit_Score")
miss <- setdiff(need_cols, names(raw))
if (length(miss)) stop("Columns not found: ", paste(miss, collapse = ", "),
                       "\nColumns present: ", paste(names(raw), collapse = ", "))
 
d <- na.omit(raw[, need_cols])
d[cat_vars] <- lapply(d[cat_vars], factor)
tier_order <- c("Basic", "Silver", "Gold", "Platinum", "Signature")
if (all(levels(d$Card_Type) %in% tier_order))
    d$Card_Type <- factor(d$Card_Type, levels = tier_order[tier_order %in% levels(d$Card_Type)])
d[money_vars] <- lapply(d[money_vars], log1p)   # monetary variables on log(1 + x) scale
if (nzchar(Sys.getenv("QUICK_TEST"))) {      # fast check: run on an 8,000-row subsample
    d <- d[sample(nrow(d), min(8000, nrow(d))), ]
    message("QUICK_TEST is set: using an 8,000-row subsample (results will differ from the full run).")
}
cat("\nRows used:", nrow(d), "| Card types:", paste(levels(d$Card_Type), collapse = ", "), "\n")
 
# Common 70/30 split and CV folds, reused by every model below
n <- nrow(d)
train_idx <- sample(n, floor(TRAIN_FRAC * n))
train <- d[train_idx, ]; test <- d[-train_idx, ]
folds <- sample(rep(seq_len(K_FOLDS), length.out = nrow(train)))
cat("Train n =", nrow(train), "| Test n =", nrow(test), "| CV folds =", K_FOLDS, "\n")
 
rmse <- function(y, p) sqrt(mean((y - p)^2))
mae  <- function(y, p) mean(abs(y - p))
r2   <- function(y, p) 1 - sum((y - p)^2) / sum((y - mean(y))^2)
metrics <- function(y, p) c(RMSE = rmse(y, p), MAE = mae(y, p), R2 = r2(y, p))
 
 
## B. Section 4.7: linear regression vs trees vs forest vs XGBoost --------------
# Predictors (all outcomes). EXCLUDED because they define or duplicate spending:
#   the seven category spends (they sum exactly to Monthly_Spending),
#   Avg_Transaction_Value (spending / transactions), Reward_Points_Earned
#   (a fixed share of spending by card tier), and Outstanding/Statement balance
#   and Payment_Amount (downstream consequences of spending, r > 0.80).
base_x <- c("Age", "Gender", "Occupation", "Card_Type", "Annual_Income", "Credit_Limit",
            "Card_Age_Months", "Monthly_Transactions", "EMI_Count",
            "International_Transactions", "Mobile_App_Logins", "Payment_Ratio",
            "Credit_Utilization", "Cash_Advance_Amount", "Reward_Points_Redeemed")
specs <- list(
    "Monthly_Spending (log)" = list(y = "Monthly_Spending", x = c(base_x, "Credit_Score")),
    "Credit_Score"           = list(y = "Credit_Score",     x = c(base_x, "Monthly_Spending"))
)
 
# Numeric design matrix for XGBoost (factor levels are fixed, so columns match)
mm <- function(df) {
    X <- model.matrix(~ . - 1, data = df)
    colnames(X) <- make.names(colnames(X))
    X
}
 
# Hyper-parameters tuned once on the training set (OOB error for the forest,
# 5-fold CV with early stopping for XGBoost)
tune_hp <- function(tr, yname, xnames) {
    dd <- tr[, c(yname, xnames)]; frm <- as.formula(paste(yname, "~ ."))
    p <- length(xnames); mtry_grid <- sort(unique(pmin(p, c(3, 5, 8))))
    oob <- sapply(mtry_grid, function(m)
        ranger(frm, dd, num.trees = 200, mtry = m, min.node.size = 5,
               num.threads = NT, seed = SEED, verbose = FALSE)$prediction.error)
    dtr <- xgb.DMatrix(mm(dd[, xnames]), label = dd[[yname]])
    best <- list(rmse = Inf)
    for (md in c(3, 4, 6)) {
        pr <- list(objective = "reg:squarederror", eval_metric = "rmse", eta = 0.1,
                   max_depth = md, min_child_weight = 5, subsample = 0.8,
                   colsample_bytree = 0.8, nthread = NT)
        set.seed(SEED)
        cv <- xgb.cv(params = pr, data = dtr, nrounds = 600, nfold = 5,
                     early_stopping_rounds = 25, verbose = 0)
        lg <- as.data.frame(cv$evaluation_log); i <- which.min(lg$test_rmse_mean)
        if (lg$test_rmse_mean[i] < best$rmse)
            best <- list(rmse = lg$test_rmse_mean[i], params = pr, nrounds = i)
    }
    list(mtry = mtry_grid[which.min(oob)], rf_trees = 300,
         xgb_params = best$params, xgb_nrounds = best$nrounds)
}
 
# Fit all four models and return prediction closures
fit_all <- function(tr, yname, xnames, hp) {
    frm <- as.formula(paste(yname, "~ .")); dd <- tr[, c(yname, xnames)]
    m_lm   <- lm(frm, dd)
    m_tree <- rpart(frm, dd, method = "anova",
                    control = rpart.control(cp = 0.0005, minsplit = 60, xval = 5))
    m_tree <- prune(m_tree, cp = m_tree$cptable[which.min(m_tree$cptable[, "xerror"]), "CP"])
    m_rf   <- ranger(frm, dd, num.trees = hp$rf_trees, mtry = hp$mtry, min.node.size = 5,
                     num.threads = NT, seed = SEED, verbose = FALSE)
    set.seed(SEED)
    m_xgb  <- xgb.train(params = hp$xgb_params, nrounds = hp$xgb_nrounds, verbose = 0,
                        data = xgb.DMatrix(mm(dd[, xnames]), label = dd[[yname]]))
    list(pred = list(
        "Linear regression" = function(nd) as.numeric(predict(m_lm, nd)),
        "Regression tree"   = function(nd) as.numeric(predict(m_tree, nd)),
        "Random forest"     = function(nd) predict(m_rf, nd, num.threads = NT)$predictions,
        "XGBoost"           = function(nd) as.numeric(predict(m_xgb, mm(nd[, xnames])))),
        lm = m_lm, tree = m_tree, rf = m_rf, xgb = m_xgb)
}
 
cv_all <- list(); test_all <- list(); fits <- list(); hps <- list()
for (oc in names(specs)) {
    yname <- specs[[oc]]$y; xnames <- specs[[oc]]$x
    cat("\n=== Section 4.7 outcome:", oc, "===\n")
    hp <- tune_hp(train, yname, xnames); hps[[oc]] <- hp
    cat("Tuned: forest mtry =", hp$mtry, "| XGBoost max_depth =", hp$xgb_params$max_depth,
        "| nrounds =", hp$xgb_nrounds, "\n")
    
    # 5-fold cross-validation on the training set
    cv_rows <- list()
    for (k in seq_len(K_FOLDS)) {
        f <- fit_all(train[folds != k, ], yname, xnames, hp); hold <- train[folds == k, ]
        for (m in names(f$pred))
            cv_rows[[length(cv_rows) + 1]] <- data.frame(
                outcome = oc, fold = k, model = m, t(metrics(hold[[yname]], f$pred[[m]](hold))))
    }
    cv_all[[oc]] <- do.call(rbind, cv_rows)
    
    # Final fit on the full training set, evaluated on the 30% holdout
    fits[[oc]] <- fit_all(train, yname, xnames, hp)
    rows <- lapply(names(fits[[oc]]$pred), function(m) {
        p <- fits[[oc]]$pred[[m]](test); y <- test[[yname]]
        out <- data.frame(outcome = oc, model = m, t(metrics(y, p)))
        out$MAE_original <- if (yname == "Monthly_Spending") mae(expm1(y), expm1(p)) else NA
        out
    })
    test_all[[oc]] <- do.call(rbind, rows)
}
cv_tab <- do.call(rbind, cv_all); test_tab <- do.call(rbind, test_all)
rownames(cv_tab) <- rownames(test_tab) <- NULL
 
cv_sum <- aggregate(cbind(RMSE, MAE, R2) ~ outcome + model, cv_tab,
                    function(v) c(mean = mean(v), sd = sd(v)))
cv_sum <- do.call(data.frame, cv_sum)
lm_rmse <- setNames(test_tab$RMSE[test_tab$model == "Linear regression"],
                    test_tab$outcome[test_tab$model == "Linear regression"])
test_tab$RMSE_vs_linear_pct <- 100 * (test_tab$RMSE / lm_rmse[test_tab$outcome] - 1)
cat("\nCross-validated metrics (mean and SD over folds):\n"); print(cv_sum, digits = 4)
cat("\nHoldout test metrics:\n"); print(test_tab, digits = 4)
#Ty Error Here. 
write.csv(cv_tab,  file.path(OUT_DIR, "tab_47_cv_folds.csv"), row.names = FALSE)
write.csv(cv_sum,  file.path(OUT_DIR, "tab_47_cv_summary.csv"), row.names = FALSE)
write.csv(test_tab, file.path(OUT_DIR, "tab_47_test_metrics.csv"), row.names = FALSE)
 
# Figure: model comparison (CV mean +/- SD bars, holdout as diamonds)
cmp <- merge(cv_sum[, c("outcome", "model", "RMSE.mean", "RMSE.sd")],
             test_tab[, c("outcome", "model", "RMSE")], by = c("outcome", "model"))
cmp$model <- factor(cmp$model, levels = names(model_cols))
fig("fig47_01_model_comparison.png",
    ggplot(cmp, aes(model, RMSE.mean, fill = model)) +
        geom_col(width = 0.65) +
        geom_errorbar(aes(ymin = RMSE.mean - RMSE.sd, ymax = RMSE.mean + RMSE.sd), width = 0.2) +
        geom_point(aes(y = RMSE), shape = 23, size = 3.5, fill = "white", colour = "black") +
        geom_text(aes(label = sprintf("%.3f", RMSE.mean)), vjust = -0.6 - 1.3, size = 3.4) +
        facet_wrap(~ outcome, scales = "free_y") +
        scale_fill_manual(values = model_cols, guide = "none") +
        scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
        labs(x = NULL, y = "RMSE (lower is better)",
             title = "Model comparison: 5-fold CV RMSE (bars, +/- 1 SD) and holdout RMSE (diamonds)") +
        theme(axis.text.x = element_text(angle = 15, hjust = 1)), w = 11, h = 6)
 
# Figure: predicted vs observed for the best holdout model of each outcome
pv <- do.call(rbind, lapply(names(specs), function(oc) {
    tt <- test_tab[test_tab$outcome == oc, ]; best <- tt$model[which.min(tt$RMSE)]
    idx <- sample(nrow(test), min(3000, nrow(test)))
    data.frame(outcome = sprintf("%s: %s (R2 = %.3f)", oc, best, max(tt$R2)),
               obs = test[[specs[[oc]]$y]][idx],
               pred = fits[[oc]]$pred[[best]](test[idx, ]))
}))
fig("fig47_02_pred_vs_observed.png",
    ggplot(pv, aes(pred, obs)) + geom_point(alpha = 0.2, size = 0.7, colour = "#1C7C85") +
        geom_abline(colour = "#C8553D", linewidth = 0.8) + facet_wrap(~ outcome, scales = "free") +
        labs(x = "Predicted", y = "Observed", title = "Best model: predicted vs observed (holdout, 3,000 customers)"),
    w = 11, h = 5.5)
 
# Variable importance: forest permutation, XGBoost gain, tree improvement (all as % of total)
agg_to_var <- function(cols, vars) {
    vapply(cols, function(cn) { hit <- vars[startsWith(cn, make.names(vars))]
    if (length(hit)) hit[which.max(nchar(hit))] else cn }, "")
}
imp_all <- list(); top_vars <- list()
for (oc in names(specs)) {
    yname <- specs[[oc]]$y; xnames <- specs[[oc]]$x
    dd <- train[, c(yname, xnames)]
    rf_imp <- ranger(as.formula(paste(yname, "~ .")), dd, num.trees = 300, mtry = hps[[oc]]$mtry,
                     min.node.size = 5, importance = "permutation", num.threads = NT, seed = SEED, verbose = FALSE)
    ri <- pmax(rf_imp$variable.importance, 0)
    xi <- xgb.importance(model = fits[[oc]]$xgb)
    xg <- tapply(xi$Gain, agg_to_var(xi$Feature, xnames), sum)
    ti <- fits[[oc]]$tree$variable.importance
    one <- function(v, lab) { vals <- as.numeric(v[xnames]); vals[is.na(vals)] <- 0
    data.frame(outcome = oc, variable = xnames, model = lab, pct = 100 * vals / sum(vals)) }
    imp_all[[oc]] <- rbind(one(ri, "Random forest"), one(xg, "XGBoost"), one(ti, "Regression tree"))
    top_vars[[oc]] <- names(sort(ri, decreasing = TRUE))
}
imp <- do.call(rbind, imp_all); rownames(imp) <- NULL
write.csv(imp, file.path(OUT_DIR, "tab_47_variable_importance.csv"), row.names = FALSE)
cat("\nTop-5 variables (random-forest permutation importance):\n")
for (oc in names(specs)) cat(" ", oc, ":", paste(head(top_vars[[oc]], 5), collapse = ", "), "\n")
 
top10 <- do.call(rbind, lapply(names(specs), function(oc) {
    keep <- head(top_vars[[oc]], 8); s <- imp[imp$outcome == oc & imp$variable %in% keep, ]
    s$variable <- factor(s$variable, levels = rev(keep)); s }))
top10$model <- factor(top10$model, levels = c("Regression tree", "Random forest", "XGBoost"))
fig("fig47_03_variable_importance.png",
    ggplot(top10, aes(variable, pct, fill = model)) +
        geom_col(position = position_dodge(width = 0.8), width = 0.75) + coord_flip() +
        facet_wrap(~ outcome, scales = "free") +
        scale_fill_manual(values = model_cols[levels(top10$model)], name = NULL) +
        labs(x = NULL, y = "Share of total importance (%)",
             title = "Variable importance: top 8 predictors by random-forest permutation importance"),
    w = 12, h = 6)
 
# Partial dependence (computed by hand): average prediction as one predictor varies
pdp_curve <- function(pred_fun, data, var, ngrid = 20, nsamp = 1000) {
    bg <- data[sample(nrow(data), min(nsamp, nrow(data))), ]
    grid <- unique(quantile(data[[var]], seq(0.02, 0.98, length.out = ngrid)))
    data.frame(value = grid, pd = sapply(grid, function(g) { bg[[var]] <- g; mean(pred_fun(bg)) }))
}
for (oc in names(specs)) {
    xnames <- specs[[oc]]$x
    nums <- top_vars[[oc]][vapply(top_vars[[oc]], function(v) is.numeric(d[[v]]), logical(1))]
    pv_vars <- head(nums, 4)
    pd <- do.call(rbind, lapply(pv_vars, function(v)
        do.call(rbind, lapply(names(model_cols), function(m) {
            cur <- pdp_curve(fits[[oc]]$pred[[m]], train[, c(specs[[oc]]$y, xnames)], v)
            cbind(variable = v, model = m, cur) }))))
    pd$model <- factor(pd$model, levels = names(model_cols))
    tag <- gsub("[^A-Za-z]+", "_", sub(" \\(log\\)", "", oc))
    fig(sprintf("fig47_04_pdp_%s.png", tag),
        ggplot(pd, aes(value, pd, colour = model, linetype = model)) + geom_line(linewidth = 0.9) +
            facet_wrap(~ variable, scales = "free") +
            scale_colour_manual(values = model_cols, name = NULL) +
            scale_linetype_discrete(name = NULL) +
            labs(x = "Predictor value (monetary variables on log scale)",
                 y = paste("Average predicted", oc),
                 title = paste("Partial dependence:", oc)), w = 11, h = 7)
    write.csv(pd, file.path(OUT_DIR, sprintf("tab_47_pdp_%s.csv", tag)), row.names = FALSE)
}
 
# Pruned regression tree for spending
ms <- "Monthly_Spending (log)"
save_png("fig47_05_regression_tree_spending.png", {
    if (requireNamespace("rpart.plot", quietly = TRUE)) {
        rpart.plot::rpart.plot(fits[[ms]]$tree, type = 4, extra = 101, roundint = FALSE,
                               box.palette = "BuGn", cex = 0.7,
                               main = "Pruned regression tree: log Monthly_Spending")
    } else {                                   # base-R fallback if rpart.plot is missing
        plot(fits[[ms]]$tree, uniform = TRUE, margin = 0.05,
             main = "Pruned regression tree: log Monthly_Spending")
        text(fits[[ms]]$tree, use.n = TRUE, cex = 0.6)
    }
}, w = 2000, h = 1300)
writeLines(capture.output(print(fits[[ms]]$tree)), file.path(OUT_DIR, "tab_47_tree_text.txt"))
cat("Pruned tree leaves (spending):", sum(fits[[ms]]$tree$frame$var == "<leaf>"), "\n")
 
 
## C. Section 4.6: discriminant analysis for Card_Type ---------------------------
# Predictors for classification. Reward_Points_Earned is EXCLUDED: with spending in
# the model its fixed reward rate per tier reveals the label (100% accuracy).
lda_vars <- c("Credit_Limit", "Monthly_Spending", "Monthly_Transactions", "Payment_Ratio",
              "Credit_Utilization", "Credit_Score", "International_Transactions",
              "Mobile_App_Logins", "Cash_Advance_Amount")
cat("\n=== Section 4.6: discriminant analysis ===\n")
 
# C1. Assumptions within card types (Section 3.1 of the proposal)
sub_idx <- sample(nrow(d), min(N_SUB * 2, nrow(d)))
within <- lapply(levels(d$Card_Type), function(ct) {
    X <- as.matrix(d[sub_idx, lda_vars][d$Card_Type[sub_idx] == ct, ]); X <- X[, apply(X, 2, sd) > 0]
    m <- psych::mardia(X, plot = FALSE)
    d2 <- mahalanobis(X, colMeans(X), cov(X))
    list(row = data.frame(card_type = ct, n = nrow(X), mardia_skew = m$b1p, p_skew = m$p.skew,
                          mardia_kurt = m$b2p, p_kurt = m$p.kurt,
                          pct_d2_above_chisq_0.999 = 100 * mean(d2 > qchisq(0.999, ncol(X))),
                          log_det_cov = as.numeric(determinant(cov(X), logarithm = TRUE)$modulus)),
         qq = data.frame(card_type = ct, theo = qchisq(ppoints(length(d2)), ncol(X)), d2 = sort(d2)))
})
norm_tab <- do.call(rbind, lapply(within, `[[`, "row"))
cat("\nWithin-card-type multivariate normality (Mardia) and outlier share:\n")
print(norm_tab, digits = 4)
write.csv(norm_tab, file.path(OUT_DIR, "tab_46_within_cardtype_normality.csv"), row.names = FALSE)
qq <- do.call(rbind, lapply(within, `[[`, "qq")); qq$card_type <- factor(qq$card_type, levels = levels(d$Card_Type))
fig("fig46_01_chisq_qq_by_cardtype.png",
    ggplot(qq, aes(theo, d2)) + geom_point(alpha = 0.3, size = 0.8, colour = "#1C7C85") +
        geom_abline(colour = "#C8553D", linewidth = 0.8) + facet_wrap(~ card_type, scales = "free") +
        labs(x = expression(chi[p]^2 ~ "quantiles"), y = "Squared Mahalanobis distance (within class)",
             title = "Chi-square Q-Q plots within each card type"), w = 11, h = 6.5)
 
# Box's M. H0: equal covariance matrices across card types
box_m <- function(X, group) {
    X <- as.matrix(X); group <- droplevels(factor(group)); p <- ncol(X); g <- nlevels(group)
    n_i <- table(group); N <- sum(n_i)
    S_i <- lapply(levels(group), function(l) cov(X[group == l, , drop = FALSE]))
    S_p <- Reduce(`+`, Map(function(S, k) (k - 1) * S, S_i, n_i)) / (N - g)
    M <- (N - g) * log(det(S_p)) - sum(mapply(function(S, k) (k - 1) * log(det(S)), S_i, n_i))
    c1 <- (sum(1 / (n_i - 1)) - 1 / (N - g)) * (2 * p^2 + 3 * p - 1) / (6 * (p + 1) * (g - 1))
    chi <- M * (1 - c1); df <- p * (p + 1) * (g - 1) / 2
    c(M = M, chisq = chi, df = df, p.value = pchisq(chi, df, lower.tail = FALSE))
}
bm <- box_m(d[sub_idx, lda_vars], d$Card_Type[sub_idx])
cat("\nBox's M (H0: equal covariance matrices across card types):\n"); print(round(bm, 4))
 
# C2. LDA and QDA for the five card types: 5-fold CV and holdout
da <- function(method) if (method == "LDA") lda else qda
cv5 <- do.call(rbind, lapply(c("LDA", "QDA"), function(mth) do.call(rbind, lapply(seq_len(K_FOLDS), function(k) {
    tr <- train[folds != k, c("Card_Type", lda_vars)]; ho <- train[folds == k, c("Card_Type", lda_vars)]
    fit <- da(mth)(Card_Type ~ ., data = tr); pr <- predict(fit, ho)$class
    data.frame(method = mth, fold = k, accuracy = mean(pr == ho$Card_Type),
               balanced_acc = mean(diag(prop.table(table(ho$Card_Type, pr), 1))))
}))))
cat("\n5-fold CV, five card types (mean over folds):\n")
print(aggregate(cbind(accuracy, balanced_acc) ~ method, cv5, mean), digits = 4)
 
lda5 <- lda(Card_Type ~ ., data = train[, c("Card_Type", lda_vars)])
qda5 <- qda(Card_Type ~ ., data = train[, c("Card_Type", lda_vars)])
conf <- list()
cls_tab <- do.call(rbind, lapply(list(LDA = lda5, QDA = qda5), function(fit) NULL))
for (mth in c("LDA", "QDA")) {
    fit <- if (mth == "LDA") lda5 else qda5
    pr <- predict(fit, test)$class; cm <- table(Observed = test$Card_Type, Predicted = pr)
    conf[[mth]] <- cm
    cat("\n", mth, "holdout confusion matrix:\n"); print(cm)
    cat(mth, "holdout accuracy:", round(sum(diag(cm)) / sum(cm), 4),
        "| balanced accuracy:", round(mean(diag(prop.table(cm, 1))), 4), "\n")
}
cat("Majority-class baseline:", round(max(prop.table(table(test$Card_Type))), 4), "\n")
cls_res <- data.frame(
    method = c("LDA", "QDA"),
    cv_accuracy = tapply(cv5$accuracy, cv5$method, mean)[c("LDA", "QDA")],
    test_accuracy = sapply(conf, function(cm) sum(diag(cm)) / sum(cm)),
    test_balanced_acc = sapply(conf, function(cm) mean(diag(prop.table(cm, 1)))),
    baseline = max(prop.table(table(test$Card_Type))))
write.csv(cls_res, file.path(OUT_DIR, "tab_46_lda_qda_summary.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(conf$LDA), file.path(OUT_DIR, "tab_46_lda_confusion.csv"))
write.csv(data.frame(card_type = names(diag(prop.table(conf$LDA, 1))),
                     recall_LDA = diag(prop.table(conf$LDA, 1)), n = rowSums(conf$LDA)),
          file.path(OUT_DIR, "tab_46_lda_recall_by_class.csv"), row.names = FALSE)
 
cm_df <- as.data.frame(prop.table(conf$LDA, 1)); names(cm_df)[3] <- "share"
fig("fig46_02_lda_confusion_heatmap.png",
    ggplot(cm_df, aes(Predicted, Observed, fill = share)) + geom_tile(colour = "white") +
        geom_text(aes(label = sprintf("%.0f%%", 100 * share)), size = 4) +
        scale_fill_gradient(low = "#F4F8F8", high = "#1C7C85", labels = scales::percent, name = "Row share") +
        scale_y_discrete(limits = rev(levels(cm_df$Observed))) +
        labs(title = "LDA holdout confusion matrix (row percentages)"), w = 8, h = 6)
 
ld <- predict(lda5, test[sample(nrow(test), min(3000, nrow(test))), ])
ld_df <- data.frame(ld$x[, 1:2], card = ld$class)
fig("fig46_03_lda_scores.png",
    ggplot(ld_df, aes(LD1, LD2, colour = card)) + geom_point(alpha = 0.45, size = 0.9) +
        labs(title = "LDA discriminant scores (holdout subsample)", colour = NULL), w = 9, h = 6.5)
 
 
## D. Section 4.7: logistic regression for Premium card --------------------------
cat("\n=== Logistic regression: Premium (Platinum/Signature) vs Other ===\n")
prem <- function(df) factor(ifelse(df$Card_Type %in% c("Platinum", "Signature"), "Premium", "Other"),
                            levels = c("Other", "Premium"))
mu <- colMeans(train[lda_vars]); sdv <- apply(train[lda_vars], 2, sd)
prep <- function(df) { x <- as.data.frame(sweep(sweep(as.matrix(df[lda_vars]), 2, mu), 2, sdv, "/"))
x$premium <- prem(df); x }
trb <- prep(train); teb <- prep(test)
cat("Premium share: train", round(mean(trb$premium == "Premium"), 4),
    "| test", round(mean(teb$premium == "Premium"), 4), "\n")
 
auc_fn <- function(score, y01) { r <- rank(score); n1 <- sum(y01 == 1); n0 <- sum(y01 == 0)
(sum(r[y01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
roc_df <- function(score, y01, lab) { o <- order(score, decreasing = TRUE)
data.frame(fpr = c(0, cumsum(1 - y01[o]) / sum(1 - y01)), tpr = c(0, cumsum(y01[o]) / sum(y01)), model = lab) }
cls_metrics <- function(p, y01, thr = 0.5) { pr <- as.integer(p >= thr)
c(accuracy = mean(pr == y01), sensitivity = mean(pr[y01 == 1] == 1), specificity = mean(pr[y01 == 0] == 0),
  AUC = auc_fn(p, y01), Brier = mean((p - y01)^2)) }
 
warns <- character()
fit_glm <- withCallingHandlers(glm(premium ~ ., data = trb, family = binomial()),
                               warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
if (length(warns)) cat("glm warnings:", paste(unique(warns), collapse = " | "), "\n") else cat("glm fit: no warnings\n")
null_glm <- glm(premium ~ 1, data = trb, family = binomial())
lr <- as.numeric(2 * (logLik(fit_glm) - logLik(null_glm))); lr_df <- length(coef(fit_glm)) - 1
cat("Likelihood-ratio test vs null: chi2 =", round(lr, 1), "on", lr_df, "df, p =",
    format.pval(pchisq(lr, lr_df, lower.tail = FALSE), digits = 3, eps = 1e-16), "\n")
cat("McFadden pseudo-R2:", round(1 - as.numeric(logLik(fit_glm) / logLik(null_glm)), 4), "\n")
 
ci <- exp(confint.default(fit_glm)); sm <- summary(fit_glm)$coefficients
or_tab <- data.frame(predictor = rownames(sm), odds_ratio_per_SD = exp(sm[, 1]), ci_low = ci[, 1],
                     ci_high = ci[, 2], z = sm[, 3], p_value = sm[, 4])[-1, ]
or_tab <- or_tab[order(-abs(log(or_tab$odds_ratio_per_SD))), ]; rownames(or_tab) <- NULL
cat("\nOdds ratios per 1-SD increase (95% Wald CI):\n"); print(or_tab, digits = 3)
write.csv(or_tab, file.path(OUT_DIR, "tab_47_logistic_odds_ratios.csv"), row.names = FALSE)
cat("\nVIF (logistic predictors):\n"); print(round(car::vif(fit_glm), 2))
 
or_plot <- or_tab; or_plot$predictor <- factor(or_plot$predictor, levels = rev(or_plot$predictor))
fig("fig47_06_logistic_odds_ratios.png",
    ggplot(or_plot, aes(odds_ratio_per_SD, predictor)) +
        geom_vline(xintercept = 1, colour = "grey60", linetype = 2) +
        geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.25, colour = "#1C7C85") +
        geom_point(size = 2.8, colour = "#C8553D") + scale_x_log10() +
        labs(x = "Odds ratio per 1-SD increase (log scale, 95% CI)", y = NULL,
             title = "Logistic regression: odds of a Premium (Platinum/Signature) card"), w = 9, h = 5.5)
 
# Compare logistic, LDA and QDA (binary) on CV and holdout
bin_fit <- list(
    Logistic = function(tr) { f <- suppressWarnings(glm(premium ~ ., data = tr, family = binomial()))
    function(nd) as.numeric(predict(f, nd, type = "response")) },
    LDA = function(tr) { f <- lda(premium ~ ., data = tr); function(nd) predict(f, nd)$posterior[, "Premium"] },
    QDA = function(tr) { f <- qda(premium ~ ., data = tr); function(nd) predict(f, nd)$posterior[, "Premium"] })
y_tr <- as.integer(trb$premium == "Premium"); y_te <- as.integer(teb$premium == "Premium")
cv_bin <- do.call(rbind, lapply(names(bin_fit), function(m) do.call(rbind, lapply(seq_len(K_FOLDS), function(k) {
    pf <- bin_fit[[m]](trb[folds != k, ]); ho <- trb[folds == k, ]
    data.frame(model = m, fold = k, t(cls_metrics(pf(ho), as.integer(ho$premium == "Premium")))) }))))
cat("\n5-fold CV (mean over folds):\n"); print(aggregate(. ~ model, cv_bin[, -2], mean), digits = 4)
 
roc_list <- list(); test_bin <- list(); probs <- list()
for (m in names(bin_fit)) {
    pf <- bin_fit[[m]](trb); p <- pf(teb); probs[[m]] <- p
    # threshold maximizing Youden's J on the training data
    ptr <- pf(trb); cand <- seq(0.05, 0.95, by = 0.01)
    j <- sapply(cand, function(t) mean(ptr[y_tr == 1] >= t) + mean(ptr[y_tr == 0] < t) - 1)
    thr <- cand[which.max(j)]
    test_bin[[m]] <- rbind(data.frame(model = m, threshold = 0.5, t(cls_metrics(p, y_te, 0.5))),
                           data.frame(model = m, threshold = thr, t(cls_metrics(p, y_te, thr))))
    roc_list[[m]] <- roc_df(p, y_te, m)
}
test_bin <- do.call(rbind, test_bin); rownames(test_bin) <- NULL
cat("\nHoldout test performance (threshold 0.5 and Youden-optimal threshold from training):\n")
print(test_bin, digits = 4)
write.csv(test_bin, file.path(OUT_DIR, "tab_47_logistic_vs_lda_qda_test.csv"), row.names = FALSE)
write.csv(aggregate(. ~ model, cv_bin[, -2], mean), file.path(OUT_DIR, "tab_47_logistic_vs_lda_qda_cv.csv"), row.names = FALSE)
 
cm_log <- table(Observed = teb$premium, Predicted = ifelse(probs$Logistic >= 0.5, "Premium", "Other"))
cat("\nLogistic holdout confusion matrix (threshold 0.5):\n"); print(cm_log)
write.csv(as.data.frame.matrix(cm_log), file.path(OUT_DIR, "tab_47_logistic_confusion.csv"))
 
roc <- do.call(rbind, roc_list); aucs <- test_bin[test_bin$threshold == 0.5, c("model", "AUC")]
roc$model <- factor(roc$model, levels = aucs$model, labels = sprintf("%s (AUC %.3f)", aucs$model, aucs$AUC))
fig("fig47_07_roc_curves.png",
    ggplot(roc, aes(fpr, tpr, colour = model)) + geom_line(linewidth = 1) +
        geom_abline(linetype = 2, colour = "grey60") + coord_equal() +
        scale_colour_manual(values = c("#1C7C85", "#F0A53A", "#C8553D"), name = NULL) +
        labs(x = "False positive rate", y = "True positive rate",
             title = "ROC curves, Premium vs Other card (holdout)"), w = 7.5, h = 6.5)
 
# Calibration of the logistic model (fixed-width probability bins; empty bins dropped)
pp <- probs$Logistic
bin <- cut(pp, seq(0, 1, 0.1), include.lowest = TRUE)
cal <- data.frame(predicted = tapply(pp, bin, mean), observed = tapply(y_te, bin, mean), n = as.integer(table(bin)))
cal <- cal[cal$n > 0, ]
hl <- sum((cal$observed * cal$n - cal$predicted * cal$n)^2 /
              (cal$n * cal$predicted * (1 - cal$predicted)))
cat("\nHosmer-Lemeshow style statistic (", nrow(cal), "bins):", round(hl, 2), " p =",
    format.pval(pchisq(hl, nrow(cal) - 2, lower.tail = FALSE), digits = 3), "\n")
write.csv(cal, file.path(OUT_DIR, "tab_47_logistic_calibration.csv"), row.names = FALSE)
fig("fig47_08_logistic_calibration.png",
    ggplot(cal, aes(predicted, observed)) + geom_abline(linetype = 2, colour = "grey60") +
        geom_line(colour = "#1C7C85") + geom_point(aes(size = n), colour = "#C8553D") +
        labs(x = "Mean predicted probability", y = "Observed Premium rate", size = "Customers",
             title = "Logistic calibration by predicted-probability bin (holdout)"), w = 7.5, h = 6)
 
 
## E. Wrap up -------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo_4_6_4_7.txt"))
cat("\nDone. Outputs written to", normalizePath(OUT_DIR), "\n")
