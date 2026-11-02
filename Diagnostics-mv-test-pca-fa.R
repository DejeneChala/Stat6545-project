
# STAT 6545: PCA, factor analysis and normality screening (requires psych)

set.seed(123)
setwd("C:/Users/ddcdq7/OneDrive - University of Missouri/Fall2026/MVA/Project")
DATA_FILE <- "synthetic_credit_card_customer_behavior_dataset.csv"
OUT <- "PCA_FA_results"
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
stopifnot(file.exists(DATA_FILE), requireNamespace("psych", quietly=TRUE))
raw <- read.csv(DATA_FILE, check.names=TRUE)
if ("Average_Transaction_Value" %in% names(raw) &&
    !"Avg_Transaction_Value" %in% names(raw))
  raw$Avg_Transaction_Value <- raw$Average_Transaction_Value

vars <- c("Annual_Income","Credit_Limit","Monthly_Spending",
          "Avg_Transaction_Value","Online_Shopping_Spending","Grocery_Spending",
          "Fuel_Spending","Dining_Spending","Travel_Spending",
          "Entertainment_Spending","Utility_Bill_Spending",
          "Outstanding_Balance","Payment_Amount","Payment_Ratio",
          "Credit_Utilization","Cash_Advance_Amount","Credit_Score")
monetary <- setdiff(vars, c("Payment_Ratio", "Credit_Utilization", "Credit_Score"))

missing <- setdiff(vars, names(raw))
if(length(missing)) stop("Missing expected PCA/FA columns: ",paste(missing,collapse=", "))
x <- raw[,vars,drop=FALSE]
if(!all(vapply(x,is.numeric,logical(1)))) stop("PCA/FA columns must be numeric")
keep_rows <- apply(x,1,function(z) all(is.finite(z)))
x <- x[keep_rows,,drop=FALSE]
row_id <- which(keep_rows)
if(nrow(x)<100) stop("Not enough complete finite observations.")
if(any(as.matrix(x[,monetary,drop=FALSE]) < 0))
  stop("Negative monetary values: log1p is not applicable; inspect data.")
# Diagnostic decision: skewed monetary variables -> log1p; remaining variables unchanged.
prepared <- x
prepared[,monetary] <- lapply(x[,monetary,drop=FALSE],log1p)
# Audit constant and exactly/nearly linearly redundant predictors.
positive_sd <- vapply(prepared,function(z) is.finite(sd(z)) && sd(z)>0,logical(1))
z0 <- scale(prepared[,positive_sd,drop=FALSE])
q <- qr(z0,tol=1e-8)
kept <- colnames(z0)[sort(q$pivot[seq_len(q$rank)])]
write.csv(data.frame(Variable=vars,Log1p=vars %in% monetary,
                     Retained=vars %in% kept),file.path(OUT,"Preparation_Audit.csv"),row.names=FALSE)
if(length(kept)<3) stop("Too few independent columns.")
x <- x[,kept,drop=FALSE]
prepared <- prepared[,kept,drop=FALSE]
Z <- scale(prepared)
n <- nrow(Z); p <- ncol(Z)
write.csv(data.frame(Original_Row=row_id),file.path(OUT,"Included_Rows.csv"),row.names=FALSE)
write.csv(data.frame(Variable=kept,Original_Mean=colMeans(x),
                     Prepared_Mean=colMeans(prepared),Prepared_SD=apply(prepared,2,sd)),
          file.path(OUT,"Scaling_Parameters.csv"),row.names=FALSE)

# Output helpers: PDF and PNG use identical plotting functions.
plot_both <- function(name,fun,w=9,h=7) {
  pdf(file.path(OUT,paste0(name,".pdf")),width=w,height=h,useDingbats=FALSE)
  tryCatch(fun(),finally=dev.off())
  png(file.path(OUT,paste0(name,".png")),width=w,height=h,units="in",res=300)
  tryCatch(fun(),finally=dev.off())
}
heat <- function(M,main,lim=1) {
  M <- as.matrix(M); nr <- nrow(M); nc <- ncol(M)
  pal <- colorRampPalette(c("#2166AC","white","#B2182B"))(101)
  par(mar=c(9,16,3,2))
  plot(0,type="n",xlim=c(.5,nc+.5),ylim=c(.5,nr+.5),
       axes=FALSE,xlab="",ylab="",main=main)
  for(i in seq_len(nr)) for(j in seq_len(nc)) {
    v <- M[i,j]; yy <- nr-i+1
    ix <- if(is.finite(v)) 1+round(100*(pmin(lim,pmax(-lim,v))+lim)/(2*lim)) else NA
    rect(j-.5,yy-.5,j+.5,yy+.5,col=if(is.na(ix)) "grey90" else pal[ix],border="white")
  }
  axis(1,seq_len(nc),gsub("_"," ",colnames(M)),las=2,cex.axis=.75)
  axis(2,seq_len(nr),rev(gsub("_"," ",rownames(M))),las=1,cex.axis=.75)
}
write_table <- function(d,name) write.csv(d,file.path(OUT,paste0(name,".csv")),row.names=FALSE)
# PCA: comparison uses exactly the same observations and transformed variables.
# Correlation PCA is PRIMARY because original variables had different scales.
pca_corr <- prcomp(prepared,center=TRUE,scale.=TRUE)
pca_cov <- prcomp(prepared,center=TRUE,scale.=FALSE)
pca_raw_cov <- prcomp(x,center=TRUE,scale.=FALSE) # original-scale sensitivity
summary_pca <- function(f) {
  eig <- f$sdev^2
  data.frame(Component=paste0("PC",seq_along(eig)),Eigenvalue=eig,
             Variance_Percent=100*eig/sum(eig),
             Cumulative_Percent=100*cumsum(eig)/sum(eig))
}
vc <- summary_pca(pca_corr)
write_table(vc,"Variance_Explained")
write_table(summary_pca(pca_cov),"Covariance_PCA_Transformed")
write_table(summary_pca(pca_raw_cov),"Covariance_PCA_Original")
k <- which(vc$Cumulative_Percent>=80)[1]
V <- pca_corr$rotation
cor_pc <- sweep(V,2,pca_corr$sdev,"*")
write.csv(V,file.path(OUT,"PCA_Loadings.csv"))
write.csv(cor_pc,file.path(OUT,"Variable_PC_Correlations.csv"))
write.csv(100*V^2,file.path(OUT,"Variable_Contributions_Percent.csv"))
write_table(data.frame(Original_Row=row_id,pca_corr$x),"PCA_Scores")
write_table(data.frame(N=n,Variables=p,Retained_80pct=k,
                       Retained_Kaiser=sum(vc$Eigenvalue>1),
                       Cumulative_Variance_Percent=vc$Cumulative_Percent[k]),"PCA_Retention")
saveRDS(pca_corr,file.path(OUT,"PCA_Model.rds"))
plot_both("Scree_Plot",function(){
  par(mar=c(5,5,3,2))
  plot(seq_len(p),vc$Eigenvalue,type="b",pch=19,xlab="Principal component",
       ylab="Eigenvalue",main="Correlation-based PCA")
  abline(h=1,lty=2,col="red")
})
plot_both("Cumulative_Variance",function(){
  par(mar=c(5,5,3,2))
  plot(seq_len(p),vc$Cumulative_Percent,type="b",pch=19,ylim=c(0,100),
       xlab="Number of components",ylab="Cumulative variance (%)")
  abline(h=80,lty=2,col="red")
  abline(v=k,lty=3)})
plot_both("PCA_Loadings_Heatmap",function(){
  heat(cor_pc[,seq_len(k),drop=FALSE],"Variable-PC correlations")
},w=9,h=8)
plot_both("PCA_PC1_PC2_Contributions",function(){
  par(mfrow=c(1,2),mar=c(5,13,3,1))
  for(j in 1:2){
    v <- sort(100*V[,j]^2)
    barplot(v,horiz=TRUE,names.arg=gsub("_"," ",names(v)),las=1,
            cex.names=.7,col="#2166AC",xlab="Contribution (%)",
            main=paste("PC",j))
    abline(v=100/p,lty=2,col="red")
  }
},w=14,h=8)
plot_both("PCA_PC1_PC2",function(){
  set.seed(123); ii <- sample.int(n,min(n,3000))
  plot(pca_corr$x[ii,1:2],pch=19,cex=.35,
       xlab="PC1",ylab="PC2",main="Correlation PCA (sampled scores)")
})
plot_both("PCA_Covariance_Comparison",function(){
  par(mar=c(5,5,3,2))
  a <- summary_pca(pca_cov)$Cumulative_Percent
  b <- summary_pca(pca_raw_cov)$Cumulative_Percent
  plot(seq_len(p),vc$Cumulative_Percent,type="l",lwd=2,ylim=c(0,100),
       xlab="Number of components",ylab="Cumulative variance (%)")
  lines(seq_len(p),a,lty=2,lwd=2)
  lines(seq_len(p),b,lty=3,lwd=2)
  abline(h=80,col="red",lty=3)
  legend("bottomright",legend=c("Prepared correlation","Prepared covariance",
                                "Original covariance"),lty=1:3,lwd=2,bty="n",cex=.8)})

# EFA: ML inference is approximate because the diagnostics rejected normality.
# Assess factorability, candidate fits, and common-factor residuals.
R <- cor(Z)
KMO <- psych::KMO(R)
bart <- psych::cortest.bartlett(R,n=n)
write_table(data.frame(KMO=KMO$MSA,Bartlett_Chi_Square=bart$chisq,
                       df=bart$df,p_value=bart$p.value),"Factorability")
write_table(data.frame(Variable=names(KMO$MSAi),MSA=as.numeric(KMO$MSAi)),
            "Variable_MSA")
candidates <- seq_len(min(6L,p-1L))
candidates <- candidates[((p-candidates)^2-p-candidates)/2>0]
cmp <- list()
for(m in candidates){
  f <- tryCatch(suppressWarnings(psych::fa(R,nfactors=m,n.obs=n,fm="ml",rotate="none")),error=function(e) NULL)
  if(is.null(f)) next
  ok <- all(is.finite(f$uniquenesses)) && all(f$uniquenesses>0) &&
    all(f$uniquenesses<=1)
  cmp[[length(cmp)+1]] <- data.frame(Factors=m,df=f$dof,
                                     Chi_square=f$STATISTIC,p_value=f$PVAL,BIC=f$BIC,Admissible=ok)
}
if(!length(cmp)) stop("All ML candidate FA models failed.")
cmp <- do.call(rbind,cmp)
write_table(cmp,"Factor_Count_Comparison")
eligible <- which(cmp$Admissible & is.finite(cmp$BIC))
if(!length(eligible)) stop("No admissible EFA fit; see factor comparison.")
best <- cmp$Factors[eligible[which.min(cmp$BIC[eligible])]]
fit <- psych::fa(R,nfactors=best,n.obs=n,fm="ml",rotate="varimax")
L <- unclass(fit$loadings)
colnames(L) <- paste0("ML",seq_len(ncol(L)))
write.csv(L,file.path(OUT,"Varimax_Factor_Loadings.csv"))
write_table(data.frame(Variable=rownames(L),Communality=fit$communality,
                       Uniqueness=fit$uniquenesses),"Communalities_Uniquenesses")
reproduced <- L%*%t(L)+diag(fit$uniquenesses)
res <- R-reproduced
write.csv(res,file.path(OUT,"Residual_Correlations.csv"))
off <- res[lower.tri(res)]
metrics <- data.frame(Factors=best,Chi_square=fit$STATISTIC,df=fit$dof,
                      p_value=fit$PVAL,BIC=fit$BIC,RMSR=sqrt(mean(off^2)),
                      TLI=if(length(fit$TLI)) fit$TLI else NA_real_,
                      RMSEA=if(length(fit$RMSEA)) fit$RMSEA[1] else NA_real_,
                      RMSEA_Lower=if(length(fit$RMSEA)>1) fit$RMSEA[2] else NA_real_,
                      RMSEA_Upper=if(length(fit$RMSEA)>2) fit$RMSEA[3] else NA_real_,
                      Max_Absolute_Residual=max(abs(off)),
                      Residual_Fraction_Above_005=mean(abs(off)>.05),
                      Near_Boundary_Uniquenesses=sum(fit$uniquenesses<.05))
write_table(metrics,"FA_Fit_Statistics")
if(!is.null(fit$Vaccounted)) write.csv(fit$Vaccounted,
                                       file.path(OUT,"Factor_Variance_Accounting.csv"))
saveRDS(fit,file.path(OUT,"Factor_Model.rds"))
plot_both("Factor_Count_BIC",function(){
  plot(cmp$Factors,cmp$BIC,type="b",pch=19,xlab="Number of factors",
       ylab="BIC (smaller is better)")
  abline(v=best,col="red",lty=2)})
plot_both("FA_Loadings_Heatmap",function(){
  heat(L,"Varimax factor loadings")
},w=10,h=8)
plot_both("FA_Communalities_Uniquenesses",function(){
  par(mar=c(5,15,3,2))
  H <- rbind(Communality=fit$communality,Uniqueness=fit$uniquenesses)
  barplot(H,horiz=TRUE,names.arg=gsub("_"," ",rownames(L)),las=1,
          cex.names=.75,col=c("#2166AC","#CCCCCC"),xlim=c(0,1),
          xlab="Proportion of variable variance")
  legend("bottomright",legend=rownames(H),fill=c("#2166AC","#CCCCCC"),
         cex=.8,bty="n")
},w=10,h=8)
plot_both("FA_Residual_Correlation_Heatmap",function(){
  S <- res; diag(S) <- NA
  heat(S,"Observed - reproduced correlations",lim=max(.05,max(abs(off))))
},w=11,h=10)
capture.output(list(PCA_Retention=read.csv(file.path(OUT,"PCA_Retention.csv")),
                    EFA_Fit=metrics,KMO=KMO$MSA,
                    Normality_Note="Prior six-variable diagnostics: marginal tests reject normality;",
                    Caution="ML EFA chi-square and RMSEA assume approximate multivariate normality."),
               file=file.path(OUT,"Analysis_Summary.txt"))
cat("Finished. Outputs:",normalizePath(OUT),"\n")

DATA_FILE <- "C:/Users/ddcdq7/OneDrive - University of Missouri/Fall2026/MVA/Project/synthetic_credit_card_customer_behavior_dataset.csv"
out <- file.path(dirname(DATA_FILE), "Normality_Cleaning")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(out, "figures")
res_dir <- file.path(out, "results")
dir.create(fig_dir, showWarnings = FALSE)
dir.create(res_dir, showWarnings = FALSE)
EXTRA_DIAGNOSTICS <- TRUE # Produce both original diagnostic PDFs AND report figures
set.seed(123)
alpha <- .05
B <- 1999L
vars <- c("Age", "Annual_Income", "Credit_Limit", "Monthly_Spending",
          "Credit_Utilization", "Credit_Score")
money <- c("Annual_Income", "Credit_Limit", "Monthly_Spending")
save_csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE)
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
money <- intersect(money, names(X))
if (any(as.matrix(X[money]) < 0)) stop("Review negative monetary values before log1p.")
prepared <- X
prepared[money] <- lapply(prepared[money], log1p)
save_csv(data.frame(Original_Row = rows, prepared), "Prepared_Selected_Data.csv")

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
  pc <- prcomp(Z)$x[, c(1L, p), drop = FALSE]
  # PC directions are fitted from these data: rQ is descriptive, not calibrated.
  save_csv(data.frame(PC = colnames(pc), rQ = apply(pc, 2, function(y)
    cor(sort(y), qnorm(ppoints(length(y), a = .5))))), paste0(label, "_PC_rQ.csv"))
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
  pdf(file.path(out, paste0(label, "_Diagnostics.pdf")), width = 11, height = 8.5)
  on.exit(dev.off())
  par(mfrow = c(2, 4), mar = c(4, 4, 3, 1), cex.main = .8)
  for (v in names(x)) {
    hist(x[[v]], breaks = "FD", main = v, xlab = "", col = "grey80")
    qqnorm(x[[v]], main = v, pch = 16, cex = .2); qqline(x[[v]], col = 2)
  }
  par(mfrow = c(2, 2))
  for (j in 1:2) {
    hist(pc[,j], breaks = "FD", main = colnames(pc)[j], xlab = "Score")
    qqnorm(pc[,j], main = colnames(pc)[j], pch = 16, cex = .2)
    qqline(pc[,j], col = 2)
  }
  par(mfrow = c(1, 1))
  pairs(Z[plot_rows, ], pch = 16, cex = .25,
        col = ifelse(flags$Flag_99pct[plot_rows], "red", "grey50"))
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1), cex.main = .7)
  for (v in sets) {
    A <- Z[, v, drop = FALSE]
    ds <- mahalanobis(A, colMeans(A), cov(A))
    plot(qchisq(ppoints(n, a = .5), length(v)), sort(ds), pch = 16, cex = .2,
      main = paste(v, collapse = " / "), xlab = "Chi-square quantiles", ylab = "Ordered D-squared")
    abline(0, 1, col = 2)
  }
  }
  invisible(list(tests = tests, pc = pc, flags = flags, iqr = iqr_flag, radial = radial))}
original <- screen(X, "Original")
transformed <- screen(prepared, "Log1p")
raw_flags <- original$flags
prepared_flags <- transformed$flags
cat("Results saved in:", out, "\n")

save_figure <- function(name, width, height, draw) {
  for (ext in c("pdf", "png")) {
    path <- file.path(fig_dir, paste0(name, ".", ext))
    if (ext == "pdf") pdf(path, width = width, height = height,
                          useDingbats = FALSE)
    else png(path, width = width, height = height, units = "in", res = 300)
    tryCatch(draw(), finally = dev.off())  }}
blue <- "#195D9C"
red <- "#B83B36"
save_figure("normal_qq", 13, 9, function() {
  par(mfrow = c(ceiling(ncol(X)/2), 4), mar = c(4, 4, 3, 1),
      mgp = c(2.4, .7, 0), cex.main = .9, cex.axis = .8, cex.lab = .85)
  q_all <- qnorm(ppoints(nrow(X), a = .5))
  for (v in names(X)) for (scale_name in c("Original", "Prepared")) {
    y <- if (scale_name == "Original") X[[v]] else prepared[[v]]
    plot(q_all, sort(y), pch = 16, cex = .22, col = blue,
         main = paste(gsub("_", " ", v), scale_name, sep = "\n"),
         xlab = "Normal quantile", ylab = "Sample quantile")
    qqline(y, col = red, lwd = 1.3) }})
save_figure("mahalanobis_qq", 10, 4.2, function() {
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1), mgp = c(2.8, .8, 0))
  theoretical <- qchisq(ppoints(nrow(X), a = .5), ncol(X))
  for (label in c("Original", "Prepared")) {
    distances <- if (label == "Original") raw_flags$D2 else prepared_flags$D2
    plot(theoretical, sort(distances), pch = 16, cex = .25, col = blue,
         main = paste(label, "scale"),
         xlab = paste0("Chi-square quantile (df = ", ncol(X), ")"),
         ylab = "Ordered squared Mahalanobis distance")
    abline(0, 1, col = red, lwd = 1.3)  }})

report_table <- function(tab, name, caption, label, digits = 4) {
  write.csv(tab, file.path(res_dir, paste0(name, ".csv")), row.names = FALSE)
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
  writeLines(lines, file.path(res_dir, paste0(name, ".tex")))}
a <- original$tests; b <- transformed$tests
normality <- data.frame(Variable = a$Variable,`Original rQ` = a$rQ, `Prepared rQ` = b$rQ,
  `Original SW adj. p` = a$SW_Bonferroni, `Prepared SW adj. p` = b$SW_Bonferroni, check.names = FALSE)
report_table(normality, "Normality_Comparison", "Normality diagnostics before and after preparation.",
             "tab:normality")
# Full test files in the root also contain Monte Carlo rQ p-values.
normality_decisions <- data.frame(Variable = a$Variable,
  Original_rQ_reject = a$rQ_Bonferroni < alpha,
  Prepared_rQ_reject = b$rQ_Bonferroni < alpha,
  Original_SW_reject = a$SW_Bonferroni < alpha,
  Prepared_SW_reject = b$SW_Bonferroni < alpha)
write.csv(normality_decisions, file.path(res_dir, "Normality_Decisions.csv"), row.names = FALSE)
pc_rq <- function(pc) apply(pc, 2, function(y)
  cor(sort(y), qnorm(ppoints(length(y), a = .5))))
joint <- data.frame(Diagnostic = c("First PC rQ", "Last PC rQ", "Central 50 percent: full vector"),
  Original = c(pc_rq(original$pc), tail(original$radial$Proportion_inside_50pct, 1)),
  Prepared = c(pc_rq(transformed$pc), tail(transformed$radial$Proportion_inside_50pct, 1)))
if (ncol(X) > 3) joint <- rbind(joint, data.frame(
  Diagnostic = "Central 50 percent: first three variables",
  Original = original$radial$Proportion_inside_50pct[nrow(original$radial)-1],
  Prepared = transformed$radial$Proportion_inside_50pct[nrow(transformed$radial)-1]))
report_table(joint, "Multivariate_Comparison",
             "Joint diagnostics; central-region entries are proportions.", "tab:joint")
ni <- c(colSums(original$iqr), sum(rowSums(original$iqr) > 0))
nt <- c(colSums(transformed$iqr), sum(rowSums(transformed$iqr) > 0))
iqr_summary <- data.frame(Variable = c(names(X), "At least one variable"),
  `Original n` = ni, `Original percent` = 100*ni/nrow(X),
  `Prepared n` = nt, `Prepared percent` = 100*nt/nrow(X), check.names = FALSE)
report_table(iqr_summary, "IQR_Comparison", "Univariate outlier screening; counts overlap across variables.","tab:iqr", digits = 2)
mo <- colSums(raw_flags[c("Flag_99pct", "Flag_Bonferroni")])
mt <- colSums(prepared_flags[c("Flag_99pct", "Flag_Bonferroni")])
mahal <- data.frame(Rule = c("99 percent", "Bonferroni"),
  Cutoff = qchisq(c(.99, 1-alpha/nrow(X)), ncol(X)),
  `Original n` = mo, `Original percent` = 100*mo/nrow(X),
  `Prepared n` = mt, `Prepared percent` = 100*mt/nrow(X), check.names = FALSE)
report_table(mahal, "Mahalanobis_Comparison", "Approximate multivariate outlier screening.",
             "tab:mahalanobis", digits = 3)
audit_summary <- data.frame(Item = c("Input rows", "Retained rows", "Excluded rows",
    "Missing cells", "Nonfinite nonmissing cells", "Rows in duplicate records",
    "Rows with duplicate IDs", "Selected constant columns"),
  Count = c(nrow(d), sum(ok), sum(!ok), sum(audit$Missing),
    sum(audit$Nonfinite_nonmissing), sum(dup), sum(dup_id), sum(!keep)))
write.csv(audit_summary, file.path(res_dir, "Audit_Summary.csv"), row.names = FALSE)
write.csv(prepared_flags[order(-prepared_flags$D2), ],
          file.path(res_dir, "Prepared_Outliers_Ranked.csv"), row.names = FALSE)
writeLines(c("All complete observations are used for plots and outlier screening.",
  paste("Normality test sample size:", length(ii)), paste("Normal simulations:", B),
  paste("Minimum Monte Carlo p-value:", 1/(B+1)),
  "PC rQ is descriptive. Chi-square thresholds are approximate screens.",
  "No outliers or duplicate rows were automatically deleted.",
  "Only selected monetary variables were transformed; see money in the script.",
  "Repeat checks for each analysis variable set and within relevant groups."),
  file.path(res_dir, "Analysis_Notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(res_dir, "R_Session_Info.txt"))
cat("\nReport figures:", fig_dir, "\nReport tables:",res_dir, "\n")
