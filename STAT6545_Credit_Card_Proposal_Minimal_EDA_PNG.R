############################################################
# STAT 6545 Project Proposal
# Credit Card Customer Data
#   1. Correlation matrix
#   2. Histograms of selected informative variables
#
############################################################

# Install once if needed:
# install.packages("corrplot")

library(corrplot)

############################################################
# 1. Read the data
############################################################
# A1
setwd("C:/Users/user/OneDrive - University of Missouri/Fall2026/MVA/Project")
# A2
setwd(r"(C:\Users\tjgfx9\OneDrive - University of Missouri\Masters\Stat6545-project)")
# A3

# A4

############################################################
# Credit Card Customer Dataset
# Enhanced Scatterplot Matrix and Histograms
############################################################
credit_data <- read.csv("synthetic_credit_card_customer_behavior_dataset.csv")
############################################################
# 2. Keep numeric variables and remove ID if it is numeric
############################################################

num_data <- credit_data[, sapply(credit_data, is.numeric), drop = FALSE]

id_names <- c("Customer_ID", "Customer ID", "Customer.ID")
num_data <- num_data[, !names(num_data) %in% id_names, drop = FALSE]

############################################################
# 3. VARIABLES KEPT FOR THE HISTOGRAMS
#
# These are the variables that produced informative/readable
# histograms in the preliminary output.
#
# Variables whose histograms were visually uninformative because
# almost all observations were compressed near zero on the plotted
# scale are intentionally omitted from the proposal-stage figures.
############################################################

hist_vars <- c(  "Age",  "Annual_Income",  "Credit_Limit",  "Card_Age_Months",  "Monthly_Spending",
  "Monthly_Transactions",  "Payment_Ratio",  "Credit_Utilization",  "EMI_Count",  "International_Transactions",
  "Mobile_App_Login",  "Credit_Score")

# Keep only names that actually exist in the data file
hist_vars <- intersect(hist_vars, names(num_data))
hist_data <- num_data[, hist_vars, drop = FALSE]
cat("Variables used for histograms:\n")
print(hist_vars)

############################################################
# 4. VARIABLES EXCLUDED FROM THE HISTOGRAM FIGURE
# These were visually uninformative in the preliminary histogram
# output and are not needed for the proposal-stage EDA figure.
############################################################

excluded_hist_vars <- c("Avg_Transaction_Value",  "Online_Shopping_Spending",  "Grocery_Spending",
  "Fuel_Spending",  "Dining_Spending",  "Travel_Spending",  "Entertainment_Spending",  "Utility_Bill_Spending",
  "Outstanding_Balance",  "Statement_Balance",  "Payment_Amount",  "Cash_Advance_Amount",  "Reward_Points_Earned",
  "Reward_Points_Redeemed")

excluded_hist_vars <- intersect(excluded_hist_vars, names(num_data))
cat("\nVariables excluded from the proposal histogram figure:\n")
print(excluded_hist_vars)

############################################################
# 5. CORRELATION MATRIX
#
############################################################
cor_vars <- c(  "Age",  "Annual_Income",  "Credit_Limit",  "Monthly_Spending",  "Monthly_Transactions",
  "Credit_Utilization",  "Outstanding_Balance",  "Credit_Score")
cor_vars <- intersect(cor_vars, names(num_data))
cor_data <- num_data[, cor_vars, drop = FALSE]
cor_matrix <- cor(  cor_data,  use = "pairwise.complete.obs",  method = "pearson")

cat("\nCorrelation matrix:\n")
print(round(cor_matrix, 3))
# Save the correlation matrix as PNG
png(  filename = "credit_card_correlation_matrix.png",  width = 2000,  height = 1700,  res = 220)

corrplot( cor_matrix,  method = "color",  type = "upper",  order = "original",  addCoef.col = "black",
  number.cex = 0.70,  tl.cex = 0.85,  tl.col = "black",  tl.srt = 45,  diag = TRUE,  mar = c(0, 0, 2, 0),
  title = "Credit Card Customer Correlation Matrix")

dev.off()

############################################################
# 6. COMBINED HISTOGRAM 
############################################################

png( filename = "credit_card_histograms_selected.png",  width = 2600,  height = 1900,  res = 200)
par(  mfrow = c(3, 4),  mar = c(4.0, 4.0, 3.2, 1.0),  oma = c(1, 1, 3, 1),  cex.axis = 0.85,
  cex.lab = 0.95,  cex.main = 0.95)
for (v in hist_vars) {  x <- hist_data[[v]]
  x <- x[is.finite(x)]
  hist( x,breaks = "FD", main = gsub("_", " ", v),xlab = gsub("_", " ", v), ylab = "Frequency", border = "white")}
mtext( "Histograms of Selected Credit Card Customer Variables", outer = TRUE,side = 3,line = 1,cex = 1.25,font = 2)


# the following are new commits, that performs the nultivariate normality tests and the test for the presence of outlier in the data.





# new commit --- Dejene------


############################################################
# 7. MULTIVARIATE NORMALITY AND OUTLIER DETECTION
# Source R codes for STAT 6545 - Chapter 4
############################################################

# Close the histogram PNG opened above
dev.off()
library(rcompanion)
library(lattice)
dir.create("Normality_Outlier_Results", showWarnings = FALSE)
out_dir <- "Normality_Outlier_Results"

############################################################
# 7.1 Select variables for multivariate analysis
############################################################
mvn_vars <- c(
  "Age",
  "Annual_Income",
  "Credit_Limit",
  "Monthly_Spending",
  "Credit_Utilization",
  "Credit_Score")

X0 <- num_data[, mvn_vars, drop = FALSE]
# Remove incomplete observations
valid_rows <- complete.cases(X0) &
  rowSums(!is.finite(as.matrix(X0))) == 0
X <- X0[valid_rows, , drop = FALSE]

# Standardize to avoid numerical scaling problems
Z <- as.data.frame(scale(X));n <- nrow(Z);p <- ncol(Z)
S <- cov(Z)
if (qr(S)$rank < p)
  stop("Selected variables have a singular covariance matrix.")

cat("Sample size:", n, "\n")
cat("Number of variables:", p, "\n")

############################################################
# 7.2 Univariate normal Q-Q plots and rQ
############################################################
png(file.path(out_dir, "Normal_QQ_Plots.png"),
    width = 2400, height = 1600, res = 200)

par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
rQ <- numeric(p)
for (j in seq_len(p)) {
    x <- Z[[j]]
    qq_plot <- qqnorm(x, main = names(Z)[j],    pch = 19, cex = 0.3  )
    qqline(x, col = 2, qtype = 5, lwd = 2)
    rQ[j] <- cor(qq_plot$x, qq_plot$y)
}

dev.off()
rQ_results <- data.frame( Variable = names(Z),  rQ = rQ)
print(rQ_results)
write.csv(rQ_results, file.path(out_dir, "Normal_QQ_Correlations.csv"),
          row.names = FALSE)

############################################################
# 7.3 Eigenvalue and eigenvector analysis
############################################################
Eigen_Analysis <- eigen(cov(Z))
Eigen_Analysis$values
Eigen_Analysis$vectors
# Linear combinations corresponding to the
# largest and smallest eigenvalues
v1 <- as.vector(  as.matrix(Z) %*% Eigen_Analysis$vectors[, 1])

vp <- as.vector(  as.matrix(Z) %*% Eigen_Analysis$vectors[, p])

############################################################
# 7.4 Normality of PCA linear combinations
############################################################
png(file.path(out_dir, "PCA_Normality.png"),
    width = 2000, height = 1700, res = 200)
par(mfrow = c(2, 2))
plotNormalHistogram(v1, main = "Largest Eigenvalue")
qq_plot1 <- qqnorm(v1, main = "Q-Q Plot: v1")
qqline(v1, col = 2, qtype = 5)
plotNormalHistogram(vp, main = "Smallest Eigenvalue")
qq_plotp <- qqnorm(vp, main = "Q-Q Plot: vp")
qqline(vp, col = 2, qtype = 5)
dev.off()
cat("rQ for v1:",    cor(qq_plot1$x, qq_plot1$y), "\n")

cat("rQ for vp:",   cor(qq_plotp$x, qq_plotp$y), "\n")

############################################################
# 7.5 Bivariate normality: Scatterplot matrix
############################################################

set.seed(6545)

# Subsample only for readable scatterplots
plot_rows <- sample(seq_len(n), min(2000, n))
png(file.path(out_dir, "Scatterplot_Matrix.png"),
    width = 2400, height = 2400, res = 200)
pairs(Z[plot_rows, , drop = FALSE],
      pch = 19, cex = 0.3)
dev.off()

############################################################
# 7.6 Bivariate normality: Chi-square plots
############################################################
if (!requireNamespace("heplots", quietly = TRUE))
  install.packages("heplots")
library(heplots)
# All pairs of selected variables
variable_pairs <- combn(names(Z), 2, simplify = FALSE)
for (vars in variable_pairs) {
    filename <- paste0(    "ChiSquare_", vars[1], "_", vars[2], ".png"  )
  
  png(file.path(out_dir, filename),
      width = 1500, height = 1200, res = 200)
  
  cqplot(    Z[, vars, drop = FALSE],    main = paste("Chi-square Plot:",  vars[1], "and", vars[2]))
    dev.off()
}

############################################################
# 7.7 Trivariate normality: Chi-square plot
############################################################
three_vars <- c("Annual_Income",  "Credit_Limit", "Monthly_Spending")

png(file.path(out_dir, "Trivariate_ChiSquare.png"),
    width = 1500, height = 1200, res = 200)

cqplot(  Z[, three_vars, drop = FALSE],  main = "Chi-square Plot: Three Variables")

dev.off()

############################################################
# 7.8 Full multivariate normality: Chi-square plot
############################################################

png(file.path(out_dir, "Multivariate_ChiSquare.png"),
    width = 1500, height = 1200, res = 200)

cqplot(  Z,  main = "Chi-square Plot: All Selected Variables")

dev.off()

############################################################
# 8. MULTIVARIATE OUTLIER DETECTION
############################################################

# Squared Mahalanobis distances
D2 <- mahalanobis(  Z,  center = colMeans(Z),  cov = cov(Z))

# Chi-square cutoff (99% reference region)
alpha <- 0.01
cutoff <- qchisq(1 - alpha, df = p)
outliers <- which(D2 > cutoff)
cat("\nChi-square cutoff:", cutoff, "\n")
cat("Number of potential outliers:",
    length(outliers), "\n")

############################################################
# 8.1 Chi-square Q-Q plot of Mahalanobis distances
############################################################

theoretical <- qchisq(  (seq_len(n) - 0.5) / n,  df = p)

observed <- sort(D2)
png(file.path(out_dir, "Mahalanobis_ChiSquare.png"),
    width = 1600, height = 1300, res = 200)

plot(
  theoretical, observed,
  pch = 19, cex = 0.4,
  xlab = "Theoretical Chi-square Quantiles",
  ylab = "Ordered Squared Mahalanobis Distances",
  main = "Chi-square Q-Q Plot")

abline(0, 1, col = "red", lwd = 2)

dev.off()

############################################################
# 8.2 Mahalanobis distance outlier plot
############################################################

png(file.path(out_dir, "Multivariate_Outliers.png"),
    width = 1800, height = 1200, res = 200)

plot(
  D2,
  pch = 19, cex = 0.35,
  col = ifelse(D2 > cutoff, "red", "#0B1F3A"),
  xlab = "Observation",
  ylab = "Squared Mahalanobis Distance",
  main = "Multivariate Outlier Detection")

abline(h = cutoff, col = "red", lwd = 2, lty = 2)

dev.off()

############################################################
# 8.3 Identify and save potential outliers
############################################################

outlier_results <- data.frame(
  Original_Row = which(valid_rows),
  Mahalanobis_D2 = D2,
  Potential_Outlier = D2 > cutoff)

outlier_results <- outlier_results[
  order(-outlier_results$Mahalanobis_D2),
]

write.csv(outlier_results,
  file.path(out_dir, "Mahalanobis_Outliers.csv"),
  row.names = FALSE)

print(head(outlier_results, 20))

############################################################
# 8.4 Scatterplot matrix for visual outlier detection
# Same method used in Chapter 4 course code
############################################################

splom(Z[plot_rows, , drop = FALSE])

# For interactive point identification in RStudio:
# trellis.focus("panel", 1, 1)
# clicked.points <- panel.link.splom(
#   pch = 19, col = "green"
# )
# clicked.points

############################################################
# FINAL SUMMARY
############################################################

cat("\nMultivariate analysis completed.\n")
cat("Sample size:", n, "\n")
cat("Number of variables:", p, "\n")
cat("Potential outliers:", length(outliers), "\n")
cat("Outlier percentage:",
    round(100 * length(outliers) / n, 2), "%\n")



############################################################
#MANOVA
############################################################

Split.Cov=function(df,column){
  #Separate the Data Frame into separate groups
  #Calculate covariance matrices for separate groups
  groups=unique(df[[column]])
  g=length(groups)
  group_data=vector("list",g)
  group_cov=vector("list",g)
  for (i in 1:g){
    group_data[[i]]=df %>%
      filter(.data[[column]]==groups[i]) %>%
      select(where(is.numeric)) %>%
      as.matrix
    group_cov[[i]]=cov(group_data[[i]])
  }
  names(group_data)=as.character(groups)
  names(group_cov)=as.character(groups)
  return(list(X=group_data,Cov=group_cov))
}






