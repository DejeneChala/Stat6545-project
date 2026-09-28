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

dev.off()



