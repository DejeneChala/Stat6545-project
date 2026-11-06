setwd(r"(C:\Users\tjgfx9\OneDrive - University of Missouri\Masters\Stat6545-project)")
credit_data <- read.csv("synthetic_credit_card_customer_behavior_dataset.csv")
attach(credit_data)
Occn=model.matrix(~Occupation-1)[,-1]
Xs=cbind(Credit_Score,Occn,Annual_Income,Age)
Y=cbind(Credit_Limit,Monthly_Spending)
YFull=lm(Y~Xs)
summary(YFull)
