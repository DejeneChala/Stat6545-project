setwd(r"(C:\Users\tjgfx9\OneDrive - University of Missouri\Masters\Stat6545-project)") #Change as needed. 
credit_data <- read.csv("synthetic_credit_card_customer_behavior_dataset.csv")
attach(credit_data)
library(car)
#Set up dummies for categorical variables
Occn=model.matrix(~Occupation-1)[,-1]
Ctype=model.matrix(~Card_Type-1)[,-1]
GMat=model.matrix(~Gender-1)[,-1]


#Reset the output directory to prevent clogging up the GitHub repository with large output files.
setwd(r"(C:\Users\tjgfx9\OneDrive - University of Missouri\Masters\Stat 6545\Project\Results)")

sink('LinearRegressionResults1.txt',split=TRUE)
X1=cbind(Credit_Score,Occn,Annual_Income,Age)
Y1=cbind(Credit_Limit,Monthly_Spending)
dat=data.frame(Y1,X1)
Y1Full=lm(cbind(Credit_Limit,Monthly_Spending)~.-Credit_Limit-Monthly_Spending,data=dat)
summary(Y1Full)
Z01=dat[1,]
pred1=predict(Y1Full,Z01)
E=t(Y1Full$residuals)%*%Y1Full$residuals
Z=as.matrix(model.matrix(Y1Full))
n=nrow(Y1)
m=ncol(Y1)
r=ncol(Z)-1
Z01=c(1,unlist(Z01[,colnames(Z)[-1]]))
c_prod=data.matrix(solve(crossprod(Z)))
temp=t(Z01)%*%c_prod%*%Z01
temp=c(temp)
alpha=0.05
conf=sqrt(m*(n-r-1)*qf(1-alpha,m,n-r-m)/(n-r-m))
print("Confidence Intervals for the Mean")
print(paste(z0=Z01))
for (i in 1:m){
    Lower=pred1[i]-conf*sqrt(temp*E[i,i]/(n-r-1))
    Upper=pred1[i]+conf*sqrt(temp*E[i,i]/(n-r-1))
    print(paste("Response variable",colnames(pred1)[i]))
    print(c(Lower,Upper))
}
print("Prediction Intervals for a new observation")
for (i in 1:m){
    Lower=pred1[i]-conf*sqrt((1+temp)*E[i,i]/(n-r-1))
    Upper=pred1[i]+conf*sqrt((1+temp)*E[i,i]/(n-r-1))
    print(paste("Response variable",colnames(pred1)[i]))
    print(c(Lower,Upper))
}




#Do a box cox transformation to see what to do? log1p??????
pt1=powerTransform(Y1Full)
pt1
Y1Trans=bcPower(Y1,lambda=pt1$lambda)
colnames(Y1Trans)=colnames(Y1)
datT=data.frame(Y1Trans,X1)
Y1TransFull=lm(cbind(Credit_Limit,Monthly_Spending)~.-Credit_Limit-Monthly_Spending,data=datT)
summary(Y1TransFull)

Z01=datT[1,]
pred1=predict(Y1TransFull,Z01)
E=t(Y1TransFull$residuals)%*%Y1TransFull$residuals
Z=as.matrix(model.matrix(Y1TransFull))
n=nrow(Y1)
m=ncol(Y1)
r=ncol(Z)-1
Z01=c(1,unlist(Z01[,colnames(Z)[-1]]))
c_prod=data.matrix(solve(crossprod(Z)))
temp=t(Z01)%*%c_prod%*%Z01
temp=c(temp)
alpha=0.05
conf=sqrt(m*(n-r-1)*qf(1-alpha,m,n-r-m)/(n-r-m))
print("Confidence Intervals for the Mean")
print(paste(z0=Z01))
for (i in 1:m){
    Lower=pred1[i]-conf*sqrt(temp*E[i,i]/(n-r-1))
    Upper=pred1[i]+conf*sqrt(temp*E[i,i]/(n-r-1))
    print(paste("Response variable",colnames(pred1)[i]))
    print(c(Lower,Upper))
}
print("Prediction Intervals for a new observation")
for (i in 1:m){
    Lower=pred1[i]-conf*sqrt((1+temp)*E[i,i]/(n-r-1))
    Upper=pred1[i]+conf*sqrt((1+temp)*E[i,i]/(n-r-1))
    print(paste("Response variable",colnames(pred1)[i]))
    print(c(Lower,Upper))
}
sink()


#A second Linear Model. Spending by category.
sink('LinearRegressionResults2.txt',split=TRUE)
X2=cbind(Age,Annual_Income,Credit_Limit,Card_Age_Months,GMat,Occn,Ctype)
Y2=cbind(Online_Shopping_Spending,Grocery_Spending,Fuel_Spending,
            Dining_Spending,Travel_Spending,Entertainment_Spending,Utility_Bill_Spending
)
dat2=data.frame(Y2,X2)
Y2Full=lm(cbind(Online_Shopping_Spending,Grocery_Spending,Fuel_Spending,
  Dining_Spending,Travel_Spending,Entertainment_Spending,Utility_Bill_Spending)
  ~.-Online_Shopping_Spending-Grocery_Spending-Fuel_Spending-Dining_Spending
  -Travel_Spending-Entertainment_Spending-Utility_Bill_Spending,data=dat2
)
summary(Y2Full)
Z01=dat2[1,]
pred1=predict(Y2Full,Z01)
E=t(Y2Full$residuals)%*%Y2Full$residuals
Z=as.matrix(model.matrix(Y2Full))
n=nrow(Y2)
m=ncol(Y2)
r=ncol(Z)-1
Z01=c(1,unlist(Z01[,colnames(Z)[-1]]))
c_prod=data.matrix(solve(crossprod(Z)))
temp=t(Z01)%*%c_prod%*%Z01
temp=c(temp)
alpha=0.05
conf=sqrt(m*(n-r-1)*qf(1-alpha,m,n-r-m)/(n-r-m))
print("Confidence Intervals for the Mean")
print(paste(z0=Z01))
for (i in 1:m){
    Lower=pred1[i]-conf*sqrt(temp*E[i,i]/(n-r-1))
    Upper=pred1[i]+conf*sqrt(temp*E[i,i]/(n-r-1))
    print(paste("Response variable",colnames(pred1)[i]))
    print(c(Lower,Upper))
}
print("Prediction Intervals for a new observation")
for (i in 1:m){
    Lower=pred1[i]-conf*sqrt((1+temp)*E[i,i]/(n-r-1))
    Upper=pred1[i]+conf*sqrt((1+temp)*E[i,i]/(n-r-1))
    print(paste("Response variable",colnames(pred1)[i]))
    print(c(Lower,Upper))
}



#Do a box cox transformation to see what to do? log1p??????
pt2=powerTransform(Y2Full)
pt2
Y2Trans=bcPower(Y2,lambda=pt2$lambda)
colnames(Y2Trans)=colnames(Y2)
datT2=data.frame(Y2Trans,X2)
Y2TransFull=lm(cbind(Online_Shopping_Spending,Grocery_Spending,Fuel_Spending,
  Dining_Spending,Travel_Spending,Entertainment_Spending,Utility_Bill_Spending)
  ~.-Online_Shopping_Spending-Grocery_Spending-Fuel_Spending-Dining_Spending
  -Travel_Spending-Entertainment_Spending-Utility_Bill_Spending,data=datT2
)
summary(Y2TransFull)

Z01=datT2[1,]
pred1=predict(Y2TransFull,Z01)
E=t(Y2TransFull$residuals)%*%Y2TransFull$residuals
Z=as.matrix(model.matrix(Y2TransFull))
n=nrow(Y2)
m=ncol(Y2)
r=ncol(Z)-1
Z01=c(1,unlist(Z01[,colnames(Z)[-1]]))
c_prod=data.matrix(solve(crossprod(Z)))
temp=t(Z01)%*%c_prod%*%Z01
temp=c(temp)
alpha=0.05
conf=sqrt(m*(n-r-1)*qf(1-alpha,m,n-r-m)/(n-r-m))
print("Confidence Intervals for the Mean")
print(paste(z0=Z01))
for (i in 1:m){
    Lower=pred1[i]-conf*sqrt(temp*E[i,i]/(n-r-1))
    Upper=pred1[i]+conf*sqrt(temp*E[i,i]/(n-r-1))
    print(paste("Response variable",colnames(pred1)[i]))
    print(c(Lower,Upper))
}
print("Prediction Intervals for a new observation")
for (i in 1:m){
    Lower=pred1[i]-conf*sqrt((1+temp)*E[i,i]/(n-r-1))
    Upper=pred1[i]+conf*sqrt((1+temp)*E[i,i]/(n-r-1))
    print(paste("Response variable",colnames(pred1)[i]))
    print(c(Lower,Upper))
}
sink()




detach(credit_data)
