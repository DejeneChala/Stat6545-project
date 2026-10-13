#Stat 6545 - Project - Exploratory
setwd(r'(C:\Users\tjgfx9\OneDrive - University of Missouri\Masters\Stat6545-project)')
credit_data=read.csv('synthetic_credit_card_customer_behavior_dataset.csv')
credit_data=as.data.frame(credit_data)
ggplot(
  data=credit_data,
  aes(x=Credit_Limit,y=Statement_Balance,colour=Occupation)
)+geom_point()

ggplot(
  data=credit_data,
  aes(x=Credit_Limit,y=Monthly_Spending,colour=Card_Type)
)+geom_point()


ggplot(
  data=credit_data,
  aes(x=Annual_Income/12,y=Monthly_Spending,colour=Occupation))+
  geom_point()+
  geom_abline(slope=1,intercept=0,colour='black')


ggplot(
  data=credit_data,
  aes(x=Annual_Income/12,y=Monthly_Spending,colour=Gender))+
  geom_point()+
  geom_abline(slope=1,intercept=0,colour='black')

ggplot(
  data=credit_data,
  aes(x=Annual_Income/12,y=Monthly_Spending,colour=Card_Type))+
  geom_point()+
  geom_abline(slope=1,intercept=0,colour='black')

P=ggplot(
  data=credit_data
)
P+geom_point(aes(x=Annual_Income,y=Monthly_Spending,colour=Card_Type))

P+geom_histogram(aes(x=Monthly_Spending,colour=Card_Type))

P+geom_point(aes(x=Annual_Income/12,y=Monthly_Spending,colour=Credit_Score))

P+geom_point(aes(x=Credit_Score,y=Entertainment_Spending))

P+geom_point(aes(x=Credit_Score,y=EMI_Count,alpha=0.0001,colour=NA))


