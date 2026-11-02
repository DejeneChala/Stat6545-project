############################################################
#MANOVA
############################################################

Split.Cov=function(df,column){
  #Separate the Data Frame into separate groups
  #Calculate covariance matrices for separate groups
  #Calculate Spool for all
  groups=unique(df[[column]])
  g=length(groups)
  group_data=vector("list",g)
  group_cov=vector("list",g)
  ns=vector('list',g)
  for (i in 1:g){
    group_data[[i]]=df %>%                    #Keep only Numeric
      filter(.data[[column]]==groups[i]) %>%
      select(where(is.numeric)) %>%
      as.matrix
    group_cov[[i]]=cov(group_data[[i]])
    ns[[i]]=nrow(group_data[[i]])
  }
  p=dim(group_data[[1]])[2]
  Single_Spool=matrix(0,ncol=p,nrow=p)
  n=0
  for (j in 1:g){
    Single_Spool=Single_Spool+(ns[[j]]-1)*group_cov[[j]]
    n=n+ns[[j]]
  }
  Single_Spool=(1/(n-g))*Single_Spool
  names(group_data)=as.character(groups)
  names(group_cov)=as.character(groups)
  return(list(X=group_data,Cov=group_cov,Spool=Single_Spool,Ns=ns))
}

BoxM.Test=function(Cov,Spool,Ns){
  #Performs Box's M Test to Test H0 = Sig1=Sig2=....
  #Ha = Otherwise.
  #Returns a p-value. 
  p=dim(Cov[[1]])[1]
  g=length(Cov)
  lambda=1
  ls=0
  rsDenom=0
  DetSpool=det(Spool)
  for (i in 1:g){
    holder=(det(Cov[[i]])/DetSpool)**((Ns[[i]]-1)/2)
    lambda=lambda*holder
    ls=ls+1/(Ns[[i]]-1)
    rsDenom=rsDenom+(Ns[[i]]-1)
  }
  M=-2*log(lambda)
  frac=(2*(p**2)+3*p-1)/(6*(p+1)*(g-1))
  rs=1/rsDenom
  u=(ls-rs)*frac
  C=(1-u)*M
  df=p*((p+1)*(g-1)/2)
  pval=pchisq(C,df,lower.tail = F)
  return(pval)
}
#Keep only these variables for MANOVA.
card_kept=
  subset(credit_data,select=c(Age,Gender,Occupation,Card_Type,Annual_Income,
                    Credit_Limit,Monthly_Spending,Outstanding_Balance,
                    Statement_Balance,Payment_Amount,Credit_Score)
)


#Set Up for MANOVA
CatVars=c('Gender','Occupation','Card_Type')
NumVars=c('Age','Annual_Income','Credit_Limit','Monthly_Spending',
          'Outstanding_Balance','Statement_Balance','Payment_Amount',
          'Credit_Score'
)
set.seed(42)
Chosen=sample(1:nrow(card_kept),100)
card_kept_Chosen_Nums=card_kept[Chosen,NumVars]
card_kept_Chosen_Nums=log1p(card_kept_Chosen_Nums)
sub_card_kept=append(card_kept_Chosen_Nums,card_kept[Chosen,CatVars])
sub_card_kept=as.data.frame(card_kept)
log_sub_card_kept_Nums=log1p(sub_card_kept[NumVars])
log_sub_card_kept=append(log_sub_card_kept_Nums,sub_card_kept[CatVars])
log_sub_card_kept=as.data.frame(log_sub_card_kept)
class(log_sub_card_kept)
#Do the Testing.
for (i in 1:length(CatVars)){
  cat("Testing for",CatVars[i],"\n")
  credit_test=Split.Cov(log_sub_card_kept,CatVars[i])
  Cov=credit_test$Cov
  Spool=credit_test$Spool
  Ns=credit_test$Ns
  pval=BoxM.Test(Cov,Spool,Ns)
  print(paste('The p-value for', CatVars[i], 'is:', pval))
  print('**************************************************')
}
