############################################################
#MANOVA
############################################################
setwd(r"(C:\Users\tjgfx9\OneDrive - University of Missouri\Masters\Stat6545-project)")
credit_data <- read.csv("synthetic_credit_card_customer_behavior_dataset.csv")
Split.Cov=function(df,column,vars){
  #Separate the Data Frame into separate groups
  #Calculate covariance matrices for separate groups
  #Calculate Spool for all
  groups=unique(df[[column]])
  total=dim(df)[1]
  g=length(groups)
  group_data=vector("list",g)
  group_cov=vector("list",g)
  ns=vector('list',g)
  for (i in 1:g){
    group_data[[i]]=df %>%                    #Keep only Numeric
      filter(.data[[column]]==groups[i]) %>%
      select(all_of(vars)) %>%
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
  return(list(X=group_data,Cov=group_cov,Spool=Single_Spool,Ns=ns,total=n))
}

BoxM.Test=function(Cov,Spool,Ns){
  #Performs Box's M Test to Test H0 = Sig1=Sig2=....
  #Ha = Otherwise.
  #Returns a p-value. 
  #Could have used the biotools package, but I wanted to do it from scratch so that
  #I could also get the covariance matrices and Spool for later use in Bonferroni CI's.
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
Chosen=sample(1:nrow(card_kept),100)    #Start with a random sample of 100 rows
card_kept_chosen=card_kept[Chosen,]     #Here is Chosen Line
#Keep only the numeric and categorical variables in separate matrices.
card_kept_chosen_nums=(card_kept_chosen %>%
  select(all_of(NumVars)) %>%
  as.matrix
)
card_kept_chosen_cats=(card_kept_chosen %>%
  select(all_of(CatVars)) %>%
  as.data.frame
)
#Log transform the numeric variables to reduce skewness. (Dejene confirmed)
card_kept_chosen_nums=log1p(card_kept_chosen_nums)   
card_kept_chosen_nums=as.data.frame(card_kept_chosen_nums)
log_sub_card_kept=cbind(card_kept_chosen_nums,card_kept_chosen_cats)
#Do Box's M Test to see if the covariance matrices are equal across groups for each categorical variable.
for (i in 1:length(CatVars)){
  credit_test=Split.Cov(log_sub_card_kept,CatVars[i],NumVars)
  Cov=credit_test$Cov
  Spool=credit_test$Spool
  Ns=credit_test$Ns
  pval=BoxM.Test(Cov,Spool,Ns)
  cat('Box M Test for:',CatVars[i])
  cat('\n','p-value =',pval)
  cat('\n\n\n\n\n')
}

#One Way MANOVA (and ANOVA)for each categorical variable.
for (i in 1:length(CatVars)){
  cat('MANOVA Results for:',CatVars[i])
  cat('\n\n')
  manova_formula=as.formula(paste('cbind(',paste(NumVars,collapse=','),')~',CatVars[i]))
  manova_model=manova(manova_formula,data=log_sub_card_kept)
  cat('\n')
  print(summary(manova_model,test='Wilks'))
  print(summary(manova_model,test='Hotelling-Lawley'))
  print(summary(manova_model,test='Roy'))
  print(summary(manova_model,test='Pillai'))
  cat('ANOVA Results for:',CatVars[i])
  cat('\n')
  for (j in 1:length(NumVars)){
    anova_formula=as.formula(paste(NumVars[j],'~',CatVars[i]))
    anova_model=aov(anova_formula,data=log_sub_card_kept)
    cat('ANOVA Results for:',NumVars[j])
    cat('\n')
    print(summary(anova_model))
  }
  cat('********************************************************')
  cat('\n\n\n\n\n')
}



#One Way Bonferroni CI's
p=length(NumVars)
alpha=0.05
for (l in 1:length(CatVars)){
  SplitIt=Split.Cov(log_sub_card_kept,CatVars[l],NumVars)
  Spool=SplitIt$Spool
  Ns=SplitIt$Ns
  n=SplitIt$total
  g=length(unique(log_sub_card_kept[[CatVars[l]]]))
  df=(n)-g
  m=p*g*(g-1)/2
  tcrit=qt(1-alpha/(2*m),df)
  X=SplitIt$X
  stopifnot(identical(colnames(X[[1]]),NumVars))
  xmean=list()
  for (i in 1:g){
    xmean[[i]]=colMeans(X[[i]])
  }
  for (i in 1:p){
    for (j in 1:(g-1)){
      for (k in (j+1):g){
        Lower=(xmean[[j]][i]-xmean[[k]][i]-tcrit*sqrt(((1/Ns[[j]])+(1/Ns[[k]]))*Spool[i,i]))
        Upper=(xmean[[j]][i]-xmean[[k]][i]+tcrit*sqrt(((1/Ns[[j]])+(1/Ns[[k]]))*Spool[i,i]))
        if(Lower*Upper>0){
        cat('Response',NumVars[i],'| Compare',CatVars[l],':',names(X)[j],'-',names(X)[k],'\n')
        cat('Lower =',Lower,'Upper =',Upper)
        cat('\n****************************')
        cat('\n\n\n\n\n')
        }
      }
    }
  }
}




#One Way CI's
p=length(NumVars)
alpha=0.05
for (l in 1:length(CatVars)){
  SplitIt=Split.Cov(log_sub_card_kept,CatVars[l],NumVars)
  Spool=SplitIt$Spool
  Ns=SplitIt$Ns
  n=SplitIt$total
  g=length(unique(log_sub_card_kept[[CatVars[l]]]))
  df=(n)-g
  tcrit=qt(1-alpha/2,df)
  X=SplitIt$X
  stopifnot(identical(colnames(X[[1]]),NumVars))
  xmean=list()
  for (i in 1:g){
    xmean[[i]]=colMeans(X[[i]])
  }
  for (i in 1:p){
    for (j in 1:(g-1)){
      for (k in (j+1):g){
        Lower=(xmean[[j]][i]-xmean[[k]][i]-tcrit*sqrt(((1/Ns[[j]])+(1/Ns[[k]]))*Spool[i,i]))
        Upper=(xmean[[j]][i]-xmean[[k]][i]+tcrit*sqrt(((1/Ns[[j]])+(1/Ns[[k]]))*Spool[i,i]))
        if(Lower*Upper>0){
        cat('Response',NumVars[i],'| Compare',CatVars[l],':',names(X)[j],'-',names(X)[k],'\n')
        cat('Lower =',Lower,'Upper =',Upper)
        cat('\n****************************')
        cat('\n\n\n\n\n')
        }
      }
    }
  }
}



#Two Way MANOVA (and ANOVA) for each pair of categorical variables.
for (i in 1:(length(CatVars)-1)){
  for (j in (i+1):length(CatVars)){
    cat('MANOVA Results for:',CatVars[i],'and',CatVars[j])
    cat('\n\n')
    manova_formula=as.formula(paste('cbind(',paste(NumVars,collapse=','),')~',CatVars[i],'*',CatVars[j]))
    manova_model=manova(manova_formula,data=log_sub_card_kept)
    cat('\n')
    print(summary(manova_model,test='Wilks'))
    print(summary(manova_model,test='Hotelling-Lawley'))
    print(summary(manova_model,test='Roy'))
    print(summary(manova_model,test='Pillai'))
    cat('ANOVA Results for:',CatVars[i],'and',CatVars[j])
    cat('\n')
    for (k in 1:length(NumVars)){
      anova_formula=as.formula(paste(NumVars[k],'~',CatVars[i],'*',CatVars[j]))
      anova_model=aov(anova_formula,data=log_sub_card_kept)
      cat('ANOVA Results for:',NumVars[k])
      cat('\n')
      print(summary(anova_model))
    }
    cat('********************************************************')
    cat('\n\n\n\n\n')
  }
}


#In Post, add in Bonferroni CI's and CI's for two way MANOVA.
