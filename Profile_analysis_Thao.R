############################################################
# PROFILE ANALYSIS                                    -- Thao

#
# Profile analysis = MANOVA asked in three ordered steps (slides 466-470):
#   Stage 1  PARALLEL?    H01: C mu_1 = C mu_2 = ... = C mu_g
#   Stage 2  COINCIDENT?  H02: 1'mu_1 = 1'mu_2 = ... = 1'mu_g  (only if Stage 1 not rejected)
#   Stage 3  LEVEL?       H03: C mu = 0                         (only if Stage 2 not rejected)
############################################################
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)
while (sink.number() > 0) sink()          # close any sink left open by an earlier error

setwd(r"(C:\Users\songt\Documents\STATProject)")   # Adjust your working directory as needed.
credit_data <- read.csv("synthetic_credit_card_customer_behavior_dataset.csv")

out_dir = "profile_results"
dir.create(out_dir, showWarnings = FALSE)
alpha = 0.05

############################################################
# 0. Helper functions
############################################################
# Tyler's Split.Cov (unchanged)
Split.Cov=function(df,column,vars){
  groups=unique(df[[column]])
  g=length(groups)
  group_data=vector("list",g); group_cov=vector("list",g); ns=vector('list',g)
  for (i in 1:g){
    group_data[[i]]=df %>% filter(.data[[column]]==groups[i]) %>%
      select(all_of(vars)) %>% as.matrix
    group_cov[[i]]=cov(group_data[[i]])
    ns[[i]]=nrow(group_data[[i]])
  }
  p=dim(group_data[[1]])[2]
  Single_Spool=matrix(0,ncol=p,nrow=p); n=0
  for (j in 1:g){
    Single_Spool=Single_Spool+(ns[[j]]-1)*group_cov[[j]]
    n=n+ns[[j]]
  }
  Single_Spool=(1/(n-g))*Single_Spool
  names(group_data)=as.character(groups); names(group_cov)=as.character(groups)
  return(list(X=group_data,Cov=group_cov,Spool=Single_Spool,Ns=ns,total=n))
}

# Box's M on the log scale (slides 415 and 419). Works for n = 50,000.
BoxM.Test=function(Cov,Spool,Ns){
  p=dim(Cov[[1]])[1]; g=length(Cov); n=unlist(Ns)
  logdet=function(A) as.numeric(determinant(A,logarithm=TRUE)$modulus)
  M=sum(n-1)*logdet(Spool)-sum(sapply(1:g,function(i)(n[i]-1)*logdet(Cov[[i]])))
  u=(sum(1/(n-1))-1/sum(n-1))*(2*p^2+3*p-1)/(6*(p+1)*(g-1))
  C=(1-u)*M
  df=p*(p+1)*(g-1)/2
  return(c(C=C,df=df,p.value=pchisq(C,df,lower.tail=FALSE)))
}

############################################################
# 1. Variables and data
############################################################
# All responses must be in SIMILAR UNITS (slide 465). From Tyler's 8
# MANOVA variables we keep the 5 money variables and leave out Age (years),
# Card_Age_Months (months) and Credit_Score (score points).
# Annual_Income is divided by 12 so it is a monthly amount like the others.
# After log1p, the difference between two neighbouring variables is a
# log-ratio (e.g. log Limit - log Income = credit limit relative to income).
ProfVars  = c('Annual_Income','Credit_Limit','Monthly_Spending',
              'Outstanding_Balance','Payment_Amount')
ShortName = c('Income','Limit','Spend','Balance','Payment')
CatVars   = c('Gender','Occupation','Card_Type')

USE_SAMPLE = FALSE      # TRUE = Tyler's 1,000-row random sample (set.seed(42))
dat = credit_data
if (USE_SAMPLE) { set.seed(42); dat = dat[sample(1:nrow(dat),1000),] }

dat$Annual_Income = dat$Annual_Income/12
X  = log1p(as.matrix(dat[,ProfVars]))          # same transform as Tyler
colnames(X) = ShortName
dat$Card_Type = factor(dat$Card_Type, levels=c('Basic','Silver','Gold','Platinum','Signature'))
dat$Gender = factor(dat$Gender); dat$Occupation = factor(dat$Occupation)

p = ncol(X)
N = nrow(X)

############################################################
# 2. Contrast matrix C ((p-1) x p), slide 466
############################################################
C = matrix(0, p-1, p)
for (i in 1:(p-1)) { C[i,i] = -1; C[i,i+1] = 1 }
SegNames = paste0(ShortName[-1],'_minus_',ShortName[-p])   # the p-1 profile segments
dimnames(C) = list(SegNames, ShortName)

CX  = X %*% t(C)                 # Stage 1 responses: segment scores Cx
Tot = rowSums(X)                 # Stage 2 response:  total 1'x
pdat = data.frame(X, CX, Total = Tot, dat[,CatVars])

############################################################
# 3. Box's M on the segment scores Cx (assumption for Stage 1)
############################################################
sink(file.path(out_dir,'1_BoxM_on_contrasts.txt'), split=TRUE)
cat('Box\'s M test on the segment scores Cx, n =', N, '\n')
cat('H0: Sigma_1 = ... = Sigma_g  (equal covariance of Cx across groups)\n\n')
boxm = data.frame()
for (f in CatVars){
  st = Split.Cov(pdat, f, SegNames)
  r  = BoxM.Test(st$Cov, st$Spool, st$Ns)
  cat(f, ': C =', round(r['C'],2), ' df =', r['df'], ' p-value =', signif(r['p.value'],4), '\n')
  boxm = rbind(boxm, data.frame(Factor=f, C=r['C'], df=r['df'], p.value=r['p.value']))
}
cat('\nIf H0 is rejected, base conclusions on Pillai\'s trace (slide 432).\n')
sink()
write.csv(boxm, file.path(out_dir,'1_BoxM_on_contrasts.csv'), row.names=FALSE)

############################################################
# 4. Two-group profile analysis for Gender (exact slide formulas, pp. 467-470)
############################################################
sink(file.path(out_dir,'2_Gender_two_group_profile.txt'), split=TRUE)
st = Split.Cov(as.data.frame(cbind(X, Gender=as.character(dat$Gender))) %>%
                 mutate(across(all_of(ShortName), as.numeric)), 'Gender', ShortName)
n1 = st$Ns[[1]]; n2 = st$Ns[[2]]; Sp = st$Spool
xbar1 = colMeans(st$X[[1]]); xbar2 = colMeans(st$X[[2]])
cat('Groups:', names(st$X)[1], '(n1 =', n1, ') vs', names(st$X)[2], '(n2 =', n2, ')\n\n')

# Stage 1: parallel profiles
CX_mean = C %*% (xbar1 - xbar2)
CX_cov  = C %*% Sp %*% t(C)
T2_1 = as.numeric(t(CX_mean) %*% solve((1/n1+1/n2)*CX_cov) %*% CX_mean)
F_1  = ((n1+n2-p)/((n1+n2-2)*(p-1)))*T2_1          # note the brackets
p_1  = pf(F_1, p-1, n1+n2-p, lower.tail=FALSE)
cat('Stage 1 (parallel):   T2 =', round(T2_1,4), ' F =', round(F_1,4),
    ' df = (', p-1, ',', n1+n2-p, ')  p-value =', signif(p_1,4), '\n')

# Stage 2: coincident profiles
T2_2 = (sum(xbar1 - xbar2)/sqrt((1/n1+1/n2)*sum(Sp)))^2
p_2  = pf(T2_2, 1, n1+n2-2, lower.tail=FALSE)
cat('Stage 2 (coincident): T2 =', round(T2_2,4), ' df = ( 1 ,', n1+n2-2,
    ')  p-value =', signif(p_2,4), '\n')

# Stage 3: level profile (uses S from all n1 + n2 observations)
S  = cov(X); xbar = colMeans(X)
Cx = C %*% xbar
T2_3 = as.numeric((n1+n2) * t(Cx) %*% solve(C %*% S %*% t(C)) %*% Cx)
F_3  = ((n1+n2-p+1)/((n1+n2-1)*(p-1)))*T2_3
p_3  = pf(F_3, p-1, n1+n2-p+1, lower.tail=FALSE)
cat('Stage 3 (level):      T2 =', round(T2_3,4), ' F =', round(F_3,4),
    ' df = (', p-1, ',', n1+n2-p+1, ')  p-value =', signif(p_3,4), '\n\n')

if (p_1 < alpha) {
  cat('Parallel profiles rejected: stop at Stage 1 (Stages 2-3 shown for completeness).\n')
} else if (p_2 < alpha) {
  cat('Profiles parallel but NOT coincident: one group lies above the other.\n')
} else if (p_3 < alpha) {
  cat('Profiles parallel and coincident, but the common profile is NOT level.\n')
} else cat('Profiles parallel, coincident and level.\n')
sink()
write.csv(data.frame(Stage=c('1 Parallel','2 Coincident','3 Level'),
                     T2=c(T2_1,T2_2,T2_3), F=c(F_1,T2_2,F_3),
                     df1=c(p-1,1,p-1), df2=c(n1+n2-p,n1+n2-2,n1+n2-p+1),
                     p.value=c(p_1,p_2,p_3)),
          file.path(out_dir,'2_Gender_two_group_profile.csv'), row.names=FALSE)

# Figure for the report: male - female mean profile with Bonferroni 95% CIs
# (the two profiles overlap in the profile plot, so this shows the differences directly)
d_mf = xbar1 - xbar2                                   # M - F (Split.Cov order)
se_mf = sqrt((1/n1 + 1/n2) * diag(Sp))
t_mf  = qt(1 - alpha/(2*p), n1 + n2 - 2)               # Bonferroni over the p variables
draw_gender_diff = function(){
  par(mfrow = c(1,1), mar = c(3, 4.5, 2.5, 1))
  plot(1:p, d_mf, pch = 19, ylim = range(d_mf - t_mf*se_mf, d_mf + t_mf*se_mf, 0),
       xaxt = 'n', xlab = '', ylab = expression(bar(x)[M] - bar(x)[F]),
       main = 'Male minus female mean profile, Bonferroni 95% intervals')
  arrows(1:p, d_mf - t_mf*se_mf, 1:p, d_mf + t_mf*se_mf, angle = 90, code = 3, length = 0.05)
  abline(h = 0, lty = 2, col = 'grey40')
  axis(1, at = 1:p, labels = ShortName)
}
pdf(file.path(out_dir, 'gender_difference.pdf'), width = 7, height = 4.2)
draw_gender_diff(); dev.off()
png(file.path(out_dir, 'gender_difference.png'), width = 7, height = 4.2, units = 'in', res = 300)
draw_gender_diff(); dev.off()

############################################################
# 5. g-group profile analysis for every factor

#    Stage 1: one-way MANOVA on Cx (4 tests) + ANOVA for each segment
#    Stage 2: one-way ANOVA on the total 1'x
#    Stage 3: level test (same for every factor, computed in Section 4)
############################################################
sink(file.path(out_dir,'3_OneWay_profile.txt'), split=TRUE)
multi_tab = data.frame(); seg_tab = data.frame()
for (f in CatVars){
  cat('==================== Profile analysis for:', f, '====================\n\n')
  fit = manova(as.formula(paste('cbind(',paste(SegNames,collapse=','),') ~',f)), data=pdat)
  cat('Stage 1 - Parallel profiles (MANOVA on Cx)\n')
  for (tt in c('Wilks','Hotelling-Lawley','Roy','Pillai')){
    s = summary(fit, test=tt); print(s)
    multi_tab = rbind(multi_tab, data.frame(Factor=f, Test=tt, Statistic=s$stats[1,2],
                                            approx.F=s$stats[1,3], num.Df=s$stats[1,4], den.Df=s$stats[1,5],
                                            p.value=s$stats[1,6]))
  }
  cat('Follow-up ANOVAs: which segments are not parallel?\n')
  for (v in c(SegNames,'Total')){
    a = summary(aov(as.formula(paste(v,'~',f)), data=pdat))[[1]]
    eta2 = a[1,'Sum Sq']/(a[1,'Sum Sq']+a[2,'Sum Sq'])        # effect size
    stage = ifelse(v=='Total','2 Coincident (1\'x)','1 Parallel (segment)')
    cat(sprintf('  %-22s F = %10.3f  p = %10.4g  partial eta^2 = %.4f\n',
                v, a[1,'F value'], a[1,'Pr(>F)'], eta2))
    seg_tab = rbind(seg_tab, data.frame(Factor=f, Stage=stage, Response=v,
                                        Df=a[1,'Df'], Res.Df=a[2,'Df'], F=a[1,'F value'],
                                        p.value=a[1,'Pr(>F)'], partial.eta2=eta2))
  }
  cat('\n\n')
}
cat('Stage 3 - Level profile (all customers): T2 =', round(T2_3,2),
    ' p-value =', signif(p_3,4), '\n')
sink()
write.csv(multi_tab, file.path(out_dir,'3_OneWay_MANOVA_on_contrasts.csv'), row.names=FALSE)
write.csv(seg_tab,   file.path(out_dir,'3_OneWay_segment_ANOVAs.csv'),     row.names=FALSE)

############################################################
# 6. Bonferroni CIs for the segments and the total (slides 440-442)
#    Version A: m = q g(g-1)/2  (simultaneous over all responses and pairs)
#    Version B: m = g(g-1)/2    (simultaneous only within one response)
#    Only CIs that do not contain 0 are printed; all are saved to CSV.
############################################################
sink(file.path(out_dir,'4_OneWay_Bonferroni_CIs.txt'), split=TRUE)
ci_tab = data.frame()
for (version in c('A: all responses','B: per response')){
  cat('================ Bonferroni CIs, version', version, '================\n')
  cat('Only CIs that dont contain zero are printed\n\n')
  for (f in CatVars){
    for (resp in list(SegNames, 'Total')){
      SplitIt = Split.Cov(pdat, f, resp)
      Spool = as.matrix(SplitIt$Spool); Ns = SplitIt$Ns; n = SplitIt$total
      Xg = SplitIt$X; g = length(Xg); q = length(resp)
      m = if (startsWith(version,'A')) q*g*(g-1)/2 else g*(g-1)/2
      tcrit = qt(1-alpha/(2*m), n-g)
      xmean = lapply(Xg, colMeans)
      for (i in 1:q) for (j in 1:(g-1)) for (k in (j+1):g){
        d  = xmean[[j]][i] - xmean[[k]][i]
        hw = tcrit*sqrt((1/Ns[[j]]+1/Ns[[k]])*Spool[i,i])
        ci_tab = rbind(ci_tab, data.frame(Version=version, Factor=f, Response=resp[i],
                                          Group1=names(Xg)[j], Group2=names(Xg)[k], Difference=d,
                                          Lower=d-hw, Upper=d+hw, Excludes0=(d-hw)*(d+hw)>0))
        if ((d-hw)*(d+hw) > 0){
          cat('Response', resp[i], '| Compare', f, ':', names(Xg)[j], '-', names(Xg)[k], '\n')
          cat('Lower =', d-hw, 'Upper =', d+hw, '\n****************************\n')
        }
      }
    }
  }
  cat('\n')
}
sink()
write.csv(ci_tab, file.path(out_dir,'4_OneWay_Bonferroni_CIs.csv'), row.names=FALSE)

############################################################
# 7. Two-way profile analysis for each pair of factors (as in Tyler's MANOVA)
#    Stage 1: MANOVA on Cx ~ A*B   (A:B row: does the shape effect of A depend on B?)
#    Stage 2: ANOVA on 1'x ~ A*B
#    R uses sequential (Type I) SS, so each model is fitted in BOTH orders;
#    the second-entered factor is the one adjusted for the other.
############################################################
sink(file.path(out_dir,'5_TwoWay_profile.txt'), split=TRUE)
two_tab = data.frame()
for (i in 1:(length(CatVars)-1)) for (j in (i+1):length(CatVars)){
  for (ord in list(c(CatVars[i],CatVars[j]), c(CatVars[j],CatVars[i]))){
    rhs = paste(ord, collapse='*')
    cat('==================== Two-way profile analysis:', rhs, '====================\n')
    fit = manova(as.formula(paste('cbind(',paste(SegNames,collapse=','),') ~',rhs)), data=pdat)
    cat('Stage 1 - Parallel profiles (MANOVA on Cx)\n')
    for (tt in c('Wilks','Pillai')){
      s = summary(fit, test=tt); print(s)
      k = nrow(s$stats)-1
      two_tab = rbind(two_tab, data.frame(Model=rhs, Stage='1 Parallel', Test=tt,
                                          Term=trimws(rownames(s$stats))[1:k], Statistic=s$stats[1:k,2],
                                          approx.F=s$stats[1:k,3], num.Df=s$stats[1:k,4],
                                          den.Df=s$stats[1:k,5], p.value=s$stats[1:k,6]))
    }
    cat('Stage 2 - Coincident profiles (ANOVA on 1\'x)\n')
    a = summary(aov(as.formula(paste('Total ~',rhs)), data=pdat)); print(a)
    a = a[[1]]; k = nrow(a)-1
    two_tab = rbind(two_tab, data.frame(Model=rhs, Stage='2 Coincident', Test='ANOVA',
                                        Term=trimws(rownames(a))[1:k], Statistic=a[1:k,'F value'],
                                        approx.F=a[1:k,'F value'], num.Df=a[1:k,'Df'],
                                        den.Df=a[k+1,'Df'], p.value=a[1:k,'Pr(>F)']))
    cat('\n')
  }
}
sink()
write.csv(two_tab, file.path(out_dir,'5_TwoWay_profile.csv'), row.names=FALSE)

############################################################
# 8. Profile plots  
############################################################
means_tab = data.frame()
for (f in CatVars){
  M = aggregate(X, by=list(Group=dat[[f]]), FUN=mean)
  means_tab = rbind(means_tab, cbind(Factor=f, M, n=as.vector(table(dat[[f]])[as.character(M$Group)])))
}
write.csv(means_tab, file.path(out_dir,'6_group_mean_profiles.csv'), row.names=FALSE)

draw_profiles = function(){
  par(mfrow=c(1,3), mar=c(3,4,3,1))
  for (f in CatVars){
    M = means_tab[means_tab$Factor==f,]
    G = as.character(M$Group); Y = as.matrix(M[,ShortName])
    matplot(1:p, t(Y), type='o', pch=1:length(G), lty=1:length(G), col=1:length(G),
            xaxt='n', xlab='', ylab='Mean log(1 + x)', main=paste('Profiles by',f))
    axis(1, at=1:p, labels=ShortName)
    legend('topright', legend=G, pch=1:length(G), lty=1:length(G),
           col=1:length(G), bty='n', cex=0.8)
  }
  par(mfrow=c(1,1))
}
png(file.path(out_dir,'7_profile_plots.png'), width=15, height=5, units='in', res=300)
draw_profiles(); dev.off()
pdf(file.path(out_dir,'7_profile_plots.pdf'), width=15, height=5)
draw_profiles(); dev.off()
try(draw_profiles())          # on screen; enlarge the Plots pane if it complains

cat('\nAll results saved in:', normalizePath(out_dir), '\n')
print(list.files(out_dir))
