#Loading Packages
library(readxl)
library(dplyr)
library(DescTools)
library(haven)
library(MASS)
library(Mediana)
library(tidyr)

#Loading Excel data
v_cov <- read_excel("Book1 1.xlsx")
View(v_cov)

# Antigen-A/-B neutralinzing antibody titers 
sub<- v_cov%>%pivot_wider(id_cols= SUBJECTID, names_from=PARAMS, values_from= c(BASE,AVAL),values_fill = 0)

# Drop incomplete observations
sub<-sub %>% dplyr::select(SUBJECTID,BASE_A,AVAL_A,BASE_B,AVAL_B) %>% drop_na()

# Compute means
means<-sub %>% summarise(PRE_A=mean(BASE_A,na.rm = T),POST_A=mean(AVAL_A,na.rm = T),PRE_B=mean(BASE_B,na.rm = T),POST_B=mean(AVAL_B,na.rm = T))
means<-as.numeric(means)

# Compute variance-covariance matrix
cc_2<-cov(sub[,2:5],use ="complete.obs")
# Compute correlation matrix
cc<-cor(sub[,2:5],use ="complete.obs")

#Number of simulations
R=10000
#Number of evaluable subjects per arm
N<-542

h<-rep(0,4)
h_all<-0

# Set seed
set.seed(31291591)

#Means for HA (assumed GMT ratio HA/OA=1.25, RSV-A and -B at D31)
meansHA<-means+c(0,log10(1.25),0,log10(1.25))

for (i in 1:R){
  # Individual immunogenicity data are modelled as a multivariate (4-dimensional) normal on the log-10 scale:
  # 1. Anti-RSV-A titer at D1, 2. Anti-RSV-A titer at D31, 3. Anti-RSV-B titer at D1, 4. Anti-RSV-B titer at D31
  # Repeat for 3 arms
  Control<-mvrnorm(N,mu=means,Sigma=cc_2)
  AR<-mvrnorm(N,mu=means,Sigma=cc_2)
  # use means for HA
  HA<-mvrnorm(N,mu=meansHA,Sigma=cc_2)
  
  # Compute fold-increase on log-10 scale (difference) and flag if >= 4
  Control_2<-tibble(PRE_A=Control[,1],POST_A=Control[,2],PRE_B=Control[,3],POST_B=Control[,4])
  Control_2<- Control_2 %>% mutate(FI_A=POST_A-PRE_A, FI_B=POST_B-PRE_B, FI_4A=as.numeric(FI_A>=log10(4)),FI_4B=as.numeric(FI_B>=log10(4))) %>% dplyr::select(POST_A,POST_B,FI_4A,FI_4B)
  
  AR_2<-tibble(PRE_A=AR[,1],POST_A=AR[,2],PRE_B=AR[,3],POST_B=AR[,4])
  AR_2<-AR_2 %>% mutate(FI_A=POST_A-PRE_A, FI_B=POST_B-PRE_B, FI_4A=as.numeric(FI_A>=log10(4)),FI_4B=as.numeric(FI_B>=log10(4))) %>% dplyr::select(POST_A,POST_B,FI_4A,FI_4B)
  
  HA_2<-tibble(PRE_A=HA[,1],POST_A=HA[,2],PRE_B=HA[,3],POST_B=HA[,4])
  HA_2<-HA_2 %>% mutate(FI_A=POST_A-PRE_A, FI_B=POST_B-PRE_B, FI_4A=as.numeric(FI_A>=log10(4)),FI_4B=as.numeric(FI_B>=log10(4))) %>% dplyr::select(POST_A,POST_B,FI_4A,FI_4B)
  
  # T-test for non-inferiority of GMT ratio: AIR vs Control on RSV-A and RSV-B
  GMT_A_AR<-t.test(AR_2$POST_A, Control_2$POST_A,alternative = "greater",mu=-0.176)$p.value
  GMT_B_AR<-t.test(AR_2$POST_B, Control_2$POST_B,alternative = "greater",mu=-0.176)$p.value
  
  # Compute 2-sided 95% and 97.5% CI for SRR difference by MN method: AIR vs Control on RSV-A
  # p-value for SRR difference: if LL at 95% <-0.1 then p=1 (hypothesis can't be rejected both at 1.25% and 2.5% alpha), else if LL at 97.5% <-0.1 then p=0.02 (hypothesis can be rejected at 2.5%, not at 1.25%), else p=0 (hp can be rejected also at 1.25% alpha)
  if(BinomDiffCI(sum(AR_2$FI_4A),N,sum(Control_2$FI_4A),N,conf.level = .95,method = "mn")[2]< -0.1){
    FI4_A_AR<-1 
  } else if (BinomDiffCI(sum(AR_2$FI_4A),N,sum(Control_2$FI_4A),N,conf.level = .975,method = "mn")[2]< -0.1){
    FI4_A_AR<-0.02 
  } else {
    FI4_A_AR<-0
  }
  
  # Compute 2-sided 95% and 97.5% CI for SRR difference by MN method: AIR vs Control on RSV-B
  if(BinomDiffCI(sum(AR_2$FI_4B),N,sum(Control_2$FI_4B),N,conf.level = .95,method = "mn")[2]< -0.1){
    FI4_B_AR<-1 
  } else if (BinomDiffCI(sum(AR_2$FI_4B),N,sum(Control_2$FI_4B),N,conf.level = .975,method = "mn")[2]< -0.1){
    FI4_B_AR<-0.02 
  } else {
    FI4_B_AR<-0
  }
  
  # Repeat for HA vs Control
  GMT_A_HA<-t.test(HA_2$POST_A, Control_2$POST_A,alternative = "greater",mu=-0.176)$p.value
  GMT_B_HA<-t.test(HA_2$POST_B, Control_2$POST_B,alternative = "greater",mu=-0.176)$p.value
  
  if(BinomDiffCI(sum(HA_2$FI_4A),N,sum(Control_2$FI_4A),N,conf.level = .95,method = "mn")[2]< -0.1){
    FI4_A_HA<-1 
  } else if (BinomDiffCI(sum(HA_2$FI_4A),N,sum(Control_2$FI_4A),N,conf.level = .975,method = "mn")[2]< -0.1){
    FI4_A_HA<-0.02 
  } else {
    FI4_A_HA<-0
  }
  
  if(BinomDiffCI(sum(HA_2$FI_4B),N,sum(Control_2$FI_4B),N,conf.level = .95,method = "mn")[2]< -0.1){
    FI4_B_HA<-1 
  } else if (BinomDiffCI(sum(HA_2$FI_4B),N,sum(Control_2$FI_4B),N,conf.level = .975,method = "mn")[2]< -0.1){
    FI4_B_HA<-0.02 
  } else {
    FI4_B_HA<-0
  }
  
  # Overall p-value is max between p of SRR difference and p of GMT ratio (co-primary endpoints)
  AR_A<-max(GMT_A_AR,FI4_A_AR)
  AR_B<-max(GMT_B_AR,FI4_B_AR)
  HA_A<-max(GMT_A_HA,FI4_A_HA)
  HA_B<-max(GMT_B_HA,FI4_B_HA)
  
  # Set graphical testing procedure
  # Raw p-values 
  rawp<-c(HA_A,AR_A,HA_B,AR_B)
  # Initial proportion of alpha allocated
  weight=c(0.75,0.25,0,0)
  
  # Adjusted p-values by graphical procedure
  adjp=AdjustPvalues(rawp,
                     proc = "FallbackAdj",
                     par=parameters(weight=weight))
  
  adjp<-round(adjp, 4)
  # Identify rejected hypotheses
  h<-h+as.numeric(adjp<=0.025)
  h_all<-h_all+all(adjp<=0.025)
}

h<-h/R
h_all<-h_all/R

h
h_all