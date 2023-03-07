##############################
# Title: 699 Project 3 R code
# Date created: 3/6/2023

# Date last edited: 3/6/2023
#     by: Jessica Aldous
##############################

#Effect in control: 0%-15%

#Delta: 15%

#type 1 error: 0.05
#type 2 error: 0.2


#Hypothesis we are testing:
#H0: Delta=0 
#Ha: Delta =/= 0

#general calculations
library(pwr)
power.prop.test(p1=0.15, p2=0.30, 
                sig.level = 0.05,power=0.80, 
                alternative = "two.sided")
#121 per group, 242 overall

source("trial_sim_bernoulli.R")
#####
#Initial n caluclation as a starting point
target_size <- 0.05;
target_power <- 0.80;
assumed_pA = 0.15;
assumed_pB = 0.3;
assumed_delta <- assumed_pB - assumed_pA;
assumed_sigmasq_A <- assumed_pA * (1 - assumed_pA);
assumed_sigmasq_B <- assumed_pB * (1 - assumed_pB);
ratio_n_B_n_A <- 1;
required_n_A <-
  (qnorm(1 - target_size/2) +
     qnorm(target_power))**2 /
  assumed_delta**2 *
  (assumed_sigmasq_A + assumed_sigmasq_B /
     ratio_n_B_n_A);
required_n_B <- ratio_n_B_n_A * required_n_A;
alpha_fin=(0.05-alpha_int)*(1-alpha_int);

check.pwr<-trial_sim_2stg(10000,118/2,118/2,assumed_pA,assumed_pB,0.01)
mean(check.pwr$decision)#power 0.67 (too small)

check.alp<-trial_sim_2stg(10000,118/2,118/2,assumed_pA,assumed_pA,0.01)
mean(check.alp$decision)# 0.03 (too small)
#####
#####
#what do we need to get a power of 0.8?
for (i in 59:200){
  check.pwr<-trial_sim_2stg(1000,i,i,assumed_pA,assumed_pB,0.01)
  S_t <-mean(check.pwr$decision)#power 0.92 (too small)
  if (S_t>=0.8){
    print(
      list(n_stg=i,
           S_t=S_t,
           CI=S_t +qnorm(c(0.025,0.975)) * sqrt(S_t * (1-S_t) / 1000)
           )
    )
    break
  }else{
    print(paste(i,"0.8 not reached:",S_t))
  }
  }
check.pwr<-trial_sim_2stg(10000,118,118,assumed_pA,assumed_pB,0.01)
S_t <-mean(check.pwr$decision)#power 0.92 (too small)
S_t +qnorm(c(0.025,0.975)) * sqrt(S_t * (1-S_t) / 10000) #little underpowered

check.alp<-trial_sim_2stg(10000,118,118,assumed_pA,assumed_pA,0.01)
S_t <-mean(check.alp$decision)#alpha=0.1139
S_t +qnorm(c(0.025,0.975)) * sqrt(S_t * (1-S_t) / 10000) #type 1 error is too low
#####
