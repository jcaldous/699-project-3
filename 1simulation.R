##############################
# Title: 699 Project 3 R code
# Date created: 3/6/2023

# Date last edited: 3/9/2023
#     by: Jessica Aldous
#~*~#: Notes to be resolved
##############################

#Effect in control: 0%-15%

#Delta under alternative: 15%

#type 1 error: 0.05
#type 2 error: 0.2


#Hypothesis we are testing:
#H0: Delta=0 
#Ha: Delta =/=0 

######Initial n calculation as a starting point#####
type_1_error <- 0.05; #type 1 error
target_power <- 0.80; #power
assumed_pA = 0.15; #proportion of SOC 
assumed_pB = assumed_pA+0.15; #Porpotion of ND314
assumed_delta <- assumed_pB - assumed_pA; #difference in proportion
assumed_sigmasq_A <- assumed_pA * (1 - assumed_pA); 
assumed_sigmasq_B <- assumed_pB * (1 - assumed_pB); 
ratio_n_B_n_A <- 1; #balanced arms
required_n_A <-
  (qnorm(1 - type_1_error/2) + 
     qnorm(target_power))**2 /
  assumed_delta**2 *
  (assumed_sigmasq_A + assumed_sigmasq_B /
     ratio_n_B_n_A); #general sample size for arm A
required_n_B <- ratio_n_B_n_A * required_n_A; #general sample size for arm B
alpha_int=0.02 #rejection threshold at stage 1 (<0.05)
alpha_fin=(0.05-alpha_int)/(1-alpha_int); #rejection threshold at stage 2 (<0.05)
#~*~#alpha_fin needs some tweaking, working on that rn

#####Number of simulations needed calculations######
#Number of simulations needed calculations
moe=0.01 #margin of error

#conservative: p=0.5
S_con= (qnorm(0.975)/moe)**2 *(.5)*(.5)
S_con
#9604
#Under null: p=0.05
moe=0.001
S_H0= (qnorm(0.975)/moe)**2 *(.05)*(.95)
S_H0 #~2e5
(qnorm(0.975)/(moe/2))**2 *(.05)*(.95)# smaller moe since t1e is so small
#1825; 7298
#Under alternative: p=0.8
S_HA= (qnorm(0.975)/moe)**2 *(.8)*(.2)
S_HA
#6147
#Use 2e5 to be conservative 
nsim=2e5
#####Check initial starting point power and t1e#####
source("trial_sim_2stg.R")
check.pwr<-trial_sim_2stg(nsim,118/2,118/2,assumed_pA,assumed_pB,alpha_int,alpha_fin) #simulate under Ha
mean(check.pwr$decision)#power 0.748 (too small)

check.alp<-trial_sim_2stg(nsim,118/2,118/2,assumed_pA,assumed_pA,alpha_int,alpha_fin) #simulate under H0
mean(check.alp$decision)# 0.047 (a little too small)

####what n do we need to get a power of 0.8?#####
dat<-data.frame(nsim=numeric(0),n_stg=numeric(0),S_t=numeric(0))
for (i in 59:118){
    check.pwr<-trial_sim_2stg(nsim,i,i,assumed_pA,assumed_pB,alpha_int,alpha_fin)
    S_t <-mean(check.pwr$decision)
    if (S_t>=0.8){
      dat_add<-data.frame(nsim=j,n_stg=i,S_t=S_t)
      dat<-rbind(dat,dat_add)
      # print(
      #   list(n_stg=i,
      #        S_t=S_t,
      #        CI=S_t +qnorm(c(0.025,0.975)) * sqrt(S_t * (1-S_t) / nsim)
      #        )
      # )
      break
    }
  }

#####what alpha_fin do we need to get a type 1 of about 0.05?#####
t1e.dat<-data.frame(nsim=numeric(0),n_stg=numeric(0),S_t=numeric(0))

for(j in 1:1e3){
  start_sim = Sys.time()
  sim<-trial_sim_2stg(nsim,66,66,assumed_pA,assumed_pA,alpha_int,alpha_fin) #run sim
  new.alpha.fin<-quantile(sim$pvalues_stg2_quant, probs=c(alpha_fin), na.rm=T) #get new alpha
  new.decision<-(check.alp$pvalues_stg1<alpha_int)*I(check.alp$pvalues_stg1 < alpha_int)+ 
    (check.alp$pvalues_stg2<new.alpha.fin)*(1-I(check.alp$pvalues_stg1 < alpha_int));
      check.t1e<-mean(new.decision, na.rm = T) #check t1e
      dat_add<-data.frame(nsim=j,alpha.fin=new.alpha.fin,S_t=check.t1e)%>%mutate(flag=(S_t>0.045 & S_t<0.05));
      t1e.dat<-rbind(t1e.dat,dat_add)
      # print(
      #   list(n_stg=i,
      #        S_t=S_t,
      #        CI=S_t +qnorm(c(0.025,0.975)) * sqrt(S_t * (1-S_t) / nsim)
      #        )
      # )
      if (j %% 100 == 0){
        print(paste("nsim",j/1e3 *100,"%"))
      }
}
hist(t1e.dat[which(t1e.dat$flag==1),"alpha.fin"])
hist(t1e.dat[which(t1e.dat$flag==1),"S_t"])

summary(t1e.dat[which(t1e.dat$flag==1),"S_t"])
summary(t1e.dat[,"S_t"])#mean t1e is a little high

summary(t1e.dat[which(t1e.dat$flag==1),"alpha.fin"]);#mean alpha.final adjusted
summary(t1e.dat[,"alpha.fin"]);#mean alpha.final adjusted

mean(t1e.dat[which(t1e.dat$flag==1),"S_t"])+qnorm(c(0.025,0.975)) * sqrt(var(t1e.dat[which(t1e.dat$flag==1),"S_t"])) #95% CI of t1e is a little high

#set alpha.fin to mean from above
new.alpha.fin<-round(median(t1e.dat[which(t1e.dat$flag==1),"alpha.fin"]),3); #median alpha.final adjusted: 0.0376
sim<-trial_sim_2stg(nsim,66,66,assumed_pA,assumed_pA,alpha_int,new.alpha.fin) #run sim
check.t1e<-mean(sim$decision, na.rm = T) #check t1e
round(mean(check.t1e)+qnorm(c(0.025,0.5,0.975)) * sqrt(check.t1e*(1-check.t1e)/nsim),4) #95% CI of t1e is a little high
#looks good

#####What n_stage do we need to achieve:t1e~0.05, pwr~0.8, pr_stop>0.40 ?#####
n.new<-data.frame(nsim=numeric(0),n_stg=numeric(0),S_t=numeric(0),t1e=numeric(0),
                  pr_stop=numeric(0),flag=numeric(0))
for (i in 66:70){
    check.pwr<-trial_sim_2stg(nsim,i,i,assumed_pA,assumed_pB,alpha_int,new.alpha.fin)
    check.t1e<-trial_sim_2stg(nsim,i,i,assumed_pA,assumed_pA,alpha_int,new.alpha.fin)
    S_t <-mean(check.pwr$decision)#power     
    t1e <-mean(check.t1e$decision)#type 1 error
    pr_stop<-mean(check.pwr$stopped_int)#power at int
    if (S_t>=0.8 & (t1e<0.050 & t1e>0.048) & (pr_stop>=0.4)){
      dat_add<-data.frame(nsim=j,n_stg=i,S_t=S_t,t1e=t1e,pr_stop=pr_stop)%>%
        mutate(flag=(S_t>=0.8 & (t1e<0.050 & t1e>0.045) & (pr_stop>=0.4)))
      n.new<-rbind(n.new,dat_add)
      # print(
      #   list(n_stg=i,
      #        S_t=S_t,
      #        CI=S_t +qnorm(c(0.025,0.975)) * sqrt(S_t * (1-S_t) / nsim)
      #        )
      # )
       break
     }
  }
(n.new$S_t)
(n.new$t1e);
(n.new$pr_stop);
(n.new$n_stg);#68

n_stg<-(n.new$n_stg)

check.pwr<-trial_sim_2stg(nsim,n_stg,n_stg,assumed_pA,assumed_pB,alpha_int,new.alpha.fin)
S_t <-mean(check.pwr$decision)#power     
pr_st<-mean(check.pwr$stopped_int)#power at int

S_t+qnorm(c(0.5,0.025,0.975)) * sqrt(S_t * (1-S_t) / nsim)
pr_st+qnorm(c(0.5,0.025,0.975)) * sqrt(pr_st * (1-pr_st) / nsim)

#Expected N:
Expect_HA<- (n_stg*2)* pr_st + (n_stg*4)*(1-pr_st)
sd_Ha<-(n_stg*2) * sqrt(pr_st*(1-pr_st)/nsim)
round(Expect_HA+qnorm(c(0.025,0.5,0.975))*sd_Ha,0)

#
check.alp<-trial_sim_2stg(nsim,n_stg,n_stg,assumed_pA,assumed_pA,alpha_int,new.alpha.fin)
t1e <-mean(check.alp$decision, na.rm = T)#t1e
pr_stop<-mean(check.alp$stopped_int, na.rm = T)#t1e at int

t1e+qnorm(c(0.5, 0.025,0.975)) * sqrt(t1e * (1-t1e) / nsim)
pr_stop+qnorm(c(0.5, 0.025,0.975)) * sqrt(pr_stop * (1-pr_stop) / nsim)

#Expected N:
Expect_H0<- (n_stg*2)* pr_stop + (n_stg*4)*(1-pr_stop)
sd_H0<-(n_stg*2) * sqrt(pr_stop*(1-pr_stop)/nsim)

round(Expect_H0+qnorm(c(0.025,0.5,0.975))*sd_H0,0)

#####FINAL VALUES#####
alpha_int #interim alpha threshold: 0.02
new.alpha.fin #final alpha threshold: 0.0375

n_stg #n/4
4*n_stg #max n: 272, blocks of size 2,4,8,34

#randomization: run the following to get a list of 136 randomized assignments exactly 
blockrand(n=136, num.levels = 2, block.sizes = c(1:4), levels=c("SOC","EXP"))
#run it again if advancing to stage 2

#####Different scenarios:Null's#####
#1----
pA=0.10
check.alp<-trial_sim_2stg(nsim,n_stg,n_stg,pA,pA,alpha_int,new.alpha.fin)
t1e<-mean(check.alp$decision)
t1e+qnorm(c(0.5,0.025,0.975)) * sqrt(t1e * (1-t1e) / nsim)

#2----
pA=0.05
check.alp<-trial_sim_2stg(nsim,n_stg,n_stg,pA,pA,alpha_int,new.alpha.fin)
t1e<-mean(check.alp$decision)
t1e+qnorm(c(0.5,0.025,0.975)) * sqrt(t1e * (1-t1e) / nsim)

#####Different scenarios:Non Null's#####
#1----
pA=0.05
pB=0.20
check.pwr<-trial_sim_2stg(nsim,n_stg,n_stg,pA,pB,alpha_int,new.alpha.fin)
power<-mean(check.pwr$decision)

power+qnorm(c(0.5,0.025,0.975)) * sqrt(power * (1-power) / nsim)
#2----
pA=0.10
pB=0.25
check.pwr<-trial_sim_2stg(nsim,n_stg,n_stg,pA,pB,alpha_int,new.alpha.fin)
power<-mean(check.pwr$decision)

power+qnorm(c(0.5,0.025,0.975)) * sqrt(power * (1-power) / nsim)

#####Different scenarios:Additional#####