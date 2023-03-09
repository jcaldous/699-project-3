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
S_H0= (qnorm(0.975)/moe)**2 *(.05)*(.95)
S_H0
#1825
#Under alternative: p=0.8
S_HA= (qnorm(0.975)/moe)**2 *(.8)*(.2)
S_HA
#6147
#Use 10,000 to be conservative 
nsim=10000
#####Check initial starting point power and t1e#####
check.pwr<-trial_sim_2stg(nsim,118/2,118/2,assumed_pA,assumed_pB,alpha_int,alpha_fin) #simulate under Ha
mean(check.pwr$decision)#power 0.67 (too small)

check.alp<-trial_sim_2stg(nsim,118/2,118/2,assumed_pA,assumed_pA,alpha_int,alpha_fin) #simulate under H0
mean(check.alp$decision)# 0.03 (too small)
####what n do we need to get a power of 0.8?#####
dat<-data.frame(nsim=numeric(0),n_stg=numeric(0),S_t=numeric(0))
for(j in 1:100){
  start_sim = Sys.time()
for (i in 59:200){
  check.pwr<-trial_sim_2stg(nsim,i,i,assumed_pA,assumed_pB,alpha_int,alpha_fin)
  S_t <-mean(check.pwr$decision)#power 0.92 (too small)
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
  print(paste("nsim",j,Sys.time() - start_sim))
}
summary(dat$n_stg)#82.8
mean(dat$S_t)+qnorm(c(0.025,0.975)) * sqrt(var(dat$S_t)) #looks okay
#####what alpha_fin do we need to get a type 1 of about 0.05?#####
t1e.dat<-data.frame(nsim=numeric(0),n_stg=numeric(0),S_t=numeric(0))

for(j in 1:1e4){
  start_sim = Sys.time()
  sim<-trial_sim_2stg(nsim,59,59,assumed_pA,assumed_pA,alpha_int,alpha_fin) #run sim
  new.alpha.fin<-quantile(sim$pvalues_stg2_quant, probs=c(alpha_fin), na.rm=T) #get new alpha
  new.decision<-(check.alp$pvalues_stg1<alpha_int)*I(check.alp$pvalues_stg1 < alpha_int)+ 
    (check.alp$pvalues_stg2<new.alpha.fin)*(1-I(check.alp$pvalues_stg1 < alpha_int));
      check.t1e<-mean(new.decision, na.rm = T) #check t1e
      dat_add<-data.frame(nsim=j,alpha.fin=new.alpha.fin,S_t=check.t1e)%>%mutate(flag=(S_t>0.048 & S_t<0.052));
      t1e.dat<-rbind(t1e.dat,dat_add)
      # print(
      #   list(n_stg=i,
      #        S_t=S_t,
      #        CI=S_t +qnorm(c(0.025,0.975)) * sqrt(S_t * (1-S_t) / nsim)
      #        )
      # )
  print(paste("nsim",j,Sys.time() - start_sim))
}
hist(t1e.dat$alpha.fin)
summary(t1e.dat$S_t)
summary(t1e.dat$alpha.fin); #mean alpha.final adjusted
mean(t1e.dat$S_t); #mean t1e is about right
mean(t1e.dat$S_t)+qnorm(c(0.025,0.975)) * sqrt(var(t1e.dat$S_t)) #95% CI of t1e

#set alpha.fin to median from above
new.alpha.fin<-round(median(t1e.dat$alpha.fin),4); #median alpha.final adjusted: 0.0376

check.alp<-trial_sim_2stg(nsim,59,59,assumed_pA,assumed_pA,alpha_int,new.alpha.fin)
S_t <-mean(check.alp$decision)#alpha=0.0499 # okay
S_t +qnorm(c(0.025,0.975)) * sqrt(S_t * (1-S_t) / nsim) #c(0.046, 0.054)

#####What n_stage do we need to achieve:t1e~0.05, pwr~0.8, pr_stop>0.40 ?#####
n.new<-data.frame(nsim=numeric(0),n_stg=numeric(0),S_t=numeric(0),t1e=numeric(0),
                  pr_stop=numeric(0))
for(j in 1:100){
  start_sim = Sys.time()
  for (i in 59:82){
    check.pwr<-trial_sim_2stg(nsim,i,i,assumed_pA,assumed_pB,alpha_int,new.alpha.fin)
    check.t1e<-trial_sim_2stg(nsim,i,i,assumed_pA,assumed_pA,alpha_int,new.alpha.fin)
    S_t <-mean(check.pwr$decision)#power     
    t1e <-mean(check.t1e$decision)#type 1 error
    pr_stop<-mean(check.pwr$stopped_int)#power at int
    if (S_t>=0.8 & (t1e<0.050 & t1e>0.045) & (pr_stop>=0.4)){
      dat_add<-data.frame(nsim=j,n_stg=i,S_t=S_t,t1e=t1e,pr_stop=pr_stop)
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
  print(paste("nsim",j,Sys.time() - start_sim))
}

ggplot(n.new, aes(x=n_stg,y=S_t))+geom_smooth()+geom_point()
summary(n.new$S_t)#range (0.8-0.86)
summary(n.new$t1e); #range (0.044-0.0500)
summary(n.new$pr_stop); #range (0.40-0.43)

summary(n.new$n_stg); #range (63-74)
mean(n.new$S_t)+qnorm(c(0.025,0.975)) * sqrt(var(n.new$S_t)) #95% CI of power 0.79-0.82
mean(n.new$t1e)+qnorm(c(0.025,0.975)) * sqrt(var(n.new$t1e)) #95% CI of t1e 0.046-0.051
mean(n.new$pr_stop)+qnorm(c(0.025,0.975)) * sqrt(var(n.new$pr_stop)) #95% CI of pr_stop 0.39-0.44

#####take the mean (rounded): 68#####
check.alp<-trial_sim_2stg(nsim,68,68,assumed_pA,assumed_pA,alpha_int,new.alpha.fin)
S_t <-mean(check.alp$decision)#alpha=0.0497  #reasonable
S_t +qnorm(c(0.025,0.975)) * sqrt(S_t * (1-S_t) / nsim) #c(0.042, 0.051)(reasonable)

check.pwr<-trial_sim_2stg(nsim,68,68,assumed_pA,assumed_pB,alpha_int,new.alpha.fin)
mean(check.pwr$decision)#reasonable: 0.8348
mean(check.pwr$decision) +qnorm(c(0.025,0.975)) * sqrt(mean(check.pwr$decision) * (1-mean(check.pwr$decision)) / nsim) 
mean(check.pwr$stopped_int)#reasonable: 0.4212

#####FINAL VALUES#####
alpha_int #interim alpha threshold: 0.02
new.alpha.fin #final alpha threshold: 0.0375

n_stg=68 #n/4
4*n_stg #max n: 272, blocks of size 2,4,8,34
#####Different scenarios:Null's#####


#####Different scenarios:Non Null's#####

#####Different scenarios:Additional#####