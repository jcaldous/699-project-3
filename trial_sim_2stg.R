trial_sim_2stg <- function(n_sim, #number of simulations
                                n_A, #number in arm 1 stage 1 (ie 1/2 of arm 1)
                                n_B, #number in arm 2 stage 1 (ie 1/2 of arm 2)
                                pA,  #probability in arm 1
                                pB,  #probability in arm 2
                                alpha_int, #interim alpha level
                                alpha_fin, #final alpha level
                                seed = sample(.Machine$integer.max, 1)) {
  set.seed(seed)
  start_sim = Sys.time(); #time it takes to run the simulation
  stg1_data_A <- 
    rbinom(n_sim, n_A, pA) / n_A; #stage 1, arm A
  stg1_data_B <-
    rbinom(n_sim, n_B, pB) / n_B;#stage 1, arm B
  stg1_sigma<- sqrt((stg1_data_A*(1-stg1_data_A)) / (n_A) +
                      (stg1_data_B*(1-stg1_data_B))/(n_B)); #all stage 1
  # estimated_sigmasq_A <- stg1_data_A * (1 -  stg1_data_A);
  # estimated_sigmasq_B <- stg1_data_B * (1 -  stg1_data_B);
  # t_stat <-
  #   (stg1_data_A - stg1_data_B) / 
  #   sqrt(estimated_sigmasq_A / n_A + estimated_sigmasq_B / n_B);
  z_stat <-
    (stg1_data_A - stg1_data_B) / stg1_sigma; #stg 1 test stat (Z)
  # df_approx <- 
  #   (estimated_sigmasq_A/n_A + estimated_sigmasq_B/n_B)^2 / 
  #   ((estimated_sigmasq_A/n_A)^2/(n_A-1) + 
  #      (estimated_sigmasq_B/n_B)^2/(n_B-1));
  # pvalues_stg1_t = 2 * pt(q = abs(t_stat), 
  #                         df = df_approx, 
  #                         lower.tail = FALSE);
  pvalues_stg1_z = 2 * pnorm(q = abs(z_stat), 
                          lower.tail = FALSE)
  
  stg2_data_B <- (stg1_data_B*n_B+rbinom(n_sim, n_B, pB))/(2*n_B);#if there is a stage 2, arm A
  stg2_data_A <- (stg1_data_A*n_A+rbinom(n_sim, n_A, pA))/(2*n_A);#if there is a stage 2, arm A
  stg2_sigma<- sqrt((stg2_data_A*(1-stg2_data_A)) / (2*n_A) +
                      (stg2_data_B*(1-stg2_data_B))/(2*n_B)); #all stage 1
  # estimated_sigmasq_A2 <- stg2_data_A * (1 -  stg2_data_A);
  # estimated_sigmasq_B2 <- stg2_data_B * (1 -  stg2_data_B);
    # t_stat_2 <-
    #   (stg2_data_A - stg2_data_B) / 
    #   sqrt(estimated_sigmasq_A2 / (2*n_A) + estimated_sigmasq_B2 / (2*n_B));
    # df_approx_2 <- 
    #   (estimated_sigmasq_A2/(2*n_A) + estimated_sigmasq_B2/(2*n_B))^2 / 
    #   ((estimated_sigmasq_A2/(2*n_A))^2/((2*n_A)-1) + 
    #      (estimated_sigmasq_B2/(2*n_B))^2/((2*n_B)-1));
  z_stat_2 <-
    (stg2_data_A - stg2_data_B) / stg2_sigma
    # pvalues_stg2_t = 2 * pt(q = abs(t_stat_2), 
    #                         df = df_approx_2, 
    #                         lower.tail = FALSE)
    pvalues_stg2_z = 2 * pnorm(q = abs(z_stat_2), 
                             lower.tail = FALSE)
    pvalues_stg2_quant= ifelse(pvalues_stg1_z < alpha_int, NA,pvalues_stg2_z)
    
    stopped_int=pvalues_stg1_z < alpha_int;
    pvalues_stg1=pvalues_stg1_z;
    pvalues_stg2=pvalues_stg2_z;
    pvalues_stg2_quant=pvalues_stg2_quant;
    decision_1=pvalues_stg1_z<alpha_int;
    decision_2=pvalues_stg2_z<alpha_fin;
    pvalues_z = pvalues_stg2_z*I(pvalues_stg1_z < alpha_int)+ 
      pvalues_stg1_z*(1-I(pvalues_stg1_z < alpha_int));
    dec_z = (z_stat)*stopped_int+ 
      (z_stat_2)*(1-stopped_int);
    decision = (pvalues_stg1<alpha_int)*I(pvalues_stg1 < alpha_int)+ 
      (pvalues_stg2<alpha_fin)*(1-I(pvalues_stg1 < alpha_int));
    
    return(list(alphas=c(alpha_int,alpha_fin),
                stopped_int=stopped_int,
                pvalues_stg1=pvalues_stg1,
                pvalues_stg2=pvalues_stg2,
                pvalues_stg2_quant=pvalues_stg2_quant,
                decision_1=decision_1,
                decision_2=decision_2,
                pvalues_z = pvalues_z,
                dec_z = dec_z,
                decision = decision,
                # decision= (pvalues_stg1_t<alpha_int)*stopped_int+ 
                #   (pvalues_stg2_t<alpha_fin)*(1-stopped_int),
                runtime = Sys.time() - start_sim
        )
    )
  }
#   return(list(stopped_int=pvalues_t < alpha_int,
#               pvalues_t = 2 * pt(q = abs(t_stat), 
#                                  df = df_approx, 
#                                  lower.tail = FALSE),
#               runtime = Sys.time() - start_sim))
# }
