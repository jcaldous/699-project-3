trial_sim_2stg <- function(n_sim, #number of simulations
                                n_A, #number in arm 1 stage 1 (ie 1/2 of arm 1)
                                n_B, #number in arm 2 stage 1 (ie 1/2 of arm 2)
                                pA,  #probability in arm 1
                                pB,  #probability in arm 2
                                alpha_int, #interim alpha level
                                alpha_fin,
                                seed = sample(.Machine$integer.max, 1)) {
  set.seed(seed)
  start_sim = Sys.time();
  stg1_data_A <- 
    rbinom(n_sim, n_A, pA) / n_A;
  stg1_data_B <-
    rbinom(n_sim, n_B, pB) / n_B;
  estimated_sigmasq_A <- stg1_data_A * (1 -  stg1_data_A);
  estimated_sigmasq_B <- stg1_data_B * (1 -  stg1_data_B);
  t_stat <-
    (stg1_data_A - stg1_data_B) / 
    sqrt(estimated_sigmasq_A / n_A + estimated_sigmasq_B / n_B);
  df_approx <- 
    (estimated_sigmasq_A/n_A + estimated_sigmasq_B/n_B)^2 / 
    ((estimated_sigmasq_A/n_A)^2/(n_A-1) + 
       (estimated_sigmasq_B/n_B)^2/(n_B-1));
  pvalues_stg1_t = 2 * pt(q = abs(t_stat), 
                          df = df_approx, 
                          lower.tail = FALSE)
  
  stg2_data_B <- (stg1_data_B*n_B+rbinom(n_sim, n_B, pB))/(2*n_B);
  stg2_data_A <- (stg1_data_A*n_A+rbinom(n_sim, n_A, pA))/(2*n_A);
  estimated_sigmasq_A2 <- stg2_data_A * (1 -  stg2_data_A);
  estimated_sigmasq_B2 <- stg2_data_B * (1 -  stg2_data_B);
    t_stat_2 <-
      (stg2_data_A - stg2_data_B) / 
      sqrt(estimated_sigmasq_A2 / (2*n_A) + estimated_sigmasq_B2 / (2*n_B));
    df_approx_2 <- 
      (estimated_sigmasq_A2/(2*n_A) + estimated_sigmasq_B2/(2*n_B))^2 / 
      ((estimated_sigmasq_A2/(2*n_A))^2/((2*n_A)-1) + 
         (estimated_sigmasq_B2/(2*n_B))^2/((2*n_B)-1));
    pvalues_stg2_t = 2 * pt(q = abs(t_stat_2), 
                            df = df_approx_2, 
                            lower.tail = FALSE)
  
    return(list(alphas=c(alpha_int,alpha_fin),
                stopped_int=pvalues_stg1_t < alpha_int,
                pvalues_stg1=pvalues_stg1_t,
                pvalues_stg2=pvalues_stg2_t,
                decision_1=pvalues_stg1_t<alpha_int,
                decision_2=pvalues_stg2_t<alpha_fin,
                pvalues_t = pvalues_stg2_t*I(pvalues_stg1_t < alpha_int)+ 
                  pvalues_stg1_t*(1-I(pvalues_stg1_t < alpha_int)),
                dec_t = t_stat*I(pvalues_stg1_t < alpha_int)+ 
                  t_stat_2*(1-I(pvalues_stg1_t < alpha_int)),
                decision= (pvalues_stg1_t<alpha_int)*stopped_int+ 
                  (pvalues_stg2_t<alpha_fin)*(1-stopped_int),
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
