trial_sim_bernoulli <- function(n_sim, 
                                n_A, 
                                n_B, 
                                pA,
                                pB,
                                seed = sample(.Machine$integer.max, 1)) {
  set.seed(seed)
  start_sim = Sys.time();
  all_data_A <- 
    rbinom(n_sim, n_A, pA) / n_A;
  all_data_B <-
    rbinom(n_sim, n_B, pB) / n_B;
  estimated_sigmasq_A <- all_data_A * (1 -  all_data_A);
  estimated_sigmasq_B <- all_data_B * (1 -  all_data_B);
  t_stat <-
    (all_data_A - all_data_B) / 
    sqrt(estimated_sigmasq_A / n_A + estimated_sigmasq_B / n_B);
  df_approx <- 
    (estimated_sigmasq_A/n_A + estimated_sigmasq_B/n_B)^2 / 
    ((estimated_sigmasq_A/n_A)^2/(n_A-1) + 
       (estimated_sigmasq_B/n_B)^2/(n_B-1));
  list(pvalues_t = 2 * pt(q = abs(t_stat), 
                          df = df_approx, 
                          lower.tail = FALSE), 
       runtime = Sys.time() - start_sim);  
}