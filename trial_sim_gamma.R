trial_sim_gamma <- function(n_sim, 
                            n_A, 
                            n_B, 
                            true_delta, 
                            true_sigmasq_A, 
                            true_sigmasq_B, 
                            seed = sample(.Machine$integer.max, 1)) {
  set.seed(seed)
  start_sim = Sys.time();
  all_data_A <- 
    matrix(rgamma(n_sim * n_A, 
                  16/true_sigmasq_A, 
                  4/true_sigmasq_A),
           nrow = n_A) - 4;
  all_data_B <-
    matrix(rgamma(n_sim * n_B, 
                  (1+true_delta)^2/true_sigmasq_B, 
                  (1+true_delta)/true_sigmasq_B),
           nrow = n_B) - 1;
  estimated_sigmasq_A <- apply(all_data_A, 2, var);
  estimated_sigmasq_B <- apply(all_data_B, 2, var);
  z_stat <- 
    (colMeans(all_data_A) - colMeans(all_data_B)) / 
    sqrt(true_sigmasq_A/n_A + true_sigmasq_B/n_B);
  t_stat <-
    z_stat * sqrt(true_sigmasq_A/n_A + true_sigmasq_B/n_B) / 
    sqrt(estimated_sigmasq_A/n_A + estimated_sigmasq_B/n_B);
  #Satterthwaite approximation for degrees of freedom of ttest
  df_approx <- 
    (estimated_sigmasq_A/n_A + estimated_sigmasq_B/n_B)^2 / 
    ((estimated_sigmasq_A/n_A)^2/(n_A-1) + 
       (estimated_sigmasq_B/n_B)^2/(n_B-1));
  list(pvalues_t = 2 * pt(q = abs(t_stat), 
                          df = df_approx, 
                          lower.tail = FALSE), 
       pvalues_z = 2 * pnorm(q = abs(z_stat), 
                             lower.tail = FALSE), 
       runtime = Sys.time() - start_sim);  
}