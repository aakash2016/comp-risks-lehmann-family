# ----------------------------------------------------
# Licensed under the MIT License.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %

# Common Functions
G <- function(t, lambda){
  if (member == "weibull") {
    return(1 - exp(-(t**lambda)))
  } else if (member == "gompertz") {
    return(1 - exp(-(exp(t * lambda) - 1)))
  } else if (member == "lomax") {
    return((t * lambda / (1 + t * lambda)))
  } else {
    print("member not found")
  }
}

# dervative of function G
G_dash <- function(t, lambda){
  if (member == "weibull") {
    return(lambda * (t ** (lambda - 1)) * exp(-(t ** lambda)))
  } else if (member == "gompertz") {
    return(lambda * exp(1 + lambda * t - exp(lambda * t)))
  } else if (member == "lomax") {
    return(lambda / (1 + t * lambda) ** 2)
  } else {
    print("member not found")
  }
}

# Survival function
S <- function(t, lambda){
  return (1 - G(t, lambda))
}

# density function
f_pdf <- function(t, lambda) {
  return(G_dash(t, lambda))
}

# Num fails for each cause
set_num_fails <- function(data) {
  m = table(data$cause)
  
  m_0 <<- unname(m[1]) # number of failures due to both causes
  m_1 <<- unname(m[2]) # number of failures due to first cause
  m_2 <<- unname(m[3]) # number of failures due to second cause
  m_c <<- unname(m[4]) # number of censored subjects
}

# same as above
get_num_fails <- function(data) {
  m = table(data$cause)
  
  m_0 <- unname(m[1]) # number of failures due to both causes
  m_1 <- unname(m[2]) # number of failures due to first cause
  m_2 <- unname(m[3]) # number of failures due to second cause
  m_c <- unname(m[4]) # number of censored subjects
  
  return (c(m_0,m_1,m_2,m_c))
}

# Save the results in csv file
write_data <- function(n_iter, 
                       MLE_est,
                       aph_0_est, aph_1_est, aph_2_est, lam_est, 
                       ci_aph0_l, ci_aph0_u, ci_aph1_l, ci_aph1_u,
                       ci_aph2_l, ci_aph2_u, ci_lam_l, ci_lam_u,
                       convergence) {
  
  cols = c("MLE", 
           "alpha_0", "alpha_1", "alpha_2", "lambda",
           "convergence")
  
  df = data.frame(matrix(nrow = n_iter, ncol = length(cols)))
  colnames(df) = cols
  
  df$MLE = MLE_est
  df$alpha_0 = aph_0_est
  df$alpha_1 = aph_1_est
  df$alpha_2 = aph_2_est
  df$lambda = lam_est
  df$CI_alpha_0_l = ci_aph0_l
  df$CI_alpha_0_u = ci_aph0_u
  df$CI_alpha_1_l = ci_aph1_l
  df$CI_alpha_1_u = ci_aph1_u
  df$CI_alpha_2_l = ci_aph2_l
  df$CI_alpha_2_u = ci_aph2_u
  df$CI_lambda_l = ci_lam_l
  df$CI_lambda_u = ci_lam_u
  df$convergence = convergence
  
  return(df)
}

# Save the results of 2 model setup - model selection
write_model_select_2 <- function(n_iter, 
                                 a_mle, b_mle,
                                 aph0_est_a, aph0_est_b,
                                 aph1_est_a, aph1_est_b,
                                 aph2_est_a, aph2_est_b,
                                 lam_est_a, lam_est_b) {
  
  cols = c("mle_a",
           "aph0_a", "aph1_a", "aph2_a", "lam_a", 
           "mle_b",
           "aph0_b", "aph1_b", "aph2_b", "lam_b")
  
  df = data.frame(matrix(nrow = n_iter, ncol = length(cols)))
  colnames(df) = cols
  
  df$mle_a = a_mle
  df$aph0_a = aph0_est_a
  df$aph1_a = aph1_est_a
  df$aph2_a = aph2_est_a
  df$lam_a = lam_est_a
  
  df$mle_b = b_mle
  df$aph0_b = aph0_est_b
  df$aph1_b = aph1_est_b
  df$aph2_b = aph2_est_b
  df$lam_b = lam_est_b
  
  return(df)
}

# Save the results of 3 model setup - model selection
write_model_select_3 <- function(n_iter, 
                                 weibull_mle, gompertz_mle, lomax_mle,
                                 aph0_est_w, aph0_est_g, aph0_est_l,
                                 aph1_est_w, aph1_est_g, aph1_est_l,
                                 aph2_est_w, aph2_est_g, aph2_est_l,
                                 lam_est_w, lam_est_g, lam_est_l) {
  
  cols = c("mle_weibull",
           "aph0_w", "aph1_w", "aph2_w", "lam_w", 
           "mle_gompertz",
           "aph0_g", "aph1_g", "aph2_g", "lam_g", 
           "mle_lomax",
           "aph0_l", "aph1_l", "aph2_l", "lam_l")
  
  df = data.frame(matrix(nrow = n_iter, ncol = length(cols)))
  colnames(df) = cols
  
  
  df$mle_weibull = weibull_mle
  df$aph0_w = aph0_est_w
  df$aph1_w = aph1_est_w
  df$aph2_w = aph2_est_w
  df$lam_w = lam_est_w
  
  df$mle_gompertz = gompertz_mle
  df$aph0_g = aph0_est_g
  df$aph1_g = aph1_est_g
  df$aph2_g = aph2_est_g
  df$lam_g = lam_est_g
  
  df$mle_lomax = lomax_mle
  df$aph0_l = aph0_est_l
  df$aph1_l = aph1_est_l
  df$aph2_l = aph2_est_l
  df$lam_l = lam_est_l
  
  return(df)
}
