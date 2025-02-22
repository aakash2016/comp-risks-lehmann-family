# ----------------------------------------------------
# Licensed under the MIT License.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Simulation for general lehmann family using Bootstrapping ##
# Percentile Bootstrapping for CIs #

rm(list=ls())
options(warn=-1)

source("scripts/utils.R")
library(VGAM)
library(foreach)
library(doParallel)
n_cores = detectCores()
cls <- makeCluster(n_cores)
registerDoParallel(cls)

member <<- "lomax" # specify the member of the Family
n_dpts = 400 # sample size
censoring = TRUE # whether to have censoring or not
tau = 1.00 # censoring factor

##
B = 250
n_iter = 500 # Number of monte carlo replications
true_param = c(0.85, 0.57, 0.74, 0.69)

# Estimates
aph_0_est = aph_1_est = aph_2_est = lam_est = MLE_est = rep(0,n_iter)
B_estimates_aph0 = B_estimates_aph1 = B_estimates_aph2 = B_estimates_lam = rep(0,B)

# CIs
ci_aph0_l = ci_aph0_u = ci_aph1_l = ci_aph1_u = rep(0,n_iter)
ci_aph2_l = ci_aph2_u = ci_lam_l = ci_lam_u = rep(0,n_iter)

# Average Length OF CIs.
AL_aph0 = AL_aph1 = AL_aph2 = AL_lam = 0

# Coverage Percentage of CIs.
cp_aph0 = cp_aph1 = cp_aph2 = cp_lam = 0

# Convergence of MLEs
convergence = rep(0,n_iter)

# Log Likelihood function
LL <- function(x,lifetime_df){
  sum_ = 0
  sum_j = 0
  
  for(i in 1:nrow(lifetime_df)) {
    row <- lifetime_df[i, ]
    t = row$lifetime
    sum_j = sum_j + log(S(t, x))
    
    if (row$cause != 3) {
      sum_ = sum_ + log(f_pdf(t, x) / S(t, x))
    }
  }
  
  SS <<- sum_j
  sum_ = sum_ - (m_0 + m_1 + m_2) * log(-sum_j)
  return(sum_)
}

# Simulation Dataset
get_dataset <- function(x){
  # initialise a dataframe
  cols = c("cause", "lifetime")
  df = data.frame(matrix(nrow = n_dpts, ncol = length(cols)))
  
  # assign column names
  colnames(df) = cols
  
  if (member == "weibull") {
    u1 <- rweibull(n_dpts, shape = x[4], scale = 1 / (x[2] ** (1 / x[4])))
    u2 <- rweibull(n_dpts, shape = x[4], scale = 1 / (x[3] ** (1 / x[4])))
    u3 <- rweibull(n_dpts, shape = x[4], scale = 1 / (x[1] ** (1 / x[4])))
  } else if (member == "gompertz") {
    u1 <- rgompertz(n_dpts, scale = x[4], shape = x[2] * x[4])
    u2 <- rgompertz(n_dpts, scale = x[4], shape = x[3] * x[4])
    u3 <- rgompertz(n_dpts, scale = x[4], shape = x[1] * x[4])
  } else if (member == "lomax") {
    u1 <- rlomax(n_dpts, scale = 1/x[4], shape3.q = x[2])
    u2 <- rlomax(n_dpts, scale = 1/x[4], shape3.q = x[3])
    u3 <- rlomax(n_dpts, scale = 1/x[4], shape3.q = x[1])
  } else {
    cat("member not found")
  }
  
  u_bind = cbind(u1, u2, u3)
  
  df$cause = apply(u_bind, 1, which.min)
  df[df$cause == 3, ]$cause = 0 # tie case
  
  df$lifetime = apply(u_bind, 1, min)
  
  # add some cases of censoring
  if (censoring) {
    df[df$lifetime > tau, ]$cause = 3
    df[df$lifetime > tau, ]$lifetime = tau  
  }
  
  # set num fails in each cause
  set_num_fails(df)
  
  return(df)
}

# Bootstrap Confidence Interval
boot_ci = function(x){
  ci_boot = unname(quantile(x, c(0.05, 0.95), na.rm = TRUE))
  return(ci_boot)
}

i = 1
mle_error = 0
start.time <- Sys.time()
while (i <= n_iter) {
  # track progress
  if(i %% 1 == 0){
    cat("iteration num: ", i)
    cat("\n")
  }
  
  # dataset
  df = get_dataset(true_param)
  cat("# fails:", c(m_0, m_1, m_2, m_c))
  cat("\n")
  
  # MLE
  MLE = try(
    {constrOptim(true_param[4],
                 LL,
                 lifetime_df = df,
                 grad = NULL,
                 control = list(fnscale=-1),
                 ui = 1,
                 ci = 0)},
    silent = TRUE
  )
  
  # if error then skip
  if(inherits(MLE, 'try-error')) {
    mle_error = mle_error + 1
    next
  }
  
  # log the values
  MLE_est[i] = MLE$value
  lam_est[i] = lam = MLE$par
  aph_0_est[i] = a0 = -m_0 / SS
  aph_1_est[i] = a1 = -m_1 / SS
  aph_2_est[i] = a2 = -m_2 / SS
  
  params = c(a0, a1, a2, lam)
  cat("estimated generating parameters:", params)
  cat("\n")
  
  # Divergence Cases
  if (abs(a0 - true_param[1]) > 10) {
    next
  }
  
  results_ci <- foreach(1:B,
                        .combine = 'rbind', 
                        .errorhandling = 'remove',
                        .packages = 'VGAM') %dopar% {
                          
                          b_df = get_dataset(params)
                          
                          # Bootstrap MLE
                          B_MLE = try(
                            {constrOptim(params[4],
                                         LL,
                                         lifetime_df = b_df,
                                         grad = NULL,
                                         control = list(fnscale=-1),
                                         ui = 1,
                                         ci = 0)},
                            silent = TRUE
                          )
                          
                          b0 = -m_0 / SS
                          b1 = -m_1 / SS
                          b2 = -m_2 / SS
                          
                          if (abs(b0 - params[1]) > 10) {
                            return(c(NA, NA, NA, NA))
                          }
                          
                          # if error then skip
                          if(inherits(B_MLE, 'try-error')) {
                            c(NA, NA, NA, NA)
                          } else {
                            c(b0, b1, b2, B_MLE$par)
                          }
                        }
  
  cat("\n")
  
  CI_aph0 = boot_ci(unname(results_ci[,1]))
  CI_aph1 = boot_ci(unname(results_ci[,2]))
  CI_aph2 = boot_ci(unname(results_ci[,3]))
  CI_lam = boot_ci(unname(results_ci[,4]))
  
  # finding the 95% CI
  ci_aph0_l[i] = CI_aph0[1]
  ci_aph0_u[i] = CI_aph0[2]
  ci_aph1_l[i] = CI_aph1[1]
  ci_aph1_u[i] = CI_aph1[2]
  ci_aph2_l[i] = CI_aph2[1]
  ci_aph2_u[i] = CI_aph2[2]
  ci_lam_l[i] = CI_lam[1]
  ci_lam_u[i] = CI_lam[2]
  
  i = i + 1
}

stopCluster(cls)

end.time <- Sys.time()
time.taken <- end.time - start.time

cat("time taken for simulation: ", time.taken)

## Average Estimate
AE_aph0 = sum(aph_0_est)/n_iter
AE_aph1 = sum(aph_1_est)/n_iter  
AE_aph2 = sum(aph_2_est)/n_iter
AE_lam = sum(lam_est)/n_iter

## MSE
MSE_aph0 = sum((aph_0_est - true_param[1])**2)/n_iter
MSE_aph1 = sum((aph_1_est - true_param[2])**2)/n_iter
MSE_aph2 = sum((aph_2_est - true_param[3])**2)/n_iter
MSE_lam = sum((lam_est - true_param[4])**2)/n_iter

## BIAS
BIAS_aph0 = AE_aph0 - true_param[1]
BIAS_aph1 = AE_aph1 - true_param[2]
BIAS_aph2 = AE_aph2 - true_param[3]
BIAS_lam = AE_lam - true_param[4]

## Empirical Variance
empvar_aph0 = MSE_aph0 - BIAS_aph0 ** 2
empvar_aph1 = MSE_aph1 - BIAS_aph1 ** 2
empvar_aph2 = MSE_aph2 - BIAS_aph2 ** 2
empvar_lam = MSE_lam - BIAS_lam ** 2

## CP and AL
for (i in 1:n_iter) {
  AL_aph0 = AL_aph0 + ci_aph0_u[i] - ci_aph0_l[i]
  AL_aph1 = AL_aph1 + ci_aph1_u[i] - ci_aph1_l[i]
  AL_aph2 = AL_aph2 + ci_aph2_u[i] - ci_aph2_l[i]
  AL_lam = AL_lam + ci_lam_u[i] - ci_lam_l[i]
  
  if ((ci_aph0_l[i] <= true_param[1]) && (ci_aph0_u[i] >= true_param[1])) {
    cp_aph0 = cp_aph0 + 1
  }
  
  if((ci_aph1_l[i] <= true_param[2]) && (ci_aph1_u[i] >= true_param[2])){
    cp_aph1 = cp_aph1 + 1
  }
  
  if((ci_aph2_l[i] <= true_param[3]) && (ci_aph2_u[i] >= true_param[3])){
    cp_aph2 = cp_aph2 + 1
  }
  
  if((ci_lam_l[i] <= true_param[4]) && (ci_lam_u[i] >= true_param[4])){
    cp_lam = cp_lam + 1
  }
}

metrics <- c(AL_aph0, AL_aph1, AL_aph2, AL_lam, cp_aph0, cp_aph1, cp_aph2, cp_lam)
metrics <- metrics / rep(n_iter, 8)
print(metrics)
