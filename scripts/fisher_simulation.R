# ----------------------------------------------------
# Licensed under the MIT License.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Simulation for general lehmann family using FIM ##

rm(list=ls())
options(warn=-1)

source("scripts/utils.R")
library(VGAM) # for rgompertz, rlomax
library(pracma) # for inverse of matrix

h = 0.001 # finite difference
member <<- "lomax" # specify the member of the Family
n_dpts = 200 # sample size
n_iter = 1000 # Number of monte carlo replications
true_param = c(0.85, 0.57, 0.74, 0.69) # generating parameter
censoring = TRUE # whether to have censoring or not
tau = 1.73 # censoring factor

# Estimates
aph_0_est = aph_1_est = aph_2_est = lam_est = MLE_est = rep(0, n_iter)

# CIs
ci_aph0_l = ci_aph0_u = ci_aph1_l = ci_aph1_u = rep(0, n_iter)
ci_aph2_l = ci_aph2_u = ci_lam_l = ci_lam_u = rep(0, n_iter)

# Average Length OF CIs.
AL_aph0 = AL_aph1 = AL_aph2 = AL_lam = 0 

# Coverage Percentage of CIs.
cp_aph0 = cp_aph1 = cp_aph2 = cp_lam = 0

# Convergence of MLEs
convergence = rep(0, n_iter)

# Log Likelihood function
LL <- function(x){
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

func1 <- function(lambda) {
  sum_j = 0
  for(i in 1:nrow(lifetime_df)) {
    t <- lifetime_df[i, ]$lifetime
    sum_j = sum_j + log(S(t, lambda))
  }
  return(sum_j)
}

derivative_1 <- function(lambda) {
  d = (func1(lambda + 0.5 * h) - func1(lambda - 0.5 * h)) / h
  return(d)
}

func2 <- function(x) {
  alpha_0 <- x[1]
  alpha_1 <- x[2]
  alpha_2 <- x[3]
  lambda <- x[4]
  sum_ = 0
  sum_j = 0
  for(i in 1:nrow(lifetime_df)) {
    row <- lifetime_df[i, ]
    
    t = row$lifetime
    sum_j = sum_j + log(S(t, lambda))
    
    if (row$cause != 3) {
      sum_ = sum_ + log(f_pdf(t, lambda) / S(t, lambda))
    }
  }
  sum_ = sum_ + sum_j * (alpha_0 + alpha_1 + alpha_2)
  return(sum_)
}

derivative_2 <- function(x) {
  dd = (
    func2(c(x[1],x[2],x[3],x[4]+h)) + 
      func2(c(x[1],x[2],x[3],x[4]-h)) - 
      2 * func2(x)) / h ** 2
  return(dd)
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

i = 1
mle_error = 0
hessian_error = 0
start.time <- Sys.time()
while (i <= n_iter) {
  # track progress
  if(i %% 1 == 0){
    cat("iteration num: ", i)
    cat("\n")
  }
  
  # dataset
  df = get_dataset(true_param)
  print(c(m_0, m_1, m_2, m_c))
  
  lifetime_df <<- df # make df global
  
  # MLE
  MLE = try(
    {constrOptim(true_param[4],
                 LL,
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
  print(params)
  
  # Divergence Cases
  if (abs(a0 - true_param[1]) > 100) {
    next
  }
  
  p = derivative_1(lam)
  dd = derivative_2(params)
  hessMat = matrix( 
    c(
      -m_0 / a0 ** 2, 0, 0, p,
      0, -m_1 / a1 ** 2, 0, p,
      0, 0, -m_2 / a2 ** 2, p,
      p, p, p, dd
    ),
    nrow = 4,
    ncol = 4,
    byrow = TRUE
  )
  
  # convergence 0 indicates successful completion 
  # 1 indicates that the iteration limit maxit had been reached.
  convergence[i] = conv = MLE$convergence
  
  # covariance matrix is inverse of negative of hessian
  CV = inv(-hessMat)
  
  # if diagonal element less than 0
  if (any(diag(CV) < 0) | all(is.na(CV))) {
    hessian_error = hessian_error + 1
    print('hessian error')
    next
  }
  
  # finding the 95% CI
  ci_aph0_l[i] = a0 - 1.96 * sqrt(CV[1, 1])
  ci_aph0_u[i] = a0 + 1.96 * sqrt(CV[1, 1])
  ci_aph1_l[i] = a1 - 1.96 * sqrt(CV[2, 2])
  ci_aph1_u[i] = a1 + 1.96 * sqrt(CV[2, 2])
  ci_aph2_l[i] = a2 - 1.96 * sqrt(CV[3, 3])
  ci_aph2_u[i] = a2 + 1.96 * sqrt(CV[3, 3])
  ci_lam_l[i] = lam - 1.96 * sqrt(CV[4, 4])
  ci_lam_u[i] = lam + 1.96 * sqrt(CV[4, 4])
  
  i = i + 1
}

end.time <- Sys.time()
time.taken <- end.time - start.time

cat("time taken for simulation: ", time.taken)

## Average Estimate
AE_aph0 = sum(aph_0_est)/n_iter
AE_aph1 = sum(aph_1_est)/n_iter
AE_aph2 = sum(aph_2_est)/n_iter
AE_lam = sum(lam_est)/n_iter

# MSE
MSE_aph0 = sum((aph_0_est - true_param[1])**2)/n_iter
MSE_aph1 = sum((aph_1_est - true_param[2])**2)/n_iter
MSE_aph2 = sum((aph_2_est - true_param[3])**2)/n_iter
MSE_lam = sum((lam_est - true_param[4])**2)/n_iter

# BIAS
BIAS_aph0 = AE_aph0 - true_param[1]
BIAS_aph1 = AE_aph1 - true_param[2]
BIAS_aph2 = AE_aph2 - true_param[3]
BIAS_lam = AE_lam - true_param[4]

# Empirical Variance
empvar_aph0 = MSE_aph0 - BIAS_aph0 ** 2
empvar_aph1 = MSE_aph1 - BIAS_aph1 ** 2
empvar_aph2 = MSE_aph2 - BIAS_aph2 ** 2
empvar_lam = MSE_lam - BIAS_lam ** 2

# CP and AL
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
