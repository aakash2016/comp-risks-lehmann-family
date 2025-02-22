# ----------------------------------------------------
# Licensed under the MIT License.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Sensitivity Analysis for general lehmann family ##

rm(list=ls())
options(warn=-1)

source("scripts/utils.R")
library(VGAM) # for rgompertz, rlomax

member <<- "weibull" # specify the member of the Family
n_dpts = 100 # sample size
n_iter = 10 # Number of monte carlo replications
true_param = c(0.85, 0.57, 0.74, 0.69) # generating parameter
censoring = FALSE # whether to have censoring or not
tau = 1.73 # censoring factor

perc_25 = 0.05 # found using 10k samples generated
perc_75 = 0.53
t_range = seq(perc_25, perc_75, by=0.25*(perc_75-perc_25))

# unit reliability
relative_ur_M1 = relative_ur_M2 = matrix(nrow=n_iter,ncol=5)

# Estimates
aph_0_est = aph_1_est = aph_2_est = lam_est = MLE_est = rep(0, n_iter)

# Convergence of MLEs
convergence = rep(0, n_iter)

# reliability of a unit
unit_reliability <- function(t, params) {
  return(exp(-(params[1] + params[2] + params[3]) * (t ** params[4])))
}
true_ur =  unit_reliability(t_range, true_param)

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

# Log Likelihood function M2: Independent We
LL_M2 <- function(x){
  sum_1 = 0
  sum_j = 0
  
  for(i in 1:nrow(lifetime_df)) {
    row <- lifetime_df[i, ]
    t = row$lifetime
    sum_j = sum_j + t ^ x
    
    if (row$cause != 3) {
      sum_1 = sum_1 + log(t)
    }
  }
  
  sum_ = (x - 1) * sum_1 + 
    (m_1 + m_2) * log(x / sum_j)
  
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
  
  # log the values Bivariate MLE
  MLE_est[i] = MLE$value
  lam_est[i] = lam = MLE$par
  aph_0_est[i] = a0 = -m_0 / SS
  aph_1_est[i] = a1 = -m_1 / SS
  aph_2_est[i] = a2 = -m_2 / SS
  
  params = c(a0, a1, a2, lam)
  print(params)
  
  ur_M1 = unit_reliability(t_range, params)
  relative_ur_M1[i,] = abs(true_ur - ur_M1) / true_ur
  
  # modify dataset and randomly split tites 
  lifetime_df[lifetime_df$cause == 0, ]$cause = (
    1 + as.integer(runif(length(lifetime_df[lifetime_df$cause == 0, ]$cause)) > 0.5)
  )
  
  lifetime_df <<- lifetime_df
  m = table(lifetime_df$cause)
  m_0 <<- 0
  m_1 <<- unname(m[1])
  m_2 <<- unname(m[2])
  print(c(m_0, m_1, m_2, m_c))
  
  # MLE: M2
  MLE_M2 = try(
    {constrOptim(true_param[4],
                 LL_M2,
                 grad = NULL,
                 control = list(fnscale=-1),
                 ui = 1,
                 ci = 0)},
    silent = TRUE
  )
  
  # log the values Bivariate MLE
  lam_M2 = MLE_M2$par
  a1_M2 = m_1 / sum(df$lifetime ** lam_M2)
  a2_M2 = m_2 / sum(df$lifetime ** lam_M2)
  
  # Divergence Cases
  if (abs(a0 - true_param[1]) > 100) {
    next
  }
  
  # if error then skip
  if(inherits(MLE, 'try-error') | inherits(MLE_M2, 'try-error')) {
    mle_error = mle_error + 1
    next
  }
  
  # t calculation
  # lms = unname(quantile(lifetime_df$lifetime, c(0.25, 0.75)))
  ur_M2 = unit_reliability(t_range, c(0, a1_M2, a2_M2, lam_M2))
  relative_ur_M2[i,] = abs(true_ur - ur_M2) / true_ur
  
  print(c(0, a1_M2, a2_M2, lam_M2))
  
  
  i = i + 1
}

end.time <- Sys.time()
time.taken <- end.time - start.time

cat("time taken for simulation: ", time.taken)

cat(t_range)
cat(colMeans(relative_ur_M1))
cat(colMeans(relative_ur_M2))
