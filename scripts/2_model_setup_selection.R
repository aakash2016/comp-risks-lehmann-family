# ----------------------------------------------------
# Licensed under the MIT License.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Model Selection in a 2-Model Setup

rm(list=ls())
options(warn=-1)

source("scripts/utils.R")
library(VGAM)

member = "gompertz" # specify the member of the Family
n_dpts = 150 # sample size
censoring = FALSE # whether to have censoring or not
tau = 2.0 # censoring factor

parent = "gompertz" # model selection
child = "weibull"
n_iter = 1000 # Number of monte carlo replications
true_param = c(1.81, 1.25, 1.42, 0.96)

# Estimates       
a_mle = b_mle = rep(0, n_iter)
aph0_est_a = aph1_est_a = aph2_est_a = lam_est_a = rep(0,n_iter)
aph0_est_b = aph1_est_b = aph2_est_b = lam_est_b = rep(0,n_iter)

pa = pb = 0 # probability if model selection

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

# Simulated Dataset
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

# Get Estimates using Optim
findMLE <- function() {
  MLE = try(
    {constrOptim(true_param[4],
                 LL,
                 grad = NULL,
                 control = list(fnscale=-1),
                 ui = 1,
                 ci = 0)},
    silent = TRUE
  )
  return(MLE)
}

i = 1
next_count = 0
start.time <- Sys.time()
while (i <= n_iter) {
  # track progress
  if(i %% 1 == 0){
    cat("iteration num: ", i)
    cat("\n")
  }
  
  # dataset
  member <<- parent
  df = get_dataset(true_param)
  lifetime_df <<- df # make df global
  
  la = lb = NULL
  
  # fit different models to generated data
  member <<- parent
  mle_a <- findMLE()
  print(inherits(mle_a, 'try-error'))
  if (!inherits(mle_a, 'try-error')) {
    la <- mle_a
  }
  
  # check if MLE doesnt exist
  if (is.null(la)) {
    next_count = next_count + 1
    next
  }
  
  a_mle[i] = ma = la$value
  aph0_est_a[i] = a0 = -m_0 / SS
  aph1_est_a[i] = a1 = -m_1 / SS
  aph2_est_a[i] = a2 = -m_2 / SS
  lam_est_a[i] = la$par
  cat("a",c(a0, a1, a2, la$par))
  cat('\n')
  
  ##
  member <<- child
  mle_b <- findMLE()
  print(inherits(mle_b, 'try-error'))
  if (!inherits(mle_b, 'try-error')) {
    lb <- mle_b
  }
  
  # check if MLE doesnt exist
  if (is.null(lb)) {
    next_count = next_count + 1
    next
  }
  
  b_mle[i] = mb = lb$value
  aph0_est_b[i] = a0 = -m_0 / SS
  aph1_est_b[i] = a1 = -m_1 / SS
  aph2_est_b[i] = a2 = -m_2 / SS
  lam_est_b[i] = lb$par
  cat("b",c(a0, a1, a2, lb$par))
  cat('\n')
  
  idx_max = which.max(c(ma, mb))
  if (idx_max == 1) {
    pa = pa + (1 / n_iter)
  } else {
    pb = pb + (1 / n_iter)
  }
  
  i = i + 1
}

end.time <- Sys.time()
time.taken <- end.time - start.time
cat("time taken for study: ", time.taken)

pest = c(pa, pb)
print(pest)
