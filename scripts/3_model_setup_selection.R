# ----------------------------------------------------
# Licensed under the MIT License.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# Model Selection in a 3-Model Setup

rm(list=ls())
options(warn=-1)

source("scripts/utils.R")
library(VGAM)

member = "gompertz" # specify the member of the Family
n_dpts = 300 # sample size
censoring = FALSE # whether to have censoring or not
tau = 2.0 # censoring factor

parent = "gompertz" # model selection
n_iter = 1000 # Number of monte carlo replications
true_param = c(1.63, 1.16, 1.39, 1.55)

# Estimates
weibull_mle = gompertz_mle = lomax_mle = rep(0, n_iter)
aph0_est_w = aph1_est_w = aph2_est_w = lam_est_w = rep(0,n_iter)
aph0_est_g = aph1_est_g = aph2_est_g = lam_est_g = rep(0,n_iter)
aph0_est_l = aph1_est_l = aph2_est_l = lam_est_l = rep(0,n_iter)

pw = pg = pl = 0 # probability if model selection

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
  
  lw = lg = ll <- NULL
  
  # fit different models to generated data
  # 1. WEIBULL
  member <<- "weibull"
  mle_w <- findMLE()
  print(inherits(mle_w, 'try-error'))
  if (!inherits(mle_w, 'try-error')) {
    lw <- mle_w
  }
  
  # check if MLE doesnt exist
  if (is.null(lw)) {
    next_count = next_count + 1
    next
  }
  
  weibull_mle[i] = mw = lw$value
  aph0_est_w[i] = -m_0 / SS
  aph1_est_w[i] = -m_1 / SS
  aph2_est_w[i] = -m_2 / SS
  lam_est_w[i] = lw$par
  
  # 2. GOMPERTZ
  member <<- "gompertz"
  mle_g <- findMLE()
  print(inherits(mle_g, 'try-error'))
  if (!inherits(mle_g, 'try-error')) {
    lg <- mle_g
  }

  # check if MLE doesnt exist
  if (is.null(lg)) {
    next_count = next_count + 1
    next
  }

  gompertz_mle[i] = mg = lg$value
  aph0_est_g[i] = -m_0 / SS
  aph1_est_g[i] = -m_1 / SS
  aph2_est_g[i] = -m_2 / SS
  lam_est_g[i] = lg$par

  # 3. LOMAX
  member <<- "lomax"
  mle_l <- findMLE()
  print(inherits(mle_l, 'try-error'))
  if (!inherits(mle_l, 'try-error')) {
    ll <- mle_l
  }

  # check if MLE doesnt exist
  if (is.null(ll)) {
    next_count = next_count + 1
    next
  }

  lomax_mle[i] = ml = ll$value
  aph0_est_l[i] = -m_0 / SS
  aph1_est_l[i] = -m_1 / SS
  aph2_est_l[i] = -m_2 / SS
  lam_est_l[i] = ll$par
  
  cat(lw$par, lg$par, ll$par)
  cat('\n')

  idx_max = which.max(c(mw, mg, ml))
  if (idx_max == 1) {
    pw = pw + (1 / n_iter)
  } else if (idx_max == 2) {
    pg = pg + (1 / n_iter)
  } else {
    pl = pl + (1 / n_iter)
  }

  i = i + 1
}

end.time <- Sys.time()
time.taken <- end.time - start.time
cat("time taken for study: ", time.taken)

pest = c(pw, pg, pl)
print(pest)
