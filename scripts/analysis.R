#$%
# Checking the saved results

rm(list=ls())
options(warn=-1)

n_iter = 500
true_param = c(0.85, 0.57, 0.74, 0.69)

df = read.csv("results/project_2/bootstrapping/lomax_400_0_para_simulation_.csv")

df = df[, -1]
colnames(df)[7: 14] <- c("ci_aph0_l",  "ci_aph0_u", "ci_aph1_l",
                         "ci_aph1_u",  "ci_aph2_l", "ci_aph2_u",  
                         "ci_lam_l", 'ci_lam_u')

## Average Estimate
AE_aph0 = sum(df$alpha_0)/n_iter
AE_aph1 = sum(df$alpha_1)/n_iter
AE_aph2 = sum(df$alpha_2)/n_iter
AE_lam = sum(df$lambda)/n_iter

## Relative MSE
MSE_aph0 = (sum((df$alpha_0 - true_param[1])**2)/n_iter) / true_param[1] **2
MSE_aph1 = (sum((df$alpha_1 - true_param[2])**2)/n_iter) / true_param[2] **2
MSE_aph2 = (sum((df$alpha_2 - true_param[3])**2)/n_iter) / true_param[3] **2
MSE_lam = (sum((df$lambda - true_param[4])**2)/n_iter) / true_param[4] **2

## Relative BIAS
BIAS_aph0 = (AE_aph0 - true_param[1]) / true_param[1]
BIAS_aph1 = (AE_aph1 - true_param[2]) / true_param[2]
BIAS_aph2 = (AE_aph2 - true_param[3]) / true_param[3]
BIAS_lam = (AE_lam - true_param[4]) / true_param[4]

# Average Length OF CIs.
AL_aph0 = AL_aph1 = AL_aph2 = AL_lam = 0 

# Coverage Percentage of CIs.
cp_aph0 = cp_aph1 = cp_aph2 = cp_lam = 0

# CP and AL
for (i in 1:n_iter) {
  AL_aph0 = AL_aph0 + df$ci_aph0_u[i] - df$ci_aph0_l[i]
  AL_aph1 = AL_aph1 + df$ci_aph1_u[i] - df$ci_aph1_l[i]
  AL_aph2 = AL_aph2 + df$ci_aph2_u[i] - df$ci_aph2_l[i]
  AL_lam = AL_lam + df$ci_lam_u[i] - df$ci_lam_l[i]
  
  if ((df$ci_aph0_l[i] <= true_param[1]) && (df$ci_aph0_u[i] >= true_param[1])) {
    cp_aph0 = cp_aph0 + 1
  }
  
  if((df$ci_aph1_l[i] <= true_param[2]) && (df$ci_aph1_u[i] >= true_param[2])){
    cp_aph1 = cp_aph1 + 1
  }
  
  if((df$ci_aph2_l[i] <= true_param[3]) && (df$ci_aph2_u[i] >= true_param[3])){
    cp_aph2 = cp_aph2 + 1
  }
  
  if((df$ci_lam_l[i] <= true_param[4]) && (df$ci_lam_u[i] >= true_param[4])){
    cp_lam = cp_lam + 1
  }
}

metrics <- c(AL_aph0, AL_aph1, AL_aph2, AL_lam, cp_aph0, cp_aph1, cp_aph2, cp_lam)
metrics <- metrics / rep(n_iter, 8)

cat("AL and CP", metrics)
cat("relative bias", c(BIAS_aph0, BIAS_aph1, BIAS_aph2, BIAS_lam))
cat("relative MSE", c(MSE_aph0, MSE_aph1, MSE_aph2, MSE_lam))
# End ##
