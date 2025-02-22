# ----------------------------------------------------
# Licensed under the MIT License.
# Written by Aakash Agrawal (akash.gagrawal@gmail.com)
# ----------------------------------------------------
# %
# data analysis
# Diabetic Retinopathy data #

rm(list=ls())
options(warn=-1)
library(dplyr)
library(survival)
library(VGAM)
library(ggplot2)

source("scripts/utils.R")
member <<- "weibull" # specify the member of the Family
data = "diabetic" # diabetic or cencer
n_dpts = 50000 # sample size
norm_ = 365 # norm for diabetic dataset
rt = 1

# diabetic retinopathy dataset
if (data == "diabetic") {
  df = read.csv('datasets/data.csv', header = FALSE)
  lifetime_df = rename(df, c('id'='V1', 'lifetime'='V2', 'cause'='V3'))
  lifetime_df$lifetime  <- lifetime_df$lifetime / norm_ # lifetime in years
  lifetime_df$ind = 1
  lifetime_df$ind[lifetime_df$cause == 3] = 0
} else {
  print("dataset not found")
}

# check for bimodality in the histogram
hist(lifetime_df$lifetime, nclass=50)
set_num_fails(lifetime_df)

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

# MLE
MLE = optim(0.4, LL, control = list(fnscale=-1), method = "Nelder-Mead")
params = c(-m_0 / SS, -m_1 / SS, -m_2 / SS, MLE$par)

cat(MLE$value)
cat(params)

SweMinXY <- function(t) {
  return(exp(-(params[1] + params[2] + params[3])*(t**params[4])))
}

## Kaplan-Meier fit -- Non Parametric
t = seq(0.0, 4.55, 0.005)
Xdf = lifetime_df[lifetime_df$cause %in% c(0,1,3), ] %>% select('lifetime', 'ind')
Ydf = lifetime_df[lifetime_df$cause %in% c(0,2,3), ] %>% select('lifetime', 'ind')

KMfit_X = survfit(Surv(lifetime, ind) ~ 1, se.fit=F, data=Xdf)
KMfit_Y = survfit(Surv(lifetime, ind) ~ 1, se.fit=F, data=Ydf)
KMfit_MinXY = survfit(Surv(lifetime, ind) ~ 1, se.fit=F, data=lifetime_df)

# Convert the KMfit_MinXY to a data frame for ggplot
km_data <- data.frame(
  time = KMfit_MinXY$time,
  survival = KMfit_MinXY$surv
)

# Create the ggplot
ggplot(km_data, aes(x = time, y = survival)) +
  geom_step() + # Kaplan-Meier survival curve as a step plot
  stat_function(fun = SweMinXY, aes(color = "SweMinXY"), linetype = "solid") +
  labs(
    x = expression("T (in years)"),
    y = "Estimated Survival Probability",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.position = "none"
  )

# Profile Likelihood plot for different members
x <- seq(1.03, 2.2, 0.001)
ll_values <- LL(x)

# Create a data frame for ggplot
df <- data.frame(lambda = x, LL = ll_values)

# Create the plot using ggplot
ggplot(df, aes(x = lambda, y = LL)) +
geom_line(color = "blue") +
labs(
  x = expression(paste(lambda)),
  y = "Log-Likelihood (LL)",
  title = sprintf("Profile likelihood plot for Bivariate %s", member)
) +
theme(
  axis.title.x = element_text(margin = margin(t = 10)),
  axis.title.y = element_text(margin = margin(r = 10)),
  plot.title = element_text(hjust = 0.5)
)
