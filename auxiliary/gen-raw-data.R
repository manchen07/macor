library(dplyr)
library(tidyr)
library(MASS)
library(mvtnorm)
library(metafor)
library(purrr)

# --------------------------------
# RAW DATA GENERATION 
# --------------------------------

set.seed(20200317)
k <-  100 # number of studies 
N <- 200  # N 

# generating sample sizes 

chisq <- rchisq(n = k, df = 3) 
samplesize <- round((N/2) * ((chisq - 3)/(sqrt(6)))  + N, digits = 0)

# population correlation matrix 

R <- matrix(cbind(1,   .18, .06,
                  .18, 1, .11,
                  .06, .11, 1), nrow = 3)

SD <- matrix(cbind(2, 0, 0, 
                   0, 2, 0,
                   0, 0, 7.67), nrow = 3)

# Sigma and Mean 

Sigma <- SD %*% R %*% t(SD) #covariance structure of observed variables
Mu_ctr <- c(10, 10, 50.73)

# for loop - generating control group individual data 

raw_dat <- list()

for (i in 1:k) {
  N <- samplesize[i] # Sample Size 
  
  # Control Group
  y_c <- rmvnorm(n = N, mean = Mu_ctr, sigma = Sigma) %>%
    data.frame() %>%
    mutate(group = "control",
           N = N)
  
  es <- 0.5 # this can be changed - we may use a vector of effect sizes (but here I used a single effect size)
  y_c_mean <- colMeans(y_c[,1:3])
  sd_c <- apply(y_c[,1:3], 2, sd)
  
  s_p1 <- sqrt(((N - 1)*sd_c[1]  + (N - 1)*sd_c[1])/( 2*N - 2))  # pooled sd for verbal fluency 
  s_p2 <- sqrt(((N - 1)*sd_c[2] + (N - 1)*sd_c[2])/( 2*N - 2))  # pooled sd for stroop color word
  s_p3 <- sqrt(((N - 1)*sd_c[3] + (N - 1)*sd_c[3])/( 2*N - 2))  # pooled sd for internalizing
  
  # Treatment Group  
  Mu_trt <- c(((es * s_p1) + y_c_mean[1]), ((es * s_p2) + y_c_mean[2]), ((es * s_p3) + y_c_mean[3]))
  
  y_t <- rmvnorm(n = N, mean = Mu_trt, sigma = Sigma) %>%
    data.frame() %>%
    mutate(group = "intervention",
           N = N)
  
  dat <- rbind(y_c, y_t)
  colnames(dat)[1:3] <- c("Verbal", "Stroop", "Internalzing")
  
  # binary variable 
  raw_dat[[i]] <- dat %>%
    mutate(Grp = ifelse(group == "control", 0, 1))

}

# Introducing missingness at variable-level
# a function to introduce missingness on either Verbal or Stroop variable
intr_miss <- function(x) {
  if(runif(1, 0, 1) <.4) { x[,c("Verbal")] <- NA}
  else{
    if(runif(1, 0, 1) <.4) {x[,c("Stroop")] <- NA}}
  return(x)
}

raw_dat_miss <- purrr::map(raw_dat, intr_miss) # please use --> raw_dat_miss

save(raw_dat_miss, file = "data/raw_dat.RData", compress = TRUE)