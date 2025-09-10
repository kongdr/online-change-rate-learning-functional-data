###############################
# Covariance Derivative Estimation

rm(list=ls())

# Required libraries
library(fda)
library(Matrix)
library(foreach)
library(mgcv)
library(pracma)
source("./Spline_Cov1-st_fun.R", encoding = 'UTF-8')

# 1. Function definitions

# True mean
fun_mu <- function(t){ 5*sin(2*pi*t) }
# True mean derivative
fun_mu1 <- function(t){ 5*2*pi*cos(2*pi*t) }
# Eigenfunctions for simulation
fun_phi <- function(t){
  phi <- matrix(0, length(t), Mpc)
  phi[,1] <- sqrt(2) * cos(2*pi*t)
  phi[,2] <- sqrt(2) * sin(2*pi*t)
  phi[,3] <- sqrt(2) * cos(4*pi*t)
  phi[,4] <- sqrt(2) * sin(4*pi*t)
  return(phi)
}

# Data generator
gene_data <- function(mk, njk){
  t <- lapply(1:mk, function(i) runif(njk[i], a, b))
  e <- lapply(1:mk, function(i) rnorm(njk[i], 0, 0.5))
  mu <- lapply(t, fun_mu)
  phi <- lapply(t, fun_phi)
  kesi <- c()
  for(i in 1:Mpc){
    kesi <- cbind(kesi, c(rnorm(mk, 0, sqrt(lam[i]))))
  }
  y <- lapply(1:mk, function(i) mu[[i]] + 
                rowSums(matrix(rep(kesi[i,], each=njk[i]), ncol=Mpc) * phi[[i]]) + e[[i]])
  mylist <- list(t, y)
  names(mylist) <- c('t','y')
  return(mylist)
}

# 2. Parameters

Mpc <- 4
lam <- ((1:Mpc)+1)^(-2)
a <- 0; b <- 1
EV1 <- 100; EV2 <- 50
eval_mu <- seq(a, b, length.out = EV1)
mu_true <- fun_mu(eval_mu)
mu_true1 <- fun_mu1(eval_mu)
eval_gam_vec <- seq(a, b, length.out = EV2)

# True covariance derivative
fun_phi1 <- function(t){
  phi1 <- matrix(0, length(t), Mpc)
  phi1[,1] <- -2*pi*sqrt(2) * sin(2*pi*t)
  phi1[,2] <- 2*pi*sqrt(2) * cos(2*pi*t)
  phi1[,3] <- -4*pi*sqrt(2) * sin(4*pi*t)
  phi1[,4] <- 4*pi*sqrt(2) * cos(4*pi*t)
  return(phi1)
}
gam1_true <- c()
for(i in 1:EV2)
  for(j in 1:EV2){
    gam1_true[(i-1)*EV2+j] <- sum(lam * fun_phi1(eval_gam_vec[i]) * fun_phi(eval_gam_vec[j]))
  }
gam1_truem <- matrix(gam1_true, EV2, EV2)

#--------------dense------------------
Kmax <- 80
sub.streams <- c(1, seq(20, Kmax, 20)) # At which K to evaluate and record
set.seed(2024)
njk_mean <- 25; njk_std <- 2
mk <- rep(3, Kmax); mk[1] <- 10
njk <- sapply(1:(2*5*Kmax), function(i){max(round(rnorm(1, njk_mean, njk_std)), 2)})
njk <- njk[which(njk<=30 & njk>=20)]
njk <- njk[1:sum(mk)]

#-------------sparse------------------
#Kmax <- 40 #
#set.seed(2024)
#sds <- sample(1:2^20,R)
#mk_mean <- 18; mk_std <- 3
#mk <- ceiling(rnorm(Kmax, mk_mean, mk_std))
#mk[1] <- 40
#njk_mean <- 8; njk_std <- 2
#njk <- sapply(1:(2*mk_mean*Kmax),function(i){max(round(rnorm(1,njk_mean,njk_std)),2)})
#njk <- njk[which(njk<=11 & njk>=5)]
#njk <- njk[1:sum(mk)]

# 3. Main simulation

set.seed(12345)
x <- c(); y <- c(); id <- c()
mfull <- 0
time_cov <- c()
rss_cov <- c()
rss2i <- c()
cov_deriv_estimates <- list()
last_subject_id <- 0

# Loop through batches
for(K in 1:Kmax) {
  # 1. Determine number of new subjects for this batch
  num_new_subjects <- mk[K]
  start_idx <- mfull + 1
  end_idx <- mfull + num_new_subjects
  current_njk_subset <- njk[start_idx:end_idx]
  
  # 2. Generate and accumulate data
  data <- gene_data(num_new_subjects, current_njk_subset)
  x <- c(x, unlist(data$t))
  y <- c(y, unlist(data$y))
  
  # 3. Assign subject ids
  new_ids <- rep((last_subject_id + 1):(last_subject_id + num_new_subjects), times = current_njk_subset)
  id <- c(id, new_ids)
  
  # 4. Update counters
  mfull <- mfull + num_new_subjects
  last_subject_id <- last_subject_id + num_new_subjects
  
  # 5. If K is in sub.streams, perform covariance derivative estimation
  if (K %in% sub.streams) {
    t0 <- Sys.time()
    estimation_results <- estimate_cov_deriv_FULLY_ADAPTIVE_SPLINE_GAM_V19_5(
      id = id, 
      x = x, 
      y = y, 
      eval_grid = eval_gam_vec
    )
    estimated_cov_deriv_matrix <- estimation_results
    t1 <- Sys.time()
    time_cov <- c(time_cov, as.numeric(difftime(t1, t0, units = 'secs')))
    rss_cov <- c(rss_cov, mean((estimated_cov_deriv_matrix - gam1_truem)^2, na.rm = TRUE))
    cov_deriv_estimates[[as.character(K)]] <- estimated_cov_deriv_matrix
    rss2i <- c(rss2i, mean(((estimated_cov_deriv_matrix - gam1_truem)^2)[3:48, 3:48]))
    cat(sprintf("--- Finished covariance derivative estimation at K = %d, RSS = %.6f ---\n", K, tail(rss_cov, 1)))
  }
}
