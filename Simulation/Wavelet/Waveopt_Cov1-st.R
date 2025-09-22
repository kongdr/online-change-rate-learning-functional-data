##############################################
# Covariance Derivative Estimation
rm(list=ls())

# Required libraries
library(wavethresh)
library(pracma)
library(dplyr)
library(mgcv)
library(data.table)
source("./Waveopt_Cov1-st_fun.R", encoding = 'UTF-8')

# 1. Function definitions

# True mean function
fun_mu <- function(t){ 5*sin(2*pi*t) }
# True mean derivative
fun_mu1 <- function(t){ 5*2*pi*cos(2*pi*t) }
# Eigenfunctions for simulation
fun_phi <- function(t){
  phi <- matrix(0,length(t),Mpc)
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
EV2 <- 50
eval_gam_vec <- seq(a, b, length.out = EV2)

# True covariance derivative
fun_phi1 <- function(t){
  phi1 <- matrix(0,length(t),Mpc)
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

Kmax <- 60
sub.streams <- c(1, seq(20, Kmax, 20))

## --- Dense setup ---
njk_mean <- 25; njk_std <- 2
mk <- rep(3, Kmax); mk[1] <- 10
njk <- sapply(1:(2*5*Kmax), function(i){max(round(rnorm(1, njk_mean, njk_std)),2)})
njk <- njk[which(njk<=30 & njk>=20)]
njk <- njk[1:sum(mk)]

## --- Sparse setup ---
# mk_mean <- 18; mk_std <- 3
# mk <- ceiling(rnorm(Kmax, mk_mean, mk_std)); mk[1] <- 40
# njk_mean <- 8; njk_std <- 2
# njk <- sapply(1:(2*mk_mean*Kmax), function(i){max(round(rnorm(1, njk_mean, njk_std)),2)})
# njk <- njk[which(njk<=11 & njk>=5)]
# njk <- njk[1:sum(mk)]

set.seed(12345)

# 3. Main simulation

mfull <- 0
time_cov <- c()
rss_cov <- c()
rss2i <- c()
cov_deriv_estimates <- list()
x_batch <- c()
y_batch <- c()
njk_cumulative <- c()

for(K in 1:Kmax) {
  # Accumulate new data
  njk_batch <- njk[(mfull + 1):(mfull + mk[K])]
  data_batch <- gene_data(mk[K], njk_batch)
  x_batch <- c(x_batch, unlist(data_batch$t))
  y_batch <- c(y_batch, unlist(data_batch$y))
  njk_cumulative <- c(njk_cumulative, njk_batch)
  rm(data_batch)
  mfull <- mfull + mk[K]
  
  # When K is in sub.streams, estimate covariance derivative
  if (K %in% sub.streams) {
    curve_indices <- rep(1:mfull, times = njk_cumulative)
    t_list <- split(x_batch, curve_indices)
    y_list <- split(y_batch, curve_indices)
    data_for_fpca <- list(t = t_list, y = y_list)
    
    cat(sprintf("--- Covariance derivative estimation at K = %d, total curves: %d ---\n", K, mfull))
    t0 <- Sys.time()
    # You need to ensure this function is sourced/loaded
    estimated_cov_deriv_matrix <- estimate_cov_deriv_WAVELET_GAM_V25_1(
      data_list = data_for_fpca,
      eval_grid = eval_gam_vec
    )
    t1 <- Sys.time()
    time_cov <- c(time_cov, as.numeric(difftime(t1, t0, units = 'secs')))
    rss_cov <- c(rss_cov, mean((estimated_cov_deriv_matrix - gam1_truem)^2, na.rm = TRUE))
    cov_deriv_estimates[[as.character(K)]] <- estimated_cov_deriv_matrix
    rss2i <- c(rss2i, mean(((estimated_cov_deriv_matrix - gam1_truem)^2)[3:48, 3:48]))
    cat(sprintf("--- Finished K = %d, RSS = %.6f ---\n", K, tail(rss_cov, 1)))
  }
}


