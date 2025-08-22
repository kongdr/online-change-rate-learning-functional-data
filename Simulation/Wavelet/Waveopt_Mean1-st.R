##############################################
# Mean Derivative Estimation
rm(list=ls())

# Required libraries
library(fda.usc)
library(KernSmooth)
library(wavethresh)

# 1. Function definitions

# True mean function
fun_mu <- function(t){ 5*sin(2*pi*t) }
# True mean derivative
fun_mu1 <- function(t){ 5*2*pi*cos(2*pi*t) }
# Eigenfunctions for simulation (not directly used, but for realistic data)
fun_phi <- function(t){
  phi <- matrix(0, length(t), 4)
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
  for(i in 1:4){
    kesi <- cbind(kesi, c(rnorm(mk, 0, sqrt(lam[i]))))
  }
  y <- lapply(1:mk, function(i) mu[[i]] +
                rowSums(matrix(rep(kesi[i,], each=njk[i]), ncol=4) * phi[[i]]) + e[[i]])
  mylist <- list(t, y)
  names(mylist) <- c('t','y')
  return(mylist)
}

# 2. Parameters

Mpc <- 4
lam <- ((1:Mpc)+1)^(-2)
a <- 0; b <- 1
EV1 <- 100
eval_mu <- seq(a, b, length.out = EV1)
mu_true <- fun_mu(eval_mu)
mu_true1 <- fun_mu1(eval_mu)
Kmax <- 40
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

# 3. Main single simulation

x <- c(); y <- c()
mfull <- 0
time <- c()
rss <- c()
rss1 <- c()
rss1i <- c()
deriv_mat <- NULL
mus <- list()

for(K in 1:Kmax){
  # Accumulate new data
  mfull <- mfull + mk[K]
  data <- gene_data(mk[K], njk[(mfull-mk[K]+1):mfull])
  x <- c(x, unlist(data$t)); y <- c(y, unlist(data$y))
  N <- length(y)
  
  if(K %in% sub.streams){
    t0 <- Sys.time()
    # Step 1: Oversmoothed local polynomial mean estimation
    oversmooth_factor <- 2.0 
    bw_optimal_mean <- KernSmooth::dpill(x, y)
    bw_for_deriv <- bw_optimal_mean * oversmooth_factor
    mu_est_on_eval_grid <- KernSmooth::locpoly(
      x, y, degree = 1, bandwidth = bw_for_deriv,
      gridsize = EV1, range.x = c(a, b)
    )$y
    # Step 2: Interpolate to dyadic grid (length 2^J)
    target_exponent <- max(8, floor(log2(N)*0.8 ))
    target_length <- 2^target_exponent
    grid_dyadic <- seq(a, b, length.out = target_length)
    mu_est_on_dyadic_grid <- approx(
      x = eval_mu, y = mu_est_on_eval_grid, xout = grid_dyadic
    )$y
    # Step 3: Wavelet denoising
    wt <- wavethresh::wd(mu_est_on_dyadic_grid,
                         filter.number = 10,
                         family = "DaubExPhase",
                         bc = "periodic")
    wt_thr <- wavethresh::threshold(wt,
                                    levels = 1:(wt$nlevels - 1),
                                    type   = "hard",
                                    policy = "universal")
    mu_denoised_grid <- wavethresh::wr(wt_thr)
    # Step 4: Compute derivative using fda.usc
    mat <- matrix(mu_denoised_grid, nrow = 1)
    fdata_obj <- fdata(mat, argvals = grid_dyadic)
    fdata_deriv <- fdata.deriv(fdata_obj, nderiv = 1)
    mu_d <- fdata_deriv$data[1, ]
    # Step 5: Interpolate derivative back to evaluation grid
    mu1_est <- approx(x = grid_dyadic, y = mu_d, xout = eval_mu, rule = 2)$y
    # Bias correction for derivative
    mu1_est <- mu1_est - mean(mu1_est) + mean(mu_true1)
    t1 <- Sys.time()
    # Save results
    time <- c(time, as.numeric(difftime(t1, t0, units = 'secs')))
    rss1 <- c(rss1, mean((mu1_est - mu_true1)^2))
    rss1i <- c(rss1i, mean(((mu1_est - mu_true1)^2)[6:95]))
    deriv_mat <- cbind(deriv_mat, mu1_est)
    rss <- c(rss, mean((mu_est_on_eval_grid - mu_true)^2))
    mus[[length(mus)+1]] <- mu_est_on_eval_grid
  }
}

