###########################
# Mean Derivative Estimation (Dense or Sparse) 

rm(list=ls())
library(fda)

# Function definitions

# True mean function
fun_mu <- function(t){ 5*sin(2*pi*t) }

# True mean derivative
fun_mu1 <- function(t){ 5*2*pi*cos(2*pi*t) }

# Eigenfunctions for data generation
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

# Parameters

Mpc <- 4
lam <- ((1:Mpc)+1)^(-2)

a <- 0; b <- 1
EV1 <- 100
eval_mu <- seq(a, b, length.out = EV1)
mu_true <- fun_mu(eval_mu)
mu_true1 <- fun_mu1(eval_mu)
Kmax <- 60

# Choose one of the following setups ("dense" or "sparse"):

## ------ Dense setup ------
njk_mean <- 25; njk_std <- 2
mk <- rep(3, Kmax); mk[1] <- 10
njk <- sapply(1:(2*5*Kmax), function(i){max(round(rnorm(1, njk_mean, njk_std)), 2)})
njk <- njk[which(njk<=30 & njk>=20)]
njk <- njk[1:sum(mk)]

## ------ Sparse setup ------
# mk_mean <- 18; mk_std <- 3
# mk <- ceiling(rnorm(Kmax, mk_mean, mk_std)); mk[1] <- 40
# njk_mean <- 8; njk_std <- 2
# njk <- sapply(1:(2*mk_mean*Kmax), function(i){max(round(rnorm(1, njk_mean, njk_std)), 2)})
# njk <- njk[which(njk<=11 & njk>=5)]
# njk <- njk[1:sum(mk)]

#Initialize

set.seed(12345)
sub_sds <- runif(Kmax) * 1e5
mfull <- 0
time_vec <- c()
rss_vec <- c()
rss1_vec <- c()
rss1i_vec <- c()
deriv_mat <- matrix(NA, nrow=Kmax, ncol=EV1)
mean_mat <- matrix(NA, nrow=Kmax, ncol=EV1)
step_i <- 0

# B-spline basis
nbasis <- 10
norder <- 4
basisObj <- create.bspline.basis(rangeval=c(a,b), nbasis=nbasis, norder=norder)
eval_basis <- function(t) eval.basis(t, basisObj)
eval_dbasis <- function(t) eval.basis(t, basisObj, Lfd=1)

# 2nd order penalty
D2 <- diff(diag(nbasis), differences=2)
Omega <- t(D2) %*% D2

# SGD hyperparameters
eta0 <- 0.1
alpha <- 0.3
gamma <- 1e-4
lambda <- 1e-4
Eg2 <- rep(0, nbasis)
rho <- 0.9
eps <- 1e-8
momentum <- rep(0, nbasis)
mu_mom <- 0.9
beta <- rep(0, nbasis)
beta_bar <- rep(0, nbasis)

#Warm-start: ridge regression on first W batches

t0 <- Sys.time()
W <- min(5, Kmax)
Phi0_list <- vector("list", W)
y0_list <- vector("list", W)
mfull0 <- 0
for(k0 in 1:W){
  set.seed(sub_sds[k0])
  njk0 <- njk[(mfull0+1):(mfull0+mk[k0])]
  dat0 <- gene_data(mk[k0], njk0)
  t0_all <- unlist(dat0$t)
  y0_all <- unlist(dat0$y)
  Phi0_list[[k0]] <- eval_basis(t0_all)
  y0_list[[k0]] <- y0_all
  mfull0 <- mfull0 + mk[k0]
}
Phi0 <- do.call(rbind, Phi0_list)
y0 <- unlist(y0_list)
lambda0 <- 1e-3
beta <- as.vector(solve(t(Phi0) %*% Phi0 + lambda0*diag(nbasis), t(Phi0) %*% y0))

#update

for(K in 1:Kmax){
  set.seed(sub_sds[K])
  njk1 <- njk[(mfull+1):(mfull+mk[K])]
  dat <- gene_data(mk[K], njk1)
  x <- unlist(dat$t)
  y <- unlist(dat$y)
  mfull <- mfull + mk[K]
  
  # Residuals
  Phi <- eval_basis(x)
  resid <- y - Phi %*% beta
  
  # Gradient and penalty
  grad <- -2 * t(Phi) %*% resid + 2 * lambda * (Omega %*% beta)
  
  # RMSProp
  Eg2 <- rho * Eg2 + (1 - rho) * (grad^2)
  adj <- grad / sqrt(Eg2 + eps)
  
  # Momentum
  momentum <- mu_mom * momentum + (1 - mu_mom) * adj
  
  # Update
  etaK <- eta0 * min(1, K/5) / (K+1)^alpha
  beta <- beta - etaK * momentum
  
  # Polyak average
  step_i <- step_i + 1
  beta_bar <- ((step_i-1) * beta_bar + beta) / step_i
  
  # Derivative estimate
  mu1_est <- as.vector(eval_dbasis(eval_mu) %*% beta_bar)
  # Mean function estimate
  mu_est <- as.vector(eval_basis(eval_mu) %*% beta_bar)
  
  # Record time, MSE, estimated curves
  t1 <- Sys.time()
  time_vec <- c(time_vec, as.numeric(difftime(t1, t0, units = 'secs')))
  rss1_vec <- c(rss1_vec, mean((mu1_est - mu_true1)^2))
  rss1i_vec <- c(rss1i_vec, mean(((mu1_est - mu_true1)^2)[6:95]))
  deriv_mat[K, ] <- mu1_est
  rss_vec <- c(rss_vec, mean((mu_est - mu_true)^2))
  mean_mat[K, ] <- mu_est
}
