################# dense #################
rm(list=ls())

# Function definitions
fun_mu <- function(t){ 5*sin(2*pi*t) }
fun_mu1 <- function(t){ 5*2*pi*cos(2*pi*t)}
fun_phi <- function(t){
  phi <- matrix(0,length(t),Mpc)
  phi[,1] <- sqrt(2) * cos(2*pi*t)
  phi[,2] <- sqrt(2) * sin(2*pi*t)
  phi[,3] <- sqrt(2) * cos(4*pi*t)
  phi[,4] <- sqrt(2) * sin(4*pi*t)
  return(phi)
}
gene_data <- function(mk, njk){
  t <- lapply(1:mk, function(i) runif(njk[i],a,b))
  e <- lapply(1:mk, function(i) rnorm(njk[i],0,0.5))
  mu <- lapply(t, fun_mu)
  phi <- lapply(t, fun_phi)
  kesi <- c()
  for(i in 1:Mpc){
    kesi <- cbind(kesi,c(rnorm(mk,0,sqrt(lam[i]))))
  }
  y <- lapply(1:mk, function(i) mu[[i]] + 
                rowSums(matrix(rep(kesi[i,],each=njk[i]),ncol=Mpc)
                        * phi[[i]]) + e[[i]])
  mylist <- list(t,y)
  names(mylist) <- c('t','y')
  return(mylist)
}

# Simulation parameters
Mpc <- 4
lam <- ((1:Mpc)+1)^(-2)
a <- 0; b <- 1
EV1 <- 100; EV2 <- 50
eval_mu <- seq(a,b,length.out = EV1)
mu_true1 <- fun_mu1(eval_mu)
eval_gam_vec <- seq(a,b,length.out = EV2)
eval_gam_mat <- cbind(rep(eval_gam_vec,each=EV2), rep(eval_gam_vec,EV2))
Kmax <- 40

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

# True derivative covariance function
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
    gam1_true[(i-1)*EV2+j] <- sum(lam * fun_phi1(eval_gam_vec[i])* fun_phi(eval_gam_vec[j]))
  }
gam1_truem <- matrix(gam1_true,EV2,EV2)

library(fda)
nbasis_cov   <- 5
norder_cov   <- 4
basisObj_cov <- create.bspline.basis(rangeval=c(a,b), nbasis=nbasis_cov, norder=norder_cov)
eval_basis_cov  <- function(t) eval.basis(t, basisObj_cov)
D2_cov    <- diff(diag(nbasis_cov), differences=2)
Omega_single <- t(D2_cov) %*% D2_cov
Omega_cov <- kronecker(diag(nbasis_cov), Omega_single) + kronecker(Omega_single, diag(nbasis_cov))

# SGD hyperparameters
eta0_cov   <- 0.05
alpha_cov <- 0.3
gamma_cov <- 1e-4
lambda_cov <- 1e-5
Eg2_cov      <- rep(0, nbasis_cov^2)
rho_cov <- 0.9
eps_cov <- 1e-8
momentum_cov <- rep(0, nbasis_cov^2)
mu_mom_cov <- 0.9

set.seed(12345)
mfull     <- 0
cov_time  <- c()
gam1s <- list() # Store estimated matrices
rss_gam1 <- c()
rssi_gam1 <- c()

# Warm-start: use first W batches for ridge regression initialization
W <- min(5, Kmax)
s_pairs_warmstart <- c()
t_pairs_warmstart <- c()
resid_prod_warmstart <- c()
mfull0 <- 0
for(k0 in 1:W){
  njk0 <- njk[(mfull0+1):(mfull0+mk[k0])]
  dat0 <- gene_data(mk[k0], njk0)
  x_list0 <- dat0$t
  x_unlist0 <- unlist(x_list0)
  mu_est0 <- approx(eval_mu, fun_mu(eval_mu), xout=x_unlist0, method = 'linear')$y
  resid_unlist0 <- unlist(dat0$y) - mu_est0
  resid_start_idx <- 1
  for (i in 1:mk[k0]) {
    n_i <- njk0[i]
    if (n_i < 2) {
      resid_start_idx <- resid_start_idx + n_i
      next
    }
    current_resids <- resid_unlist0[resid_start_idx:(resid_start_idx + n_i - 1)]
    current_times <- x_list0[[i]]
    indices <- combn(1:n_i, 2)
    s_pairs_warmstart <- c(s_pairs_warmstart, current_times[indices[1,]])
    t_pairs_warmstart <- c(t_pairs_warmstart, current_times[indices[2,]])
    resid_prod_warmstart <- c(resid_prod_warmstart, current_resids[indices[1,]] * current_resids[indices[2,]])
    resid_start_idx <- resid_start_idx + n_i
  }
  mfull0 <- mfull0 + mk[k0]
}

# Ridge regression for initialization
lambda0_cov <- 1e-4
Phi_s0 <- eval_basis_cov(s_pairs_warmstart)
Phi_t0 <- eval_basis_cov(t_pairs_warmstart)
Psi0 <- t(sapply(1:length(s_pairs_warmstart), function(i) kronecker(Phi_s0[i,], Phi_t0[i,])))
theta <- as.vector(solve(t(Psi0) %*% Psi0 + lambda0_cov * Omega_cov, t(Psi0) %*% resid_prod_warmstart))
theta_bar <- theta
step_i    <- 1

# Online update for each batch
for(K in 1:Kmax){
  t_start_online <- Sys.time()
  njk1 <- njk[(mfull+1):(mfull+mk[K])]
  dat  <- gene_data(mk[K], njk1)
  x_list <- dat$t
  y_list <- dat$y
  x_unlist <- unlist(x_list)
  y_unlist <- unlist(y_list)
  mfull <- mfull + mk[K]
  mu_est_K <- approx(eval_mu, fun_mu(eval_mu), xout=x_unlist, method = 'linear')$y
  resid_unlist <- y_unlist - mu_est_K
  
  s_pairs <- c(); t_pairs <- c(); resid_prod <- c()
  resid_start_idx <- 1
  for (i in 1:mk[K]) {
    n_i <- njk1[i]
    if (n_i < 2) {
      resid_start_idx <- resid_start_idx + n_i
      next
    }
    current_resids <- resid_unlist[resid_start_idx:(resid_start_idx + n_i - 1)]
    current_times <- x_list[[i]]
    indices <- combn(1:n_i, 2)
    s_pairs <- c(s_pairs, current_times[indices[1,]])
    t_pairs <- c(t_pairs, current_times[indices[2,]])
    resid_prod <- c(resid_prod, current_resids[indices[1,]] * current_resids[indices[2,]])
    resid_start_idx <- resid_start_idx + n_i
  }
  
  if(length(s_pairs) == 0) next
  
  Phi_s <- eval_basis_cov(s_pairs)
  Phi_t <- eval_basis_cov(t_pairs)
  Psi <- t(sapply(1:length(s_pairs), function(i) kronecker(Phi_s[i,], Phi_t[i,])))
  gamma_hat <- as.vector(Psi %*% theta)
  grad_cov <- -2 * t(Psi) %*% (resid_prod - gamma_hat)  + 2 * lambda_cov * (Omega_cov %*% theta)
  Eg2_cov   <- rho_cov * Eg2_cov + (1 - rho_cov) * (grad_cov^2)
  adj_cov   <- grad_cov / sqrt(Eg2_cov + eps_cov)
  momentum_cov <- mu_mom_cov * momentum_cov + (1 - mu_mom_cov) * adj_cov
  etaK_cov <- eta0_cov / ((1 + gamma_cov * K)^alpha_cov)
  theta <- theta - etaK_cov * momentum_cov
  theta_bar <- ((step_i-1) * theta_bar + theta) / step_i
  step_i   <- step_i + 1
  
  eval_dbasis_s <- eval.basis(eval_gam_vec, basisObj_cov, Lfd=1)
  eval_basis_t  <- eval.basis(eval_gam_vec, basisObj_cov, Lfd=0)
  gam1_est <- matrix(NA, nrow=EV2, ncol=EV2)
  for(i in 1:EV2){
    for(j in 1:EV2){
      Psi_st <- kronecker(eval_dbasis_s[i,], eval_basis_t[j,])
      gam1_est[i,j] <- as.vector(Psi_st %*% theta_bar)
    }
  }
  t_end_online <- Sys.time()
  cov_time <- c(cov_time, difftime(t_end_online, t_start_online, units = 'secs'))
  gam1s <- c(gam1s,list(gam1_est))
  rss_gam1 <- c(rss_gam1, mean((gam1_est - gam1_truem)^2))
  rssi_gam1 <- c(rssi_gam1, mean(((gam1_est - gam1_truem)^2)[3:48,3:48]))
}
