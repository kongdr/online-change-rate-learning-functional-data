#################### sparse - online ####################
rm(list=ls())

source("./basic functions.R", encoding = 'UTF-8')

# --- Parameters and Functions ---
Mpc <- 10
lam <- 0.4*(1:Mpc)^(-2)
a <- 0; b <- 1
EV1 <- 100; EV2 <- 50
eval_mu <- seq(a,b,length.out = EV1)
eval_gam_vec <- seq(a,b,length.out = EV2)
eval_gam_mat <- cbind(rep(eval_gam_vec,each=EV2), rep(eval_gam_vec,EV2))
fun_mu <- function(t){ 5*sin(2*pi*t) }
fun_phi <- function(t){
  phi <- matrix(0,length(t),Mpc)
  phi[,1] <- sqrt(2) * cos(2*pi*t)
  phi[,2] <- sqrt(2) * sin(2*pi*t)
  phi[,3] <- sqrt(2) * cos(4*pi*t)
  phi[,4] <- sqrt(2) * sin(4*pi*t)
  return(phi)
}
G <- 0.9
Kmax <- 20
sub.streams <- c(1,seq(10,Kmax,10))
mk_mean <- 18; mk_std <- 3
mk <- ceiling(rnorm(Kmax, mk_mean, mk_std)); mk[1] <- 40
njk_mean <- 8; njk_std <- 2
njk <- sapply(1:(2*mk_mean*Kmax),function(i){max(round(rnorm(1,njk_mean,njk_std)),2)})
njk <- njk[which(njk<=11 & njk>=5)]
njk <- njk[1:sum(mk)]

# --- --- ---
set.seed(12345)
N <- 0; mfull <- 0; N_gam <- 0
sigma_gam5 <- 0
res_sigma_gam1 <- list()
res_sigma_gam1$centroids <- rep(0, 3)
res_sigma_gam1$P <- array(0, dim = c(3,3,EV2^2,3))
res_sigma_gam1$q <- array(0, dim = c(3,EV2^2,3))
res_sigma_gam2 <- res_sigma_gam1
res_theta_gam <- list()
res_theta_gam$centroids <- rep(0, 3)
res_theta_gam$P <- array(0, dim = c(15,15,EV2^2,3))
res_theta_gam$q <- array(0, dim = c(15,EV2^2,3))
dens <- list(); deris <- list(); s5 <- c()

for(K in 1:Kmax){
  set.seed(12345 + K)
  njk1 <- njk[(mfull+1):(mfull+mk[K])]
  data <- gene_data(mk[K], njk1)
  x <- unlist(data$t); y <- unlist(data$y); z <- as.vector(unlist(data$z))
  rm(data)
  NK <- length(y); N <- N + NK; mfull <- mfull + mk[K]
  mu <- fun_mu(eval_mu)
  mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
  res_gam_data <- gene_gam_data(NK, njk1, mk[K], x, y, mu_est)
  u <- res_gam_data$u; v <- res_gam_data$v
  NK_gam <- sum(njk1 * (njk1 - 1))
  N_gam <- N_gam + NK_gam

  # --- theta ---
  h_theta_gam <- G * N_gam^(-1/10)
  res_theta_gam <- online_LQuar4(u, v, eval_gam_mat, h_theta_gam, 3, res_theta_gam, N_gam, NK_gam, 2)
  gam_thi_deri <- sapply(1:EV2^2, function(i){
    6*(sum((solve(res_theta_gam$P[,,i,1]+
                    diag(1e-12,15)) %*%
              matrix(res_theta_gam$q[,i,1],15,1))[c(7,8)]))
  })
  gam_thi_deri <- matrix(gam_thi_deri, EV2, EV2)
  den <- matrix(res_theta_gam$P[1,1, ,1], EV2, EV2)
  dens <- c(dens, list(den)); deris <- c(deris, list(gam_thi_deri))
  
  # --- sigma ---
  h_sigma_gam <- G * N_gam^(-1/6)
  res_sigma_gam1 <- online_LL(u, v, eval_gam_mat, h_sigma_gam,
                              3, res_sigma_gam1, N_gam, NK_gam, 2)
  gam <- sapply(1:EV2^2, function(i){
    (solve(res_sigma_gam1$P[,,i,1]+diag(1e-12,3)) %*% matrix(res_sigma_gam1$q[,i,1],3,1))[1]
  })
  r <- (v - sapply(1:NK_gam, function(i){
    gam[which.min(abs(u[i,1]-eval_gam_mat[,1])+
                  abs(u[i,2]-eval_gam_mat[,2]))]
  }))^2
  r[r>(mean(r)+5*sqrt(var(r)))] <- mean(r)
  r[r<(mean(r)-5*sqrt(var(r)))] <- mean(r)
  res_sigma_gam2 <- online_LL(u, r, eval_gam_mat, h_sigma_gam,
                              3, res_sigma_gam2, N_gam, NK_gam, 2)
  r <- sapply(1:EV2^2, function(i){
    (solve(res_sigma_gam2$P[,,i,1]+diag(1e-12,3)) %*% matrix(res_sigma_gam2$q[,i,1],3,1))[1]
  })
  r <- matrix(r,EV2,EV2)
  gam <- matrix(gam, EV2, EV2)
  sigma_eps <- 0.5 - mean(diag(gam)) # Or use a fixed value
  V1 <- mean(r) + 4*mean(diag(gam))*sigma_eps + 2*sigma_eps^2 + mean(diag(r))
  sigma_gam5 <- (N_gam-NK_gam)/N_gam * sigma_gam5 + NK_gam/N_gam * V1
  s5 <- c(s5, sigma_gam5)
}

# --- Calculate C parameter ---
th <- c()
for(k in 1:length(dens)){
  den <- dens[[k]]
  deri <- deris[[k]]
  th <- c(th, mean(den[6:45,6:45] *deri[6:45,6:45]^2))
}
C <- (27*5.444444*2*1.285714* s5 / th)^(1/8)
