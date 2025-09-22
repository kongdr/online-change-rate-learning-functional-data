################# sparse - online #################
rm(list=ls())
source("./basic functions.R", encoding = 'UTF-8')

fun_mu <- function(t){ 5*sin(2*pi*t) }
fun_mu1 <- function(t){ 5*2*pi*cos(2*pi*t)}
fun_phi <- function(t){
  phi <- matrix(0,length(t),4)
  phi[,1] <- sqrt(2) * cos(2*pi*t)
  phi[,2] <- sqrt(2) * sin(2*pi*t)
  phi[,3] <- sqrt(2) * cos(4*pi*t)
  phi[,4] <- sqrt(2) * sin(4*pi*t)
  return(phi)
}

Mpc <- 4
lam <- ((1:Mpc)+1)^(-2)
a <- 0; b <- 1
EV1 <- 100; EV2 <- 50
eval_mu <- seq(a,b,length.out = EV1)
mu_true1 <- fun_mu1(eval_mu)
G <- 0.8
Kmax <- 20
sub.streams <- c(1,seq(10,Kmax,10))
njk_mean <- 8; njk_std <- 2
mk <- ceiling(rnorm(Kmax, 8, 2)); mk[1] <- 12
njk <- sapply(1:(2*8*Kmax),function(i){max(round(rnorm(1,njk_mean,njk_std)),2)})
njk <- njk[which(njk<=12 & njk>=4)]
njk <- njk[1:sum(mk)]

L1 <- 3 

N <- 0; mfull <- 0
res_theta_mu <- list()
res_theta_mu$centroids <- rep(0, L1)
res_theta_mu$P <- array(0, dim = c(5,5,EV1,L1))
res_theta_mu$q <- array(0, dim = c(5,EV1,L1))
res_sigma_mu1 <- list()
res_sigma_mu1$centroids <- rep(0, L1)
res_sigma_mu1$P <- array(0, dim = c(2,2,EV1,L1))
res_sigma_mu1$q <- array(0, dim = c(2,EV1,L1))
res_sigma_mu2 <- res_sigma_mu1
res_den1 <- list()
res_den1$centroids <- rep(0, L1)
res_den1$P <- array(0, dim = c(2,2,EV1,L1))
res_den1$q <- array(0, dim = c(2,EV1,L1))
res_mu <- list()
res_mu$centroids <- rep(0, L1)
res_mu$P <- array(0, dim = c(3,3,EV1,L1))
res_mu$q <- array(0, dim = c(3,EV1,L1))
r01 <- list()
r01$centroids <- rep(0, L1)
r01$P <- array(0, dim = c(2,2,EV1,L1))
r01$q <- array(0, dim = c(2,EV1,L1))
r10 <- list()
r10$centroids <- rep(0, L1)
r10$P <- array(0, dim = c(3,3,EV1,L1))
r10$q <- array(0, dim = c(3,EV1,L1))
r12 <- list()
r12$centroids <- rep(0, L1)
r12$P <- array(0, dim = c(4,4,EV1,L1))
r12$q <- array(0, dim = c(4,EV1,L1))

time <- c(); rss1 <- c(); h1 <- c()
sigma <- c(); theta<-c(); mus <- c()
sigma_mu1 <- c(); rss1i <- c()

set.seed(12345)
for(K in 1:Kmax){
  set.seed(12345 + K)
  njk1 <- njk[(mfull+1):(mfull+mk[K])]
  data <- gene_data(mk[K], njk1)
  x <- unlist(data$t); y <- unlist(data$y)
  NK <- length(y); N <- N + NK; mfull <- mfull + mk[K]
  
  t0 <- Sys.time()
  h_eps <- 0.33*mfull^(-3/20)
  eps_hat <- c()
  for(i in 1:mk[K]){
    ysm <- batch_LL(data$t[[i]], data$y[[i]], eval_mu, h_eps, njk1[i], 1)$est
    ysm <- approx(eval_mu, ysm, xout=data$t[[i]], method = 'linear')$y
    eps_hat <- c(eps_hat, data$y[[i]] - ysm)
  }
  idx <- which((eps_hat < mean(eps_hat)+3*sqrt(var(eps_hat)))&
               (eps_hat > mean(eps_hat)-3*sqrt(var(eps_hat))))
  eps_hat <- eps_hat[idx]
  sigma_eps <- mean(eps_hat^2)
  rm(data)
  
  h_theta_mu <- 0.55 * N^(-1/9)
  res_theta_mu <- online_LQuar4(x, y, eval_mu, h_theta_mu, L1, res_theta_mu, N, NK, 1)
  mu_thi_deri <- sapply(1:EV1, function(i){
    6*(solve(res_theta_mu$P[,,i,1]+diag(1e-12,5)) %*% matrix(res_theta_mu$q[,i,1],5,1))[4]
  })
  theta_mu <-  mean(res_theta_mu$P[1,1, ,1] * mu_thi_deri^2)
  den <- sapply(1:EV1, function(i){res_theta_mu$P[1,1,i,1]})
  
  h_sigma_mu <- G * N^(-1/5)
  res_sigma_mu1 <- online_LL(x, y, eval_mu, h_sigma_mu, L1, res_sigma_mu1, N, NK, 1)
  mu <- sapply(1:EV1, function(i){
    (solve(res_sigma_mu1$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_sigma_mu1$q[,i,1],2,1))[1]
  })
  res_den1 <- online_LL(den, y, eval_mu, h_sigma_mu, L1, res_den1, N, NK, 1)
  den1 <- sapply(1:EV1, function(i){res_den1$P[1,1,i,1]})
  mu_est<-approx(eval_mu, mu, xout=x, method = 'linear')$y
  r <- (y-mu_est)^2
  res_sigma_mu2 <- online_LL(x, r, eval_mu, h_sigma_mu, L1, res_sigma_mu2, N, NK, 1)
  rr <- sapply(1:EV1, function(i){
    (solve(res_sigma_mu2$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_sigma_mu2$q[,i,1],2,1))[1]
  })
  h_sigma_mu1 <- G * N^(-1/7)
  r10 <- online_LQuad2(x, r, eval_mu, h_sigma_mu1, L1, r10, N, NK, 1)
  r11 <- sapply(1:EV1, function(i){
    (solve(r10$P[,,i,1]+diag(1e-12,3)) %*% matrix(r10$q[,i,1],3,1))[2]
  })
  h_sigma_mu2 <- G * N^(-1/9)
  r12 <- online_LCub3(x, r, eval_mu, h_sigma_mu2, L1, r12, N, NK, 1)
  r2 <- sapply(1:EV1, function(i){
    2*(solve(r12$P[,,i,1]+diag(1e-12,4)) %*% matrix(r12$q[,i,1],4,1))[3]
  })
  r00 <- r-sigma_eps
  r01 <- online_LL(x, r00, eval_mu, h_sigma_mu, L1, r01, N, NK, 1)
  r0 <- sapply(1:EV1, function(i){
    (solve(r01$P[,,i,1]+diag(1e-12,2)) %*% matrix(r01$q[,i,1],2,1))[1]
  })
  s2 <- mean((r0*(den1)^2/den+r2*den+2*r11*den1)[6:95])
  if (s2<0) s2 <- abs(s2)
  sigma_mu <- 2.142857*mean(rr)+s2
  
  h_mu1 <- min((27 *5.444444* sigma_mu / theta_mu)^(1/7) * N^(-1/7), 1)
  res_mu <- online_LQuad2(x, y, eval_mu, h_mu1, L1, res_mu, N, NK, 1)
  mu1 <- sapply(1:EV1, function(i){
    (solve(res_mu$P[,,i,1]+diag(1e-12,3)) %*% matrix(res_mu$q[,i,1],3,1))[2]
  })
  t1 <- Sys.time()
  time <- c(time,difftime(t1,t0,units = 'secs'))
  rss1 <- c(rss1, mean((mu1 - mu_true1)^2))
  rss1i <-c(rss1i, mean(((mu1 - mu_true1)^2)[6:95])) 
  h1 <- c(h1, h_mu1)
  sigma <- c(sigma,sigma_eps)
  theta <-c(theta, theta_mu)
  mus <- cbind(mus, mu1)
  sigma_mu1 <- c(sigma_mu1,sigma_mu)
}
