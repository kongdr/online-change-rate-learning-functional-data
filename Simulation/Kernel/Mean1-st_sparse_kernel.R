################# sparse - kernel #################

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
eval_gam_vec <- seq(a,b,length.out = EV2)
eval_gam_mat <- cbind(rep(eval_gam_vec,each=EV2), rep(eval_gam_vec,EV2))
G <- 0.8
Kmax <- 20
sub.streams <- c(1,seq(10,Kmax,10)) 
mk <- ceiling(rnorm(Kmax, 8, 2)); mk[1] <- 12 
njk <- sapply(1:(2*8*Kmax),function(i){max(round(rnorm(1,njk_mean,njk_std)),2)})
njk <- njk[which(njk<=12 & njk>=4)]
njk <- njk[1:sum(mk)]

x <- c(); y <- c()
N <- 0; mfull <- 0
time <- c(); rss1 <- c(); h1 <- c()
sigma <- c(); mus <- c(); theta <- c()
sigma_mu1 <- c(); rss1i <- c()

for(K in 1:Kmax){
  set.seed(12345 + K)
  mfull <- mfull + mk[K]
  njk1 <- njk[1:mfull]
  data <- gene_data(mk[K], njk[(mfull-mk[K]+1):mfull])
  x <- c(x, unlist(data$t)); y <- c(y, unlist(data$y))
  N <- length(y)
  
  if(K %in% sub.streams){
    t0 <- Sys.time()
    # Estimate error variance
    h_eps <- 0.33*mfull^(-3/20)
    eps_hat <- c()
    for(i in 1:mfull){
      idx <- (sum(njk1[1:(i-1)])*(i>1)+1):sum(njk1[1:i])
      x1 <- x[idx]; y1 <- y[idx]
      ysm <- batch_LL(x1, y1, eval_mu, h_eps, njk1[i], 1)$est
      ysm <- approx(eval_mu, ysm, xout=x1, method = 'linear')$y
      eps_hat <- c(eps_hat, y1 - ysm)
    }
    idx <- which((eps_hat < mean(eps_hat)+3*sqrt(var(eps_hat)))&
                 (eps_hat > mean(eps_hat)-3*sqrt(var(eps_hat))))
    eps_hat <- eps_hat[idx]
    sigma_eps <- mean(eps_hat^2)
    # Estimate theta
    h_theta_mu <- 1.2*N^(-1/9)
    theta_mu <- batch_LQuar4(x, y, eval_mu, h_theta_mu, N, 1)
    # Estimate mean and its MSE
    h_sigma_mu <- G * N^(-1/5)
    res <- batch_LL(x, y, eval_mu, h_sigma_mu, N, 1)
    mu <- res$est
    den <- res$den
    res1 <- batch_LL(den, y, eval_mu, h_sigma_mu, N, 1)
    den1 <- res1$den
    mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
    r <- (y-mu_est)^2
    rr <- batch_LL(x, r, eval_mu, h_sigma_mu, N, 1)$est
    r0 <- r-sigma_eps
    r0 <- batch_LL(x, r0, eval_mu, h_sigma_mu, N, 1)$est
    h_sigma_mu1 <- G * N^(-1/7)
    r1v <- batch_LQuad2(x, r, eval_mu, h_sigma_mu1, N, 1)$est
    h_theta_mu2 <- G * N^(-1/9)
    r2 <- batch_LCub3(x, r, eval_mu, h_theta_mu2, N, 1)$sec_der
    s2 <- mean((r0*(den1)^2/den+r2*den+2*r1v*den1))
    if(s2 < 0) s2 <- abs(s2)
    sigma_mu <- 2.142857*mean(rr) + s2
    theta1 <-  mean((theta_mu$den[15:85] * (theta_mu$thi_der)[15:85]^2))
    h_mu1 <- (27*5.444444* sigma_mu /theta1)^(1/7) * N^(-1/7)
    mu1 <- batch_LQuad2(x, y, eval_mu, h_mu1, N, 1)$est
    t1 <- Sys.time()
    time <- c(time,difftime(t1,t0,units = 'secs'))
    rss1 <- c(rss1, mean((mu1 - mu_true1)^2))
    rss1i <- c(rss1i, mean(((mu1 - mu_true1)^2)[6:95]))
    h1 <- c(h1, h_mu1)
    sigma <- c(sigma, sigma_eps)
    mus <- cbind(mus, mu1)
    theta <- c(theta, theta_mu$theta)
    sigma_mu1 <- c(sigma_mu1,sigma_mu)
  }
}
