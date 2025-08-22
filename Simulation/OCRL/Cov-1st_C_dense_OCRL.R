################## dense - online##################
rm(list=ls())

source("./basic functions.R", encoding = 'UTF-8')

# --- Parameters & Functions ---
Mpc <- 4
lam <- ((1:Mpc)+1)^(-2)
a <- 0; b <- 1
EV1 <- 100; EV2 <- 50
eval_mu <- seq(a, b, length.out = EV1)
eval_gam_vec <- seq(a, b, length.out = EV2)
eval_gam_mat <- cbind(rep(eval_gam_vec, each=EV2), rep(eval_gam_vec, EV2))
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
Kmax <- 300
sub.streams <- c(1,seq(20,Kmax,20))
mk <- rep(3, Kmax); mk[1] <- 10
njk_mean <- 20; njk_std <- 2
njk <- sapply(1:(2*5*Kmax),function(i){max(round(rnorm(1,njk_mean,njk_std)),2)})
njk <- njk[which(njk<=25 & njk>=15)]
njk <- njk[1:sum(mk)]

#=============

N <- 0; mfull <- 0; N_gam <- 0; N_3 <- 0; N_4 <- 0
N_gam1 <- c(); N_3c <- c(); N_4c <- c()
E1s5 <- c(); E2s5 <- c(); E3s5 <- c()
dens <- list(); deris <- list()
time_band <- c()

# --- Initialize online estimator objects ---
res_sigma_gam1 <- list()
res_sigma_gam1$centroids <- rep(0, 3)
res_sigma_gam1$P <- array(0, dim = c(3,3,EV2^2,3))
res_sigma_gam1$q <- array(0, dim = c(3,EV2^2,3))
res_sigma_gam2 <- res_sigma_gam1
res_theta_gam <- list()
res_theta_gam$centroids <- rep(0, 3)
res_theta_gam$P <- array(0, dim = c(15,15,EV2^2,3))
res_theta_gam$q <- array(0, dim = c(15,EV2^2,3))

for(K in 1:Kmax){
  set.seed(12345 + K)
  njk1 <- njk[(mfull+1):(mfull+mk[K])]
  data <- gene_data(mk[K], njk1)
  x <- unlist(data$t); y <- unlist(data$y)
  NK <- length(y); N <- N + NK; mfull <- mfull + mk[K]
  mu <- fun_mu(eval_mu)
  mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
  res_gam_data <- gene_gam_data(NK, njk1, mk[K], x, y, mu_est)
  u <- res_gam_data$u; v <- res_gam_data$v
  NK_gam <- sum(njk1 * (njk1 - 1))
  N_gam <- N_gam + NK_gam
  N_3 <- N_3+sum(njk1 * (njk1 - 1)*(njk1 - 2))
  N_4 <- N_4+sum(njk1 * (njk1 - 1)*(njk1 - 2)*(njk1 - 3))
  N_gam1 <- c(N_gam1,N_gam)
  N_3c <- c(N_3c,N_3)
  N_4c <- c(N_4c,N_4)

  t0<-Sys.time()
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
  dens<-c(dens,list(den)); deris<-c(deris,list(gam_thi_deri))
  
  # --- sigma ---
  h_sigma_gam <- G * N_gam^(-1/6)
  res_sigma_gam1 <- online_LL(u, v, eval_gam_mat, h_sigma_gam,
                              3, res_sigma_gam1, N_gam, NK_gam, 2)
  gam <- sapply(1:EV2^2, function(i){
    (solve(res_sigma_gam1$P[,,i,1]+diag(1e-12,3)) %*%matrix(res_sigma_gam1$q[,i,1],3,1))[1]
  })
  den <- sapply(1:EV2^2, function(i){res_sigma_gam1$P[1,1,i,1]})
  denm <- matrix(den,EV2, EV2)
  
  rr <- v - sapply(1:NK_gam, function(i){
    gam[which.min(abs(u[i,1]-eval_gam_mat[,1])+
                    abs(u[i,2]-eval_gam_mat[,2]))]
  })
  r <- rr^2
  r[r>(mean(r)+5*sqrt(var(r)))] <- mean(r)
  r[r<(mean(r)-5*sqrt(var(r)))] <- mean(r)
  res_sigma_gam2 <- online_LL(u, r, eval_gam_mat, h_sigma_gam,
                              3, res_sigma_gam2, N_gam, NK_gam, 2)
  r0 <- sapply(1:EV2^2, function(i){
    (solve(res_sigma_gam2$P[,,i,1]+diag(1e-12,3)) %*%matrix(res_sigma_gam2$q[,i,1],3,1))[1]
  })
  r <- matrix(r0,EV2, EV2)
  gamm <- matrix(gam, EV2, EV2)
  sigma_eps <- 0.5 # demo value
  E1_new <- mean(r) + 4*mean(diag(gamm))*sigma_eps + 2*sigma_eps^2+ mean(diag(r))
  E2_new <- mean(r*sqrt(denm)) + 2*mean(diag(gamm*sqrt(denm)))*sigma_eps + mean(diag(r*sqrt(denm))) 
  h_den1 <- 0.5 * N^(-1/5)
  # den1, den2, V3, r1, r2, E3_new computation as batch above
  den1 <- rep(1, EV2^2); den2 <- rep(1, EV2^2); V3 <- mean(r) - mean(gam^2)
  h_r <- G * N_gam^(-1/8)
  r1 <- rep(1, EV2^2); r2 <- rep(1, EV2^2)
  E3_new <- V3*mean(den1^2)+mean(r1*den2)+mean(r2^2)
  
  E1s5 <- c(E1s5, E1_new)
  E2s5 <- c(E2s5, E2_new)
  E3s5 <- c(E3s5, E3_new)
  t1<-Sys.time()
  time_band <- c(time_band, difftime(t1,t0,units = 'secs'))
}

# --- Calculate C parameter ---
theta1 <- sapply(seq_along(dens), function(k){
  den <- dens[[k]]
  deri <- deris[[k]]
  mean(den[6:45,6:45] * deri[6:45,6:45]^2)
})

f <- function(x,a,b,c,d){a*x^8+b*x^4+c*x+d}
C <- numeric(length(theta1))
for(k in seq_along(theta1)){
  C[k] <- uniroot(f, c(0,10),
                  a = 0.1836735*theta1[k],
                  b = -27*(N_4c[1]/N_gam1[1]^(3/2))*E3s5[k],
                  c = -27*2*2.142857*(N_3c[1]/N_gam1[1]^(9/8))*E2s5[k],
                  d = -27*2*1.285714*E1s5[k])$root
}
