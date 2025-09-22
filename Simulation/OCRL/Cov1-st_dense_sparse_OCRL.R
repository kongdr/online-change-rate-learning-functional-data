################## dense & sparse - online ##################
rm(list=ls())

source("./basic functions.R", encoding = 'UTF-8')

# --- Parameters ---
Mpc <- 4
lam <- ((1:Mpc)+1)^(-2)
EV1 <- 100; EV2 <- 50
a <- 0; b <- 1
eval_mu <- seq(a,b,length.out = EV1)
eval_gam_vec <- seq(a,b,length.out = EV2)
eval_gam_mat <- cbind(rep(eval_gam_vec,each=EV2), rep(eval_gam_vec,EV2))
fun_phi <- function(t){
  phi <- matrix(0,length(t),Mpc)
  phi[,1] <- sqrt(2) * cos(2*pi*t)
  phi[,2] <- sqrt(2) * sin(2*pi*t)
  phi[,3] <- sqrt(2) * cos(4*pi*t)
  phi[,4] <- sqrt(2) * sin(4*pi*t)
  return(phi)
}
fun_phi1 <- function(t){
  phi1 <- matrix(0,length(t),Mpc)
  phi1[,1] <- -2*pi*sqrt(2) * sin(2*pi*t)
  phi1[,2] <- 2*pi*sqrt(2) * cos(2*pi*t)
  phi1[,3] <- -4*pi*sqrt(2) * sin(4*pi*t)
  phi1[,4] <- 4*pi*sqrt(2) * cos(4*pi*t)
  return(phi1)
}
G <- 0.9
Kmax <- 500
L2 <- 20 # Number of local polynomials for online estimator

#dense common parameters
mk <- rep(3, Kmax); mk[1] <- 10
njk_mean <- 20; njk_std <- 2
njk <- sapply(1:(2*5*Kmax),function(i){max(round(rnorm(1,njk_mean,njk_std)),2)})
njk <- njk[which(njk<=25 & njk>=15)]
njk <- njk[1:sum(mk)]

# sparse common parameters
mk_mean <- 18; mk_std <- 3
mk <- ceiling(rnorm(Kmax, mk_mean, mk_std)); mk[1] <- 40
njk_mean <- 8; njk_std <- 2
njk <- sapply(1:(2*mk_mean*Kmax),function(i){max(round(rnorm(1,njk_mean,njk_std)),2)})
njk <- njk[which(njk<=11 & njk>=5)]
njk <- njk[1:sum(mk)]

# --- Calculate true derivative covariance matrix ---
gam1_true <- c()
for(i in 1:EV2)
  for(j in 1:EV2){
    gam1_true[(i-1)*EV2+j] <- sum(lam * fun_phi1(eval_gam_vec[i]) * fun_phi(eval_gam_vec[j]))
  }
gam1_truem <- matrix(gam1_true,EV2,EV2)

# --- --- ---
N <- 0; mfull <- 0; N_gam <- 0; N_gam1 <- c()
C_old <- 1.2 # Fixed C for demonstration
res_gam <- list()
res_gam$centroids <- rep(0, L2)
res_gam$P <- array(0, dim = c(6,6,EV2^2,L2))
res_gam$q <- array(0, dim = c(6,EV2^2,L2)
rss2 <- c(); h2 <- c(); rss2i <- c(); time <- c(); gam1s <- list()

for(K in 1:Kmax){
  set.seed(12345 + K)
  njk1 <- njk[(mfull+1):(mfull+mk[K])]
  data <- gene_data(mk[K], njk1)
  x <- unlist(data$t); y <- unlist(data$y)
  NK <- length(y); N <- N + NK; mfull <- mfull + mk[K]
  mu <- rep(0, length(eval_mu)) # Use zero mean for demonstration
  mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
  res_gam_data <- gene_gam_data(NK, njk1, mk[K], x, y, mu_est)
  u <- res_gam_data$u; v <- res_gam_data$v
  NK_gam <- sum(njk1 * (njk1 - 1))
  N_gam <- N_gam + NK_gam
  N_gam1 <- c(N_gam1,N_gam)
  
  t0 <- Sys.time()
  if(K<=300){
    C_gam <- C
  }else{
    C_gam <- C[300]
  }
  h_gam <- C_gam * N_gam^(-1/8)
  res_gam <- online_LQuad2(u, v, eval_gam_mat, h_gam, L2, res_gam, N_gam, NK_gam, 2)
  gam1 <- sapply(1:EV2^2, function(i){
    (solve(res_gam$P[,,i,1]+diag(1e-12,6)) %*% matrix(res_gam$q[,i,1],6,1))[2]
  })
  t1 <- Sys.time()
  #save
  time <- c(time, difftime(t1,t0,units = 'secs'))
  rss2 <- c(rss2, mean((gam1 - gam1_true)^2))
  gam1m <- matrix(gam1,EV2,EV2)
  rss2i <- c(rss2i, mean(((gam1m - gam1_truem)^2)[3:48,3:48]))
  h2 <- c(h2, h_gam)
  gam1s <- c(gam1s,list(gam1m))
}

