################## dense - kernel##################
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
gam_true <- c()
for(i in 1:EV2)
  for(j in 1:EV2){
    gam_true[(i-1)*EV2+j] <- sum(lam * fun_phi(eval_gam_vec[i]) * fun_phi(eval_gam_vec[j]))
  }
G <- 0.9
Kmax <- 300
sub.streams <- c(1,seq(20,Kmax,20))
mk <- rep(3, Kmax); mk[1] <- 10
njk_mean <- 20; njk_std <- 2
njk <- sapply(1:(2*5*Kmax),function(i){max(round(rnorm(1,njk_mean,njk_std)),2)})
njk <- njk[which(njk<=25 & njk>=15)]
njk <- njk[1:sum(mk)]

# --- ---- ---
x <- c(); y <- c()
N <- 0; mfull <- 0; N_gam <- 0
E1s <- c(); E2s <- c(); E3s <- c()
N_gam1 <- c(); N_3 <- c(); N_4 <- c()
time_band <- c()
deris <- list(); dens <- list()

for(K in 1:Kmax){
  set.seed(12345 + K)
  mfull <- mfull + mk[K]
  njk1 <- njk[1:mfull]
  data <- gene_data(mk[K], njk[(mfull-mk[K]+1):mfull])
  x <- c(x, unlist(data$t)); y <- c(y, unlist(data$y))
  N <- length(y)
  
  if(K %in% sub.streams){
    mu <- fun_mu(eval_mu)
    sigma_eps <- 0.5 # demo value
    mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
    res_gam_data <- gene_gam_data(N, njk1, mfull, x, y, mu_est)
    u <- res_gam_data$u; v <- res_gam_data$v
    N_gam <- sum(njk1 * (njk1 - 1))
    N_gam1 <- c(N_gam1,N_gam)
    N_3 <- c(N_3,sum(njk1 * (njk1 - 1)*(njk1 - 2)))
    N_4 <- c(N_4,sum(njk1 * (njk1 - 1)*(njk1 - 2)*(njk1 - 3)))
    
    t0<-Sys.time()
    # --- theta ---
    h_theta_gam <- G * N_gam^(-1/10)
    res <- batch_LQuar4(u, v, eval_gam_mat, h_theta_gam, N_gam, 2)
    dens <- c(dens,list(matrix(res$den,EV2,EV2)))
    deris <- c(deris,list(matrix(res$thi_der,EV2,EV2)))
    # --- sigma ---
    h_sigma_gam <- G * N_gam^(-1/6)
    res2 <- batch_LL(u, v, eval_gam_mat, h_sigma_gam, N_gam, 2)
    gam <- res2$est
    denm <- matrix(res2$den, EV2, EV2)
    rr <- v - sapply(1:N_gam, function(i){
      gam[which.min(abs(u[i,1]-eval_gam_mat[,1])+
                      abs(u[i,2]-eval_gam_mat[,2]))]
    })
    r <- rr^2
    r[r>(mean(r)+5*sqrt(var(r)))] <- mean(r)
    r[r<(mean(r)-5*sqrt(var(r)))] <- mean(r)
    r0 <- batch_LL(u, r, eval_gam_mat, h_sigma_gam, N_gam, 2)$est
    r <- matrix(r0, EV2, EV2)
    gamm <- matrix(gam, EV2, EV2)
    E1 <- mean(r) + 4*mean(diag(gamm))*sigma_eps + 2*sigma_eps^2+ mean(diag(r))
    E2 <- mean(r*sqrt(denm)) + 2*mean(diag(gamm*sqrt(denm)))*sigma_eps + mean(diag(r*sqrt(denm))) 
    h_den1 <- 0.5 * N^(-1/5)
    den1 <- batch_LL(sqrt(res2$den), v, eval_gam_mat[,1], h_den1, N, 1)$den
    den2 <- batch_LL(res2$den, v, eval_gam_mat[,1], h_den1, N, 1)$den
    V3 <- mean(r) - mean(gam^2)
    h_r <- G * N_gam^(-1/8)
    r1 <- batch_LQuad2(u, r0-gam^2, eval_gam_mat, h_r, N_gam, 2)$est
    r2 <- batch_LQuad2(u, rr, eval_gam_mat, h_r, N_gam, 2)$est
    E3 <- V3*mean(den1^2)+mean(r1*den2)+mean(r2^2)
    t1<-Sys.time()
    time_band <- c(time_band, difftime(t1,t0,units = 'secs'))
    E1s<-c(E1s,E1); E2s<-c(E2s,E2); E3s<-c(E3s,E3)
  }
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
                  b = -27*(N_4[1]/N_gam1[1]^(3/2))*E3s[k],
                  c = -27*2*2.142857*(N_3[1]/N_gam1[1]^(9/8))*E2s[k],
                  d = -27*2*1.285714*E1s[k])$root
}
