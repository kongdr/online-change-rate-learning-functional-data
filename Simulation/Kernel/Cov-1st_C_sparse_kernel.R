#################### sparse - kernel ####################
rm(list=ls())

source("./basic functions.R", encoding = 'UTF-8')

# --- Parameters and Functions ---
Mpc <- 4
lam <- ((1:Mpc)+1)^(-2)
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
gam_true <- c()
for(i in 1:EV2)
  for(j in 1:EV2){
    gam_true[(i-1)*EV2+j] <- sum(lam * fun_phi(eval_gam_vec[i]) * fun_phi(eval_gam_vec[j]))
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
x <- c(); y <- c(); z <- c()
N <- 0; mfull <- 0; N_gam <- 0
s <- c()
dens <- list(); deris <- list()

for(K in 1:Kmax){
  set.seed(12345 + K)
  mfull <- mfull + mk[K]
  njk1 <- njk[1:mfull]
  data <- gene_data(mk[K], njk[(mfull-mk[K]+1):mfull])
  x <- c(x, unlist(data$t)); y <- c(y, unlist(data$y)); z <- c(z, as.vector(unlist(data$z)))
  rm(data)
  N <- length(y)
  
  if(K %in% sub.streams){
    mu <- fun_mu(eval_mu)
    sigma_mu <- 0.5 # Or use a fixed value
    mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
    res_gam_data <- gene_gam_data(N, njk1, mfull, x, y, mu_est)
    u <- res_gam_data$u; v <- res_gam_data$v
    N_gam <- sum(njk1 * (njk1 - 1))
    
    # --- theta ---
    h_theta_gam <- G * N_gam^(-1/10)
    res <- batch_LQuar4(u, v, eval_gam_mat, h_theta_gam, N_gam, 2)
    dens <- c(dens, list(matrix(res$den, EV2, EV2)))
    deris <- c(deris, list(matrix(res$thi_der, EV2, EV2)))
    # --- sigma ---
    h_sigma_gam <- G * N_gam^(-1/6)
    gam <- batch_LL(u, v, eval_gam_mat, h_sigma_gam, N_gam, 2)$est
    r <- (v - sapply(1:N_gam, function(i){
      gam[which.min(abs(u[i,1]-eval_gam_mat[,1])+
                    abs(u[i,2]-eval_gam_mat[,2]))]
    }))^2
    r[r>(mean(r)+5*sqrt(var(r)))] <- mean(r)
    r[r<(mean(r)-5*sqrt(var(r)))] <- mean(r)
    r <- batch_LL(u, r, eval_gam_mat, h_sigma_gam, N_gam, 2)$est
    r <- matrix(r, EV2, EV2)
    gam <- matrix(gam, EV2, EV2)
    sigma_eps <- sigma_mu - mean(diag(gam))
    sigma_gam <- mean(r) + 4*mean(diag(gam))*sigma_eps + 2*sigma_eps^2 + mean(diag(r)) 
    s <- c(s, sigma_gam)
  }
}

# --- Calculate C parameter ---
th <- c()
for(k in 1:length(dens)){
  den <- dens[[k]]
  deri <- deris[[k]]
  th <- c(th, mean(den[5:45,5:45] * deri[5:45,5:45]^2))
}
C <- (27*5.444444*2*1.285714* s / th)^(1/8)

