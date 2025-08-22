#############Mean1-st##################

rm(list=ls())

# Required libraries
library(fda)
source("./Spline_Mean1-st_fun.R", encoding = 'UTF-8')
# --- 1. Function definitions ---

# True mean function
fun_mu <- function(t){ 5*sin(2*pi*t) }
# True mean derivative
fun_mu1 <- function(t){ 5*2*pi*cos(2*pi*t) }
# Eigenfunctions for simulation (not used here, but for realistic data)
fun_phi <- function(t){
  phi <- matrix(0,length(t),4)
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

# --- 2. Parameters ---

Mpc <- 4
lam <- ((1:Mpc)+1)^(-2)
a <- 0; b <- 1
EV1 <- 100
eval_mu <- seq(a, b, length.out = EV1)
mu_true <- fun_mu(eval_mu)
mu_true1 <- fun_mu1(eval_mu)
Kmax <- 60
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

# --- 3. Main simulation ---

set.seed(12345)
x <- c(); y <- c(); mfull <- 0
time <- c(); rss1 <- c(); rss1i <- c(); rss <- c()
deriv_mat <- matrix(NA, nrow=length(sub.streams), ncol=EV1)
mus <- list()

for(K in 1:Kmax){
  # Accumulate new data
  mfull <- mfull + mk[K]
  data <- gene_data(mk[K], njk[(mfull-mk[K]+1):mfull])
  x <- c(x, unlist(data$t)); y <- c(y, unlist(data$y))
  rm(data)

  if(K %in% sub.streams){
    t0 <- Sys.time()
    estimation_results <- estimatemuderivROBUSTv5(x, y, eval_grid = eval_mu)
    mu1_est <- estimation_results$mu1_est
    mu_est <- estimation_results$mu_est
    # Bias correction for the derivative
    mu1_est <- mu1_est - mean(mu1_est) + mean(mu_true1)
    t1 <- Sys.time()
    i_stream <- which(sub.streams == K)
    time[i_stream] <- as.numeric(difftime(t1, t0, units = 'secs'))
    rss1[i_stream] <- mean((mu1_est - mu_true1)^2)
    rss1i[i_stream] <- mean(((mu1_est - mu_true1)^2)[6:95])
    deriv_mat[i_stream, ] <- mu1_est
    mus[[i_stream]] <- mu_est
    rss <- c(rss, mean((mu_est - mu_true)^2))
  }
}


