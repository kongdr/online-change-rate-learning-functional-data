rm(list=ls())
### basic functions
source("./basic functions.R", encoding = 'UTF-8')


#### common parameters
{
  a <- 0; b <- 1
  EV1 <- 100; EV2 <- 50
  eval_mu <- seq(a,b,length.out = EV1)
  eval_gam_vec <- seq(a,b,length.out = EV2)
  eval_gam_mat <- cbind(rep(eval_gam_vec,each=EV2), rep(eval_gam_vec,EV2))
  G <- 0.8
  #Kmax <- 20
  #sub.streams <- c(1,seq(10,Kmax,10))
  R <- 1#112
}


###################### online ########################
### mean
{
  # initialize
  {
    N <- 0; mfull <- 0
    
    L1=5
    res_theta_mu <- list()
    res_theta_mu$centroids <- rep(0, L1)
    res_theta_mu$P <- array(0, dim = c(5,5,EV1,L1))
    res_theta_mu$q <- array(0, dim = c(5,EV1,L1))
    res_sigma_mu1 <- list()
    res_sigma_mu1$centroids <- rep(0, L1)
    res_sigma_mu1$P <- array(0, dim = c(2,2,EV1,L1))
    res_sigma_mu1$q <- array(0, dim = c(2,EV1,L1))
    res_sigma_mu2 <- res_sigma_mu1
    res_mu <- list()
    res_mu$centroids <- rep(0, L1)
    res_mu$P <- array(0, dim = c(3,3,EV1,L1))
    res_mu$q <- array(0, dim = c(3,EV1,L1))
    
    time <- c(); h1 <- c(); mus1 <- c()
    sigma_mu1 <- c();theta_mu1 <- c()
    start <- 0
    theta_mu <- 0; sigma_mu <- 0; sigma_eps <- 0
    h_old <- 1
  }
  # update

  for(year in 2013:2018){#data each year
    load(paste(year,'sparse_fdda.Rdata',sep=''))
    Kmax <- length(datafull$block) + start
    
  for(K in (start+1):Kmax){
    
    print(paste('K=',K))
    
    t0 <- Sys.time()
    # generate data
    {
      #i <- K - start
      #data_ <- data[[i]]
      data_ <- datafull$block[[K-start]]
      x <- unlist(data_$t); y <- unlist(data_$y)
      NK <- length(y); mk <- length(data_$y)
      njk1 <- sapply(1:mk, function(i_mk){length(data_$t[[i_mk]])})
      N <- N + NK; mfull <- mfull + mk
    }
    
    # theta
    {
      h_theta_mu <- G* N^(-1/9)
      res_theta_mu <- online_LQuar4(x, y, eval_mu, h_theta_mu, L1, 
                                    res_theta_mu, N, NK, 1)
      mu_thi_deri <- sapply(1:EV1, function(i){
        6*(solve(res_theta_mu$P[,,i,1]+diag(1e-12,5)) %*% matrix(res_theta_mu$q[,i,1],5,1))[4]
      })
      theta_mu <-  mean(res_theta_mu$P[1,1,6:95 ,1] * mu_thi_deri[6:95]^2)
    }
    
    # sigma
    {
      h_sigma_mu <- G * N^(-1/5)
      res_sigma_mu1 <- online_LL(x, y, eval_mu, h_sigma_mu, L1, res_sigma_mu1, N, NK, 1)
      mu <- sapply(1:EV1, function(i){
        (solve(res_sigma_mu1$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_sigma_mu1$q[,i,1],2,1))[1]
      })
      mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
      r <- (y-mu_est)^2
      res_sigma_mu2 <- online_LL(x, r, eval_mu, h_sigma_mu, L1, res_sigma_mu2, N, NK, 1)
      r <- sapply(1:EV1, function(i){
        (solve(res_sigma_mu2$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_sigma_mu2$q[,i,1],2,1))[1]
      })
      sigma_mu <- (N-NK)/N * sigma_mu + NK/N * mean(r)
    }
    
    # main regression
    {
      h_mu1 <- min((27 *5.444444* sigma_mu / theta_mu)^(1/7) * N^(-1/7),h_old)
      res_mu <- online_LQuad2(x, y, eval_mu, h_mu1, L1, res_mu, N, NK, 1)
      mu1 <- sapply(1:EV1, function(i){
        (solve(res_mu$P[,,i,1]+diag(1e-12,3)) %*% matrix(res_mu$q[,i,1],3,1))[2]
      })
      h_old <- h_mu1
    }
    t1 <- Sys.time()
    
    # combine
    {
      time <- c(time,difftime(t1,t0,units = 'secs'))
      h1 <- c(h1, h_mu1)
      sigma_mu1 <- c(sigma_mu1,sigma_mu)
      theta_mu1 <- c(theta_mu1,theta_mu)
      mus1 <- cbind(mus1, mu1)
    }
  }
  start <- Kmax
  }
  save(mus1, h1, time,sigma_mu1,theta_mu1, file='sparse_res_mean1_online.Rdata')
}


###################### batch ########################
sub.streams <- c(214,579,944,1310,1675,1889)#579,1310,1889
### mean
{
  # initialize
  {
    x <- c(); y <- c()
    mus1 <- c(); h1 <- c(); time <- c()
    sigma_mu1 <- c();theta_mu1 <- c()
    N <- 0; mfull <- 0; njk <- c(); start <- 0
    #theta1 <- 0; sigma_mu <- 0
  }
  
  # update
 
  for(year in 2013:2018){
    load(paste(year,'sparse_fdda.Rdata',sep=''))
    Kmax <- length(datafull$block) + start
    for(K in (start+1):Kmax){
    
    t0 <- Sys.time()
    
    # generate data
    {
      data_ <- datafull$block[[K-start]]
      x <- c(x, unlist(data_$t))
      y <- c(y, unlist(data_$y))
      mk <- length(data_$y)
      njk <- c(njk,sapply(1:mk, function(i_mk){length(data_$t[[i_mk]])}))
      mfull <- mfull + mk
      njk1 <- njk[1:mfull]
      N <- length(y)
    }
    
    if(K%in%sub.streams){
      print(paste('K=', K))
      # mean function estimation
      {
        h_theta_mu <- G * N^(-1/9)
        theta_mu <- batch_LQuar4(x, y, eval_mu, h_theta_mu, N, 1) 
        h_sigma_mu <- G * N^(-1/5)
        mu <- batch_LL(x, y, eval_mu, h_sigma_mu, N, 1)$est
        mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
        r <- (y-mu_est)^2
        r <- batch_LL(x, r, eval_mu, h_sigma_mu, N, 1)$est
        sigma_mu <- mean(r)
        theta1 <-  mean((theta_mu$den[6:95] * (theta_mu$thi_der)[6:95]^2))# [6:95]
        h_mu1 <- (27*5.444444 * sigma_mu /theta1)^(1/7) * N^(-1/7)#0.65
        mu1 <- batch_LQuad2(x, y, eval_mu, h_mu1, N, 1)$est
        
      }
      t1 <- Sys.time()
      # save
      {
        mus1 <- cbind(mus1, mu1)
        h1 <- c(h1,h_mu1)
        sigma_mu1 <- c(sigma_mu1,sigma_mu)
        theta_mu1 <- c(theta_mu1,theta1)
        time <- c(time,difftime(t1,t0,units = 'secs'))
      }
    }
    }
    start <- Kmax
  }
  save(mus1, h1, time,sigma_mu1,theta_mu1, file='sparse_res_mean1_batch.Rdata')
}

