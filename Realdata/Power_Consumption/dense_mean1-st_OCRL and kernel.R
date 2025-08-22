
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
    res_sigma_mu2$centroids <- rep(0, L1)
    res_sigma_mu2$P <- array(0, dim = c(2,2,EV1,L1))
    res_sigma_mu2$q <- array(0, dim = c(2,EV1,L1))
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
    
    time <- c(); h1 <- c(); mus1 <- c()
    sigma_mu1 <- c();theta_mu1 <- c();sigma <- c()
    start <- 0
    theta_mu <- 0; sigma_mu <- 0; sigma_eps <- 0
    h_old <- 1
  }
  # update
  for(year in 2013:2018){#data each year 
  load(paste(year,'dense_fdda.Rdata',sep=''))
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
    
    # sigma_eps
    {
      h_eps <- 0.33*mfull^(-3/20)
      eps_hat <- c()
      for(i in 1:mk){
        ysm <- batch_LL(data_$t[[i]], data_$y[[i]], eval_mu, h_eps, njk1[i], 1)$est
        ysm <- approx(eval_mu, ysm, xout=data_$t[[i]], method = 'linear')$y
        eps_hat <- c(eps_hat, data_$y[[i]] - ysm)
      }
      sigma_eps <- mean(eps_hat^2)
    }
    
    
    # theta
    {
      h_theta_mu <- G * N^(-1/9)
      res_theta_mu <- online_LQuar4(x, y, eval_mu, h_theta_mu, L1, 
                                    res_theta_mu, N, NK, 1)
      mu_thi_deri <- sapply(1:EV1, function(i){
        6*(solve(res_theta_mu$P[,,i,1]+diag(1e-12,5)) %*% matrix(res_theta_mu$q[,i,1],5,1))[4]
      })
      #theta_mu <-  mean(res_theta_mu$P[1,1, ,1] * mu_thi_deri^2)
      theta_mu <- (N-NK)/N * theta_mu + NK/N*(mean(res_theta_mu$P[1,1,8:92 ,1] * mu_thi_deri[8:92]^2)) 
      den <- sapply(1:EV1, function(i){res_theta_mu$P[1,1,i,1]})
    }
    
    # sigma
    {
      h_sigma_mu <- G * N^(-1/5)
      res_sigma_mu1 <- online_LL(x, y, eval_mu, h_sigma_mu, L1, res_sigma_mu1, N, NK, 1)
      mu <- sapply(1:EV1, function(i){
        (solve(res_sigma_mu1$P[,,i,1]+diag(1e-12,2)) %*% matrix(res_sigma_mu1$q[,i,1],2,1))[1]
      })
      #den <- sapply(1:EV1, function(i){res_sigma_mu1$P[1,1,i,1]})
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
      s2 <- mean((r0*(den1)^2/den+r2*den+2*r11*den1)[6:95])#
      if (s2>0){
        s2 <- s2
      }else{
        s2 <- abs(s2)
      }
      sigma_mu <- 2.142857*mean(rr)+s2
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
      sigma <- c(sigma,sigma_eps)
      sigma_mu1 <- c(sigma_mu1,sigma_mu)
      theta_mu1 <- c(theta_mu1,theta_mu)
      mus1 <- cbind(mus1, mu1)
    }
  }
  start <- Kmax
  }
  save(mus1, h1, time, sigma,sigma_mu1,theta_mu1, file='dense_res_mean1_online.Rdata')
}


###################### batch ########################
sub.streams <- c(214)#,1310,1889
### mean
{
  # initialize
  {
    x <- c(); y <- c()
    mus1 <- c(); h1 <- c(); sigma <- c(); time <- c()
    sigma_mu1 <- c();theta_mu1 <- c()
    N <- 0; mfull <- 0; njk <- c(); njk1 <- c();start <- 0
  }
  
  # update
  #for(year in 2013:2018){
  for(year in 2013){
    load(paste(year,'dense_fdda.Rdata',sep=''))
    #Kmax <- length(datafull$block) + start
    Kmax <- 214
  for(K in (start+1):Kmax){
    
    #print(paste('year=', year))
    
    t0 <- Sys.time()
    
    # generate data
    {
      data_ <- datafull$block[[K-start]]
      x <- c(x, unlist(data_$t))
      y <- c(y, unlist(data_$y))
      mk <- length(data_$y)
      #x <- c(x, unlist(sapply(1:length(data_), function(i){unlist(data_[[i]]$t)})))
      #y <- c(y, unlist(sapply(1:length(data_), function(i){unlist(data_[[i]]$y)})))
      njk <- c(njk,sapply(1:mk, function(i_mk){length(data_$t[[i_mk]])}))
      mfull <- mfull + mk
      njk1 <- njk[1:mfull]
      N <- length(y)
      #mk <- sum(sapply(1:length(data_), function(i){length(data_[[i]]$t)}))
      #mfull <- mfull + mk
      #njk1 <- c(njk1, unlist(sapply(1:length(data_), function(i){
      #  sapply(1:length(data_[[i]]$t), 
      #         function(j){length(data_[[i]]$t[[j]])})})))
    }
    
    if(K%in%sub.streams){
      print(paste('K=', K))
      # sigma_eps
      {
        h_eps <- 0.33*mfull^(-3/20)
        eps_hat <- c()
        for(i in 1:mfull){
          idx <- (sum(njk[1:(i-1)])*(i>1)+1):sum(njk[1:i])
          x1 <- x[idx]; y1 <- y[idx]
          ysm <- batch_LL(x1, y1, eval_mu, h_eps, njk1[i], 1)$est
          ysm <- approx(eval_mu, ysm, xout=x1, method = 'linear')$y
          eps_hat <- c(eps_hat, y1 - ysm)
        }
        # idx <- which((eps_hat < mean(eps_hat)+3*sqrt(var(eps_hat)))&
        #                (eps_hat > mean(eps_hat)-3*sqrt(var(eps_hat))))
        # eps_hat <- eps_hat[idx]
        sigma_eps <- mean(eps_hat^2)  
      }
      
      # theta, sigma, main regression
      {
        #if(year==1989){
        # theta
        h_theta_mu <- 0.6 * N^(-1/9)
        theta_mu <- batch_LQuar4(x, y, eval_mu, h_theta_mu, N, 1) 
        # sigma_mu
        h_sigma_mu <- G * N^(-1/5)
        res <- batch_LL(x, y, eval_mu, h_sigma_mu, N, 1)
        mu <- res$est
        den <- res$den
        h_sigma_mu <- G * N^(-1/5)
        res1 <- batch_LL(den, y, eval_mu, h_sigma_mu, N, 1)
        den1 <- res1$den
        mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
        r <- (y-mu_est)^2
        rr <- batch_LL(x, r, eval_mu, h_sigma_mu, N, 1)$est
        r0 <- r-sigma_eps
        r0 <- batch_LL(x, r0, eval_mu, h_sigma_mu, N, 1)$est
        h_sigma_mu1 <- G * N^(-1/7)
        r1 <- batch_LQuad2(x, r, eval_mu, h_sigma_mu1, N, 1)$est
        h_theta_mu2 <- G * N^(-1/9)
        r2 <- batch_LCub3(x, r, eval_mu, h_theta_mu2, N, 1)$sec_der
        h_theta_mu2 <- G * N^(-1/7)
        s2 <- mean((r0*(den1)^2/den+r2*den+2*r1*den1)[6:95])
        if(s2>0){
          s2 <- s2 
        }else{
          s2 <- abs(s2)
        }
        sigma_mu <- 2.142857*mean(rr) + s2
        #}
        theta1 <-  mean((theta_mu$den[6:90] * (theta_mu$thi_der)[6:90]^2))# [6:95]
        h_mu1 <- (27*5.444444* sigma_mu /theta1)^(1/7) * N^(-1/7)
        mu1 <- batch_LQuad2(x, y, eval_mu, h_mu1, N, 1)$est  
      }
      t1 <- Sys.time()
      {
        mus1 <- cbind(mus1, mu1)
        h1 <- c(h1,h_mu1)
        sigma <- c(sigma,sigma_eps)
        sigma_mu1 <- c(sigma_mu1,sigma_mu)
        theta_mu1 <- c(theta_mu1,theta1)
        time <- c(time,difftime(t1,t0,units = 'secs'))
      }
    }
    # save
  }
    start <- Kmax
  }
  save(mus1, h1, time,sigma,sigma_mu1,theta_mu1,  file='dense_res_mean1_batch.Rdata')
}



#####baseline
###################### batch ########################
sub.streams <- c(214,579,944,1310,1675,1889)#579,1310,1889
### mean
{
  # initialize
  {
    x <- c(); y <- c()
    mus1 <- c(); h1 <- c(); sigma <- c(); time <- c()
    sigma_mu1 <- c();theta_mu1 <- c()
    N <- 0; mfull <- 0; njk <- c(); njk1 <- c();start <- 0
  }
  
  # update
  for(year in 2013:2018){
    load(paste(year,'dense_fdda.Rdata',sep=''))
    Kmax <- length(datafull$block) + start
    for(K in (start+1):Kmax){
      
      #print(paste('year=', year))
      
      t0 <- Sys.time()
      
      # generate data
      {
        data_ <- datafull$block[[K-start]]
        x <- c(x, unlist(data_$t))
        y <- c(y, unlist(data_$y))
        mk <- length(data_$y)
        #x <- c(x, unlist(sapply(1:length(data_), function(i){unlist(data_[[i]]$t)})))
        #y <- c(y, unlist(sapply(1:length(data_), function(i){unlist(data_[[i]]$y)})))
        njk <- c(njk,sapply(1:mk, function(i_mk){length(data_$t[[i_mk]])}))
        mfull <- mfull + mk
        njk1 <- njk[1:mfull]
        N <- length(y)
        #mk <- sum(sapply(1:length(data_), function(i){length(data_[[i]]$t)}))
        #mfull <- mfull + mk
        #njk1 <- c(njk1, unlist(sapply(1:length(data_), function(i){
        #  sapply(1:length(data_[[i]]$t), 
        #         function(j){length(data_[[i]]$t[[j]])})})))
      }
      
      if(K%in%sub.streams){
        print(paste('K=', K))
        # sigma_eps
        {
          h_eps <- 0.33*mfull^(-3/20)
          eps_hat <- c()
          for(i in 1:mfull){
            idx <- (sum(njk[1:(i-1)])*(i>1)+1):sum(njk[1:i])
            x1 <- x[idx]; y1 <- y[idx]
            ysm <- batch_LL(x1, y1, eval_mu, h_eps, njk1[i], 1)$est
            ysm <- approx(eval_mu, ysm, xout=x1, method = 'linear')$y
            eps_hat <- c(eps_hat, y1 - ysm)
          }
          sigma_eps <- mean(eps_hat^2)  
        }
        
        # theta, sigma, main regression
        {
          h_theta_mu <- G * N^(-1/9)
          theta_mu <- batch_LQuar4(x, y, eval_mu, h_theta_mu, N, 1) 
          # sigma_mu
          h_sigma_mu <- G * N^(-1/5)
          res <- batch_LL(x, y, eval_mu, h_sigma_mu, N, 1)
          mu <- res$est
          den <- res$den
          h_sigma_mu <- G * N^(-1/5)
          res1 <- batch_LL(den, y, eval_mu, h_sigma_mu, N, 1)
          den1 <- res1$den
          mu_est <- approx(eval_mu, mu, xout=x, method = 'linear')$y
          r <- (y-mu_est)^2
          rr <- batch_LL(x, r, eval_mu, h_sigma_mu, N, 1)$est
          r0 <- r-sigma_eps
          r0 <- batch_LL(x, r0, eval_mu, h_sigma_mu, N, 1)$est
          h_sigma_mu1 <- G * N^(-1/7)
          r1 <- batch_LQuad2(x, r, eval_mu, h_sigma_mu1, N, 1)$est
          h_theta_mu2 <- G * N^(-1/9)
          r2 <- batch_LCub3(x, r, eval_mu, h_theta_mu2, N, 1)$sec_der
          h_theta_mu2 <- G * N^(-1/7)
          s2 <- mean((r0*(den1)^2/den+r2*den+2*r1*den1)[6:95])
          if(s2>0){
            s2 <- s2 
          }else{
            s2 <- abs(s2)
          }
          sigma_mu <- 2.142857*mean(rr) + s2
          theta1 <-  mean((theta_mu$den[11:90] * (theta_mu$thi_der)[11:90]^2))# [6:95]
          h_mu1 <- (27*5.444444* sigma_mu /theta1)^(1/7) * N^(-1/7)
          mu1 <- batch_LQuad2(x, y, eval_mu, h_mu1, N, 1)$est  
        }
        t1 <- Sys.time()
        {
          mus1 <- cbind(mus1, mu1)
          h1 <- c(h1,h_mu1)
          sigma <- c(sigma,sigma_eps)
          sigma_mu1 <- c(sigma_mu1,sigma_mu)
          theta_mu1 <- c(theta_mu1,theta1)
          time <- c(time,difftime(t1,t0,units = 'secs'))
        }
      }
      # save
    }
    start <- Kmax
  }
  save(mus1, h1, time,sigma,sigma_mu1,theta_mu1, file='dense_res_mean1_batch_baseline.Rdata')
}

