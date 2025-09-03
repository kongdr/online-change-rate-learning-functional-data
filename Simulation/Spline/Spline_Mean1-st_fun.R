
estimate_mu_deriv_CLASSIC <- function(x, y, eval_grid, 
                                      n_basis = 150, 
                                      k_folds = 8, 
                                      lambda_grid_size = 25) {

  clean_data <- na.omit(data.frame(x = x, y = y))
  x_clean <- clean_data$x
  y_clean <- clean_data$y
  n_obs <- length(x_clean)

  #  step 1: B spline
  full_range <- range(c(x_clean, eval_grid))
  
  bspline_basis <- create.bspline.basis(
    rangeval = full_range, 
    nbasis = n_basis,
    norder = 5
  )
  
  #  step 2: K-foldCV-lambda
  lambda_grid <- 10^seq(-3, -15, length.out = lambda_grid_size)

  set.seed(999) 
  folds_idx <- sample(1:k_folds, n_obs, replace = TRUE)
  
  cv_errors_matrix <- foreach(k = 1:k_folds, .combine = 'rbind', .packages = 'fda') %dopar% {
    x_train <- x_clean[folds_idx != k]; y_train <- y_clean[folds_idx != k]
    x_test <- x_clean[folds_idx == k]; y_test <- y_clean[folds_idx == k]
    
    fold_errors <- numeric(length(lambda_grid))
    
    for (i in 1:length(lambda_grid)) {
      lambda_val <- lambda_grid[i]
      fdPar_cv <- fdPar(fdobj = bspline_basis, Lfdobj = 3, lambda = lambda_val)
      smooth_cv <- smooth.basis(argvals = x_train, y = y_train, fdParobj = fdPar_cv)
      
      y_pred <- eval.fd(x_test, smooth_cv$fd)
      fold_errors[i] <- mean((y_test - y_pred)^2, na.rm = TRUE)
    }
    fold_errors
  }
  
  mean_cv_errors <- colMeans(cv_errors_matrix, na.rm = TRUE)
  best_lambda <- lambda_grid[which.min(mean_cv_errors)]
  
  #  step 3: fitting
  final_fdPar <- fdPar(fdobj = bspline_basis, Lfdobj = 3, lambda = best_lambda)
  final_smooth <- smooth.basis(argvals = x_clean, y = y_clean, fdParobj = final_fdPar)
  
  mu_fd <- final_smooth$fd
  mu_deriv_fd <- deriv.fd(mu_fd, Lfdobj = 1)
  mu1_est <- as.vector(eval.fd(evalarg = eval_grid, fdobj = mu_deriv_fd))
  mu_est <- as.vector(eval.fd(evalarg = eval_grid, fdobj = mu_fd))
  
  return(list(mu1_est = mu1_est, mu_est = mu_est))
}
