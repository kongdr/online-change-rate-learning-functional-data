estimatemuderivROBUSTv5 <- function(x, y, eval_grid, 
                                    k_folds = 5,
                                    nbasis_grid = c(8, 12, 16, 20, 24),
                                    lambda_grid = 10^seq(0, -10, length.out = 15),
                                    Lfd_order = 3) {
  
  
  n_obs <- length(x)
  param_grid <- expand.grid(nbasis = nbasis_grid, lambda = lambda_grid)
  n_combinations <- nrow(param_grid)
  
  
  set.seed(999)
  folds_idx <- sample(1:k_folds, n_obs, replace = TRUE)
  cv_errors_matrix <- matrix(NA, nrow = k_folds, ncol = n_combinations)
  full_range <- range(c(x, eval_grid))
  
  for (k in 1:k_folds) {
    x_train <- x[folds_idx != k]; y_train <- y[folds_idx != k]
    x_test <- x[folds_idx == k]; y_test <- y[folds_idx == k]
    
    for (j in 1:n_combinations) {
      nbasis_val <- param_grid$nbasis[j]
      lambda_val <- param_grid$lambda[j]
      
      bspline_basis_cv <- try(create.bspline.basis(rangeval = full_range, nbasis = nbasis_val, norder = 5), silent = TRUE)
      if (inherits(bspline_basis_cv, "try-error")) { next }
      
      fdPar_cv <- fdPar(fdobj = bspline_basis_cv, Lfdobj = Lfd_order, lambda = lambda_val)
      smooth_cv_fit <- try(smooth.basis(argvals = x_train, y = y_train, fdParobj = fdPar_cv), silent = TRUE)
      
      if (!inherits(smooth_cv_fit, "try-error")) {
        y_pred <- eval.fd(x_test, smooth_cv_fit$fd)
        cv_errors_matrix[k, j] <- mean((y_test - y_pred)^2, na.rm = TRUE)
      }
    }
  }
  
  mean_cv_errors <- colMeans(cv_errors_matrix, na.rm = TRUE)
  se_cv_errors <- apply(cv_errors_matrix, 2, function(col) sd(col, na.rm=T) / sqrt(sum(!is.na(col))))
  
  
  provisional_best_idx <- which.min(mean_cv_errors)
  min_error <- mean_cv_errors[provisional_best_idx]
  se_at_min <- se_cv_errors[provisional_best_idx]
  
  threshold <- min_error + se_at_min
  
  candidate_indices <- which(mean_cv_errors <= threshold)
  
  candidate_params <- param_grid[candidate_indices, ]
  sorted_candidates <- candidate_params[order(candidate_params$nbasis, -candidate_params$lambda), ]
  
  final_best_params <- sorted_candidates[1, ]
  best_nbasis <- final_best_params$nbasis
  best_lambda <- final_best_params$lambda
  
  final_bspline_basis <- create.bspline.basis(rangeval = full_range, nbasis = best_nbasis, norder = 5)
  final_fdPar <- fdPar(fdobj = final_bspline_basis, Lfdobj = Lfd_order, lambda = best_lambda)
  final_smooth <- smooth.basis(argvals = x, y = y, fdParobj = final_fdPar)
  
  mu_fd <- final_smooth$fd
  mu_deriv_fd <- deriv.fd(mu_fd, Lfdobj = 1)
  
  mu1est <- as.vector(eval.fd(evalarg = eval_grid, fdobj = mu_deriv_fd))
  muest <- as.vector(eval.fd(evalarg = eval_grid, fdobj = mu_fd))
  
  return(list(
    mu1_est = mu1est,
    mu_est = muest,
    debug_info = list(
      best_lambda = best_lambda,
      best_nbasis = best_nbasis,
      Lfd_order_used = Lfd_order,
      min_cv_error_1SE = mean_cv_errors[which(param_grid$nbasis == best_nbasis & param_grid$lambda == best_lambda)],
      df = final_smooth$df,
      k_folds = k_folds
    )
  ))
}
