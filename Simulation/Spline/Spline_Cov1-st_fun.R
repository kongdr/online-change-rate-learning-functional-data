estimate_cov_deriv_FULLY_ADAPTIVE_SPLINE_GAM_V19_6 <- function(id, x, y, eval_grid,
                                                               n_basis_mean = 80,
                                                               lambda_grid_size = 25, 
                                                               k_folds = 10,
                                                           
                                                               gam_k_base = 15,
                                                               gam_k_power = 0.1,
                                                               gam_k_max = 40,
                                                               seed = 12345) {
  
  cat("Starting Efficient Fully Adaptive Hybrid Spline-GAM Refinement Strategy v19.6 (High Performance Version)...\n")
  set.seed(seed)
  
  # --- Stage 1: Purely Adaptive Spline Fitting (Optimized Version V19.3) ---
  cat("--- STAGE 1: Purely Adaptive Spline Fitting ---\n")
  
  row_kronecker_product <- function(X, Y) {
    if (nrow(X) != nrow(Y)) stop("Matrices X and Y must have the same number of rows.")
    n_row <- nrow(X)
    n_col_X <- ncol(X)
    n_col_Y <- ncol(Y)
    # Expand X and Y matrices, then perform element-wise multiplication
    X_expanded <- X[, rep(1:n_col_X, each = n_col_Y), drop = FALSE]
    Y_expanded <- Y[, rep(1:n_col_Y, times = n_col_X), drop = FALSE]
    return(X_expanded * Y_expanded)
  }
  
  cat("--> Steps 1 & 2: GCV Mean Estimation and Raw Covariance Generation...\n")
  raw_range <- range(c(x, eval_grid)); range_span <- raw_range[2] - raw_range[1]
  full_range <- c(raw_range[1] - 0.001 * range_span, raw_range[2] + 0.001 * range_span)
  mean_basis <- create.bspline.basis(full_range, n_basis_mean, norder = 5)
  mean_lambda_grid <- 10^seq(-2, -12, length.out = 20)
  gcv_scores <- sapply(mean_lambda_grid, function(lambda) { smooth_obj <- smooth.basis(x, y, fdPar(mean_basis, 3, lambda)); if (is.null(smooth_obj$gcv)) 1e9 else smooth_obj$gcv })
  best_mean_lambda <- mean_lambda_grid[which.min(gcv_scores)]
  mean_fdPar <- fdPar(mean_basis, Lfdobj = 3, lambda = best_mean_lambda)
  mean_smooth <- smooth.basis(x, y, mean_fdPar); mu_est_func <- mean_smooth$fd
  residuals <- y - eval.fd(x, mu_est_func)
  data <- data.frame(id = id, t = x, res = residuals); unique_ids <- unique(id)
  num_curves <- length(unique_ids)
  s_raw_full <- c(); t_raw_full <- c(); cov_raw_full <- c()
  for (uid in unique_ids) {
    curve_data <- data[data$id == uid, ]; n_curve <- nrow(curve_data); if (n_curve < 2) next
    idx_pairs <- expand.grid(i = 1:n_curve, j = 1:n_curve)
    s_raw_full <- c(s_raw_full, curve_data$t[idx_pairs$i]); t_raw_full <- c(t_raw_full, curve_data$t[idx_pairs$j])
    cov_raw_full <- c(cov_raw_full, curve_data$res[idx_pairs$i] * curve_data$res[idx_pairs$j])
  }
  cat("--> Step 3: Intelligent Adaptive Data Compression...\n")
  num_points <- length(cov_raw_full); n_bins <- round(25 + (num_points/10000)^(1/3) * 15); n_bins <- max(25, min(100, n_bins))
  cat(sprintf("    - Original point cloud: %d. Adaptive bin count: %d x %d\n", num_points, n_bins, n_bins))
  bin_breaks <- seq(full_range[1], full_range[2], length.out = n_bins + 1)
  binned_df <- data.frame(s_idx = cut(s_raw_full, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE),
                          t_idx = cut(t_raw_full, breaks = bin_breaks, include.lowest = TRUE, labels = FALSE),
                          cov = cov_raw_full)
  binned_df_clean <- na.omit(binned_df); binned_agg <- aggregate(cov ~ s_idx + t_idx, data = binned_df_clean, FUN = mean)
  bin_centers <- (bin_breaks[-1] + bin_breaks[-(n_bins + 1)]) / 2
  s_model_data_orig <- bin_centers[binned_agg$s_idx]; t_model_data_orig <- bin_centers[binned_agg$t_idx]; cov_model_data_orig <- binned_agg$cov
  cat(sprintf("    - Data compression complete. Number of points used for weight generation: %d\n", length(cov_model_data_orig)))
  
  cat("--> Steps 4, 4.5, 4.6: Model Setup, Weight Generation, and Data Purification...\n")
  n_basis_cov <- round(12 + (length(cov_model_data_orig)/100)^(1/3) * 5); n_basis_cov <- max(12, min(45, n_basis_cov))
  s_basis <- create.bspline.basis(full_range, n_basis_cov, norder = 4); t_basis <- create.bspline.basis(full_range, n_basis_cov, norder = 4)
  P_s_pen <- Matrix(bsplinepen(s_basis, Lfdobj = 2), sparse = TRUE); P_t_pen <- Matrix(bsplinepen(t_basis, Lfdobj = 2), sparse = TRUE)
  P_s_ident <- Diagonal(n_basis_cov); P_t_ident <- Diagonal(n_basis_cov)
  # [Optimization] Adjust lambda search range to better adapt to large datasets
  log_lambda_min <- -8; log_lambda_max <- -3.0 - log10(length(cov_model_data_orig))/4; log_lambda_max <- min(-2.5, log_lambda_max)
  lambda_grid <- 10^seq(log_lambda_max, log_lambda_min, length.out = lambda_grid_size)
  lambda_init <- lambda_grid[floor(length(lambda_grid) * 0.6)]
  B_tensor_full_orig <- row_kronecker_product(eval.basis(s_model_data_orig, s_basis), eval.basis(t_model_data_orig, t_basis))
  BtB_full_orig <- crossprod(B_tensor_full_orig); Bty_full_orig <- crossprod(B_tensor_full_orig, cov_model_data_orig)
  R_init <- lambda_init * (kronecker(P_s_pen, P_t_ident) + kronecker(P_s_ident, P_t_pen))
  LHS_init <- BtB_full_orig + R_init
  init_coefs_matrix <- matrix(as.vector(solve(LHS_init, Bty_full_orig)), nrow = n_basis_cov)
  C_init_bifd <- bifd(coef = init_coefs_matrix, sbasisobj = s_basis, tbasisobj = t_basis)
  n_model_pts <- length(s_model_data_orig)
  s_rough_vec <- numeric(n_model_pts); t_rough_vec <- numeric(n_model_pts)
  for(i in 1:n_model_pts) {
    s_rough_vec[i] <- eval.bifd(s_model_data_orig[i], t_model_data_orig[i], C_init_bifd, sLfd = 2, tLfd = 0)
    t_rough_vec[i] <- eval.bifd(s_model_data_orig[i], t_model_data_orig[i], C_init_bifd, sLfd = 0, tLfd = 2)
  }
  point_weights_orig <- abs(s_rough_vec) + abs(t_rough_vec)
  valid_idx <- !is.na(point_weights_orig)
  if(sum(!valid_idx) > 0) { cat(sprintf("    - Warning: Detected and removed %d data points with NA weights due to boundary effects.\n", sum(!valid_idx))) }
  s_model_data <- s_model_data_orig[valid_idx]; t_model_data <- t_model_data_orig[valid_idx]
  cov_model_data <- cov_model_data_orig[valid_idx]; point_weights_raw <- point_weights_orig[valid_idx]
  point_weights <- pmax(point_weights_raw, quantile(point_weights_raw, 0.1)); point_weights <- point_weights / mean(point_weights)
  cat(sprintf("    - Data purification complete. Final number of points used for modeling: %d\n", length(cov_model_data)))
  
  cat("--> Step 5: Weighted CV and Final Fitting...\n")
  set.seed(seed); folds_idx <- sample(1:k_folds, length(cov_model_data), replace = TRUE)
  
  # [Optimization] Efficient CV: Precompute basis matrices for all points, extract via indices in loop
  B_s_full <- eval.basis(s_model_data, s_basis)
  B_t_full <- eval.basis(t_model_data, t_basis)
  B_tensor_full_unweighted <- row_kronecker_product(B_s_full, B_t_full)
  
  # [Optimization] Explicitly load packages in foreach to enhance parallel robustness
  cv_errors_list <- foreach(k = 1:k_folds, .packages=c("Matrix")) %dopar% {
    train_idx <- which(folds_idx != k); test_idx <- which(folds_idx == k)
    if (length(test_idx) == 0 || length(train_idx) == 0) return(NA)
    
    # [Optimization] Directly obtain training/test data and basis matrices via indices, no need to recalculate
    cov_train <- cov_model_data[train_idx]; w_train <- point_weights[train_idx]
    cov_test <- cov_model_data[test_idx]
    
    B_tensor_train_unweighted_k <- B_tensor_full_unweighted[train_idx, , drop=FALSE]
    B_tensor_test_k <- B_tensor_full_unweighted[test_idx, , drop=FALSE]
    
    W_train_mat <- Diagonal(x = w_train)
    B_tensor_train <- W_train_mat %*% B_tensor_train_unweighted_k
    cov_train_weighted <- w_train * cov_train
    
    BtB_train <- crossprod(B_tensor_train, B_tensor_train_unweighted_k) # More efficient crossprod(W*B, B)
    Bty_train <- crossprod(B_tensor_train, cov_train) 
    
    fold_errors <- numeric(length(lambda_grid))
    
    # Define penalty matrix components
    P_kron_sum <- kronecker(P_s_pen, P_t_ident) + kronecker(P_s_ident, P_t_pen)
    
    for (i in 1:length(lambda_grid)) {
      lambda_val <- lambda_grid[i]
      R <- lambda_val * P_kron_sum
      LHS <- BtB_train + R
      # Use `try` to avoid interruption due to singular matrices
      coefs <- try(as.vector(solve(LHS, Bty_train)), silent = TRUE)
      if (inherits(coefs, "try-error")) {
        fold_errors[i] <- 1e9 # Assign a large error if solving fails
        next
      }
      cov_pred <- as.vector(B_tensor_test_k %*% coefs)
      fold_errors[i] <- mean((cov_test - cov_pred)^2, na.rm = TRUE)
    }
    fold_errors
  }
  cv_errors_matrix <- do.call(rbind, cv_errors_list); mean_cv_errors <- colMeans(cv_errors_matrix, na.rm = TRUE)
  best_lambda <- if (all(is.nan(mean_cv_errors)) || all(mean_cv_errors > 1e8)) {
    lambda_grid[floor(length(lambda_grid) / 2)] 
  } else {
    lambda_grid[which.min(mean_cv_errors)]
  }
  
  W_full_mat <- Diagonal(x = point_weights)
  B_tensor_full <- W_full_mat %*% B_tensor_full_unweighted
  BtB_full <- crossprod(B_tensor_full, B_tensor_full_unweighted) # More efficient
  Bty_full <- crossprod(B_tensor_full, cov_model_data)
  
  R_final <- best_lambda * (kronecker(P_s_pen, P_t_ident) + kronecker(P_s_ident, P_t_pen))
  LHS_final <- BtB_full + R_final; final_coefs_matrix <- matrix(as.vector(solve(LHS_final, Bty_full)), nrow = n_basis_cov)
  C_est_bifd <- bifd(coef = final_coefs_matrix, sbasisobj = s_basis, tbasisobj = t_basis)
  cat(sprintf("--- STAGE 1 Complete. Optimal lambda = %s ---\n", format(best_lambda, scientific=T)))
  
  # --- Stage 2: Fully Adaptive GAM Final Refinement ---
  cat("--- STAGE 2: Fully Adaptive GAM Polishing ---\n")
  
  # 5.5a. Evaluate spline fitting results on a fine intermediate grid
  J_dyadic <- 7 
  refine_grid_size <- 2^J_dyadic
  refine_grid_points <- seq(full_range[1], full_range[2], length.out = refine_grid_size)
  spline_output_matrix <- eval.bifd(refine_grid_points, refine_grid_points, C_est_bifd)
  
  # 5.5b. Prepare input data for GAM
  gam_input_df <- data.frame(
    s = rep(refine_grid_points, times = refine_grid_size),
    t = rep(refine_grid_points, each = refine_grid_size),
    v = as.vector(spline_output_matrix)
  )
  
  # 5.5c. (Core improvement) Dynamically calculate basis dimension k for GAM
  k_val_gam <- floor(gam_k_base * (num_curves / 10)^(gam_k_power))
  k_val_gam <- min(k_val_gam, gam_k_max)
  k_val_gam <- max(k_val_gam, 10)
  cat(sprintf("--> Step 5.5: Dynamically calculated GAM refinement basis dimension k = %d\n", k_val_gam))
  
  # 5.5d. Run lightweight GAM for final smoothing
  final_smoother <- mgcv::gam(v ~ te(s, t, k = k_val_gam), data = gam_input_df, method = "REML")
  
  # --- Step 6: Calculate derivatives on the GAM-refined surface ---
  cat("--> Step 6: Calculating derivatives on the final refined surface...\n")
  grid_for_pred <- expand.grid(s = eval_grid, t = eval_grid)
  final_smooth_vector <- predict(final_smoother, newdata = grid_for_pred)
  final_cov_matrix <- matrix(final_smooth_vector, nrow = length(eval_grid), ncol = length(eval_grid))
  final_cov_matrix <- (final_cov_matrix + t(final_cov_matrix)) / 2
  h <- eval_grid[2] - eval_grid[1]
  grad <- pracma::gradient(final_cov_matrix, h1 = h, h2 = h)
  final_deriv_matrix <- grad[[1]]
  
  cat("--> Covariance partial derivative estimation complete.\n")
  return(final_deriv_matrix)
}
