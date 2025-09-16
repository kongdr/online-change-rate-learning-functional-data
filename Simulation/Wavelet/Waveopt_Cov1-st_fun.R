# install.packages(c("mgcv", "dplyr", "wavethresh", "pracma", "data.table"))

estimate_cov_deriv_WAVELET_GAM_V25_1 <- function(
    data_list,
    eval_grid,
    J_dyadic = 10,
    gam_max_points = 25000,
    max_pairs_per_curve = 3000,
    gam_k_base = 15,
    gam_k_max = 25,
    mean_J_dyadic = 4,
    gam_k_power=0.1,
    filter.number = 6,
    family = "DaubLeAsymm",
    gam_method = "REML", 
    wavelet_threshold_type = "hard",
    wavelet_threshold_policy = "universal",
    verbose = TRUE
) {
  
  all_t <- unlist(lapply(data_list$t, function(x) x[!is.na(x)]))
  all_y <- unlist(lapply(data_list$y, function(x) x[!is.na(x)]))
  
  if (length(all_t) != length(all_y)) stop("Lengths of all_t and all_y do not match after NA removal.")
  if (length(unique(all_t)) < 4) stop("Not enough unique time points (< 4) to estimate mean function.")
  
  wavelet_mean_smoother <- function(x, y, J, filt.num, fam, thresh_type, thresh_policy) {
    grid_len <- 2^J
    t_range <- range(x, na.rm = TRUE)
    t_breaks <- seq(t_range[1], t_range[2], length.out = grid_len + 1)
    t_mids <- (t_breaks[-1] + t_breaks[-(grid_len + 1)]) / 2
    
    bins <- .bincode(x, t_breaks, include.lowest = TRUE)
    y_binned_means <- tapply(y, bins, mean, na.rm = TRUE)
    
    y_grid <- rep(NA, grid_len) 
    valid_bins <- as.integer(names(y_binned_means))
    y_grid[valid_bins] <- y_binned_means
    
    is_na_grid <- is.na(y_grid)
    if(any(is_na_grid)) {
      non_na_indices <- which(!is_na_grid)
      if (length(non_na_indices) > 1) {
        interp_result <- stats::approx(x = non_na_indices, y = y_grid[non_na_indices], xout = 1:grid_len, method = "linear", rule = 2)
        y_grid <- interp_result$y
      } else if (length(non_na_indices) == 1) { 
        y_grid <- rep(y_grid[non_na_indices], grid_len)
      } else {
        y_grid <- rep(mean(y, na.rm = TRUE), grid_len)
      }
    }
    
    wt_mean <- wavethresh::wd(y_grid, filter.number = filt.num, family = fam)
    wt_mean_thr <- wavethresh::threshold(wt_mean, type = thresh_type, policy = thresh_policy)
    y_smoothed_grid <- wavethresh::wr(wt_mean_thr)
    
    mu_hat_func_inner <- stats::splinefun(t_mids, y_smoothed_grid)
    
    return(mu_hat_func_inner)
  }
  
  mu_hat_func <- wavelet_mean_smoother(
    x = all_t, y = all_y, J = mean_J_dyadic,
    filt.num = filter.number, fam = family,
    thresh_type = wavelet_threshold_type, thresh_policy = wavelet_threshold_policy
  )
  
  num_curves <- length(data_list$t)
  if (num_curves == 0) stop("data_list is empty.")
  
  cov_list <- lapply(1:num_curves, function(i) {
    curve_t_orig <- data_list$t[[i]]; curve_y_orig <- data_list$y[[i]]
    valid_indices <- !is.na(curve_t_orig) & !is.na(curve_y_orig)
    curve_t <- curve_t_orig[valid_indices]; curve_y <- curve_y_orig[valid_indices]
    n_obs <- length(curve_t)
    if (n_obs < 2) return(NULL)
    
    residuals <- curve_y - mu_hat_func(curve_t)
    if (any(!is.finite(residuals))) { warning(paste("Curve", i, "has non-finite residuals. Skipping.")); return(NULL) }
    
    total_possible_pairs <- n_obs * (n_obs - 1) / 2
    if (total_possible_pairs > max_pairs_per_curve) {
      all_idx_pairs <- expand.grid(s_idx = 1:n_obs, t_idx = 1:n_obs)
      all_idx_pairs <- all_idx_pairs[all_idx_pairs$s_idx < all_idx_pairs$t_idx, ]
      sample_rows <- sample(nrow(all_idx_pairs), min(max_pairs_per_curve, nrow(all_idx_pairs)))
      s_idx_sampled <- all_idx_pairs$s_idx[sample_rows]; t_idx_sampled <- all_idx_pairs$t_idx[sample_rows]
    } else {
      pair_indices <- which(lower.tri(matrix(0, n_obs, n_obs)), arr.ind = TRUE)
      s_idx_sampled <- pair_indices[, 2]; t_idx_sampled <- pair_indices[, 1]
    }
    
    s_vals <- curve_t[s_idx_sampled]; t_vals <- curve_t[t_idx_sampled]
    v_vals <- residuals[s_idx_sampled] * residuals[t_idx_sampled]
    
    finite_pairs <- is.finite(s_vals) & is.finite(t_vals) & is.finite(v_vals)
    if(sum(finite_pairs) == 0) return(NULL)
    
    data.frame(s = s_vals[finite_pairs], t = t_vals[finite_pairs], v = v_vals[finite_pairs])
  })
  
  cov_data_df <- data.table::rbindlist(cov_list)
  if (is.null(cov_data_df) || nrow(cov_data_df) == 0) stop("No valid covariance pairs generated.")
  
  grid_size <- 2^J_dyadic
  grid_range <- range(eval_grid, na.rm = TRUE)
  s_breaks <- seq(grid_range[1], grid_range[2], length.out = grid_size + 1)
  t_breaks <- s_breaks
  s_mids <- (s_breaks[-1] + s_breaks[-(grid_size+1)]) / 2; t_mids <- s_mids
  
  binned_cov <- cov_data_df %>%
    dplyr::mutate(
      s_bin = .bincode(s, s_breaks, include.lowest = TRUE),
      t_bin = .bincode(t, t_breaks, include.lowest = TRUE)
    ) %>%
    dplyr::filter(!is.na(s_bin) & !is.na(t_bin)) %>%
    dplyr::group_by(s_bin, t_bin) %>%
    dplyr::summarise(v_mean = mean(v, na.rm = TRUE), .groups = 'drop')
  
  raw_cov_matrix <- matrix(0, nrow = grid_size, ncol = grid_size)
  raw_cov_matrix[as.matrix(binned_cov[, c("s_bin", "t_bin")])] <- binned_cov$v_mean
  raw_cov_matrix <- (raw_cov_matrix + t(raw_cov_matrix)) / 2
  raw_cov_matrix[is.nan(raw_cov_matrix)] <- 0
  
  wt <- wavethresh::imwd(raw_cov_matrix, filter.number = filter.number, family = family)
  wt_thr <- wavethresh::threshold(wt, type = wavelet_threshold_type, policy = wavelet_threshold_policy)
  cov_wavelet_denoised_matrix <- wavethresh::imwr(wt_thr)
  cov_wavelet_denoised_matrix <- (cov_wavelet_denoised_matrix + t(cov_wavelet_denoised_matrix)) / 2
  
  
  k_val <- floor(gam_k_base * (num_curves / 10)^(gam_k_power))
  k_val <- min(k_val, gam_k_max); k_val <- max(k_val, 5)
  k_val_s <- min(k_val, length(unique(s_mids)) - 1); k_val_t <- min(k_val, length(unique(t_mids)) - 1)
  k_val_s <- max(k_val_s, 3); k_val_t <- max(k_val_t, 3) 
  

  gam_input_df <- data.frame(
    s = rep(s_mids, times = grid_size),
    t = rep(t_mids, each = grid_size),
    v = as.vector(cov_wavelet_denoised_matrix)
  )
  gam_input_df <- na.omit(gam_input_df)
  if (nrow(gam_input_df) == 0) stop("GAM input data is empty after NA removal.")
  
  if (nrow(gam_input_df) > gam_max_points) {
    gam_input_df <- gam_input_df[sample(nrow(gam_input_df), gam_max_points, replace = FALSE), ]
  }
  
  final_smoother <- mgcv::gam(v ~ te(s, t, k = c(k_val_s, k_val_t)), data = gam_input_df, method = gam_method)
  
  grid_for_pred <- expand.grid(s = eval_grid, t = eval_grid)
  final_smooth_vector <- predict(final_smoother, newdata = grid_for_pred)
  final_cov_matrix <- matrix(final_smooth_vector, nrow = length(eval_grid), ncol = length(eval_grid))
  final_cov_matrix <- (final_cov_matrix + t(final_cov_matrix)) / 2
  
  h <- eval_grid[2] - eval_grid[1]
  if (is.na(h) || h <= 0) stop("eval_grid must be equally spaced and have length > 1.")
  
  grad <- pracma::gradient(final_cov_matrix, h1 = h, h2 = h)
  final_deriv_matrix <- grad[[1]]
  
  return(final_deriv_matrix)
}
