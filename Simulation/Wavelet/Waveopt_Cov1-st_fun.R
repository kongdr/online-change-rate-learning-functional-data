# library(mgcv)
# library(dplyr)
# library(wavethresh)
# library(pracma)

estimate_cov_deriv_WAVELET_GAM_V25_1 <- function(
    data_list,
    eval_grid,
    J_dyadic = 8,
    filter.number = 6,
    family = "DaubLeAsymm",
    gam_k_base = 15,
    gam_k_power = 0.1,
    gam_k_max = 25,
    gam_method = "REML", 
    gam_downsample_ratio = 0.2,
    wavelet_threshold_type = "hard",
    wavelet_threshold_policy = "universal",
    max_pairs_per_curve = 1000,
    
    verbose = TRUE
) {
  all_t <- unlist(lapply(data_list$t, function(x) x[!is.na(x)]))
  all_y <- unlist(lapply(data_list$y, function(x) x[!is.na(x)]))
  
  if (length(all_t) != length(all_y)) {
    stop("Lengths of all_t and all_y do not match after NA removal.")
  }
  if (length(unique(all_t)) < 4) { ### <<< robustness
    stop("Not enough unique time points (< 4) in data_list to estimate mean function.")
  }
  
  safe_mean_smoother <- function(x, y){
    mean_obj <- NULL
    tryCatch({
      mean_obj <- smooth.spline(x = x, y = y, tol = 1e-6, df = min(length(unique(x)) - 1, 7))
    }, error = function(e){
      message("smooth.spline failed: ", e$message, ". Falling back to gam for mean estimation.")
      df_gam <- data.frame(x=x, y=y)
      k_gam_mean <- min(length(unique(x)) -1, 7)
      if(k_gam_mean < 3) {
        warning("Not enough unique x values for gam mean smoothing. Returning constant mean.")
        return(function(t) rep(mean(y, na.rm=TRUE), length(t)))
      }
      mean_obj <- mgcv::gam(y ~ s(x, k = k_gam_mean), data = df_gam, method = "REML")
    })
    
    if ("smooth.spline" %in% class(mean_obj)) {
      mu_hat_func_inner <- function(t) predict(mean_obj, x = t)$y
    } else if ("gam" %in% class(mean_obj)) {
      mu_hat_func_inner <- function(t) predict(mean_obj, newdata = data.frame(x = t))
    } else {
      stop("Failed to estimate mean function.")
    }
    return(mu_hat_func_inner)
  }
  
  mu_hat_func <- safe_mean_smoother(x = all_t, y = all_y)
  
  num_curves <- length(data_list$t)
  if (num_curves == 0) stop("data_list is empty.")
  
  cov_list <- lapply(1:num_curves, function(i) {
    curve_t_orig <- data_list$t[[i]]
    curve_y_orig <- data_list$y[[i]]
    
    valid_indices_curve <- !is.na(curve_t_orig) & !is.na(curve_y_orig)
    curve_t <- curve_t_orig[valid_indices_curve]
    curve_y <- curve_y_orig[valid_indices_curve]
    
    n_obs <- length(curve_t)
    if (n_obs < 2) return(NULL)
    
    residuals <- curve_y - mu_hat_func(curve_t)
    if (any(!is.finite(residuals))) {
      warning(paste("Curve", i, "has non-finite residuals. Skipping."))
      return(NULL)
    }
    
    total_possible_pairs <- n_obs * (n_obs - 1)
    
    if (total_possible_pairs > max_pairs_per_curve) {
      all_idx_pairs <- expand.grid(s_idx = 1:n_obs, t_idx = 1:n_obs)
      all_idx_pairs <- all_idx_pairs[all_idx_pairs$s_idx != all_idx_pairs$t_idx, ]
      
      num_sample_pairs <- min(max_pairs_per_curve, nrow(all_idx_pairs))
      sample_rows <- sample(1:nrow(all_idx_pairs), num_sample_pairs, replace = FALSE)
      
      s_idx_sampled <- all_idx_pairs$s_idx[sample_rows]
      t_idx_sampled <- all_idx_pairs$t_idx[sample_rows]
    } else {
      s_idx_sampled <- rep(1:n_obs, each = n_obs)
      t_idx_sampled <- rep(1:n_obs, times = n_obs)
      valid_pairs <- s_idx_sampled != t_idx_sampled
      s_idx_sampled <- s_idx_sampled[valid_pairs]
      t_idx_sampled <- t_idx_sampled[valid_pairs]
    }
    
    s_vals <- curve_t[s_idx_sampled]
    t_vals <- curve_t[t_idx_sampled]
    v_vals <- residuals[s_idx_sampled] * residuals[t_idx_sampled]
    
    finite_pairs <- is.finite(s_vals) & is.finite(t_vals) & is.finite(v_vals)
    if(sum(finite_pairs) == 0) return(NULL)
    
    data.frame(s = s_vals[finite_pairs], t = t_vals[finite_pairs], v = v_vals[finite_pairs])
  })
  
  cov_data_df <- do.call(rbind, cov_list)
  if (is.null(cov_data_df) || nrow(cov_data_df) == 0) {
    stop("No valid covariance pairs generated.")
  }
  
  if (requireNamespace("data.table", quietly = TRUE)) {
    cov_dt <- data.table::as.data.table(cov_data_df)
    cov_dt[, s_sorted := pmin(s, t)]
    cov_dt[, t_sorted := pmax(s, t)]
    cov_data_df_unique <- as.data.frame(unique(cov_dt, by = c("s_sorted", "t_sorted"))[, c("s", "t", "v")])
  } else { # use dplyr
    cov_data_df_unique <- cov_data_df %>%
      dplyr::mutate(s_sorted = pmin(s, t), t_sorted = pmax(s, t)) %>%
      dplyr::distinct(s_sorted, t_sorted, .keep_all = TRUE) %>%
      dplyr::select(-s_sorted, -t_sorted)
  }
  grid_size <- 2^J_dyadic
  s_breaks <- seq(min(eval_grid), max(eval_grid), length.out = grid_size + 1)
  t_breaks <- s_breaks
  s_mids <- (s_breaks[-1] + s_breaks[-(grid_size+1)]) / 2
  t_mids <- s_mids
  
  binned_cov <- cov_data_df_unique %>%
    dplyr::mutate(s_bin = .bincode(s, s_breaks, include.lowest = TRUE),
                  t_bin = .bincode(t, t_breaks, include.lowest = TRUE)) %>%
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
  k_val <- min(k_val, gam_k_max)
  k_val <- max(k_val, 5)

  k_val_s <- min(k_val, length(unique(s_mids)) - 1)
  k_val_t <- min(k_val, length(unique(t_mids)) - 1)
  k_val_s <- max(k_val_s, 3) # k at least 3
  k_val_t <- max(k_val_t, 3) 
  
  gam_input_df <- data.frame(s = rep(s_mids, times = grid_size),
                             t = rep(t_mids, each = grid_size),
                             v = as.vector(cov_wavelet_denoised_matrix))
  gam_input_df <- na.omit(gam_input_df)
  if (nrow(gam_input_df) == 0) stop("GAM input data is empty after NA removal.")
  
  if (gam_downsample_ratio < 1 && nrow(gam_input_df) > 500) {
    sample_size <- max(500, floor(nrow(gam_input_df) * gam_downsample_ratio))
    gam_input_df <- gam_input_df[sample(1:nrow(gam_input_df), sample_size, replace = FALSE), ]
    if (verbose) cat(paste("--> step 3", nrow(gam_input_df), "row.\n"))
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
