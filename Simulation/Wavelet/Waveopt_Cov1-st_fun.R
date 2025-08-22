
# install.packages(c("wavethresh", "mgcv", "pracma", "dplyr"))

estimate_cov_deriv_WAVELET_GAM_V25_1 <- function(
    data_list, 
    eval_grid,
    J_dyadic = 8,  
    filter.number = 6,
    family = "DaubLeAsymm",
    gam_k_base = 15,    
    gam_k_power = 0.1, 
    gam_k_max = 40,     
    verbose = TRUE
) {
  
  all_t <- unlist(data_list$t); all_y <- unlist(data_list$y)
  mean_smoother <- smooth.spline(x = all_t, y = all_y)
  mu_hat_func <- function(t) predict(mean_smoother, x = t)$y
  
  num_curves <- length(data_list$t)
  cov_list <- lapply(1:num_curves, function(i) {
    curve_t <- data_list$t[[i]]; curve_y <- data_list$y[[i]]; n_obs <- length(curve_t)
    if (n_obs < 2) return(NULL)
    residuals <- curve_y - mu_hat_func(curve_t)
    idx_pairs <- t(combn(1:n_obs, 2)); s <- curve_t[idx_pairs[, 1]]; t <- curve_t[idx_pairs[, 2]]
    v <- residuals[idx_pairs[, 1]] * residuals[idx_pairs[, 2]]
    return(data.frame(s = c(s, t), t = c(t, s), v = c(v, v)))
  })
  cov_data_df <- do.call(rbind, cov_list)
  if (verbose) cat(paste("--> step0: nrow(cov_data_df), "raw covariance (from", num_curves, "curves).\n"))
  
  grid_size <- 2^J_dyadic
  if (verbose) cat(paste("--> step1: start", grid_size, "x", grid_size, "binning...\n"))
  s_breaks <- seq(min(cov_data_df$s, na.rm=T), max(cov_data_df$s, na.rm=T), length.out = grid_size + 1)
  t_breaks <- seq(min(cov_data_df$t, na.rm=T), max(cov_data_df$t, na.rm=T), length.out = grid_size + 1)
  s_mids <- (s_breaks[-1] + s_breaks[-(grid_size+1)])/2
  t_mids <- (t_breaks[-1] + t_breaks[-(grid_size+1)])/2
  binned_cov <- cov_data_df %>% dplyr::mutate(s_bin = .bincode(s, s_breaks, include.lowest = TRUE), t_bin = .bincode(t, t_breaks, include.lowest = TRUE)) %>% dplyr::filter(!is.na(s_bin) & !is.na(t_bin)) %>% dplyr::group_by(s_bin, t_bin) %>% dplyr::summarise(v_mean = mean(v, na.rm = TRUE))
  raw_cov_matrix <- matrix(0, nrow = grid_size, ncol = grid_size)
  raw_cov_matrix[cbind(binned_cov$s_bin, binned_cov$t_bin)] <- binned_cov$v_mean
  raw_cov_matrix <- (raw_cov_matrix + t(raw_cov_matrix)) / 2
  if (verbose) cat("--> step1: binning finished. \n")
  
  wt <- wavethresh::imwd(raw_cov_matrix, filter.number = filter.number, family = family)
  wt_thr <- wavethresh::threshold(wt, type = "hard", policy = "universal")
  cov_wavelet_denoised_matrix <- wavethresh::imwr(wt_thr)
  
  k_val <- floor(gam_k_base * (num_curves / 10)^(gam_k_power))
  k_val <- min(k_val, gam_k_max)
  k_val <- max(k_val, 10) 
  
  gam_input_df <- data.frame(s = rep(s_mids, times = grid_size), t = rep(t_mids, each = grid_size), v = as.vector(cov_wavelet_denoised_matrix))
  
  final_smoother <- mgcv::gam(v ~ te(s, t, k = k_val), data = gam_input_df, method = "REML")
  
  grid_for_pred <- expand.grid(s = eval_grid, t = eval_grid)
  final_smooth_vector <- predict(final_smoother, newdata = grid_for_pred)
  final_cov_matrix <- matrix(final_smooth_vector, nrow = length(eval_grid), ncol = length(eval_grid))
  final_cov_matrix <- (final_cov_matrix + t(final_cov_matrix)) / 2

  h <- eval_grid[2] - eval_grid[1]
  grad <- pracma::gradient(final_cov_matrix, h1 = h, h2 = h)
  final_deriv_matrix <- grad[[1]]
  
  return(final_deriv_matrix)
}
