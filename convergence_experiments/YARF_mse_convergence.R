

#' 
#' @param YARF_oob_matrix       A matrix of all out-of-bag estimates
#' @param num_trees             Number of trees to include in MSE calculation
#' 
YARF_oob_mse = function(YARF_oob_matrix, num_trees = NULL) {
  if (is.null(num_trees)){
    num_trees = ncol(YARF_oob_matrix)
  }
  
  mse = array(NA, num_trees)
  for (j in 1:num_trees){
    mse[j] = mean((y - rowMeans(y_oob_matrix[, 1:j, drop = FALSE], na.rm = TRUE))^2, na.rm = TRUE)
  }
  return(mse)
}



#' 
#' @param YARF_mod        A YARF model object
#' @param min_trees       The minimum number of trees to be sampled for MSE calculation
#' @param max_trees       The maximum number of trees to be sampled for MSE calculation
#' 
#' @param delta           Interval spacing for sampling        
#' @param n_R             Number of samples for each number of trees
#' 
YARF_oob_mse_table = function(YARF_mod, min_trees = 100, max_trees = 1000, delta = 25, n_R = 20) {
  num_trees = YARF_mod$num_trees
  if (max_trees >= num_trees) stop(paste("max_trees must be less than the number of trees used in YARF_mod: ", num_trees))
  
  num_ms = floor((max_trees - min_trees) / delta)
  m_grid = min_trees + (delta * (0 : num_ms))
  
  y_oob_matrix = YARF::YARF_all_oob_results_matrix(YARF_mod)
  
  regression_df = matrix(NA, nrow = length(m_grid) * n_R, ncol = 2)
  colnames(regression_df) = c("inv_num_trees", "MSE_estimate")
  
  for (j in 1:length(m_grid)) {
    for (i in 1:n_R) {
      sample_indices = sample(1:ncol(y_oob_matrix), m_grid[j], replace = FALSE)
      oob = y_oob_matrix[, sample_indices]
      regression_df[(j-1) * n_R + i, ] = c(1 / m_grid[j], mean((YARF_mod$y - rowMeans(oob, na.rm = TRUE))^2, na.rm = TRUE))
    }
  }
  return(data.frame(regression_df))
}


#' 
#' @param yarf_mod      A YARF model object
#' @param mse_df        An MSE table for regression
#' @param beta          Quantile
#' @param alpha         
#'
YARF_optim_num_trees = function(yarf_mod = NULL, mse_df = NULL, beta = 0.99, alpha = 0.05){
  if (is.null(mse_df)){
    mse_df = YARF_oob_mse_table(yarf_mod)
  }
  
  mse_ols_mod = lm(MSE_estimate ~ inv_num_trees, data = mse_df)
  a_hat = as.numeric(mse_ols_mod$coefficients[1])
  b_hat = as.numeric(mse_ols_mod$coefficients[2])
  ols_rmse = summary(mse_ols_mod)$sigma
  sum_ms = sum((mse_df[, 1] - mean(mse_df[, 1]))^2)
  
  for (i in 1:nrow(mse_df)){
    sigma_hat = ols_rmse * 
      sqrt(1 + (1 / nrow(mse_df)) + ((mse_df[i, 1] - mean(mse_df[, 1]))^2 / sum_ms))
    lower_bound = b_hat / ((alpha * a_hat) - (sigma_hat * qnorm(beta)))
    if (mse_df[i, 1]^-1 >= lower_bound){
      m_star = ceiling(lower_bound)
      break
    }
  }
  return(list(m_star = m_star, a_hat = a_hat, b_hat = b_hat))
}





# library(dplyr)
# mse_df = YARF_oob_MSE_table(yarf_mod)
# mse_df = mse_df %>% 
#   mutate(num_trees = inv_num_trees^-1)
# ggplot(mse_df) +
#   geom_point(aes(x = num_trees, y = MSE_estimate)) +
#   scale_x_log10() +
#   geom_hline(yintercept = 11.22, col = "blue") +
#   geom_hline(yintercept = 10.688, col = "red")













