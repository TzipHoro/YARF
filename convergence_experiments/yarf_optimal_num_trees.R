#' YARF Convergence
#' 
#' @param yarf_mod            A YARF model object.
#' @param num_sample_trees    Number of trees in each resampling.
#' @param mse_start           Minimum number of trees to start calculating MSE.
#' @param mse_by              
#' @param Nsim                Number of times to resample.
#' @param alpha               
#' 
#' 
yarf_optimal_num_trees = function(yarf_mod, num_sample_trees = 2000, mse_start = 100, mse_by = 25, Nsim = 100, alpha = 0.01) {
  assertClass(yarf_mod, "YARF")
  
  if (yarf_mod$num_trees <= num_sample_trees) {
    stop("num_trees must be greater than num_sample_trees")
  } else {
    message("Calculating OOB results...")
    oob_matrix = YARF_all_oob_results_matrix(yarf_mod)
  }
  
  # find all mse's for each sampling
  trees_to_sample = 1 : num_sample_trees
  num_trees = seq(mse_start, num_sample_trees, by = mse_by)
  regression_results = matrix(NA, nrow = Nsim, ncol = 2)
  
  message("Beginning mstar calculation...")
  for (n in 1:Nsim){
    indices = sample(trees_to_sample, num_sample_trees, replace = FALSE)
    oob_sample = oob_matrix[, indices]
    
    mse = array(NA, length(num_trees))
    for (j in 1:length(num_trees)){
      mse[j] = mean((yarf_mod$y - rowMeans(oob_sample[, 1:num_trees[j], drop = FALSE], na.rm = TRUE))^2, na.rm = TRUE)
    }
    
    regress = lm(mse ~ I(1 / num_trees))
    regression_results[n, ] = coef(regress)
    
    if (n %% 10 == 0) {
      message("   Finding mstar: ", round((n/Nsim) * 100), "% complete")
    }
  }
  
  mstars = regression_results[, 2] / (alpha * regression_results[, 1])
  round(mean(mstars))
}



options(java.parameters = "-Xmx5g")
library(YARF)
library(MASS); data(Boston)
X = Boston[, 1 : 13]; y = Boston[, 14]

# with 10000 trees
yarf_mod = YARF(X, y, num_trees = 10000)
m_star_50sims = yarf_optimal_num_trees(yarf_mod = yarf_mod, Nsim = 50)
m_star_100sims = yarf_optimal_num_trees(yarf_mod = yarf_mod)
m_star_200sims = yarf_optimal_num_trees(yarf_mod = yarf_mod, Nsim = 200)
m_star_300sims = yarf_optimal_num_trees(yarf_mod = yarf_mod, Nsim = 300)

# with fewer trees
yarf_mod1000 = YARF(X, y, num_trees = 1000)
m_star_1000_100sims = yarf_optimal_num_trees(yarf_mod1000, num_sample_trees = 500)
m_star_1000_200sims = yarf_optimal_num_trees(yarf_mod1000, num_sample_trees = 500, Nsim = 200)

# with num_trees close to num_sample_trees
m_star_1000_950trees = yarf_optimal_num_trees(yarf_mod1000, num_sample_trees = 950, mse_start = 5)
m_star_1000_990trees = yarf_optimal_num_trees(yarf_mod1000, num_sample_trees = 990)






