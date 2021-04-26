

mse_array = function(y, oob_matrix, m_trees_array){
  
  mse = array(NA, length(m_trees_array))
  
  for (j in 1:length(m_trees_array)){
    mse[j] = mean((y - rowMeans(oob_matrix[, 1:m_trees_array[j], drop = FALSE], na.rm = TRUE))^2, na.rm = TRUE)
  }
  
  return(mse)
}




#' Locate the optimal number of trees for a YARF model to converge
#' 
#' @param yarf_mod              The yarf model object.
#' @param num_sample_trees      Number of trees to sample for each MSE calculation. If \code{yarf_mod$wait == TRUE}, 
#'                              this parameter is ignored.
#' @param m_0                   Minimum number of trees to begin calculation MSE from.
#' @param delta_m               Interval of trees between MSE calculations.
#' @param n_sim                 Number of simulations for re-sampling.
#' @param alpha_star            Closeness parameter for convergence.
#' 
#' 
locate_optimal_num_trees = function(yarf_mod, num_sample_trees = NULL, m_0 = 100, delta_m = 20, n_sim = 100, alpha_star = 0.01){
  assertClass(yarf_mod, "YARF")
  
  
  if (yarf_mod$wait == TRUE){
    
    if (is.null(num_sample_trees)) {
      num_sample_trees = 1000
    }
    
    if (num_sample_trees >= yarf_mod$num_trees) {
      warning("Not enough trees to sample \nUsing default num_sample_trees = num_trees - 100\n")
      num_sample_trees = yarf_mod$num_trees - 100
    }
    
    trees_to_sample = 1 : num_sample_trees
    num_trees = seq(m_0, num_sample_trees, by = delta_m)
    
    mse_matrix = matrix(NA, nrow = length(num_trees), ncol = n_sim)
    
    # get all oob results for mse calculation
    oob_matrix = YARF_all_oob_results_matrix(yarf_mod)
    
    cat("Model has total trees =", yarf_mod$num_trees, "\nBeginning m_star calculation...\n")
    
    # resampling
    for (i in 1:n_sim){
      
      indices = sample(trees_to_sample, num_sample_trees, replace = FALSE)
      oob_sample = oob_matrix[, indices]
      
      mse_matrix[, i] = mse_array(yarf_mod$y, oob_sample, num_trees)
      
    }
    
    mse = rowMeans(mse_matrix)
    
    # model the mse to find m_star
    mse_mod = lm(mse ~ I(1 / num_trees))
    coeff = coef(mse_mod)
    
    alpha_m = array(NA, length(num_trees))
    
    for (i in 1:length(num_trees)){
      
      alpha_m[i] = coeff[2] / (num_trees[i] * coeff[1])
      Sys.sleep(.05)
      
      # if you are not within alpha_star of convergence, continue letting the user know
      if (alpha_m[i] > alpha_star) {
        cat("Testing m =", paste0(num_trees[i], ","), "MSE =", paste0(round(mse[i], 5), ","), "alpha_m =", alpha_m[i], "\n")
        m_star = i
      }
    }
    
    if (is.na(num_trees[m_star + 1])){
      stop("Failed to converge with specified parameters: alpha_star = ", alpha_star)
    }
    
    # point of convergence
    m_star = num_trees[m_star + 1]
    cnvg_line = mse[which(num_trees == m_star)]
    
    # convergence plot
    plot = ggplot(data.frame(mse, num_trees)) +
      geom_point(aes(num_trees, mse)) +
      geom_smooth(aes(num_trees, mse), formula = y ~ I(1/x), method = "lm", se = FALSE) +
      geom_hline(aes(yintercept = cnvg_line), linetype = "dashed") +
      geom_point(aes(x = m_star, y = cnvg_line, col = "pink")) +
      theme(legend.position = "none")
    print(plot)
    
    # save results
    yarf_mod$MSE_plot <<- plot
    yarf_mod$MSEs <<- data.frame(cbind(m = num_trees, mse = mse))
    yarf_mod$optimal_num_trees_raw_data <<- data.frame(cbind(m = num_trees, mse = mse, alpha = alpha_m))
    yarf_mod$m_star <<- list(m_star = m_star, a = coeff[[1]], b = coeff[[2]])
    
    cat("\nm_star ~", m_star)
    
  }
  
  else {
    
    num_trees_completed = .jcall(yarf_mod$java_YARF, "I", "progress")
    
    while (num_trees_completed < 750) {
      Sys.sleep(1)
      num_trees_completed = .jcall(yarf_mod$java_YARF, "I", "progress")
    }
    
    # starting parameters
    oob_matrix = YARF_all_oob_results_matrix(yarf_mod)
    num_sample_trees = dim(oob_matrix)[2] - 100
    
    trees_to_sample = 1 : num_sample_trees
    num_trees = seq(m_0, num_sample_trees, by = delta_m)
    
    mse_matrix = matrix(NA, nrow = length(num_trees), ncol = n_sim)
    
    while (TRUE){
      
      cat(num_trees_completed, "trees completed. \nCalculating m_star...\n")
      
      # resampling
      for (i in 1:n_sim){
        
        indices = sample(trees_to_sample, num_sample_trees, replace = FALSE)
        oob_sample = oob_matrix[, indices]
        
        mse_matrix[, i] = mse_array(yarf_mod$y, oob_sample, num_trees)
        
      }
      
      mse = rowMeans(mse_matrix)
      
      # model the mse to find m_star
      mse_mod = lm(mse ~ I(1 / num_trees))
      coeff = coef(mse_mod)
      
      alpha_m = array(NA, length(num_trees))
      
      for (i in 1:length(num_trees)){
        
        alpha_m[i] = coeff[2] / (num_trees[i] * coeff[1])
        
        # if you are not within alpha_star of convergence, continue letting the user know
        if (alpha_m[i] > alpha_star) {
          cat("Testing m =", paste0(num_trees[i], ","), "MSE =", paste0(round(mse[i], 5), ","), "alpha_m =", alpha_m[i], "\n")
          m_star = i
        }
      }
      
      # if m_star was found record it
      if (!is.na(num_trees[m_star + 1])){
        
        m_star = num_trees[m_star + 1]
        cnvg_line = mse[which(num_trees == m_star)]
        
        # convergence plot
        plot = ggplot(data.frame(mse, num_trees)) +
          geom_point(aes(num_trees, mse)) +
          geom_smooth(aes(num_trees, mse), formula = y ~ I(1/x), method = "lm", se = FALSE) +
          geom_hline(aes(yintercept = cnvg_line), linetype = "dashed") +
          geom_point(aes(x = m_star, y = cnvg_line, col = "pink")) +
          theme(legend.position = "none")
        print(plot)
        
        # save results
        yarf_mod$MSE_plot <<- plot
        yarf_mod$MSEs <<- data.frame(cbind(m = num_trees, mse = mse))
        yarf_mod$optimal_num_trees_raw_data <<- data.frame(cbind(m = num_trees, mse = mse, alpha = alpha_m))
        yarf_mod$m_star <<- list(m_star = m_star, a = coeff[[1]], b = coeff[[2]])
        
        cat("\nm_star ~", m_star, "\n\n")
        
        break
      }
      
      else {
        
        # add trees
        num_trees_completed = .jcall(yarf_mod$java_YARF, "I", "progress")
        
        oob_matrix = YARF_all_oob_results_matrix(yarf_mod)
        num_sample_trees = dim(oob_matrix)[2] 
        
        trees_to_sample = 1 : num_sample_trees
        num_trees = c(num_trees, seq(num_trees[length(num_trees)] + delta_m, (num_sample_trees), by = delta_m))
        
        mse_matrix = matrix(NA, nrow = length(num_trees), ncol = n_sim)
        
      }
      
      
    }
    
    YARF_stop(yarf_mod = yarf_mod)
    
    num_trees_completed = .jcall(yarf_mod$java_YARF, "I", "progress")
    cat("YARF model stopped at", num_trees_completed, "trees\n")
    
  }
  
}





#############################################################

options(java.parameters = "-Xmx10000m")
pacman::p_load(YARF)
set_YARF_num_cores(7)

library(MASS); data(Boston)
X = Boston[, 1 : 13]; y = Boston[, 14]

yarf_mod = YARF(X, y, num_trees = 5000, wait = F)
locate_optimal_num_trees(yarf_mod, alpha_star = .0075)

yarf_mod_wait = YARF(X, y, num_trees = 1000)
locate_optimal_num_trees(yarf_mod_wait)


