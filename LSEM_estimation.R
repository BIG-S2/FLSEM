LSEM_estimation <- function(Y, X, Z, K_gram, beta_true=NA, gamma_true=NA, divide_x=FALSE){
  
#Y: the outcome
#X: the functional variable
#Z: the scalar outcome
#beta_true: the true functional coefficient, to calculate MSE in simulation
#gamma_true: the true scalar coefficient, to calculate MSE in simulation
#divide_x: use divide and conquer in the exposure model estimation
  
  ##################
  
  pn_all <- dim(Z)[2]
  n <- dim(Z)[1]
  lent <- dim(X)[2]
  

  sample_for_x <- sample(n, 50)
  
  ###################### variable selection for the functional variable #####
  if(divide_x==FALSE){
      est_X <- functional_GroupL0_gradient( X[sample_for_x, ],Z[sample_for_x, ], K_gram, lambda_range, sp_range, max_iter, nbasis=min(30,lent))
      C_fun_all <- t(do.call(cbind, est_X$final_results$C_function_est))
    }else{
      est_X <- functional_divide_conquer(X, Z, locations_1, locations_2, locations, lambda_range, sp_range, 
        n_split_num=2, overlap_length=10, num_core=30, max_iter=50, nbasis=30)
      C_fun_all <- est_X$C_fun_all
    }
  
  
  ind_for_X <- est_X$idx
  
  Z_for_X <- Z[, which(ind_for_X==TRUE)]

  
  ##########################################################################
  C_sig_est <- crossprod( solve( crossprod(Z_for_X) + 1e-4* diag( dim(Z_for_X)[2] ), diag( dim(Z_for_X)[2] )), crossprod(Z_for_X, X) )
  C_fun_all <- matrix(0, nrow=pn, ncol=lent)
  C_fun_all[which(ind_for_X==TRUE), ] <- C_sig_est
  
  X_fitted <- Z_for_X%*%C_sig_est  #for results_linear
  
  Sigma_mat <- X_fitted%*%K_gram%*%t(X_fitted)/(lent^2)
  Sigma_mat_raw <- X%*%K_gram%*%t(X)/(lent^2)
  
  final_est <- PFLR_RKHS_estimate(Y, X_fitted, Z, K_gram, Sigma_mat,locations, lambda_range=1e-2, sp_range, init_gamma, max_iter)
  est_raw <- PFLR_RKHS_estimate(Y, X, Z, K_gram, Sigma_mat_raw,locations, lambda_range=1e-2, sp_range, init_gamma, max_iter)
  
  
  final_est$MSE_beta_raw <- mean( (est_raw$final_results$beta_est - beta_true)^2 )  
  final_est$MSE_gamma_raw <- sum( (as.numeric(est_raw$final_results$gamma_est) - gamma_true)^2 )  
  
  final_est$MSE_gamma <- sum( ( as.numeric(final_est$final_results$gamma_est) - gamma_true)^2 )
  final_est$MSE_beta <- mean( (final_est$final_results$beta_est - beta_true)^2 )   
  
  final_est$ind_for_X <- ind_for_X
  final_est$ind_for_Y  <- final_est$final_results$idx
  final_est$ind_for_Y_raw  <- est_raw$final_results$idx
  final_est$C_fun_all <- C_fun_all
  
  return(final_est)
  
}