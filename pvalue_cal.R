
pvalue_cal <- function(Y, X,Z_for_Y, K_gram){
 
  projection_P_U <- diag(n) - Z_for_Y%*%tcrossprod( solve( crossprod(Z_for_Y), diag(dim(Z_for_Y)[2]) ), Z_for_Y )
    
  y_new <- projection_P_U%*%Y



  sigma_est <-mean(y_new^2)

  lent <- dim(X)[2]

  Sigma_mat <- X%*%K_gram%*%t(X)/(lent^2)
  
  S_n  <- t(y_new)%*%Sigma_mat%*%y_new / (n* sigma_est)
  
  Rn <- projection_P_U%*%X%*%K_gram%*%t(X)%*%projection_P_U/(n*lent^2)
  Rn_square <- Rn %*% Rn
  
  tr_Rn <- sum(diag(Rn))
  tr_Rn_square <-sum(diag(Rn_square))
  
  kappa <- tr_Rn_square / tr_Rn
  zeta <- tr_Rn^2 / tr_Rn_square 
  
  S_n_normed <- S_n / kappa
  
  p_value <- 1- pchisq(S_n_normed, df= zeta )
  
  test_results <- c( S_n,  kappa, zeta, S_n_normed, p_value )

  names(test_results ) <- c('statistic', 'para_kappa', 'para_zeta', 'statistic_standarized', 'pvalue')
  
  return(test_results)
}




