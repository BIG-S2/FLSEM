library(Matrix)
functional_GroupL0_gradient <- function(z, x, K_gram, lambda_range, sp_range, max_iter=200,nbasis=NULL, position_fixed=NULL){
  
  nobs <- dim(z)[1]
  mobs <- dim(z)[2]
  nmobs <- nobs*mobs
  
  pn<- dim(x)[2]
  
  vz <- as.vector(t(z))
  
  
  
  id.basis <- sample(c(1:mobs),nbasis,replace=FALSE)
  id.basis <- sort(id.basis)
  ############## generate final id.basis ######
  K_gram_select <- K_gram[, id.basis ]
  pen_mat <- K_gram[id.basis, id.basis ]

  pen_mat_all <-as.matrix( bdiag(replicate(pn,pen_mat,simplify=FALSE))  )
  
  xmat_all_matrix <- NULL
  
  for(p in 1:pn){
    xmat_all_matrix <- cbind(xmat_all_matrix,  kronecker(x[, p], K_gram_select) )
  }
  
  sp_length <- length(sp_range)
  lambda_length <- length(lambda_range)
  
  init_C<- rep( list(rep(0, nbasis)), pn )
  
  results.sp <- NULL
  HBIC.sp <- NULL 
  
  for(i in 1:sp_length){
    sp_tmp <- sp_range[i]
    
    for(j in 1:lambda_length){
      
      lambda_tmp <- lambda_range[j]
      
      
      active_idx <- NULL
      iter <- 0
      
      update_C <- init_C
      
      while(iter < max_iter){
        #print(iter)
        iter=iter+1
        pre_active_idx<- active_idx

        update_C_vector <- do.call(c, update_C)

        gradient <- -crossprod( xmat_all_matrix, vz - crossprod( t(xmat_all_matrix), update_C_vector ) )/(nmobs) + lambda_tmp* crossprod(pen_mat_all, update_C_vector )

        eta <- 0.03   
        update_C_vector_next  <- update_C_vector  - eta *  gradient 
        update_C_next <- lapply(seq_len(pn), function(x)update_C_vector_next[((x-1)*nbasis+1):(x*nbasis)]   )
        rm(list=c("update_C_vector", "update_C_vector_next"))
        gc()
        
        arr <- as.numeric( lapply(update_C_next , function(x){sum(x^2)}) )
        thre <- sort(arr, decreasing = T)[sp_tmp]
        
        
        active_idx <-  unique( c( which( arr >= thre ),position_fixed) )
        inactive_idx <- c(1:pn)[-active_idx]
        
        idx <- logical(pn)
        idx[active_idx] <- TRUE

        update_C_next[inactive_idx ]=rep(list(rep(0, nbasis)), length(inactive_idx) )

        C_diff <- mean( unlist( lapply(seq_len(pn),function(x){sum( (update_C_next[[x]]- update_C[[x]])^2) }  )   )  )
        
        if( C_diff < 1e-3  ){break}


        update_C  <- update_C_next 
        rm(update_C_next)
        gc()
      }

     
     res_vz <- vz - crossprod( t(xmat_all_matrix), do.call(c, update_C) )

      
    HBIC.tmp <-  log(sum(res_vz^2) ) +  8*sp_tmp*nbasis*log(log(nmobs))*log(pn)/(nmobs)
    
    results.tmp <- list(  C_est= update_C,  idx=idx, sp=sp_tmp,
                          lambda =lambda_tmp, iter=iter, hbic=HBIC.tmp )  #fitted_vz=fitted_vz, 
    results.sp <- c(results.sp, list(results.tmp))
    HBIC.sp <- c(HBIC.sp, HBIC.tmp)

    }
    
    
  }
  
  sp_lambda_position <- which( HBIC.sp==min(HBIC.sp))[1]
  

  final_results <- results.sp[[sp_lambda_position]]
  
  final_results$C_function_est <- lapply(final_results$C_est, function(x) K_gram_select%*%x )

  fitted_value <- x%*% t(do.call(cbind, final_results$C_function_est))
  
  list( final_results=final_results, results.sp.all=results.sp, fitted_value=fitted_value, idx=final_results$idx )
  
}