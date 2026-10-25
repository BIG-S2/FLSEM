library(Matrix)
functional_divide_conquer <- function(z, x, locations_1, locations_2, locations, lambda_range, sp_range, n_split_num, overlap_length, num_core, max_iter=200, nbasis=NULL){
  
  nobs <- dim(z)[1]
  mobs <- dim(z)[2]
  nmobs <- nobs*mobs
  
  pn<- dim(x)[2]
  
  vz <- as.vector(t(z))


  subject_groups<- split(1:nobs, rep(1:n_split_num, each = ceiling(n/n_split_num)))


  # Create a sample matrix of dimension 6x8 (replace with your actual matrix)
matrix_rows <- length(locations_1)
matrix_cols <- length(locations_2)
sample_matrix <- matrix(1:(matrix_rows * matrix_cols), nrow = matrix_rows, ncol = matrix_cols)

# Define the dimensions of the submatrices
submatrix_rows <- 2*overlap_length
submatrix_cols <- 2*overlap_length
overlap <- overlap_length # Number of points to overlap on each side

# Initialize a list to store the submatrices
submatrices <- list()

# Loop through the original matrix to extract submatrices
for (i in 1:(matrix_rows /(overlap)- 1)){
    for (j in 1:(matrix_cols /(overlap) - 1)){
        submatrix <- sample_matrix[((i-1)*overlap+1):((i-1)*overlap + submatrix_rows ), ((j-1)*overlap+1):((j-1)*overlap+ submatrix_cols )]
        submatrices[[length(submatrices) + 1]] <- submatrix
    }
}

##### parallel estimation ########
cl <- makeCluster(num_core,type="FORK")
registerDoParallel(cl)


# Execute the function in parallel
results <- foreach(i=1:n_split_num, .combine = rbind) %:%
  foreach(j = 1:length(submatrices), .combine = rbind) %dopar%
  {
    x_sub <- x[subject_groups[[i]], ]
    z_sub <- z[subject_groups[[i]], as.vector(submatrices[[j]] ) ]

    locations_tmp <- locations[as.vector(submatrices[[j]] ), ]
    dim_tmp<- dim(locations_tmp)[1]

    kp <- 1

    K_gram_tmp <- outer( c(1:dim_tmp), c(1:dim_tmp), function(x, y) exp(-( (locations_tmp [x,1]-locations_tmp [y,1])^2 +(locations_tmp [x,2]-locations_tmp [y,2])^2  )/kp)  )


    result <- functional_GroupL0_gradient(z_sub, x_sub, K_gram_tmp, lambda_range, sp_range, max_iter=200, nbasis=nbasis)$final_results
    return(list(list(i, j, result) ) ) 
  }

# Stop the parallel cluster
stopCluster(cl)

###### merge the results #########
#average_fitted_z <- vector("list", length=length(submatrices))

C_fun_recover_single <- function(px_ind){

  C_fun_recover <- matrix(0, nrow=mobs, ncol=length(submatrices))
  count_num <-  matrix(0, nrow=mobs, ncol=length(submatrices))


for (j in 1:length(submatrices)) {
  sub_data <- subset(results, sapply(results, function(x) x[[2]] == j))
  
  #average_f <- mean(sapply(sub_data, function(x) x[[3]]$fitted_vz ))
  average_C_fun_temp <- apply(sapply(sub_data, function(x) x[[3]]$C_function_est[[px_ind]] ),1,mean )
  
  #average_fitted_z[[j]] <- list(j = j, fitted_z = average_f)

  C_fun_recover[ as.vector(submatrices[[j]] ),j] = average_C_fun_temp
  count_num[ as.vector(submatrices[[j]] ),j] <- 1
}

count_non_zeros <- rowSums(count_num)

C_fun_final <- rowSums(C_fun_recover)/count_non_zeros

return(C_fun_final)

}

cl <- makeCluster(num_core)
registerDoParallel(cl)
C_fun_all <- foreach(j = 1:pn, .combine = rbind) %dopar%
  {
   C_fun_est <- C_fun_recover_single(j)
   }
stopCluster(cl)  


fitted_vz <- x%*% C_fun_all 

  
  
  list( C_fun_all=C_fun_all, results=results )
  
}