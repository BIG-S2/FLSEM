library(MASS)
Gen_data_new <- function(n, pn,   locations, rho2, cor, B, k=10, nu=2, rho1=0,eigen_center=1, sd_noise=1){
  # n: sample size
  # pn: number of scalar variables
  # k: number of basis
  # locations: observation points of the functional variable
  # nu: smoothness of the functional variable
  # rho1: correlation of the FPC scores and the scalar variables
  # rho2: correlation of the scalar variables
  # eigen_center: to generate disordered representative basis
  # sd_noise: sd of the error
  # Cor: correlation between the errors
  # B: signal of the functional coefficient
  lent <- length(locations)
  
  ######################################
  ## xi1-xi4 and  z ####################
  ######################################
  pn4=pn+4;   ## pn is the number of scalar variates
  sigma=matrix(0,pn4,pn4)  ## generate a zeros matrix
  sigma[1,1]=16; sigma[2,2]=16*2^(-nu); sigma[3,3]=16*3^(-nu); sigma[4,4]=16*4^(-nu);
  
  ## the four leading FPC scores correlate with z
  for(i in 1:4){
    for(j in 5:pn4){
      sigma[i,j]=rho1^(abs(i-j+4)+1)*sqrt(sigma[i,i]);  
      sigma[j,i]=sigma[i,j];
    }
  }
  
  for(i in 5:pn4){
    for(j in 5:pn4){
      sigma[i,j]=rho2^(abs(i-j));} 
  }
  
  xi_z=mvrnorm(n=n,rep(0,pn4),sigma)  
  xi_14=xi_z[,1:4]; z=xi_z[,5:pn4] 
  
  
  
  #####################################
  ## generate  xi and x ###############
  #####################################
  phi=function(k,t){
    if(k%%2==0)  return(sqrt(2)*sin((k-1)*pi*t))  else  return(sqrt(2)*cos(k*pi*t));
  }
  
  
  xi=apply(matrix(1:k,1,k), 2 , function(x)rnorm( n,0,1*( abs( (x-eigen_center))+1 )^(-nu/2) ))   ### n*k
  
  
  j_trans=as.vector( matrix(rep(c(1:k),lent),lent,k,byrow=TRUE) )
  basis= matrix( mapply( function(x, t){ phi(x,t)}, j_trans, rep(locations,k) ), k, lent, byrow=TRUE)
  
  
  ### generate error distribution #######
  
  errors <- mvrnorm(n, rep(0, 3),  matrix(c(1,0, cor,0,0.64,0.8*cor, cor,0.8*cor, sd_noise), 3,3)  ) #matrix(c(16,0, 4*cor,0,4,2*cor, 4*cor,2*cor, 1), 3,3)
  xi_1 <- errors[, 1]
  xi_2 <- errors[, 2]  #rnorm(n, 0, sqrt(4))
  epsilon <- errors[, 3]
  
  #### generate betas ########################
  beta_1 <- 2*locations^2
  beta_2 <- 3*cos( pi*locations /2  ) 
  beta_3 <- 0.5*sqrt(2)*sin(pi*locations /2) 
  beta_4 <- exp(-locations)
  beta_5 <- 2 + 1*locations
  
  C_all <- rbind(beta_1, beta_2, beta_3, beta_4, beta_5)
  
  #######################################
  
  x_aux <- z[, 1:5]%*%C_all
  x <-x_aux +  cbind(xi_1, xi_2)%*%basis[1:2, ]
  
  
  ##### generate beta ###########
   beta<- B*exp(-locations)
    
  
  
  
  ### generate response #####
  gamma <- c(7, 0, 0, 0, 0, 5.5, 4, 3.5, 5, 4.5, rep(0,pn-10) )  # for high-dimensional 
  #gamma <- c(2, 1, 0, 0, 0, 5.5, 4, 3.5, 5, 4.5, rep(0,pn-10) ) # for pn 20
  
  y <- as.numeric( x%*%as.matrix(beta ))/lent + as.numeric( z%*%gamma ) + epsilon 
  
  K_gram <- t(basis)%*%basis
  
  
  list(y=y, x=x, z=z, x_aux =x_aux, beta_true=beta, C_all=C_all, gamma=gamma, K_gram=K_gram, idx_true_X = c(rep(1,5), rep(0, pn-5)), 
       idx_true_Y =as.numeric(gamma >0) )
  
}