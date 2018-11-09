#### BIOSTAT615 Final Project
#### Simulation in R 

##### A R function of "Complete Pivoting LU Decomposition" should be added for getting sigma

simulation = function(p,r,u,d,N){
  ## Generate Beta
  Gamma = matrix(runif(u*r),nrow=r)
  #Gamma = scale(Gamma, center=T, scale = T)
  GammaU = svd(Gamma)$u
  eta = matrix(runif(u*d),nrow=u)
  etaU = svd(eta)$u
  B = matrix(runif(d*p),nrow=d)
  
  beta = (GammaU%*%etaU%*%B)
  beta = beta / norm(beta,type='f')
  
 ## Generate Sigma
  library(Matrix)
  omega =  matrix(-0.9 , ncol = u, nrow = u)
  for(i in 1:u){
    for(j in 1:u){
      omega[i,j] = omega[i,j]^abs(i-j)}}
  
  omega0 = matrix(-0.5 , ncol = r-u, nrow = r-u)
  for(i in 1:(r-u)){
    for(j in 1:(r-u)){
      omega0[i,j] = omega[i,j]^abs(i-j)}}
  omega0 = 5*omega0
  
  sigma = GammaU %*% omega %*% t(GammaU)
  
  ## Generate X and Y
  
  X = scale(matrix(rnorm(p*N), nrow=p),center = T)
  
  Z = scale(matrix(rnorm(r*N), nrow=r),center = T)
  c_sigma = t(chol(sigma)) ## Factor L in the Cholesky decomposition
  error = c_sigma %*% Z
  Y = (beta %*% X) + error
}
