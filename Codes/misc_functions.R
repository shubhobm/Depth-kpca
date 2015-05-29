# Kdepth program
## function to generate from multivariate normal
my.mvrnorm = function(n, mu, Sigma){
  p = length(mu)
  # compute square root of covariance matrix
  eo=eigen(Sigma, symmetric=TRUE)
  sigma.sqrt=eo$vec%*%diag(eo$val^0.5)%*%t(eo$vec)
  
  # generate random normals from runif by box-muller transform
  rnorm.vec = sqrt(-2*log(runif(n*p)))*cos(2*pi*runif(n*p))
  
  # generate sample matrix
  sample.matrix = matrix(rep(mu, n), nrow=n, byrow=T) +
    matrix(rnorm.vec, nrow=n, ncol=p)%*%sigma.sqrt
  return(sample.matrix)
}

ones = function(m,n){
  matrix(1, nrow=m, ncol=n)
}

Kern = function(a,b,sigma) exp(-sum((a - b)^2) / sigma)

# SLOW... writwe in C
kdepth.SP = function(x, xx, sigma){
  n = nrow(xx)
  K.mat = matrix(0, nrow=n, ncol=n)
  
  ## initialize kernel matrix
  for(i in 1:nrow(K.mat)){
    for(j in 1:i){
      K.mat[i,j] = Kern(xx[i,], xx[j,], sigma)
      K.mat[j,i] = K.mat[i,j]
    }
  }
  
  # calculate depth values for rows of x
  nn = nrow(x)
  dep.vec = rep(0, nn)
  for(i in 1:nn){
    ix = x[i,]
    
    xi = rep(0, n)
    for(j in 1:n) xi[j] = Kern(ix, xx[j,], sigma)
    
    del = sqrt(Kern(ix, ix, sigma) + diag(K.mat) - 2*xi)
    z = ifelse(del==0, 0, 1/del)
    
    iK.mat = Kern(ix, ix, sigma) + K.mat
    for(j in 1:nrow(K.mat)){
      for(k in 1:j){
        iK.mat[j,k] = iK.mat[j,k] - xi[j] - xi[k]
        iK.mat[k,j] = K.mat[j,k]
      }
    }
    dep.vec[i] = 1 - sqrt(t(z) %*% iK.mat %*% z) / n
  }
  
  dep.vec
}