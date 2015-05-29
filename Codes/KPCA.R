#'Kernel PCA
#'
#'TODO: should have similar interface as princomp from base R
#'
#'@param x A matrix of dim nxp.
#'@param kernel A kernel function, gaussian by default.
#'@param nfeatures number of features to select.
#'@param ... other parameters needed for kernel function.
#'
#'@examples
#'out <- KPCA(as.matrix(iris[,-5]),nfeatures=2)
#'out$proj
#'@author: Abhirup Mallik

gausskern <- function(x, y=x, sigma = 1){
  as.numeric(exp(-crossprod(x-y)/(2*sigma^2)))
}

kernmat <- function(kern, x, y=x,...){
  n <- NROW(x)
  m <- NROW(y)
  ret <- matrix(NA,nrow=n,ncol=m)
  for(i in 1:n){
    for(j in 1:m){
     ret[i,j] <- kern(x[i,],y[j,],...) 
    }
  }
  ret
}

KPCA <- function(x, kern = gausskern, nfeatures=1, ...){
  n <- NROW(x)
  km <- kernmat(kern,x,...)
  kc <- t(t(km - colSums(km)/n) -  rowSums(km)/n) + sum(km)/n^2
  res <- eigen(kc/n,symmetric=TRUE)
  ret <- list()
  ret$pc <- pc <- t(t(res$vectors[,1:nfeatures])/sqrt(res$values[1:nfeatures]))
  ret$eigenval <- res$values[1:nfeatures]
  ret$proj <- kc %*% pc
  ret$x <- x
  ret
}