## kpcaLocantore function
## author : sbm
## Same syntax as kpca in package kernlab
## v1: only matrix inputs supported
## Matrix Interface
kpcaLocantore = function(x, kernel = "rbfdot", kpar = list(sigma = 0.1),
	features = 0, th = 1e-4, delta=0.001, is.depth=FALSE, na.action = na.omit,
	...)
{
  x <- na.action(x)
  x <- as.matrix(x)
  m <- nrow(x)
  p <- ncol(x)
  ret <- list()
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  # determine whether to do depth kpca
  if(is.depth){
	mult <- kdepth.SP(x, x, kernel)
	mult <- exp(-mult)
  }
  else{
	mult <- 1
  }

  ## center data matrix by spatial median
  sp <- spatial.median(x, kernel, delta, deps=mult)
#  x <- x - matrix(sp$mu, m, p, byrow=TRUE)
  km <- kernelMatrix(kernel,x)
  ww <- matrix(sp$w, m, m, byrow=F)
  km.ww = km %*% ww
  kw <- km + t(ww) %*% km %*% ww - t(km.ww) - km.ww
  
  kdsqrt <- mult/sqrt(diag(kw))
  km <- outer(kdsqrt, kdsqrt) * kw

  ## center kernel matrix
  kc <- t(t(km - colSums(km)/m) -  rowSums(km)/m) + sum(km)/m^2

  ## compute eigenvectors
  res <- eigen(kc/m,symmetric=TRUE)
  
  if(features == 0)
    features <- sum(res$values > th)
#  else 
#    if(res$values[features] < th)
#      warning(paste("eigenvalues of the kernel matrix are below threshold!"))
 
  ret$pcv <- t(t(res$vectors[,1:features])/sqrt(res$values[1:features]))
  ret$eig <- res$values[1:features]
  names(ret$eig) <- paste("Comp.", 1:features, sep = "")
  ret$rotated <- kc %*% ret$pcv
  ret$kernelf <- kernel
  ret$xmatrix <- x
  ret$is.depth = is.depth
  ret$kdsqrt = kdsqrt
  ret$centerwts = sp$w
  class(ret) = "kpcaLocantore"
  ret
}

## computes the spatial median
spatial.median <- function(x, kernel, delta=1e-3, deps=1)
{
    dime = dim(x)
    n=dime[1]
    p=dime[2]
    delta1=delta*sqrt(p)
    mu0=apply(x,2,median)
    h=delta1+1
    tt=0
    if(length(deps)>1){
      deps2 = outer(deps,deps)
    }
    else{
      deps2 = 1
    }

    while(h>delta1)
    {
        tt=tt+1
        TT=matrix(mu0,n,p,byrow=TRUE)
        U=kernelMatrix(kernel, (x - TT))
        w=sqrt(apply(U,1,sum))
        w0=median(w)
        ep=delta*w0

        z=(w<=ep)
        w[z]=ep
        w[!z]=1/w[!z]
        w=w/sum(w)
        x1=x
        for(i in 1:n)
            x1[i,]=w[i]*x[i,]
        mu=apply(x1,2,sum)
        h=sqrt(sum((mu-mu0)^2))
        mu0=mu
    }
    out=list(mu=mu0,ep=ep, w=w)
    out
}

kdepth.SP = function(xx, x, kernel="rbfdot", kpar = list(sigma = 0.1)){
  n = nrow(x)
  
  ## initialize kernel matrix
  if(!is(kernel,"kernel"))
    {
      if(is(kernel,"function")) kernel <- deparse(substitute(kernel))
      kernel <- do.call(kernel, kpar)
    }
  if(!is(kernel,"kernel")) stop("kernel must inherit from class `kernel'")

  km.x <- kernelMatrix(kernel, x)
  km.xx <- kernelMatrix(kernel, xx)
  km.xx.x <- kernelMatrix(kernel, xx, x)

  # calculate depth values for rows of xx
  nn = nrow(xx)
  dep.vec = rep(0, nn)
  
  for(i in 1:nn){
    xi = km.xx.x[i,]
    del = sqrt(km.xx[i,i] + diag(km.x) - 2*xi)
    z = (del==0)
    del[z] = 0
    del[!z] = 1/del[!z]
    
    km.i = km.xx[i,i] + km.x - outer(xi, xi, function(x,y) x+y)
    dep.vec[i] = 1 - sqrt(t(del) %*% km.i %*% del) / n
  }
  
  dep.vec
}

## function to reconstruct pre-image given a kpca model
## unstable... looking for new algo
reconGood = function(model, newdata, k=10){

	if(class(model)=="kpcaLocantore"){
		xdata = model$xmatrix
	
		# get kernel matrix and center accordingly
		km = kernelMatrix(kernel=model$kernelf, newdata, xdata)
		km.w = as.numeric(km %*% model$centerwts)
		km = t(t(km - km.w) * model$kdsqrt)
	}
	else {
		xdata = xmatrix(model)
		km = kernelMatrix(kernel=kernelf(model), newdata, xdata)
	}
	kc <- t(t(km - colMeans(km)) -  rowMeans(km)) + mean(km)

	z1 = newdata
	xmean = matrix(1, ncol=nrow(xdata), nrow=1) %*% xdata /nrow(xdata)
	for(i in 1:nrow(newdata)){
		isort = sort(kc[i,], index.return=TRUE)
		iX = xdata[isort$ix[1:k],]
		id = isort$x[1:k]^2
		isvd = svd(iX)
		id0 = apply((isvd$d * isvd$u)^2, 2, mean)
		z1[i,] = -0.5 * t((isvd$u * 1/isvd$d) %*% t(isvd$v)) %*% (id-id0)
		z1[i,] = z1[i,] + xmean
	}
	z1
}

reconGood1 = function(model, newdata, k=10){

	if(class(model)=="kpcaLocantore"){
		xdata = model$xmatrix
		kern = model$kernelf
		P = model$pcv

		km = kernelMatrix(kern, newdata, xdata)
		km.w = as.numeric(km %*% model$centerwts)
		km = t(t(km - km.w) * model$kdsqrt)
	}
	else {
		xdata = xmatrix(model)
		kern = kernelf(model)
		P = pcv(model)

		km = kernelMatrix(kern, newdata, xdata)
	}
	kx = kernelMatrix(kern, xdata)
	kmc <- t(t(km - colMeans(km)) -  rowMeans(kx)) + mean(km)
	kxc <- t(t(kx - colMeans(kx)) -  rowMeans(kx)) + mean(kx)

	# get 3 parts of the sum
	s1 = apply((kmc %*% P)^2, 1, sum)
	s2 = diag(kxc)
	s3 = kmc %*% P %*% t(P)
	dist = outer(s1, s2, FUN="+") - 2*s3

	z1 = newdata
	xmean = matrix(1, ncol=nrow(xdata), nrow=1) %*% xdata /nrow(xdata)
	for(i in 1:nrow(newdata)){
		isort = sort(kmc[i,], index.return=TRUE)
		iX = xdata[isort$ix[1:k],]
		id = isort$x[1:k]^2
		isvd = svd(iX)
		id0 = apply((isvd$d * isvd$u)^2, 2, mean)
		z1[i,] = -0.5 * t((isvd$u * 1/isvd$d) %*% t(isvd$v)) %*% (id-id0)
		z1[i,] = z1[i,] + xmean
	}
	z1
}


#s = z1
#sm = matrix(s[1,], nrow=480, byrow=F)
#plot(imagematrix(sm))

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

## function to get best gaussian kernel for data
## uses sigma = .5 times mean knn distance
tau = function(x, k=5){
	dist = get.knn(x, k=k)
	mean(dist$nn.dist)/2	
}

