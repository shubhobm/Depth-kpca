## source functions required for depth-based kernel PCA
## author : SBM

## kpcaLocantore: function for spherical and depth kernel PCA
## Matrix Interface
## is.depth = FALSE (default) -> spherical (Locantore) PCA
## is.depth = TRUE -> depth PCA

kpcaLocantore = function(x, kernel = "rbfdot", kpar = list(sigma = 0.1),
	features = 0, th = 1e-4, delta=1e-3, is.depth=FALSE, na.action = na.omit,
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

  # center data matrix by spatial median
  # TODO: we can also do depth-based spatial median: but need theory
  # for that which hasn't been worked out yet
  sp <- spatial.median(x, kernel, delta, deps=mult)
#  x <- x - matrix(sp$mu, m, p, byrow=TRUE)

  ## obtain kernel matrix
  km <- kernelMatrix(kernel,x)
  ww <- matrix(sp$w, m, m, byrow=F)
  km.ww = km %*% ww
  kw <- km + t(ww) %*% km %*% ww - t(km.ww) - km.ww
  
  kdsqrt <- mult/sqrt(diag(kw))
  km <- outer(kdsqrt, kdsqrt) * kw

  # center kernel matrix
  kc <- t(t(km - colSums(km)/m) -  rowSums(km)/m) + sum(km)/m^2

  # compute eigenvectors
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

## spatial.median: computes the spatial median given data x
spatial.median <- function(x, kernel, delta=1e-3, deps=1)
{
    dime = dim(x)
    n=dime[1]
    p=dime[2]
    delta1=delta*sqrt(p)
    mu0=apply(x,2,median)
    h=delta1+1
    tt=0

# setting deps=mult and multiplying U below with deps2 will give depth-weighted spatial median:
# worth exploring
#    if(length(deps)>1){
#      deps2 = outer(deps,deps)
#    }
#    else{
#      deps2 = 1
#    }

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

## kdepth.SP: computes kernelized spatial depth for points in xx given data x and a kernel
kdepth.SP = function(xx, x, kernel="rbfdot", kpar = list(sigma = 0.1)){
  n = nrow(x)
  
  # initialize kernel matrix
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

## reconKwok: function to reconstruct pre-image given a kpca model,
## by Kwok and Tsan method: stable
reconKwok = function(model, newdata, k=10){

	# get kernel matrices, based on what type of pca we are doing
	if(class(model)=="kpcaLocantore"){
		xdata = model$xmatrix
		kern = model$kernelf
		P = model$pcv
		eig = model$eig
		
		# get kernel matrix and center accordingly
		kx = kernelMatrix(kern, newdata, xdata)
		km = kernelMatrix(kern, xdata)
		kmbar = as.numeric(km %*% model$centerwts)
		kx = t(t(kx - kmbar) * model$kdsqrt)
		km = t(t(km - kmbar) * outer(model$kdsqrt, model$kdsqrt)) # optional (*)
	}
	else{
		xdata = xmatrix(model)
		kern = kernelf(model)
		P = pcv(model)
		eig = eig(model)
		kx = kernelMatrix(kern, newdata, xdata)
		km = kernelMatrix(kern, xdata)
		kmbar = colMeans(km)
	}
	kx = t(kx)
	kmbar = colMeans(km) # optional if km is scaled above (*)
	kxPlus = kx + kmbar
	kxMinus = kx - kmbar

	# compute summands
	N = nrow(xdata)
	onex = rep(1, N)
	H = diag(onex) - onex %*% t(onex)/N

	M = P %*% (t(P) * eig)
	HtMH = t(H) %*% M %*% H
	s1 = diag(t(kxPlus) %*% HtMH %*% kxMinus) + mean(kmbar)
	s2 = diag(km)
	s3 = km %*% HtMH %*% kxMinus + kmbar
	
	# get distances
	df2 = (outer(s1,s2, "+") - t(2*s3))
	d2 = -kpar(kern)$sigma * log(1-df2/2)

	# reconstruct
	z1 = newdata
	Hxdata = H %*% xdata
	xmean = matrix(1, ncol=N, nrow=1) %*% xdata / N
	for(i in 1:nrow(newdata)){
		isort = sort((d2[i,]), index.return=TRUE)
		iX = Hxdata[isort$ix[1:k],]
		id = isort$x[1:k]
		isvd = svd(iX)
		id0 = apply((isvd$u * isvd$d)^2, 1, sum)
		z1[i,] = -0.5 * t((isvd$u * 1/isvd$d) %*% t(isvd$v)) %*% (id-id0)
		z1[i,] = z1[i,] + xmean
	}
	z1
}

## reconMika: function to reconstruct pre-image given a kpca model,
## by Mika et al method: unstable
reconMika = function(model, newdata, ep=1e-7, maxit=1e3){

	# get kernel matrices, based on what type of pca we are doing
	if(class(model)=="kpcaLocantore"){
		xdata = model$xmatrix
		kern = model$kernelf
		P = model$pcv

		km = kernelMatrix(kern, newdata, xdata)
		km.w = as.numeric(km %*% model$centerwts)
		km = t(t(km - km.w) * model$kdsqrt)
	}
	else{
		xdata = xmatrix(model)
		kern = kernelf(model)
		P = pcv(model)
		km = kernelMatrix(kern, newdata, xdata)

	}
	gamma = km %*% P %*% t(P)

	# reconstruct iteratively
	delta = 1e3
	z0 = newdata
	iter=0
	while(delta>ep & iter<maxit){
		kmg = kernelMatrix(kern, z0, xdata)
		if(class(model)=="kpcaLocantore"){
			kmg.w = as.numeric(kmg %*% model$centerwts)
			kmg = t(t(kmg - kmg.w) * model$kdsqrt)
		}
	
		km.gamma = kmg * gamma
		numerat = km.gamma %*% xdata
		denom = apply(km.gamma, 1, sum)
		z1 = numerat/as.numeric(denom)
		diff = abs((z1-z0)/z0)
		delta = mean(diff, na.rm=T)
		iter = iter+1
		z0=z1
	}
	z1
}

#s = z1
#sm = matrix(s, nrow=480, byrow=F)
#plot(imagematrix(sm))

## my.mvrnorm: function to generate from multivariate normal
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

## tau: function to get best gaussian kernel for data by knn method
## uses sigma = .5 times mean knn distance
tau = function(x, k=5){
	dist = get.knn(x, k=k)
	mean(dist$nn.dist)/2	
}

