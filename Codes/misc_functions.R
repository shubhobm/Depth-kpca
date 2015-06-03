## kpcaLocantore function
## author : sbm
## Same syntax as kpca in package kernlab
## v1: only matrix inputs supported
############ START OF CLASS DEFINITION ###################
########################################################

##kpcaLocantore object
setClass("prc", representation(pcv = "matrix",
                               eig = "vector",
                               kernelf = "kfunction",
                               kpar = "list",
                               xmatrix = "input",
                               kcall = "ANY",
                               terms = "ANY"),contains="VIRTUAL")
#accessor functions 
if(!isGeneric("pcv")){
  if (is.function("pcv"))
    fun <- pcv
  else fun <- function(object) standardGeneric("pcv")
  setGeneric("pcv", fun)
}
setMethod("pcv", "prc", function(object) object@pcv)
setGeneric("pcv<-", function(x, value) standardGeneric("pcv<-"))
setReplaceMethod("pcv", "prc", function(x, value) {
  x@pcv <- value
  x
})

if(!isGeneric("eig")){
  if (is.function("eig"))
    fun <- eig
  else fun <- function(object) standardGeneric("eig")
  setGeneric("eig", fun)
}
setMethod("eig", "prc", function(object) object@eig)
setGeneric("eig<-", function(x, value) standardGeneric("eig<-"))
setReplaceMethod("eig", "prc", function(x, value) {
  x@eig <- value
  x
})

if(!isGeneric("kernelf")){
  if (is.function("kernelf"))
    fun <- kernelf
  else fun <- function(object) standardGeneric("kernelf")
  setGeneric("kernelf", fun)
}
setMethod("kernelf", "prc", function(object) object@kernelf)
setGeneric("kernelf<-", function(x, value) standardGeneric("kernelf<-"))
setReplaceMethod("kernelf", "prc", function(x, value) {
  x@kernelf <- value
  x
})

if(!isGeneric("xmatrix")){
  if (is.function("xmatrix"))
    fun <- xmatrix
  else fun <- function(object) standardGeneric("xmatrix")
  setGeneric("xmatrix", fun)
}
setMethod("xmatrix", "prc", function(object) object@xmatrix)
setGeneric("xmatrix<-", function(x, value) standardGeneric("xmatrix<-"))
setReplaceMethod("xmatrix", "prc", function(x, value) {
  x@xmatrix <- value
  x
})

if(!isGeneric("kcall")){
  if (is.function("kcall"))
    fun <- kcall
  else fun <- function(object) standardGeneric("kcall")
  setGeneric("kcall", fun)
}
setMethod("kcall", "prc", function(object) object@kcall)
setGeneric("kcall<-", function(x, value) standardGeneric("kcall<-"))
setReplaceMethod("kcall", "prc", function(x, value) {
  x@kcall <- value
  x
})


setMethod("terms", "prc", function(x, ...) x@terms)
setGeneric("terms<-", function(x, value) standardGeneric("terms<-"))
setReplaceMethod("terms", "prc", function(x, value) {
  x@terms <- value
  x
})

setClass("kpcaLocantore", representation(rotated = "matrix"),contains="prc")
#accessor functions 

if(!isGeneric("rotated")){
  if (is.function("rotated"))
    fun <- rotated
  else fun <- function(object) standardGeneric("rotated")
  setGeneric("rotated", fun)
}
setMethod("rotated", "kpcaLocantore", function(object) object@rotated)
setGeneric("rotated<-", function(x, value) standardGeneric("rotated<-"))
setReplaceMethod("rotated", "kpcaLocantore", function(x, value) {
  x@rotated <- value
  x
})

setGeneric("kpcaLocantore",function(x, ...) standardGeneric("kpcaLocantore"))
setMethod("kpcaLocantore", signature(x = "formula"),
function(x,  data = NULL, na.action = na.omit, ...)
{
    mt <- terms(x, data = data)
    if(attr(mt, "response") > 0) stop("response not allowed in formula")
    attr(mt, "intercept") <- 0
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    mf$formula <- mf$x
    mf$... <- NULL
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    na.act <- attr(mf, "na.action")
    Terms <- attr(mf, "terms")
    x <- model.matrix(mt, mf)
    res <- kpcaLocantore(x, ...)
    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("kpcaLocantore")
    kcall(res) <- cl
    attr(Terms,"intercept") <- 0
    terms(res) <- Terms
    if(!is.null(na.act)) 
        n.action(res) <- na.act
  
    return(res)
  })

############ END OF CLASS DEFINITION ###################
########################################################

## Matrix Interface
setMethod("kpcaLocantore",signature(x="matrix"),
          function(x, kernel = "rbfdot", kpar = list(sigma = 0.1), features = 0, th = 1e-4, delta=0.001, is.depth=FALSE, na.action = na.omit, ...)
{
  x <- na.action(x)
  x <- as.matrix(x)
  m <- nrow(x)
  p <- ncol(x)
  ret <- new("kpcaLocantore")
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
  kw <- km + t(ww) %*% km %*% ww - t(ww) %*% km - km %*% ww
  
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
 
  pcv(ret) <- t(t(res$vectors[,1:features])/sqrt(res$values[1:features]))
  eig(ret) <- res$values[1:features]
  names(eig(ret)) <- paste("Comp.", 1:features, sep = "")
  rotated(ret) <- kc %*% pcv(ret)
  kcall(ret) <- match.call()
  kernelf(ret) <- kernel
  xmatrix(ret) <- x
  return(ret)
})

## computes the spatial median
spatial.median <- function(x, kernel, delta, deps=1)
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

