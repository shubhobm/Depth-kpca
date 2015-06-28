#' Simulation Experiment 2
#' Framework of Yang et al.
#' @author Abhirup Mallik
#' 


angdist <- function(x,y){
  theta <- crossprod(x,y)/(sqrt(crossprod(x)*crossprod(y)))
  acos(theta)
}

vecdist <- function(x,y){
  sqrt(crossprod(x,y))
}

set.seed(12345)
n <- 101
x <- runif(n)
r <- rnorm(n)
y <- 2*sin(2*pi*x)+0.5+0.1*r
yout <- y
#number of outliers
numout <- 30
outliers <- rnorm(numout,mean = 5,sd=1)
#randomly replace points in x by outliers
indexout <- sample.int(n,size=numout)
yout[indexout] <- outliers
data <- cbind(x,yout)
#Distance between eigen vectors of KPCA and other PCA

#Cosine distance:
#Kernel Choise
k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), data))/.0005))
#KPCA
pcnum <- 1
mod.kpc <- kpca(data, features=pcnum, kernel=k)
mod.kpcLoc <- kpcaLocantore(data, features=pcnum, kernel=k)
mod.kpcDep <- kpcaLocantore(data, features=pcnum, is.depth=TRUE, kernel=k)

## Study for increasing number of outliers from 0 to 50

theta.Loc <- numeric(51)
theta.Dep <- numeric(51)
dist.Loc <- numeric(51)
dist.Dep <- numeric(51)

for(i in 1:51){
  numout <- i -1
  outliers <- rnorm(numout,mean = 5,sd=1)
  #randomly replace points in x by outliers
  indexout <- sample.int(n,size=numout)
  yout[indexout] <- outliers
  data <- cbind(x,yout)
  #Distance between eigen vectors of KPCA and other PCA
  
  #Cosine distance:
  #Kernel Choise
  k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), data))/.0005))
  #KPCA
  pcnum <- 1
  mod.kpc <- kpca(data, features=pcnum, kernel=k)
  mod.kpcLoc <- kpcaLocantore(data, features=pcnum, kernel=k)
  mod.kpcDep <- kpcaLocantore(data, features=pcnum, is.depth=TRUE, kernel=k)
  theta.Loc[i] <- angdist(pcv(mod.kpc),mod.kpcLoc$pcv)
  theta.Dep[i] <- angdist(pcv(mod.kpc),mod.kpcDep$pcv)
  dist.Loc[i] <- vecdist(pcv(mod.kpc),mod.kpcLoc$pcv)
  dist.Dep[i] <- vecdist(pcv(mod.kpc),mod.kpcDep$pcv)
}

