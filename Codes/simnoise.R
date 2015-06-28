#' Simulation Setup
#'@author Abhirup Mallik


# Simulation framework of Mika et al
ntrain <- 100
ntest <- 33
sigma <- 0.25
set.seed(12345)
#Error vectors
err_pc <- numeric(11)
err_kpc <- numeric(11)
err_kpcLoc <- numeric(11)
err_kpcDep <- numeric(11)


mu <- matrix(runif(11*10,min=-1,max=1),nrow=11,ncol=10)

for(i in 1:11){
  x <- matrix(rnorm(10*(ntrain + ntest),mean=mu[i,],sd=sigma),ncol=10)
  x.test <- x[-(1:100),]
  x.train <- x[1:100,]
  #PCA
  res <- prcomp(x.train, center = TRUE, scale = FALSE)
  test.res <- predict(res,x.test)
  pcnum <- 2
  res.trunc <- test.res[,1:pcnum] %*% t(res$rotation[,1:pcnum])
  
  if(res$scale != FALSE){
    res.trunc <- scale(res.trunc, center = FALSE , scale=1/res$scale)
  }
  if(res$center != FALSE){
    res.trunc <- scale(res.trunc, center = -1 * res$center, scale=FALSE)
  }
  #Kernel Choise
  k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), x.train))/.0005))
  #KPCA
  mod.kpc <- kpca(x.train, features=pcnum, kernel=k)
  mod.kpcLoc <- kpcaLocantore(x.train, features=pcnum, kernel=k)
  mod.kpcDep <- kpcaLocantore(x.train, features=pcnum, is.depth=TRUE, kernel=k)
  #Reconstruct
  s_kpc.Kw = reconKwok(mod.kpc, x.test, k=5)
  s_kpcLoc.Kw <- reconKwok(mod.kpcLoc, x.test, k=5)
  s_kpcDep.Kw <- reconKwok(mod.kpcDep, x.test, k=5)
  #Mean squared Dist
  err_pc[i] <- norm(scale(res.trunc,center=mu[i,],scale = FALSE),type="F")
  err_kpc[i] <- norm(scale(s_kpc.Kw,center=mu[i,],scale = FALSE),type="F")
  err_kpcLoc[i] <- norm(scale(s_kpcLoc.Kw,center=mu[i,],scale = FALSE),type="F")
  err_kpcDep[i] <- norm(scale(s_kpcDep.Kw,center=mu[i,],scale = FALSE),type="F")
}

err <- rbind(err_pc,err_kpc,err_kpcLoc,err_kpcDep)
apply(err,1,mean)



# 
# #x <- matrix(rnorm(10*(ntrain + ntest),sd=sigma),ncol=10)
# x.test <- x[-(1:100),]
# x.train <- x[1:100,]
# 
# #Add noise to Training data
# 
# #Gaussian noise
# x.train.g <- x.train + matrix(rnorm(10*ntrain,sd=sigma/10),ncol=10)
# #salt & Pepper noise
# p <- 0.1
# noiseindex <- sample.int(10*ntrain,floor(10*ntrain*p/2))
# tmp <- as.numeric(x.train)
# tmp[noiseindex] <- -sign(tmp[noiseindex])*tmp[noiseindex]
# x.train.s <- matrix(tmp,nrow=ntrain,ncol = 10)
# 
# #PCA
# res <- prcomp(x.train, center = TRUE, scale = FALSE)
# test.res <- predict(res,x.test)
# pcnum <- 2
# res.trunc <- test.res[,1:pcnum] %*% t(res$rotation[,1:pcnum])
# 
# if(res$scale != FALSE){
#   res.trunc <- scale(res.trunc, center = FALSE , scale=1/res$scale)
# }
# if(res$center != FALSE){
#   res.trunc <- scale(res.trunc, center = -1 * res$center, scale=FALSE)
# }
# 
# k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), x.train))/.0005))
# system.time(mod.kpc <- kpca(x.train, features=pcnum, kernel=k))
# system.time(mod.kpcLoc <- kpcaLocantore(x.train, features=pcnum, kernel=k))
# system.time(mod.kpcDep <- kpcaLocantore(x.train, features=pcnum, is.depth=TRUE, kernel=k))
# # reconstruct
# 
# s_kpc.Kw = reconKwok(mod.kpc, x.test, k=5)
# s_kpcLoc.Kw <- reconKwok(mod.kpcLoc, x.test, k=5)
# s_kpcDep.Kw <- reconKwok(mod.kpcDep, x.test, k=5)
# 
# #s_kpc.mika <- reconMika(mod.kpc, x.test)
# #s_kpcLoc.mika <- reconMika(mod.kpcLoc, x.test)
# #s_kpcDep.mika <- reconMika(mod.kpcDep, x.test)
# 
# #Find the distance
# diff0 <- norm(x.test - res.trunc,type = "F")
# 
# 
# diff_kpc.Kw <- norm(x.test - s_kpc.Kw,type = "F")
# diff_kpcLoc.Kw <- norm(x.test - s_kpcLoc.Kw,type = "F")
# diff_kpcDep.Kw <- norm(x.test - s_kpcDep.Kw,type = "F")
# 
# #diff_kpc.mika <- norm(x.test - s_kpc.mika,type = "F")
# #diff_kpcLoc.mika <- norm(x.test - s_kpcLoc.mika,type = "F")
# #diff_kpcDep.mika <- norm(x.test - s_kpcDep.mika,type = "F")
# 
# result <- data.frame(kpc.Kw = diff_kpc.Kw/diff0,
#                      kpcLoc.Kw = diff_kpcLoc.Kw/diff0,
#                      kpcDep.Kw = diff_kpcDep.Kw/diff0
# #                     kpc.mika = diff_kpc.mika/diff0,
# #                     kpcLoc.mika = diff_kpcLoc.mika/diff0,
# #                     kpcDep.mika = diff_kpcDep.mika/diff0
#                      )
# result
# 
# ##PCA on Noisy Dataset:
# 
# res <- prcomp(x.train.g, center = TRUE, scale = FALSE)
# test.res <- predict(res,x.test)
# pcnum <- 2
# res.trunc <- test.res[,1:pcnum] %*% t(res$rotation[,1:pcnum])
# 
# if(res$scale != FALSE){
#   res.trunc <- scale(res.trunc, center = FALSE , scale=1/res$scale)
# }
# if(res$center != FALSE){
#   res.trunc <- scale(res.trunc, center = -1 * res$center, scale=FALSE)
# }
# k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), x.train))/.0005))
# 
# mod.kpc <- kpca(x.train, features=pcnum, kernel=k)
# mod.kpcLoc <- kpcaLocantore(x.train, features=pcnum, kernel=k)
# mod.kpcDep <- kpcaLocantore(x.train, features=pcnum, is.depth=TRUE, kernel=k)

# Why is every method performing worse than linear PCA?