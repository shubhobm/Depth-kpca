### Image reconstruction using kpca, uses block kpca
rm(list=ls())
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/KDPCA")

library(kernlab)
library(png)
library(ripa)
library(EBImage)
library(FNN)
source("misc_functions1.R")

## read in image
img = readPNG(paste0("images/aug4.43.1.png"))
#t = resize(t, 200, 200)
#t = imagematrix(t)

img0 = matrix(as.numeric(img), nrow=30, byrow=F)

img.vec = as.numeric(img)
#samps = sample(1:length(img.vec), length(img.vec)*1e-3)
#img.vec[samps] = sample(0:1, length(samps), replace=T)
#set.seed(20150604)
samp = sample(1:length(img.vec), .05*length(img.vec))
img.vec[samp] = rnorm(length(samp))*.09

#img.vec = img.vec + rnorm(length(img.vec))*.09

img.mat = matrix(img.vec, nrow=30, byrow=F)
dim(img.mat)

k = rbfdot(1/(2*tau(img.mat)))
k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), img.mat))/.0005))

system.time(mod.kpc <- kpca(img.mat, features=20, kernel=k))
system.time(mod.kpcLoc <- kpcaLocantore(img.mat, features=20, kernel=k))
system.time(mod.kpcDep <- kpcaLocantore(img.mat, features=20, is.depth=TRUE, kernel=k))

# plot all results for one some image sample
par(mfrow=c(2,2))
sm = matrix(img.vec , nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconKwok(mod.kpc, img0, k=5)
sm = matrix(s, nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconKwok(mod.kpcLoc, img0, k=5)
sm = matrix(s, nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconKwok(mod.kpcDep, img0, k=5)
sm = matrix(s, nrow=480, byrow=F)
plot(imagematrix(sm))
par(mfrow=c(1,1))

par(mfrow=c(2,2))
sm = matrix(img.vec , nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconMika(mod.kpc, img0)
sm = matrix(s, nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconMika(mod.kpcLoc, img0)
sm = matrix(s, nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconMika(mod.kpcDep, img0)
sm = matrix(s, nrow=480, byrow=F)
plot(imagematrix(sm))
par(mfrow=c(1,1))