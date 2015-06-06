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

img0 = matrix(as.numeric(img), nrow=60, byrow=F)

img[201:240, 301:340] = 1
img[301:340, 201:240] = 1
img[1:40, 101:140] = 1

img.vec = as.numeric(img)
#set.seed(20150604)
#samp = sample(1:length(img.vec), .05*length(img.vec))
#img.vec[samp] = runif(length(samp))

img.mat = matrix(img.vec, nrow=60, byrow=F)
dim(img.mat)
#img.mat = scale(img.mat, center=F, scale=F)

k = rbfdot(1)
system.time(mod.kpc <- kpca(img.mat, features=3, kernel=k))
system.time(mod.kpcLoc <- kpcaLocantore(img.mat, features=3, kernel=k))
system.time(mod.kpcDep <- kpcaLocantore(img.mat, features=3, is.depth=TRUE, kernel=k))

# plot all results for one some image sample
par(mfrow=c(2,2))
sm = matrix(img.mat, nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconGood(mod.kpc, img0, k=50)
sm = matrix(s, nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconGood(mod.kpcLoc, img0, k=50)
sm = matrix(s, nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconGood(mod.kpcDep, img0, k=50)
sm = matrix(s, nrow=480, byrow=F)
plot(imagematrix(sm))
par(mfrow=c(1,1))