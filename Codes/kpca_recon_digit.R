### Image reconstruction using kpca
### uses uale face data, one sample only: B01
rm(list=ls())
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/KDPCA")

library(kernlab)
library(pixmap)
library(ripa)
library(EBImage)
library(FNN)
source("misc_functions1.R")

## read in images: extended yale database B
## first 20 images of 3 persons each: 1,2 and 28: chosen so that facial features are very different
## 7th image is corrupted with 25% speckle noise for each person
a = 50; b = 50
ns=ceiling(.5*a)
set.seed(13062015)

digitdata = read.csv("C:/data/train_MNIST.csv")
digitdata[,-1] = digitdata[,-1]/255

i = 1
sm = matrix(as.numeric(digitdata[i,-1], nrow=28, byrow=T))
plot(imagematrix(sm))

## choose kernel
#k = rbfdot(1/(2*tau(img.mat)))
(k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), img.mat))/.0005)))
#k = rbfdot(.0005)

## Change number of features and compare plots
## basis retrieval and noise removal

mod.kpc = list()
mod.kpcLoc = list()
mod.kpcDep = list()
for(nf in 1:10){
	mod.kpc[[nf]] <- kpca(img.mat, features=nf, kernel=k)
	mod.kpcLoc[[nf]] <- kpcaLocantore(img.mat, features=nf, kernel=k)
	mod.kpcDep[[nf]] <- kpcaLocantore(img.mat, features=nf, is.depth=TRUE, kernel=k)
}

###### plot all results for one some image sample
### Reconstruction by Kwok and Tsan method
### img.mat = denoising, img.mat0 = basis retrieval
defaultPar = par()
par(mfcol=c(3,7), mai=rep(0,4))
for(nf in 2*(1:5)){
	s = reconKwok(mod.kpc[[nf]], img.mat[i:(i+1),], k=2)
	sm = matrix(s[1,], nrow=a, byrow=F)
	plot(imagematrix(sm))

	s = reconKwok(mod.kpcLoc[[nf]], img.mat[i:(i+1),], k=2)
	sm = matrix(s[1,], nrow=a, byrow=F)
	plot(imagematrix(sm))

	s = reconKwok(mod.kpcDep[[nf]], img.mat[i:(i+1),], k=2)
	sm = matrix(s[1,], nrow=a, byrow=F)
	plot(imagematrix(sm))
}
sm = matrix(img.mat0[i,], nrow=a, byrow=F)
plot(imagematrix(sm))
plot(imagematrix(sm))
plot(imagematrix(sm))

sm = matrix(img.mat[i,], nrow=a, byrow=F)
plot(imagematrix(sm))
plot(imagematrix(sm))
plot(imagematrix(sm))
par(defaultPar)
