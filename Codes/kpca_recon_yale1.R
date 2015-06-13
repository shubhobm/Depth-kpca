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
n = 20
a = 192; b = 168
img.mat = matrix(0, nrow=3*n, ncol=a*b)
ns=ceiling(.5*a)
set.seed(12062015)

img.list = list.files("../CroppedYale/yaleB01")
for(i in 4:(n+3)){ 
	t = read.pnm(paste0("../CroppedYale/yaleB01/",img.list[[i]]))
	t = imagematrix(getChannels(t))
	t = resize(t, w=a, h=b)

	# speckle noise for sample 10
	if(i==10){
		samps = sample(1:min(a,b), 2*ns, replace=T)
		t1 = t
		t[samps[1:ns], samps[(ns+1):(2*ns)]] = sample(c(0,1), ns^2, replace=T)
	}
	img.mat[i-3,] = as.numeric(t)
}
img.mat0 = img.mat
img.mat0[7,] = as.numeric(t1)

img.list = list.files("../CroppedYale/yaleB02")
for(i in 4:(n+3)){ 
	t = read.pnm(paste0("../CroppedYale/yaleB02/",img.list[[i]]))
	t = imagematrix(getChannels(t))
	t = resize(t, w=a, h=b)

	# speckle noise for sample 10
	if(i==10){
		samps = sample(1:min(a,b), 2*ns, replace=T)
		t1 = t
		t[samps[1:ns], samps[(ns+1):(2*ns)]] = sample(c(0,1), ns^2, replace=T)
	}
	img.mat[n+i-3,] = as.numeric(t)
}
img.mat0[n+7,] = as.numeric(t1)

img.list = list.files("../CroppedYale/yaleB28")
for(i in 4:(n+3)){ 
	t = read.pnm(paste0("../CroppedYale/yaleB28/",img.list[[i]]))
	t = imagematrix(getChannels(t))
	t = resize(t, w=a, h=b)

	# speckle noise for sample 10
	if(i==10){
		samps = sample(1:min(a,b), 2*ns, replace=T)
		t1 = t
		t[samps[1:ns], samps[(ns+1):(2*ns)]] = sample(c(0,1), ns^2, replace=T)
	}
	img.mat[2*n+i-3,] = as.numeric(t)
}
img.mat0[2*n+7,] = as.numeric(t1)

i = 7 # change i to 7, 27 or 47 for 3 different persons
sm = matrix(img.mat[i,], nrow=a, byrow=F)
plot(imagematrix(sm))

## choose kernel
#k = rbfdot(1/(2*tau(img.mat)))
k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), img.mat))/.0005))
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

# plot first 2 PCs
cols = c(rep("red",n), rep("blue",n), rep("green",n))
par(mfrow=c(1,3))
plot(rotated(mod.kpc[[2]]), col=cols, pch=19)
plot(mod.kpcLoc[[2]]$rotated, col=cols, pch=19)
plot(mod.kpcDep[[2]]$rotated, col=cols, pch=19)
par(mfrow=c(1,1))

###### plot all results for one some image sample
### Reconstruction by Kwok and Tsan method
### img.mat = denoising, img.mat0 = basis retrieval
defaultPar = par()
par(mfcol=c(3,7), mai=rep(0,4))
for(nf in 2*(1:5)){
	s = reconKwok(mod.kpc[[nf]], img.mat[i:(i+1),], k=5)
	sm = matrix(s[1,], nrow=a, byrow=F)
	plot(imagematrix(sm))

	s = reconKwok(mod.kpcLoc[[nf]], img.mat[i:(i+1),], k=5)
	sm = matrix(s[1,], nrow=a, byrow=F)
	plot(imagematrix(sm))

	s = reconKwok(mod.kpcDep[[nf]], img.mat[i:(i+1),], k=5)
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
