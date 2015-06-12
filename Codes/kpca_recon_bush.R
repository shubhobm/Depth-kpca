### Image reconstruction using kpca
rm(list=ls())
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/KDPCA")

library(kernlab)
library(jpeg)
library(ripa)
library(EBImage)
library(FNN)
source("misc_functions1.R")

## read in images
img.list = list.files("all")
n = 40
a=100; b=100

s=1:n
ns=10
img.mat = matrix(0, nrow=length(s), ncol=a*b)
for(i in 1:length(s)){ 
	t = readJPEG(paste0("all/",img.list[[i]]), native=F)[,,1]
	t = resize(t, w=a, h=b)
	t = imagematrix(t)

	# speckle noise for sample 10
	if(i==10){
		samps = c(11:20, 31:40)
		t1 = t
#		t[111:140, 111:140] = 0
		t[samps[1:ns], samps[(ns+1):(2*ns)]] = 1

	}

	img.mat[i,] = as.numeric(t)
}

# store normal data
img.mat0 = img.mat
img.mat0[10,] = as.numeric(t1)

#img.mat = img.mat0

k = rbfdot(1/(2*tau(img.mat)))
k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), img.mat))/.0005))
k = rbfdot(.0005)

system.time(mod.kpc <- kpca(img.mat, features=9, kernel=k))
system.time(mod.kpcLoc <- kpcaLocantore(img.mat, features=9, kernel=k))
system.time(mod.kpcDep <- kpcaLocantore(img.mat, features=9, is.depth=TRUE, kernel=k))

# plot all results for one some image sample
i = 10
par(mfrow=c(2,2))
sm = matrix(img.mat[i,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconKwok(mod.kpc, img.mat0[i:(i+1),], k=5)
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconKwok(mod.kpcLoc, img.mat0[i:(i+1),], k=5)
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconKwok(mod.kpcDep, img.mat0[i:(i+1),], k=5)
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))
par(mfrow=c(1,1))

i = 10
par(mfrow=c(2,2))
sm = matrix(img.mat[i,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconMika(mod.kpc, img.mat0[i:(i+1),])
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconMika(mod.kpcLoc, img.mat0[i:(i+1),])
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconMika(mod.kpcDep, img.mat0[i:(i+1),])
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))
par(mfrow=c(1,1))