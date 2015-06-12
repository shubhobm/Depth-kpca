### Image reconstruction using kpca, uses block kpca
rm(list=ls())
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/KDPCA")

library(kernlab)
library(pixmap)
library(ripa)
library(EBImage)
library(FNN)
source("misc_functions1.R")

## read in image
img.list = list.files("yaleB01")
a = 100; b = 100
ns=20

n = length(img.list)
img.mat = matrix(0, nrow=n, ncol=a*b)
for(i in 4:n){ 
	t = read.pnm(paste0("yaleB01/",img.list[[i]]))
	t = imagematrix(getChannels(t))
	t = resize(t, w=a, h=b)

	# speckle noise for sample 10
	if(i==10){
		samps = sample(1:a, 2*ns)
		t1 = t
		t[11:30, 11:30] = 0
#		t[71:90, 41:60] = 0
#		t[71:100, 41:70] = runif(100)
#		t[samps[1:ns], samps[(ns+1):(2*ns)]] = 1
	}
	img.mat[i,] = as.numeric(t)
}


img.mat0 = img.mat
img.mat0[10,] = as.numeric(t1)

#img.mat = img.mat0

k = rbfdot(1/(2*tau(img.mat)))
k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), img.mat))/.0005))
k = rbfdot(.0005)

## Change number of features and compare plots
## basis reconstruction and noise removal
## img.mat0 and features - 9 and 10 shows difference well
system.time(mod.kpc <- kpca(img.mat, features=9, kernel=k))
system.time(mod.kpcLoc <- kpcaLocantore(img.mat, features=9, kernel=k))
system.time(mod.kpcDep <- kpcaLocantore(img.mat, features=9, is.depth=TRUE, kernel=k))

# plot all results for one some image sample
i = 10
par(mfrow=c(2,2))
sm = matrix(img.mat[i,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconKwok(mod.kpc, img.mat0[i:(i+1),], k=10)
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconKwok(mod.kpcLoc, img.mat0[i:(i+1),], k=10)
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconKwok(mod.kpcDep, img.mat0[i:(i+1),], k=10)
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))
par(mfrow=c(1,1))

i = 10
par(mfrow=c(2,2))
sm = matrix(img.mat[i,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconMika(mod.kpc, img.mat[i:(i+1),])
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconMika(mod.kpcLoc, img.mat[i:(i+1),])
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconMika(mod.kpcDep, img.mat[i:(i+1),])
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))
par(mfrow=c(1,1))