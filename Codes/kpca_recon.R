### Image reconstruction using kpca
rm(list=ls())
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/KDPCA")

library(kernlab)
library(png)
library(ripa)
source("misc_functions1.R")

## read in images
img.list = list.files("images")
n = 93
img.mat = matrix(0, nrow=n, ncol=480*512)

for(i in 1:n){ 
	t = readPNG(paste0("images/",img.list[[i]]))
	t = imagematrix(t)

	# speckle noise for sample 10
	if(i==10){
		samps = c(sample(1:nrow(t), 50), sample(1:ncol(t), 50))
		t1 = t
		t[201:250, 201:250] = runif(2500)	
	}

	img.mat[i,] = as.numeric(t)
}

# store normal data
img.mat0 = img.mat
img.mat0[10,] = as.numeric(t1)

#img.mat = img.mat0
#k = rbfdot(1/2*ta)
system.time(mod.kpc <- kpca(img.mat, features=10))
system.time(mod.kpcLoc <- kpcaLocantore(img.mat, features=10))
system.time(mod.kpcDep <- kpcaLocantore(img.mat, features=10, is.depth=TRUE))

# plot all results for one some image sample
i = 10
par(mfrow=c(2,2))
sm = matrix(img.mat[i,], nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconGood(mod.kpc, img.mat0[i:(i+1),], k=20)
sm = matrix(s[1,], nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconGood(mod.kpcLoc, img.mat0[i:(i+1),], k=20)
sm = matrix(s[1,], nrow=480, byrow=F)
plot(imagematrix(sm))

s = reconGood(mod.kpcDep, img.mat0[i:(i+1),], k=20)
sm = matrix(s[1,], nrow=480, byrow=F)
plot(imagematrix(sm))
par(mfrow=c(1,1))