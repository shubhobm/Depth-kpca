### Image reconstruction using kpca
rm(list=ls())
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/KDPCA")

library(kernlab)
library(png)
library(ripa)
library(EBImage)
source("misc_functions1.R")

## read in images
img.list = list.files("images")
n = 93
a=480; b=512

s=sample(1:n, 30)

s=1:50
img.mat = matrix(0, nrow=length(s), ncol=a*b)
for(i in 1:length(s)){ 
	t = readPNG(paste0("images/",img.list[[s[i]]]))
	t = resize(t, w=a)
	t = imagematrix(t)

	# speckle noise for sample 10
	if(i==10){
		samps = c(sample(1:nrow(t), 50), sample(1:ncol(t), 50))
		t1 = t
		t[101:150, 101:150] = runif(2500)
		t[321:370, 221:270] = runif(2500)	
	}

	img.mat[i,] = as.numeric(t)
}

# store normal data
img.mat0 = img.mat
img.mat0[10,] = as.numeric(t1)

#img.mat = img.mat0

k = rbfdot(.0005)
system.time(mod.kpc <- kpca(img.mat, features=3, kernel=k))
system.time(mod.kpcLoc <- kpcaLocantore(img.mat, features=3, kernel=k))
system.time(mod.kpcDep <- kpcaLocantore(img.mat, features=3, is.depth=TRUE, kernel=k))

# plot all results for one some image sample
i = 10
par(mfrow=c(2,2))
sm = matrix(img.mat0[i,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconGood1(mod.kpc, img.mat0[i:(i+1),], k=5)
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconGood1(mod.kpcLoc, img.mat0[i:(i+1),], k=5)
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))

s = reconGood1(mod.kpcDep, img.mat0[i:(i+1),], k=5)
sm = matrix(s[1,], nrow=a, byrow=F)
plot(imagematrix(sm))
par(mfrow=c(1,1))