rm(list=ls())
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/KDPCA/")

library(kernlab)
source("misc_functions.R")

## function to get best gaussian kernel for data
tau = function(x){
	n = nrow(x)
	mat = matrix(0, n, n)
	for(i in 1:n){
		for(j in 1:i){
			mat[i,j] = sum((x[i,] - x[j,])^2)
			mat[j,i] = mat[i,j]
		}
	}
	mean(mat)
}

## replicate the simulations of Huang et al, Neurocomputing 2011, 74:18, 3921-3930
## Paper name: An Iterative Algorithm for Robust Kernel Principal Component Analysis

# generate data
n = 150
set.seed(20150603)
x = seq(-1,1,length=n)
y = -.3*x^2 + .1*rnorm(n)
data = cbind(x,y)

# perform 3 types of PCA
k = polydot(2)
mod.kpc = kpca(data, features=1, kernel=k)
mod.kpcLoc = kpcaLocantore(data, features=1, kernel=k)
mod.kpcDep = kpcaLocantore(data, features=1, is.depth=TRUE, kernel=k)

# standardize PCs
pc.kpc = scale(pcv(mod.kpc))
pc.kpcLoc = scale(pcv(mod.kpcLoc))
pc.kpcDep = scale(pcv(mod.kpcDep))

# dot product: measure of similarity between two vectors... more is good
abs(mean(pc.kpc * pc.kpcLoc))
abs(mean(pc.kpc * pc.kpcDep))

## Bring in contaminations now
## randomly choose samples, change their y values
nc = 10
yc = c(rnorm(nc/2)*1.5 - .5, rnorm(nc/2)*1.5 + .5)
dataf = data
dataf[sample(1:n, nc),2] = yc

kf = k
cmod.kpc = kpca(dataf, features=1, kernel=kf)
cmod.kpcLoc = kpcaLocantore(dataf, features=1, kernel=kf)
cmod.kpcDep = kpcaLocantore(dataf, features=1, is.depth=TRUE, kernel=kf)

abs(mean(pc.kpc * scale(pcv(cmod.kpc))))
abs(mean(pc.kpc * scale(pcv(cmod.kpcLoc))))
abs(mean(pc.kpc * scale(pcv(cmod.kpcDep))))
