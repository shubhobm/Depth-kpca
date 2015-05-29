rm(list=ls())
setwd("C:/Study/My projects/Depth-kpca/Codes")
source("misc_functions.R")

# Bivariate normal mixture
set.seed(120214)
sig = matrix(c(1,.9,.9,1), nrow=2)
sig2 = matrix(c(1,-.9,-.9,1), nrow=2)
X1 = my.mvrnorm(500, mu=c(-2,2), Sigma=sig)
X3 = my.mvrnorm(500, mu=c(2,-2), Sigma=sig)
X = rbind(X1,X3)

# make grid of points
mingrid=-5
maxgrid=5
res=.2

pts = seq(mingrid, maxgrid, by=res)
lengrid = length(pts)
xcoord = rep(pts, rep(lengrid,lengrid))
ycoord = rep(pts, lengrid)
xygrid = cbind(xcoord,ycoord)
rm(xcoord,ycoord)

d = kdepth.SP(X, X, 1)
