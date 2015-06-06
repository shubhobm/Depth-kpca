rm(list=ls())
setwd("C:/Study/My projects/Depth-kpca/Codes")
source("misc_functions.R")
library(kernlab)

X = matrix(rnorm(400), ncol=2)
normz = rowSums(X^2)
label = ifelse(normz > .7, 1, 0)
plot(X, pch=19, col=ifelse(label==0, "red", "blue"))

par(mfrow=c(1,2))
mod.kpca = kpca(X, features=2, kernel="rbfdot",
            kpar=list(sigma=100))
plot(rotated(mod.kpca), pch=19, col=ifelse(label==0, "red", "blue"))
mod.kpcaLoc = kpcaLocantore(X, features=2, kernel="rbfdot",
            kpar=list(sigma=100))
plot(rotated(mod.kpcaLoc), pch=19, col=ifelse(label==0, "red", "blue"))
par(mfrow=c(1,1))

######### Sample comparison
data(iris)
test <- sample(1:150,20)

system.time(kpcLoc <- kpcaLocantore(~.,data=iris[-test,-5],kernel="rbfdot",
            kpar=list(sigma=2),features=2)
)
system.time(kpcDep <- kpcaLocantore(~.,data=iris[-test,-5],kernel="rbfdot",
            kpar=list(sigma=2),features=2, is.depth=TRUE)
)
system.time(kpc <- kpca(~.,data=iris[-test,-5],kernel="rbfdot",
            kpar=list(sigma=2),features=2)
)

#plot the data projection on the components
par(mfrow=c(1,3))
plot(rotated(kpc),col=as.integer(iris[-test,5]), pch=19,
     xlab="1st Principal Component",ylab="2nd Principal Component")
plot(rotated(kpcLoc),col=as.integer(iris[-test,5]), pch=19,
     xlab="1st Principal Component",ylab="2nd Principal Component")
plot(rotated(kpcDep),col=as.integer(iris[-test,5]), pch=19,
     xlab="1st Principal Component",ylab="2nd Principal Component")
par(mfrow=c(1,1))

# Bivariate normal mixture
set.seed(120214)
sig = matrix(c(1,.9,.9,1), nrow=2)
sig2 = matrix(c(1,-.9,-.9,1), nrow=2)
X1 = my.mvrnorm(500, mu=c(-2,-2), Sigma=sig2)
X3 = my.mvrnorm(500, mu=c(2,-2), Sigma=sig)
X = rbind(X1,X3)
plot(X)

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

system.time(d <- kdepth.SP(xygrid, X))

