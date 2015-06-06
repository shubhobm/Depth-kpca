<<<<<<< HEAD
rm(list=ls())
=======
Enter file contents hererm(list=ls())
>>>>>>> 5746fd20528722df13ac0ddb999a60c88c029c27
setwd("\\\\dfs.com/root/Dept-Decision/Dept-Users/Majumdar/KDPCA/")

library(kernlab)
source("misc_functions.R")

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

