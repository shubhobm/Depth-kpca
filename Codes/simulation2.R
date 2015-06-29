#' Simulation Experiment 2
#' Framework of Yang et al.
#' @author Abhirup Mallik
#' 


angdist <- function(x,y){
  theta <- crossprod(x,y)/(sqrt(crossprod(x)*crossprod(y)))
  acos(theta)
}

vecdist <- function(x,y){
  sqrt(crossprod(x,y))
}

set.seed(12345)
n <- 101
x <- runif(n)
r <- rnorm(n)
y <- 2*sin(2*pi*x)+0.5+0.1*r
yout <- y
#number of outliers
numout <- 30
outliers <- rnorm(numout,mean = 5,sd=1)
#randomly replace points in x by outliers
indexout <- sample.int(n,size=numout)
yout[indexout] <- outliers
data <- cbind(x,yout)
#Distance between eigen vectors of KPCA and other PCA

#Cosine distance:
#Kernel Choise
k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), data))/.0005))
#KPCA
pcnum <- 1
mod.kpc <- kpca(data, features=pcnum, kernel=k)
mod.kpcLoc <- kpcaLocantore(data, features=pcnum, kernel=k)
mod.kpcDep <- kpcaLocantore(data, features=pcnum, is.depth=TRUE, kernel=k)

## Study for increasing number of outliers from 0 to 50

theta.Loc <- numeric(51)
theta.Dep <- numeric(51)
dist.Loc <- numeric(51)
dist.Dep <- numeric(51)

for(i in 1:51){
  numout <- i -1
  outliers <- rnorm(numout,mean = 5,sd=1)
  #randomly replace points in x by outliers
  indexout <- sample.int(n,size=numout)
  yout[indexout] <- outliers
  data <- cbind(x,yout)
  #Distance between eigen vectors of KPCA and other PCA
  
  #Cosine distance:
  #Kernel Choise
  k = rbfdot(- 1/mean(log(kernelMatrix(rbfdot(.0005), data))/.0005))
  #KPCA
  pcnum <- 1
  mod.kpc <- kpca(data, features=pcnum, kernel=k)
  mod.kpcLoc <- kpcaLocantore(data, features=pcnum, kernel=k)
  mod.kpcDep <- kpcaLocantore(data, features=pcnum, is.depth=TRUE, kernel=k)
  theta.Loc[i] <- angdist(pcv(mod.kpc),mod.kpcLoc$pcv)
  theta.Dep[i] <- angdist(pcv(mod.kpc),mod.kpcDep$pcv)
  dist.Loc[i] <- vecdist(pcv(mod.kpc),mod.kpcLoc$pcv)
  dist.Dep[i] <- vecdist(pcv(mod.kpc),mod.kpcDep$pcv)
}

dat_theta <- data.frame(num = out_num[1:20], dep = theta.Dep[1:20], loc = theta.Loc[1:20])
dat_theta_melt <- melt(dat_theta, id = "num")
ggplot(data=dat_theta_melt,
aes(x=num, y=value, colour=variable)) +
geom_line() + labs(title = "Anglular Similarity") + xlab("number of outliers") + ylab("Angle")

dat_dist <- data.frame(num = out_num[1:20], dep = dist.Dep[1:20], loc=dist.Loc[1:20])
dat_dist_melt <- melt(dat_dist, id = "num")
ggplot(data=dat_dist_melt,
       aes(x=num, y=value, colour=variable)) +
  geom_line()+ labs(title = "Distance") + xlab("number of outliers") + ylab("Distance")


require(tikzDevice)
require(gridExtra)
plot1 <- ggplot(data=dat_theta_melt,
aes(x=num, y=value, colour=variable)) +
geom_line() + labs(title = "Anglular Similarity") + xlab("number of outliers") + ylab("Angle")
plot2 <- ggplot(data=dat_dist_melt,
aes(x=num, y=value, colour=variable)) +
geom_line()+ labs(title = "Distance") + xlab("number of outliers") + ylab("Distance")

tikz(file = "simulation2plot.tex",width = 7, height = 4, standAlone = TRUE)
grid.arrange(plot1, plot2, ncol=2)
dev.off()

