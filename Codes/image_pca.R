#' Reconstructing Images using PCA
#' Examples are shown for PCA 
#' TODO: example of reconstruction using KPCA, depth-pca and other pca
#' @author: Abhirup Mallik

library(png)
lenna.orig <- readPNG("Codes/lenna.png")
res <- prcomp(lenna.orig, center = TRUE, scale = FALSE)
res.kpca <- KPCA(lenna.orig,nfeatures = 20)
pcnum <- 20 #Use first 20 principal components
lenna.trunc <- res$x[,1:pcnum] %*% t(res$rotation[,1:pcnum])
#how to recover image in kpca?
#lenna.trunc.kpc <- res.kpca$proj%*%t(res.kpca$pc)%*%res.kpca$x
if(res$scale != FALSE){
  lenna.trunc <- scale(lenna.trunc, center = FALSE , scale=1/res$scale)
}
if(res$center != FALSE){
  lenna.trunc <- scale(lenna.trunc, center = -1 * res$center, scale=FALSE)
}

r1 <- writePNG(lenna.orig,target="lenna_orig.png")
r <- writePNG(lenna.trunc,target="lenna_recon.png")
#r3 <- writePNG(lenna.trunc.kpc,target="lenna_recon_kpc.png")
