# Depth-kpca

Contains codes, background papers and manuscript for depth-based robust KERNEL principal component estimation.

## A fast and accurate denoising technique for high-resolution image streams using data depth-based robust kernel PCA

### Abstract
Classical kernel principal component analysis (PCA) has been shown to be quite sensitive to the presence of outliers, so many robust versions of kernel PCA exists in the literature. However, they suffer from lack of scalability with increasing amount of data and less accuracy compared to the classical version. Here we propose to use multivariate ranks obtained using data depth calculated in kernel spaces to get a robust version of kernel PCA applicable in big data scenarios, for example, denoising of large images. Given the kernel matrix, the extra calculation required here on top of classical kernel PCA does not grow with number of features, hence being more scalable. We provide theoretical properties related to the method, as well as demonstrate its effectiveness through simulations and a real data application. We also propose a scheme to outline how to use our method to check a continuous stream of images for noise, and denoise them if necessary.

### Owner
Subho Majumdar <majum010@umn.edu>

Abhirup Mallik <malli066@umn.edu>
