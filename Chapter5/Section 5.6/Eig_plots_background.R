sec_path = 'Rcode/Chapter5/Section 5.6/'
setwd(paste0(SLEDbook_path,sec_path))

library(sp)
library(viridis)
library(classInt)

source('pointSimSyst.R')
source('corModels.R')
source('distGeoAni.R')

#-------------------------------------------------------------------------------
#     Investigate spatial patterning of eigenvectors for systematic points
#-------------------------------------------------------------------------------

#create systematic grid of points
xy = pointSimSyst(nrow = 40, ncol = 40, lower.x.lim = 0, upper.x.lim = 1, 
    lower.y.lim = 0, upper.y.lim = 1) 
   
#create distance matrix 
dismat = distGeoAni(xy$x, xy$y, xy$x, xy$y, rotate = 90, range = .5, minorp = 1)
#create covariance matrix
covmat = corModelExponential(dismat) + diag(rep(0.000001, times = 40*40))
#eigen decomposition
eig_covmat = eigen(covmat)

#make a plot of eigenvalue, and map eigenvectors spatially
# make coordinates for image plot
x = (1:40)/40 - .05
y = x
# pick a set of 25 eigenvectors to plot
evecset = c(1, 2, 3, 4, 5,
	6, 8, 10, 12, 14, 
	20, 30, 40, 50, 60, 
	100, 200, 300, 400, 500,
	1596, 1597, 1598, 1599, 1600 )
# number of color classes
nbrks = 20
file_name = "Eigval_bg"
pdf(file = paste0(file_name,'.pdf'), width = 10, height = 10)
  oldpar = par(mar = c(5,5,1,1))
  plot(sqrt(eig_covmat$values), pch = 19, cex = 1, cex.axis = 1.5, cex.lab = 2,
  xlab = 'index', ylab = 'Square Root of Eigenvalue')
  par(oldpar)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

file_name = "Eigvec25_bg"
pdf(file = paste0(file_name,'.pdf'), width = 8, height = 8)
  oldpar = par(mar = c(2,2,3,1))
  layout(matrix(1:25, nrow = 5, byrow = TRUE))
  for(i in 1:length(evecset)) {
    z = eig_covmat$vector[,evecset[i]]*sqrt(eig_covmat$values[evecset[i]])
    z = matrix(z, nrow = 40)
    brks = quantile(z, probs = (0:nbrks)/nbrks)
    cramp = viridis(nbrks)
    image(x, y, z, breaks = brks, col = cramp, 
      main = paste0("Eigenvector ",evecset[i]),
      cex.main = 1.2, cex.axis = .8, xlab = '', ylab = '')
  #	points(xy[,1],xy[,2], pch = 19, cex = 2)
  }
  oldpar = par(mar = c(2,2,3,1))
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
