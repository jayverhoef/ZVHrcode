sec_path = 'Rcode/Chapter3/Section 3.4/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(scatterplot3d)
library(sf)
library(spdep)



data(SO4obs)
SO4data = st_drop_geometry(SO4obs)
SO4data$xcoord = st_coordinates(SO4obs)[,1]
SO4data$ycoord = st_coordinates(SO4obs)[,2]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  so4-nn-scatterplot
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# nearest neighbor (NN) scatterplot
file_name = 'figures/so4-nn-scatterplot'

np = length(SO4data[,1])
storeResults = NULL
# loop through each record and find index of nearest neighbor by Euclidean distance
for(i in 1:np) {
  j = which((1:np)!=i)[which.min(sqrt(
    (SO4data$xcoord[i] - SO4data$xcoord[(1:np)!=i])^2 +
    (SO4data$ycoord[i] - SO4data$ycoord[(1:np)!=i])^2))]
  #i is site index, j is NN index, then values for i and j
  storeResults = rbind(storeResults,
    c(i,j,SO4data[i,'SO4'], SO4data[j,'SO4']))
}
# get absolute difference between site and nearest neighbor
d = abs(storeResults[,3] - storeResults[,4])
# create standard z-values for absolute differences
d.z = (d - mean(d))/sqrt(var(d))

# nearest neighbor scatter plot
pdf(file = paste0(file_name,'.pdf'))

  old.par = par(mar = c(5,5,1,1))
  plot(storeResults[,3:4], pch = 1, 
    xlab = 'SO4 (kg/ha) at Originating Site',
    ylab = 'SO4 (kg/ha) at Nearest Neighbor',
    cex.lab = 2, cex.axis = 1.5)
  lines(c(min(storeResults[,3]),max(storeResults[,3])),
    c(min(storeResults[,3]),max(storeResults[,3])),
    lty = 2, lwd = 2)
  points(storeResults[d.z > 3,3:4], pch = 19)
  par(old.par)
  
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  so4-outlier-sites
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# make a map from nearest neighbor scatterplots
file_name = 'figures/so4-outlier-sites'

# find which NN scatterplots have standardized z-values greater than 3
spatialOutliers = storeResults[d.z > 3,]
colnames(spatialOutliers) = c('Orig_indx','NN_indx','Orig_val','NN_val')
data(USboundary)
states = c("Virginia", "West Virginia", "North Carolina", "Pennsylvania", "Tennessee", "Maryland", "Kentucky", "New York", "Ohio", "Indiana")
# make a plot using arrows from originating site to NN with standardized z-value
# greater than 3
pdf(paste0(file_name,'.pdf'))
plot(st_geometry(USboundary[USboundary$STATE_NAME %in% states,]))
  for(i in 1:length(spatialOutliers[,1])) {
    outlier_number = i
    #plot pairs of points
    points(SO4data$xcoord[spatialOutliers[outlier_number,1]],
      SO4data$ycoord[spatialOutliers[outlier_number,1]], 
      pch = 19, col = 'red', cex = 1.2)
    points(SO4data$xcoord[spatialOutliers[outlier_number,2]],
      SO4data$ycoord[spatialOutliers[outlier_number,2]], 
      pch = 19, col = 'red', cex = 1.2)
  }
  # outliers
  points(SO4data$xcoord[146], SO4data$ycoord[146], 
		pch = 1,  cex = 2.5, lwd = 1.5)
  points(SO4data$xcoord[153], SO4data$ycoord[153], 
		pch = 1,  cex = 2.5, lwd = 1.5)
  points(SO4data$xcoord[173], SO4data$ycoord[173], 
		pch = 1,  cex = 2.5, lwd = 1.5)
  for(i in 1:length(spatialOutliers[,1])) {
    outlier_number = i
    #plot arrows connecting pairs of point, originating toward NN
    arrows(SO4data$xcoord[spatialOutliers[outlier_number,1]],
      SO4data$ycoord[spatialOutliers[outlier_number,1]],
        SO4data$xcoord[spatialOutliers[outlier_number,2]],
      SO4data$ycoord[spatialOutliers[outlier_number,2]],
      length = .09, angle = 20)
  }
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

cor(storeResults[, 3:4], method = "pearson")
# [1] 0.828628
cor(storeResults[, 3:4], method = "spearman")
# [1] 0.8447826

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  seals-nn-scatterplot
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = 'figures/seals-nn-scatterplot'
data(sealPolys)

#     We want to do Moran's I and Geary's c on observed data only

# get a list of all polygons that have observed values
ind_samp = which(!is.na(sealPolys$Estimate))
# subset the SpatialPolygonsDataFrame to those with observed values
sealPolys_obs = sealPolys[ind_samp,]
# use poly2nb to create neighborhood list for each polygon
Nlist_obs = poly2nb(sealPolys_obs, snap = 2000)
# some polygons have no neighbors, so remove them as well
sealPolys_obs = sealPolys_obs[which(unlist(
	lapply(Nlist_obs,function(x){all(x!=0)}))),]
# again use poly2nb to create neighborhood list for each polygon
Nlist_obs = poly2nb(sealPolys_obs, snap = 2000)
attr(Nlist_obs,'polyid') = as.factor(as.character(sealPolys_obs$polyid))
attr(Nlist_obs,'stockid') = as.factor(as.character(sealPolys_obs$stockid))

DF_obs = st_drop_geometry(sealPolys_obs)
storeResults = NULL
for(i in 1:length(Nlist_obs))
 storeResults = rbind(storeResults,
  c(i, DF_obs[i,'Estimate'], mean(DF_obs[Nlist_obs[[i]], 'Estimate'])))
  
# get absolute difference between site and nearest neighbor
d = abs(storeResults[,2] - storeResults[,3])
# create standard z-values for absolute differences
d.z = (d - mean(d))/sqrt(var(d))

pdf(paste0(file_name,'.pdf'))
  old.par = par(mar = c(5,5,1,1))
  plot(storeResults[,2:3], pch = 1, 
    xlab = 'Seal Trend at Originating Site',
    ylab = 'Seal Trend at Mean of Neighbors',
    cex.lab = 2, cex.axis = 1.5, xlim = c(-.6,.85), ylim = c(-.6,.85))
  lines(c(min(storeResults[,2]),max(storeResults[,2])),
    c(min(storeResults[,2]),max(storeResults[,2])),
    lty = 2, lwd = 2)
  points(storeResults[d.z > 3,2:3], pch = 19)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
