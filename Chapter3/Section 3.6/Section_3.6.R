sec_path = 'Rcode/Chapter3/Section 3.6/figures/'
setwd(paste0(SLEDbook_path,sec_path))

################################################################################
#-------------------------------------------------------------------------------
#                 Nearest Neighbor Correlation for SO4 data
#-------------------------------------------------------------------------------
################################################################################
# Pearson's and Spearman's correlations
library(ZVHdata)
data(SO4obs)

np = length(SO4obs@coords[,1])
storeResults = NULL
# loop through each record and find index of nearest neighbor by Euclidean distance
for(i in 1:np) {
  j = which((1:np)!=i)[which.min(sqrt(
    (SO4obs@coords[i,1] - SO4obs@coords[(1:np)!=i,1])^2 +
    (SO4obs@coords[i,2] - SO4obs@coords[(1:np)!=i,2])^2))]
  #i is site index, j is NN index, then values for i and j
  storeResults = rbind(storeResults,
    c(i,j,SO4obs@data[i,'SO4'], SO4obs@data[j,'SO4']))
}
cor(storeResults[,3:4])[1,2]
cor(storeResults[,3:4], method = 'spearman')[1,2]

################################################################################
#-------------------------------------------------------------------------------
#                            so4-acf-svgm-3x2-clean
#-------------------------------------------------------------------------------
################################################################################
# semivariograms and acf plots (Figure 3.7)
# Need to run script Section_3.4.R to generate SpatialOutliers object
# uncomment below and source as necessary
# source(paste0(SLEDbook_path,'Rcode/Chapter3/Section 3.4/Section_3.4.R'))
file_name = "so4-acf-svgm-3x2-clean"
spatialOutliers
# we declare sites with indexes 153, 146, and 173 as spatial outliers
# just to make sure, here is a plot
data(USboundary)
states = c("Virginia", "West Virginia", "North Carolina", "Pennsylvania", "Tennessee", "Maryland", "Kentucky", "New York", "Ohio", "Indiana")
# make a plot using arrows from originating site to NN with standardized z-value
# greater than 3
plot(USboundary[USboundary$STATE_NAME %in% states,])
plot(SO4obs[c(146,153,173),], add = TRUE, pch = 19, cex = 2, col = 'red')
# remove the outlers
SO4clean = SO4obs[!(1:length(SO4obs) %in% c(146,153,173)),]
#use distance function, which is lower triangular, for all distances
# turn into a vector
ds = as.vector(dist(SO4clean@coords))
# get all mean-corrected cross-products, subset to lower triangular, then vector
cp = outer(sqrt(SO4clean@data$SO4) - mean(sqrt(SO4clean@data$SO4)),
  sqrt(SO4clean@data$SO4) - mean(sqrt(SO4clean@data$SO4)))[
    lower.tri(outer(1:length(SO4clean),1:length(SO4clean)))]
# get all 0.5*squared differences, subset to lower triangular, then vector
sv = (0.5*(outer(sqrt(SO4clean@data$SO4), rep(1, times = length(SO4clean))) -
  outer(rep(1, times = length(SO4clean)), sqrt(SO4clean@data$SO4)))^2)[
    lower.tri(outer(1:length(SO4clean),1:length(SO4clean)))]
#2500 kilometers is approximately 1/2 maximum distance
# truncate vectors to shorter distances
ind = ds/1000 < 2500
ds = ds[ind]/1000
cp = cp[ind]
sv = sv[ind]
## Make the data into 15 equally spaced bins. 
bins <- cut(ds,  breaks = seq(0, 2500, by = 2500/15))
bins <- factor(bins, levels = levels(bins), 
  labels = format(round(tapply(ds, bins, mean), 0), 0))
# Compute covariogram and semivariograms using bin means
sample_cvg <- tapply(cp, bins, mean)
sample_svg <- tapply(sv, bins, mean)
sample_n <- tapply(cp, bins, length)

# make plots
mtxtadj = -.12
cex_lab = 2.5
cex_axis = 1.8
pdf(paste0(file_name,'.pdf'), width = 11, height = 14)
  layout(matrix(1:6, nrow = 3, byrow = TRUE))
  old.par = par(mar = c(5,5,5,2))
  plot(ds, cp, pch = 19, cex = .5, 
    ylab = 'Autocovariance', xlab = 'Distance (km)', cex.lab = cex_lab, 
    cex.axis = cex_axis)
  mtext('A', adj = mtxtadj, cex = 3)
  plot(ds, sv, pch = 19, cex = .5, 
    ylab = 'Semivariogram', xlab = 'Distance (km)', cex.lab = cex_lab, 
    cex.axis = cex_axis)
  mtext('B', adj = mtxtadj, cex = 3)
  boxplot(cp ~ bins, ylab = 'Autocovariance', xlab = 'Distance (km)',
    cex.lab = cex_lab, cex.axis = cex_axis)
  mtext('C', adj = mtxtadj, cex = 3)
  boxplot(sv ~ bins, ylab = 'Semivariogram', xlab = 'Distance (km)',
    cex.lab = cex_lab, cex.axis = cex_axis)
  mtext('D', adj = mtxtadj, cex = 3)
  plot(as.numeric(names(sample_cvg)), sample_cvg, pch = 19, 
    ylab = 'Autocovariance', xlab = 'Distance (km)',
    cex.lab = cex_lab, cex.axis = cex_axis, 
    cex = 8*sqrt(sample_n)/max(sqrt(sample_n)))
  mtext('E', adj = mtxtadj, cex = 3)
  plot(as.numeric(names(sample_svg)), sample_svg, pch = 19,
    ylab = 'Semivariogram', xlab = 'Distance (km)',
    cex.lab = cex_lab, cex.axis = cex_axis, 
    cex = 8*sqrt(sample_n)/max(sqrt(sample_n)))
  mtext('F', adj = mtxtadj, cex = 3)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
system(paste0('convert ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\' ','-flatten ','\'',SLEDbook_path,
  sec_path,file_name,'.png','\''))

################################################################################
#-------------------------------------------------------------------------------
#                            so4-omnisvgm-quad-resids
#-------------------------------------------------------------------------------
################################################################################
# semivariograms and acf plots on residuals from quadratic surface (Figure 3.8)
file_name = "so4-omnisvgm-quad-resids"

dat = data.frame(SO4sqrt = sqrt(SO4clean@data$SO4),
  x = SO4clean@coords[,1]/1000, y = SO4clean@coords[,2]/1000)
res = resid(lm(SO4sqrt ~ poly(x, y, degree = 2, raw = TRUE), data = dat))
# get all mean-corrected cross-products, subset to lower triangular, then vector
cpr = outer(res - mean(res), res - mean(res))[
    lower.tri(outer(1:length(SO4clean),1:length(SO4clean)))]
# get all 0.5*squared differences, subset to lower triangular, then vector
svr = (0.5*(outer(res, rep(1, times = length(SO4clean))) -
  outer(rep(1, times = length(SO4clean)), res))^2)[
    lower.tri(outer(1:length(SO4clean),1:length(SO4clean)))]
cpr = cpr[ind]
svr = svr[ind]

# Compute covariogram and semivariograms using bin means
sample_cvgr <- tapply(cpr, bins, mean)
sample_svgr <- tapply(svr, bins, mean)

# make plots
mtxtadj = -.16
pdf(paste0(file_name,'.pdf'), width = 14, height = 7)
  layout(matrix(1:2, nrow = 1, byrow = TRUE))
  old.par = par(mar = c(5,5,5,2))
  plot(as.numeric(names(sample_cvgr)), sample_cvgr, pch = 19,
    ylab = 'Autocovariance', xlab = 'Distance (km)',
    cex.lab = 2, cex.axis = 1.5, cex = 5*sample_n/max(sample_n))
  mtext('A', adj = mtxtadj, cex = 3)
  plot(as.numeric(names(sample_svgr)), sample_svgr, pch = 19, 
    ylab = 'Semivariogram', xlab = 'Distance (km)', ylim = c(0,0.8),
    cex.lab = 2, cex.axis = 1.5, cex = 5*sample_n/max(sample_n))
  mtext('B', adj = mtxtadj, cex = 3)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

################################################################################
#-------------------------------------------------------------------------------
#                            polar-partition
#-------------------------------------------------------------------------------
################################################################################
# Create a polar partition for binning semivariograms (Figure 3.9)
file_name = "polar-partition"

library(shape)
pdf(paste0(file_name,'.pdf'))
  old.par = par(mar = c(5,5,1,1))
  xpts <- c(-4,4)
  ypts <- c(-4,4)
  par(pty="s")
  plot(xpts,ypts,type="n",xlab="x-coordinate of lag",ylab="y-coordinate of lag",
    cex.axis = 1.5, cex.lab = 2)
  plotcircle(r = 1, mid = c(0,0), from = -pi/8, to = 7*pi/8, lwd = 3)
  plotcircle(r = 2, mid = c(0,0), from = -pi/8, to = 7*pi/8, lwd = 3)
  plotcircle(r = 3, mid = c(0,0), from = -pi/8, to = 7*pi/8, lwd = 3)
  plotcircle(r = 4, mid = c(0,0), from = -pi/8, to = 7*pi/8, lwd = 3)
  lines(c(0,3.695518), c(0,sqrt(16-3.695518^2)), lwd = 3)
  lines(c(0,1.530734), c(0,sqrt(16-1.530734^2)), lwd = 3)
  lines(c(0,-1.530734), c(0,sqrt(16-1.530734^2)), lwd = 3)
  lines(c(-3.695518,3.695518), c(sqrt(16-3.695518^2), -sqrt(16-3.695518^2)),
    lwd = 3)
  lines(c(0,0),c(-4.2,4.2), lty = 2, lwd = 2)
  lines(c(-4.2,4.2),c(0,0), lty = 2, lwd = 2)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

################################################################################
#-------------------------------------------------------------------------------
#                            so4-dirsvgm-quad-resids
#-------------------------------------------------------------------------------
################################################################################
# directional semivariogram on residuals from quadratic surface (Figure 3.10)
file_name = "so4-dirsvgm-quad-resids"
library(gstat)
vgm_resid_dir <- variogram(res ~ 1, loc=~x+y, data=data.frame(res = res, 
  x = dat$x, y = dat$y), alpha=c(0,45,90,135), cutoff = 2500, width = 2500/15)
 
pdf(paste0(file_name,'.pdf'), height = 8, width = 11)
  old_par = par(mar = c(5,5,3,1))
  plot(vgm_resid_dir$dist, vgm_resid_dir$gamma, type = 'n', xlab = 'Distance (km)',
    ylab = 'Semivariogram', cex.lab = 2, cex.axis = 1.5)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], pch = 1, 
    cex = 10*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], pch = 2, 
    cex = 10*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], pch = 5, 
    cex = 10*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], pch = 22, 
    cex = 10*vgm_resid_dir$np/max(vgm_resid_dir$np))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], lwd = 2)
  legend(100, 1.1, legend=(c('N-S','NE-SW','E-W','SE-NW')), 
    pch = c(1,2,5,22), cex = 2)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

################################################################################
#-------------------------------------------------------------------------------
#                            so4-smoothed-svgm
#-------------------------------------------------------------------------------
################################################################################
# directional semivariogram on residuals from quadratic surface (Figure 3.10)
file_name = "so4-smoothed-svgm"

# let theta be vector of 3 parameters such that: exp(theta[1]) = nugget, 
# exp(theta[2]) = partial sill, exp(theta[3]) = range.
# Exponential semivariogram function
svgm_exp = function(theta, x)
{
  exp(theta[1]) + exp(theta[2])*(1 - exp(-x/exp(theta[3])))
}

#creata a function to evaluation Cressie's weighted least squares for
# a given semivariogram function and sample semivariogram
cwls = function(theta, svgm_fun, dists, sample_svgm, sample_n)
{
  sum(sample_n*(sample_svgm - svgm_fun(theta, dists))^2/svgm_fun(theta, dists)^2)
}  
# minimize cwls using optim
opt_out = optim(log(c(.2, .6, 500)), cwls, svgm_fun = svgm_exp, dists = dists, 
  sample_svgm = sample_svgr, sample_n = sample_n)
# nugget
exp(opt_out$par[1])
# partial sill
exp(opt_out$par[2])
# range
exp(opt_out$par[3])
# plot the sample semivariogram and it fit by Cressie's weighted least squares
pdf(paste0(file_name,'.pdf'), height = 8, width = 11)
  old.par = par(mar = c(5,5,1,1))
  plot(as.numeric(names(sample_svgr)), sample_svgr, pch = 19, 
    ylab = 'Semivariogram', xlab = 'Distance (km)', ylim = c(0,0.8),
    cex.lab = 2, cex.axis = 1.5, cex = 5*sample_n/max(sample_n),
    xlim = c(0,2500))
  lines(1:2500, svgm_exp(opt_out$par, 1:2500), lwd = 3)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
