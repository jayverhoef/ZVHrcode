sec_path = 'Rcode/Chapter4/Section 4.5/'
setwd(paste0(SLEDbook_path,sec_path))

# attach data library
library(ZVHdata)
library(sf)
library(akima)
library(viridis)
library(classInt)
library(colorspace)
library(gstat)
source('addBreakColorLegend.R')

# load data for graphics and analysis
data(SO4obs)

SO4_data = unlist(as.vector(st_drop_geometry(SO4obs)))

# from Section 3.6, remove the outlers and use sqrt of response
SO4clean = SO4_data[!(1:length(SO4_data) %in% c(146,153,173))]
xy = st_coordinates(SO4obs)
xy = xy[!(1:length(SO4_data) %in% c(146,153,173)),]
DF = data.frame(z = sqrt(SO4clean), easting = xy[,1], northing = xy[,2])

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        Wire-frame plots of planar, quadratic, and cubic surfaces
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/so4-poly-surfaces-1x3"

# Planar surface
# fit the linear model on the coordinates
fit_1 <- lm(z ~ poly(easting, northing, degree = 1), data = DF)
# create new data.frame with fitted values
data <- data.frame(u = DF$easting , v = DF$northing, z = fit_1$fitted.values)
# interpolate the surface on a 40 x 40 grid
surf_1 <- with(data, interp(u,v,z))
 
# Quadratic surface
# fit the linear model on the coordinates
fit_2 <- lm(z ~ poly(easting, northing, degree = 2), data = DF)
# create new data.frame with fitted values
data <- data.frame(u = DF$easting , v = DF$northing, z = fit_2$fitted.values)
# interpolate the surface on a 40 x 40 grid
surf_2 <- with(data, interp(u,v,z))

# Cubic surface
# fit the linear model on the coordinates
fit_3 <- lm(z ~ poly(easting, northing, degree = 3), data = DF)
# create new data.frame with fitted values
data <- data.frame(u = DF$easting , v = DF$northing, z = fit_3$fitted.values)
# interpolate the surface on a 40 x 40 grid
surf_3 <- with(data, interp(u,v,z))

pdf(paste0(file_name,'.pdf'), height = 5, width = 15)
  layout(matrix(1:3, nrow = 1))
  persp(surf_1$x, surf_1$y, surf_1$z*1000000, phi = 25, theta = 315, expand = .4, 
    scale = FALSE, xlab = 'Easting', ylab = 'Northing', zlab = 'SO4',
    cex.lab = 1.8)
  text(-.37, .2, 'A              Planar', cex = 3, pos = 4)
  persp(surf_2$x, surf_2$y, surf_2$z*1000000, phi = 25, theta = 315, expand = .4, 
    scale = FALSE, xlab = 'Easting', ylab = 'Northing', zlab = 'SO4',
    cex.lab = 1.8)
  text(-.37, .2, 'B             Quadratic', cex = 3, pos = 4)
  persp(surf_3$x, surf_3$y, surf_3$z*1000000, phi = 25, theta = 315, expand = .4, 
    scale = FALSE, xlab = 'Easting', ylab = 'Northing', zlab = 'SO4',
    cex.lab = 1.8)
  text(-.37, .2, 'C              Cubic', cex = 3, pos = 4)
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
#        Leverage quadratic surface
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# standardize variables for stability of inverse
xc = (DF$easting - mean(DF$easting))/sqrt(var(DF$easting))
yc = (DF$northing - mean(DF$northing))/sqrt(var(DF$northing))

# quadratic surface
# create design matrix
X = model.matrix(~ xc + yc + I(xc^2) + I(yc^2) + xc*yc, data = DF)
# compute leverage
levrg = diag(X %*% solve(t(X) %*% X, t(X)))
# make a boxplot of leverage values
boxplot(levrg)
# find out which site have highest leverage
data(USboundary)
plot(st_geometry(USboundary))
points(DF[which(levrg >.15),c('easting','northing')], pch = 19, 
	cex = 2, col = 'red')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        Leverage Cubic Surface
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# create design matrix
X = model.matrix(~ xc + yc + I(xc^2) + I(yc^2) + xc*yc + I(xc^3) + I(yc^3) + 
  I(xc^2)*yc + xc*I(yc^2), data = DF)
# compute leverage
levrg = diag(X %*% solve(t(X) %*% X, t(X)))
# make a boxplot of leverage values
boxplot(levrg)
# find out which site have highest leverage
data(USboundary)
plot(st_geometry(USboundary))
points(DF[which(levrg >.4),c('easting','northing')], add = TRUE, 
	pch = 19, cex = 2, col = 'red')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        F-tests for order of polynomial
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

fit_4 <- lm(z ~ poly(easting, northing, degree = 4), data = DF)
fit_5 <- lm(z ~ poly(easting, northing, degree = 5), data = DF)

AOV_1_2 = anova(fit_1, fit_2)
AOV_2_3 = anova(fit_2, fit_3)
AOV_3_4 = anova(fit_3, fit_4)
AOV_4_5 = anova(fit_4, fit_5)

tab = data.frame(Surface = c('Planar', 'Quadratic', 'Cubic', 'Quartic',
    'Quintic'), 
  R2 = c(summary(fit_1)$r.squared, summary(fit_2)$r.squared, 
    summary(fit_3)$r.squared, summary(fit_4)$r.squared, 
    summary(fit_5)$r.squared),
  ResDF = c(summary(fit_1)$df[2], summary(fit_2)$df[2], summary(fit_3)$df[2],
    summary(fit_4)$df[2], summary(fit_5)$df[2]),
  F = c(NA, AOV_1_2$F[2], AOV_2_3$F[2], AOV_3_4$F[2], AOV_4_5$F[2]),
  Pvalue = c(NA, AOV_1_2[['Pr(>F)']][2], AOV_2_3[['Pr(>F)']][2], 
    AOV_3_4[['Pr(>F)']][2], AOV_4_5[['Pr(>F)']][2])
)
tab

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        Plot Residuals Along Coordinates
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# plot residuals along coordinates
r2 = residuals(fit_2)
plot(DF$northing, r2)
plot(DF$easting, r2)

# plot residuals along coordinates
r3 = residuals(fit_3)
plot(DF$northing, r3)
plot(DF$easting, r3)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        Plot Residuals Spatially
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# plot residuals spatially from quadratic surface
cip = classIntervals(r2, n = 6, style = 'fisher')
palp = viridis(6)
cip_colors = findColours(cip, palp)
  old.par = par(mar = c(0,0,5,0))
  layout(matrix(1:2, nrow = 1, byrow = TRUE), widths = c(3,1))
  plot(st_geometry(USboundary), border = 'black')
  plot(SO4obs, col = cip_colors, add = TRUE, pch = 19, cex = 1.5)
  par(mar = c(0,0,0,0))
  plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
  addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
    breaks = cip$brks, colors = palp, cex = 1.5)
  par(old.par)

# plot residuals spatially from cubic surface
ip = classIntervals(r3, n = 6, style = 'fisher')
palp = viridis(6)
cip_colors = findColours(cip, palp)
source('addBreakColorLegend.R')
  old.par = par(mar = c(0,0,5,0))
  layout(matrix(1:2, nrow = 1, byrow = TRUE), widths = c(3,1))
  plot(st_geometry(USboundary), border = 'black')
  plot(SO4obs, col = cip_colors, add = TRUE, pch = 19, cex = 1.5)
  par(mar = c(0,0,0,0))
  plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
  addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
    breaks = cip$brks, colors = palp, cex = 1.5)
  par(old.par)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        Directional Semivariogram for Cubic Residuals
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/so4-dirsvgm-cubic-resids"
vgm_resid_dir <- variogram(res ~ 1, loc=~x+y, data=data.frame(res = r3, 
  x = DF$easting/1000, y = DF$northing/1000), alpha=c(0,45,90,135), 
  cutoff = 2500, width = 2500/15)
 
pdf(paste0(file_name,'.pdf'), height = 8, width = 11)
  cexA = 4
  old_par = par(mar = c(5,5,3,1))
  plot(vgm_resid_dir$dist, vgm_resid_dir$gamma, type = 'n', 
		xlab = 'Distance (km)', ylab = 'Semivariogram', cex.lab = 2, 
		cex.axis = 1.5)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], pch = 1, 
    cex = cexA*sqrt(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'np'])/
			max(sqrt(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'np'])))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], pch = 2, 
    cex = cexA*sqrt(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'np'])/
			max(sqrt(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'np'])))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], pch = 5, 
    cex = cexA*sqrt(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'np'])/
			max(sqrt(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'np'])))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], pch = 22, 
    cex = cexA*sqrt(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'np'])/
			max(sqrt(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'np'])))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], lwd = 2)
  legend(100, 0.52, legend=(c('N-S','NE-SW','E-W','SE-NW')), 
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

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        Spatial Tests
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# spatial test for isotropy
#library(spTest)
# the spTest is archived.  We have taken these functions from spTest and
# included them as source files.  The spTest function was created by
# Zachary Weller, email = wellerz@stat.colostate.edu
source('GuanTestGrid.R')
source('GuanTestUnif.R')
source('GuanTestGrid_help.R')
source('GuanTestUnif_help.R')
source('MaityTest.R')
source('MaityTest_help.R')
source('htestIso_class.R')
# quadratic
# transform the coordinates by scaling them to be much smaller, 
# and translate them so they have lowest value at zero
dt_res<- cbind(DF[,c('easting','northing')]/100000, r2)
dt_res[,1] = dt_res[,1] - min(dt_res[,1])
dt_res[,2] = dt_res[,2] - min(dt_res[,2])

min(dt_res[,1])
max(dt_res[,1])
min(dt_res[,2])
max(dt_res[,2])


GuanTestUnif(as.matrix(dt_res), 
  lagmat = rbind(c(1, 0), c(0, 1), c(1,1),c(-1, 1)), 
  A = rbind(c(1, -1, 0, 0), c(0, 0, 1, -1)), df = 2, 
  xlims = c(0, 45), 
  ylims = c(0, 30), 
  grid.spacing = c(1,1),
  window.dims = c(15,15))

set.seed(10001)
MaityTest(as.matrix(dt_res), 
  lagmat = rbind(c(1, 0), c(0, 1), c(1, 1),c(-1, 1)), 
  A = rbind(c(1, -1, 0, 0), c(0, 0, 1, -1)), df = 2, 
  xlims = c(0, 45), 
  ylims = c(0, 30), 
  block.dims = c(15,10))

# cubic
# transform the coordinates by scaling them to be much smaller, 
# and translate them so they have lowest value at zero
dt_res<- cbind(DF[,c('easting','northing')]/100000, r3)
dt_res[,1] = dt_res[,1] - min(dt_res[,1])
dt_res[,2] = dt_res[,2] - min(dt_res[,2])

min(dt_res[,1])
max(dt_res[,1])
min(dt_res[,2])
max(dt_res[,2])

GuanTestUnif(as.matrix(dt_res), 
  lagmat = rbind(c(1, 0), c(0, 1), c(1,1),c(-1, 1)), 
  A = rbind(c(1, -1, 0, 0), c(0, 0, 1, -1)), df = 2, 
  xlims = c(0, 45), 
  ylims = c(0, 30), 
  grid.spacing = c(1,1),
  window.dims = c(15,15))

set.seed(10002)
MaityTest(as.matrix(dt_res), 
  lagmat = rbind(c(1, 0), c(0, 1), c(1, 1),c(-1, 1)), 
  A = rbind(c(1, -1, 0, 0), c(0, 0, 1, -1)), df = 2, 
  xlims = c(0, 45), 
  ylims = c(0, 30), 
  block.dims = c(15,10))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#    Omnidirectional Semivariogram  and Autocovariance for Cubic Residuals
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

vgm_resid_omni <- variogram(res ~ 1, loc=~x+y, data=data.frame(res = r3, 
  x = DF$easting/1000, y = DF$northing/1000), cutoff = 2500, 
  width = 2500/15)

cvgm_resid_omni <- variogram(res ~ 1, loc=~x+y, data=data.frame(res = r3, 
  x = DF$easting/1000, y = DF$northing/1000), cutoff = 2500, 
  width = 2500/15, covariogram = TRUE)

file_name = "figures/so4-omnisvgm-cubic-resids"
pdf(paste0(file_name,'.pdf'), height = 7, width = 15)
  layout(matrix(1:2, nrow = 1))
    old_par = par(mar = c(5,5,5,1))
    plot(cvgm_resid_omni$dist, cvgm_resid_omni$gamma, xlab = 'Distance (km)',
      ylab = 'Autocovariance Function', cex.lab = 2, cex.axis = 1.5, pch = 19,
      cex = 5*sqrt(cvgm_resid_omni$np)/max(sqrt(cvgm_resid_omni$np)))
    mtext('A', adj = -.15, padj = -.4, cex = 3)
    plot(vgm_resid_omni$dist, vgm_resid_omni$gamma, xlab = 'Distance (km)',
      ylab = 'Semivariogram', cex.lab = 2, cex.axis = 1.5, pch = 19, 
      ylim = c(0, max(vgm_resid_omni$gamma)),
      cex = 5*sqrt(vgm_resid_omni$np)/max(sqrt(vgm_resid_omni$np)))
    mtext('B', adj = -.15, padj = -.4, cex = 3)
    par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
