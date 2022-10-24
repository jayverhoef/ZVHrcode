sec_path = 'Rcode/Chapter4/Section 4.4/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(sf)
data(SO4obs)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Figure corr-not-trend 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/corr-not-trend"

# Create plot of simulated data on transect with strong positive correlation and no trend (Figure 4.3)
pdf(paste0(file_name,'.pdf'), height = 7, width = 9)

  y <- matrix(0,20,1)
  set.seed(1016)
  for(i in 2:20){
    y[i,1] <- 0.95*y[i-1,1]+rnorm(1,0,0.1)
  }
  i <- c(1:20)
  old.par = par(mar = c(4,4.5,.5,.5))
  plot(i, y, pch = 19, cex = 2.5, cex.lab = 2, cex.axis = 1.5)
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
#          VIF original and centered trend surface analysis for SO4 data 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# VIF original and centered trend surface analysis for SO4 data

# from Section 3.6, remove the outlers
SO4clean = SO4obs[!(1:length(SO4obs) %in% c(146,153,173)),]
xy = st_coordinates(SO4clean)
str(xy)

# create polynomials as data.frame
xy = data.frame(x = xy[,1], y = xy[,2], x2 = xy[,1]^2, y2 = xy[,2]^2, 
  xy = xy[,1]*xy[,2], x3 = xy[,1]^3, y3 = xy[,2]^3, x2y = xy[,1]^2*xy[,2],
  xy2 = xy[,1]*xy[,2]^2)
# create centered polynomials as data.frame
xy_c = cbind(xy[,1] - mean(xy[,1]), xy[,2] - mean(xy[,2]))
xy_c = data.frame(x = xy_c[,1], y = xy_c[,2], x2 = xy_c[,1]^2, y2 = xy_c[,2]^2, 
  xy = xy_c[,1]*xy_c[,2], x3 = xy_c[,1]^3, y3 = xy_c[,2]^3, x2y = xy_c[,1]^2*xy_c[,2],
  xy2 = xy_c[,1]*xy_c[,2]^2)

# VIF for uncentered quadratic model
VIF_x = 1/(1 - summary(lm(x ~ y + x2 + y2 + xy, data = xy))$r.squared)
VIF_y = 1/(1 - summary(lm(y ~ x + x2 + y2 + xy, data = xy))$r.squared)
VIF_x2 = 1/(1 - summary(lm(x2 ~ x + y + y2 + xy, data = xy))$r.squared)
VIF_y2 = 1/(1 - summary(lm(y2 ~ x + y + x2 + xy, data = xy))$r.squared)
VIF_xy = 1/(1 - summary(lm(xy ~ x + y + x2 + y2, data = xy))$r.squared)
c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy)
min(c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy))
max(c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy))

# VIF for centered quadratic model
VIF_x = 1/(1 - summary(lm(x ~ y + x2 + y2 + xy, data = xy_c))$r.squared)
VIF_y = 1/(1 - summary(lm(y ~ x + x2 + y2 + xy, data = xy_c))$r.squared)
VIF_x2 = 1/(1 - summary(lm(x2 ~ x + y + y2 + xy, data = xy_c))$r.squared)
VIF_y2 = 1/(1 - summary(lm(y2 ~ x + y + x2 + xy, data = xy_c))$r.squared)
VIF_xy = 1/(1 - summary(lm(xy ~ x + y + x2 + y2, data = xy_c))$r.squared)
c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy)
min(c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy))
max(c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy))

# VIF for uncentered cubic model
VIF_x = 1/(1 - summary(lm(x ~ y + x2 + y2 + xy + x3 + y3 + x2y + xy2, 
  data = xy))$r.squared)
VIF_y = 1/(1 - summary(lm(y ~ x + x2 + y2 + xy + x3 + y3 + x2y + xy2, 
  data = xy))$r.squared)
VIF_x2 = 1/(1 - summary(lm(x2 ~ x + y + y2 + xy + x3 + y3 + x2y + xy2, 
  data = xy))$r.squared)
VIF_y2 = 1/(1 - summary(lm(y2 ~ x + y + x2 + xy + x3 + y3 + x2y + xy2, 
  data = xy))$r.squared)
VIF_xy = 1/(1 - summary(lm(xy ~ x + y + x2 + y2 + x3 + y3 + x2y + xy2, 
  data = xy))$r.squared)
VIF_x3 = 1/(1 - summary(lm(x3 ~ x + y + x2 + y2 + xy + y3 + x2y + xy2, 
  data = xy))$r.squared)
VIF_y3 = 1/(1 - summary(lm(y3 ~ x + y + x2 + y2 + xy + x3 + x2y + xy2, 
  data = xy))$r.squared)
VIF_x2y = 1/(1 - summary(lm(x2y ~ x + y + x2 + y2 + xy + x3 + y3 + xy2, 
  data = xy))$r.squared)
VIF_xy2 = 1/(1 - summary(lm(xy2 ~ x + y + x2 + y2 + xy + x3 + y3 + x2y, 
  data = xy))$r.squared)
c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy, VIF_x3, VIF_y3, VIF_x2y, VIF_xy2)
min(c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy, VIF_x3, VIF_y3, VIF_x2y, VIF_xy2))
max(c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy, VIF_x3, VIF_y3, VIF_x2y, VIF_xy2))

# VIF for centered cubic model
VIF_x = 1/(1 - summary(lm(x ~ y + x2 + y2 + xy + x3 + y3 + x2y + xy2, 
  data = xy_c))$r.squared)
VIF_y = 1/(1 - summary(lm(y ~ x + x2 + y2 + xy + x3 + y3 + x2y + xy2, 
  data = xy_c))$r.squared)
VIF_x2 = 1/(1 - summary(lm(x2 ~ x + y + y2 + xy + x3 + y3 + x2y + xy2, 
  data = xy_c))$r.squared)
VIF_y2 = 1/(1 - summary(lm(y2 ~ x + y + x2 + xy + x3 + y3 + x2y + xy2, 
  data = xy_c))$r.squared)
VIF_xy = 1/(1 - summary(lm(xy ~ x + y + x2 + y2 + x3 + y3 + x2y + xy2, 
  data = xy_c))$r.squared)
VIF_x3 = 1/(1 - summary(lm(x3 ~ x + y + x2 + y2 + xy + y3 + x2y + xy2, 
  data = xy_c))$r.squared)
VIF_y3 = 1/(1 - summary(lm(y3 ~ x + y + x2 + y2 + xy + x3 + x2y + xy2, 
  data = xy_c))$r.squared)
VIF_x2y = 1/(1 - summary(lm(x2y ~ x + y + x2 + y2 + xy + x3 + y3 + xy2, 
  data = xy_c))$r.squared)
VIF_xy2 = 1/(1 - summary(lm(xy2 ~ x + y + x2 + y2 + xy + x3 + y3 + x2y, 
  data = xy_c))$r.squared)
c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy, VIF_x3, VIF_y3, VIF_x2y, VIF_xy2)
min(c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy, VIF_x3, VIF_y3, VIF_x2y, VIF_xy2))
max(c(VIF_x, VIF_y, VIF_x2, VIF_y2, VIF_xy, VIF_x3, VIF_y3, VIF_x2y, VIF_xy2))

