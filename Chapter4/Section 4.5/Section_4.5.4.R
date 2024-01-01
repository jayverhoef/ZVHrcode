sec_path = 'Rcode/Chapter4/Section 4.5/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(sf)
library(gstat)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           MOSS ANOVA Table
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# load data for graphics and analysis
data(MOSSobs)

# pull out the data.frame
MOSS_data = st_drop_geometry(MOSSobs)
# change year to a factor variable
MOSS_data$year = as.factor(as.character(MOSS_data$year))

# fit some models using OLS
lmfit1 = lm(log(Pb) ~ year + sideroad + log(dist2road) + log(dist2road)*year, 
  data = MOSS_data)
summary(lmfit1)
lmfit2 = lm(log(Pb) ~ year + sideroad + log(dist2road) + log(dist2road)*sideroad, 
  data = MOSS_data)
summary(lmfit2)
lmfit3 = lm(log(Pb) ~ year + sideroad + log(dist2road) + log(dist2road)*sideroad +
  log(dist2road)*year, data = MOSS_data)
summary(lmfit3)

# choose lmfit2 model and create ANOVA table
aovOut0 = summary(aov(log(Pb) ~ 1, data = MOSS_data))
aovOut2 = summary(aov(log(Pb) ~ year + sideroad + log(dist2road) + 
  log(dist2road)*sideroad, data = MOSS_data))
aovTab = rbind(aovOut2[[1]],aovOut0[[1]])
aovTab[dim(aovTab)[1],3:dim(aovTab)[2]] = NA
aovTab = cbind(c('Year', 
  'Side-of-Road',
  'log(Distance-to-Road)',
  'Side-of-Road times log(Distance-to-Road)', 
  'Error',
  'Corrected Total'), aovTab)
# create ANOVA table to paste into latex
library(xtable)
  print(
    xtable(aovTab, 
      align = c('l',rep('l', times = length(aovTab[1,]))),
      digits = c(0, 0, 0, 2, 2, 2, 4),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
  )

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Coefficients Table
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# create coefficients table to paste into latex
lmTab = summary(lmfit2)$coefficients[,1:2]
  print(
    xtable(lmTab, 
      align = c('l',rep('l', times = length(lmTab[1,]))),
      digits = c(0, 3, 3),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
  )

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Graph of the Fitted Model
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# create graphic of fitted model
file_name = "figures/Moss_Heavy_Pb_fit"
pdf(file = paste0(file_name,'.pdf'), width = 10, height = 10)
  xaxs = rep(0, times = dim(MOSS_data)[1])
  xaxs[MOSS_data$sideroad == 'N'] = -1
  xaxs[MOSS_data$sideroad == 'S'] = 1
  xaxs = xaxs*log(MOSS_data$dist2road)
  old.par = par(mar = c(5,5,1,1))
  plot(xaxs, log(MOSS_data$Pb), xlab = 'Log Distance (ln m) From Road',
    ylab = 'Log Lead Concentration (ln Pb mg/kg)', pch = 19, type = 'n',
    cex = 1.5, cex.lab = 2, cex.axis = 1.5,
    xaxt = 'n', col = 'white', ylim = c(-.25,7.5))
  axis(1, at = c(-10,-5,0,5,10), label = c('10','5','0','5','10'), 
		cex.axis = 1.5)
  points(xaxs[MOSS_data$year == '2001'], 
    log(MOSS_data$Pb)[MOSS_data$year == '2001'], 
    pch = 19, cex = 1.5)
  points(xaxs[MOSS_data$year == '2006'], 
    log(MOSS_data$Pb)[MOSS_data$year == '2006'], 
    pch = 3, cex = 2, lwd = 2)
  text(-9,6,label = 'North', cex = 3)
  text(9,6,label = 'South', cex = 3)
  # 2001 and north
  lines(-(0:11), lmTab[1,1] + lmTab[4,1]*(0:11), lwd = 2)
  # 2006 and north
  lines(-(0:11), lmTab[1,1] + lmTab[2,1] + lmTab[4,1]*(0:11), lty = 2, lwd = 2)
  # 2001 and south
  lines((0:11), lmTab[1,1] + lmTab[3,1] + -(lmTab[4,1] + lmTab[5,1])*-(0:11), lwd = 2)
  # 2006 and south
  lines((0:11), lmTab[1,1] + lmTab[2,1] + lmTab[3,1] + 
    -(lmTab[4,1] + lmTab[5,1])*-(0:11), lty = 2, lwd = 2)
  legend(-2.7, 3, legend = c('2001', '2006'),
    pch = c(19, 3), cex = 2.5)
  legend(-3.75, 1.3, legend = c('2001','2006'), lty = c(1,2),
    cex = 2.5, lwd = 2)
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
#           Semivariograms of Raw data and Residuals
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# add residuals to data.frame
MOSS_data$resids = resid(lmfit2)
# get the coordinates from the SpatialPointsDataFrame
xy = st_coordinates(MOSSobs)
# add the coordinates to the data.frame
colnames(xy) = c('x','y')
MOSS_data$x = xy[,'x']
MOSS_data$y = xy[,'y']

file_name = "figures/Moss_Heavy_Pb_OLSvariograms"
# create graphic of various variograms using gstat
pdf(file = paste0(file_name,'.pdf'), width = 10, height = 10)

  cexAB = 4
  cexCD = 5
  layout(matrix(1:4, nrow = 2, byrow = TRUE))

  vgm_raw <- variogram(log(Pb) ~ 1, loc= ~x+y, data=MOSS_data, cutoff = 15000, width = 1000)
  vgm_raw$dist = vgm_raw$dist/1000
  old_par = par(mar = c(5,5,3,1))
  plot(vgm_raw$dist, vgm_raw$gamma, 
		cex = cexAB*sqrt(vgm_raw$np)/max(sqrt(vgm_raw$np)),
    pch = 19, xlab = 'Distance (km)', ylab = 'Semivariogram', cex.axis = 1.5, cex.lab = 2)
  mtext('A', adj = -.17, cex = 3)
  par(old_par)

  vgm_resid <- variogram(resids~1, loc= ~x+y, data=MOSS_data, cutoff = 15000, width = 1000)
  vgm_resid$dist = vgm_resid$dist/1000
  old_par = par(mar = c(5,5,3,1))
  plot(vgm_resid$dist, vgm_resid$gamma, 
		cex = cexAB*sqrt(vgm_resid$np)/max(sqrt(vgm_resid$np)),
    pch = 19, xlab = 'Distance (km)', ylab = 'Semivariogram', cex.axis = 1.5, cex.lab = 2)
  mtext('B', adj = -.17, cex = 3)
  par(old_par)

  vgm_resid_dir <- variogram(log(Pb)~1, loc=~x+y, data=MOSS_data,
       alpha=c(0,45,90,135), cutoff = 15000, width = 1000)
  vgm_resid_dir$dist = vgm_resid_dir$dist/1000
  vgm_resid_dir = vgm_resid_dir[vgm_resid_dir$dist < 35,]
  old_par = par(mar = c(5,5,3,1))
  plot(vgm_resid_dir$dist, vgm_resid_dir$gamma, type = 'n', xlab = 'Distance (km)',
    ylab = 'Semivariogram', cex.lab = 2, cex.axis = 1.5)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], pch = 1, 
    cex = cexCD*sqrt(vgm_resid_dir$np)/max(sqrt(vgm_resid_dir$np)))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], pch = 2, 
    cex = cexCD*sqrt(vgm_resid_dir$np)/max(sqrt(vgm_resid_dir$np)))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], pch = 5, 
    cex = cexCD*sqrt(vgm_resid_dir$np)/max(sqrt(vgm_resid_dir$np)))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], pch = 22, 
    cex = cexCD*sqrt(vgm_resid_dir$np)/max(sqrt(vgm_resid_dir$np)))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], lwd = 2)
  legend(.1, 2.75, legend=(c(0,45,90,135)), pch = c(1,2,5,22), cex = 1.5)
  mtext('C', adj = -.17, cex = 3)
  par(old_par)

  vgm_resid_dir <- variogram(resids~1, loc=~x+y, data=MOSS_data,
       alpha=c(0,45,90,135), cutoff = 15000, width = 1000)
  vgm_resid_dir$dist = vgm_resid_dir$dist/1000
  vgm_resid_dir = vgm_resid_dir[vgm_resid_dir$dist < 35,]
  old_par = par(mar = c(5,5,3,1))
  plot(vgm_resid_dir$dist, vgm_resid_dir$gamma, type = 'n', xlab = 'Distance (km)',
    ylab = 'Semivariogram', cex.lab = 2, cex.axis = 1.5)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], pch = 1, 
    cex = cexCD*sqrt(vgm_resid_dir$np)/max(sqrt(vgm_resid_dir$np)))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 0,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], pch = 2, 
    cex = cexCD*sqrt(vgm_resid_dir$np)/max(sqrt(vgm_resid_dir$np)))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 45,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], pch = 5, 
    cex = cexCD*sqrt(vgm_resid_dir$np)/max(sqrt(vgm_resid_dir$np)))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 90,'gamma'], lwd = 2)
  points(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], pch = 22, 
    cex = cexCD*sqrt(vgm_resid_dir$np)/max(sqrt(vgm_resid_dir$np)))
  lines(vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'dist'],
    vgm_resid_dir[vgm_resid_dir$dir.hor == 135,'gamma'], lwd = 2)
  #legend(.1,.26, legend=(c(0,45,90,135)), pch = c(1,2,5,22), cex = 1.4)
  mtext('D', adj = -.17, cex = 3)
  par(old_par)

dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

