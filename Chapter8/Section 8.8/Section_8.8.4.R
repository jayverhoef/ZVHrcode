sec_path = 'Rcode/Chapter8/Section 8.8/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Get the Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# attach data library
library(ZVHdata)
library(sf)
library(viridis)
library(classInt)
library(colorspace)
library(gstat)
library(spmodel)
library(stringr)

# load data for graphics and analysis
data(MOSSobs)
data(MOSSpreds)
data(CAKRboundary)

# transform some of the variables
DF = data.frame(st_drop_geometry(MOSSobs), 
	easting = st_coordinates(MOSSobs)[,1]/1e+3,
	northing = st_coordinates(MOSSobs)[,2]/1e+3)
DF$year = as.factor(DF$year)
DF$field_dup = as.factor(DF$field_dup)
DF$dist2road = log(DF$dist2road)
DF$Pb = log(DF$Pb)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Sample Design with Autocorrelation Points
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

DF01 = DF[DF$year == '2001',]
indAC = str_detect(DF01$sample, 'AC')

file_name = 'figures/Moss_autocorr_samples'
pdf(paste0(file_name,'.pdf'), width = 16, height = 8)

layout(matrix(1:2, nrow = 1))

	par(mar = c(0,0,1,0))
	plot(st_geometry(CAKRboundary))
	rect(-420000, 1985000, -410000, 1995000, lwd = 2, col = 'grey80')
	plot(st_geometry(MOSSobs[MOSSobs$year == 2001,]), add = TRUE, pch = 19)
	text(-437000, 2011000, 'A', cex = 4)

	par(mar = c(5,5,3,1))
	plot(DF01[!indAC, c('easting','northing')], pch = 3, cex = 3,
		xlim = c(-420,-410), ylim = c(1985,1995), cex.lab = 2, cex.axis = 1.5,
		xlab = 'Northing', ylab = 'Easting', lwd = 2)
	points(DF01[indAC,c('easting','northing')], pch = 4, cex = 3, lwd = 2)
	mtext('B', adj = -.2, cex = 4, padj = .4)

	layout(1)

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
#                  Fit a 3-way Model with Exponential Covariance
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# all terms
PbOut = splm(Pb ~ year*dist2road*sideroad, data = DF,
	xcoord = 'easting', ycoord = 'northing', spcov_type = 'exponential',
	partition_factor = ~ year,
	random = ~ (1 | sample) + (1 |sample:field_dup),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(PbOut)

# remove main effect sideroad

PbOut = splm(Pb ~ year + dist2road + year:dist2road + year:sideroad +
	dist2road:sideroad + year:dist2road:sideroad, data = DF,
	xcoord = 'easting', ycoord = 'northing', spcov_type = 'exponential',
	partition_factor = ~ year,
	random = ~ (1 | sample) + (1 |sample:field_dup),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(PbOut)

# remove year:sideroad

PbOut = splm(Pb ~ year + dist2road + year:dist2road +
	dist2road:sideroad + year:dist2road:sideroad, data = DF,
	xcoord = 'easting', ycoord = 'northing', spcov_type = 'exponential',
	partition_factor = ~ year,
	random = ~ (1 | sample) + (1 |sample:field_dup),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(PbOut)

# remove year:dist2road:sideroad

PbOut = splm(Pb ~ year + dist2road + year:dist2road +
	dist2road:sideroad, data = DF,
	xcoord = 'easting', ycoord = 'northing', spcov_type = 'exponential',
	partition_factor = ~ year,
	random = ~ (1 | sample) + (1 |sample:field_dup),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(PbOut)

# remove year:dist2road

PbOutFinal = splm(Pb ~ year + dist2road + dist2road:sideroad, data = DF,
	xcoord = 'easting', ycoord = 'northing', spcov_type = 'exponential',
	partition_factor = ~ year,
	random = ~ (1 | sample) + (1 |sample:field_dup),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(PbOutFinal)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#    Model Diagnostics Using LOOCV
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

loocv_out = loocv(PbOutFinal, cv_predict = TRUE, se.fit = TRUE)
hist((loocv_out$cv_predict - DF$Pb)/loocv_out$se.fit) 
indOutlier = abs((loocv_out$cv_predict - DF$Pb)/loocv_out$se.fit) > 4

#DFclean = DF[!indOutlier,]
#PbOutFinal = splm(Pb ~ year + dist2road + dist2road:sideroad, data = DFclean,
#	xcoord = 'easting', ycoord = 'northing', spcov_type = 'exponential',
#	partition_factor = ~ year,
#	random = ~ (1 | sample) + (1 |sample:field_dup),
#	control = list(reltol = 1e-12), estmethod = 'reml')
#summary(PbOutFinal)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#    Partition the variability among fixed effects and variance components
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

(1 - 0.812)*0.2016/(0.2016 + 0.0640 + 0.0267 + 0.0028)
(1 - 0.812)*0.064/(0.2016 + 0.0640 + 0.0267 + 0.0028)
(1 - 0.812)*0.0267/(0.2016 + 0.0640 + 0.0267 + 0.0028)
(1 - 0.812)*0.0028/(0.2016 + 0.0640 + 0.0267 + 0.0028)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Variance Estimate Bases on Lab Replications Only
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

summary(lm(Pb ~ sample,
	rbind(DF[DF$lab_rep == 2,c('sample','Pb')],
	DF[which(DF$lab_rep == 2) - 1, c('sample','Pb')])
))
0.0531^2

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                            Graph of Fitted Model
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# there will be 4 fitted submodels.  The 3-way interaction specified a different
# slope for each year/side-of-road combination, while the year*sideroad interaction
# provides for 4 intercept.  Let's compute the 4 intercepts and slopes for the full model:
int_2001N = PbOut$coefficient$fixed['(Intercept)'] 
int_2001S = sum(PbOut$coefficient$fixed[c('(Intercept)','sideroadS')]) 
int_2006N = sum(PbOut$coefficient$fixed[c('(Intercept)','year2006')]) 
int_2006S = sum(PbOut$coefficient$fixed[c('(Intercept)','sideroadS','year2006',
	'year2006:sideroadS')]) 
slp_2001N = PbOut$coefficient$fixed[c('dist2road')] 
slp_2001S = sum(PbOut$coefficient$fixed[c('dist2road','dist2road:sideroadS')])
slp_2006N = sum(PbOut$coefficient$fixed[c('dist2road','year2006:dist2road')]) 
slp_2006S = sum(PbOut$coefficient$fixed[c('dist2road','dist2road:sideroadS',
	'year2006:dist2road','year2006:dist2road:sideroadS')])

# slopes and intercepts for parsimonious model

int_2001N = PbOutFinal$coefficient$fixed['(Intercept)'] 
int_2001S = PbOutFinal$coefficient$fixed['(Intercept)'] 
int_2006N = sum(PbOutFinal$coefficient$fixed[c('(Intercept)','year2006')]) 
int_2006S = sum(PbOutFinal$coefficient$fixed[c('(Intercept)','year2006')])
slp_2001N = PbOutFinal$coefficient$fixed[c('dist2road')] 
slp_2001S = sum(PbOutFinal$coefficient$fixed[c('dist2road','dist2road:sideroadS')])
slp_2006N = PbOutFinal$coefficient$fixed[c('dist2road')] 
slp_2006S = sum(PbOutFinal$coefficient$fixed[c('dist2road','dist2road:sideroadS')])


file_name = 'figures/Moss_modelfit'
pdf(paste0(file_name,'.pdf'), width = 10, height = 10)

  xaxs = rep(0, times = dim(DF)[1])
  xaxs[DF$sideroad == 'N'] = -1
  xaxs[DF$sideroad == 'S'] = 1
  xaxs = xaxs*DF$dist2road
  old.par = par(mar = c(5,5,1,1))
  plot(xaxs, DF$Pb, xlab = 'Log Distance (ln m) From Road',
    ylab = 'Log Lead Concentration (ln Pb mg/kg)', pch = 19, type = 'n',
    cex = 1.5, cex.lab = 2, cex.axis = 1.5,
    xaxt = 'n', col = 'white', ylim = c(-.25,7.8))
  axis(1, at = c(-10,-5,0,5,10), label = c('10','5','0','5','10'), 
		cex.axis = 1.5)
  points(xaxs[DF$year == '2001'], 
    DF$Pb[DF$year == '2001'], 
    pch = 19, cex = 1.5)
  points(xaxs[DF$year == '2006'], 
    DF$Pb[DF$year == '2006'], 
    pch = 3, cex = 2, lwd = 2)
  text(-9,6,label = 'North', cex = 3)
  text(9,6,label = 'South', cex = 3)
  lines(c(0,-11),c(int_2001N, int_2001N + slp_2001N*(11)), lwd = 3, lty = 1)
  lines(c(0,-11),c(int_2006N, int_2006N + slp_2006N*(11)), lwd = 3, lty = 2)
  lines(c(0,11),c(int_2001S, int_2001S + slp_2001S*(11)), lwd = 3, lty = 1)
  lines(c(0,11),c(int_2006S, int_2006S + slp_2006S*(11)), lwd = 3, lty = 2)
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

