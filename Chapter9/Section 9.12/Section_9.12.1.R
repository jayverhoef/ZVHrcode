sec_path = 'Rcode/Chapter9/Section 9.12/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Get the Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# attach data library
library(ZVHdata)
library(sf)
library(viridis)
library(classInt)
library(spmodel)
library(xtable)
#library(colorspace)
# load data for graphics and analysis
data(SO4obs)


# from Section 3.6, remove the outlers and use sqrt of response
SO4clean = SO4obs[!(1:dim(SO4obs)[1] %in% c(146,153,173)),]
xy = st_coordinates(SO4clean)
# change spatial coordinates to 1000km units, rather than meters
# we will be making polynomials on the coordinates, and such large values can
# cause computer overflows and loss of precision
DF = data.frame(y = sqrt(SO4clean$SO4), easting = xy[,1]/1e+6, 
	northing = xy[,2]/1e+6)


################################################################################
#-------------------------------------------------------------------------------
#         Fit Models as Series of Polynomials, Order 0 to 5
#         Covariance Structure is Independence and Exponential Model
#         Spatial Fitting Methods are MLE and REMLE
#-------------------------------------------------------------------------------
################################################################################

# Independence Models
lm_0 = splm(y ~ 1, data = DF, estmethod = 'ml',
	xcoord = easting, ycoord = northing, spcov_type = "none")
lm_1 = splm(y ~ poly(easting, northing, degree = 1, raw = TRUE), data = DF, 
	spcov_type = "none", estmethod = 'ml')
lm_2 = splm(y ~ poly(easting, northing, degree = 2, raw = TRUE), data = DF, 
	spcov_type = "none", estmethod = 'ml')
lm_3 = splm(y ~ poly(easting, northing, degree = 3, raw = TRUE), data = DF, 
	spcov_type = "none", estmethod = 'ml')
lm_4 = splm(y ~ poly(easting, northing, degree = 4, raw = TRUE), data = DF, 
	spcov_type = "none", estmethod = 'ml')
lm_5 = splm(y ~ poly(easting, northing, degree = 5, raw = TRUE), data = DF, 
	spcov_type = "none", estmethod = 'ml')

# Spatial Exponential Models Estimated with ML and REMLE
sp_ml_0 = splm(y ~ 1, data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'ml')
sp_ml_1 = splm(y ~ poly(easting, northing, degree = 1, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'ml')
sp_ml_2 = splm(y ~ poly(easting, northing, degree = 2, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'ml')
sp_ml_3 = splm(y ~ poly(easting, northing, degree = 3, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'ml')
sp_ml_4 = splm(y ~ poly(easting, northing, degree = 4, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'ml')
sp_ml_5 = splm(y ~ poly(easting, northing, degree = 5, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'ml')

# Spatial Exponential Models Estimated with ML and REMLE
sp_reml_0 = splm(y ~ 1, data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'reml')
sp_reml_1 = splm(y ~ poly(easting, northing, degree = 1, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'reml')
sp_reml_2 = splm(y ~ poly(easting, northing, degree = 2, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'reml')
sp_reml_3 = splm(y ~ poly(easting, northing, degree = 3, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'reml')
sp_reml_4 = splm(y ~ poly(easting, northing, degree = 4, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'reml')
sp_reml_5 = splm(y ~ poly(easting, northing, degree = 5, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, spcov_type = "exponential",
	estmethod = 'reml')


################################################################################
#-------------------------------------------------------------------------------
#          Leave-one-out Cross-validation (LOOCV)
#-------------------------------------------------------------------------------
################################################################################

loocv_lm_rmspe = c(
	sqrt(loocv(lm_0)),
	sqrt(loocv(lm_1)),
	sqrt(loocv(lm_2)),
	sqrt(loocv(lm_3)),
	sqrt(loocv(lm_4)),
	sqrt(loocv(lm_5))
)

loocv_spml_rmspe = c(
	sqrt(loocv(sp_ml_0)),
	sqrt(loocv(sp_ml_1)),
	sqrt(loocv(sp_ml_2)),
	sqrt(loocv(sp_ml_3)),
	sqrt(loocv(sp_ml_4)),
	sqrt(loocv(sp_ml_5))
)

loocv_spreml_rmspe = c(
	sqrt(loocv(sp_reml_0)),
	sqrt(loocv(sp_reml_1)),
	sqrt(loocv(sp_reml_2)),
	sqrt(loocv(sp_reml_3)),
	sqrt(loocv(sp_reml_4)),
	sqrt(loocv(sp_reml_5))
)

# a function to compute prediction interval coverage using LOOCV
PIC = function(spmod, alpha = 0.1)
{
	loocvout = loocv(spmod, cv_predict = TRUE, se.fit = TRUE)
	mean(
		loocvout$cv_predict - qnorm(1-alpha/2)*loocvout$se.fit <= 
				spmod$obdata[,as.character(spmod$formula[[2]])] &
			spmod$obdata[,as.character(spmod$formula[[2]])] <= 
				loocvout$cv_predict + qnorm(1-alpha/2)*loocvout$se.fit)
}	

loocv_lm_pic90 = c(
	PIC(lm_0),
	PIC(lm_1),
	PIC(lm_2),
	PIC(lm_3),
	PIC(lm_4),
	PIC(lm_5)
)

loocv_spml_pic90 = c(
	PIC(sp_ml_0),
	PIC(sp_ml_1),
	PIC(sp_ml_2),
	PIC(sp_ml_3),
	PIC(sp_ml_4),
	PIC(sp_ml_5)
)

loocv_spreml_pic90 = c(
	PIC(sp_reml_0),
	PIC(sp_reml_1),
	PIC(sp_reml_2),
	PIC(sp_reml_3),
	PIC(sp_reml_4),
	PIC(sp_reml_5)
)

################################################################################
#-------------------------------------------------------------------------------
#                       N-fold Cross-validation
#-------------------------------------------------------------------------------
################################################################################

# a function for n-fold cross-validation where the data set has been
# appended with a column of sequential integers, beginning with 1,
# that indicates the grouping for creating test data sets.  It returns
# the mean-squared prediction error, and prediction interval coverage
# where the alpha-level can be specified
Nfold = function(formula, spcov_type, estmethod, data, group_col_name, 
	alpha = 0.1)
{
	# hold the true and predicted results
	true_pred = NULL
	for(i in 1:ngrp) {
		DFtemp = data
		DFtemp[data[,group_col_name] == i, as.character(formula[[2]])] = NA
		splmtemp = splm(formula, data = DFtemp, 
			xcoord = "easting", ycoord = "northing", spcov_type = spcov_type,
			estmethod = estmethod)
		predout = predict(splmtemp, interval = 'prediction', level = 1 - alpha)
		predout = data.frame(predout, data.frame(
			obsval = data[data[,group_col_name] == i, as.character(formula[[2]])]))
		true_pred = rbind(true_pred, predout)
	}
	c(mspe = mean((true_pred$fit - true_pred$obsval)^2),
		pic90 = mean(true_pred$lwr <= true_pred$obsval & 
			true_pred$obsval <=true_pred$upr) )
}

# how many groups for N-fold cross-validation
ngrp = 3

# formulas to use as we loop through the models
formulas = c(formula(y ~ 1), 
	y ~ poly(easting, northing, degree = 1, raw = TRUE),
	y ~ poly(easting, northing, degree = 2, raw = TRUE),
	y ~ poly(easting, northing, degree = 3, raw = TRUE),
	y ~ poly(easting, northing, degree = 4, raw = TRUE),
	y ~ poly(easting, northing, degree = 5, raw = TRUE) )

# set seed so it is reproducible
set.seed(101)

all_lm_results = NULL
all_spml_results = NULL
all_spreml_results = NULL
for(j in 1:10) {
	# randomly assign observations into 3 groups
	DFgrp = DF
	DFgrp$grp = as.factor(
		sample(
			rep(1:ngrp, times = ceiling(dim(DF)[1]/ngrp))[1:dim(DF)[1]], 
				size = dim(DF)[1])
	)
	# independence models
	lm_results = matrix(nrow = 6, ncol = 3)
	for(i in 1:6) {
		Nfout = Nfold(formula = formulas[[i]], 
			spcov_type = 'none', estmethod = 'ml', 
			data = DFgrp, group_col_name = 'grp')
		lm_results[i,1] = i - 1
		lm_results[i,2] = Nfout['mspe']
		lm_results[i,3] = Nfout['pic90']
	}
	# spatial mle models	
	spml_results = matrix(nrow = 6, ncol = 3)
	for(i in 1:6) {
		Nfout = Nfold(formula = formulas[[i]], 
			spcov_type = 'exponential', estmethod = 'ml', 
			data = DFgrp, group_col_name = 'grp')
		spml_results[i,1] = i - 1
		spml_results[i,2] = Nfout['mspe']
		spml_results[i,3] = Nfout['pic90']
	}
	# spatial remle models
	spreml_results = matrix(nrow = 6, ncol = 3)
	for(i in 1:6) {
		Nfout = Nfold(formula = formulas[[i]], 
			spcov_type = 'exponential', estmethod = 'reml', 
			data = DFgrp, group_col_name = 'grp')
		spreml_results[i,1] = i - 1
		spreml_results[i,2] = Nfout['mspe']
		spreml_results[i,3] = Nfout['pic90']
	}
	
	all_lm_results = rbind(all_lm_results, lm_results)
	all_spml_results = rbind(all_spml_results, spml_results)
	all_spreml_results = rbind(all_spreml_results, spreml_results)

}

# average after doing 3-fold cross-validation 10 times
nfold_lm_results = rbind(
	apply(all_lm_results[all_lm_results[,1] == 0,],2,mean),
	apply(all_lm_results[all_lm_results[,1] == 1,],2,mean),
	apply(all_lm_results[all_lm_results[,1] == 2,],2,mean),
	apply(all_lm_results[all_lm_results[,1] == 3,],2,mean),
	apply(all_lm_results[all_lm_results[,1] == 4,],2,mean),
	apply(all_lm_results[all_lm_results[,1] == 5,],2,mean)
)
nfold_spml_results = rbind(
	apply(all_spml_results[all_spml_results[,1] == 0,],2,mean),
	apply(all_spml_results[all_spml_results[,1] == 1,],2,mean),
	apply(all_spml_results[all_spml_results[,1] == 2,],2,mean),
	apply(all_spml_results[all_spml_results[,1] == 3,],2,mean),
	apply(all_spml_results[all_spml_results[,1] == 4,],2,mean),
	apply(all_spml_results[all_spml_results[,1] == 5,],2,mean)
)
nfold_spreml_results = rbind(
	apply(all_spreml_results[all_spreml_results[,1] == 0,],2,mean),
	apply(all_spreml_results[all_spreml_results[,1] == 1,],2,mean),
	apply(all_spreml_results[all_spreml_results[,1] == 2,],2,mean),
	apply(all_spreml_results[all_spreml_results[,1] == 3,],2,mean),
	apply(all_spreml_results[all_spreml_results[,1] == 4,],2,mean),
	apply(all_spreml_results[all_spreml_results[,1] == 5,],2,mean)
)
nfold_lm_results[,2] = sqrt(nfold_lm_results[,2])
nfold_spml_results[,2] = sqrt(nfold_spml_results[,2])
nfold_spreml_results[,2] = sqrt(nfold_spreml_results[,2])


file_name = "figures/SO4_crossval"

pdf(paste0(file_name,'.pdf'), width = 12, height = 12)
	padj = -.5
	adj = -.16
	cex_lab = 2.5
	cex_axis = 2
	cex_main = 3
	cex_sym = 4
	cex_sym_lm = 3
	layout(matrix(1:4, ncol = 2, byrow = TRUE))
	old.par = par(mar = c(5,5,5,3))	
	# A
	up = max(loocv_lm_rmspe, loocv_spml_rmspe, loocv_spreml_rmspe)
	lo = min(loocv_lm_rmspe, loocv_spml_rmspe, loocv_spreml_rmspe)
	plot(0:5, loocv_lm_rmspe, ylim = c(lo, .9), pch = 19, cex = cex_sym_lm,
		xlab = 'Order of Polynomial', ylab = 'RMSPE', cex.lab = cex_lab, cex.axis = cex_axis,
		main = 'LOOCV RMSPE', cex.main = cex_main)
	lines(0:5, loocv_lm_rmspe, lwd = 3)
	points(0:5, loocv_spml_rmspe, pch = 22, cex = cex_sym)
	lines(0:5, loocv_spml_rmspe, lty = 2, lwd = 3)
	points(0:5, loocv_spreml_rmspe, pch = 24, cex = cex_sym)
	lines(0:5, loocv_spreml_rmspe, lty = 3, lwd = 3)
	legend(2.3, .902, legend = c('Indep Model','Spatial MLE', 'Spatial REMLE'), 
		lty = c(1,2,3), lwd = 2, pch = c(19,22,24), cex = 1.7)
	mtext('A', adj = adj, cex = 4, padj = padj)
	# B
	up = max(loocv_lm_pic90, loocv_spml_pic90, loocv_spreml_pic90)
	lo = min(loocv_lm_pic90, loocv_spml_pic90, loocv_spreml_pic90)
	plot(0:5, loocv_lm_pic90, ylim = c(lo, 1.003*up), xlim = c(-.1, 5.1), 
		pch = 19, 
		cex = cex_sym_lm, xlab = 'Order of Polynomial', ylab = 'Interval Coverage', 
		cex.lab = cex_lab, cex.axis = cex_axis, main = 'LOOCV PIC90', cex.main = cex_main)
	points(0:5, loocv_spml_pic90, pch = 22, cex = cex_sym)
	points(0:5, loocv_spreml_pic90, pch = 24, cex = cex_sym)
	lines(c(-.2,5.2),c(0.9,0.9), lty = 2, lwd = 2)
	legend(2.6, .962, legend = c('Indep Model','Spatial MLE', 'Spatial REMLE'), 
		pch = c(19,22,24), cex = 1.7)
	mtext('B', adj = adj, cex = 4, padj = padj)
	# C
	up = max(nfold_lm_results[,2], nfold_spml_results[,2], 
		nfold_spreml_results[,2])
	lo = min(nfold_lm_results[,2], nfold_spml_results[,2], 
		nfold_spreml_results[,2])
	plot(nfold_lm_results[,1], nfold_lm_results[,2], ylim = c(lo, .9), 
		pch = 19, cex = cex_sym_lm, xlab = 'Order of Polynomial', ylab = 'RMSPE', 
		cex.lab = cex_lab, cex.axis = cex_axis,
		main = 'N-FOLD RMSPE', cex.main = cex_main)
	lines(nfold_lm_results[,1], nfold_lm_results[,2], lwd = 3)
	points(0:5, nfold_spml_results[,2], pch = 22, cex = cex_sym)
	lines(0:5, nfold_spml_results[,2], lty = 2, lwd = 3)
	points(0:5, nfold_spreml_results[,2], pch = 24, cex = cex_sym)
	lines(0:5, nfold_spreml_results[,2], lty = 3, lwd = 3)
	legend(2.3, .902, legend = c('Indep Model','Spatial MLE', 'Spatial REMLE'), 
		lty = c(1,2,3), lwd = 2, pch = c(19,22,24), cex = 1.7)
	mtext('C', adj = adj, cex = 4, padj = padj)
	# D
	up = max(nfold_lm_results[,3], nfold_spml_results[,3], 
		nfold_spreml_results[,3])
	lo = min(nfold_lm_results[,3], nfold_spml_results[,3], 
		nfold_spreml_results[,3])
	plot(0:5, nfold_lm_results[,3], ylim = c(lo, up), 
		pch = 19, cex = cex_sym_lm, xlab = 'Order of Polynomial', 
		ylab = 'Interval Coverage', 
		cex.lab = cex_lab, cex.axis = cex_axis,
		main = 'N-FOLD PIC90', cex.main = cex_main)
	points(0:5, nfold_spml_results[,3], pch = 22, cex = cex_sym)
	points(0:5, nfold_spreml_results[,3], pch = 24, cex = cex_sym)
	lines(c(-.2,5.2),c(0.9,0.9), lty = 2, lwd = 2)
	legend(2.6, .944, legend = c('Indep Model','Spatial MLE', 'Spatial REMLE'), 
		pch = c(19,22,24), cex = 1.7)
	mtext('D', adj = adj, cex = 4, padj = padj)

	par(old.par)
	
	layout(1)
	
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
#            Choosing Among Spatial Autocovariance Models
#-------------------------------------------------------------------------------
################################################################################

comp_cormodels = NULL
cormodel_list = c('exponential', 'spherical','gaussian','circular','pentaspherical',
	'wave','jbessel','gravity','rquad','magnetic','matern','cauchy','pexponential')
DFgrp = DF
set.seed(103)
ngrp = 3
DFgrp$grp = as.factor(
	sample(
		rep(1:ngrp, times = ceiling(dim(DF)[1]/ngrp))[1:dim(DF)[1]], 
			size = dim(DF)[1])
)
j = 1
for(j in 1:length(cormodel_list)) {
	splmfit = splm(y ~ 1, 
		data = DFgrp, xcoord = easting, ycoord = northing, 
		spcov_type = cormodel_list[j], estmethod = 'reml')
	Nfout = Nfold(y ~ 1, spcov_type = cormodel_list[j], estmethod = 'reml', 
		data = DFgrp, group_col_name = 'grp', alpha = 0.1)
loocv(splmfit, cv_predict = TRUE, se.fit = TRUE)
	comp_cormodels = rbind(comp_cormodels,
		data.frame(
			cormod = cormodel_list[j], 
			m2ll= -2*logLik(splmfit),
			AIC = AIC(splmfit),
			loocv = sqrt(loocv(splmfit)), 
			cv3fold = sqrt(Nfout['mspe']),
			PIC90_Lo = PIC(splmfit),
			PIC90_Nf = Nfout['pic90']  )
	)
}
comp_cormodels
comp_cormodels$BIC = NA
comp_cormodels[1:10,'BIC'] = comp_cormodels[1:10,'m2ll'] + 3*log(dim(DFgrp)[1])
comp_cormodels[11:13,'BIC'] = comp_cormodels[11:13,'m2ll'] + 4*log(dim(DFgrp)[1])
comp_cormodels = comp_cormodels[,
	c('cormod','m2ll','AIC','BIC','loocv','cv3fold','PIC90_Lo','PIC90_Nf')]

print(
    xtable(comp_cormodels, 
      align = c('l',rep('l', times = length(comp_cormodels[1,]))),
      digits = c(0,0,rep(1, times = 3), rep(3, times = 4)),
      caption = 'Metrics for model selection',
      label = 'tab:ModSelMetrics'
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

################################################################################
#-------------------------------------------------------------------------------
#       Mapping Spatial versus High-Order-Polynomial-Indpendence Models
#-------------------------------------------------------------------------------
################################################################################

data(USboundary)
st_bbox(USboundary)[3] - st_bbox(USboundary)[1]
st_bbox(USboundary)[4] - st_bbox(USboundary)[2]
# spacing
spacing = (st_bbox(USboundary)[3] - st_bbox(USboundary)[1])/100
xpredlocs = st_bbox(USboundary)[1] + spacing*(1:100 - 0.5)
ypredlocs = st_bbox(USboundary)[2] + spacing*(1:100 - 0.5)
# create a data frame of predictions covering all of US
preds = data.frame(x = xpredlocs %x% rep(1, times = 100), 
	y = rep(1, times = 100) %x% ypredlocs)
# turn it into a spatial points data frame, 
# give it the same projection as US polygons
preds = st_as_sf(preds, coords = 1:2, crs = st_crs(USboundary))
# clip the points to be within US borders
preds_sub = st_filter(preds, USboundary)

plot(preds_sub, pch = 19, cex = .3)
# check the result
plot(USboundary, add = TRUE, col = 'transparent', lwd = 5)
plot(SO4clean, add = TRUE, pch = 19, cex = .8, col = 'blue')
# get the coordinates of the clipped points
DFpred = as.data.frame(st_coordinates(preds_sub))
colnames(DFpred) = c('easting', 'northing')
# change coordinates to 1000 km
DFpred = DFpred/1e+6

#         Fit some models
sp_sph = splm(y~1, data = DF, spcov_type = 'spherical',
			xcoord = "easting", ycoord = "northing")
sp_gau = splm(y~1, data = DF, spcov_type = 'gaussian',
			xcoord = "easting", ycoord = "northing")
sp_sph3 = splm(y ~ poly(easting, northing, degree = 3, raw = TRUE), data = DF,
	spcov_type = 'spherical', xcoord = "easting", ycoord = "northing")
# predict points in DFpred using independence model with 4th order polynomial
pred_ind = predict(lm_4, DFpred, se.fit = TRUE)
# predict points with constant mean spherical model
pred_sp1 = predict(sp_sph, DFpred, se.fit = TRUE)
# predict points with constant mean gaussian model
pred_sp2 = predict(sp_gau, DFpred, se.fit = TRUE)
# fit 3rd order polynomial with spherical model
pred_sp3 = predict(sp_sph3, DFpred, se.fit = TRUE)
# predict points with constant mean exponential model
pred_sp4 = predict(sp_reml_0, DFpred, se.fit = TRUE)

file_name = "figures/SO4_Prediction_Maps"

pdf(paste0(file_name,'.pdf'), width = 16, height = 24)
brks_cex = 2.3
printF = "1.2"
mtext.cex = 4.8
map_mar = c(0,4,7,2)
mtext.y = 3100000
mtext.x = -2400000
leg_right = 0.6
pch_num = 15

brks_pred = quantile(c(pred_ind$fit, pred_sp1$fit, 
	pred_sp2$fit, pred_sp3$fit, pred_sp4$fit), probs = seq(0, 1, .1))
brks_se = quantile(c(pred_ind$se.fit, pred_sp1$se.fit, 
	pred_sp2$se.fit, pred_sp3$se.fit, pred_sp4$se.fit), probs = seq(0, 1, .125))

# make a color map of predictions
cip = classIntervals(pred_ind$fit, style = 'fixed', fixedBreaks= brks_pred)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = map_mar, bg = 'gray80')
source('addBreakColorLegend.R')
layout(matrix(1:20, ncol = 4, byrow = TRUE), 
	widths = rep(c(3,.4),times = 10))
plot(st_coordinates(preds_sub), col = cip_colors, pch = pch_num, cex = 1.5, 
	bty = 'n', xaxt = 'n', yaxt = 'n')
plot(USboundary, add = TRUE, border = 'black', col = 'transparent')
text(mtext.x, mtext.y, 'A', cex = mtext.cex)
old.par = par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

# make a color map of prediction standard errors
cip = classIntervals(pred_ind$se.fit, style = 'fixed', fixedBreaks= brks_se)
palp = viridis(7, option = 'C')
cip_colors = findColours(cip, palp)
old.par = par(mar = map_mar)
plot(st_coordinates(preds_sub), col = cip_colors, pch = pch_num, cex = 1.5, 
	bty = 'n', xaxt = 'n', yaxt = 'n')
plot(USboundary, add = TRUE, border = 'black', col = 'transparent')
points(st_coordinates(SO4clean), pch = 19)
text(mtext.x, mtext.y, 'B', cex = mtext.cex)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

#         Spatial Models

# spherical model, constant mean

cip = classIntervals(pred_sp1$fit, style = 'fixed', fixedBreaks= brks_pred)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = map_mar)
plot(st_coordinates(preds_sub), col = cip_colors, pch = pch_num, cex = 1.5, 
	bty = 'n', xaxt = 'n', yaxt = 'n')
plot(USboundary, add = TRUE, border = 'black', col = 'transparent')
text(mtext.x, mtext.y, 'C', cex = mtext.cex)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

cip = classIntervals(pred_sp1$se.fit, style = 'fixed', fixedBreaks= brks_se)
palp = viridis(7, option = 'C')
cip_colors = findColours(cip, palp)
old.par = par(mar = map_mar)
plot(st_coordinates(preds_sub), col = cip_colors, pch = pch_num, cex = 1.5, 
	bty = 'n', xaxt = 'n', yaxt = 'n')
plot(USboundary, add = TRUE, border = 'black', col = 'transparent')
points(coordinates(SO4clean), pch = 19)
text(mtext.x, mtext.y, 'D', cex = mtext.cex)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

# exponential model, constant mean

cip = classIntervals(pred_sp4$fit, style = 'fixed', fixedBreaks= brks_pred)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = map_mar)
plot(st_coordinates(preds_sub), col = cip_colors, pch = pch_num, cex = 1.5, 
	bty = 'n', xaxt = 'n', yaxt = 'n')
plot(USboundary, add = TRUE, border = 'black', col = 'transparent')
text(mtext.x, mtext.y, 'E', cex = mtext.cex)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

cip = classIntervals(pred_sp4$se.fit, style = 'fixed', fixedBreaks= brks_se)
palp = viridis(7, option = 'C')
cip_colors = findColours(cip, palp)
old.par = par(mar = map_mar)
plot(st_coordinates(preds_sub), col = cip_colors, pch = pch_num, cex = 1.5, 
	bty = 'n', xaxt = 'n', yaxt = 'n')
plot(USboundary, add = TRUE, border = 'black', col = 'transparent')
points(coordinates(SO4clean), pch = 19)
text(mtext.x, mtext.y, 'F', cex = mtext.cex)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

# 3rd order polynomial, spherical model

cip = classIntervals(pred_sp3$fit, style = 'fixed', fixedBreaks= brks_pred)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = map_mar)
plot(st_coordinates(preds_sub), col = cip_colors, pch = pch_num, cex = 1.5, 
	bty = 'n', xaxt = 'n', yaxt = 'n')
plot(USboundary, add = TRUE, border = 'black', col = 'transparent')
text(mtext.x, mtext.y, 'G', cex = mtext.cex)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

cip = classIntervals(pred_sp3$se.fit, style = 'fixed', fixedBreaks= brks_se)
palp = viridis(7, option = 'C')
cip_colors = findColours(cip, palp)
old.par = par(mar = map_mar)
plot(st_coordinates(preds_sub), col = cip_colors, pch = pch_num, cex = 1.5, 
	bty = 'n', xaxt = 'n', yaxt = 'n')
plot(USboundary, add = TRUE, border = 'black', col = 'transparent')
points(coordinates(SO4clean), pch = 19)
text(mtext.x, mtext.y, 'H', cex = mtext.cex)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)
  
# constant mean, Gaussian model

cip = classIntervals(pred_sp2$fit, style = 'fixed', fixedBreaks= brks_pred)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = map_mar)
plot(st_coordinates(preds_sub), col = cip_colors, pch = pch_num, cex = 1.5, 
	bty = 'n', xaxt = 'n', yaxt = 'n')
plot(USboundary, add = TRUE, border = 'black', col = 'transparent')
text(mtext.x, mtext.y, 'I', cex = mtext.cex)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

cip = classIntervals(pred_sp2$se.fit, style = 'fixed', fixedBreaks= brks_se)
palp = viridis(7, option = 'C')
cip_colors = findColours(cip, palp)
old.par = par(mar = map_mar)
plot(st_coordinates(preds_sub), col = cip_colors, pch = pch_num, cex = 1.5, 
	bty = 'n', xaxt = 'n', yaxt = 'n')
plot(USboundary, add = TRUE, border = 'black', col = 'transparent')
points(coordinates(SO4clean), pch = 19)
text(mtext.x, mtext.y, 'J', cex = mtext.cex)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

	par(old.par)
	
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))



