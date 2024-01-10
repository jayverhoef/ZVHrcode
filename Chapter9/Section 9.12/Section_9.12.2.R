sec_path = 'Rcode/Chapter9/Section 9.12/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Get the Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# attach data library
library(ZVHdata)
library(sp)
library(spdep)
library(sf)
library(viridis)
library(classInt)
library(colorspace)
library(spmodel)
source('addBreakColorLegend.R')

# load data for graphics and analysis
data(sealPolys)
seals_sf = st_as_sf(sealPolys)
seals_sf$stockname = as.factor(as.character(seals_sf$stockname))

################################################################################
#-------------------------------------------------------------------------------
#                    Create Neighborhood Matrices
#-------------------------------------------------------------------------------
################################################################################

nTot = length(sealPolys)
nObs = sum(!is.na(sealPolys$Estimate))
nMiss = nTot - nObs

# a function to create matrix from lists and numbers of neighbors
Neighmat <- function(adj, num, n)
{
	N.mat <- matrix( 0, nrow = n, ncol = n )
	k <- 1
	for (i in 1:n){
		if(num[i] > 0) {
				N.mat[i,adj[[i]]] <- 1
			}
		}
	N.mat
}

# use spdep to find touching neighors, and then add rest manually
Nlist = poly2nb(sealPolys, snap = 2000)
Nlist[[79]] = as.integer(c(211, 463))
Nlist[[211]] = as.integer(c(Nlist[[211]]), 79)
Nlist[[463]] = as.integer(c(Nlist[[463]]), 79)
Nlist[[130]] = as.integer(302)
Nlist[[302]] = as.integer(c(Nlist[[302]],130))
Nlist[[325]] = as.integer(c(326, 353))
Nlist[[326]] = as.integer(c(Nlist[[326]],325))
Nlist[[353]] = as.integer(c(Nlist[[353]],325))
Nlist[[435]] = as.integer(437)
Nlist[[437]] = as.integer(c(Nlist[[437]],435))
Nlist[[436]] = as.integer(c(86, 88))
Nlist[[86]] = as.integer(c(Nlist[[86]],436))
Nlist[[88]] = as.integer(c(Nlist[[88]],436))
Nlist[[437]] = as.integer(87)
Nlist[[87]] = as.integer(c(Nlist[[87]],437))
Nlist[[438]] = as.integer(436)
Nlist[[436]] = as.integer(c(Nlist[[436]],438))
Nlist[[439]] = as.integer(346)
Nlist[[346]] = as.integer(c(Nlist[[346]],439))
Nlist[[443]] = as.integer(281)
Nlist[[281]] = as.integer(c(Nlist[[281]],443))
Nlist[[463]] = as.integer(79)
attr(Nlist,'polyid') = as.factor(as.character(sealPolys@data$polyid))
attr(Nlist,'stockid') = as.factor(as.character(sealPolys@data$stockid))
num = lapply(Nlist, function(x) length(x))
num = unlist(num)
Nmat = Neighmat(Nlist, num, length(num))
Nmat1 = pmax(Nmat,t(Nmat))
Nmat2 = (((Nmat1 %*% Nmat1 > 0)*1 + Nmat1) > 0)*1 - diag(dim(Nmat1)[1])
Nmat3 = (((Nmat2 %*% Nmat1 > 0)*1 + Nmat2) > 0)*1 - diag(dim(Nmat1)[1])
Nmat4 = (((Nmat3 %*% Nmat1 > 0)*1 + Nmat3) > 0)*1 - diag(dim(Nmat1)[1])
Nmat5 = (((Nmat4 %*% Nmat1 > 0)*1 + Nmat4) > 0)*1 - diag(dim(Nmat1)[1])
Nmat6 = (((Nmat5 %*% Nmat1 > 0)*1 + Nmat5) > 0)*1 - diag(dim(Nmat1)[1])


################################################################################
#-------------------------------------------------------------------------------
#                  Prediction
#-------------------------------------------------------------------------------
################################################################################


lmfit_X = splm(Estimate ~ stockname, data = seals_sf, estmethod = 'reml', 
	control = list(reltol = 1e-7), spcov_type = 'none')
summary(lmfit_X)
spfit_m = spautor(Estimate ~ 1, data = seals_sf, 
	estmethod = 'reml', spcov_type = 'car', row_st = TRUE)
summary(spfit_m)
spfit_X = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'reml', spcov_type = 'car', row_st = TRUE)
summary(spfit_m)

seals_pred = seals_sf
seals_pred[,'Estimate_se'] = NA
seals_pred_m = seals_sf
seals_pred_m[,'Estimate_se'] = NA
predout = predict(spfit_X, se.fit = TRUE)
predout_m = predict(spfit_m, se.fit = TRUE)
ind = is.na(seals_pred$Estimate)
seals_pred[ind,'Estimate'] = predout$fit
seals_pred[ind,'Estimate_se'] = predout$se.fit
seals_pred_m[ind,'Estimate'] = predout_m$fit
seals_pred_m[ind,'Estimate_se'] = predout_m$se.fit
loocvout = loocv(spfit_X, cv_predict = TRUE, se.fit = TRUE)
loocvout_m = loocv(spfit_m, cv_predict = TRUE, se.fit = TRUE)
seals_smooth = seals_pred
seals_smooth[!ind,'Estimate'] = loocvout$cv_predict
seals_smooth[!ind,'Estimate_se'] = loocvout$se.fit
seals_smooth_m = seals_pred_m
seals_smooth_m[!ind,'Estimate'] = loocvout_m$cv_predict
seals_smooth_m[!ind,'Estimate_se'] = loocvout_m$se.fit

file_name = 'figures/seal_predhist'
pdf(paste0(file_name,'.pdf'), width = 8, height = 8)

	layout(matrix(1:2, nrow = 2))
	padj = -.5
	adj = -.15
	cex.mtext = 3
	par(mar = c(5,5,4,1))
	hist(seals_sf$Estimate,  breaks = seq(-0.6, 1, by = .04),
		xlim = c(-.42, .5), col = rgb(.8, .8, .8, .3), ylim = c(0,150),
		main = '', cex.lab = 2, cex.axis = 1.5, xlab = 'log(Trend)')
	hist(seals_smooth$Estimate,  breaks = seq(-0.6, 1, by = .02),
		xlim = c(-.5, .5), add = TRUE, col = rgb(.1, .1, .1, .4))
	legend(x = .1, y = 150, legend =c('Raw Values','LOOCV Predictions'),
		pch = 15, col = c(rgb(.7, .7, .7, .5),rgb(.1, .1, .1, .6)),
		cex = 1.4)
	mtext('A', adj = adj, cex = cex.mtext, padj = padj)


	h1out = hist(predict(spfit_m), breaks = seq(-0.15, 0.25, by = .01), plot = FALSE)
	h2out = hist(predict(spfit_X), breaks = seq(-0.15, 0.25, by = .01), plot = FALSE)
	h3out = hist(predict(lmfit_X), breaks = seq(-0.15, 0.25, by = .01), plot = FALSE)

	par(mar = c(5,5,4,1), bg = 'white')
	lwd_bars = 5
	offset = .003
	padj = -.5
	adj = -.15
	cex.mtext = 3
	plot(h1out$mids, h1out$counts, type = 'n', bty = 'n',
		cex.lab = 2, cex.axis = 1.5, ylab = 'Frequency', xlab = 'log(Trend)',
		xaxt = 'n')
	for(i in 1:40) {
		if(h1out$counts[i] > 0)
			lines(rep(h1out$mids[i],2), c(0,h1out$counts[i]), lwd = lwd_bars,
				col = 'gray50', lend = 2)
	}
	for(i in 1:40) {
		if(h2out$counts[i] > 0)
			lines(rep(h2out$mids[i]-offset, 2), c(0,h2out$counts[i]), lwd = lwd_bars,
				col = 'gray80', lend = 2)
	}
	for(i in 1:40) {
		if(h3out$counts[i] > 0)
			lines(rep(h3out$mids[i]+offset, 2), c(0,h3out$counts[i]), lwd = lwd_bars,
				col = 'black', lend = 2)
	}
	lines(c(-.145,.245), c(-.1,-.1))
	axis(1, at = c(-.1,-.05, 0,.05,.1,.15, .2), cex.axis = 1.5)
	d1out = density(predict(spfit_m), bw = .02)
	lines(d1out$x, d1out$y, col = 'gray50', lwd = 4)
	d2out = density(predict(spfit_X), bw = .02)
	lines(d2out$x, d2out$y, col = 'gray80', lwd = 4)
	legend(x = .08, y = 44, legend =c('Indep. Mean','Autocorr. Mean', 
		'Autocorr. Stock'),
		lwd = 5, col = c('black','gray50','gray80'),
		cex = 1.4)
	mtext('B', adj = adj, cex = cex.mtext, padj = padj)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))




spfit_X_norowst = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'reml', spcov_type = 'car', row_st = FALSE)
spfit_X_Nmat1 = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'reml', W = Nmat1, spcov_type = 'car', row_st = TRUE)
spfit_X_Nmat4 = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'reml', W = Nmat4, spcov_type = 'car', row_st = TRUE)
loocv(spfit_X)
loocv(spfit_X_norowst)
loocv(spfit_X_Nmat1)
loocv(lmfit_X)
loocv(spfit_m)




file_name = 'figures/seal_maps'
	pdf(paste0(file_name,'.pdf'), width = 8, height = 12)

	layout(matrix(1:6, ncol = 2, byrow = TRUE))
	par(mar = c(0,0,0,0), bg = 'gray80')
	brks_pred = quantile(seals_pred$Estimate, probs = seq(0, 1, .125))
	brks_se = quantile(seals_pred[ind,]$Estimate_se, probs = seq(0, 1, .125))
	brks_smoo = quantile(seals_smooth$Estimate, probs = seq(0, 1, .125))
	brks_smoo_se = quantile(seals_smooth$Estimate_se, probs = seq(0, 1, .125))
	brks_smoo_m = quantile(seals_smooth_m$Estimate, probs = seq(0, 1, .125))
	brks_smoo_m_se = quantile(seals_smooth_m$Estimate_se, probs = seq(0, 1, .125))

	# make a color map of predictions

	# A

	cip = classIntervals(seals_sf$Estimate, style = 'fixed', 
		fixedBreaks= brks_pred)
	palp = viridis(8)
	cip_colors = findColours(cip, palp)
	plot(st_geometry(seals_sf), col = cip_colors)
	polycent = st_centroid(seals_pred)
	ind = is.na(seals_sf$Estimate)
	cip = classIntervals(seals_pred[ind,]$Estimate, style = 'fixed', 
		fixedBreaks = brks_pred)
	cip_colors = findColours(cip, palp)
	plot(st_geometry(polycent[ind,]), add = TRUE, pch = 19, cex = 1.3, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_pred, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'A', cex = 4)

	# B

	cip = classIntervals(seals_pred[ind,]$Estimate_se, style = 'fixed', 
		fixedBreaks = brks_se)
	palp = plasma(8)
	cip_colors = findColours(cip, palp)
	plot(st_geometry(seals_sf))
	plot(st_geometry(polycent[ind,]), add = TRUE, pch = 19, cex = 1.3, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_se, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'B', cex = 4)

	# C
	cip = classIntervals(seals_smooth$Estimate, style = 'fixed', 
		fixedBreaks= brks_smoo)
	palp = viridis(8)
	cip_colors = findColours(cip, palp)
	# make a color map of prediction standard errors
	plot(st_geometry(seals_smooth), col = cip_colors)
	polycent = st_centroid(seals_smooth)
	plot(st_geometry(polycent), add = TRUE, pch = 19, cex = 1, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_smoo, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'C', cex = 4)

	# D
	cip = classIntervals(seals_smooth$Estimate_se, style = 'fixed', 
		fixedBreaks= brks_smoo_se)
	palp = plasma(8)
	cip_colors = findColours(cip, palp)
	# make a color map of prediction standard errors
	plot(st_geometry(seals_smooth), col = cip_colors)
	polycent = st_centroid(seals_smooth)
	plot(st_geometry(polycent), add = TRUE, pch = 19, cex = 1, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_smoo_se, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'D', cex = 4)

	# E
	cip = classIntervals(seals_smooth_m$Estimate, style = 'fixed', 
		fixedBreaks= brks_smoo_m)
	palp = viridis(8)
	cip_colors = findColours(cip, palp)
	# make a color map of prediction standard errors
	plot(st_geometry(seals_smooth_m), col = cip_colors)
	polycent = st_centroid(seals_smooth_m)
	plot(st_geometry(polycent), add = TRUE, pch = 19, cex = 1, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_smoo_m, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'E', cex = 4)

	# F
	cip = classIntervals(seals_smooth_m$Estimate_se, style = 'fixed', 
		fixedBreaks= brks_smoo_m_se)
	palp = plasma(8)
	cip_colors = findColours(cip, palp)
	# make a color map of prediction standard errors
	plot(st_geometry(seals_smooth), col = cip_colors)
	polycent = st_centroid(seals_smooth)
	plot(st_geometry(polycent), add = TRUE, pch = 19, cex = 1, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_smoo_m_se, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'F', cex = 4)

layout(1)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))


