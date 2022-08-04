sec_path = 'Rcode/Chapter8/Section 8.8/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Get the Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# attach data library
library(ZVHdata)
library(sp)
library(viridis)
library(classInt)
library(spmodel)
library(numDeriv)
library(xtable)
#library(colorspace)
# load data for graphics and analysis
data(SO4obs)


# from Section 3.6, remove the outlers and use sqrt of response
SO4clean = SO4obs[!(1:length(SO4obs) %in% c(146,153,173)),]
xy = coordinates(SO4clean)
# change spatial coordinates to 1000km units, rather than meters
# we will be making polynomials on the coordinates, and such large values can
# cause computer overflows and loss of precision
DF = data.frame(y = sqrt(SO4clean@data$SO4), easting = xy[,1]/1e+6, 
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

# loglikelihood
m2ll_indep = c(-2*logLik(lm_0), -2*logLik(lm_1), -2*logLik(lm_2),
	-2*logLik(lm_3), -2*logLik(lm_4), -2*logLik(lm_5))
m2ll_spatial = c(-2*logLik(sp_ml_0), -2*logLik(sp_ml_1), -2*logLik(sp_ml_2),
	-2*logLik(sp_ml_3), -2*logLik(sp_ml_4), -2*logLik(sp_ml_5))

# AIC
AIC_indep = c(AIC(lm_0), AIC(lm_1), AIC(lm_2), 
	AIC(lm_3), AIC(lm_4), AIC(lm_5))
AIC_spatial = c(AIC(sp_ml_0), AIC(sp_ml_1), AIC(sp_ml_2), 
	AIC(sp_ml_3), AIC(sp_ml_4), AIC(sp_ml_5))

# BIC
BIC_indep = 
	c(-2*logLik(lm_0) + (length(coef(lm_0)) + 1)*log(dim(DF)[1]), 
	-2*logLik(lm_1) + (length(coef(lm_1)) + 1)*log(dim(DF)[1]), 
	-2*logLik(lm_2) + (length(coef(lm_2)) + 1)*log(dim(DF)[1]), 
	-2*logLik(lm_3) + (length(coef(lm_3)) + 1)*log(dim(DF)[1]),  
	-2*logLik(lm_4) + (length(coef(lm_4)) + 1)*log(dim(DF)[1]), 
	-2*logLik(lm_5) + (length(coef(lm_5)) + 1)*log(dim(DF)[1]) )
BIC_spatial = 
	c(-2*logLik(sp_ml_0) + (length(coef(sp_ml_0)) + 3)*log(dim(DF)[1]), 
	-2*logLik(sp_ml_1) + (length(coef(sp_ml_1)) + 3)*log(dim(DF)[1]),  
	-2*logLik(sp_ml_2) + (length(coef(sp_ml_2)) + 3)*log(dim(DF)[1]), 
	-2*logLik(sp_ml_3) + (length(coef(sp_ml_3)) + 3)*log(dim(DF)[1]),  
	-2*logLik(sp_ml_4) + (length(coef(sp_ml_4)) + 3)*log(dim(DF)[1]),  
	-2*logLik(sp_ml_5) + (length(coef(sp_ml_5)) + 3)*log(dim(DF)[1]) )
		
# plot them
file_name = "SO4_AIC"
pdf(paste0(file_name,'.pdf'), width = 9, height = 9)
	old.par = par(mar = c(5,5,1,1))
	plot(0:5, AIC_indep, ylim = c(230, 500), type = 'l', lwd = 3, 
		xlab = 'Order of Polynomial', ylab = '-2loglikelihood or AIC or BIC', cex.axis = 1.5, cex.lab = 2)
	points(0:5, AIC_indep, pch = 19, cex = 3)
	lines(0:5, m2ll_indep,lty = 1, lwd = 3)
	points(0:5, m2ll_indep, pch = 17, cex = 3)
	lines(0:5, BIC_indep,lty = 1, lwd = 3)
	points(0:5, BIC_indep, pch = 15, cex = 3)
	lines(0:5, AIC_spatial,lty = 2, lwd = 3)
	points(0:5, AIC_spatial, pch = 21, cex = 3)
	lines(0:5, m2ll_spatial,lty = 2, lwd = 3)
	points(0:5, m2ll_spatial, pch = 24, cex = 3)
	lines(0:5, BIC_spatial,lty = 2, lwd = 3)
	points(0:5, BIC_spatial, pch = 22, cex = 3)
	legend(2.9, 500, legend = 
		c('Indep m2LL','Indep AIC','Indep BIC',
		'Spatial m2LL', 'Spatial AIC', 'Spatial BIC'), lty = c(1,1,1,2,2,2), 
		lwd = 2, pch = c(17, 19, 15, 24, 21, 22), cex = 2)
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
#        Visualize Profiled Likelihood for Autocorrelation Parameters
#-------------------------------------------------------------------------------
################################################################################

#set grids sequences and dimensions for all plots
seqs = seq(-1,1,by = 0.05)
grid_dim = length(seqs)

# Model A
corModel = 'circular'

# inital covariance values using REMLE
splmfitA = splm(y ~ 1, 
	data = DF, xcoord = easting, ycoord = northing, 
	spcov_type = corModel, estmethod = 'reml')
summary(splmfitA)
# create a sequence of values to try for range and partial sill
de_seqA = log(coef(splmfitA, type = "spcov")['de']) + 3*seqs
range_seqA = log(coef(splmfitA, type = "spcov")['range']) + 3*seqs
# compute minus 2 times loglikelihood on a grid around REMLE estimates
# empty matrix to store results
zA = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		splmout = splm(
			y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, estmethod = 'reml',
			spcov_initial = spcov_initial(spcov_type = corModel, 
				de = exp(de_seqA[i]), range = exp(range_seqA[j]), 
				ie = coef(splmfitA, type = "spcov")['ie'],
				known = c('de','range','ie'))
		)	
		zA[i,j] = logLik(splmout)
	}
}

	nbrks = 20
	brks = quantile(zA, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(de_seqA, range_seqA, zA, breaks = brks, col = cramp,
		cex.main = 2, xlab = 'log(partial sill)', ylab = 'log(range)',
		cex.axis = 1.5, cex.lab = 2)
	points(log(coef(splmfitA, type = "spcov")['de']), 
		log(coef(splmfitA, type = "spcov")['range']), 
		pch = 19, cex = 2, col = 'white')

# Model B 
corModel = 'spherical'

# inital covariance values using REMLE
splmfitB = splm(y ~ 1, 
	data = DF, xcoord = easting, ycoord = northing, 
	spcov_type = corModel, estmethod = 'reml')
summary(splmfitB)

# create a sequence of values to try for range and partial sill

de_seqB = log(coef(splmfitB, type = "spcov")['de']) + 3*seqs
range_seqB = log(coef(splmfitB, type = "spcov")['range']) + 3*seqs
# compute minus 2 times loglikelihood on a grid around REMLE estimates
# empty matrix to store results
zB = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		splmout = splm(
			y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, estmethod = 'reml',
			spcov_initial = spcov_initial(spcov_type = corModel, 
				de = exp(de_seqB[i]), range = exp(range_seqB[j]), 
				ie = coef(splmfitB, type = "spcov")['ie'],
				known = c('de','range','ie'))
		)	
		zB[i,j] = logLik(splmout)
	}
}


# Model C
corModel = 'gaussian'

# inital covariance values using REMLE
splmfitC = splm(y ~ 1, 
	data = DF, xcoord = easting, ycoord = northing, 
	spcov_type = corModel, estmethod = 'reml')

# create a sequence of values to try for range and partial sill
de_seqC = log(coef(splmfitC, type = "spcov")['de']) + 3*seqs
range_seqC = log(coef(splmfitC, type = "spcov")['range']) + seqs
# compute minus 2 times loglikelihood on a grid around REMLE estimates
# empty matrix to store results
zC = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		splmout = splm(
			y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, estmethod = 'reml',
			spcov_initial = spcov_initial(spcov_type = corModel, 
				de = exp(de_seqC[i]), range = exp(range_seqC[j]), 
				ie = coef(splmfitC, type = "spcov")['ie'],
				known = c('de','range','ie'))
		)	
		zC[i,j] = logLik(splmout)
	}
}


# Model D
corModel = 'gravity'

# inital covariance values using REMLE
splmfitD = splm(y ~ 1, 
	data = DF, xcoord = easting, ycoord = northing, 
	spcov_type = corModel, estmethod = 'reml')

# create a sequence of values to try for range and partial sill
de_seqD = log(coef(splmfitD, type = "spcov")['de']) + 3*seqs
range_seqD = log(coef(splmfitD, type = "spcov")['range']) + 1.5*seqs
# compute minus 2 times loglikelihood on a grid around REMLE estimates
# empty matrix to store results
zD = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		splmout = splm(
			y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, estmethod = 'reml',
			spcov_initial = spcov_initial(spcov_type = corModel, 
				de = exp(de_seqD[i]), range = exp(range_seqD[j]), 
				ie = coef(splmfitD, type = "spcov")['ie'],
				known = c('de','range','ie'))
		)	
		zD[i,j] = logLik(splmout)
	}
}

file_name = "SO4_Viz_m2LL_covParms"

#pdf(paste0(file_name,'.pdf'), width = 12.5, height = 12.5)
tiff(paste0(file_name,'.tiff'), width = 720, height = 720)

	padj = -.5
	adj = -.15
	dotcol = 'black'
	layout(matrix(1:4, ncol = 2, byrow = TRUE))
	old.par = par(mar = c(5,5,5,1))

	nbrks = 20
	brks = quantile(zA, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(de_seqA, range_seqA, zA, breaks = brks, col = cramp,
		cex.main = 2, xlab = 'log(partial sill)', ylab = 'log(range)',
		cex.axis = 1.5, cex.lab = 2)
	points(log(coef(splmfitA, type = "spcov")['de']), 
		log(coef(splmfitA, type = "spcov")['range']), 
		pch = 19, cex = 2, col = dotcol)
	mtext('A', adj = adj, cex = 3, padj = padj)

	brks = quantile(zB, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(de_seqB, range_seqB, zB, breaks = brks, col = cramp,
		cex.main = 2, xlab = 'log(partial sill)', ylab = 'log(range)',
		cex.axis = 1.5, cex.lab = 2)
	points(log(coef(splmfitB, type = "spcov")['de']), 
		log(coef(splmfitB, type = "spcov")['range']), 
		pch = 19, cex = 2, col = dotcol)
	mtext('B', adj = adj, cex = 3, padj = padj)

	brks = quantile(zC, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(de_seqC, range_seqC, zC, breaks = brks, col = cramp,
		cex.main = 2, xlab = 'log(partial sill)', ylab = 'log(range)',
		cex.axis = 1.5, cex.lab = 2)
	points(log(coef(splmfitC, type = "spcov")['de']), 
		log(coef(splmfitC, type = "spcov")['range']), 
		pch = 19, cex = 2, col = dotcol)
	mtext('C', adj = adj, cex = 3, padj = padj)

	brks = quantile(zD, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(de_seqD, range_seqD, zD, breaks = brks, col = cramp,
		cex.main = 2, xlab = 'log(partial sill)', ylab = 'log(range)',
		cex.axis = 1.5, cex.lab = 2)
	points(log(coef(splmfitD, type = "spcov")['de']), 
		log(coef(splmfitD, type = "spcov")['range']), 
		pch = 19, cex = 2, col = dotcol)
	mtext('D', adj = adj, cex = 3, padj = padj)

	par(old.par)
	
dev.off()

system(paste0('tiff2pdf -o','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.tiff','\''))

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

# create a sequence of values to try for range and partial sill

de_seqE = log(coef(splmfitB, type = "spcov")['de']) + seqs
ie_seqE = log(coef(splmfitB, type = "spcov")['ie']) + seqs
# compute minus 2 times loglikelihood on a grid around REMLE estimates
# empty matrix to store results
zE = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		splmout = splm(
			y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, estmethod = 'reml',
			spcov_initial = spcov_initial(spcov_type = 'spherical', 
				de = exp(de_seqE[i]), ie = exp(ie_seqE[j]), 
				range = coef(splmfitB, type = "spcov")['range'],
				known = c('de','range','ie'))
		)	
		zE[i,j] = logLik(splmout)
	}
}

file_name = "SO4_psillvsnugget"

#pdf(paste0(file_name,'.pdf'), width = 8, height = 8)
tiff(paste0(file_name,'.tiff'), width = 480, height = 480)

	old.par = par(mar = c(5,5,1,1))
	
	brks = quantile(zE, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(de_seqE, ie_seqE, zE, breaks = brks, col = cramp,
		cex.main = 2, xlab = 'log(partial sill)', ylab = 'log(nugget)',
		cex.axis = 1.5, cex.lab = 2)
	points(log(coef(splmfitB, type = "spcov")['de']), 
		log(coef(splmfitB, type = "spcov")['ie']), 
		pch = 19, cex = 2, col = dotcol)

	par(old.par)
	
dev.off()

system(paste0('tiff2pdf -o','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.tiff','\''))

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

################################################################################
#-------------------------------------------------------------------------------
# Hessian, Fisher Information, and Covariance Matrix for Covariance Parameters
#-------------------------------------------------------------------------------
################################################################################

asy_se = data.frame(NA, nrow = 7, ncol = 10)

spcov_list = c('exponential', 'spherical','gaussian','gravity',
	'rquad','magnetic','circular')
for(i in 1:7) {
	# pick a model and an estimation method
	spcov_mod = spcov_list[i]
	pvec_len = 3
	estmeth = 'reml'

	#fit the model to obtain the (RE)MLE
	ModelFit = splm(y ~ 1, 
		data = DF, xcoord = easting, ycoord = northing, 
		spcov_type = spcov_mod, estmethod = estmeth,
		control = list(reltol = 1e-7))
	summary(ModelFit)
	theta_reml = coef(ModelFit, type = 'spcov')[1:pvec_len]
	theta = theta_reml

	# ---------------------- Numerical Hessian Approach ----------------------------

	# create a function that returns the log-likelihood at any value of the
	# covariance parameters
	REML_parms_only = function(theta)
	{
		#initialize covariance parameters and hold them constant
		if(pvec_len == 3) {
			spini = spcov_initial = spcov_initial(spcov_type = spcov_mod, de = theta[1], 
				ie = theta[2], range = theta[3], known = c('de','ie','range'))
		} else {
			spini = spcov_initial = spcov_initial(spcov_type = spcov_mod, de = theta[1], 
				ie = theta[2], range = theta[3], extra = theta[4],
					known = c('de','ie','range','extra'))
		}
		# return the log-likelihood for the specified theta value
		ModelFit = splm(y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, 
			estmethod = estmeth, spcov_initial = spini )
			logLik(ModelFit)
	}
	# use a numerical method to find the Hessian matrix
	# observed Fisher Information is the negative Hessian evaluated at (RE)MLE
	FishInf = -hessian(REML_parms_only, theta_reml)
	# asymptotic covariance matrix
	vcov_theta_comp = solve(FishInf)

	# ---------------------- Analytical Approach ----------------------------

	# get the matrix of all pairwise distances
	D = as.matrix(dist(cbind(DF$easting,DF$northing)))
	# create a diagonal matrix with dimension D
	II = diag(dim(D)[1])
	# create a matrix of all 1's with dimension D
	allones = matrix(1, nrow = dim(D)[1], ncol = dim(D)[1])
	# for models with restricted range
	Dind = D*(D < theta['range'])
	# create the spatial correlation matrix
	if(spcov_mod == 'exponential') R = exp(-D/theta['range'])
	if(spcov_mod == 'spherical') R = (allones - 1.5*Dind/theta['range'] + 
		0.5*(Dind/theta['range'])^3)
	if(spcov_mod == 'gaussian') R = exp(-D/(theta['range'])^2)
	if(spcov_mod == 'gravity') R = (allones + (D/theta['range'])^2)^(-0.5)
	if(spcov_mod == 'rquad') R = (allones + (D/theta['range'])^2)^(-1)
	if(spcov_mod == 'magnetic') R = (allones + (D/theta['range'])^2)^(-1.5)
	if(spcov_mod == 'circular')  R = 2/pi*(acos(Dind/theta['range']) - 
		Dind/theta['range']*sqrt(1 - Dind^2/theta['range']^2))

	# create the covariance matrix and take its inverse
	V = theta['de']*R + theta['ie']*II
	Vi = solve(V)

	# create the derivative matrix
	if(spcov_mod == 'exponential') A = theta['de']*D*R/theta['range']^2
	if(spcov_mod == 'spherical') A = -theta['de']*
		1.5*(Dind^3 - Dind*theta['range']^2)/(theta['range']^4)
	if(spcov_mod == 'gaussian') A = theta['de']*2*D^2*R/theta['range']^3
	if(spcov_mod == 'gravity') A = theta['de']*D^2/(theta['range']^3*
		(D^2/theta['range']^2 + 1)^(1.5))
	if(spcov_mod == 'rquad') A = theta['de']*(2*D^2*theta['range'])/((theta['range']^2 +
		(D^2))^(2))
	if(spcov_mod == 'magnetic') A = theta['de']*3*D^2/(theta['range']^3*
		(D^2/theta['range']^2 + 1)^(2.5))
	if(spcov_mod == 'circular') A = theta['de']*2/pi*2*Dind/theta['range']^2*
		sqrt(1 - Dind^2/theta['range']^2)

	# compute Fisher Information element by element using trace formula
	FI = matrix(rep(NA, times = 9), nrow = 3)
	FI[1,1] = sum(diag(Vi %*% R %*% Vi %*% R))
	FI[1,2] = FI[2,1] = sum(diag(Vi %*% R %*% Vi %*% II))
	FI[1,3] = FI[3,1] = sum(diag(Vi %*% R %*% Vi %*% A))
	FI[2,2] = sum(diag(Vi %*% II %*% Vi %*% II))
	FI[2,3] = FI[3,2] = sum(diag(Vi %*% II %*% Vi %*% A))
	FI[3,3] = sum(diag(Vi %*% A %*% Vi %*% A))
	FI = 0.5*FI

	# get the variance covariance matrix of theta from Fisher Information
	vcov_theta_ana = solve(FI)

	asy_se[i,1] = spcov_mod
	asy_se[i,2:4] = theta
	asy_se[i,5:7] = sqrt(diag(vcov_theta_ana))
	asy_se[i,8:10] = sqrt(diag(vcov_theta_comp))
}

print(
    xtable(asy_se, 
      align = c('l',rep('l', times = length(asy_se[1,]))),
      digits = c(0,0,rep(3, times = 9)),
      caption = 'Aymptotic Standard Errors for Covariance Parameters',
      label = 'tab:AsySE'
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# show the asymptotic correlation matrix for the last model used, the
# circular model
print(
	xtable(
		diag(1/sqrt(diag(vcov_theta_ana))) %*% vcov_theta_ana %*% 
				diag(1/sqrt(diag(vcov_theta_ana))),
     align = c('l',rep('l', times = 3)),
      digits = c(0,rep(4, times = 3)),
      caption = 'Aymptotic correlation matrix for covariance parameters',
      label = 'tab:AsySE'
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
#                            The Matern Model
#-------------------------------------------------------------------------------
################################################################################

seqs = seq(-1,1,by = 0.05)
grid_dim = length(seqs)


corModel = 'exponential'

# inital covariance values using REMLE
splmfitA = splm(y ~ 1, 
	data = DF, xcoord = easting, ycoord = northing, 
	spcov_type = corModel, estmethod = 'reml',
	spcov_initial(spcov_type = corModel, 
				de = exp(.1), 
				range = exp(1)))
summary(splmfitA)
# create a sequence of values to try for range and partial sill
de_seqA = log(coef(splmfitA, type = "spcov")['de']) + 3*seqs
range_seqA = log(coef(splmfitA, type = "spcov")['range']) + 3*seqs
# compute minus 2 times loglikelihood on a grid around REMLE estimates
# empty matrix to store results
zA = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		splmout = splm(
			y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, estmethod = 'reml',
			spcov_initial = spcov_initial(spcov_type = corModel, 
				de = exp(de_seqA[i]), range = exp(range_seqA[j]), 
				ie = coef(splmfitA, type = "spcov")['ie'],
				known = c('de','range','ie'))
		)	
		zA[i,j] = logLik(splmout)
	}
}

	nbrks = 20
	brks = quantile(zA, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(de_seqA, range_seqA, zA, breaks = brks, col = cramp,
		cex.main = 2, xlab = 'log(partial sill)', ylab = 'log(range)',
		cex.axis = 1.5, cex.lab = 2)
	points(log(coef(splmfitA, type = "spcov")['de']), 
		log(coef(splmfitA, type = "spcov")['range']), 
		pch = 19, cex = 2, col = 'white')
















loocv(ModelFit)
seqs = seq(-1,1,by = 0.05)
grid_dim = length(seqs)

seqs = seq(-1,1,by = 0.1)
range_seqs = exp(log(coef(MaternFit, type = "spcov")['range']) + seqs)
seqs_len = length(seqs)
store_range = rep(NA, times = seqs_len)
for(i in 1:seqs_len) {
		splmout = splm(
			y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, estmethod = 'reml',
			spcov_initial = spcov_initial(spcov_type = 'matern', 
#				de = coef(MaternFit, type = "spcov")['de'], 
				range = range_seqs[i], 
#				ie = coef(MaternFit, type = "spcov")['ie'],
#				extra = coef(MaternFit, type = "spcov")['extra'],
				known = c('range')),
				control = list(reltol = 1e-7)
		)	
	store_range[i] = logLik(splmout)
	}
old.par = par(mar = c(5,5,1,1))
plot(range_seqs,store_range, pch = 1, cex = 2, xlab = 'Range Parameter Values', 
	ylab = 'Profile Likelihood Value', cex.lab = 2, cex.axis = 1.5)
points(coef(MaternFit, type = "spcov")['range'], logLik(MaternFit), pch = 3,
	cex = 3, lwd = 3)
par(old.par)

range_seqs[which(store_range == max(store_range))]

seqs = seq(-1,1,by = 0.1)
de_seqs = exp(log(coef(MaternFit, type = "spcov")['de']) + seqs)
seqs_len = length(seqs)
store_de = rep(NA, times = seqs_len)
for(i in 1:seqs_len) {
		splmout = splm(
			y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, estmethod = 'reml',
			spcov_initial = spcov_initial(spcov_type = 'matern', 
				de = de_seqs[i], 
#				range = coef(MaternFit, type = "spcov")['range'], 
#				ie = coef(MaternFit, type = "spcov")['ie'],
#				extra = coef(MaternFit, type = "spcov")['extra'],
				known = c('de')),
				control = list(reltol = 1e-7)
		)	
	store_de[i] = logLik(splmout)
	}
old.par = par(mar = c(5,5,1,1))
plot(de_seqs,store_de, pch = 1, cex = 2, xlab = 'Partial Sill', 
	ylab = 'Profile Likelihood', cex.lab = 2, cex.axis = 1.5)
points(coef(MaternFit, type = "spcov")['de'], logLik(MaternFit), pch = 3,
	cex = 3, lwd = 3)
par(old.par)

de_seqs[which(store_de == max(store_de))]

seqs = seq(-1,1,by = 0.1)
extra_seqs = exp(1.6*seqs)
seqs_len = length(seqs)
store_extra = rep(NA, times = seqs_len)
for(i in 1:seqs_len) {
		splmout = splm(
			y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, estmethod = 'reml',
			spcov_initial = spcov_initial(spcov_type = 'matern', 
#				de = coef(MaternFit, type = "spcov")['de'], 
#				range = coef(MaternFit, type = "spcov")['range'], 
#				ie = coef(MaternFit, type = "spcov")['ie'],
				extra = extra_seqs[i],
				known = c('extra'))
		)	
	store_extra[i] = logLik(splmout)
	}
old.par = par(mar = c(5,5,1,1))
plot(extra_seqs,store_extra, pch = 1, cex = 2, xlab = 'Smoothness Parameter', 
	ylab = 'Profile Likelihood', cex.lab = 2, cex.axis = 1.5)
points(coef(MaternFit, type = "spcov")['extra'], logLik(MaternFit), pch = 3,
	cex = 3, lwd = 3)
par(old.par)

seqs = seq(-1,1,by = 0.1)
ie_seqs = exp(log(coef(MaternFit, type = "spcov")['ie']) + seqs)
seqs_len = length(seqs)
store_ie = rep(NA, times = seqs_len)
i = 6
for(i in 1:seqs_len) {
		splmout = splm(
			y ~ 1, 
			data = DF, xcoord = easting, ycoord = northing, estmethod = 'reml',
			spcov_initial = spcov_initial(spcov_type = 'matern', 
				ie = ie_seqs[i], known = c('ie'))
		)	
	store_ie[i] = logLik(splmout)
	}
old.par = par(mar = c(5,5,1,1))
plot(ie_seqs,store_ie, pch = 1, cex = 2, xlab = 'Nugget', 
	ylab = 'Profile Likelihood Value', cex.lab = 2, cex.axis = 1.5)
points(coef(MaternFit, type = "spcov")['ie'], logLik(MaternFit), pch = 3,
	cex = 3, lwd = 3)
par(old.par)

de_seqs[which(store_de == max(store_de))]

# profile likelihood for range parameter

file_name = "Matern_SO4"

pdf(paste0(file_name,'.pdf'), width = 9, height = 9)

	layout(matrix(1:4, ncol = 2, byrow = TRUE))

	padj = -.6
	adj = -.25
	old.par = par(mar = c(5,5,4,1))
	plot(de_seqs,store_de, pch = 1, cex = 2, xlab = 'Partial Sill', 
		ylab = 'Profile Likelihood', cex.lab = 2, cex.axis = 1.5)
	points(coef(MaternFit, type = "spcov")['de'], logLik(MaternFit), pch = 3,
		cex = 3, lwd = 3)
	mtext('A', adj = adj, cex = 3, padj = padj)
	plot(ie_seqs,store_ie, pch = 1, cex = 2, xlab = 'Nugget', 
		ylab = 'Profile Likelihood Value', cex.lab = 2, cex.axis = 1.5)
	points(coef(MaternFit, type = "spcov")['ie'], logLik(MaternFit), pch = 3,
		cex = 3, lwd = 3)
	mtext('B', adj = adj, cex = 3, padj = padj)
	plot(range_seqs,store_range, pch = 1, cex = 2, xlab = 'Range', 
		ylab = 'Profile Likelihood', cex.lab = 2, cex.axis = 1.5)
	points(coef(MaternFit, type = "spcov")['range'], logLik(MaternFit), pch = 3,
		cex = 3, lwd = 3)
	mtext('C', adj = adj, cex = 3, padj = padj)
	plot(extra_seqs,store_extra, pch = 1, cex = 2, xlab = 'Smoothness', 
		ylab = 'Profile Likelihood', cex.lab = 2, cex.axis = 1.5)
	points(coef(MaternFit, type = "spcov")['extra'], logLik(MaternFit), pch = 3,
		cex = 3, lwd = 3)
	mtext('D', adj = adj, cex = 3, padj = padj)
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
#       Mapping Spatial versus High-Order-Polynomial-Indpendence Models
#-------------------------------------------------------------------------------
################################################################################

data(USboundary)
USboundary@bbox[1,2] - USboundary@bbox[1,1]
USboundary@bbox[2,2] - USboundary@bbox[2,1]
# spacing
spacing = (USboundary@bbox[1,2] - USboundary@bbox[1,1])/100
xpredlocs = USboundary@bbox[1,1] + spacing*(1:100 - 0.5)
ypredlocs = USboundary@bbox[2,1] + spacing*(1:100 - 0.5)
# create a data frame of predictions covering all of US
preds = data.frame(x = xpredlocs %x% rep(1, times = 100), 
	y = rep(1, times = 100) %x% ypredlocs)
# turn it into a spatial points data frame
coordinates(preds) <- ~ x + y
# give it the same projection as US polygons
preds@proj4string = USboundary@proj4string
# clip the points to be within US borders
preds_sub = preds[USboundary]
# check the result
plot(USboundary)
plot(preds_sub, add = TRUE)
# get the coordinates of the clipped points
DFpred = as.data.frame(preds_sub@coords)
colnames(DFpred) = c('easting', 'northing')
# change coordinates to 1000 km
DFpred = DFpred/1e+6

#         Independence Model

# predict points in DFpred using independence model with 3th order polynomial
pred_ind = predict(lm_3, DFpred, se.fit = TRUE)
# predict points with constant mean spherical model
pred_sp1 = predict(splmfitB, DFpred, se.fit = TRUE)
# predict points with constant mean gaussian model
pred_sp2 = predict(splmfitC, DFpred, se.fit = TRUE)
# fit 3rd order polynomial with exponential model
# inital covariance values using REMLE
ExpoFit = splm(y ~ poly(easting, northing, degree = 3, raw = TRUE), 
	data = DF, xcoord = easting, ycoord = northing, 
	spcov_type = 'exponential', estmethod = 'reml',
	control = list(reltol = 1e-7))
summary(ExpoFit)
pred_sp3 = predict(ExpoFit, DFpred, se.fit = TRUE)

file_name = "SO4_Prediction_Maps"

pdf(paste0(file_name,'.pdf'), width = 17, height = 24)
brks_cex = 3
printF = "1.3"

brks_pred = quantile(c(pred_ind$fit, pred_sp1$fit, 
	pred_sp2$fit, pred_sp3$fit), probs = seq(0, 1, .1))
brks_se = quantile(c(pred_ind$se.fit, pred_sp1$se.fit, 
	pred_sp2$se.fit, pred_sp3$se.fit), probs = seq(0, 1, .1))
# make a color map of predictions
cip = classIntervals(pred_ind$fit, style = 'fixed', fixedBreaks= brks_pred)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
source('addBreakColorLegend.R')
layout(matrix(1:16, ncol = 4, byrow = TRUE), 
	widths = c(3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
text(-2350000, 3200000, 'A', cex = 6)
old.par = par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

# make a color map of prediction standard errors
cip = classIntervals(pred_ind$se.fit, style = 'fixed', fixedBreaks= brks_se)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
points(coordinates(SO4clean), pch = 19)
text(-2350000, 3200000, 'B', cex = 6)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

#         Spatial Models

cip = classIntervals(pred_sp1$fit, style = 'fixed', fixedBreaks= brks_pred)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
text(-2350000, 3200000, 'C', cex = 6)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

cip = classIntervals(pred_sp1$se.fit, style = 'fixed', fixedBreaks= brks_se)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
points(coordinates(SO4clean), pch = 19)
text(-2350000, 3200000, 'D', cex = 6)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

cip = classIntervals(pred_sp2$fit, style = 'fixed', fixedBreaks= brks_pred)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
text(-2350000, 3200000, 'E', cex = 6)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

cip = classIntervals(pred_sp2$se.fit, style = 'fixed', fixedBreaks= brks_se)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
points(coordinates(SO4clean), pch = 19)
text(-2350000, 3200000, 'F', cex = 6)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

cip = classIntervals(pred_sp3$fit, style = 'fixed', fixedBreaks= brks_pred)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
text(-2350000, 3200000, 'G', cex = 6)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = printF)

cip = classIntervals(pred_sp3$se.fit, style = 'fixed', fixedBreaks= brks_se)
palp = viridis(9)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
points(coordinates(SO4clean), pch = 19)
text(-2350000, 3200000, 'H', cex = 6)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
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
