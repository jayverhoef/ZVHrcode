# Set a path as a working directory
sec_path = 'Rcode/Chapter8/Section 8.6/'
setwd(paste0(SLEDbook_path,sec_path))

# This code computes entries for the Tables 8.5, 8.6, and 8.7 for the "exp" fitted model and all three true models 

library(xtable)
library(spmodel)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#             Code for Tables 8.5, 8.6, and 8.7
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 12
locxy <- data.frame(xcoord = rep(1:K,K), ycoord = rep(1:K, each = K))
n <- K^2

range_sph = 2.88
range_exp = 1.44
range_gau = sqrt(1.44)

# Compute distances between sites on the grid
d = as.matrix(dist(locxy))

# true covariance parameters for simulation
spcov_params_sph <- spcov_params("spherical", de = 1, ie = 0.000001, 
	range = range_sph)
spcov_params_exp <- spcov_params("exponential", de = 1, ie = 0.000001, 
	range = range_exp)
spcov_params_gau <- spcov_params("gaussian", de = 1, ie = 0.000001, 
	range = range_gau)
simy = sprnorm(spcov_params_val, mean = 0, data = locxy, 
	xcoord = xcoord, ycoord = ycoord)
	
cor_sph = function(h, range) {1-1.5*(h/range)+0.5*(h/range)^3} 
cor_exp = function(h, range) {exp(-h/range)}
cor_gau = function(h, range) {exp(-(h/range)^2)}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                       Simulations
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

set.seed(507)
lengthCIbeta = matrix(0, nrow = 12, ncol = 3)
CI90 = matrix(0, nrow = 12, ncol = 3)
biassig2 = matrix(0, nrow = 9, ncol = 3)
biasUDC = matrix(0, nrow = 9, ncol = 3)
MSEsig2 = matrix(0, nrow = 9, ncol = 3)
MSEUDC = matrix(0, nrow = 9, ncol = 3)
niter = 1000
for(k in 1:niter){
	cat("\r", "iteration: ", k)
	# simulate data, all true regression coefficients are zero
	# use sprnorm from the spmodel package
	simsph = sprnorm(spcov_params_sph, mean = 0, data = locxy, 
		xcoord = xcoord, ycoord = ycoord)
	simexp = sprnorm(spcov_params_exp, mean = 0, data = locxy, 
		xcoord = xcoord, ycoord = ycoord)
	simgau = sprnorm(spcov_params_gau, mean = 0, data = locxy, 
		xcoord = xcoord, ycoord = ycoord)
	# create data.frame for spmodel
	DF = data.frame(simsph = simsph, simexp = simexp, simgau = simgau,
		xcoord = locxy[,1], ycoord = locxy[,2])
	# ------- ML ESTIMATION
	# fit model using spmodel, holding nugget effect at 0.000001
	# in the name, true model is first, fitted model is second
	# sph = spherical, exp = exponential, gau = gaussian, unc = uncorrelated
	               # constant mean
	sphsph_const = splm(simsph ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	sphexp_const = splm(simsph ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	sphgau_const = splm(simsph ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	sphunc_const = lm(simsph ~ 1, data = DF)
	expsph_const = splm(simexp ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	expexp_const = splm(simexp ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	expgau_const = splm(simexp ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	expunc_const = lm(simexp ~ 1, data = DF)
	gausph_const = splm(simgau ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	gauexp_const = splm(simgau ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	gaugau_const = splm(simgau ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	gauunc_const = lm(simgau ~ 1, data = DF)
	               # row effects
	sphsph_row = splm(simsph ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	sphexp_row = splm(simsph ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	sphgau_row = splm(simsph ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	sphunc_row = lm(simsph ~ I(as.factor(xcoord)), data = DF)
	expsph_row = splm(simexp ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	expexp_row = splm(simexp ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	expgau_row = splm(simexp ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	expunc_row = lm(simexp ~ I(as.factor(xcoord)), data = DF)
	gausph_row = splm(simgau ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	gauexp_row = splm(simgau ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	gaugau_row = splm(simgau ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	gauunc_row = lm(simgau ~ I(as.factor(xcoord)), data = DF)
	               # row and column effects
	sphsph_rowcol = splm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)),
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	sphexp_rowcol = splm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	sphgau_rowcol = splm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	sphunc_rowcol = lm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF)
	expsph_rowcol = splm(simexp ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	expexp_rowcol = splm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	expgau_rowcol = splm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	expunc_rowcol = lm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF)
	gausph_rowcol = splm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	gauexp_rowcol = splm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	gaugau_rowcol = splm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'ml')
	gauunc_rowcol = lm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF)

	# ------- REML ESTIMATION
	# fit model using spmodel, holding nugget effect at 0.000001
	# in the name, true model is first, fitted model is second
	# sph = spherical, exp = exponential, gau = gaussian, unc = uncorrelated
	               # constant mean
	sphsph_const_REML = splm(simsph ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphexp_const_REML = splm(simsph ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphgau_const_REML = splm(simsph ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphunc_const_REML = lm(simsph ~ 1, data = DF)
	expsph_const_REML = splm(simexp ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expexp_const_REML = splm(simexp ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expgau_const_REML = splm(simexp ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expunc_const_REML = lm(simexp ~ 1, data = DF)
	gausph_const_REML = splm(simgau ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gauexp_const_REML = splm(simgau ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gaugau_const_REML = splm(simgau ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gauunc_const_REML = lm(simgau ~ 1, data = DF)
	               # row effects
	sphsph_row_REML = splm(simsph ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphexp_row_REML = splm(simsph ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphgau_row_REML = splm(simsph ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphunc_row_REML = lm(simsph ~ I(as.factor(xcoord)), data = DF)
	expsph_row_REML = splm(simexp ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expexp_row_REML = splm(simexp ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expgau_row_REML = splm(simexp ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expunc_row_REML = lm(simexp ~ I(as.factor(xcoord)), data = DF)
	gausph_row_REML = splm(simgau ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gauexp_row_REML = splm(simgau ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gaugau_row_REML = splm(simgau ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gauunc_row_REML = lm(simgau ~ I(as.factor(xcoord)), data = DF)
	               # row and column effects
	sphsph_rowcol_REML = splm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)),
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphexp_rowcol_REML = splm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)),
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphgau_rowcol_REML = splm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)),
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphunc_rowcol_REML = lm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF)
	expsph_rowcol_REML = splm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expexp_rowcol_REML = splm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expgau_rowcol_REML = splm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expunc_rowcol_REML = lm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF)
	gausph_rowcol_REML = splm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gauexp_rowcol_REML = splm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gaugau_rowcol_REML = splm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gauunc_rowcol_REML = lm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), 
		data = DF)

	#
	#                           estimation summaries
	#
	# length of confidence intervals for beta
			#constant model
	lengthCIbeta[1,1] = lengthCIbeta[1,1] + mean(summary(sphsph_const_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[1,2] = lengthCIbeta[1,2] + mean(summary(expsph_const_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[1,3] = lengthCIbeta[1,3] + mean(summary(gausph_const_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[2,1] = lengthCIbeta[2,1] + mean(summary(sphexp_const_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[2,2] = lengthCIbeta[2,2] + mean(summary(expexp_const_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[2,3] = lengthCIbeta[2,3] + mean(summary(gauexp_const_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[3,1] = lengthCIbeta[3,1] + mean(summary(sphgau_const_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[3,2] = lengthCIbeta[3,2] + mean(summary(expgau_const_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[3,3] = lengthCIbeta[3,3] + mean(summary(gaugau_const_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[4,1] = lengthCIbeta[4,1] + mean(summary(sphunc_const_REML)$
		coefficients[1,'Std. Error']*2*1.645)
	lengthCIbeta[4,2] = lengthCIbeta[4,2] + mean(summary(expunc_const_REML)$
		coefficients[1,'Std. Error']*2*1.645)
	lengthCIbeta[4,3] = lengthCIbeta[4,3] + mean(summary(gauunc_const_REML)$
		coefficients[1,'Std. Error']*2*1.645)
			#row model
	lengthCIbeta[5,1] = lengthCIbeta[5,1] + mean(summary(sphsph_row_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[5,2] = lengthCIbeta[5,2] + mean(summary(expsph_row_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[5,3] = lengthCIbeta[5,3] + mean(summary(gausph_row_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[6,1] = lengthCIbeta[6,1] + mean(summary(sphexp_row_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[6,2] = lengthCIbeta[6,2] + mean(summary(expexp_row_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[6,3] = lengthCIbeta[6,3] + mean(summary(gauexp_row_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[7,1] = lengthCIbeta[7,1] + mean(summary(sphgau_row_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[7,2] = lengthCIbeta[7,2] + mean(summary(expgau_row_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[7,3] = lengthCIbeta[7,3] + mean(summary(gaugau_row_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[8,1] = lengthCIbeta[8,1] + mean(summary(sphunc_row_REML)$
		coefficients[,'Std. Error']*2*1.645)
	lengthCIbeta[8,2] = lengthCIbeta[8,2] + mean(summary(expunc_row_REML)$
		coefficients[,'Std. Error']*2*1.645)
	lengthCIbeta[8,3] = lengthCIbeta[8,3] + mean(summary(gauunc_row_REML)$
		coefficients[,'Std. Error']*2*1.645)
			#row and column model
	lengthCIbeta[9,1] = lengthCIbeta[9,1] + mean(summary(sphsph_rowcol_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[9,2] = lengthCIbeta[9,2] + mean(summary(expsph_rowcol_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[9,3] = lengthCIbeta[9,3] + mean(summary(gausph_rowcol_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[10,1] = lengthCIbeta[10,1] + mean(summary(sphexp_rowcol_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[10,2] = lengthCIbeta[10,2] + mean(summary(expexp_rowcol_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[10,3] = lengthCIbeta[10,3] + mean(summary(gauexp_rowcol_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[11,1] = lengthCIbeta[11,1] + mean(summary(sphgau_rowcol_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[11,2] = lengthCIbeta[11,2] + mean(summary(expgau_rowcol_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[11,3] = lengthCIbeta[11,3] + mean(summary(gaugau_rowcol_REML)$
		coefficients$fixed[,'Std_Error']*2*1.645)
	lengthCIbeta[12,1] = lengthCIbeta[12,1] + mean(summary(sphunc_rowcol_REML)$
		coefficients[,'Std. Error']*2*1.645)
	lengthCIbeta[12,2] = lengthCIbeta[12,2] + mean(summary(expunc_rowcol_REML)$
		coefficients[,'Std. Error']*2*1.645)
	lengthCIbeta[12,3] = lengthCIbeta[12,3] + mean(summary(gauunc_rowcol_REML)$
		coefficients[,'Std. Error']*2*1.645)
	# 90% confidence interval coverage
			#constant model
	CI90[1,1] = CI90[1,1] + 
			mean(summary(sphsph_const_REML)$coefficients$fixed$p > 0.1)
	CI90[1,2] = CI90[1,2] + 
			mean(summary(expsph_const_REML)$coefficients$fixed$p > 0.1)
	CI90[1,3] = CI90[1,3] + 
			mean(summary(gausph_const_REML)$coefficients$fixed$p > 0.1)
	CI90[2,1] = CI90[2,1] + 
			mean(summary(sphexp_const_REML)$coefficients$fixed$p > 0.1)
	CI90[2,2] = CI90[2,2] + 
			mean(summary(expexp_const_REML)$coefficients$fixed$p > 0.1)
	CI90[2,3] = CI90[2,3] + 
			mean(summary(gauexp_const_REML)$coefficients$fixed$p > 0.1)
	CI90[3,1] = CI90[3,1] + 
			mean(summary(sphgau_const_REML)$coefficients$fixed$p > 0.1)
	CI90[3,2] = CI90[3,2] + 
			mean(summary(expgau_const_REML)$coefficients$fixed$p > 0.1)
	CI90[3,3] = CI90[3,3] + 
			mean(summary(gaugau_const_REML)$coefficients$fixed$p > 0.1)
	CI90[4,1] = CI90[4,1] + 
			mean(summary(sphunc_const_REML)$coefficients[1,'Pr(>|t|)'] > 0.1)
	CI90[4,2] = CI90[4,2] + 
			mean(summary(expunc_const_REML)$coefficients[1,'Pr(>|t|)'] > 0.1)
	CI90[4,3] = CI90[4,3] + 
			mean(summary(gauunc_const_REML)$coefficients[1,'Pr(>|t|)'] > 0.1)
			#row model
	CI90[5,1] = CI90[5,1] + 
			mean(summary(sphsph_row_REML)$coefficients$fixed$p > 0.1)
	CI90[5,2] = CI90[5,2] + 
			mean(summary(expsph_row_REML)$coefficients$fixed$p > 0.1)
	CI90[5,3] = CI90[5,3] + 
			mean(summary(gausph_row_REML)$coefficients$fixed$p > 0.1)
	CI90[6,1] = CI90[6,1] + 
			mean(summary(sphexp_row_REML)$coefficients$fixed$p > 0.1)
	CI90[6,2] = CI90[6,2] + 
			mean(summary(expexp_row_REML)$coefficients$fixed$p > 0.1)
	CI90[6,3] = CI90[6,3] + 
			mean(summary(gauexp_row_REML)$coefficients$fixed$p > 0.1)
	CI90[7,1] = CI90[7,1] + 
			mean(summary(sphgau_row_REML)$coefficients$fixed$p > 0.1)
	CI90[7,2] = CI90[7,2] + 
			mean(summary(expgau_row_REML)$coefficients$fixed$p > 0.1)
	CI90[7,3] = CI90[7,3] + 
			mean(summary(gaugau_row_REML)$coefficients$fixed$p > 0.1)
	CI90[8,1] = CI90[8,1] + 
			mean(summary(sphunc_row_REML)$coefficients[,'Pr(>|t|)'] > 0.1)
	CI90[8,2] = CI90[8,2] + 
			mean(summary(expunc_row_REML)$coefficients[,'Pr(>|t|)'] > 0.1)
	CI90[8,3] = CI90[8,3] + 
			mean(summary(gauunc_row_REML)$coefficients[,'Pr(>|t|)'] > 0.1)
			#rowcol model
	CI90[9,1] = CI90[9,1] + 
			mean(summary(sphsph_rowcol_REML)$coefficients$fixed$p > 0.1)
	CI90[9,2] = CI90[9,2] + 
			mean(summary(expsph_rowcol_REML)$coefficients$fixed$p > 0.1)
	CI90[9,3] = CI90[9,3] + 
			mean(summary(gausph_rowcol_REML)$coefficients$fixed$p > 0.1)
	CI90[10,1] = CI90[10,1] + 
			mean(summary(sphexp_rowcol_REML)$coefficients$fixed$p > 0.1)
	CI90[10,2] = CI90[10,2] + 
			mean(summary(expexp_rowcol_REML)$coefficients$fixed$p > 0.1)
	CI90[10,3] = CI90[10,3] + 
			mean(summary(gauexp_rowcol_REML)$coefficients$fixed$p > 0.1)
	CI90[11,1] = CI90[11,1] + 
			mean(summary(sphgau_rowcol_REML)$coefficients$fixed$p > 0.1)
	CI90[11,2] = CI90[11,2] + 
			mean(summary(expgau_rowcol_REML)$coefficients$fixed$p > 0.1)
	CI90[11,3] = CI90[11,3] + 
			mean(summary(gaugau_rowcol_REML)$coefficients$fixed$p > 0.1)
	CI90[12,1] = CI90[12,1] + 
			mean(summary(sphunc_rowcol_REML)$coefficients[,'Pr(>|t|)'] > 0.1)
	CI90[12,2] = CI90[12,2] + 
			mean(summary(expunc_rowcol_REML)$coefficients[,'Pr(>|t|)'] > 0.1)
	CI90[12,3] = CI90[12,3] + 
			mean(summary(gauunc_rowcol_REML)$coefficients[,'Pr(>|t|)'] > 0.1)

	# bias in estimation of overall variance parameter
			#constant model
	biassig2[1,1] = biassig2[1,1] + 
		coef(sphsph_const, type = 'spcov')['de'] - 1
	biassig2[1,2] = biassig2[1,2] + 
		coef(expsph_const, type = 'spcov')['de'] - 1
	biassig2[1,3] = biassig2[1,3] + 
		coef(gausph_const, type = 'spcov')['de'] - 1
	biassig2[2,1] = biassig2[2,1] + 
		coef(sphexp_const, type = 'spcov')['de'] - 1
	biassig2[2,2] = biassig2[2,2] + 
		coef(expexp_const, type = 'spcov')['de'] - 1
	biassig2[2,3] = biassig2[2,3] + 
		coef(gauexp_const, type = 'spcov')['de'] - 1
	biassig2[3,1] = biassig2[3,1] + 
		coef(sphgau_const, type = 'spcov')['de'] - 1
	biassig2[3,2] = biassig2[3,2] + 
		coef(expgau_const, type = 'spcov')['de'] - 1
	biassig2[3,3] = biassig2[3,3] + 
		coef(gaugau_const, type = 'spcov')['de'] - 1
			#row model
	biassig2[4,1] = biassig2[4,1] + 
		coef(sphsph_row, type = 'spcov')['de'] - 1
	biassig2[4,2] = biassig2[4,2] + 
		coef(expsph_row, type = 'spcov')['de'] - 1
	biassig2[4,3] = biassig2[4,3] + 
		coef(gausph_row, type = 'spcov')['de'] - 1
	biassig2[5,1] = biassig2[5,1] + 
		coef(sphexp_row, type = 'spcov')['de'] - 1
	biassig2[5,2] = biassig2[5,2] + 
		coef(expexp_row, type = 'spcov')['de'] - 1
	biassig2[5,3] = biassig2[5,3] + 
		coef(gauexp_row, type = 'spcov')['de'] - 1
	biassig2[6,1] = biassig2[6,1] + 
		coef(sphgau_row, type = 'spcov')['de'] - 1
	biassig2[6,2] = biassig2[6,2] + 
		coef(expgau_row, type = 'spcov')['de'] - 1
	biassig2[6,3] = biassig2[6,3] + 
		coef(gaugau_row, type = 'spcov')['de'] - 1
			#row and column model
	biassig2[7,1] = biassig2[7,1] + 
		coef(sphsph_rowcol, type = 'spcov')['de'] - 1
	biassig2[7,2] = biassig2[7,2] + 
		coef(expsph_rowcol, type = 'spcov')['de'] - 1
	biassig2[7,3] = biassig2[7,3] + 
		coef(gausph_rowcol, type = 'spcov')['de'] - 1
	biassig2[8,1] = biassig2[8,1] + 
		coef(sphexp_rowcol, type = 'spcov')['de'] - 1
	biassig2[8,2] = biassig2[8,2] + 
		coef(expexp_rowcol, type = 'spcov')['de'] - 1
	biassig2[8,3] = biassig2[8,3] + 
		coef(gauexp_rowcol, type = 'spcov')['de'] - 1
	biassig2[9,1] = biassig2[9,1] + 
		coef(sphgau_rowcol, type = 'spcov')['de'] - 1
	biassig2[9,2] = biassig2[9,2] + 
		coef(expgau_rowcol, type = 'spcov')['de'] - 1
	biassig2[9,3] = biassig2[9,3] + 
		coef(gaugau_rowcol, type = 'spcov')['de'] - 1

	# bias in estimation of unit-distance correlation (UDC)
			#constant model
	biasUDC[1,1] = biasUDC[1,1] + 
		cor_sph(1, coef(sphsph_const, type = 'spcov')['range']) - 0.5
	biasUDC[1,2] = biasUDC[1,2] + 
		cor_sph(1, coef(expsph_const, type = 'spcov')['range']) - 0.5
	biasUDC[1,3] = biasUDC[1,3] + 
		cor_sph(1, coef(gausph_const, type = 'spcov')['range']) - 0.5
	biasUDC[2,1] = biasUDC[2,1] + 
		cor_exp(1, coef(sphexp_const, type = 'spcov')['range']) - 0.5
	biasUDC[2,2] = biasUDC[2,2] + 
		cor_exp(1, coef(expexp_const, type = 'spcov')['range']) - 0.5
	biasUDC[2,3] = biasUDC[2,3] + 
		cor_exp(1, coef(gauexp_const, type = 'spcov')['range']) - 0.5
	biasUDC[3,1] = biasUDC[3,1] + 
		cor_gau(1, coef(sphgau_const, type = 'spcov')['range']) - 0.5
	biasUDC[3,2] = biasUDC[3,2] + 
		cor_gau(1, coef(expgau_const, type = 'spcov')['range']) - 0.5
	biasUDC[3,3] = biasUDC[3,3] + 
		cor_gau(1, coef(gaugau_const, type = 'spcov')['range']) - 0.5
			#row model
	biasUDC[4,1] = biasUDC[4,1] + 
		cor_sph(1, coef(sphsph_row, type = 'spcov')['range']) - 0.5
	biasUDC[4,2] = biasUDC[4,2] + 
		cor_sph(1, coef(expsph_row, type = 'spcov')['range']) - 0.5
	biasUDC[4,3] = biasUDC[4,3] + 
		cor_sph(1, coef(gausph_row, type = 'spcov')['range']) - 0.5
	biasUDC[5,1] = biasUDC[5,1] + 
		cor_exp(1, coef(sphexp_row, type = 'spcov')['range']) - 0.5
	biasUDC[5,2] = biasUDC[5,2] + 
		cor_exp(1, coef(expexp_row, type = 'spcov')['range']) - 0.5
	biasUDC[5,3] = biasUDC[5,3] + 
		cor_exp(1, coef(gauexp_row, type = 'spcov')['range']) - 0.5
	biasUDC[6,1] = biasUDC[6,1] + 
		cor_gau(1, coef(sphgau_row, type = 'spcov')['range']) - 0.5
	biasUDC[6,2] = biasUDC[6,2] + 
		cor_gau(1, coef(expgau_row, type = 'spcov')['range']) - 0.5
	biasUDC[6,3] = biasUDC[6,3] + 
		cor_gau(1, coef(gaugau_row, type = 'spcov')['range']) - 0.5
			#row and column model
	biasUDC[7,1] = biasUDC[7,1] + 
		cor_sph(1, coef(sphsph_rowcol, type = 'spcov')['range']) - 0.5
	biasUDC[7,2] = biasUDC[7,2] + 
		cor_sph(1, coef(expsph_rowcol, type = 'spcov')['range']) - 0.5
	biasUDC[7,3] = biasUDC[7,3] + 
		cor_sph(1, coef(gausph_rowcol, type = 'spcov')['range']) - 0.5
	biasUDC[8,1] = biasUDC[8,1] + 
		cor_exp(1, coef(sphexp_rowcol, type = 'spcov')['range']) - 0.5
	biasUDC[8,2] = biasUDC[8,2] + 
		cor_exp(1, coef(expexp_rowcol, type = 'spcov')['range']) - 0.5
	biasUDC[8,3] = biasUDC[8,3] + 
		cor_exp(1, coef(gauexp_rowcol, type = 'spcov')['range']) - 0.5
	biasUDC[9,1] = biasUDC[9,1] + 
		cor_gau(1, coef(sphgau_rowcol, type = 'spcov')['range']) - 0.5
	biasUDC[9,2] = biasUDC[9,2] + 
		cor_gau(1, coef(expgau_rowcol, type = 'spcov')['range']) - 0.5
	biasUDC[9,3] = biasUDC[9,3] + 
		cor_gau(1, coef(gaugau_rowcol, type = 'spcov')['range']) - 0.5

	# MSE in estimation of overall variance parameter
			#constant model
	MSEsig2[1,1] = MSEsig2[1,1] + 
		(coef(sphsph_const, type = 'spcov')['de'] - 1)^2
	MSEsig2[1,2] = MSEsig2[1,2] + 
		(coef(expsph_const, type = 'spcov')['de'] - 1)^2
	MSEsig2[1,3] = MSEsig2[1,3] + 
		(coef(gausph_const, type = 'spcov')['de'] - 1)^2
	MSEsig2[2,1] = MSEsig2[2,1] + 
		(coef(sphexp_const, type = 'spcov')['de'] - 1)^2
	MSEsig2[2,2] = MSEsig2[2,2] + 
		(coef(expexp_const, type = 'spcov')['de'] - 1)^2
	MSEsig2[2,3] = MSEsig2[2,3] + 
		(coef(gauexp_const, type = 'spcov')['de'] - 1)^2
	MSEsig2[3,1] = MSEsig2[3,1] + 
		(coef(sphgau_const, type = 'spcov')['de'] - 1)^2
	MSEsig2[3,2] = MSEsig2[3,2] + 
		(coef(expgau_const, type = 'spcov')['de'] - 1)^2
	MSEsig2[3,3] = MSEsig2[3,3] + 
		(coef(gaugau_const, type = 'spcov')['de'] - 1)^2
			#row model
	MSEsig2[4,1] = MSEsig2[4,1] + 
		(coef(sphsph_row, type = 'spcov')['de'] - 1)^2
	MSEsig2[4,2] = MSEsig2[4,2] + 
		(coef(expsph_row, type = 'spcov')['de'] - 1)^2
	MSEsig2[4,3] = MSEsig2[4,3] + 
		(coef(gausph_row, type = 'spcov')['de'] - 1)^2
	MSEsig2[5,1] = MSEsig2[5,1] + 
		(coef(sphexp_row, type = 'spcov')['de'] - 1)^2
	MSEsig2[5,2] = MSEsig2[5,2] + 
		(coef(expexp_row, type = 'spcov')['de'] - 1)^2
	MSEsig2[5,3] = MSEsig2[5,3] + 
		(coef(gauexp_row, type = 'spcov')['de'] - 1)^2
	MSEsig2[6,1] = MSEsig2[6,1] + 
		(coef(sphgau_row, type = 'spcov')['de'] - 1)^2
	MSEsig2[6,2] = MSEsig2[6,2] + 
		(coef(expgau_row, type = 'spcov')['de'] - 1)^2
	MSEsig2[6,3] = MSEsig2[6,3] + 
		(coef(gaugau_row, type = 'spcov')['de'] - 1)^2
			#row and column model
	MSEsig2[7,1] = MSEsig2[7,1] + 
		(coef(sphsph_rowcol, type = 'spcov')['de'] - 1)^2
	MSEsig2[7,2] = MSEsig2[7,2] + 
		(coef(expsph_rowcol, type = 'spcov')['de'] - 1)^2
	MSEsig2[7,3] = MSEsig2[7,3] + 
		(coef(gausph_rowcol, type = 'spcov')['de'] - 1)^2
	MSEsig2[8,1] = MSEsig2[8,1] + 
		(coef(sphexp_rowcol, type = 'spcov')['de'] - 1)^2
	MSEsig2[8,2] = MSEsig2[8,2] + 
		(coef(expexp_rowcol, type = 'spcov')['de'] - 1)^2
	MSEsig2[8,3] = MSEsig2[8,3] + 
		(coef(gauexp_rowcol, type = 'spcov')['de'] - 1)^2
	MSEsig2[9,1] = MSEsig2[9,1] + 
		(coef(sphgau_rowcol, type = 'spcov')['de'] - 1)^2
	MSEsig2[9,2] = MSEsig2[9,2] + 
		(coef(expgau_rowcol, type = 'spcov')['de'] - 1)^2
	MSEsig2[9,3] = MSEsig2[9,3] + 
		(coef(gaugau_rowcol, type = 'spcov')['de'] - 1)^2
		
	# MSE in estimation of unit-distance correlation (UDC)
			#constant model
	MSEUDC[1,1] = MSEUDC[1,1] + 
		(cor_sph(1, coef(sphsph_const, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[1,2] = MSEUDC[1,2] + 
		(cor_sph(1, coef(expsph_const, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[1,3] = MSEUDC[1,3] + 
		(cor_sph(1, coef(gausph_const, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[2,1] = MSEUDC[2,1] + 
		(cor_exp(1, coef(sphexp_const, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[2,2] = MSEUDC[2,2] + 
		(cor_exp(1, coef(expexp_const, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[2,3] = MSEUDC[2,3] + 
		(cor_exp(1, coef(gauexp_const, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[3,1] = MSEUDC[3,1] + 
		(cor_gau(1, coef(sphgau_const, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[3,2] = MSEUDC[3,2] + 
		(cor_gau(1, coef(expgau_const, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[3,3] = MSEUDC[3,3] + 
		(cor_gau(1, coef(gaugau_const, type = 'spcov')['range']) - 0.5)^2
			#row model
	MSEUDC[4,1] = MSEUDC[4,1] + 
		(cor_sph(1, coef(sphsph_row, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[4,2] = MSEUDC[4,2] + 
		(cor_sph(1, coef(expsph_row, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[4,3] = MSEUDC[4,3] + 
		(cor_sph(1, coef(gausph_row, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[5,1] = MSEUDC[5,1] + 
		(cor_exp(1, coef(sphexp_row, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[5,2] = MSEUDC[5,2] + 
		(cor_exp(1, coef(expexp_row, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[5,3] = MSEUDC[5,3] + 
		(cor_exp(1, coef(gauexp_row, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[6,1] = MSEUDC[6,1] + 
		(cor_gau(1, coef(sphgau_row, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[6,2] = MSEUDC[6,2] + 
		(cor_gau(1, coef(expgau_row, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[6,3] = MSEUDC[6,3] + 
		(cor_gau(1, coef(gaugau_row, type = 'spcov')['range']) - 0.5)^2
			#row and column model
	MSEUDC[7,1] = MSEUDC[7,1] + 
		(cor_sph(1, coef(sphsph_rowcol, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[7,2] = MSEUDC[7,2] + 
		(cor_sph(1, coef(expsph_rowcol, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[7,3] = MSEUDC[7,3] + 
		(cor_sph(1, coef(gausph_rowcol, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[8,1] = MSEUDC[8,1] + 
		(cor_exp(1, coef(sphexp_rowcol, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[8,2] = MSEUDC[8,2] + 
		(cor_exp(1, coef(expexp_rowcol, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[8,3] = MSEUDC[8,3] + 
		(cor_exp(1, coef(gauexp_rowcol, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[9,1] = MSEUDC[9,1] + 
		(cor_gau(1, coef(sphgau_rowcol, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[9,2] = MSEUDC[9,2] + 
		(cor_gau(1, coef(expgau_rowcol, type = 'spcov')['range']) - 0.5)^2
	MSEUDC[9,3] = MSEUDC[9,3] + 
		(cor_gau(1, coef(gaugau_rowcol, type = 'spcov')['range']) - 0.5)^2

}

Table1 = cbind(lengthCIbeta/niter,
	CI90/niter)

Table2 = cbind(biassig2/niter,
	biasUDC/niter)

Table3 = cbind(MSEsig2/niter,
	MSEUDC/niter)

print(
    xtable(Table1, 
      align = c('l',rep('l', times = length(Table1[1,]))),
      digits = c(0,rep(3, times = 6)),
      caption = 'Table 8.5'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

print(
    xtable(Table2, 
      align = c('l',rep('l', times = length(Table2[1,]))),
      digits = c(0,rep(3, times = 6)),
      caption = 'Table 8.6'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

print(
    xtable(Table3, 
      align = c('l',rep('l', times = length(Table3[1,]))),
      digits = c(0,rep(3, times = 6)),
      caption = 'Table 8.7'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
