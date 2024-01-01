sec_path = 'Rcode/Chapter9/Section 9.9/'
setwd(paste0(SLEDbook_path,sec_path))

library(spmodel)
library(xtable)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Kriging in Practice: ML vs REML Predictions
#              Table 9.6 with nugget and RMSPE
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# Code for comparing quality of predictions based on ML and REML estimation 
# (Table 9.6)
# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
# use an even number for K so that prediction locations are between observed
# locations
K <- 12
	
# create data.frames for observed and prediction locations
obsDF <- data.frame(x = matrix(rep(1:K,K),ncol=1,byrow=T),
	y = matrix(rep(1:K,each=K),ncol=1,byrow=T))
predDF = data.frame(x = c((K + 1)/2, (K + 1)/2), y = c((K + 1)/2, 0.5))
xyall = rbind(obsDF[,1:2], predDF[,1:2])
n <- K^2

# theta^d, where d is distance, is equal to an exponential model with more
# typical form: e^(-d/alpha), where alpha = -1/log(theta).  So, 
# theta = 0.3 implies alpha = 0.8305835
spcov_params_val <- spcov_params("exponential", de = 0.9, ie = 0.1, 
	range = -1/log(0.3))

set.seed(2222)

#number of simulations
nsim = 5000
# store ml performance at theta = 0.3
ml03_rmspe = matrix(NA, nrow = nsim, ncol = 2)
ml03_width = matrix(NA, nrow = nsim, ncol = 2)
ml03_cover = matrix(NA, nrow = nsim, ncol = 2)
reml03_rmspe = matrix(NA, nrow = nsim, ncol = 2)
reml03_width = matrix(NA, nrow = nsim, ncol = 2)
reml03_cover = matrix(NA, nrow = nsim, ncol = 2)
true03_rmspe = matrix(NA, nrow = nsim, ncol = 2)
true03_width = matrix(NA, nrow = nsim, ncol = 2)
true03_cover = matrix(NA, nrow = nsim, ncol = 2)
for(kk in 1:nsim) {
	cat("\r", "iteration: ", kk)
	# simulate errors
	errsim = sprnorm(spcov_params_val, data = xyall,
		xcoord = x, ycoord = y)
	# constant but unknown mean = 0, simulated values at observed and prediction
	# locations
	z = errsim[1:n]
	zp = errsim[(n + 1):(n + 2)]
	obsDF$z = z

	# ml prediction at both locations
	ml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_type = "exponential", 
		estmethod = 'ml')
	ml_pred = predict(ml_out, predDF, se.fit = TRUE)
	# reml prediction at both locations
	reml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_type = "exponential", 
		estmethod = 'reml')
	reml_pred = predict(reml_out, predDF, se.fit = TRUE)
	#true
	true = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial("exponential", ie = 0.1, de = 0.9, 
			range = -1/log(0.3), known = c("given")))
	true_pred = predict(true, predDF, se.fit = TRUE)
	# rmspe
	ml03_rmspe[kk,] = (ml_pred$fit - zp)^2
	# width of 90% interval using t-distribution
	p = 1
	ml03_width[kk,] = 2*ml_pred$se.fit*qt(.90, n - p)
	# coverage
	ml03_cover[kk,] = ml_pred$fit - ml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < ml_pred$fit + ml_pred$se.fit*qt(.90, n - p)
	# rmspe
	reml03_rmspe[kk,] = (reml_pred$fit - zp)^2
	# width of 90% interval using t-distribution
	reml03_width[kk,] = 2*reml_pred$se.fit*qt(.90, n - p)
	# coverage
	reml03_cover[kk,] = reml_pred$fit - reml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < reml_pred$fit + reml_pred$se.fit*qt(.90,n - p)
	# rmspe
	true03_rmspe[kk,] = (true_pred$fit - zp)^2
	# width of 90% interval using t-distribution
	true03_width[kk,] = 2*true_pred$se.fit*qt(.90, n - p)
	# coverage
	true03_cover[kk,] = true_pred$fit - true_pred$se.fit*qt(.90,n - p) < zp & 
		zp < true_pred$fit + true_pred$se.fit*qt(.90,n - p)
}
sqrt(apply(ml03_rmspe,2,mean))
sqrt(apply(reml03_rmspe,2,mean))
sqrt(apply(true03_rmspe,2,mean))
apply(ml03_width,2,mean)
apply(reml03_width,2,mean)
apply(true03_width,2,mean)
apply(ml03_cover,2,mean)
apply(reml03_cover,2,mean)
apply(true03_cover,2,mean)



#theta^d, where d is distance, is equal to an exponential model with more
# typical form e^(-d/alpha), where alpha = -1/log(theta).  So, 
# theta = 0.7 implies alpha = 2.803673
spcov_params_val <- spcov_params("exponential", de = 0.9, ie = 0.1, 
	range = -1/log(0.7))

#number of simulations
nsim = 5000
# store ml performance at theta = 0.3
ml07_rmspe = matrix(NA, nrow = nsim, ncol = 2)
ml07_width = matrix(NA, nrow = nsim, ncol = 2)
ml07_cover = matrix(NA, nrow = nsim, ncol = 2)
reml07_rmspe = matrix(NA, nrow = nsim, ncol = 2)
reml07_width = matrix(NA, nrow = nsim, ncol = 2)
reml07_cover = matrix(NA, nrow = nsim, ncol = 2)
true07_rmspe = matrix(NA, nrow = nsim, ncol = 2)
true07_width = matrix(NA, nrow = nsim, ncol = 2)
true07_cover = matrix(NA, nrow = nsim, ncol = 2)
for(kk in 1:nsim) {
	cat("\r", "iteration: ", kk)
	# simulate errors
	errsim = sprnorm(spcov_params_val, data = xyall,
		xcoord = x, ycoord = y)
	# constant but unknown mean = 0, simulated values at observed and prediction
	# locations
	z = errsim[1:n]
	zp = errsim[(n + 1):(n + 2)]
	obsDF$z = z

	# ml prediction at both locations
	ml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_type = "exponential", 
		estmethod = 'ml')
	ml_pred = predict(ml_out, predDF, se.fit = TRUE)
	# reml prediction at both locations
	reml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_type = "exponential", 
		estmethod = 'reml')
	reml_pred = predict(reml_out, predDF, se.fit = TRUE)
	#true
	true = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial("exponential", ie = 0.1, de = 0.9, 
			range = -1/log(0.7), known = c("given")))
	true_pred = predict(true, predDF, se.fit = TRUE)
	# rmspe
	ml07_rmspe[kk,] = (ml_pred$fit - zp)^2
	# width of 90% interval using t-distribution
	ml07_width[kk,] = 2*ml_pred$se.fit*qt(.90, n - p)
	# coverage
	ml07_cover[kk,] = ml_pred$fit - ml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < ml_pred$fit + ml_pred$se.fit*qt(.90,n - p)
	reml07_rmspe[kk,] = (reml_pred$fit - zp)^2
	# width of 90% interval using t-distribution
	reml07_width[kk,] = 2*reml_pred$se.fit*qt(.90, n - p)
	# coverage
	reml07_cover[kk,] = reml_pred$fit - reml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < reml_pred$fit + reml_pred$se.fit*qt(.90,n - p)
	# rmspe
	true07_rmspe[kk,] = (true_pred$fit - zp)^2
	# width of 90% interval using t-distribution
	true07_width[kk,] = 2*true_pred$se.fit*qt(.90, n - p)
	# coverage
	true07_cover[kk,] = true_pred$fit - true_pred$se.fit*qt(.90,n - p) < zp & 
		zp < true_pred$fit + true_pred$se.fit*qt(.90,n - p)
}
sqrt(apply(ml07_rmspe,2,mean))
sqrt(apply(reml07_rmspe,2,mean))
sqrt(apply(true07_rmspe,2,mean))
apply(ml07_width,2,mean)
apply(reml07_width,2,mean)
apply(true07_width,2,mean)
apply(ml07_cover,2,mean)
apply(reml07_cover,2,mean)
apply(true07_cover,2,mean)


# create the table
ml_vs_repl_pred = rbind(
	c(sqrt(apply(ml03_rmspe,2,mean)), sqrt(apply(ml07_rmspe,2,mean))),
	c(sqrt(apply(reml03_rmspe,2,mean)), sqrt(apply(reml07_rmspe,2,mean))),
	c(sqrt(apply(true03_rmspe,2,mean)), sqrt(apply(true07_rmspe,2,mean))),
	c(apply(ml03_width,2,mean), apply(ml07_width,2,mean)),
  c(apply(reml03_width,2,mean), apply(reml07_width,2,mean)),
  c(apply(true03_width,2,mean), apply(true07_width,2,mean)),
	c(apply(ml03_cover,2,mean), apply(ml07_cover,2,mean)),
	c(apply(reml03_cover,2,mean), apply(reml07_cover,2,mean)),
	c(apply(true03_cover,2,mean), apply(true07_cover,2,mean))
)

print(
    xtable(ml_vs_repl_pred, 
      align = c('l',rep('l', times = length(ml_vs_repl_pred[1,]))),
      digits = c(0, rep(3, times = 4))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#             Robustness of kriging to autocorrelation models
#                          Table 9.7
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
#                       Simulations
#-------------------------------------------------------------------------------

set.seed(1507)
bias = matrix(0, nrow = 12, ncol = 3)
RMSPE = matrix(0, nrow = 12, ncol = 3)
PIlength = matrix(0, nrow = 12, ncol = 3)
coverage = matrix(0, nrow = 12, ncol = 3)
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
	indNA = sort(sample(1:n, 10))
	simsph_p = DF[indNA,'simsph']
	simexp_p = DF[indNA,'simexp']
	simgau_p = DF[indNA,'simgau']
	DF[indNA,'simsph'] = NA
	DF[indNA,'simexp'] = NA
	DF[indNA,'simgau'] = NA
	# ------- ML ESTIMATION
	# fit model using spmodel, holding nugget effect at 0.000001
	# in the name, true model if first, fitted model is second
	# sph = spherical, exp = exponential, gau = gaussian, unc = uncorrelated
	
	               # constant mean
	sphsph_const = splm(simsph ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphsph_const_pred = predict(sphsph_const, interval = 'prediction', 
		level = 0.90)
	sphexp_const = splm(simsph ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphexp_const_pred = predict(sphexp_const, interval = 'prediction', 
		level = 0.90)
	sphgau_const = splm(simsph ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphgau_const_pred = predict(sphgau_const, interval = 'prediction', 
		level = 0.90)
	sphunc_const = splm(simsph ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_type = 'none', estmethod = 'reml')
	sphunc_const_pred = predict(sphunc_const, interval = 'prediction', 
		level = 0.90)
		
	expsph_const = splm(simexp ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expsph_const_pred = predict(expsph_const, interval = 'prediction', 
		level = 0.90)
	expexp_const = splm(simexp ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expexp_const_pred = predict(expexp_const, interval = 'prediction', 
		level = 0.90)
	expgau_const = splm(simexp ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expgau_const_pred = predict(expgau_const, interval = 'prediction', 
		level = 0.90)
	expunc_const = splm(simexp ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_type = 'none', estmethod = 'reml')
	expunc_const_pred = predict(expunc_const, interval = 'prediction', 
		level = 0.90)

	gausph_const = splm(simgau ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gausph_const_pred = predict(gausph_const, interval = 'prediction', 
		level = 0.90)
	gauexp_const = splm(simgau ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gauexp_const_pred = predict(gauexp_const, interval = 'prediction', 
		level = 0.90)
	gaugau_const = splm(simgau ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gaugau_const_pred = predict(gaugau_const, interval = 'prediction', 
		level = 0.90)
	gauunc_const = splm(simgau ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_type = 'none', estmethod = 'reml')
	gauunc_const_pred = predict(gauunc_const, interval = 'prediction', 
		level = 0.90)


	               # row effects
	sphsph_row = splm(simsph ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphsph_row_pred = predict(sphsph_row, interval = 'prediction', 
		level = 0.90)
	sphexp_row = splm(simsph ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphexp_row_pred = predict(sphexp_row, interval = 'prediction', 
		level = 0.90)
	sphgau_row = splm(simsph ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphgau_row_pred = predict(sphgau_row, interval = 'prediction', 
		level = 0.90)
	sphunc_row = splm(simsph ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_type = 'none', estmethod = 'reml')
	sphunc_row_pred = predict(sphunc_row, interval = 'prediction', 
		level = 0.90)
		
	expsph_row = splm(simexp ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expsph_row_pred = predict(expsph_row, interval = 'prediction', 
		level = 0.90)
	expexp_row = splm(simexp ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expexp_row_pred = predict(expexp_row, interval = 'prediction', 
		level = 0.90)
	expgau_row = splm(simexp ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expgau_row_pred = predict(expgau_row, interval = 'prediction', 
		level = 0.90)
	expunc_row = splm(simexp ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_type = 'none', estmethod = 'reml')
	expunc_row_pred = predict(expunc_row, interval = 'prediction', 
		level = 0.90)

	gausph_row = splm(simgau ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gausph_row_pred = predict(gausph_row, interval = 'prediction', 
		level = 0.90)
	gauexp_row = splm(simgau ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gauexp_row_pred = predict(gauexp_row, interval = 'prediction', 
		level = 0.90)
	gaugau_row = splm(simgau ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gaugau_row_pred = predict(gaugau_row, interval = 'prediction', 
		level = 0.90)
	gauunc_row = splm(simgau ~ I(as.factor(xcoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_type = 'none', estmethod = 'reml')
	gauunc_row_pred = predict(gauunc_row, interval = 'prediction', 
		level = 0.90)
	
	               # row and column effects
	               
	               # constant mean
	sphsph_rowcol = splm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphsph_rowcol_pred = predict(sphsph_rowcol, interval = 'prediction', 
		level = 0.90)
	sphexp_rowcol = splm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphexp_rowcol_pred = predict(sphexp_rowcol, interval = 'prediction', 
		level = 0.90)
	sphgau_rowcol = splm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	sphgau_rowcol_pred = predict(sphgau_rowcol, interval = 'prediction', 
		level = 0.90)
	sphunc_rowcol = splm(simsph ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_type = 'none', estmethod = 'reml')
	sphunc_rowcol_pred = predict(sphunc_rowcol, interval = 'prediction', 
		level = 0.90)
		
	expsph_rowcol = splm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expsph_rowcol_pred = predict(expsph_rowcol, interval = 'prediction', 
		level = 0.90)
	expexp_rowcol = splm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expexp_rowcol_pred = predict(expexp_rowcol, interval = 'prediction', 
		level = 0.90)
	expgau_rowcol = splm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	expgau_rowcol_pred = predict(expgau_rowcol, interval = 'prediction', 
		level = 0.90)
	expunc_rowcol = splm(simexp ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_type = 'none', estmethod = 'reml')
	expunc_rowcol_pred = predict(expunc_rowcol, interval = 'prediction', 
		level = 0.90)

	gausph_rowcol = splm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('spherical', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gausph_rowcol_pred = predict(gausph_rowcol, interval = 'prediction', 
		level = 0.90)
	gauexp_rowcol = splm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gauexp_rowcol_pred = predict(gauexp_rowcol, interval = 'prediction', 
		level = 0.90)
	gaugau_rowcol = splm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('gaussian', ie = 0.000001, known = c('ie')),
		estmethod = 'reml')
	gaugau_rowcol_pred = predict(gaugau_rowcol, interval = 'prediction', 
		level = 0.90)
	gauunc_rowcol = splm(simgau ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_type = 'none', estmethod = 'reml')
	gauunc_rowcol_pred = predict(gauunc_rowcol, interval = 'prediction', 
		level = 0.90)

	#
	#                           prediction summaries
	#
	# bias
			#constant model
	bias[1,1] = bias[1,1] + mean(sphsph_const_pred[,'fit'] - simsph_p)
	bias[2,1] = bias[2,1] + mean(sphexp_const_pred[,'fit'] - simsph_p)
	bias[3,1] = bias[3,1] + mean(sphgau_const_pred[,'fit'] - simsph_p)
	bias[4,1] = bias[4,1] + mean(sphunc_const_pred[,'fit'] - simsph_p)
	bias[1,2] = bias[1,2] + mean(expsph_const_pred[,'fit'] - simexp_p)
	bias[2,2] = bias[2,2] + mean(expexp_const_pred[,'fit'] - simexp_p)
	bias[3,2] = bias[3,2] + mean(expgau_const_pred[,'fit'] - simexp_p)
	bias[4,2] = bias[4,2] + mean(expunc_const_pred[,'fit'] - simexp_p)
	bias[1,3] = bias[1,3] + mean(gausph_const_pred[,'fit'] - simgau_p)
	bias[2,3] = bias[2,3] + mean(gauexp_const_pred[,'fit'] - simgau_p)
	bias[3,3] = bias[3,3] + mean(gaugau_const_pred[,'fit'] - simgau_p)
	bias[4,3] = bias[4,3] + mean(gauunc_const_pred[,'fit'] - simgau_p)
			#row model
	bias[5,1] = bias[5,1] + mean(sphsph_row_pred[,'fit'] - simsph_p)
	bias[6,1] = bias[6,1] + mean(sphexp_row_pred[,'fit'] - simsph_p)
	bias[7,1] = bias[7,1] + mean(sphgau_row_pred[,'fit'] - simsph_p)
	bias[8,1] = bias[8,1] + mean(sphunc_row_pred[,'fit'] - simsph_p)
	bias[5,2] = bias[5,2] + mean(expsph_row_pred[,'fit'] - simexp_p)
	bias[6,2] = bias[6,2] + mean(expexp_row_pred[,'fit'] - simexp_p)
	bias[7,2] = bias[7,2] + mean(expgau_row_pred[,'fit'] - simexp_p)
	bias[8,2] = bias[8,2] + mean(expunc_row_pred[,'fit'] - simexp_p)
	bias[5,3] = bias[5,3] + mean(gausph_row_pred[,'fit'] - simgau_p)
	bias[6,3] = bias[6,3] + mean(gauexp_row_pred[,'fit'] - simgau_p)
	bias[7,3] = bias[7,3] + mean(gaugau_row_pred[,'fit'] - simgau_p)
	bias[8,3] = bias[8,3] + mean(gauunc_row_pred[,'fit'] - simgau_p)
			#row + column model
	bias[9,1] = bias[9,1] + mean(sphsph_rowcol_pred[,'fit'] - simsph_p)
	bias[10,1] = bias[10,1] + mean(sphexp_rowcol_pred[,'fit'] - simsph_p)
	bias[11,1] = bias[11,1] + mean(sphgau_rowcol_pred[,'fit'] - simsph_p)
	bias[12,1] = bias[12,1] + mean(sphunc_rowcol_pred[,'fit'] - simsph_p)
	bias[9,2] = bias[9,2] + mean(expsph_rowcol_pred[,'fit'] - simexp_p)
	bias[10,2] = bias[10,2] + mean(expexp_rowcol_pred[,'fit'] - simexp_p)
	bias[11,2] = bias[11,2] + mean(expgau_rowcol_pred[,'fit'] - simexp_p)
	bias[12,2] = bias[12,2] + mean(expunc_rowcol_pred[,'fit'] - simexp_p)
	bias[9,3] = bias[9,3] + mean(gausph_rowcol_pred[,'fit'] - simgau_p)
	bias[10,3] = bias[10,3] + mean(gauexp_rowcol_pred[,'fit'] - simgau_p)
	bias[11,3] = bias[11,3] + mean(gaugau_rowcol_pred[,'fit'] - simgau_p)
	bias[12,3] = bias[12,3] + mean(gauunc_rowcol_pred[,'fit'] - simgau_p)
			
			
	# RMSPE
			#constant model
	RMSPE[1,1] = RMSPE[1,1] + mean((sphsph_const_pred[,'fit'] - simsph_p)^2)/10
	RMSPE[2,1] = RMSPE[2,1] + mean((sphexp_const_pred[,'fit'] - simsph_p)^2)/10
	RMSPE[3,1] = RMSPE[3,1] + mean((sphgau_const_pred[,'fit'] - simsph_p)^2)/10
	RMSPE[4,1] = RMSPE[4,1] + mean((sphunc_const_pred[,'fit'] - simsph_p)^2)/10
	RMSPE[1,2] = RMSPE[1,2] + mean((expsph_const_pred[,'fit'] - simexp_p)^2)/10
	RMSPE[2,2] = RMSPE[2,2] + mean((expexp_const_pred[,'fit'] - simexp_p)^2)/10
	RMSPE[3,2] = RMSPE[3,2] + mean((expgau_const_pred[,'fit'] - simexp_p)^2)/10
	RMSPE[4,2] = RMSPE[4,2] + mean((expunc_const_pred[,'fit'] - simexp_p)^2)/10
	RMSPE[1,3] = RMSPE[1,3] + mean((gausph_const_pred[,'fit'] - simgau_p)^2)/10
	RMSPE[2,3] = RMSPE[2,3] + mean((gauexp_const_pred[,'fit'] - simgau_p)^2)/10
	RMSPE[3,3] = RMSPE[3,3] + mean((gaugau_const_pred[,'fit'] - simgau_p)^2)/10
	RMSPE[4,3] = RMSPE[4,3] + mean((gauunc_const_pred[,'fit'] - simgau_p)^2)/10
			#row model
	RMSPE[5,1] = RMSPE[5,1] + mean((sphsph_row_pred[,'fit'] - simsph_p)^2)/10
	RMSPE[6,1] = RMSPE[6,1] + mean((sphexp_row_pred[,'fit'] - simsph_p)^2)/10
	RMSPE[7,1] = RMSPE[7,1] + mean((sphgau_row_pred[,'fit'] - simsph_p)^2)/10
	RMSPE[8,1] = RMSPE[8,1] + mean((sphunc_row_pred[,'fit'] - simsph_p)^2)/10
	RMSPE[5,2] = RMSPE[5,2] + mean((expsph_row_pred[,'fit'] - simexp_p)^2)/10
	RMSPE[6,2] = RMSPE[6,2] + mean((expexp_row_pred[,'fit'] - simexp_p)^2)/10
	RMSPE[7,2] = RMSPE[7,2] + mean((expgau_row_pred[,'fit'] - simexp_p)^2)/10
	RMSPE[8,2] = RMSPE[8,2] + mean((expunc_row_pred[,'fit'] - simexp_p)^2)/10
	RMSPE[5,3] = RMSPE[5,3] + mean((gausph_row_pred[,'fit'] - simgau_p)^2)/10
	RMSPE[6,3] = RMSPE[6,3] + mean((gauexp_row_pred[,'fit'] - simgau_p)^2)/10
	RMSPE[7,3] = RMSPE[7,3] + mean((gaugau_row_pred[,'fit'] - simgau_p)^2)/10
	RMSPE[8,3] = RMSPE[8,3] + mean((gauunc_row_pred[,'fit'] - simgau_p)^2)/10
			#row + column model
	RMSPE[9,1] = RMSPE[9,1] + mean((sphsph_rowcol_pred[,'fit'] - simsph_p)^2)
	RMSPE[10,1] = RMSPE[10,1] + mean((sphexp_rowcol_pred[,'fit'] - simsph_p)^2)
	RMSPE[11,1] = RMSPE[11,1] + mean((sphgau_rowcol_pred[,'fit'] - simsph_p)^2)
	RMSPE[12,1] = RMSPE[12,1] + mean((sphunc_rowcol_pred[,'fit'] - simsph_p)^2)
	RMSPE[9,2] = RMSPE[9,2] + mean((expsph_rowcol_pred[,'fit'] - simexp_p)^2)
	RMSPE[10,2] = RMSPE[10,2] + mean((expexp_rowcol_pred[,'fit'] - simexp_p)^2)
	RMSPE[11,2] = RMSPE[11,2] + mean((expgau_rowcol_pred[,'fit'] - simexp_p)^2)
	RMSPE[12,2] = RMSPE[12,2] + mean((expunc_rowcol_pred[,'fit'] - simexp_p)^2)
	RMSPE[9,3] = RMSPE[9,3] + mean((gausph_rowcol_pred[,'fit'] - simgau_p)^2)
	RMSPE[10,3] = RMSPE[10,3] + mean((gauexp_rowcol_pred[,'fit'] - simgau_p)^2)
	RMSPE[11,3] = RMSPE[11,3] + mean((gaugau_rowcol_pred[,'fit'] - simgau_p)^2)
	RMSPE[12,3] = RMSPE[12,3] + mean((gauunc_rowcol_pred[,'fit'] - simgau_p)^2)
			

	# Prediction interval length
			#constant model
	PIlength[1,1] = PIlength[1,1] + mean(sphsph_const_pred[,'upr'] -
		sphsph_const_pred[,'lwr'])
	PIlength[2,1] = PIlength[2,1] + mean(sphexp_const_pred[,'upr'] -
		sphexp_const_pred[,'lwr'])
	PIlength[3,1] = PIlength[3,1] + mean(sphgau_const_pred[,'upr'] -
		sphgau_const_pred[,'lwr'])
	PIlength[4,1] = PIlength[4,1] + mean(sphunc_const_pred[,'upr'] -
		sphunc_const_pred[,'lwr'])
	PIlength[1,2] = PIlength[1,2] + mean(expsph_const_pred[,'upr'] -
		expsph_const_pred[,'lwr'])
	PIlength[2,2] = PIlength[2,2] + mean(expexp_const_pred[,'upr'] -
		expexp_const_pred[,'lwr'])
	PIlength[3,2] = PIlength[3,2] + mean(expgau_const_pred[,'upr'] -
		expgau_const_pred[,'lwr'])
	PIlength[4,2] = PIlength[4,2] + mean(expunc_const_pred[,'upr'] -
		expunc_const_pred[,'lwr'])
	PIlength[1,3] = PIlength[1,3] + mean(gausph_const_pred[,'upr'] -
		gausph_const_pred[,'lwr'])
	PIlength[2,3] = PIlength[2,3] + mean(gauexp_const_pred[,'upr'] -
		gauexp_const_pred[,'lwr'])
	PIlength[3,3] = PIlength[3,3] + mean(gaugau_const_pred[,'upr'] -
		gaugau_const_pred[,'lwr'])
	PIlength[4,3] = PIlength[4,3] + mean(gauunc_const_pred[,'upr'] -
		gauunc_const_pred[,'lwr'])
			#row model
	PIlength[5,1] = PIlength[5,1] + mean(sphsph_row_pred[,'upr'] -
		sphsph_row_pred[,'lwr'])
	PIlength[6,1] = PIlength[6,1] + mean(sphexp_row_pred[,'upr'] -
		sphexp_row_pred[,'lwr'])
	PIlength[7,1] = PIlength[7,1] + mean(sphgau_row_pred[,'upr'] -
		sphgau_row_pred[,'lwr'])
	PIlength[8,1] = PIlength[8,1] + mean(sphunc_row_pred[,'upr'] -
		sphunc_row_pred[,'lwr'])
	PIlength[5,2] = PIlength[5,2] + mean(expsph_row_pred[,'upr'] -
		expsph_row_pred[,'lwr'])
	PIlength[6,2] = PIlength[6,2] + mean(expexp_row_pred[,'upr'] -
		expexp_row_pred[,'lwr'])
	PIlength[7,2] = PIlength[7,2] + mean(expgau_row_pred[,'upr'] -
		expgau_row_pred[,'lwr'])
	PIlength[8,2] = PIlength[8,2] + mean(expunc_row_pred[,'upr'] -
		expunc_row_pred[,'lwr'])
	PIlength[5,3] = PIlength[5,3] + mean(gausph_row_pred[,'upr'] -
		gausph_row_pred[,'lwr'])
	PIlength[6,3] = PIlength[6,3] + mean(gauexp_row_pred[,'upr'] -
		gauexp_row_pred[,'lwr'])
	PIlength[7,3] = PIlength[7,3] + mean(gaugau_row_pred[,'upr'] -
		gaugau_row_pred[,'lwr'])
	PIlength[8,3] = PIlength[8,3] + mean(gauunc_row_pred[,'upr'] -
		gauunc_row_pred[,'lwr'])
			#row + column model
	PIlength[9,1] = PIlength[9,1] + mean(sphsph_rowcol_pred[,'upr'] -
		sphsph_rowcol_pred[,'lwr'])
	PIlength[10,1] = PIlength[10,1] + mean(sphexp_rowcol_pred[,'upr'] -
		sphexp_rowcol_pred[,'lwr'])
	PIlength[11,1] = PIlength[11,1] + mean(sphgau_rowcol_pred[,'upr'] -
		sphgau_rowcol_pred[,'lwr'])
	PIlength[12,1] = PIlength[12,1] + mean(sphunc_rowcol_pred[,'upr'] -
		sphunc_rowcol_pred[,'lwr'])
	PIlength[9,2] = PIlength[9,2] + mean(expsph_rowcol_pred[,'upr'] -
		expsph_rowcol_pred[,'lwr'])
	PIlength[10,2] = PIlength[10,2] + mean(expexp_rowcol_pred[,'upr'] -
		expexp_rowcol_pred[,'lwr'])
	PIlength[11,2] = PIlength[11,2] + mean(expgau_rowcol_pred[,'upr'] -
		expgau_rowcol_pred[,'lwr'])
	PIlength[12,2] = PIlength[12,2] + mean(expunc_rowcol_pred[,'upr'] -
		expunc_rowcol_pred[,'lwr'])
	PIlength[9,3] = PIlength[9,3] + mean(gausph_rowcol_pred[,'upr'] -
		gausph_rowcol_pred[,'lwr'])
	PIlength[10,3] = PIlength[10,3] + mean(gauexp_rowcol_pred[,'upr'] -
		gauexp_rowcol_pred[,'lwr'])
	PIlength[11,3] = PIlength[11,3] + mean(gaugau_rowcol_pred[,'upr'] -
		gaugau_rowcol_pred[,'lwr'])
	PIlength[12,3] = PIlength[12,3] + mean(gauunc_rowcol_pred[,'upr'] -
		gauunc_rowcol_pred[,'lwr'])
			
	# coverage
			#constant model
	coverage[1,1] = coverage[1,1] + mean(sphsph_const_pred[,'lwr']	< simsph_p &
		simsph_p < sphsph_const_pred[,'upr'])
	coverage[2,1] = coverage[2,1] + mean(sphexp_const_pred[,'lwr']	< simsph_p &
		simsph_p < sphexp_const_pred[,'upr'])
	coverage[3,1] = coverage[3,1] + mean(sphgau_const_pred[,'lwr']	< simsph_p &
		simsph_p < sphgau_const_pred[,'upr'])
	coverage[4,1] = coverage[4,1] + mean(sphunc_const_pred[,'lwr']	< simsph_p &
		simsph_p < sphunc_const_pred[,'upr'])
	coverage[1,2] = coverage[1,2] + mean(expsph_const_pred[,'lwr']	< simexp_p &
		simexp_p < expsph_const_pred[,'upr'])
	coverage[2,2] = coverage[2,2] + mean(expexp_const_pred[,'lwr']	< simexp_p &
		simexp_p < expexp_const_pred[,'upr'])
	coverage[3,2] = coverage[3,2] + mean(expgau_const_pred[,'lwr']	< simexp_p &
		simexp_p < expgau_const_pred[,'upr'])
	coverage[4,2] = coverage[4,2] + mean(expunc_const_pred[,'lwr']	< simexp_p &
		simexp_p < expunc_const_pred[,'upr'])
	coverage[1,3] = coverage[1,3] + mean(gausph_const_pred[,'lwr']	< simgau_p &
		simgau_p < gausph_const_pred[,'upr'])
	coverage[2,3] = coverage[2,3] + mean(gauexp_const_pred[,'lwr']	< simgau_p &
		simgau_p < gauexp_const_pred[,'upr'])
	coverage[3,3] = coverage[3,3] + mean(gaugau_const_pred[,'lwr']	< simgau_p &
		simgau_p < gaugau_const_pred[,'upr'])
	coverage[4,3] = coverage[4,3] + mean(gauunc_const_pred[,'lwr']	< simgau_p &
		simgau_p < gauunc_const_pred[,'upr'])
			#row model
	coverage[5,1] = coverage[5,1] + mean(sphsph_row_pred[,'lwr']	< simsph_p &
		simsph_p < sphsph_row_pred[,'upr'])
	coverage[6,1] = coverage[6,1] + mean(sphexp_row_pred[,'lwr']	< simsph_p &
		simsph_p < sphexp_row_pred[,'upr'])
	coverage[7,1] = coverage[7,1] + mean(sphgau_row_pred[,'lwr']	< simsph_p &
		simsph_p < sphgau_row_pred[,'upr'])
	coverage[8,1] = coverage[8,1] + mean(sphunc_row_pred[,'lwr']	< simsph_p &
		simsph_p < sphunc_row_pred[,'upr'])
	coverage[5,2] = coverage[5,2] + mean(expsph_row_pred[,'lwr']	< simexp_p &
		simexp_p < expsph_row_pred[,'upr'])
	coverage[6,2] = coverage[6,2] + mean(expexp_row_pred[,'lwr']	< simexp_p &
		simexp_p < expexp_row_pred[,'upr'])
	coverage[7,2] = coverage[7,2] + mean(expgau_row_pred[,'lwr']	< simexp_p &
		simexp_p < expgau_row_pred[,'upr'])
	coverage[8,2] = coverage[8,2] + mean(expunc_row_pred[,'lwr']	< simexp_p &
		simexp_p < expunc_row_pred[,'upr'])
	coverage[5,3] = coverage[5,3] + mean(gausph_row_pred[,'lwr']	< simgau_p &
		simgau_p < gausph_row_pred[,'upr'])
	coverage[6,3] = coverage[6,3] + mean(gauexp_row_pred[,'lwr']	< simgau_p &
		simgau_p < gauexp_row_pred[,'upr'])
	coverage[7,3] = coverage[7,3] + mean(gaugau_row_pred[,'lwr']	< simgau_p &
		simgau_p < gaugau_row_pred[,'upr'])
	coverage[8,3] = coverage[8,3] + mean(gauunc_row_pred[,'lwr']	< simgau_p &
		simgau_p < gauunc_row_pred[,'upr'])
			#row + column model
	coverage[9,1] = coverage[9,1] + mean(sphsph_rowcol_pred[,'lwr']	< simsph_p &
		simsph_p < sphsph_rowcol_pred[,'upr'])
	coverage[10,1] = coverage[10,1] + mean(sphexp_rowcol_pred[,'lwr']	< simsph_p &
		simsph_p < sphexp_rowcol_pred[,'upr'])
	coverage[11,1] = coverage[11,1] + mean(sphgau_rowcol_pred[,'lwr']	< simsph_p &
		simsph_p < sphgau_rowcol_pred[,'upr'])
	coverage[12,1] = coverage[12,1] + mean(sphunc_rowcol_pred[,'lwr']	< simsph_p &
		simsph_p < sphunc_rowcol_pred[,'upr'])
	coverage[9,2] = coverage[9,2] + mean(expsph_rowcol_pred[,'lwr']	< simexp_p &
		simexp_p < expsph_rowcol_pred[,'upr'])
	coverage[10,2] = coverage[10,2] + mean(expexp_rowcol_pred[,'lwr']	< simexp_p &
		simexp_p < expexp_rowcol_pred[,'upr'])
	coverage[11,2] = coverage[11,2] + mean(expgau_rowcol_pred[,'lwr']	< simexp_p &
		simexp_p < expgau_rowcol_pred[,'upr'])
	coverage[12,2] = coverage[12,2] + mean(expunc_rowcol_pred[,'lwr']	< simexp_p &
		simexp_p < expunc_rowcol_pred[,'upr'])
	coverage[9,3] = coverage[9,3] + mean(gausph_rowcol_pred[,'lwr']	< simgau_p &
		simgau_p < gausph_rowcol_pred[,'upr'])
	coverage[10,3] = coverage[10,3] + mean(gauexp_rowcol_pred[,'lwr']	< simgau_p &
		simgau_p < gauexp_rowcol_pred[,'upr'])
	coverage[11,3] = coverage[11,3] + mean(gaugau_rowcol_pred[,'lwr']	< simgau_p &
		simgau_p < gaugau_rowcol_pred[,'upr'])
	coverage[12,3] = coverage[12,3] + mean(gauunc_rowcol_pred[,'lwr']	< simgau_p &
		simgau_p < gauunc_rowcol_pred[,'upr'])

}

tab_bias = bias/niter

tab_RMSPE = sqrt(RMSPE/niter)

tab_PIlength = PIlength/niter

tab_coverage = coverage/niter

print(
    xtable(tab_bias, 
      align = c('l',rep('l', times = length(tab_bias[1,]))),
      digits = c(0,rep(3, times = 3)),
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
    xtable(tab_RMSPE, 
      align = c('l',rep('l', times = length(tab_RMSPE[1,]))),
      digits = c(0,rep(3, times = 3)),
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
    xtable(tab_PIlength, 
      align = c('l',rep('l', times = length(tab_PIlength[1,]))),
      digits = c(0,rep(3, times = 3)),
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
    xtable(tab_coverage, 
      align = c('l',rep('l', times = length(tab_coverage[1,]))),
      digits = c(0,rep(3, times = 3)),
      caption = 'Table 8.5'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# Table 9.7
tab_9.7 = cbind(tab_PIlength, tab_coverage)
print(
    xtable(tab_9.7, 
      align = c('l',rep('l', times = length(tab_9.7[1,]))),
      digits = c(0,rep(3, times = 6)),
      caption = 'Table 9.7'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
