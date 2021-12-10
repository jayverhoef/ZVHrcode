sec_path = 'Rcode/Chapter8/Section 8.6/'
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
library(colorspace)
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


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Create some functions that we will use
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# some useful transformations
logit = function(x) {log(x/(1 - x))}
expit = function(x) {exp(x)/(1 + exp(x))}

# some spatial autocorrelation models with unscaled range
exponential_spatial_model = function(x)
{
	exp(-x)
}
circular_spatial_model <- function(x)
{
	d <- x
	d[x > 1] <- 0
	CorMat <- 2*(acos(d) - d*sqrt(1 - d^2))/pi
	CorMat[x >= 1] <- 0
	CorMat
}
spherical_spatial_model <- function(x)
{
	CorMat <- (1 - 1.5*x + 0.5*x^3)
	CorMat[x > 1] <- 0
	CorMat
}
rationalQuad_spatial_model <- function(x)
{
	1/(1+x^2)
}
Gaussian_spatial_model <- function(x)
{
	exp(-x^2) 
}
Matern_spatial_model <- function(x, extrap)
{
	d <- x
	d[x == 0] <- 1
	CorMat <- d^extrap*besselK(d, extrap)/(2^(extrap - 1)*gamma(extrap))
	CorMat[x == 0] <- 1
	CorMat
}

#minus two times the profiled loglikelihood
m2LL = function(theta, y, X, distmat, spatial_model, MLmeth)
{
	# range parameter
	range = exp(theta[1])
	# variance due to nugget
	nugget = exp(theta[2])
	# variance due to partial sill (spatially structured)
	psill = exp(theta[3])
	n = dim(X)[1]
	p = dim(X)[2]
	cormat = psill*spatial_model(distmat/range)
	diag(cormat) = psill + nugget
	cormatinv = solve(cormat)
	cormatinvX = cormatinv %*% X
	cormatinvY = cormatinv %*% y
	bhat = solve(t(X) %*% cormatinvX, t(X) %*% cormatinvY)
	r = y - X %*% bhat
	L1 = as.numeric(determinant(cormat, logarithm = TRUE)$modulus)
	L2 = t(r) %*% (cormatinv %*% r)
	L3 = as.numeric(determinant(t(X) %*% cormatinvX, logarithm = TRUE)$modulus)
	if(MLmeth == 'MLE') return(L1 + L2 + n*log(2*pi))
	if(MLmeth == 'REMLE') return(L1 + L2 + L3 + (n - p)*log(2*pi))
}

#minus two times the profiled loglikelihood
m2LLprof = function(theta, extrap = NULL, y, X, distmat, spatial_model, MLmeth)
{
	# range parameter
	range = exp(theta[1])
	# proportion of variance due to nugget
	propnug = exp(theta[2])/(1 + exp(theta[2]))
	n = dim(X)[1]
	p = dim(X)[2]
	cormat = (1 - propnug)*spatial_model(distmat/range, extrap)
	diag(cormat) = 1
	cormatinv = solve(cormat)
	cormatinvX = cormatinv %*% X
	cormatinvY = cormatinv %*% y
	bhat = solve(t(X) %*% cormatinvX, t(X) %*% cormatinvY)
	r = y - X %*% bhat
	L1 = as.numeric(determinant(cormat, logarithm = TRUE)$modulus)
	L2 = t(r) %*% (cormatinv %*% r)
	L3 = as.numeric(determinant(t(X) %*% cormatinvX, logarithm = TRUE)$modulus)
	if(MLmeth == 'MLE') return(L1 + n*log(L2) + n + n*log(2*pi/n))
	if(MLmeth == 'REMLE') return(L1 + (n - p)*log(L2) + L3 +
		(n - p) + (n - p)*log(2*pi/(n-p)))
}

M2LLprof_parms = function(optM2LLout, y, X, distmat, spatial_model, MLmeth)
{
	range = exp(optM2LLout$par[1])
	# proportion of variance due to nugget
	propnug = expit(optM2LLout$par[2])
	n = dim(X)[1]
	p = dim(X)[2]
	cormat = (1 - propnug)*spatial_model(distmat/range)
	diag(cormat) = 1
	cormatinv = solve(cormat)
	cormatinvX = cormatinv %*% X
	cormatinvY = cormatinv %*% y
	bhat = solve(t(X) %*% cormatinvX, t(X) %*% cormatinvY)
	r = y - X %*% bhat
	if(MLmeth == 'MLE') var_est = t(r) %*% (cormatinv %*% r)/(n)
	if(MLmeth == 'REMLE') var_est = t(r) %*% (cormatinv %*% r)/(n - p)
	list(var_est = var_est, beta_est = bhat)
}

M2LL_parms = function(optM2LLout, y, X, distmat, spatial_model, MLmeth)
{
	range = exp(optM2LLout$par[1])
	# variance due to nugget
	nugget = exp(optM2LLout$par[2])
	# variance due to partial sill (spatially structured)
	psill = exp(optM2LLout$par[3])
	# proportion of sill due to nugget
	propnug = nugget/(psill + nugget)
	n = dim(X)[1]
	p = dim(X)[2]
	cormat = (1 - propnug)*spatial_model(distmat/range)
	diag(cormat) = 1
	cormatinv = solve(cormat)
	cormatinvX = cormatinv %*% X
	cormatinvY = cormatinv %*% y
	bhat = solve(t(X) %*% cormatinvX, t(X) %*% cormatinvY)
	r = y - X %*% bhat
	if(MLmeth == 'MLE') var_est = t(r) %*% (cormatinv %*% r)/(n)
	if(MLmeth == 'REMLE') var_est = t(r) %*% (cormatinv %*% r)/(n - p)
	list(var_est = var_est, beta_est = bhat)
}

LOO_crossvalidation = function(optM2LLout, M2LLprof = TRUE, y, X, distmat, 
	spatial_model)
{
	if(M2LLprof == TRUE) {
		# range parameter
		range = exp(optM2LLout$par[1])
		# proportion of variance due to nugget
		propnug = exp(optM2LLout$par[2])/(1 + exp(optM2LLout$par[2]))
		n = dim(X)[1]
		p = dim(X)[2]
		V = (1 - propnug)*spatial_model(distmat/range)
		diag(V) = 1
	}
	if(M2LLprof == FALSE) {
		range = exp(optM2LLout$par[1])
		nugget = exp(optM2LLout$par[2])
		psill = exp(optM2LLout$par[3])
		n = dim(X)[1]
		p = dim(X)[2]
		V = psill*spatial_model(distmat/range)
		diag(V) = psill + nugget
	}
  z <- y
  Vi <- solve(V)
  cdd.out <- matrix(-999.9, nrow = n, ncol = 3)
  cdd.out[,1] <- 1:n
	for(i in 1:n) {
		Vi.i <- Vi[(1:n) != i,(1:n) != i] -
			matrix(Vi[(1:n) != i,i],ncol = 1) %*%
			matrix(Vi[i,(1:n) != i],nrow = 1)/Vi[i,i]
		c.i <- matrix(V[(1:n) != i,i],ncol = 1)
		xi <- matrix(X[i,], ncol = 1)
		X.i <- X[(1:n) != i,]
		z.i <- matrix(z[(1:n) != i], ncol = 1)
		xxi <- xi - t(X.i) %*% Vi.i %*% c.i
		covb.i <- solve(t(X.i) %*% Vi.i %*% X.i)
		si <- V[i,i]  - t(c.i) %*% Vi.i %*% c.i
		lam <- t(c.i + X.i %*% covb.i %*% xxi) %*% Vi.i

		cdd.out[i,2] <- lam %*% z.i
		cdd.out[i,3] <- sqrt(si + t(xxi) %*% covb.i %*% xxi)

	}
	cdd.out <- as.data.frame(cdd.out)
	names(cdd.out) <- c("y.row","cv.pred","cv.se")
	cdd.out

}

LOO_crossvalidation_usingLM = function(formula, DF)
{
	n = dim(DF)[1]
	cdd.out = matrix(-999.9, nrow = n, ncol = 3)
	lmout = lm(formula, data = DF)
	predout = predict(lmout, se.fit = TRUE)
	cdd.out[,1] = 1:n
	cdd.out[,2] = predout$fit
	cdd.out[,3] = sqrt(predout$se.fit^2 + predout$residual.scale^2)
	cdd.out <- as.data.frame(cdd.out)
	names(cdd.out) <- c("y.row","cv.pred","cv.se")
	cdd.out
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#      Fit models using MLE and compare profiled to unprofiled likelihoods
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# create fixed effects design matrices for up to 5th order polynomial
X0 = model.matrix(y ~ 1, data = DF)
X1 = model.matrix(y ~ poly(easting, northing, degree = 1, raw = TRUE), data = DF)
X2 = model.matrix(y ~ poly(easting, northing, degree = 2, raw = TRUE), data = DF)
X3 = model.matrix(y ~ poly(easting, northing, degree = 3, raw = TRUE), data = DF)
X4 = model.matrix(y ~ poly(easting, northing, degree = 4, raw = TRUE), data = DF)
X5 = model.matrix(y ~ poly(easting, northing, degree = 5, raw = TRUE), data = DF)

# create distance matrix
distmat = as.matrix(dist(DF[,c('easting','northing')]))
# asign the spatial model
spatial_model = exponential_spatial_model
# get response variable
y = DF$y

# inital covariance values for profiled likelihood
theta = c(0, logit(.2/.6))
# minimization
optM2LLprof = optim(theta, m2LLprof, y = y, X = X2, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'MLE')
# minimized value of the minus 2 times profiled loglikelihood
optM2LLprof$value
# MLE for the range parameter
exp(optM2LLprof$par[1])
# MLE for the nugget to total sill ratio
expit(optM2LLprof$par[2])
# MLE for covariance parameters and fixed parameters
M2LLprof_parms(optM2LLprof, y = y, X = X2, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'MLE')
	
# inital covariance values for regular likelihood
theta = c(0, log(.2), log(.4))
# minimization
optM2LL = optim(theta, m2LL, y = y, X = X2, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'MLE')
# minimized value of the minus 2 times loglikelihood
optM2LL$value
# MLE for the range parameter
exp(optM2LL$par[1])
# MLE for the nugget parameter
exp(optM2LL$par[2])
# MLE for the partial sill parameter
exp(optM2LL$par[3])
# MLE for sill
exp(optM2LL$par[2]) + exp(optM2LL$par[3])
# MLE for the nugget to total sill ratio
exp(optM2LL$par[2])/(exp(optM2LL$par[2]) + exp(optM2LL$par[3]))
# MLE for covariance parameters and fixed parameters
M2LL_parms(optM2LL, y = y, X = X2, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'MLE')

#use R package spmodel
library(spmodel)
# create a covariance parameter object
cov_params_val = cov_params('exponential', de = .4, ie = .2, range = 1)
# create an initial value object from covariance parameters
cov_initial_vals = cov_initial(cov_params_val)
# use MLE to fit the model
slmmout = slmm(y ~ poly(easting, northing, degree = 2, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, cov_type = "exponential",
	cov_initial = cov_initial_vals, estmethod = "ml")
# show the output
summary(slmmout)
# check the minimized -2 loglikelihood
-2*logLik(slmmout)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#     Fit models using REMLE and compare profiled to unprofiled likelihoods
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# inital covariance values for profiled restricted likelihood
theta = c(0, logit(.2/.6))
# minimization
optM2LLprof = optim(theta, m2LLprof, y = y, X = X2, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
# minimized value of the minus 2 times profiled restricted loglikelihood
optM2LLprof$value
# REMLE for the range parameter
exp(optM2LLprof$par[1])
# REMLE for the nugget to total sill ratio
expit(optM2LLprof$par[2])
# REMLE for covariance parameters and fixed parameters
M2LLprof_parms(optM2LLprof, y = y, X = X2, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
	
# inital covariance values for restricted likelihood
theta = c(0, log(.2), log(.4))
# minimization
optM2LL = optim(theta, m2LL, y = y, X = X2, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
# minimized value of the minus 2 times restricted loglikelihood
optM2LL$value
# REMLE for the range parameter
exp(optM2LL$par[1])
# REMLE for the nugget
exp(optM2LL$par[2])
# REMLE for the partial sill
exp(optM2LL$par[3])
# REMLE for the total sill
exp(optM2LL$par[2]) + exp(optM2LL$par[3])
# REMLE for the nugget to total sill ratio
exp(optM2LL$par[2])/(exp(optM2LL$par[2]) + exp(optM2LL$par[3]))
# REMLE for covariance parameters and fixed parameters
M2LL_parms(optM2LL, y = y, X = X2, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')

#use spmodel
# create a covariance parameter object
cov_params_val = cov_params('exponential', de = .4, ie = .2, range = 1)
# create an initial value object from covariance parameters
cov_initial_vals = cov_initial(cov_params_val)
# use REMLE to fit the model
slmmout = slmm(y ~ poly(easting, northing, degree = 2, raw = TRUE), data = DF, 
	xcoord = easting, ycoord = northing, cov_type = "exponential",
	cov_initial = cov_initial_vals, estmethod = "reml")
# show the output
summary(slmmout)
# check the minimized -2 loglikelihood
-2*logLik(slmmout)


################################################################################
#-------------------------------------------------------------------------------
#                                  AIC
#-------------------------------------------------------------------------------
################################################################################

# Independence Models for up to 5th Order Polynomial
fit_0 = lm(y ~ 1, data = DF)
fit_1 = lm(y ~ poly(easting, northing, degree = 1, raw = TRUE), data = DF)
fit_2 = lm(y ~ poly(easting, northing, degree = 2, raw = TRUE), data = DF)
fit_3 = lm(y ~ poly(easting, northing, degree = 3, raw = TRUE), data = DF)
fit_4 = lm(y ~ poly(easting, northing, degree = 4, raw = TRUE), data = DF)
fit_5 = lm(y ~ poly(easting, northing, degree = 5, raw = TRUE), data = DF)

# use profiling and fit all 6 models with both MLE and REMLE
optM2LLprof_X0_MLE = optim(theta, m2LLprof, y = y, X = X0, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'MLE')
optM2LLprof_X1_MLE = optim(theta, m2LLprof, y = y, X = X1, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'MLE')
optM2LLprof_X2_MLE = optim(theta, m2LLprof, y = y, X = X2, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'MLE')
optM2LLprof_X3_MLE = optim(theta, m2LLprof, y = y, X = X3, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'MLE')
optM2LLprof_X4_MLE = optim(theta, m2LLprof, y = y, X = X4, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'MLE')
optM2LLprof_X5_MLE = optim(theta, m2LLprof, y = y, X = X5, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'MLE')
optM2LLprof_X0_REMLE = optim(theta, m2LLprof, y = y, X = X0, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
optM2LLprof_X1_REMLE = optim(theta, m2LLprof, y = y, X = X1, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
optM2LLprof_X2_REMLE = optim(theta, m2LLprof, y = y, X = X2, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
optM2LLprof_X3_REMLE = optim(theta, m2LLprof, y = y, X = X3, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
optM2LLprof_X4_REMLE = optim(theta, m2LLprof, y = y, X = X4, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
optM2LLprof_X5_REMLE = optim(theta, m2LLprof, y = y, X = X5, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')

# compare all 6 MLE models using AIC
AIC_spatial = c(
	optM2LLprof_X0_MLE$value + 2,
	optM2LLprof_X1_MLE$value + 2*dim(X1)[2],
	optM2LLprof_X2_MLE$value + 2*dim(X2)[2],
	optM2LLprof_X3_MLE$value + 2*dim(X3)[2],
	optM2LLprof_X4_MLE$value + 2*dim(X4)[2],
	optM2LLprof_X5_MLE$value + 2*dim(X5)[2]
)
# and compare to independence models
AIC_indep = c(
	AIC(fit_0),
	AIC(fit_1),
	AIC(fit_2),
	AIC(fit_3),
	AIC(fit_4),
	AIC(fit_5)
)
# plot them
file_name = "SO4_AIC"

pdf(paste0(file_name,'.pdf'), width = 8.5, height = 8.5)
	old.par = par(mar = c(5,5,1,1))
	plot(0:5, AIC_indep, ylim = c(250, 500), type = 'l', lwd = 3, 
		xlab = 'Order of Polynomial', ylab = 'AIC', cex.axis = 1.5, cex.lab = 2)
	points(0:5, AIC_indep, pch = 19, cex = 3)
	lines(0:5, AIC_spatial,lty = 2, lwd = 3)
	points(0:5, AIC_spatial, pch = 1, cex = 3)
	legend(2.5, 500, legend = c('Indep Models','Spatial Models'), lty = c(1,2), 
		lwd = 2, pch = c(19,1), cex = 2)
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
#          Leave-one-out Cross-validation (LOOCV)
#-------------------------------------------------------------------------------
################################################################################

# compare all of the models using LOO crossvalidation
LOOCV_X0_MLE = LOO_crossvalidation(optM2LLprof_X0_MLE, M2LLprof = TRUE, y = y, 
	X = X0, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X0_MLE[,2])
mean((LOOCV_X0_MLE[,2] - y)^2)
LOOCV_X1_MLE = LOO_crossvalidation(optM2LLprof_X1_MLE, M2LLprof = TRUE, y = y, 
	X = X1, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X1_MLE[,2])
mean((LOOCV_X1_MLE[,2] - y)^2)
LOOCV_X2_MLE = LOO_crossvalidation(optM2LLprof_X2_MLE, M2LLprof = TRUE, y = y, 
	X = X2, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X2_MLE[,2])
mean((LOOCV_X2_MLE[,2] - y)^2)
LOOCV_X3_MLE = LOO_crossvalidation(optM2LLprof_X3_MLE, M2LLprof = TRUE, y = y, 
	X = X3, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X3_MLE[,2])
mean((LOOCV_X3_MLE[,2] - y)^2)
LOOCV_X4_MLE = LOO_crossvalidation(optM2LLprof_X4_MLE, M2LLprof = TRUE, y = y, 
	X = X4, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X4_MLE[,2])
mean((LOOCV_X4_MLE[,2] - y)^2)
LOOCV_X5_MLE = LOO_crossvalidation(optM2LLprof_X5_MLE, M2LLprof = TRUE, y = y,
	X = X5, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X5_MLE[,2])
mean((LOOCV_X5_MLE[,2] - y)^2)

LOOCV_X0_REMLE = LOO_crossvalidation(optM2LLprof_X0_REMLE, M2LLprof = TRUE, 
	y = y, X = X0, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X0_REMLE[,2])
mean((LOOCV_X0_REMLE[,2] - y)^2)
LOOCV_X1_REMLE = LOO_crossvalidation(optM2LLprof_X1_REMLE, M2LLprof = TRUE, 
	y = y, X = X1, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X1_REMLE[,2])
mean((LOOCV_X1_REMLE[,2] - y)^2)
LOOCV_X2_REMLE = LOO_crossvalidation(optM2LLprof_X2_REMLE, M2LLprof = TRUE, 
	y = y, X = X2, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X2_REMLE[,2])
mean((LOOCV_X2_REMLE[,2] - y)^2)
LOOCV_X3_REMLE = LOO_crossvalidation(optM2LLprof_X3_REMLE, M2LLprof = TRUE, 
	y = y, X = X3, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X3_REMLE[,2])
mean((LOOCV_X3_REMLE[,2] - y)^2)
LOOCV_X4_REMLE = LOO_crossvalidation(optM2LLprof_X4_REMLE, M2LLprof = TRUE, 
	y = y, X = X4, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X4_REMLE[,2])
mean((LOOCV_X4_REMLE[,2] - y)^2)
LOOCV_X5_REMLE = LOO_crossvalidation(optM2LLprof_X5_REMLE, M2LLprof = TRUE, 
	y = y, X = X5, distmat = distmat, spatial_model = exponential_spatial_model)
plot(y,LOOCV_X5_REMLE[,2])
mean((LOOCV_X5_REMLE[,2] - y)^2)

LOOCV_X0_IND = LOO_crossvalidation_usingLM(y ~ 1, DF)
plot(y,LOOCV_X0_IND[,2])
mean((LOOCV_X0_IND[,2] - y)^2)
LOOCV_X1_IND = LOO_crossvalidation_usingLM(y ~ poly(easting, northing, 
	degree = 1, raw = TRUE), DF)
plot(y,LOOCV_X1_IND[,2])
mean((LOOCV_X1_IND[,2] - y)^2)
LOOCV_X2_IND = LOO_crossvalidation_usingLM(y ~ poly(easting, northing, 
	degree = 2, raw = TRUE), DF)
plot(y,LOOCV_X2_IND[,2])
mean((LOOCV_X2_IND[,2] - y)^2)
LOOCV_X3_IND = LOO_crossvalidation_usingLM(y ~ poly(easting, northing, 
	degree = 3, raw = TRUE), DF)
plot(y,LOOCV_X3_IND[,2])
mean((LOOCV_X3_IND[,2] - y)^2)
LOOCV_X4_IND = LOO_crossvalidation_usingLM(y ~ poly(easting, northing, 
	degree = 4, raw = TRUE), DF)
plot(y,LOOCV_X4_IND[,2])
mean((LOOCV_X4_IND[,2] - y)^2)
LOOCV_X5_IND = LOO_crossvalidation_usingLM(y ~ poly(easting, northing, 
	degree = 5, raw = TRUE), DF)
plot(y,LOOCV_X5_IND[,2])
mean((LOOCV_X5_IND[,2] - y)^2)

LOOCV_spatial_MLE = c(
	mean((LOOCV_X0_MLE[,2] - y)^2),
	mean((LOOCV_X1_MLE[,2] - y)^2),
	mean((LOOCV_X2_MLE[,2] - y)^2),
	mean((LOOCV_X3_MLE[,2] - y)^2),
	mean((LOOCV_X4_MLE[,2] - y)^2),
	mean((LOOCV_X5_MLE[,2] - y)^2)
)
LOOCV_spatial_REMLE = c(
	mean((LOOCV_X0_REMLE[,2] - y)^2),
	mean((LOOCV_X1_REMLE[,2] - y)^2),
	mean((LOOCV_X2_REMLE[,2] - y)^2),
	mean((LOOCV_X3_REMLE[,2] - y)^2),
	mean((LOOCV_X4_REMLE[,2] - y)^2),
	mean((LOOCV_X5_REMLE[,2] - y)^2)
)
LOOCV_indep = c(
	mean((LOOCV_X0_IND[,2] - y)^2),
	mean((LOOCV_X1_IND[,2] - y)^2),
	mean((LOOCV_X2_IND[,2] - y)^2),
	mean((LOOCV_X3_IND[,2] - y)^2),
	mean((LOOCV_X4_IND[,2] - y)^2),
	mean((LOOCV_X5_IND[,2] - y)^2)
)

file_name = "SO4_LOOCV"

# plot them
pdf(paste0(file_name,'.pdf'), width = 8.5, height = 8.5)
	old.par = par(mar = c(5,5,1,1))
	plot(0:5, LOOCV_indep, ylim = c(0.19, 0.56), type = 'l', lwd = 3, 
		xlab = 'Order of Polynomial', ylab = 'RMSPE using LOOCV', cex.axis = 1.5, 
		cex.lab = 2)
	points(0:5, LOOCV_indep, pch = 19, cex = 3)
	lines(0:5, LOOCV_spatial_MLE, lty = 2, lwd = 3)
	points(0:5, LOOCV_spatial_MLE, pch = 1, cex = 3)
	lines(0:5, LOOCV_spatial_REMLE, lty = 3, lwd = 3)
	points(0:5, LOOCV_spatial_REMLE, pch = 2, cex = 3)
	legend(2.4, 0.57, legend = c('Indep Models','Spatial MLE', 'Spatial REMLE'), 
		lty = c(1, 2, 3), lwd = 2, pch = c(19, 1, 2), cex = 1.7)
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
#                       N-fold Cross-validation
#-------------------------------------------------------------------------------
################################################################################

# how many groups for N-fold cross-validation
ngrp = 3
# set seed so it is reproducible
set.seed(3005)
# copy the data.frame
DFgrp = DF
# create some random numbers
ru = runif(length(y))
# create breaks for the random numbers based on number of groups
brks = quantile(ru, probs = (0:ngrp)/(ngrp))
brks[1] = brks[1] - 1e-10
# add grp column to data based on random grouping
DFgrp$grp = as.factor(as.integer(cut(ru, 
	breaks = brks)))

# N-fold cross-validation based on independence model up to 5th-order polynomial
lm_Nfold_results = NULL
for(polydeg in 0:5) {
	hold = NULL
	for(i in 1:ngrp) {
		DFtrain = DF[DFgrp$grp != i,]
		if(polydeg == 0) {lmout = lm(y ~ 1, DFtrain) } else {
			lmout = lm(y ~ poly(easting, northing, degree = polydeg, raw = TRUE), 
				DFtrain)
		}
		hold = rbind(hold,
			cbind(DF[DFgrp$grp == i, 'y'],
				predict(lmout, newdata = DF[DFgrp$grp ==i,]))
			)
	}
	lm_Nfold_results = rbind(lm_Nfold_results, data.frame(polydeg = polydeg,
		RMSPE = sqrt(mean((hold[,1] - hold[,2])^2))))
}
lm_Nfold_results

# N-fold cross-validation based on spatial linear model up to 5th-order polynomial

xcoord = 'easting'
ycoord = 'northing'
Xset = list(X0, X1, X2, X3, X4, X5)
slmm_Nfold_results = NULL
for(polydeg in 0:5) {
	Xpoly = Xset[[polydeg + 1]]
	hold = NULL
	for(i in 1:ngrp) {
			DFtrain = DF[DFgrp$grp != i,]
			DFpred = DF[DFgrp$grp == i,]
			X = as.matrix(Xpoly[DFgrp$grp != i,])
			Xp = Xpoly[DFgrp$grp == i,]
			z = y[DFgrp$grp != i]
			yp = y[DFgrp$grp == i]
			no = dim(DFtrain)[1]
			np = dim(DFpred)[1]
			dismat <- as.matrix(dist(rbind(DFtrain[,c(xcoord, ycoord)],
				DFpred[,c(xcoord, ycoord)])))
			dist_oo = dismat[1:no, 1:no]
			dist_op = dismat[1:no, (no + 1):(no + np)]
			dist_pp = dismat[(no + 1):(no + np), (no + 1):(no + np)]
			# inital covariance values for restricted likelihood
			itheta = c(0, log(.2), log(.4))
			# minimization
			optM2LL = optim(itheta, m2LL, y = z, X = X, distmat = dist_oo, 
				spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
			# REMLE for the range parameter
			range = exp(optM2LL$par[1])
			# REMLE for the nugget
			nugget = exp(optM2LL$par[2])
			# REMLE for the partial sill
			parsil = exp(optM2LL$par[3])

			covMat <- parsil*spatial_model(dist_oo/range) +
					nugget*diag(no)
			Vpred <- parsil*spatial_model(dist_op/range)
			qrlist = qr(covMat, LAPACK = TRUE)
			ViX = solve(qrlist, X)
			Viz = solve(qrlist, z)
			ViVpred = solve(qrlist, Vpred)
			XViX <- crossprod(X, ViX)
			covb <- solve(XViX)
			bhat <- covb %*% crossprod(ViX, z)
			sill <- parsil + nugget
			
			preds <- matrix(NA, nrow = np, ncol = 3)
			preds[,1] = yp
			preds[,2] <- apply(as.vector(Viz) * Vpred, 2, sum) +
				Xp %*% bhat - t(Vpred) %*% (ViX %*% bhat)	
			preds[,3] <- sqrt(rep(sill, times = np) - 
				apply(ViVpred * Vpred, 2, sum) +
				apply((covb %*% t(Xp)) * t(Xp), 2, sum) -
				2*apply((covb %*% t(Xp)) * (t(X) %*% ViVpred), 2, sum) +
				apply(((covb %*% t(ViX)) %*% Vpred) * (t(X) %*% ViVpred), 2, sum))
			hold = rbind(hold, preds)
	}
	slmm_Nfold_results = rbind(slmm_Nfold_results, data.frame(polydeg = polydeg,
		RMSPE = sqrt(mean((hold[,1] - hold[,2])^2))))
}
lm_Nfold_results
slmm_Nfold_results

file_name = "SO4_Nfold_RMSPE"

pdf(paste0(file_name,'.pdf'), width = 8.5, height = 8.5)
	old.par = par(mar = c(5,5,1,1))
	plot(lm_Nfold_results, type = 'l', ylim = c(0.45, 1.5), lwd = 3,
		cex.axis = 1.5, cex.lab = 2, xlab = 'Order of Polynomial',
		ylab = 'RMSPE Using 3-Fold CV')
	points(lm_Nfold_results, pch = 19, cex = 3)
	lines(slmm_Nfold_results, lty = 2, lwd = 3)
	points(slmm_Nfold_results, pch = 1, cex = 3)
	legend(2, 1.5, legend = c('Indep Models','Spatial Models'), lty = c(1,2), 
		lwd = 2, pch = c(19,1), cex = 2)
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
#            Choosing Among Spatial Autocorrelation Models
#-------------------------------------------------------------------------------
################################################################################

# create distance matrix
distmat = as.matrix(dist(DF[,c('easting','northing')]))
# asign the spatial model
spatial_model = exponential_spatial_model
# get response variable
y = DF$y
# inital covariance values for profiled restricted likelihood
theta = c(0, logit(.2/.6))
# minimization
m2LL_X0_expCorr = optim(theta, m2LLprof, y = y, X = X0, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
m2LL_X0_sphCorr = optim(theta, m2LLprof, y = y, X = X0, distmat = distmat, 
	spatial_model = spherical_spatial_model, MLmeth = 'REMLE')
m2LL_X0_cirCorr = optim(theta, m2LLprof, y = y, X = X0, distmat = distmat, 
	spatial_model = circular_spatial_model, MLmeth = 'REMLE')
m2LL_X0_raqCorr = optim(theta, m2LLprof, y = y, X = X0, distmat = distmat, 
	spatial_model = rationalQuad_spatial_model, MLmeth = 'REMLE')
m2LL_X0_GauCorr = optim(theta, m2LLprof, y = y, X = X0, distmat = distmat, 
	spatial_model = Gaussian_spatial_model, MLmeth = 'REMLE')
m2LL_X4_expCorr = optim(theta, m2LLprof, y = y, X = X4, distmat = distmat, 
	spatial_model = exponential_spatial_model, MLmeth = 'REMLE')
m2LL_X4_sphCorr = optim(theta, m2LLprof, y = y, X = X4, distmat = distmat, 
	spatial_model = spherical_spatial_model, MLmeth = 'REMLE')
m2LL_X4_cirCorr = optim(theta, m2LLprof, y = y, X = X4, distmat = distmat, 
	spatial_model = circular_spatial_model, MLmeth = 'REMLE')
m2LL_X4_raqCorr = optim(theta, m2LLprof, y = y, X = X4, distmat = distmat, 
	spatial_model = rationalQuad_spatial_model, MLmeth = 'REMLE')
m2LL_X4_GauCorr = optim(theta, m2LLprof, y = y, X = X4, distmat = distmat, 
	spatial_model = Gaussian_spatial_model, MLmeth = 'REMLE')

m2LL_X0_expCorr$value
m2LL_X0_sphCorr$value
m2LL_X0_cirCorr$value
m2LL_X0_raqCorr$value
m2LL_X0_GauCorr$value

m2LL_X4_expCorr$value
m2LL_X4_sphCorr$value
m2LL_X4_cirCorr$value
m2LL_X4_raqCorr$value
m2LL_X4_GauCorr$value

################################################################################
#-------------------------------------------------------------------------------
#            Circular Model Range Parameters using X0
#-------------------------------------------------------------------------------
################################################################################

m2LL_X0_cirCorr = optim(theta, m2LLprof, y = y, X = X0, distmat = distmat, 
	spatial_model = circular_spatial_model, MLmeth = 'REMLE')$par[1]
m2LL_X1_cirCorr = optim(theta, m2LLprof, y = y, X = X1, distmat = distmat, 
	spatial_model = circular_spatial_model, MLmeth = 'REMLE')$par[1]
m2LL_X2_cirCorr = optim(theta, m2LLprof, y = y, X = X2, distmat = distmat, 
	spatial_model = circular_spatial_model, MLmeth = 'REMLE')$par[1]
m2LL_X3_cirCorr = optim(theta, m2LLprof, y = y, X = X3, distmat = distmat, 
	spatial_model = circular_spatial_model, MLmeth = 'REMLE')$par[1]
m2LL_X4_cirCorr = optim(theta, m2LLprof, y = y, X = X4, distmat = distmat, 
	spatial_model = circular_spatial_model, MLmeth = 'REMLE')$par[1]
m2LL_X5_cirCorr = optim(theta, m2LLprof, y = y, X = X5, distmat = distmat, 
	spatial_model = circular_spatial_model, MLmeth = 'REMLE')$par[1]

m2LL_X0_cirCorr
m2LL_X1_cirCorr
m2LL_X2_cirCorr
m2LL_X3_cirCorr
m2LL_X4_cirCorr
m2LL_X5_cirCorr

################################################################################
#-------------------------------------------------------------------------------
#        Visualize Profiled Likelihood for Autocorrelation Parameters
#-------------------------------------------------------------------------------
################################################################################

# create distance matrix
distmat = as.matrix(dist(DF[,c('easting','northing')]))
# asign the spatial model
spatial_model = circular_spatial_model
# get response variable
y = DF$y

# inital covariance values for profiled restricted likelihood
theta = c(0, logit(.2/.6))
# minimization
optM2LLprofA = optim(theta, m2LLprof, y = y, X = X0, distmat = distmat, 
	spatial_model = spatial_model, MLmeth = 'REMLE')
# minimized value of the minus 2 times profiled restricted loglikelihood
optM2LLprofA$value
# REMLE for the range parameter on log scale
optM2LLprofA$par[1]
# REMLE for the nugget to total sill ratio on logit scale
optM2LLprofA$par[2]
	
# compute minus 2 times loglikelihood on a grid around REMLE estimates
xx = (1:40)/40 - 1/40/2
yy = xx
xxA = (xx - 0.5)*0.5 + optM2LLprofA$par[1]
yyA = (yy - 0.5)*2 + optM2LLprofA$par[2]
zA = matrix(NA, nrow = 40, ncol = 40)
for(i in 1:40) {
	for(j in 1:40) {
		zA[i,j] = m2LLprof(c(xxA[i],yyA[j]), y = y, X = X0, distmat = distmat, 
			spatial_model = spatial_model, MLmeth = 'REMLE')
	}
}

# inital covariance values for profiled restricted likelihood
theta = c(0, logit(.2/.6))
# minimization
optM2LLprofB = optim(theta, m2LLprof, y = y, X = X4, distmat = distmat, 
	spatial_model = spatial_model, MLmeth = 'REMLE')
# minimized value of the minus 2 times profiled restricted loglikelihood
optM2LLprofB$value
# REMLE for the range parameter on log scale
optM2LLprofB$par[1]
# REMLE for the nugget to total sill ratio on logit scale
optM2LLprofB$par[2]
	
# compute minus 2 times loglikelihood on a grid around REMLE estimates
xxB = (xx - 0.5)*0.5 + optM2LLprofB$par[1]
yyB = (yy - 0.5)*2 + optM2LLprofB$par[2]
zB = matrix(NA, nrow = 40, ncol = 40)
for(i in 1:40) {
	for(j in 1:40) {
		zB[i,j] = m2LLprof(c(xxB[i],yyB[j]), y = y, X = X4, distmat = distmat, 
			spatial_model = spatial_model, MLmeth = 'REMLE')
	}
}

# inital covariance values for profiled restricted likelihood
theta = c(0, logit(.2/.6))
# minimization
optM2LLprofC = optim(theta, m2LLprof, y = y, X = X2, distmat = distmat, 
	spatial_model = spatial_model, MLmeth = 'REMLE')
# minimized value of the minus 2 times profiled restricted loglikelihood
optM2LLprofC$value
# REMLE for the range parameter on log scale
optM2LLprofC$par[1]
# REMLE for the nugget to total sill ratio on logit scale
optM2LLprofC$par[2]
	
# compute minus 2 times loglikelihood on a grid around REMLE estimates
xxC = (xx - 0.5)*0.5 + optM2LLprofC$par[1]
yyC = (yy - 0.5)*2 + optM2LLprofC$par[2]
zC = matrix(NA, nrow = 40, ncol = 40)
for(i in 1:40) {
	for(j in 1:40) {
		zC[i,j] = m2LLprof(c(xxC[i],yyC[j]), y = y, X = X2, distmat = distmat, 
			spatial_model = spatial_model, MLmeth = 'REMLE')
	}
}

# inital covariance values for profiled restricted likelihood
theta = c(.9, -2.2)
# minimization
optM2LLprofD = optim(theta, m2LLprof, y = y, X = X2, distmat = distmat, 
	spatial_model = spatial_model, MLmeth = 'REMLE', method = 'BFGS')
# minimized value of the minus 2 times profiled restricted loglikelihood
optM2LLprofD$value
# REMLE for the range parameter on log scale
optM2LLprofD$par[1]
# REMLE for the nugget to total sill ratio on logit scale
optM2LLprofD$par[2]
	
# compute minus 2 times loglikelihood on a grid around REMLE estimates
xxD = (xx - 0.5)*0.5 + optM2LLprof$par[1]
yyD = (yy - 0.5)*2 + optM2LLprof$par[2]
zD = matrix(NA, nrow = 40, ncol = 40)
for(i in 1:40) {
	for(j in 1:40) {
		zD[i,j] = m2LLprof(c(xxD[i],yyD[j]), y = y, X = X2, distmat = distmat, 
			spatial_model = spatial_model, MLmeth = 'REMLE')
	}
}

file_name = "SO4_Viz_m2LL_covParms"

pdf(paste0(file_name,'.pdf'), width = 12.5, height = 12.5)

	padj = -.5
	adj = -.15
	layout(matrix(1:4, ncol = 2, byrow = TRUE))
	old.par = par(mar = c(5,5,5,1))

	nbrks = 20
	brks = quantile(zA, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(xxA, yyA, zA, breaks = brks, col = cramp, main = 'Constant Mean Model', 
		cex.main = 2, xlab = 'log(range)', ylab = 'logit(proportion)',
		cex.axis = 1.5, cex.lab = 2)
	points(optM2LLprofA$par[1], optM2LLprofA$par[2], pch = 19, cex = 2, col = 'white')
	mtext('A', adj = adj, cex = 3, padj = padj)

	brks = quantile(zB, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(xxB, yyB, zB, breaks = brks, col = cramp, main = '4th Order Polynomial Model', 
		cex.main = 2, xlab = 'log(range)', ylab = 'logit(proportion)',
		cex.axis = 1.5, cex.lab = 2)
	points(optM2LLprofB$par[1], optM2LLprofB$par[2], pch = 19, cex = 2, col = 'white')
	mtext('B', adj = adj, cex = 3, padj = padj)

	brks = quantile(zC, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(xxC, yyC, zC, breaks = brks, col = cramp, main = '2nd Order Polynomial Model', 
		cex.main = 2, xlab = 'log(range)', ylab = 'logit(proportion)',
		cex.axis = 1.5, cex.lab = 2)
	points(optM2LLprofC$par[1], optM2LLprofC$par[2], pch = 19, cex = 2, col = 'white')
	mtext('C', adj = adj, cex = 3, padj = padj)

	brks = quantile(zD, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	image(xxD, yyD, zD, breaks = brks, col = cramp, main = '2nd Order Polynomial Model', 
		cex.main = 2, xlab = 'log(range)', ylab = 'logit(proportion)',
		cex.axis = 1.5, cex.lab = 2)
	points(optM2LLprofD$par[1], optM2LLprofD$par[2], pch = 19, cex = 2, col = 'white')
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
#                            The Matern Model
#-------------------------------------------------------------------------------
################################################################################

# create distance matrix
distmat = as.matrix(dist(DF[,c('easting','northing')]))
# get response variable
y = DF$y

# compute minus 2 times loglikelihood on a grid around REMLE estimates
xx = (1:40)/40 - 1/40/2
yy = xx
xx = (xx - 0.5)*0.5 + optM2LLprof$par[1]
yy = (yy - 0.5)*0.5 + optM2LLprof$par[2]
z = matrix(NA, nrow = 40, ncol = 40)
for(i in 1:40) {
	for(j in 1:40) {
		z[i,j] = m2LLprof(c(xx[i],yy[j]), extrap = 1, y = y, X = X0, 
			distmat = distmat, spatial_model = Matern_spatial_model, MLmeth = 'REMLE')
	}
}

# check for multimodality, and that we are at the minimum
nbrks = 20
brks = quantile(z, probs = (0:nbrks)/nbrks)
cramp = viridis(nbrks)
image(xx, yy, z, breaks = brks, col = cramp, main = 'Constant Mean Model', 
	cex.main = 2, xlab = 'log(range)', ylab = 'logit(proportion)',
	cex.axis = 1.5, cex.lab = 2)
points(optM2LLprof$par[1], optM2LLprof$par[2], pch = 19, cex = 2, col = 'white')

extrap_set = c(0.2, 0.5, 0.7, 1, 1.5, 2.5, 4, 8)
hold_m2LL = hold_logrange = hold_logitprop = NULL
for(extrap in extrap_set) {
	# inital covariance values for profiled restricted likelihood
	theta = c(0.5, logit(.01/.6))
	# minimization
	optM2LLprof = optim(theta, m2LLprof, extrap = extrap, y = y, X = X0, 
		distmat = distmat, spatial_model = Matern_spatial_model, MLmeth = 'REMLE')
	# minimized value of the minus 2 times profiled restricted loglikelihood
	hold_m2LL = c(hold_m2LL, optM2LLprof$value)
	# REMLE for the range parameter on log scale
	hold_logrange = c(hold_logrange, optM2LLprof$par[1])
	# REMLE for the nugget to total sill ratio on logit scale
	hold_logitprop = c(hold_logitprop, optM2LLprof$par[2])
}

# inital covariance values for profiled restricted likelihood
theta = c(0, logit(.1/.6))
# minimization
optM2LLprof = optim(theta, m2LLprof, extrap = 1, y = y, X = X0, 
	distmat = distmat, spatial_model = Matern_spatial_model, MLmeth = 'REMLE')
# minimized value of the minus 2 times profiled restricted loglikelihood
optM2LLprof$value
# REMLE for the range parameter on log scale
optM2LLprof$par[1]
# REMLE for the nugget to total sill ratio on logit scale
optM2LLprof$par[2]

file_name = "Matern_SO4"

pdf(paste0(file_name,'.pdf'), width = 15, height = 7.5)

	layout(matrix(c(1,2,3,4,4,4), ncol = 2))

	par(mar = c(0,6,1,1))
	plot(1:8, hold_m2LL, type = 'l', lwd = 3, xaxt = 'n', xlab = '', yaxt = 'n',
		ylim = c(297,309), ylab = '-2*loglikelihood', cex.axis = 1.5, cex.lab = 2.3)
	points(1:8, hold_m2LL, pch = 19, cex = 3)
	lines(c(1,8), c(hold_m2LL[4]+ 3.84,hold_m2LL[4]+ 3.84), lty = 2, lwd = 3)
	axis(2, at = c(298, 303, 308), cex.axis = 2)
	par(mar = c(0,6,0,1))
	plot(1:8, hold_logrange, type = 'l', lwd = 3, xaxt = 'n', xlab = '', 
		yaxt = 'n', ylim = c(-2,5), ylab = 'log(range)', cex.axis = 1.5, 
		cex.lab = 2.3)
	points(1:8, hold_logrange, pch = 19, cex = 3)
	axis(2, at = c(-1, 1, 3), cex.axis = 2)
	par(mar = c(6,6,0,1))
	plot(1:8, hold_logitprop, type = 'l', lwd = 3, xaxt = 'n', yaxt = 'n',
		ylim = c(-5,0), ylab = 'logit(proportion)', cex.axis = 1.5, cex.lab = 2.3,
		xlab = 'Smoothness Parameter')
	points(1:8, hold_logitprop, pch = 19, cex = 3)
	axis(2, at = c(-5, -3, -1), cex.axis = 2)
	axis(1, at = 1:8, labels = extrap_set, cex.axis = 2)

	nbrks = 20
	brks = quantile(z, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	par(mar = c(6,6,1,1))
	image(xx, yy, z, breaks = brks, col = cramp,
		cex.main = 2, xlab = 'log(range)', ylab = 'logit(proportion)',
		cex.axis = 2, cex.lab = 2.3)
	points(optM2LLprof$par[1], optM2LLprof$par[2], pch = 19, cex = 2, col = 'white')

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

file_name = "SO4_Prediction_Maps"

pdf(paste0(file_name,'.pdf'), width = 16, height = 10)

#         Independence Model

# predict points in DFpred using independence model with 5th order polynomial
pred_IND_X4 = predict(fit_4, DFpred, se.fit = TRUE)
# make a color map of predictions
cip = classIntervals(pred_IND_X4$fit, n = 6, style = 'fisher')
palp = viridis(6)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
source('addBreakColorLegend.R')
layout(matrix(1:8, nrow = 2, byrow = TRUE), widths = c(3,1,3,1,3,1,3,1))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
text(-2350000, 3200000, 'A', cex = 6)
old.par = par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = 2)

# make a color map of prediction standard errors
cip = classIntervals(sqrt(pred_IND_X4$se.fit^2 + pred_IND_X4$residual.scale^2), 
	n = 6, style = 'fisher')
palp = viridis(6)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
text(-2350000, 3200000, 'B', cex = 6)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = 2)

#         Spatial Model

xcoord = 'easting'
ycoord = 'northing'
no = dim(DF)[1]
np = dim(DFpred)[1]
X = X0
Xp = as.matrix(rep(1, times = np))
dismat <- as.matrix(dist(rbind(DF[,c(xcoord, ycoord)],
		DFpred[,c(xcoord, ycoord)])))	
dist_oo = dismat[1:no, 1:no]
dist_op = dismat[1:no, (no + 1):(no + np)]
# inital covariance values for restricted likelihood
itheta = c(0, log(.2), log(.4))
# minimization
optM2LL = optim(itheta, m2LL, y = y, X = X0, distmat = dist_oo, 
		spatial_model = circular_spatial_model, MLmeth = 'REMLE')
# REMLE for the range parameter
range = exp(optM2LL$par[1])
# REMLE for the nugget
nugget = exp(optM2LL$par[2])
# REMLE for the partial sill
parsil = exp(optM2LL$par[3])

covMat <- parsil*spatial_model(dist_oo/range) +
		nugget*diag(no)
Vpred <- parsil*spatial_model(dist_op/range)
qrlist = qr(covMat, LAPACK = TRUE)
ViX = solve(qrlist, X)
Viz = solve(qrlist, y)
ViVpred = solve(qrlist, Vpred)
XViX <- crossprod(X, ViX)
covb <- solve(XViX)
bhat <- covb %*% crossprod(ViX, y)
sill <- parsil + nugget
			
preds <- matrix(NA, nrow = np, ncol = 2)
preds[,1] <- apply(as.vector(Viz) * Vpred, 2, sum) +
	Xp %*% bhat - t(Vpred) %*% (ViX %*% bhat)	
preds[,2] <- sqrt(rep(sill, times = np) - 
	apply(ViVpred * Vpred, 2, sum) +
	apply((covb %*% t(Xp)) * t(Xp), 2, sum) -
	2*apply((covb %*% t(Xp)) * (t(X) %*% ViVpred), 2, sum) +
	apply(((covb %*% t(ViX)) %*% Vpred) * (t(X) %*% ViVpred), 2, sum))

cip = classIntervals(preds[,1], n = 6, style = 'fisher')
palp = viridis(6)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
text(-2350000, 3200000, 'C', cex = 6)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = 2)

cip = classIntervals(preds[,2], n = 6, style = 'fisher')
palp = viridis(6)
cip_colors = findColours(cip, palp)
old.par = par(mar = c(0,0,5,0))
plot(preds_sub, col = cip_colors, pch = 19, cex = 1.5)
plot(USboundary, add = TRUE, border = 'black')
text(-2350000, 3200000, 'D', cex = 6)
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = 2)

	par(old.par)
	
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
