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
library(gstat)
library(spmodel)

# load data for graphics and analysis
data(MOSSobs)
data(MOSSpreds)

DF = data.frame(MOSSobs@data, easting = MOSSobs@coords[,1]/1e+3,
	northing = MOSSobs@coords[,2]/1e+3)
DF$year = as.factor(DF$year)
DF$field_dup = as.factor(DF$field_dup)
DF$Zn = log(DF$Zn)
DF$Pb = log(DF$Pb)
DF$dist2road = log(DF$dist2road)


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
m2LL = function(theta, y, X, distmat1, distmat2, indx1, indx2,
	spatial_model, Z1, Z2, MLmeth = 'REMLE')
{
	# range parameter
	range = exp(theta[1])
	# variance due to nugget
	nugget = exp(theta[2])
	# variance due to partial sill (spatially structured)
	psill = exp(theta[3])
	# variance due to field duplication
	sigrep = exp(theta[4])
	# variance due to laboratory replication
	siglab = exp(theta[5])
	n = dim(X)[1]
	p = dim(X)[2]
	cormat = matrix(0, nrow = n, ncol = n)
	cormat1 = psill*spatial_model(distmat1/range)
	diag(cormat1) = psill + nugget
	cormat2 = psill*spatial_model(distmat2/range)
	diag(cormat2) = psill + nugget
	cormat[indx1,indx1] = cormat1
	cormat[indx2,indx2] = cormat2
	cormat = cormat + sigrep*(Z1 %*% t(Z1)) + siglab*(Z2 %*% t(Z2))
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

M2LL_parms = function(optM2LLout, y, X, distmat1, distmat2, indx1, indx2,
	spatial_model, Z1, Z2, MLmeth = 'REMLE')
{
	theta = optM2LLout$par
	# range parameter
	range = exp(theta[1])
	# variance due to nugget
	nugget = exp(theta[2])
	# variance due to partial sill (spatially structured)
	psill = exp(theta[3])
	# variance due to field duplication
	sigrep = exp(theta[4])
	# variance due to laboratory replication
	siglab = exp(theta[5])
	n = dim(X)[1]
	p = dim(X)[2]
	cormat = matrix(0, nrow = n, ncol = n)
	cormat1 = psill*spatial_model(distmat1/range)
	diag(cormat1) = psill + nugget
	cormat2 = psill*spatial_model(distmat2/range)
	diag(cormat2) = psill + nugget
	cormat[indx1,indx1] = cormat1
	cormat[indx2,indx2] = cormat2
	cormat = cormat + sigrep*(Z1 %*% t(Z1)) + siglab*(Z2 %*% t(Z2))
	cormatinv = solve(cormat)
	cormatinvX = cormatinv %*% X
	cormatinvY = cormatinv %*% y
	bhat = solve(t(X) %*% cormatinvX, t(X) %*% cormatinvY)
	covb = solve(t(X) %*% cormatinvX)
	list(var_est = nugget + psill + sigrep + siglab, beta_est = bhat, 
		covb = covb)
}

LOO_crossvalidation = function(optM2LLout, y, X, distmat1, distmat2, 
	indx1, indx2, spatial_model, Z1, Z2)
{
  theta = optM2LLout$par
  # range parameter
  range = exp(theta[1])
  # variance due to nugget
  nugget = exp(theta[2])
  # variance due to partial sill (spatially structured)
  psill = exp(theta[3])
  # variance due to field duplication
  sigrep = exp(theta[4])
  # variance due to laboratory replication
  siglab = exp(theta[5])
  n = dim(X)[1]
  p = dim(X)[2]
  cormat = matrix(0, nrow = n, ncol = n)
  cormat1 = psill*spatial_model(distmat1/range)
  diag(cormat1) = psill + nugget
  cormat2 = psill*spatial_model(distmat2/range)
  diag(cormat2) = psill + nugget
  cormat[indx1,indx1] = cormat1
  cormat[indx2,indx2] = cormat2
  V = cormat + sigrep*(Z1 %*% t(Z1)) + siglab*(Z2 %*% t(Z2))
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

X = model.matrix(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF)

Z1 = model.matrix( ~ -1 + sample, data = DF)
Z2 = model.matrix( ~ -1 + I(as.factor(paste(sample,field_dup))), data = DF)

# create distance matrix
indx1 = 1:244
indx2 = 245:365
distmat2001 = as.matrix(dist(DF[indx1,c('easting','northing')]))
distmat2006 = as.matrix(dist(DF[indx2,c('easting','northing')]))
# asign the spatial model
spatial_model = exponential_spatial_model
# get response variable
y = DF$Pb

# inital covariance values for likelihood
lmout = lm(Pb ~ year + dist2road + sideroad, data = DF)
lmresid = resid(lmout)
DFres = data.frame(resid = lmresid, x = DF$easting, y = DF$northing)
library(gstat)
vgm_resid_dir <- variogram(resid ~ 1,  
	loc=~x+y, data = DFres[indx1,], cutoff = 14, width = 10/10)
plot(vgm_resid_dir$dist, vgm_resid_dir$gamma, xlab = 'Distance (km)',
    ylab = 'Semivariogram', cex.lab = 2, cex.axis = 1.5, type = 'l',
    xlim = c(0,15), ylim = c(0,.4))
points(vgm_resid_dir$dist, vgm_resid_dir$gamma, pch = 19, cex = 3)
vgm_resid_dir <- variogram(resid ~ 1,  
	loc=~x+y, data = DFres[indx2,], cutoff = 14, width = 10/10)
lines(vgm_resid_dir$dist, vgm_resid_dir$gamma, lty = 2)
points(vgm_resid_dir$dist, vgm_resid_dir$gamma, pch = 1, cex = 3)

#initial values based on empirical semivariogram
thetai = c(log(10), log(.03), log(.1), log(.03), log(.03))  
# minimization
optM2LL = optim(thetai, m2LL, y = y, X = X, distmat1 = distmat2001, 
	distmat2 = distmat2006, indx1 = indx1, indx2 = indx2,
	spatial_model = spherical_spatial_model, Z1 = Z1, Z2 = Z2,
	MLmeth = 'REMLE')
# minimized value of the minus 2 times profiled loglikelihood
optM2LL$value
exp(optM2LL$par)

plot((0:70)/2,  
	exp(optM2LL$par[3])*spherical_spatial_model((0:70)/2/exp(optM2LL$par[1])), 
	type = 'l', lwd = 3, ylim = c(0,.4), ylab = 'Fitted Autocorrelation Model',
	xlab = 'Distance (km)')
points(0, exp(optM2LL$par[2]) + exp(optM2LL$par[4]) + exp(optM2LL$par[5]) + 
	exp(optM2LL$par[3]), pch = 19, cex = 2)

#undebug(splm)
#undebug(spmodel:::cov_estimate_gloglik_splm)
PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing', spcov_type = 'spherical',
	partition_factor = ~ year,
#	spcov_initial = spcov_initial('spherical', range = 4, known = c('range')),
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(PbOut)
PbOut$optim$value

min(DF$dist2road)


PbOut = splmm(Pb ~ year + dist2road + sideroad + sideroad*dist2road + 
	year*dist2road + year*sideroad + year*sideroad*dist2road, data = DF,
	spcov_type = 'none')
summary(PbOut)

summary(lm(Pb ~ year + dist2road + sideroad + sideroad*dist2road + 
	year*dist2road + year*sideroad + year*sideroad*dist2road, data = DF,
	spcov_type = 'none'))

lm(Pb ~ dist2road, data = DF[DF$year == '2001' & DF$sideroad == 'N',])
lm(Pb ~ dist2road, data = DF[DF$year == '2001' & DF$sideroad == 'S',])
lm(Pb ~ dist2road, data = DF[DF$year == '2006' & DF$sideroad == 'N',])
lm(Pb ~ dist2road, data = DF[DF$year == '2006' & DF$sideroad == 'S',])

est_parms = M2LL_parms(optM2LL, y = y, X = X, distmat1 = distmat2001, 
	distmat2 = distmat2006, indx1 = indx1, indx2 = indx2,
	spatial_model = spherical_spatial_model, Z1 = Z1, Z2 = Z2,
	MLmeth = 'MLE')
est_parms$beta_est

data.frame(est = est_parms$beta_est, 
	se = sqrt(diag(est_parms$covb)),
	t = est_parms$beta_est/sqrt(diag(est_parms$covb)))

PbOut = splmm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF1,
	xcoord = 'easting', ycoord = 'northing',  spcov_type = 'spherical',
	partition_factor = ~ year,
#	spcov_initial = spcov_initial('spherical', range = 0.5, known = c('extra')),
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(PbOut)
PbOut$optim$value
# there will be 4 fitted submodels.  The 3-way interaction specified a different
# slope for each year/side-of-road combination, while the year*sideroad interaction
# provides for 4 intercept.  Let's compute the 4 intercepts and slope:
int_2001N = PbOut$coefficient$fixed[c('(Intercept)')]
int_2001S = sum(PbOut$coefficient$fixed[c('(Intercept)','sideroadS')])
int_2006N = sum(PbOut$coefficient$fixed[c('(Intercept)','year2006')])
int_2006S = sum(PbOut$coefficient$fixed[c('(Intercept)','year2006',
	'sideroadS', 'year2006:sideroadS')])
slp_2001N = PbOut$coefficient$fixed[c('dist2road')]
slp_2001S = sum(PbOut$coefficient$fixed[c('dist2road','dist2road:sideroadS')])
slp_2006N = sum(PbOut$coefficient$fixed[c('dist2road','year2006:dist2road')])
slp_2006S = sum(PbOut$coefficient$fixed[c('dist2road','dist2road:sideroadS',
	'year2006:dist2road','year2006:dist2road:sideroadS')])
MOSS_data = DF
#X11()
old.par = par(mar = c(5,5,1,1))
  plot(
    c(-(MOSS_data[MOSS_data$year == 2001 & 
        MOSS_data$sideroad == 'N','dist2road']),
      (MOSS_data[MOSS_data$year == 2001 & 
        MOSS_data$sideroad == 'S','dist2road'])),
    c((MOSS_data[MOSS_data$year == 2001 & MOSS_data$sideroad == 'N','Pb']),
      (MOSS_data[MOSS_data$year == 2001 & MOSS_data$sideroad == 'S','Pb'])),
    xlab = 'Log Distance From Road Towards South', 
    ylab = 'Log Zinc Concentration', 
    pch = 19, cex = 1.2, cex.lab = 2, cex.axis = 1.5, ylim = c(-1,9)
  )
  points(
    c(-(MOSS_data[MOSS_data$year == 2006 & 
        MOSS_data$sideroad == 'N','dist2road']),
      (MOSS_data[MOSS_data$year == 2006 & 
        MOSS_data$sideroad == 'S','dist2road'])),
    c((MOSS_data[MOSS_data$year == 2006 & MOSS_data$sideroad == 'N','Pb']),
      (MOSS_data[MOSS_data$year == 2006 & MOSS_data$sideroad == 'S','Pb'])),
    pch = 3, cex = 2, lwd = 1.2
  )
  lines(c(0,-11),c(int_2001N, int_2001N + slp_2001N*(11)), lwd = 3, col = 'red')
  lines(c(0,-11),c(int_2006N, int_2006N + slp_2006N*(11)), lwd = 3, col = 'blue')
  lines(c(0,11),c(int_2001S, int_2001S + slp_2001S*(11)), lwd = 3, col = 'red')
  lines(c(0,11),c(int_2006S, int_2006S + slp_2006S*(11)), lwd = 3, col = 'blue')
  
  legend(6, 8.1, legend = c('2001', '2006'),
    pch = c(19, 3), cex = 2.5, lty = c(1,1), col = c('red','blue'))
par(old.par)
	
DF2006 = DF[DF$year == '2006',]
Pb06 = splmm(Pb ~ dist2road + sideroad + sideroad:dist2road, data = DF2006,
	spcov_type = 'spherical', xcoord = 'easting', ycoord = 'northing',
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(Pb06)
int_2006N = Pb06$coefficient$fixed[c('(Intercept)')]
int_2006S = sum(Pb06$coefficient$fixed[c('(Intercept)','sideroadS')])
slp_2006N = Pb06$coefficient$fixed[c('dist2road')]
slp_2006S = sum(Pb06$coefficient$fixed[c('dist2road','dist2road:sideroadS')])

lm06 = lm(Pb ~ dist2road + sideroad + sideroad:dist2road, data = DF2006)
summary(lm06)
sumLm06 = summary(lm06)
int_2006N = sumLm06$coefficients[c('(Intercept)'), 'Estimate']
int_2006S = sum(sumLm06$coefficient[c('(Intercept)','sideroadS'), 'Estimate'])
slp_2006N = sumLm06$coefficient[c('dist2road'), 'Estimate']
slp_2006S = sum(sumLm06$coefficient[c('dist2road','dist2road:sideroadS'), 'Estimate'])
X11()
  plot(
    c(-(DF2006[DF2006$sideroad == 'N','dist2road']),
      (DF2006[DF2006$sideroad == 'S','dist2road'])),
    c((DF2006[DF2006$sideroad == 'N','Pb']),
      (DF2006[DF2006$sideroad == 'S','Pb'])),
    xlab = 'Log Distance From Road Towards South', 
    ylab = 'Log Lead Concentration', 
    pch = 19, cex = 1.2, cex.lab = 2, cex.axis = 1.5, ylim = c(-1,9)
  )
  lines(c(0,-11),c(int_2006N, int_2006N + slp_2006N*(11)), lwd = 3, col = 'red')
  lines(c(0,11),c(int_2006S, int_2006S + slp_2006S*(11)), lwd = 3, col = 'red')

lm06 = lm(Pb ~ dist2road, data = DF2006[DF2006$sideroad == 'N',])
summary(lm06)
sumLm06 = summary(lm06)
int_2006N = sumLm06$coefficients[c('(Intercept)'), 'Estimate']
slp_2006N = sumLm06$coefficient[c('dist2road'), 'Estimate']
Pb06 = splmm(Pb ~ dist2road, data = DF2006[DF2006$sideroad == 'N',],
	spcov_type = 'spherical', xcoord = 'easting', ycoord = 'northing',
#	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-12), estmethod = 'reml')
int_2006N = Pb06$coefficient$fixed[c('(Intercept)')]
slp_2006N = Pb06$coefficient$fixed[c('dist2road')]
summary(Pb06)
X11()
  plot(
    c((DF2006[DF2006$sideroad == 'N','dist2road'])),
    c((DF2006[DF2006$sideroad == 'N','Pb'])),
    xlab = 'Log Distance From Road Towards South', 
    ylab = 'Log Lead Concentration', 
    pch = 19, cex = 1.2, cex.lab = 2, cex.axis = 1.5, ylim = c(-1,9)
  )
  lines(c(0,11),c(int_2006N, int_2006N + slp_2006N*(11)), lwd = 3, col = 'red')




















X = model.matrix(Zn ~ year + dist2road + sideroad, data = DF)
thetai = c(log(10), log(.2), log(.01), log(.01), log(.01))  
optM2LL1 = optim(thetai, m2LL, y = y, X = X, distmat1 = distmat2001, 
	distmat2 = distmat2006, indx1 = indx1, indx2 = indx2,
	spatial_model = exponential_spatial_model, Z1 = Z1, Z2 = Z2,
	MLmeth = 'MLE')
plot((0:70)/2,  
	exp(optM2LL1$par[3])*spherical_spatial_model((0:70)/2/exp(optM2LL1$par[1])), 
	type = 'l', lwd = 3, ylim = c(0,.4), ylab = 'Fitted Autocorrelation Model',
	xlab = 'Distance (km)')
points(0, exp(optM2LL1$par[2]) + exp(optM2LL1$par[4]) + exp(optM2LL1$par[5]) + 
	exp(optM2LL1$par[3]), pch = 19, cex = 2)

est_parms1 = M2LL_parms(optM2LL1, y = y, X = X, distmat1 = distmat2001, 
	distmat2 = distmat2006, indx1 = indx1, indx2 = indx2,
	spatial_model = exponential_spatial_model, Z1 = Z1, Z2 = Z2,
	MLmeth = 'MLE')
data.frame(est = est_parms1$beta_est, 
	se = sqrt(diag(est_parms1$covb)),
	t = est_parms1$beta_est/sqrt(diag(est_parms1$covb)))



################################################################################
#-------------------------------------------------------------------------------
#          Leave-one-out Cross-validation (LOOCV)
#-------------------------------------------------------------------------------
################################################################################


# compare all of the models using LOO crossvalidation
#undebug(LOO_crossvalidation)
LOOCVout = LOO_crossvalidation(optM2LL, y = y, X = X, distmat1 = distmat2001, 
	distmat2 = distmat2006, indx1 = indx1, indx2 = indx2,
	spatial_model = spherical_spatial_model, Z1 = Z1, Z2 = Z2)
plot(y,LOOCVout[,2])
mean((LOOCVout[,2] - y)^2)
std_LOO_resids = abs(LOOCVout[,2] - y)/sqrt(var(LOOCVout[,2] - y))
plot(std_LOO_resids)
plot(DF$easting, DF$northing)
points(DF[which(std_LOO_resids > 3), c('easting','northing')],
	pch = 19, col = 'red')
DF1 = DF[which(std_LOO_resids > 3),]

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',  spcov_type = 'spherical',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-12), estmethod = 'reml')
loocv_out = loocv(PbOut, cv_fitted = TRUE)
plot(y, loocv_out$cv_fitted)
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
