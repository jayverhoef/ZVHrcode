sec_path = 'Rcode/Chapter11/Section 11.6/'
setwd(paste0(SLEDbook_path,sec_path))

# attach data library
library(ZVHdata)
library(sp)
library(viridis)
library(classInt)
library(colorspace)
library(spdep)
library(spmodel)
library(xtable)

################################################################################
#-------------------------------------------------------------------------------
#                     Caribou Example
#-------------------------------------------------------------------------------
################################################################################

# load data for graphics and analysis
data(caribouDF)

Xall = model.matrix(z ~ water*tarp, data = caribouDF)
X0 = Xall[,1]
X1 = Xall[,1:2]
X2 = Xall[,1:4]
X3 = Xall[,1:6]

# Table 11.5

# fit exponential geostatistical model
spfitEXP = splm(z ~ water*tarp, data = caribouDF, xcoord = 'x', ycoord = 'y', 
	spcov_type = 'exponential', estmethod = 'reml')
smryEXP = summary(spfitEXP)

totvar = coef(spfitEXP, type = 'spcov')['de'] + 
	coef(spfitEXP, type = 'spcov')['ie']
propnug = coef(spfitEXP, type = 'spcov')['ie']/totvar

Sigma = (1-propnug)*
	exp(-as.matrix(dist(caribouDF[,c('x','y')]))/
	coef(spfitEXP, type = 'spcov')['range']) + 
	propnug*diag(30)
Sigma_mhalf = eigen(Sigma)$vectors %*% diag(1/sqrt(eigen(Sigma)$values)) %*%
	t(eigen(Sigma)$vectors)
ystar = Sigma_mhalf %*% caribouDF$z
X0star = Sigma_mhalf %*% X0
X1star = Sigma_mhalf %*% X1
X2star = Sigma_mhalf %*% X2
X3star = Sigma_mhalf %*% X3

Px0 = X0star %*% solve(t(X0star) %*% X0star, t(X0star))
Px1 = X1star %*% solve(t(X1star) %*% X1star, t(X1star))
Px2 = X2star %*% solve(t(X2star) %*% X2star, t(X2star))
Px3 = X3star %*% solve(t(X3star) %*% X3star, t(X3star))

SS0 = t(ystar) %*% Px0 %*% ystar
SS1 = t(ystar) %*% (Px1 - Px0) %*% ystar
SS2 = t(ystar) %*% (Px2 - Px1) %*% ystar
SS3 = t(ystar) %*% (Px3 - Px2) %*% ystar
SSE = t(ystar) %*% (diag(30) - Px3) %*% ystar
SST = t(ystar) %*% ystar
MSE = SSE/24
# corrected total
SSCT = t(ystar) %*% ystar - SS0

# sums of squares
SSvec = c(SS1, SS2, SS3, SSE, SSCT)

# degrees of freedom
dfvec = c(
	sum(svd(Px1)$d > 1e-10) - 1,
	sum(svd(Px2)$d > 1e-10) - sum(svd(Px1)$d > 1e-10),
	sum(svd(Px3)$d > 1e-10) - sum(svd(Px2)$d > 1e-10),
	length(ystar) - sum(svd(Px3)$d > 1e-10),
	length(ystar) - 1
)

# mean square
MSvec = SSvec/dfvec
MSvec[5] = NA

# F-values
Fvals = rep(NA, times = length(dfvec))
Fvals[1:3] = MSvec[1:3]/MSvec[4]

# probabilities
Probvals = rep(NA, times = length(dfvec))
Probvals[1:3] = c(
	1 - pf(Fvals[1], df1 = dfvec[1], df2 = dfvec[4]),
	1 - pf(Fvals[2], df1 = dfvec[2], df2 = dfvec[4]),
	1 - pf(Fvals[3], df1 = dfvec[3], df2 = dfvec[4])
)

AOVtable = cbind(dfvec, SSvec, MSvec, Fvals, Probvals)
print(
    xtable(AOVtable, 
      align = c('l',rep('l', times = length(AOVtable[1,]))),
      digits = c(0,0,3,3,3,4)
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

R2 = sum(SSvec[1:3])/SSvec[5]
R2

ell = c(0, 1, 0, 0, 0, 0)
tval = t(ell) %*% coef(spfitEXP)/
	sqrt(MSE*t(ell) %*% solve(t(X3star) %*% X3star) %*% ell)
tval
2*(1 - pt(tval, df = 1))
1 - pchisq(tval^2, df = 1)

(0.143/0.099)^2

anova(spfitEXP)

ell = c(0, 0, -0.5, 1, 0, 0)
tval = t(ell) %*% coef(spfitEXP)/
	sqrt(MSE*t(ell) %*% solve(t(X3star) %*% X3star) %*% ell)
2*(1 - pt(tval, df = 1))
1 - pchisq(tval^2, df = 1)

ell = c(1, 0, 0, 1, 0, 0)

# Results in last paragraph of Section 11.6.3
drvlmEXPest = (ell %*% coef(spfitEXP) + mean(resid(spfitEXP)))
drvlmEXPse = sqrt(t(ell - apply(model.matrix(z ~ water*tarp, 
	data = caribouDF),2,mean)) %*% vcov(spfitEXP) %*% 
	(ell - apply(model.matrix(z ~ water*tarp, data = caribouDF),2,mean)))

lmEXPest = ell %*% coef(spfitEXP)
lmEXPse = sqrt(t(ell) %*% vcov(spfitEXP) %*% ell)
##

anova(spfitEXP)
spfitEXP1 = splm(z ~ water + tarp, data = caribouDF, xcoord = 'x', ycoord = 'y', 
	spcov_type = 'exponential', estmethod = 'ml')
summary(spfitEXP)
anova(spfitEXP, spfitEXP1)
# test as Chi-squared
1-pchisq((coef(spfit)[2]/sqrt(vcov(spfit)[2,2]))^2, df = 1)
# test as standard normal
2*(1-pnorm(coef(spfit)[2]/sqrt(vcov(spfit)[2,2])))
ell = rbind(c(0, 0, 1, 0, 0, 0), c(0, 0, 0, 1, 0, 0))

chisqval = t(ell %*% coef(spfit)) %*% solve(ell %*% vcov(spfit) %*% t(ell)) %*% 
	(ell %*% coef(spfit))
1-pchisq(chisqval, df = 2)

ell = c(1, 1, 0)

(ell %*% coef(spfitEXP) + mean(resid(spfitEXP)))/
sqrt(t(ell - apply(model.matrix(z ~ tarp, data = caribouDF),2,mean)) %*% 
	vcov(spfitEXP) %*% (ell - apply(model.matrix(z ~ tarp, data = caribouDF),2,mean)))

sqrt(t(ell - apply(model.matrix(z ~ tarp, data = caribouDF),2,mean)) %*% 
	vcov(spfitSAR) %*% (ell - apply(model.matrix(z ~ tarp, data = caribouDF),2,mean)))

bhatstar = solve(t(X3star) %*% X3star, t(X3star) %*% ystar)
resids = ystar - X3star %*% bhatstar
t(resids) %*% resids

ell = rbind(c(0, 0, 0, 0, 1, 0), c(0, 0, 0, 0, 0, 1))
anova(spfitEXP)
summary(spfitEXP)
t(ell %*% coef(spfitEXP)) %*% 
	solve(ell %*% solve(t(Xall) %*% solve(Sigma) %*% Xall) %*% t(ell)) %*%
	(ell %*% coef(spfitEXP))/(totvar)
vcov(spfitEXP)

file_name = "Caribou_means_interact"
pdf(paste0(file_name,'.pdf'), width = 6, height = 8)

	offset = .05
	par(mar = c(5,5,1,1))
	plot(c(1 - offset, 3 + offset), c(1.75, 2.45), type = 'n',
		 xaxt = 'n', xlab = '', ylab = '% Nitrogen', cex.lab = 2, cex.axis = 1.5)
	axis(1, at = c(1:3), labels = c('shade', 'clear', 'none'), 
		cex.axis = 2)
	trt = 1
	trp = 1
	points(trp - offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 19)
	lines(c(trp - offset,trp - offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	trt = 2
	trp = 2
	points(trp - offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 19)
	lines(c(trp - offset,trp - offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	trt = 3
	trp = 3
	points(trp - offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 19)
	lines(c(trp - offset,trp - offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	trt = 4
	trp = 1
	points(trp + offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 1)
	lines(c(trp + offset,trp + offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	trt = 5
	trp = 2
	points(trp + offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 1)
	lines(c(trp + offset,trp + offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	trt = 6
	trp = 3
	points(trp + offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 1)
	lines(c(trp + offset,trp + offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	lines(c((1:3) - offset), coef(spfit)[1:3], lty = 2, lwd = 3)
	lines(c((1:3) + offset), coef(spfit)[4:6], lty = 2, lwd = 3)
	legend(2,2.45, legend = c('Watered','Control'), pch = c(19,1), cex = 2)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))





lm_mains = lm(z ~ water + tarp, data = caribouDF)

car_mains = spautor(z ~ water + tarp, data = caribouDF, W = Nmat1, 'sar',
  estmethod = 'ml')

geo_mains = splm(z ~ water + tarp, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical',
	control = list(reltol = 1e-6), estmethod = 'ml')

lm_tarp = lm(z ~ tarp, data = caribouDF)

car_tarp = spautor(z ~ tarp, data = caribouDF, W = Nmat1, 'car',
  estmethod = 'ml')

geo_tarp = splm(z ~ tarp, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical',
	control = list(reltol = 1e-6), estmethod = 'ml')

lm_mean = lm(z ~ 1, data = caribouDF)

car_mean = spautor(z ~ 1, data = caribouDF, W = Nmat1, 'car',
  estmethod = 'ml')

geo_mean = splm(z ~ 1, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical',
	control = list(reltol = 1e-6), estmethod = 'ml')

#-------------------------------------------------------------------------------
#                    Testing anova() function
#-------------------------------------------------------------------------------

#behavior for single objects
anova(lm_all)
anova(car_tarp)
summary(car_tarp)

#behavior for 2 objects
anova(lm_all, lm_mains)

anova(car_all, car_mains)
car_all$optim$value
car_mains$optim$value

anova(geo_all, geo_mains)
geo_all$optim$value
geo_mains$optim$value

#behavior for multiple objects
anova(lm_all, lm_mains, lm_tarp, lm_mean)
#effect of order
anova(lm_all, lm_mean, lm_tarp, lm_mains)

anova(car_all, car_mains, car_tarp)
anova(geo_all, geo_mains, geo_tarp)

#effect of defaults
geo_all_def = splm(z ~ water*tarp, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical')
geo_mains_def = splm(z ~ water + tarp, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical')
anova(geo_all_def, geo_mains_def)
geo_all_def$optim$value
geo_mains_def$optim$value

# contrasts for geo model
colnames(model.matrix(geo_tarp))
tidy(anova(geo_tarp, L = list(c(0, 1, -1))))
#by hand
L = matrix(c(0, 1, -1),ncol = 1)
t(L) %*% geo_tarp$coefficients$fixed
sqrt(t(L) %*% geo_tarp$vcov$fixed %*% L)
tval = abs(as.numeric(t(L) %*% geo_tarp$coefficients$fixed/
	sqrt(t(L) %*% geo_tarp$vcov$fixed %*% L)))
tval
Fval = tval^2
Fval
# use t-distribution
2*(1 - pt(2.80963, 24))
# use standard normal
2*(1 - pnorm(2.80963))


# contrasts for car model
colnames(model.matrix(car_tarp))
tidy(anova(car_tarp, L = list(c(0, 1, -1))))
#by hand
L = matrix(c(0, 1, -1),ncol = 1)
t(L) %*% car_tarp$coefficients$fixed
sqrt(t(L) %*% car_tarp$vcov$fixed %*% L)
tval = abs(as.numeric(t(L) %*% car_tarp$coefficients$fixed/
	sqrt(t(L) %*% car_tarp$vcov$fixed %*% L)))
tval
Fval = tval^2
Fval
# use t-distribution
2*(1 - pt(2.80963, 24))
# use standard normal
2*(1 - pnorm(2.80963))


#-------------------------------------------------------------------------------
#
#           makeCovMat
#
#-------------------------------------------------------------------------------

#' make a CAR/SAR covariance matrix for modeling
#'
#' make a CAR/SAR covariance matrix for modeling
#'
#' @param theta covariance parameters, with overall variance parameter profiled out.
#' @param indComp an additive independent component to the model.  Default is TRUE.  
#' @param Nmat neighborhood matrix
#' @param distMat distance matrix
#' @param indSamp indicator vector for wheter location was sampled. Zero, or FALSE, indicates missing value
#' @param model either 'CAR' or 'SAR' to be the model associated with the Nmat argument
#' @param logical value on whether Nmat should be row-standardized
#' @param rhobound a vector of two elements containing bounds for rho.  This should be determined from the eigenvalues of Nmat prior to running function
#'
#' @return two times the negative log-likelihood
#'
#' @author Jay Ver Hoef jverhoef
#' @rdname m2LL
#' @export m2LL 

makeCovMat = function(theta, nN, indComp = TRUE, Nmat = NULL, 
  model = 'CAR', rowStand = TRUE, rhoBound = c(-1,1))
{
  nN = dim(distMat)[1]
  V = matrix(1, nrow = nN, ncol = nN )
  diag(V) = 0
  itheta = 0
  if(!is.null(Nmat)) {
    V = as(Nmat, 'sparseMatrix')
  }
  if(!is.null(model)) {
    itheta = itheta + 1
    rho = rhoBound[1] + .00005 + exp(theta[itheta])/
      (1 + exp(theta[itheta]))*.9999*(rhoBound[2] - rhoBound[1])
    attr(theta,'names')[itheta] = 'logitRho'
    rs = rep(1, times = nN)
    if(rowStand) rs = apply(V,1,sum)
    if(model == 'CAR')  V = diag(rs) - rho*V
    if(model == 'SAR') V = (diag(nN) - rho*(1/rs)*V) %*%
      (diag(nN) - rho*t((1/rs)*V))
  }
  if(indComp & is.null(Nmat)) {
    itheta = itheta + 1
    relEps = exp(theta[itheta])
    attr(theta,'names')[itheta] = 'relEps'
    V = V -  V %*% solve(V + diag(rep(relEps,times = nN)),V)
  }
  if(indComp & is.null(Nmat)) {
    V = diag(nN)
  }
  V
}

#-------------------------------------------------------------------------------
#
#           m2LL
#
#-------------------------------------------------------------------------------

#' two times the negative log-likelihood
#'
#' two times the negative log-likelihood
#'
#' @param theta covariance parameters, with overall variance parameter profiled out.
#' @param X design matrix for fixed effects
#' @param y vector of data for response variable
#' @param indComp an additive independent component to the model.  Default is TRUE.  
#' @param Nmat neighborhood matrix
#' @param distMat distance matrix
#' @param indSamp indicator vector for wheter location was sampled. Zero, or FALSE, indicates missing value
#' @param model either 'CAR' or 'SAR' to be the model associated with the Nmat argument
#' @param logical value on whether Nmat should be row-standardized
#' @param rhobound a vector of two elements containing bounds for rho.  This should be determined from the eigenvalues of Nmat prior to running function
#'
#' @return two times the negative log-likelihood
#'
#' @author Jay Ver Hoef jverhoef
#' @rdname m2LL
#' @export m2LL 

m2LL = function(theta, X, y, indComp = TRUE, Nmat = NULL, 
  model = 'CAR', rowStand = TRUE, rhoBound = c(-1,1), MLmeth = 'REMLE')
{
  if(any(abs(theta) > 10.1)) return(1e+32)

  ntheta = 0
  nN = length(y)
	Vi.oo = makeCovMat(theta = theta, nN = nN, indComp = indComp, Nmat = Nmat, 
		model = model, rowStand = rowStand, rhoBound = rhoBound)
  XVi = t(X) %*% Vi.oo
  covbi = XVi %*% X
  covb = solve(covbi)
  bHat = covb %*% XVi %*% y
  r = y - X %*% bHat
  n = length(y)
  p = length(X[1,])
  if(MLmeth == 'MLE') {
  m2LL = n*log(t(r) %*% Vi.oo %*% r) - 
    as.numeric(determinant(Vi.oo, logarithm = TRUE)$modulus) +
    n*(log(2*pi) + 1 - log(n))
	} else if(MLmeth == 'REMLE') {
  m2LL = (n-p)*log(t(r) %*% Vi.oo %*% r) - 
    as.numeric(determinant(Vi.oo, logarithm = TRUE)$modulus) +
    as.numeric(determinant(XVi %*% X, logarithm = TRUE)$modulus) +
    (n - p)*(log(2*pi) + 1 - log((n - p)))
	} else {return('MLmeth argument must be either MLE or REMLE')}
 
	attr(m2LL,'covParms') = theta
	as.numeric(m2LL)
}


X0 = as.matrix(model.matrix(z ~ 1, data = caribouDF))
X1 = as.matrix(model.matrix(z ~ water + tarp + water:tarp, data = caribouDF))
X2 = as.matrix(model.matrix(z ~ water + tarp, data = caribouDF))
X3 = as.matrix(model.matrix(z ~ water + tarp, data = caribouDF))
y = caribouDF$z

ntheta = 2

if(ntheta == 1) {
	# undebug(m2LL)
	# undebug(makeCovMat)
  optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
		Nmat = Nmat4, indComp = FALSE)
  theta = optOut$minimum
  m2LLargmin = optOut$objective
  } else {
	# undebug(m2LL)
	# undebug(makeCovMat)
  optOut = optim(rep(0, times = ntheta), m2LL, X = X1, y = y, 
		Nmat = Nmat1, indComp = TRUE, model = 'SAR')
  theta = optOut$par
  m2LLargmin = optOut$value
}
X = X1
Vi.oo = makeCovMat(theta, nN = length(y), Nmat = Nmat1, indComp = TRUE,
	model = 'SAR')
XVi = t(X) %*% Vi.oo
covbi = XVi %*% X
covb = solve(covbi)
bHat = covb %*% XVi %*% y
bHat
r = as.matrix(y - X %*% bHat)
n = length(y)
p = length(X[1,])
sigma = as.numeric((t(r) %*% Vi.oo %*% r)/n)
bHat_se = sqrt(sigma*diag(covb))
bHat_se
bHat/bHat_se

optOut = optim(rep(0, times = 3), m2LL, X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat4, distMat = distMat4, indComp = FALSE,
	MLmeth = 'MLE')
theta = optOut$par
m2LLargmin = optOut$value

X = X1
WMi = makeCovMat(theta, indSamp = indSamp, Nmat = Nmat4, distMat = distMat4,
	indComp = FALSE)
WMi.oo = WMi[indSamp,indSamp] 
WMi.uu = WMi[!indSamp,!indSamp]
WMi.uo = WMi[!indSamp,indSamp]
WMi.ou = WMi[indSamp,!indSamp]
Vi.oo = WMi.oo - WMi.ou %*% solve(WMi.uu, WMi.uo)
XVi = t(X) %*% Vi.oo
covbi = XVi %*% X
covb = solve(covbi)
bHat = covb %*% XVi %*% y
bHat
r = as.matrix(y - X %*% bHat)
n = length(y)
p = length(X[1,])
sigma = as.numeric((t(r) %*% Vi.oo %*% r)/n)
bHat_se = sqrt(sigma*diag(covb))
bHat_se

