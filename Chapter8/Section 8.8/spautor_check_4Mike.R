library(spmodel)
library(sf)
library(sp)
library(spdep)

# some useful transformations
logit = function(x) {log(x/(1 - x))}
expit = function(x) {exp(x)/(1 + exp(x))}

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
#' @param Nmat neighborhood matrix
#' @param model either 'CAR' or 'SAR' to be the model associated with the Nmat argument
#' @param rowStand logical value on whether Nmat should be row-standardized
#' @param rhoBound a vector of two elements containing bounds for rho.  This should be determined from the eigenvalues of Nmat prior to running function
#'
#' @return A covariance matrix
#'
#' @author Jay Ver Hoef jverhoef
#' @rdname makeCovMat
#' @export makeCovMat 

makeCovMatInv1D = function(theta,  Nmat, 
	model = 'car', rowStand = TRUE, 
  rhoBound = c(-1,1))
{
  nN = dim(Nmat)[1]
 
  rho = rhoBound[1] + .00005 + expit(theta)*.9999*(rhoBound[2] - rhoBound[1])
  attr(theta,'names') = 'logitRho'
  rs = rep(1, times = nN)
  if(rowStand) rs = apply(Nmat,1,sum)
  if(model == 'car')  Vinv = diag(rs) - rho*Nmat
  if(model == 'sar') Vinv = (diag(nN) - rho*(1/rs)*Nmat) %*%
      (diag(nN) - rho*t((1/rs)*Nmat))

  Vinv
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
#' @param Nmat neighborhood matrix
#' @param model either 'CAR' or 'SAR' to be the model associated with the Nmat argument
#' @param rowStand logical value on whether Nmat should be row-standardized
#' @param rhobound a vector of two elements containing bounds for rho.  This should be determined from the eigenvalues of Nmat prior to running function
#'
#' @return two times the profiled negative log-likelihood
#'
#' @author Jay Ver Hoef jverhoef
#' @rdname m2LL
#' @export m2LL 

m2LL1D = function(theta, X, y, Nmat = NULL, model = 'car', rowStand = TRUE, 
  rhoBound = c(-1,1), MLmeth = 'reml')
{
  indsamp = !is.na(y)
  ysamp = y[indsamp]
	Xsamp = X[indsamp,] 
  nn = length(ysamp)
  nN = length(y)

	WMi = makeCovMatInv1D(theta = theta, Nmat = Nmat, model = model, 
		rowStand = rowStand, rhoBound = rhoBound)
  WMi.oo = WMi[indsamp,indsamp] 
  WMi.uu = WMi[!indsamp,!indsamp]
  WMi.uo = WMi[!indsamp,indsamp]
  WMi.ou = WMi[indsamp,!indsamp]
  Vi.oo = WMi.oo - WMi.ou %*% solve(WMi.uu, WMi.uo)
	XVi = t(Xsamp) %*% Vi.oo
  covbi = XVi %*% Xsamp
  covb = solve(covbi)
  bHat = covb %*% XVi %*% ysamp
  r = ysamp - Xsamp %*% bHat
  n = length(ysamp)
  p = length(X[1,])
  if(MLmeth == 'ml') {
		m2LL = n*log(t(r) %*% Vi.oo %*% r) - 
			as.numeric(determinant(Vi.oo, logarithm = TRUE)$modulus) +
			n*(log(2*pi) + 1 - log(n))
	} else if(MLmeth == 'reml') {
		m2LL = (n-p)*log(t(r) %*% Vi.oo %*% r) - 
			as.numeric(determinant(Vi.oo, logarithm = TRUE)$modulus) +
			as.numeric(determinant(XVi %*% Xsamp, logarithm = TRUE)$modulus) +
			(n - p)*(log(2*pi) + 1 - log((n - p)))
	} else {return('MLmeth argument must be either ml or reml')}
 
	m2LL = as.numeric(m2LL)
}

#-------------------------------------------------------------------------------
#
#           profvar
#
#-------------------------------------------------------------------------------

#' obtain variance parameter after profiling
#'
#' obtain variance parameter after profiling
#'
#' @param theta covariance parameters, with overall variance parameter profiled out.
#' @param X design matrix for fixed effects
#' @param y vector of data for response variable
#' @param Nmat neighborhood matrix
#' @param model either 'CAR' or 'SAR' to be the model associated with the Nmat argument
#' @param rowStand logical value on whether Nmat should be row-standardized
#' @param rhobound a vector of two elements containing bounds for rho.  This should be determined from the eigenvalues of Nmat prior to running function
#'
#' @return estimated variance parameter
#'
#' @author Jay Ver Hoef jverhoef
#' @rdname profvar
#' @export profvar 

profvar = function(theta, X, y, Nmat = NULL, model = 'car', rowStand = TRUE, 
  rhoBound = c(-1,1), MLmeth = 'reml')
{
  indsamp = !is.na(y)
  ysamp = y[indsamp]
	Xsamp = X[indsamp,] 
  nn = length(ysamp)
  nN = length(y)

	WMi = makeCovMatInv1D(theta = theta, Nmat = Nmat, model = model, 
		rowStand = rowStand, rhoBound = rhoBound)
  WMi.oo = WMi[indsamp,indsamp] 
  WMi.uu = WMi[!indsamp,!indsamp]
  WMi.uo = WMi[!indsamp,indsamp]
  WMi.ou = WMi[indsamp,!indsamp]
  Vi.oo = WMi.oo - WMi.ou %*% solve(WMi.uu, WMi.uo)
	XVi = t(Xsamp) %*% Vi.oo
  covbi = XVi %*% Xsamp
  covb = solve(covbi)
  bHat = covb %*% XVi %*% ysamp
  r = ysamp - Xsamp %*% bHat
  n = length(ysamp)
  p = length(X[1,])
  if(MLmeth == 'ml') {
		sig2 = t(r) %*% Vi.oo %*% r / n
	} else if(MLmeth == 'reml') {
		sig2 = t(r) %*% Vi.oo %*% r / (n - p) 
	} else {return('MLmeth argument must be either ml or reml')}
 
	sig2
}

#-------------------------------------------------------------------------------
#
#           Prepare Data using seal data in spmodel
#
#-------------------------------------------------------------------------------

# I want to create a neighborhood structure where all polygons have neighbors
seals_sp = as(seal, 'Spatial')
Nlist = poly2nb(seals_sp, snap = 2000)
# only 1 doesn't have a neighbor, let site 1 be its neighbor
Nlist[[7]] = as.integer(1)

# a few lines of code to create W matrix from Nlist
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
num = lapply(Nlist, function(x) length(x))
num = unlist(num)
Nmat = Neighmat(Nlist, num, length(num))
Nmat1 = pmax(Nmat,t(Nmat))

# now create a data.frame with a categorical variable, just like the
# real data set
DF = seals_sp@data
DF$stock = as.factor(rep(1:4, times = 20)[1:dim(DF)[1]])

X0 = as.matrix(model.matrix( ~ 1, data = DF))
X1 = as.matrix(model.matrix( ~ stock, data = DF))
y = DF$log_trend

#-------------------------------------------------------------------------------
#
#           Fit models using my code versus spmodel
#
#-------------------------------------------------------------------------------

#### My code, pass upper and lower bounds, and response y and design matrix X

# set same options here for my code and spmodel
row_stand = TRUE
cov_mod = 'car'
est_meth = 'ml'

# take care of bounds
evals = eigen(Nmat1)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
if(row_stand == TRUE) {
	LB = -1
	UB = 1
}

# fit with my code
#undebug(m2LL1D)
optOut = optimize(m2LL1D, interval = c(-10,10), X = X0, y = y, model = cov_mod, 
	rowStand = row_stand, Nmat = Nmat1, MLmeth = est_meth, rhoBound = c(LB,UB)
)
theta = optOut$minimum
optOut$objective
LB + expit(theta)*(UB - LB)
profvar(theta, X = X0, y = y, model = cov_mod, 
	rowStand = row_stand, Nmat = Nmat1, MLmeth = est_meth, rhoBound = c(LB,UB))

#### fit with spmodel

spfit = spautor(log_trend ~ 1, data = DF, estmethod = est_meth, 
	control = list(reltol = 1e-7), W = Nmat1, spcov_type = cov_mod, 
	row_st = row_stand)
-2*logLik(spfit)
coef(spfit, type = 'spcov')
