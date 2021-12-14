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
library(spdep)

# load data for graphics and analysis
data(caribouDF)

#-------------------------------------------------------------------------------
#                    Create Neighborhood Matrices
#-------------------------------------------------------------------------------

nTot = length(sealPolys)
nObs = sum(!is.na(sealPolys$Estimate))
nMiss = nTot - nObs

# get Euclidean distance between centroids of plots
Distmat = as.matrix(dist(caribouDF[,c('x','y')]))
# create first-order neighbor matrix (rook's move) from distances
Nmat1 = (Distmat < 1.01)*1
diag(Nmat1) = 0
Nmat2 = (Nmat1 %*% Nmat1 > 0 | Nmat1 > 0)*1
diag(Nmat2) = 0

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
