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
data(sealPolys)

#-------------------------------------------------------------------------------
#                    Create Neighborhood Matrices
#-------------------------------------------------------------------------------

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
Nmat2 = (Nmat1 %*% Nmat1 > 0 | Nmat1 > 0)*1
diag(Nmat2) = 0
Nmat4 = (Nmat2 %*% Nmat2 > 0 | Nmat2 > 0)*1
diag(Nmat4) = 0

AllPolyCentroids = data.frame(x = coordinates(sealPolys)[,1], 
    y = coordinates(sealPolys)[,2], 
    stockid = as.factor(as.character(sealPolys@data$stockid)),
    polyid = as.factor(as.character(sealPolys@data$polyid)))
distMat = as.matrix(dist(AllPolyCentroids[,c('x','y')]))/1000
distMat1 = distMat*Nmat
rownames(distMat1) = attr(Nlist,'polyid')
colnames(distMat1) = attr(Nlist,'polyid')
distMat2 = distMat*Nmat2
rownames(distMat2) = attr(Nlist,'polyid')
colnames(distMat2) = attr(Nlist,'polyid')
distMat4 = distMat*Nmat4
rownames(distMat4) = attr(Nlist,'polyid')
colnames(distMat4) = attr(Nlist,'polyid')

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

makeCovMat = function(theta, indComp = TRUE, Nmat = NULL, 
  distMat = NULL, indSamp, model = 'CAR', rowStand = TRUE, 
  rhoBound = c(-1,1))
{
  nn = sum(indSamp)
  nN = length(indSamp)
  V = matrix(1, nrow = nN, ncol = nN )
  diag(V) = 0
  itheta = 0
  if(!is.null(Nmat)) {
    V = as(Nmat, 'sparseMatrix')
  }
  if(!is.null(distMat)) {
    itheta = itheta + 1
    attr(theta,'names')[itheta] = 'logDist'
    V = V * exp(-distMat/exp(theta[itheta]))
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
  if(indComp & any(c(!is.null(Nmat), !is.null(distMat)))) {
    itheta = itheta + 1
    relEps = exp(theta[itheta])
    attr(theta,'names')[itheta] = 'relEps'
    V = V -  V %*% solve(V + diag(rep(relEps,times = nN)),V)
  }
  if(indComp & !any(c(!is.null(Nmat), !is.null(distMat)))) {
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
  distMat = NULL, indSamp, model = 'CAR', rowStand = TRUE, 
  rhoBound = c(-1,1), MLmeth = 'REMLE')
{
  if(any(abs(theta) > 10.1)) return(1e+32)

  ntheta = 0
  nn = length(y)
  nN = length(indSamp)
	V = makeCovMat(theta = theta, indComp = indComp, Nmat = Nmat, distMat = distMat, 
		indSamp = indSamp, model = model, rowStand = rowStand, rhoBound = rhoBound)
  WMi = V
  WMi.oo = WMi[indSamp,indSamp] 
  WMi.uu = WMi[!indSamp,!indSamp]
  WMi.uo = WMi[!indSamp,indSamp]
  WMi.ou = WMi[indSamp,!indSamp]
  Vi.oo = WMi.oo - WMi.ou %*% solve(WMi.uu, WMi.uo)
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

DF = sealPolys@data
indSamp = !is.na(DF$Estimate)
X0 = as.matrix(model.matrix(Estimate ~ 1, data = DF))
X1 = as.matrix(model.matrix(Estimate ~ I(as.factor(stockid)), data = DF))
y = DF$Estimate[!is.na(DF$Estimate)]

theta = c(2, 0)

ntheta = 2

if(ntheta == 1) {
	# undebug(m2LL)
	# undebug(makeCovMat)
  optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
		Nmat = Nmat4, indComp = FALSE, indSamp = indSamp)
  theta = optOut$minimum
  m2LLargmin = optOut$objective
  } else {
	# undebug(m2LL)
	# undebug(makeCovMat)
  optOut = optim(rep(0, times = ntheta), m2LL, X = X1, y = y, 
		indSamp = indSamp, Nmat = Nmat4, indComp = FALSE)
  theta = optOut$par
  m2LLargmin = optOut$value
}
X = X1
WMi = makeCovMat(theta, indSamp = indSamp, Nmat = Nmat4, indComp = FALSE)
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
