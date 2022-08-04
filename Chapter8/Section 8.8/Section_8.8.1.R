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
library(colorspace)
library(spdep)
library(Matrix)
library(xtable)
library(spmodel)

# load data for graphics and analysis
data(sealPolys)
seals_sf = st_as_sf(sealPolys)
seals_sf$stockname = as.factor(as.character(seals_sf$stockname))

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
#' @param indComp an additive independent component to the model.  Default is TRUE.  
#' @param Nmat neighborhood matrix
#' @param distMat distance matrix
#' @param indSamp indicator vector for wheter location was sampled. Zero, or FALSE, indicates missing value
#' @param model either 'CAR' or 'SAR' to be the model associated with the Nmat argument
#' @param logical value on whether Nmat should be row-standardized
#' @param rhobound a vector of two elements containing bounds for rho.  This should be determined from the eigenvalues of Nmat prior to running function
#'
#' @return A covariance matrix
#'
#' @author Jay Ver Hoef jverhoef
#' @rdname makeCovMat
#' @export makeCovMat 

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
      (diag(nN) - rho*Matrix:::t((1/rs)*V))
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
  list(V = V, theta = theta)
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
	Vlist = makeCovMat(theta = theta, indComp = indComp, Nmat = Nmat, distMat = distMat, 
		indSamp = indSamp, model = model, rowStand = rowStand, rhoBound = rhoBound)
  WMi = Vlist$V
  theta = Vlist$theta
  WMi.oo = WMi[indSamp,indSamp] 
  WMi.uu = WMi[!indSamp,!indSamp]
  WMi.uo = WMi[!indSamp,indSamp]
  WMi.ou = WMi[indSamp,!indSamp]
  Vi.oo = WMi.oo - WMi.ou %*% Matrix:::solve(WMi.uu, WMi.uo)
  XVi = t(X) %*% as.matrix(Vi.oo)
  covbi = XVi %*% X
  covb = solve(covbi)
  bHat = covb %*% XVi %*% y
  r = y - X %*% bHat
  n = length(y)
  p = length(X[1,])
  if(MLmeth == 'MLE') {
  m2LL = n*log(t(r) %*% Vi.oo %*% r) - 
    as.numeric(Matrix:::determinant(Vi.oo, logarithm = TRUE)$modulus) +
    n*(log(2*pi) + 1 - log(n))
	} else if(MLmeth == 'REMLE') {
  m2LL = (n-p)*log(t(r) %*% Vi.oo %*% r) - 
    as.numeric(Matrix:::determinant(Vi.oo, logarithm = TRUE)$modulus) +
    as.numeric(determinant(XVi %*% X, logarithm = TRUE)$modulus) +
    (n - p)*(log(2*pi) + 1 - log((n - p)))
	} else {return('MLmeth argument must be either MLE or REMLE')}
 
	m2LL = as.numeric(m2LL)
	attr(m2LL,'covParms') = theta
	m2LL
}

DF = sealPolys@data
indSamp = !is.na(DF$Estimate)
X0 = as.matrix(model.matrix(Estimate ~ 1, data = DF))
X1 = as.matrix(model.matrix(Estimate ~ I(as.factor(stockid)), data = DF))
y = DF$Estimate[!is.na(DF$Estimate)]

################################################################################
#-------------------------------------------------------------------------------
#                  CAR and SAR Likelihoods
#-------------------------------------------------------------------------------
################################################################################

#reproduce parts of Figure 5 in Ver Hoef et al. 2018
#investigate CAR versus SAR, and Row-standardized versus unstandardized
evals = eigen(Nmat1)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
#undebug(m2LL)
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat1, indComp = FALSE, MLmeth = 'MLE',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_XC1 = optOut$objective

spfit_mC1R = spautor(Estimate ~ 1, data = seals_sf, spcov_type = 'car',
	W = Nmat1, row_st = TRUE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_mC1R)

spfit_XU = splm(Estimate ~ stockname, data = seals_sf, spcov_type = 'none',
	control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XU)

spfit_XC1 = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'car',
	W = Nmat1, row_st = FALSE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XC1)

spfit_XC1R = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'car',
	W = Nmat1, row_st = TRUE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XC1R)

spfit_XS1 = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'sar',
	W = Nmat1, row_st = FALSE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XS1)

spfit_XS1R = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'sar',
	W = Nmat1, row_st = TRUE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XS1R)

spfit_XC2 = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'car',
	W = Nmat2, row_st = FALSE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XC2)

spfit_XC2R = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'car',
	W = Nmat2, row_st = TRUE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XC2R)

spfit_XS2 = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'sar',
	W = Nmat2, row_st = FALSE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XS2)

spfit_XS2R = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'sar',
	W = Nmat2, row_st = TRUE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XS2R)

spfit_XC4 = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'car',
	W = Nmat4, row_st = FALSE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XC4)

spfit_XC4R = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'car',
	W = Nmat4, row_st = TRUE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XC4R)

spfit_XS4 = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'sar',
	W = Nmat4, row_st = FALSE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XS4)

spfit_XS4R = spautor(Estimate ~ stockname, data = seals_sf, spcov_type = 'sar',
	W = Nmat4, row_st = TRUE, control = list(reltol = 1e-7), estmethod = 'ml')
2*logLik(spfit_XS4R)


evals = eigen(Nmat2)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat2, indComp = FALSE, MLmeth = 'MLE',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_C2U = optOut$objective

evals = eigen(Nmat4)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'MLE',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_C4U = optOut$objective


evals = eigen(Nmat1)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
#undebug(makeCovMat)
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat1, indComp = FALSE, MLmeth = 'MLE', model = 'SAR',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_S1U = optOut$objective


evals = eigen(Nmat2)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat2, indComp = FALSE, MLmeth = 'MLE', model = 'SAR',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_S2U = optOut$objective

evals = eigen(Nmat4)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'MLE', model = 'SAR',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_S4U = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat1, indComp = FALSE, MLmeth = 'MLE')
theta = optOut$minimum
m2LLargmin_C1R = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat2, indComp = FALSE, MLmeth = 'MLE')
theta = optOut$minimum
m2LLargmin_C2R = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'MLE')
theta = optOut$minimum
m2LLargmin_C4R = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, model = 'SAR',
	indSamp = indSamp, Nmat = Nmat1, indComp = FALSE, MLmeth = 'MLE')
theta = optOut$minimum
m2LLargmin_S1R = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, model = 'SAR',
	indSamp = indSamp, Nmat = Nmat2, indComp = FALSE, MLmeth = 'MLE')
theta = optOut$par
m2LLargmin_S2R = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, model = 'SAR',
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'MLE')
theta = optOut$par
m2LLargmin_S4R = optOut$objective

file_name = 'seal_2LL'
pdf(paste0(file_name,'.pdf'), width = 8.5, height = 8.5)
	m2LLargmin = c(m2LLargmin_C1U, m2LLargmin_C2U, m2LLargmin_C4U, m2LLargmin_S1U,
		m2LLargmin_S2U, m2LLargmin_S4U, m2LLargmin_C1R, m2LLargmin_C2R,
		m2LLargmin_C4R, m2LLargmin_S1R,  m2LLargmin_S2R,  m2LLargmin_S4R)
	labs = c('C1U', 'C2U', 'C4U', 'S1U', 'S2U', 'S4U', 	
		'C1R', 'C2R', 'C4R', 'S1R', 'S2R', 'S4R')
	par(mar = c(5,5,1,1))
	plot(1:12, -m2LLargmin, pch = 19, cex = 3, ylab = '2(log-likelihood)', 
		cex.axis = 1.5, cex.lab = 2, xaxt = 'n', xlab = '')
	axis(1, at = 1:12, labels = labs, las = 2, cex.axis = 1.5)
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

# Because the fixed effects do not change, we can also use REMLE

evals = eigen(Nmat1)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat1, indComp = FALSE, MLmeth = 'REMLE',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_REML_C1U = optOut$objective

evals = eigen(Nmat2)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat2, indComp = FALSE, MLmeth = 'REMLE',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_REML_C2U = optOut$objective

evals = eigen(Nmat4)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'REMLE',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_REML_C4U = optOut$objective


evals = eigen(Nmat1)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat1, indComp = FALSE, MLmeth = 'REMLE', model = 'SAR',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_REML_S1U = optOut$objective


evals = eigen(Nmat2)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat2, indComp = FALSE, MLmeth = 'REMLE', model = 'SAR',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_REML_S2U = optOut$objective

evals = eigen(Nmat4)$values
minevals = min(evals)
maxevals = max(evals)
LB = 1/minevals
UB = 1/maxevals
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, rowStand = FALSE,
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'REMLE', model = 'SAR',
	rhoBound = c(LB,UB))
theta = optOut$minimum
m2LLargmin_REML_S4U = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat1, indComp = FALSE, MLmeth = 'REMLE')
theta = optOut$minimum
m2LLargmin_REML_C1R = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat2, indComp = FALSE, MLmeth = 'REMLE')
theta = optOut$minimum
m2LLargmin_REML_C2R = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'REMLE')
theta = optOut$minimum
m2LLargmin_REML_C4R = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, model = 'SAR',
	indSamp = indSamp, Nmat = Nmat1, indComp = FALSE, MLmeth = 'REMLE')
theta = optOut$minimum
m2LLargmin_REML_S1R = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, model = 'SAR',
	indSamp = indSamp, Nmat = Nmat2, indComp = FALSE, MLmeth = 'REMLE')
theta = optOut$par
m2LLargmin_REML_S2R = optOut$objective

optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, model = 'SAR',
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'REMLE')
theta = optOut$par
m2LLargmin_REML_S4R = optOut$objective

file_name = 'seal_2LL_REMLE'
pdf(paste0(file_name,'.pdf'), width = 8.5, height = 8.5)
	m2LLargmin_REML = c(m2LLargmin_REML_C1U, m2LLargmin_REML_C2U, m2LLargmin_REML_C4U, m2LLargmin_REML_S1U,
		m2LLargmin_REML_S2U, m2LLargmin_REML_S4U, m2LLargmin_REML_C1R, m2LLargmin_REML_C2R,
		m2LLargmin_REML_C4R, m2LLargmin_REML_S1R,  m2LLargmin_REML_S2R,  m2LLargmin_REML_S4R)
	labs = c('C1U', 'C2U', 'C4U', 'S1U', 'S2U', 'S4U', 	
		'C1R', 'C2R', 'C4R', 'S1R', 'S2R', 'S4R')
	par(mar = c(5,5,1,1))
	plot(1:12, -m2LLargmin_REML, pch = 19, cex = 3, ylab = '2(log-likelihood)', 
		cex.axis = 1.5, cex.lab = 2, xaxt = 'n', xlab = '')
	axis(1, at = 1:12, labels = labs, las = 2, cex.axis = 1.5)
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
#          Profile Likelihood for Autocorrelation Parameter
#-------------------------------------------------------------------------------
################################################################################

thetatry = (-10:70)/20
rhovals = -1 + 2*expit(thetatry)
LL2_C4R = NULL
for(i in thetatry)
	LL2_C4R = c(LL2_C4R, 
		-m2LL(i, X = X1, y = y, indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, 
			MLmeth = 'MLE'))

file_name = 'seal_rhoprofile'
pdf(paste0(file_name,'.pdf'), width = 6, height = 6)

	plot(rhovals, LL2_C4R, type = 'l', lwd = 3, ylab = '2 x loglikelihood',
		xlab = 'rho values')
	LB = -m2LLargmin_C4R - qchisq(0.95, df = 1)
	lines(c(rhovals[1], rhovals[81]), c(LB, LB),
		lty = 2, lwd = 3)
	
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

# use profile likelihood to get confidence interval on autocorrelation parameter	

minrhoindx = min(which(LL2_C4R > -m2LLargmin_C4R - qchisq(0.95, df = 1)))
maxrhoindx = max(which(LL2_C4R > -m2LLargmin_C4R - qchisq(0.95, df = 1)))

# linear interpolation for lower bound if confidence interval
mean(rhovals[minrhoindx],rhovals[minrhoindx-1])
mean(rhovals[maxrhoindx],rhovals[maxrhoindx+1])

################################################################################
#-------------------------------------------------------------------------------
#                  Estimating Fixed Effects
#-------------------------------------------------------------------------------
################################################################################

#use independence model
lm_summ_out = summary(lm(Estimate ~ I(as.factor(stockid)), data = DF))
FixEff_LM = as.data.frame(lm_summ_out$coefficients)
rownames(FixEff_LM) = 1:5
FixEff_LM = cbind(data.frame(Effect = c('Intercept', 'Stock 9', 'Stock 10', 'Stock 11', 'Stock 12')),
	FixEff_LM)

#use C4R MLE model
X = X1
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'MLE')
theta = optOut$minimum

WMi = as.matrix(makeCovMat(theta, indSamp = indSamp, Nmat = Nmat4, 
	indComp = FALSE)$V)
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
Pval = 2*(1-pt(abs(bHat/bHat_se), df = n - p))
FixEff_MLE = cbind(bHat,
	bHat_se,
	bHat/bHat_se,
	Pval)
rownames(FixEff_MLE) = 1:5
FixEff_MLE = cbind(data.frame(Effect = c('Intercept', 'Stock 9', 'Stock 10', 'Stock 11', 'Stock 12')),
	FixEff_MLE)
colnames(FixEff_MLE) = colnames(FixEff_LM)

#use C4R REMLE model
X = X1
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'REMLE')
theta = optOut$minimum

WMi = as.matrix(makeCovMat(theta, indSamp = indSamp, Nmat = Nmat4, 
	indComp = FALSE)$V)
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
Pval = 2*(1-pt(abs(bHat/bHat_se), df = n - p))
FixEff_REMLE = cbind(bHat,
	bHat_se,
	bHat/bHat_se,
	Pval)
rownames(FixEff_REMLE) = 1:5
FixEff_REMLE = cbind(data.frame(Effect = c('Intercept', 'Stock 9', 'Stock 10', 'Stock 11', 'Stock 12')),
	FixEff_REMLE)
colnames(FixEff_REMLE) = colnames(FixEff_LM)

FixEff = rbind(FixEff_LM, FixEff_MLE, FixEff_REMLE)
print(
    xtable(FixEff, 
      align = c('l',rep('l', times = length(FixEff[1,]))),
      digits = c(0,0,rep(3, times = 3),5),
      caption = 'Fitted fixed effects',
      label = 'tab:SealsFixEff'
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# Compare Stock 9 to Stock 11 using REML estimates
cont = matrix(c(0,-1, 0, 1, 0), ncol = 1)
t(cont) %*% bHat
# standard error
sqrt(t(cont) %*% (sigma*covb) %*% cont)
# t-value 
(t(cont) %*% bHat)/sqrt(t(cont) %*% (sigma*covb) %*% cont)
# Prob of t-value given null hypothesis of equality
2*(1-pt(abs((t(cont) %*% bHat)/sqrt(t(cont) %*% (sigma*covb) %*% cont)), 
	df = n - p))


################################################################################
#-------------------------------------------------------------------------------
#                  Prediction
#-------------------------------------------------------------------------------
################################################################################

X = X1
optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat4, indComp = FALSE, MLmeth = 'MLE')
theta = optOut$minimum
theta
-1 + expit(theta)*2
      
WMi = as.matrix(makeCovMat(theta, indSamp = indSamp, Nmat = Nmat4, 
	indComp = FALSE)$V)
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
Vi.oo = Vi.oo/sigma
covb = sigma*covb

Xall = as.matrix(model.matrix( ~ I(as.factor(stockid)), data = DF))
Xp = Xall[!indSamp,]
Vall = sigma*solve(WMi)
Vpred = Vall[indSamp,!indSamp]
ViVpred = Vi.oo %*% Vpred
ViX = Vi.oo %*% X1
preds <- matrix(NA, nrow = sum(!indSamp), ncol = 2)
preds[,1] <- apply(as.vector(Vi.oo %*% y) * Vpred, 2, sum) +
	Xp %*% bHat - t(Vpred) %*% ((Vi.oo %*% X1) %*% bHat)	
preds[,2] <- sqrt(diag(Vall[!indSamp, !indSamp]) - 
	apply(ViVpred * Vpred, 2, sum) +
	apply((covb %*% t(Xp)) * t(Xp), 2, sum) -
	2*apply((covb %*% t(Xp)) * (t(X1) %*% ViVpred), 2, sum) +
	apply(((covb %*% t(ViX)) %*% Vpred) * (t(X) %*% ViVpred), 2, sum))
# program it from Schabenberger/Gotway, pg. 243
pse  = sqrt(diag(Vall[!indSamp, !indSamp] - t(Vpred) %*% ViVpred + 
	(Xp - t(Vpred) %*% ViX) %*% solve(t(X1) %*% ViX) %*%
	t(Xp - t(Vpred) %*% ViX)))
cbind(preds[,2],pse)


plot(sealPolys)
Pts = coordinates(sealPolys)
Pts = Pts[!indSamp,]
points(Pts, pch = 19, col = 'red')

file_name = 'seal_preds'
pdf(paste0(file_name,'.pdf'), width = 17, height = 8.5)

	layout(matrix(c(1,2), nrow = 1))
	source('addBreakColorLegend.R')
	cip = classIntervals(preds[,1], n = 6, style = 'fisher')
	palp = viridis(6)
	cip_colors = findColours(cip, palp)
	old.par = par(mar = c(0,0,5,0))
	plot(sealPolys)
	points(Pts, pch = 19, col = cip_colors, cex = 2)
	addBreakColorLegend(xleft = 1330000, ybottom = 986649, xright = 1390000, ytop = 1201343,
		breaks = cip$brks, colors = palp, cex = 2, printFormat = "4.3")
	text(920000, 1190000, 'A', cex = 6)

	cip = classIntervals(preds[,2], n = 6, style = 'fisher')
	palp = cividis(6)
	cip_colors = findColours(cip, palp)
	old.par = par(mar = c(0,0,5,0))
	plot(sealPolys)
	points(Pts, pch = 19, col = cip_colors, cex = 2)
	addBreakColorLegend(xleft = 1330000, ybottom = 986649, xright = 1390000, ytop = 1201343,
		breaks = cip$brks, colors = palp, cex = 2, printFormat = "4.3")
	text(920000, 1190000, 'B', cex = 6)
	
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))
