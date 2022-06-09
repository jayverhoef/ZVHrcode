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
library(spdep)
library(spmodel)

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

		
D1 = sealPolys@data
# need to eliminate un-used factor levels
D1$stockname = as.factor(as.character(D1$stockname)

#undebug(spautor)
#undebug(spmodel:::cov_estimate_gloglik)
#undebug(spmodel:::norand_ar)
#undebug(spmodel:::cov_initial_search)
# try the function with just data.frame and weights matrix
spautor(Estimate ~ stockname, data = D1, W = Nmat1, 'car')
# the above statement breaks at this code
#  if (!requireNamespace("sf", quietly = TRUE)) {
#    coords <- NULL
#  } else {
#    coords <- sf::st_geometry(data)
#  }

# try the function as sp object with default weights
spautor(Estimate ~ stockname, data = sealPolys, 'car')
# breaks because of un-used factors, but lm() works
#undebug(lm)
lm(Estimate ~ stockname, data = D1)

# get rid of un-used factors
sealPolys1 = sealPolys
sealPolys1@data = D1
spautor_out1 = spautor(Estimate ~ stockname, data = sealPolys1, 'car')
summary(spautor_out1)
# Not sure why there is an extra parameter estimated here?  
# Ah, because of unconnected sites! Right, good.

# Now try it with my weights, where I forced all locations to be connected
spautor_out2 = spautor(Estimate ~ stockname, W = Nmat1, data = sealPolys1, 'car')
summary(spautor_out2)
# Works, only 3 parameters.  My weights are not row-standardized, 
# so let's change to row-standardized, to see if anything changes.
Rmat1 = Nmat1/apply(Nmat1,1,sum)
spautor_out3 = spautor(Estimate ~ stockname, W = Rmat1, data = sealPolys1, 'car')
# Error in asMethod(object) : not a positive definite matrix

# I want a model without the independence parameter, so I need to
# force it to be zero.
print(spcov_initial("car", de = .2, ie = 0, range = 0.9, known = c("ie")))
spautor_out4 = spautor(Estimate ~ stockname, data = sealPolys1, 'car', W = Nmat1, 
 spcov_initial = spcov_initial("car", de = .2, ie = 0, range = 0.9, known = c("ie")))
summary(spautor_out4)

# I want a model without the independence parameter, but using
# unconnected sites. This is probably the best model.
print(spcov_initial("car", de = .2, ie = 0, range = 0.9, known = c("ie")))
spautor_out5 = spautor(Estimate ~ stockname, data = sealPolys1, 'car', 
 spcov_initial = spcov_initial("car", de = .2, ie = 0, range = 0.9, known = c("ie")))
summary(spautor_out5)

# compare to lm().  Not a lot of autocorrelation in the model.
summary(lm(Estimate ~ stockname, data = D1))

# Try model 4, but using ML rather than REML.
print(spcov_initial("car", de = .2, ie = 0, range = 0.9, known = c("ie")))
spautor_out4a = spautor(Estimate ~ stockname, data = sealPolys1, 'car', W = Nmat1, 
 spcov_initial = spcov_initial("car", de = .2, ie = 0, range = 0.9, known = c("ie")),
 estmethod = 'ml')
summary(spautor_out4a)

# Try model 4, but using ML and weights from Nmat4, to compare to 
# Ver Hoef et al, 2018, Figure 6
#undebug(spautor)
spautor_out4b = spautor(Estimate ~ I(as.factor(stockid)), 
  data = sealPolys1, 'car', W = Nmat4, 
  spcov_initial = spcov_initial("car", de = .2, ie = 0, range = 0.9, known = c("ie")),
  estmethod = 'ml', control = list(reltol = 1e-8))
summary(spautor_out4b)
# looks pretty good.  Range parameter seems to match Ver Hoef et al, 2018,
# as well as parameter estimates in Table 3.
# Check out the loglikelihood.
spautor_out4b$optim$value
# looks like it matches Ver Hoef et al, 2018, Figure 6
# I used my own code and our results match almost exactly.

# try prediction
# undebug(spmodel:::predict.spautor)
predict(spautor_out4b, se.fit = TRUE, interval = 'prediction')
# upper and lower bounds don't see right.  See below.
# I think that we want se.prediction too?
predict(spautor_out4b, se.fit = TRUE, interval = 'confidence')
# This fails, with error
# Error in vapply(newdata_model_split, function(x) crossprod(x, vcov(object) %*%  : 
#  values must be type 'double',
# but FUN(X[[1]]) result is type 'S4'

# Here is what I think it should be for prediction
level = .95
sigma = spautor_out4b$coefficients$spcov['de']
range = spautor_out4b$coefficients$spcov['range']
W = spautor_out4b$W
indSamp = spautor_out4b$observed_index
indMiss = spautor_out4b$missing_index
# get the formula
form = spautor_out4b$formula
# get the full data
fulldata = as.data.frame(spautor_out4b$fulldata)
# get just the right-hand side of the formula to create model matrix for all data
Xall = model.matrix(form[-2], fulldata)
# get design matrices for observed and missing data
Xo = Xall[indSamp,]
Xp = Xall[indMiss,]
# get the observed data
y = spautor_out4b$model$y
# inverse covariance matrix from CAR model
Vi = (diag(apply(W,1,sum)) - range*W)/sigma
# get partitions of inverse covariance matrix
Vi.oo = Vi[indSamp,indSamp] 
Vi.uu = Vi[indMiss,indMiss]
Vi.uo = Vi[indMiss,indSamp]
Vi.ou = Vi[indSamp,indMiss]
# inverse covariance matrix for observed data
Vi.oo = Vi.oo - Vi.ou %*% solve(Vi.uu, Vi.uo)
# covariance matrix
V = Matrix:::solve(Vi)
# covariance between observed and predicted sites
Vpred = V[indSamp, indMiss]
ViVpred = as.matrix(Vi.oo %*% Vpred)
# parts for estimating fixed effects
XVi = as.matrix(t(Xo) %*% Vi.oo)
ViX = t(XVi)
covbi = XVi %*% X
covb = solve(covbi)
# estimated regression coefficients
bHat = covb %*% XVi %*% y
bHat
# raw residuals
r = as.matrix(y - X %*% bHat)
preds <- matrix(NA, nrow = length(indMiss), ncol = 4)
preds[,1] <- as.vector(apply(as.vector(Vi.oo %*% y) * Vpred, 2, sum) +
	Xp %*% bHat - t(Vpred) %*% ((Vi.oo %*% Xo) %*% bHat))	
# program prediction standard errors from Schabenberger/Gotway, pg. 243
preds[,2] =  sqrt(diag(V[indMiss, indMiss] - t(Vpred) %*% ViVpred + 
	(Xp - t(Vpred) %*% ViX) %*% solve(t(Xo) %*% ViX) %*%
	t(Xp - t(Vpred) %*% ViX)))
# a computing formula that may be faster for prediction standard errors
predcompalt <- sqrt(diag(V[indMiss, indMiss]) - 
	apply(ViVpred * Vpred, 2, sum) +
	apply((covb %*% t(Xp)) * t(Xp), 2, sum) -
	2*apply((covb %*% t(Xp)) * (t(X1) %*% ViVpred), 2, sum) +
	apply(((covb %*% t(ViX)) %*% Vpred) * (t(X) %*% ViVpred), 2, sum))
tstar <- qnorm(1 - (1 - level) / 2)
lwr <- preds[,1] - tstar * preds[,2]
upr <- preds[,1] + tstar * preds[,2]
preds[,3] = lwr
preds[,4] = upr
