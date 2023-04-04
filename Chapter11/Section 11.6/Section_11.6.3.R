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


ell = c(0, 0, -0.5, 1, 0, 0)
tval = t(ell) %*% coef(spfitEXP)/
	sqrt(MSE*t(ell) %*% solve(t(X3star) %*% X3star) %*% ell)
tval^2
1 - pchisq(tval^2, df = 1)

ell = c(1, 0, 0, 1, 0, 0)

# Results in last paragraph of Section 11.6.3
drvlmEXPest = (ell %*% coef(spfitEXP) + mean(resid(spfitEXP)))
drvlmEXPse = sqrt(t(ell - apply(model.matrix(z ~ water*tarp, 
	data = caribouDF),2,mean)) %*% vcov(spfitEXP) %*% 
	(ell - apply(model.matrix(z ~ water*tarp, data = caribouDF),2,mean)))
drvlmEXPest
drvlmEXPse

lmEXPest = ell %*% coef(spfitEXP)
lmEXPse = sqrt(t(ell) %*% vcov(spfitEXP) %*% ell)
lmEXPest
lmEXPse
##

