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

# fit exponential geostatistical model
spfitEXP = splm(z ~ water*tarp, data = caribouDF, xcoord = 'x', ycoord = 'y', 
	spcov_type = 'exponential', estmethod = 'reml')
smryEXP = summary(spfitEXP)

# model matrix for all factors
Xall = model.matrix(z ~ water*tarp, data = caribouDF)
# model matrix for mean only
X0 = Xall[,1]
# model matrix for mean and first factor
X1 = Xall[,1:2]
# model matrix for mean, first and second factors
X2 = Xall[,1:4]
# model matrix for mean and factors 1-3 (same as Xall, but renamed)
X3 = Xall[,1:6]

# total variance
totvar = coef(spfitEXP, type = 'spcov')['de'] + 
	coef(spfitEXP, type = 'spcov')['ie']
# proportion of variance that is nugget effect
propnug = coef(spfitEXP, type = 'spcov')['ie']/totvar

# R-tilda (correlation matrix, without overall variance
Sigma = (1-propnug)*
	exp(-as.matrix(dist(caribouDF[,c('x','y')]))/
	coef(spfitEXP, type = 'spcov')['range']) + 
	propnug*diag(30)
# inverse square-root of R-tilda
Sigma_mhalf = eigen(Sigma)$vectors %*% diag(1/sqrt(eigen(Sigma)$values)) %*%
	t(eigen(Sigma)$vectors)
# y-star and X-star by premultiplying by inverse square-root of R-tilda
ystar = Sigma_mhalf %*% caribouDF$z
X0star = Sigma_mhalf %*% X0
X1star = Sigma_mhalf %*% X1
X2star = Sigma_mhalf %*% X2
X3star = Sigma_mhalf %*% X3

# orthogonal projection matrices for X-star matrices
Px0 = X0star %*% solve(t(X0star) %*% X0star, t(X0star))
Px1 = X1star %*% solve(t(X1star) %*% X1star, t(X1star))
Px2 = X2star %*% solve(t(X2star) %*% X2star, t(X2star))
Px3 = X3star %*% solve(t(X3star) %*% X3star, t(X3star))

# ANOVA, as in Table 4.2 (can use standard ANOVA table because we 
# "whitened" the model by premultiplying by inverse square-root of R-tilda)
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

################################################################################
#                   Table 11.6
################################################################################

# put together the ANOVA table from the parts given above
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

# R-squared from the ANOVA table
R2 = sum(SSvec[1:3])/SSvec[5]
R2
# or use spmodel
pseudoR2(spfitEXP)

# R-squared assuming independence (compute it in approximately the same way)
lm_anova = anova(lm(z ~ water*tarp, data = caribouDF))
sum(lm_anova[1:3,2])/(29*var(caribouDF$z))
# or use built-in from lm()
summary(lm(z ~ water*tarp, data = caribouDF))$r.squared

# ANOVA for water effects (two levels)
# using a t-value on a 1-degree-of-freedom test
# create the constrast
ell = c(0, 1, 0, 0, 0, 0)
# compute the t-value
tval = t(ell) %*% coef(spfitEXP)/
	sqrt(MSE*t(ell) %*% solve(t(X3star) %*% X3star) %*% ell)
tval
# compute the probability
2*(1 - pt(tval, df = 24))

# compute the probability assuming a Chi-square (can generalize to more than
# a single degree of freedom)
1 - pchisq(tval^2, df = 1)
# or, we can use spmodel
anova(spfitEXP)[2,]

# ANOVA for tarp effects (three levels)
ell = cbind(c(0, 0, 1, 0, 0, 0), c(0, 0, 0, 1, 0, 0))
# compute the C-value
Cval = t(t(ell) %*% coef(spfitEXP)) %*%
	solve(t(ell) %*% solve(t(X3star) %*% X3star) %*% ell) %*% 
	t(ell) %*% coef(spfitEXP)/MSE
Cval
# compute the probability
1 - pchisq(Cval, df = ncol(ell))
# or, we can use spmodel
anova(spfitEXP)[3,]

# ANOVA for interaction effects (three levels)
ell = cbind(c(0, 0, 0, 0, 1, 0), c(0, 0, 0, 0, 0, 1))
# compute the C-value
Cval = t(t(ell) %*% coef(spfitEXP)) %*%
	solve(t(ell) %*% solve(t(X3star) %*% X3star) %*% ell) %*% 
	t(ell) %*% coef(spfitEXP)/MSE
Cval
# compute the probability
1 - pchisq(Cval, df = ncol(ell))
# or, we can use spmodel
anova(spfitEXP)[3,]

