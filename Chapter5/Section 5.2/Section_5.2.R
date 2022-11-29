sec_path = 'Rcode/Chapter5/Section 5.2/'
setwd(paste0(SLEDbook_path,sec_path))

library(spmodel)
library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#              Table 5.1 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Compute some results in Section 4.2
xpts <- matrix(c(1, 1, 2, 2, 5), 5, 1)
ypts <- matrix(c(4, 3, 3, 2, 4), 5, 1)
Dmat = as.matrix(dist(cbind(xpts, ypts)))
Rmat = (1 - 3*Dmat/8 + Dmat^3/128)*(Dmat < 4)
Rmat
Rlower = Rmat
Rhat = (1 - 3*Dmat/6 + Dmat^3/54)*(Dmat < 3)
Rhat
Rlower[lower.tri(Rmat)] = NA

print(
    xtable(Rlower, 
      align = c('l',rep('l', times = length(Rlower[1,]))),
      digits = c(0, rep(3, times = 5))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# Obtain entries of Table 5.1a
X <- matrix(1, 5, 1)
olswts = as.vector(solve(t(X) %*% X) %*% t(X))
glswts = as.vector(solve(t(X) %*% solve(Rmat) %*% X) %*% 
	t(X) %*% solve(Rmat))
eglswts = as.vector(solve(t(X) %*% solve(Rhat) %*% X) %*% 
	t(X) %*% solve(Rhat))
OLS_GLS_wts = data.frame(
	olswts = olswts,
	glswts = glswts,
	eglswts = eglswts
)

print(
    xtable(OLS_GLS_wts, 
      align = c('l',rep('l', times = length(OLS_GLS_wts[1,]))),
      digits = c(0, rep(3, times = 3))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# Obtain entries of Table 5.1b
y <- matrix(c(1,0,3,1,5),5,1)
muhat = olswts%*%y
mutilde = glswts%*%y
mutiltil = eglswts%*%y
mhats = data.frame(
	muhat = muhat,
	mutilde = mutilde,
	mutiltil = mutiltil,
	sigmahat2 = (1/4)*t(y-X%*%muhat)%*%(y-X%*%muhat),
	sigmatilde2 = (1/4)*t(y-X%*%mutilde)%*%solve(Rmat)%*%(y-X%*%mutilde),
	sigmatiltil2 = (1/4)*t(y-X%*%mutiltil)%*%solve(Rhat)%*%(y-X%*%mutiltil)
)

print(
    xtable(mhats, 
      align = c('l',rep('l', times = length(mhats[1,]))),
      digits = c(0, rep(3, times = 6))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# Obtain entries of Table 5.1c
varmhats = data.frame(
	varolsmuhat = t(olswts) %*% olswts,
	varolsmutilde = t(glswts) %*% glswts,
	varglsmuhat = t(olswts) %*% Rmat %*% olswts,
	varglsmutilde = t(glswts) %*% Rmat %*% glswts
)

print(
    xtable(varmhats, 
      align = c('l',rep('l', times = length(varmhats[1,]))),
      digits = c(0, rep(3, times = 4))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# Obtain entries of Table 5.1d
toyreg <- lm(y~X)
Rhalf <- eigen(Rmat)$vectors %*% diag(sqrt(eigen(Rmat)$values)) %*% 
	t(eigen(Rmat)$vectors)
RhalfiX <- solve(Rhalf)%*%X
RiX <- solve(Rmat)%*%X
P <- RiX %*% solve(t(RiX) %*% X) %*% t(RiX)
QQ <- solve(Rmat) %*% (y - X %*% mutilde)
sigmatilde2 = (1/4)*t(y - X %*% mutilde) %*% solve(Rmat) %*% 
	(y - X %*% mutilde)

hsharp <- diag(RhalfiX %*% solve(t(RhalfiX) %*% RhalfiX) %*% t(RhalfiX))
hstar <- diag(P)/diag(solve(Rmat))
tcad <- (QQ/sqrt(diag(solve(Rmat))*(1 - hstar))) %*% (1/sqrt(sigmatilde2))
tstar <- tcad*sqrt(3/(4 - tcad^2))
Dstar <- (tcad*tcad*hstar)/(1-hstar)
covratio <- ((4-tcad*tcad)/3)/(1-hstar)

casediagnostics = data.frame(
	hii = hatvalues(toyreg),
	ti = rstudent(toyreg),
	Di = cooks.distance(toyreg),
	covrati = covratio(toyreg),
	hsharp = hsharp,
	hstar = hstar,
	tstar = tstar,
	Dstar = Dstar,
	covratio = covratio
)

print(
    xtable(casediagnostics, 
      align = c('l',rep('l', times = length(casediagnostics[1,]))),
      digits = c(0, rep(4, times = 9))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                Simulation and Graph of Dish-Shaped Residuals
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

set.seed(1004)
xcoord = rnorm(50)
y = xcoord^2 + rnorm(50,0,0.5)
plot(xcoord, y)
DFsim = data.frame(xcoord = xcoord, ycoord = rep(1, times = length(y)), y = y) 

fitSim = splm(y ~ 1, data = DFsim, xcoord = 'xcoord', ycoord = 'ycoord',
	spcov_type = 'exponential')
summary(fitSim)
SigmaSimFit = coef(fitSim,type = 'spcov')['de']*
	exp(-as.matrix(dist(DFsim$xcoord))/
	coef(fitSim,type = 'spcov')['range']) + 
	coef(fitSim,type = 'spcov')['ie']*diag(50)
FitSimWts = (1/sum(solve(SigmaSimFit))) %*% 
	rep(1, times = 50) %*% solve(SigmaSimFit)
FitSimWts %*% DFsim$y


file_name = 'figures/Moss_dishsim'

pdf(paste0(file_name,'.pdf'), width = 8, height = 12)
	layout(matrix(1:2, nrow = 2))
	par(mar = c(5,5,5,1))
	plot(DFsim$xcoord, DFsim$y, pch = 19, xlab = 'x-coordinate', 
		ylab = 'Response Data', cex.lab = 2.5, cex.axis = 2, cex = 2)
	lines(c(min(DFsim$xcoord),max(DFsim$xcoord)), c(coef(fitSim),coef(fitSim)),
		lwd = 3, lty = 2)
	mtext('A', adj = -.15, cex = 4, padj = -.1)

	par(mar = c(5,5,5,1))
	plot(DFsim$xcoord, as.vector(FitSimWts), pch = 19, cex = 2, 
		xlab = 'x-coordinate', ylab = 'Weights', cex.axis = 2, cex.lab = 2.5)
	mtext('B', adj = -.15, cex = 4, padj = -.1)

	layout(1)
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))
