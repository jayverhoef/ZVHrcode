sec_path = 'Rcode/Chapter7/Section 7.2/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)
library(viridis)
library(classInt)
source('addBreakColorLegend.R')

################################################################################
#-------------------------------------------------------------------------------
#              Topologically-induced Heterogeneity
#-------------------------------------------------------------------------------
################################################################################

# Computation of marginal variances and correlations for the SAR and CAR models 
# on the 7-site example of Section 7.2 with row-standardized weights

# create W-matrix
W0 <- matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,byrow=T)
id3 <- diag(3)
ze3 <- matrix(0,3,3)
W1 <- cbind(W0,id3,ze3)
W2 <- cbind(id3,W0,id3)
W3 <- cbind(ze3,id3,W0)
W <- rbind(W1,W2,W3)
W <- W[-1,-1]
W <- W[-2,-2]
# diagonal matrix of ones
K <- diag(7)
# minimum and maximum eigenvalues for W
min(eigen(W)$values)
max(eigen(W)$values)
# parameter space for rho
1/min(eigen(W)$values)
1/max(eigen(W)$values)

# SAR covariance matrix
SigmaSAR <- solve(diag(7)-0.3*W) %*% K %*% solve(diag(7)-0.3*t(W))
# CAR covariance matrix
SigmaCAR <- solve(diag(7)-0.3*W) %*% K
# SAR correlation matrix
corrSAR <- solve(diag(sqrt(diag(SigmaSAR)))) %*% 
	SigmaSAR%*%solve(diag(sqrt(diag(SigmaSAR))))
# CAR correlation matrix
corrCAR <- solve(diag(sqrt(diag(SigmaCAR)))) %*% 
	SigmaCAR%*%solve(diag(sqrt(diag(SigmaCAR))))
	
print(
    xtable(W, 
      align = c('c',rep('c', times = length(W[1,]))),
      digits = c(0,rep(0, times = 7))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

SigmaSAR[lower.tri(SigmaSAR)] = NA
SigmaSAR[upper.tri(SigmaSAR)] = corrSAR[upper.tri(SigmaSAR)]
print(
    xtable(SigmaSAR, 
      align = c('c',rep('c', times = length(SigmaSAR[1,]))),
      digits = c(0,rep(2, times = 7))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

SigmaCAR[lower.tri(SigmaCAR)] = NA
SigmaCAR[upper.tri(SigmaCAR)] = corrCAR[upper.tri(SigmaCAR)]
print(
    xtable(SigmaCAR, 
      align = c('c',rep('c', times = length(SigmaCAR[1,]))),
      digits = c(0,rep(2, times = 7))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

################################################################################
#-------------------------------------------------------------------------------
#              Simulations of SAR and CAR
#-------------------------------------------------------------------------------
################################################################################

# create 20 x 20 grid of spatial locations
xy = expand.grid(x = seq(1, 20), y = seq(1, 20))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
Wmat = (Distmat < 1.01)*1
diag(Wmat) = 0
# create row-standardized W-matrix
Wbarmat = Wmat/apply(Wmat,1,sum)
# minimum and maximum eigenvalues for W
min(eigen(Wmat)$values)
max(eigen(Wmat)$values)
# parameter space for rho
1/min(eigen(Wmat)$values)
1/max(eigen(Wmat)$values)

# SAR binary weights covariance matrix with rho parameter at 0.25
SigmaSAR2 = solve(diag(400) - 0.25*Wmat) %*% solve(diag(400) - 0.25*t(Wmat))
# Cholesky decomposition of previous covariance matrix
SigmaSAR2.chol = chol(SigmaSAR2)
# CAR binary weights covariance matrix with rho parameter at 0.25
SigmaCAR = solve(diag(400) - 0.25*Wmat)
# Cholesky decomposition of previous covariance matrix
SigmaCAR.chol = chol(SigmaCAR)
# SAR binary weights covariance matrix with rho parameter at 0.20
SigmaSAR1 = solve(diag(400) - 0.2*Wmat)%*%solve(diag(400) - 0.2*t(Wmat))
# Cholesky decomposition of previous covariance matrix
SigmaSAR1.chol = chol(SigmaSAR1)
# K matrix for CAR model when using row-standardized weights
KWCAR = diag(1/apply(Wmat, 1, sum))
# CAR row-standardized weights covariance matrix with rho parameter at 0.95
SigmaWCAR = solve(diag(400) - 0.95*Wbarmat) %*% KWCAR
# Cholesky decomposition of previous covariance matrix
SigmaWCAR.chol = chol(SigmaWCAR)
#DTief <- diag(sqrt(as.vector(rowsum)))
#KTiefCAR <- solve(DTief)
#SigmaTiefCAR <- solve(diag(400)-0.95*WTiefmat)%*%KTiefCAR
#SigmaTiefCAR.chol <- chol(SigmaTiefCAR)
# SMA binary weights covariance matrix with rho parameter at 0.25
SigmaSMA = (0.25*Wmat + diag(400)) %*% (0.25*t(Wmat) + diag(400))
# Cholesky decomposition of previous covariance matrix
SigmaSMA.chol = chol(SigmaSMA)
set.seed(273)
# simulate i.i.d. dataset
dataiid = data.frame(as.data.frame(xy), z = rnorm(400, 0, 1))
# simulate SAR2 dataset
datasar2 = data.frame(as.data.frame(xy), z = SigmaSAR2.chol %*% rnorm(400,0,1))
# simulate CAR dataset
datacar = data.frame(as.data.frame(xy), z = SigmaCAR.chol%*%rnorm(400,0,1))
# simulate SAR1 dataset
datasar1 = data.frame(as.data.frame(xy), z = SigmaSAR1.chol%*%rnorm(400,0,1))
# simulate WCAR dataset
datawcar = data.frame(as.data.frame(xy), z = SigmaWCAR.chol%*%rnorm(400,0,1))
#data <- SigmaTiefCAR.chol%*%rnorm(400,0,1)
#datatiefcar <- matrix(data,nrow=20,byrow=T)
# simulate WCAR dataset
datasma = data.frame(as.data.frame(xy), z = SigmaSMA.chol%*%rnorm(400,0,1))

file_name = "figures/SimAutoreg"

pdf(paste0(file_name,'.pdf'), width = 8, height = 11)   

	cex_all = 3.2
	
	layout(matrix(1:12, ncol = 4, byrow = TRUE), widths = c(4,1,4,1))

	#i.i.d. simulation
	old.par = par(mar = c(1,2,3,0))
	cip = classIntervals(dataiid$z, n = 12, style = 'fisher')
	palp = viridis(12)
	cip_colors = findColours(cip, palp)
	plot(dataiid$x, dataiid$y, col = cip_colors, 
		bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = cex_all,
		main = expression(paste("Independent Normals")),
		cex.main = 1.9)
	par(mar = c(0,0,0,0))
	plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
	addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 1.3)
  
	#SAR simulation with binary weights and rho = 0.2
	old.par = par(mar = c(1,2,3,0))
	cip = classIntervals(datasar1$z, n = 12, style = 'fisher')
	palp = viridis(12)
	cip_colors = findColours(cip, palp)
	plot(datasar1$x, datasar1$y, col = cip_colors, 
		bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = cex_all,
		main = expression(paste("SAR with ", rho[SAR], " = 0.2")),
		cex.main = 1.9)
	par(mar = c(0,0,0,0))
	plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
	addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 1.3)

	#SAR simulation with binary weights and rho = 0.25
	old.par = par(mar = c(1,2,3,0))
	cip = classIntervals(datasar2$z, n = 12, style = 'fisher')
	palp = viridis(12)
	cip_colors = findColours(cip, palp)
	plot(datasar2$x, datasar2$y, col = cip_colors, 
		bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = cex_all,
		main = expression(paste("SAR with ", rho[SAR], " = 0.25")),
		cex.main = 1.9)
	par(mar = c(0,0,0,0))
	plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
	addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 1.3)

	#CAR simulation with binary weights and rho = 0.25
	old.par = par(mar = c(1,2,3,0))
	cip = classIntervals(datacar$z, n = 12, style = 'fisher')
	palp = viridis(12)
	cip_colors = findColours(cip, palp)
	plot(datacar$x, datacar$y, col = cip_colors, 
		bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = cex_all,
		main = expression(paste("CAR with ", rho[CAR], " = 0.25")),
		cex.main = 1.9)
	par(mar = c(0,0,0,0))
	plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
	addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 1.3)

	# CAR simulation with row-standardized weights and rho = 0.95
	old.par = par(mar = c(1,2,3,0))
	cip = classIntervals(datawcar$z, n = 12, style = 'fisher')
	palp = viridis(12)
	cip_colors = findColours(cip, palp)
	plot(datawcar$x, datacar$y, col = cip_colors, 
		bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = cex_all,
		main = expression(paste("WCAR with ", rho[CAR], " = 0.95")),
		cex.main = 1.9)
	par(mar = c(0,0,0,0))
	plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
	addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 1.3)

	# SMA simulation with binary weights and rho = 0.25
	old.par = par(mar = c(1,2,3,0))
	cip = classIntervals(datasma$z, n = 12, style = 'fisher')
	palp = viridis(12)
	cip_colors = findColours(cip, palp)
	plot(datasma$x, datasma$y, col = cip_colors, 
		bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = cex_all,
		main = expression(paste("SMA with ", rho[SMA], " = 0.25")),
		cex.main = 1.9)
	par(mar = c(0,0,0,0))
	plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
	addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 1.3)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

