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

file_name = "SimAutoreg"

pdf(paste0(file_name,'.pdf'), width = 8, height = 11)   

	layout(matrix(1:12, ncol = 4, byrow = TRUE), widths = c(4,1,4,1))

	#i.i.d. simulation
	old.par = par(mar = c(1,2,3,0))
	cip = classIntervals(dataiid$z, n = 12, style = 'fisher')
	palp = viridis(12)
	cip_colors = findColours(cip, palp)
	plot(dataiid$x, dataiid$y, col = cip_colors, bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = 3,
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
	plot(datasar1$x, datasar1$y, col = cip_colors, bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = 3,
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
	plot(datasar2$x, datasar2$y, col = cip_colors, bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = 3,
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
	plot(datacar$x, datacar$y, col = cip_colors, bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = 3,
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
	plot(datawcar$x, datacar$y, col = cip_colors, bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = 3,
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
	plot(datasma$x, datasma$y, col = cip_colors, bty = 'n', xaxt = 'n', yaxt = 'n', 
		pch = 15, cex = 3,
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

==========================================================================================================

# Code for Figure 7.2
# Plots of marginal correlations versus rho in SAR and CAR models
# First for a 3x3 square lattice without row standardization
W0 <- matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,byrow=T)
id3 <- diag(3)
ze3 <- matrix(0,3,3)
W1 <- cbind(W0,id3,ze3)
W2 <- cbind(id3,W0,id3)
W3 <- cbind(ze3,id3,W0)
W <- rbind(W1,W2,W3)
K <- diag(9)
# reciprocals of smallest and largest eigenvalues of this W are +-0.35355
count <- 0
nborpairs <- matrix(0,12,3)
for(i in 1:8){
 for(j in i:9){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,707,1)
margcorrCAR <- matrix(0,12,707)
margcorrSAR <- matrix(0,12,707)
for(i in 1:12){
 for(j in 1:707){
  rho[j,1] <- (j-354)/1000
  SigmaCAR <- solve(diag(9)-rho[j,1]*W)%*%K
  SigmaSAR <- solve(diag(9)-rho[j,1]*W)%*%K%*%solve(diag(9)-rho[j,1]*t(W))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
par(mfrow=c(2,2))
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
for(i in 2:12){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_CAR")
for(i in 2:12){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}

# Then for a 6x6 square lattice without row standardization
W0 <- matrix(c(0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0),nrow=6,byrow=T)
id6 <- diag(6)
ze6 <- matrix(0,6,6)
W1 <- cbind(W0,id6,ze6,ze6,ze6,ze6)
W2 <- cbind(id6,W0,id6,ze6,ze6,ze6)
W3 <- cbind(ze6,id6,W0,id6,ze6,ze6)
W4 <- cbind(ze6,ze6,id6,W0,id6,ze6)
W5 <- cbind(ze6,ze6,ze6,id6,W0,id6)
W6 <- cbind(ze6,ze6,ze6,ze6,id6,W0)
W <- rbind(W1,W2,W3,W4,W5,W6)
K <- diag(36)
# reciprocals of smallest and largest eigenvalues of this W are +-0.27748
count <- 0
nborpairs <- matrix(0,60,3)
for(i in 1:35){
 for(j in i:36){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,545,1)
margcorrCAR <- matrix(0,60,545)
margcorrSAR <- matrix(0,60,545)
for(i in 1:60){
 for(j in 1:545){
  rho[j,1] <- (j-278)/1000
  SigmaCAR <- solve(diag(36)-rho[j,1]*W)%*%K
  SigmaSAR <- solve(diag(36)-rho[j,1]*W)%*%K%*%solve(diag(36)-rho[j,1]*t(W))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
for(i in 2:60){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_CAR")
for(i in 2:60){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
========================================================================================================

# Code for Figure 7.3
# Similar to code for Figure 7.2, but with row standardization, first for the 3x3 lattice:
W0 <- matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,byrow=T)
id3 <- diag(3)
ze3 <- matrix(0,3,3)
W1 <- cbind(W0,id3,ze3)
W2 <- cbind(id3,W0,id3)
W3 <- cbind(ze3,id3,W0)
W <- rbind(W1,W2,W3)
rowsum <- W%*%matrix(1,9,1)
D <- diag(as.vector(rowsum))
Wbar <- solve(D)%*%W
K <- solve(D)
count <- 0
nborpairs <- matrix(0,12,3)
for(i in 1:8){
 for(j in i:9){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,1999,1)
margcorrCAR <- matrix(0,12,1999)
margcorrSAR <- matrix(0,12,1999)
for(i in 1:12){
 for(j in 1:1999){
  rho[j,1] <- (j-1000)/1000
  SigmaCAR <- solve(diag(9)-rho[j,1]*Wbar)%*%K
  SigmaSAR <- solve(diag(9)-rho[j,1]*Wbar)%*%K%*%solve(diag(9)-rho[j,1]*t(Wbar))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
par(mfrow=c(2,2))
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
for(i in 2:12){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_CAR")
for(i in 2:12){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}

# Then for a 6x6 square lattice:
W0 <- matrix(c(0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0),nrow=6,byrow=T)
id6 <- diag(6)
ze6 <- matrix(0,6,6)
W1 <- cbind(W0,id6,ze6,ze6,ze6,ze6)
W2 <- cbind(id6,W0,id6,ze6,ze6,ze6)
W3 <- cbind(ze6,id6,W0,id6,ze6,ze6)
W4 <- cbind(ze6,ze6,id6,W0,id6,ze6)
W5 <- cbind(ze6,ze6,ze6,id6,W0,id6)
W6 <- cbind(ze6,ze6,ze6,ze6,id6,W0)
W <- rbind(W1,W2,W3,W4,W5,W6)
rowsum <- W%*%matrix(1,36,1)
D <- diag(as.vector(rowsum))
Wbar <- solve(D)%*%W
K <- solve(D)
count <- 0
nborpairs <- matrix(0,60,3)
for(i in 1:35){
 for(j in i:36){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,1999,1)
margcorrSAR <- matrix(0,60,1999)
margcorrCAR <- matrix(0,60,1999)
for(i in 1:60){
 for(j in 1:1999){
  rho[j,1] <- (j-1000)/1000
  SigmaCAR <- solve(diag(36)-rho[j,1]*Wbar)%*%K
  SigmaSAR <- solve(diag(36)-rho[j,1]*Wbar)%*%K%*%solve(diag(36)-rho[j,1]*t(Wbar))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
for(i in 2:60){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_CAR")
for(i in 2:60){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
============================================================================================================

# Now, for Figure 7.4, do a similar thing for the lattice of the 48 contiguous United States plus DC
library(spdep)
library(maps)
library(maptools)
library(classInt)
library(RColorBrewer)

## Create an adjacency matrix for the states in the US
usa.state = map(database="state", fill=TRUE, plot=FALSE)
state.ID <- sapply(strsplit(usa.state$names, ":"), function(x) x[1])
usa.poly = map2SpatialPolygons(usa.state, IDs=state.ID)
usa.nb = poly2nb(usa.poly)
usa.adj.mat = nb2mat(usa.nb, style="B")

## Write the 0-1 adjacency matrix
W = usa.adj.mat
W[(W>0)] = 1

# First do without row standardization:
Wusa <- matrix(scan("geostatsbook/usaneighbors.txt"),ncol=49,byrow=T)
K <- diag(49)
# Reciprocals of the smallest and largest eigenvalues are -0.3489 and 0.1847, respectively
count <- 0
nborpairs <- matrix(0,109,3)
for(i in 1:48){
 for(j in i:49){
  if(Wusa[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,53,1)
margcorrSAR <- matrix(0,109,53)
margcorrCAR <- matrix(0,109,53)
for(i in 1:109){
 for(j in 1:53){
  rho[j,1] <- -0.08+(j-27)/100
  SigmaCAR <- solve(diag(49)-rho[j,1]*Wusa)%*%K
  SigmaSAR <- solve(diag(49)-rho[j,1]*Wusa)%*%K%*%solve(diag(49)-rho[j,1]*t(Wusa))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
par(mfrow=c(2,2))
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
for(i in 2:109){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_CAR")
for(i in 2:109){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}

# Now repeat, but with row standardization
Wusa <- matrix(scan("geostatsbook/usaneighbors.txt"),ncol=49,byrow=T)
rowsum <- Wusa%*%matrix(1,49,1)
D <- diag(as.vector(rowsum))
Wusabar <- solve(D)%*%Wusa
K <- solve(D)
count <- 0
nborpairs <- matrix(0,109,3)
for(i in 1:48){
 for(j in i:49){
  if(Wusa[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,239,1)
margcorrSAR <- matrix(0,109,239)
margcorrCAR <- matrix(0,109,239)
for(i in 1:109){
 for(j in 1:239){
  rho[j,1] <- -0.2+(j-120)/100
  SigmaCAR <- solve(diag(49)-rho[j,1]*Wusabar)%*%K
  SigmaSAR <- solve(diag(49)-rho[j,1]*Wusabar)%*%K%*%solve(diag(49)-rho[j,1]*t(Wusabar))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
for(i in 2:109){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
abline(v=-1,lty=3)
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_CAR")
for(i in 2:109){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
abline(v=-1,lty=3)
=============================================================================================================

# Some additional variations not included in the text
# Add diagonal neighbors to 6x6:
W0 <- matrix(c(0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0),nrow=6,byrow=T)
id6 <- matrix(c(1,1,0,0,0,0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,1,1),nrow=6,byrow=T)
ze6 <- matrix(0,6,6)
W1 <- cbind(W0,id6,ze6,ze6,ze6,ze6)
W2 <- cbind(id6,W0,id6,ze6,ze6,ze6)
W3 <- cbind(ze6,id6,W0,id6,ze6,ze6)
W4 <- cbind(ze6,ze6,id6,W0,id6,ze6)
W5 <- cbind(ze6,ze6,ze6,id6,W0,id6)
W6 <- cbind(ze6,ze6,ze6,ze6,id6,W0)
W <- rbind(W1,W2,W3,W4,W5,W6)
K <- diag(36)
count <- 0
nborpairs <- matrix(0,110,3)
for(i in 1:35){
 for(j in i:36){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,45,1)
margcorrSAR <- matrix(0,110,45)
margcorrCAR <- matrix(0,110,45)
for(i in 1:110){
 for(j in 1:45){
  rho[j,1] <- -.08+(j-23)/100
  SigmaCAR <- solve(diag(36)-rho[j,1]*W)%*%K
  SigmaSAR <- solve(diag(36)-rho[j,1]*W)%*%K%*%solve(diag(36)-rho[j,1]*t(W))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
par(mfrow=c(2,2))
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
for(i in 2:110){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_CAR")
for(i in 2:110){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}

# anisotropy, with row standardization
W0 <- matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,byrow=T)
id3 <- 0.5*diag(3)
ze3 <- matrix(0,3,3)
W1 <- cbind(W0,id3,ze3)
W2 <- cbind(id3,W0,id3)
W3 <- cbind(ze3,id3,W0)
W <- rbind(W1,W2,W3)
rowsum <- W%*%matrix(1,9,1)
D <- diag(as.vector(rowsum))
Wbar <- solve(D)%*%W
K <- solve(D)
count <- 0
nborpairs <- matrix(0,12,3)
for(i in 1:8){
 for(j in i:9){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,199,1)
margcorrCAR <- matrix(0,12,199)
margcorrSAR <- matrix(0,12,199)
for(i in 1:12){
 for(j in 1:199){
  rho[j,1] <- (j-100)/100
  SigmaCAR <- solve(diag(9)-rho[j,1]*Wbar)%*%K
  SigmaSAR <- solve(diag(9)-rho[j,1]*Wbar)%*%K%*%solve(diag(9)-rho[j,1]*t(Wbar))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
par(mfrow=c(2,2))
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l")
for(i in 2:12){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l")}
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l")
for(i in 2:12){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l")}

# Then for a 6x6 square lattice
W0 <- matrix(c(0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0),nrow=6,byrow=T)
id6 <- 0.5*diag(6)
ze6 <- matrix(0,6,6)
W1 <- cbind(W0,id6,ze6,ze6,ze6,ze6)
W2 <- cbind(id6,W0,id6,ze6,ze6,ze6)
W3 <- cbind(ze6,id6,W0,id6,ze6,ze6)
W4 <- cbind(ze6,ze6,id6,W0,id6,ze6)
W5 <- cbind(ze6,ze6,ze6,id6,W0,id6)
W6 <- cbind(ze6,ze6,ze6,ze6,id6,W0)
W <- rbind(W1,W2,W3,W4,W5,W6)
rowsum <- W%*%matrix(1,36,1)
D <- diag(as.vector(rowsum))
Wbar <- solve(D)%*%W
K <- solve(D)
count <- 0
nborpairs <- matrix(0,60,3)
for(i in 1:35){
 for(j in i:36){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,199,1)
margcorrSAR <- matrix(0,60,199)
margcorrCAR <- matrix(0,60,199)
for(i in 1:60){
 for(j in 1:199){
  rho[j,1] <- (j-100)/100
  SigmaCAR <- solve(diag(36)-rho[j,1]*Wbar)%*%K
  SigmaSAR <- solve(diag(36)-rho[j,1]*Wbar)%*%K%*%solve(diag(36)-rho[j,1]*t(Wbar))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
for(i in 2:60){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_CAR")
for(i in 2:60){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
==========================================================================================================

# Now repeat for SMA models, for Figure 7.5:
# First for a 6x6 square lattice without row standardization:
W0 <- matrix(c(0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0),nrow=6,byrow=T)
id6 <- diag(6)
ze6 <- matrix(0,6,6)
W1 <- cbind(W0,id6,ze6,ze6,ze6,ze6)
W2 <- cbind(id6,W0,id6,ze6,ze6,ze6)
W3 <- cbind(ze6,id6,W0,id6,ze6,ze6)
W4 <- cbind(ze6,ze6,id6,W0,id6,ze6)
W5 <- cbind(ze6,ze6,ze6,id6,W0,id6)
W6 <- cbind(ze6,ze6,ze6,ze6,id6,W0)
W <- rbind(W1,W2,W3,W4,W5,W6)
K <- diag(36)
# reciprocals of smallest and largest eigenvalues of this W are +-0.27748
count <- 0
nborpairs <- matrix(0,60,3)
for(i in 1:35){
 for(j in i:36){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,545,1)
margcorrSMA <- matrix(0,60,545)
for(i in 1:60){
 for(j in 1:545){
  rho[j,1] <- (j-278)/1000
  SigmaSMA <- (diag(36)+rho[j,1]*W)%*%K%*%(diag(36)+rho[j,1]*t(W))
  margcorrSMA[i,j] <- SigmaSMA[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSMA[nborpairs[i,2],nborpairs[i,2]]*SigmaSMA[nborpairs[i,3],nborpairs[i,3]])
}}
par(mfrow=c(2,2))
plot(rho,margcorrSMA[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SMA")
for(i in 2:60){
lines(rho,margcorrSMA[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}

# Then for a 6x6 square lattice with row standardization:
W0 <- matrix(c(0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0),nrow=6,byrow=T)
id6 <- diag(6)
ze6 <- matrix(0,6,6)
W1 <- cbind(W0,id6,ze6,ze6,ze6,ze6)
W2 <- cbind(id6,W0,id6,ze6,ze6,ze6)
W3 <- cbind(ze6,id6,W0,id6,ze6,ze6)
W4 <- cbind(ze6,ze6,id6,W0,id6,ze6)
W5 <- cbind(ze6,ze6,ze6,id6,W0,id6)
W6 <- cbind(ze6,ze6,ze6,ze6,id6,W0)
W <- rbind(W1,W2,W3,W4,W5,W6)
rowsum <- W%*%matrix(1,36,1)
D <- diag(as.vector(rowsum))
Wbar <- solve(D)%*%W
K <- solve(D)
count <- 0
nborpairs <- matrix(0,60,3)
for(i in 1:35){
 for(j in i:36){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,1999,1)
margcorrSMA <- matrix(0,60,1999)
for(i in 1:60){
 for(j in 1:1999){
  rho[j,1] <- (j-1000)/1000
  SigmaSMA <- (diag(36)+rho[j,1]*Wbar)%*%K%*%(diag(36)+rho[j,1]*t(Wbar))
  margcorrSMA[i,j] <- SigmaSMA[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSMA[nborpairs[i,2],nborpairs[i,2]]*SigmaSMA[nborpairs[i,3],nborpairs[i,3]])
}}
plot(rho,margcorrSMA[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SMA")
for(i in 2:60){
lines(rho,margcorrSMA[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}

# Now do a similar thing for the lattice of the 48 contiguous United States plus DC
library(spdep)
library(maps)
library(maptools)
library(classInt)
library(RColorBrewer)

## Create an adjacency matrix for the states in the US
usa.state = map(database="state", fill=TRUE, plot=FALSE)
state.ID <- sapply(strsplit(usa.state$names, ":"), function(x) x[1])
usa.poly = map2SpatialPolygons(usa.state, IDs=state.ID)
usa.nb = poly2nb(usa.poly)
usa.adj.mat = nb2mat(usa.nb, style="B")

## Write the 0-1 adjacency matrix
W = usa.adj.mat
W[(W>0)] = 1

# First do without row standardization:
Wusa <- matrix(scan("geostatsbook/usaneighbors.txt"),ncol=49,byrow=T)
K <- diag(49)
count <- 0
nborpairs <- matrix(0,109,3)
for(i in 1:48){
 for(j in i:49){
  if(Wusa[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,53,1)
margcorrSMA <- matrix(0,109,53)
for(i in 1:109){
 for(j in 1:53){
  rho[j,1] <- -0.08+(j-27)/100
  SigmaSMA <- (rho[j,1]*Wusa+diag(49))%*%K%*%(rho[j,1]*t(Wusa)+diag(49))
  margcorrSMA[i,j] <- SigmaSMA[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSMA[nborpairs[i,2],nborpairs[i,2]]*SigmaSMA[nborpairs[i,3],nborpairs[i,3]])
}}
plot(rho,margcorrSMA[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SMA")
for(i in 2:109){
lines(rho,margcorrSMA[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}

# Now repeat, but with row standardization
Wusa <- matrix(scan("geostatsbook/usaneighbors.txt"),ncol=49,byrow=T)
rowsum <- Wusa%*%matrix(1,49,1)
D <- diag(as.vector(rowsum))
Wusabar <- solve(D)%*%Wusa
K <- solve(D)
count <- 0
nborpairs <- matrix(0,109,3)
for(i in 1:48){
 for(j in i:49){
  if(Wusa[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,239,1)
margcorrSMA <- matrix(0,109,239)
for(i in 1:109){
 for(j in 1:239){
  rho[j,1] <- -0.2+(j-120)/100
  SigmaSMA <- (rho[j,1]*Wusabar+diag(49))%*%K%*%(rho[j,1]*t(Wusabar)+diag(49))
  margcorrSMA[i,j] <- SigmaSMA[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSMA[nborpairs[i,2],nborpairs[i,2]]*SigmaSMA[nborpairs[i,3],nborpairs[i,3]])
}}
plot(rho,margcorrSMA[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SMA")
for(i in 2:109){
lines(rho,margcorrSMA[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}

-------------------------------------------------------------------------------------------------------

# For Exercise 7.10a, investigate behavior of SAR model's marginal correlations for rho_SAR outside the interval within # which I-rho_CAR*W is positive definite.
# First for 3x3 lattice
W0 <- matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,byrow=T)
id3 <- diag(3)
ze3 <- matrix(0,3,3)
W1 <- cbind(W0,id3,ze3)
W2 <- cbind(id3,W0,id3)
W3 <- cbind(ze3,id3,W0)
W <- rbind(W1,W2,W3)
K <- diag(9)
# reciprocals of smallest and largest eigenvalues of this W are +-0.35355
count <- 0
nborpairs <- matrix(0,12,3)
for(i in 1:8){
 for(j in i:9){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,21,1)
margcorrSAR <- matrix(0,12,21)
for(i in 1:12){
 for(j in 1:21){
  rho[j,1] <- (j-11)/10
  SigmaSAR <- solve(diag(9)-rho[j,1]*W)%*%K%*%solve(diag(9)-rho[j,1]*t(W))
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
par(mfrow=c(2,1))
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
for(i in 2:12){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}

# Then for a 6x6 square lattice
W0 <- matrix(c(0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0),nrow=6,byrow=T)
id6 <- diag(6)
ze6 <- matrix(0,6,6)
W1 <- cbind(W0,id6,ze6,ze6,ze6,ze6)
W2 <- cbind(id6,W0,id6,ze6,ze6,ze6)
W3 <- cbind(ze6,id6,W0,id6,ze6,ze6)
W4 <- cbind(ze6,ze6,id6,W0,id6,ze6)
W5 <- cbind(ze6,ze6,ze6,id6,W0,id6)
W6 <- cbind(ze6,ze6,ze6,ze6,id6,W0)
W <- rbind(W1,W2,W3,W4,W5,W6)
K <- diag(36)
# reciprocals of smallest and largest eigenvalues of this W are +-0.27748
count <- 0
nborpairs <- matrix(0,60,3)
for(i in 1:35){
 for(j in i:36){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,21,1)
margcorrSAR <- matrix(0,60,21)
for(i in 1:60){
 for(j in 1:21){
  rho[j,1] <- (j-11)/10
  SigmaSAR <- solve(diag(36)-rho[j,1]*W)%*%K%*%solve(diag(36)-rho[j,1]*t(W))
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
for(i in 2:60){
lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
==============================================================================================================

# Now, for Exercise 7.10c, do something similar for the "autocorrelated CAR model":
# First for 3x3 lattice
W0 <- matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,byrow=T)
id3 <- diag(3)
ze3 <- matrix(0,3,3)
W1 <- cbind(W0,id3,ze3)
W2 <- cbind(id3,W0,id3)
W3 <- cbind(ze3,id3,W0)
W <- rbind(W1,W2,W3)
rowsum <- W%*%matrix(1,9,1)
D <- diag(as.vector(rowsum))
K <- solve(D)
Wac <- matrix(0,9,9)
for(i in 1:9){
 for(j in 1:9){
  Wac[i,j] <- sqrt(rowsum[j,1])/sqrt(rowsum[i,1])
}}
Wac <- Wac*W
count <- 0
nborpairs <- matrix(0,12,3)
for(i in 1:8){
 for(j in i:9){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,707,1)
margcorrCAR <- matrix(0,12,707)
margcorrSAR <- matrix(0,12,707)
for(i in 1:12){
 for(j in 1:707){
  rho[j,1] <- (j-354)/1000
  SigmaCAR <- solve(diag(9)-rho[j,1]*Wac)%*%K
  SigmaSAR <- solve(diag(9)-rho[j,1]*Wac)%*%K%*%solve(diag(9)-rho[j,1]*t(Wac))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
par(mfrow=c(2,1))
#plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
#for(i in 2:12){
#lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_CAR")
for(i in 2:12){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}

# Then for a 6x6 lattice:
W0 <- matrix(c(0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0),nrow=6,byrow=T)
id6 <- diag(6)
ze6 <- matrix(0,6,6)
W1 <- cbind(W0,id6,ze6,ze6,ze6,ze6)
W2 <- cbind(id6,W0,id6,ze6,ze6,ze6)
W3 <- cbind(ze6,id6,W0,id6,ze6,ze6)
W4 <- cbind(ze6,ze6,id6,W0,id6,ze6)
W5 <- cbind(ze6,ze6,ze6,id6,W0,id6)
W6 <- cbind(ze6,ze6,ze6,ze6,id6,W0)
W <- rbind(W1,W2,W3,W4,W5,W6)
rowsum <- W%*%matrix(1,36,1)
D <- diag(as.vector(rowsum))
Wbar <- solve(D)%*%W
K <- solve(D)
Wac <- matrix(0,36,36)
for(i in 1:36){
 for(j in 1:36){
  Wac[i,j] <- sqrt(rowsum[j,1])/sqrt(rowsum[i,1])
}}
Wac <- Wac*W
count <- 0
nborpairs <- matrix(0,60,3)
for(i in 1:35){
 for(j in i:36){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
rho <- matrix(0,545,1)
margcorrCAR <- matrix(0,60,545)
margcorrSAR <- matrix(0,60,545)
for(i in 1:60){
 for(j in 1:545){
  rho[j,1] <- (j-278)/1000
  SigmaCAR <- solve(diag(36)-rho[j,1]*Wac)%*%K
  SigmaSAR <- solve(diag(36)-rho[j,1]*Wac)%*%K%*%solve(diag(36)-rho[j,1]*t(Wac))
  margcorrCAR[i,j] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
  margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaSAR[nborpairs[i,2],nborpairs[i,2]]*SigmaSAR[nborpairs[i,3],nborpairs[i,3]])
}}
#plot(rho,margcorrSAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_SAR")
#for(i in 2:60){
#lines(rho,margcorrSAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
plot(rho,margcorrCAR[1,],ylim=c(-1.0,1.0),type="l",lwd=0.5,ylab="marginal correlation",xlab="rho_CAR")
for(i in 2:60){
lines(rho,margcorrCAR[i,],ylim=c(-1.0,1.0),type="l",lwd=0.5)}
=========================================================================================================

# Now, for Exercise 7.10b, do something similar for the unconstrained CAR model:
# First for the 3x3 lattice
W0 <- matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,byrow=T)
id3 <- diag(3)
ze3 <- matrix(0,3,3)
W1 <- cbind(W0,id3,ze3)
W2 <- cbind(id3,W0,id3)
W3 <- cbind(ze3,id3,W0)
W <- rbind(W1,W2,W3)
count <- 0
nborpairs <- matrix(0,12,3)
for(i in 1:8){
 for(j in i:9){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
margcorrCAR <- matrix(0,12,401)
C <- matrix(0,9,9)
phi <- matrix(0,401,1)
rowsum <- W%*%matrix(1,9,1)
D <- diag(as.vector(rowsum))
for(i in 1:12){
for(k in 1:401){
phi[k,1] <- (k-201)/10
Kinv <- diag(9)+abs(phi[k,1])*D
K <- solve(Kinv)
for(l in 1:9){
 for(m in 1:9){
  C[l,m] <- phi[k,1]*W[l,m]/(1+abs(phi[k,1])*rowsum[l,1])
}}
SigmaCAR <- solve(diag(9)-C)%*%K
margcorrCAR[i,k] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
}}
par(mfrow=c(2,1))
plot(phi,margcorrCAR[1,],xlim=c(-22,22),ylim=c(-1,1),type="l",lwd=0.5,ylab="marginal correlation",xlab="phi")
for(i in 2:12){
lines(phi,margcorrCAR[i,],xlim=c(-22,22),ylim=c(-1,1),type="l",lwd=0.5)}

# Then for 6x6 lattice
W0 <- matrix(c(0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,0,1,0),nrow=6,byrow=T)
id6 <- diag(6)
ze6 <- matrix(0,6,6)
W1 <- cbind(W0,id6,ze6,ze6,ze6,ze6)
W2 <- cbind(id6,W0,id6,ze6,ze6,ze6)
W3 <- cbind(ze6,id6,W0,id6,ze6,ze6)
W4 <- cbind(ze6,ze6,id6,W0,id6,ze6)
W5 <- cbind(ze6,ze6,ze6,id6,W0,id6)
W6 <- cbind(ze6,ze6,ze6,ze6,id6,W0)
W <- rbind(W1,W2,W3,W4,W5,W6)
count <- 0
nborpairs <- matrix(0,60,3)
for(i in 1:35){
 for(j in i:36){
  if(W[i,j]!=0){
count <- count+1
nborpairs[count,1] <- count
nborpairs[count,2] <- i
nborpairs[count,3] <- j
}}}
margcorrCAR <- matrix(0,60,401)
C <- matrix(0,36,36)
phi <- matrix(0,401,1)
rowsum <- W%*%matrix(1,36,1)
D <- diag(as.vector(rowsum))
for(i in 1:60){
for(k in 1:401){
phi[k,1] <- (k-201)/10
Kinv <- diag(36)+abs(phi[k,1])*D
K <- solve(Kinv)
for(l in 1:36){
 for(m in 1:36){
  C[l,m] <- phi[k,1]*W[l,m]/(1+abs(phi[k,1])*rowsum[l,1])
}}
SigmaCAR <- solve(diag(36)-C)%*%K
margcorrCAR[i,k] <- SigmaCAR[nborpairs[i,2],nborpairs[i,3]]/sqrt(SigmaCAR[nborpairs[i,2],nborpairs[i,2]]*SigmaCAR[nborpairs[i,3],nborpairs[i,3]])
}}
plot(phi,margcorrCAR[1,],xlim=c(-22,22),ylim=c(-1,1),type="l",lwd=0.5,ylab="marginal correlation",xlab="phi")
for(i in 2:60){
lines(phi,margcorrCAR[i,],xlim=c(-22,22),ylim=c(-1,1),type="l",lwd=0.5)}
==========================================================================================




