sec_path = 'Rcode/Chapter7/Section 7.7/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(xtable)

################################################################################
#-------------------------------------------------------------------------------
#                         Spatial Moving Averages
#-------------------------------------------------------------------------------
################################################################################

# ----------------------- Simple 7 site example --------------------------------

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
# reciprocals of smallest and largest eigenvalues of W 
1/min(eigen(W)$values)
1/max(eigen(W)$values)
K <- diag(7)
# SAR covariance and correlation matrices for row-standardized weights
SigmaSMA <- (diag(7) + 0.39*W)%*%K%*%(diag(7) + 0.39*t(W))
corrSigmaSMA = diag(sqrt(1/diag(SigmaSMA))) %*% 
	SigmaSMA%*% diag(sqrt(1/diag(SigmaSMA)))
# create matrix with variances on diagonal and correlations on off-diagonal
SigmaSMA[lower.tri(SigmaSMA)] = NA
SigmaSMA[upper.tri(SigmaSMA)] = corrSigmaSMA[upper.tri(SigmaSMA)]

# print matrices
print(
    xtable(SigmaSMA, 
      align = c('c',rep('c', times = length(W[1,]))),
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
#     Relationships between spatial weights and marginal correlations
#-------------------------------------------------------------------------------
################################################################################

# Next for a 6x6 square lattice without row standardization
xy = expand.grid(x = seq(1, 6), y = seq(1, 6))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 1.01)*1
diag(W) = 0
# reciprocals of smallest and largest eigenvalues of W 
1/min(eigen(W)$values)
1/max(eigen(W)$values)
K <- diag(36)
Kbar = diag(1/apply(W, 1, sum))
Wbar = W/apply(W,1,sum)

# get first-order only indexes for pairwise correlations from 6 x 6 grid (36 locations)
nborpairs6 <- NULL
for(i in 1:35){
	for(j in (i + 1):36){
		if(W[i,j]!=0) nborpairs6 = rbind(nborpairs6, c(i,j))
	}
}
rho6 = rep(NA, times = 545)
margcorrSMA = matrix(NA, dim(nborpairs6)[1], 545)
for(i in 1:dim(nborpairs6)[1]){
	for(j in 1:545){
		rho6[j] <- (j-278)/1000
		SigmaSMA <- (diag(36) + rho6[j]*W)%*%K%*%(diag(36) + rho6[j]*t(W))
		margcorrSMA[i,j] = SigmaSMA[nborpairs6[i,1],nborpairs6[i,2]]/
			sqrt(SigmaSMA[nborpairs6[i,1],nborpairs6[i,1]]*
			SigmaSMA[nborpairs6[i,2],nborpairs6[i,2]])
	}
}
rho6bar = rep(NA, times = 1999)
margcorrSMAbar = matrix(NA, dim(nborpairs6)[1], 1999)
for(i in 1:dim(nborpairs6)[1]){
	for(j in 1:1999){
		rho6bar[j] <- (j-1000)/1000
		SigmaSMA <- (diag(36) + rho6bar[j]*Wbar)%*%Kbar%*%(diag(36) + 
			rho6bar[j]*t(Wbar))
		margcorrSMAbar[i,j] = SigmaSMA[nborpairs6[i,1],nborpairs6[i,2]]/
			sqrt(SigmaSMA[nborpairs6[i,1],nborpairs6[i,1]]*
			SigmaSMA[nborpairs6[i,2],nborpairs6[i,2]])
	}
}

data(USboundary)
library(spdep)
# use spdep to create get neighbors from state polygons
usa_nb = poly2nb(USboundary)
# turn it into a matrix
usa_adj_mat = nb2mat(usa_nb, style="B")

## write the 0-1 adjacency matrix
K = diag(49)
W = usa_adj_mat
# create row-standardized matrices
Krs = diag(1/apply(W,1,sum))
Wrs = W/apply(W,1,sum)

# get the limit of parameter space for rho for binary weights
rho_min = 1/min(eigen(W)$values)
rho_max = 1/max(eigen(W)$values)
rho_min
rho_max
# get the limit of parameter space for rho for row-standardized weights
rhors_min = 1/min(eigen(Wrs)$values)
rhors_max = 1/max(eigen(Wrs)$values)
rhors_min
rhors_max

# get first-order only indexes for pairwise correlations
nborpairsusa <- NULL
for(i in 1:48){
	for(j in (i + 1):49){
		if(W[i,j]!=0) nborpairsusa = rbind(nborpairsusa, c(i,j))
	}
}

# range of correlations for rho for binary weights
rhousa = rep(NA, times = 53)
margcorrusaSMA = matrix(NA, dim(nborpairsusa)[1], 53)
for(i in 1:dim(nborpairsusa)[1]){
	for(j in 1:53){
		rhousa[j] <- -0.08+(j-27)/100
		SigmaSMA <- (diag(49)+rhousa[j]*W)%*%K%*%(diag(49)+rhousa[j]*t(W))
		margcorrusaSMA[i,j] <- SigmaSMA[nborpairsusa[i,1],nborpairsusa[i,2]]/
			sqrt(SigmaSMA[nborpairsusa[i,1],nborpairsusa[i,1]]*
			SigmaSMA[nborpairsusa[i,2],nborpairsusa[i,2]])
	}
}

# range of correlations for rho for binary weights
rhousars = rep(NA, times = 239)
margcorrusaSMArs = matrix(NA, dim(nborpairsusa)[1], 239)
for(i in 1:dim(nborpairsusa)[1]){
	for(j in 1:239){
		rhousars[j] <- -0.2+(j-120)/100
		SigmaSMA <- (diag(49)+rhousars[j]*Wrs)%*%Krs%*%(diag(49)+rhousars[j]*t(Wrs))
		margcorrusaSMArs[i,j] <- SigmaSMA[nborpairsusa[i,1],nborpairsusa[i,2]]/
			sqrt(SigmaSMA[nborpairsusa[i,1],nborpairsusa[i,1]]*
			SigmaSMA[nborpairsusa[i,2],nborpairsusa[i,2]])
	}
}


file_name = "figures/MargCorrSMA"

pdf(paste0(file_name,'.pdf'), width = 11, height = 11)   

	layout(matrix(1:4, ncol = 2, byrow = TRUE))

	old.par = par(mar = c(5,5,5,1))
	plot(rho6, margcorrSMA[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[SMA])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = '6 x 6 grid, binary weights', cex.main = 2)
	for(i in 2:dim(nborpairs6)[1]){
		lines(rho6,margcorrSMA[i,], lwd=1.0)}
	plot(rho6bar, margcorrSMAbar[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[SMA])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = '6 x 6 grid, row standardized', cex.main = 2)
	for(i in 2:dim(nborpairs6)[1]){
		lines(rho6bar,margcorrSMAbar[i,], lwd=1.0)}
	plot(rhousa, margcorrusaSMA[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[SMAR])),
		ylab = 'Marginal Correlation',cex.axis = 1.5, cex.lab = 2, 
		main = 'US States, binary weights', cex.main = 2)
	for(i in 2:dim(nborpairsusa)[1]){
		lines(rhousa,margcorrusaSMA[i,], lwd=1.0)}
	plot(rhousars, margcorrusaSMArs[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[SMA])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = 'US States, row-standardized', cex.main = 2)
	for(i in 2:dim(nborpairsusa)[1]){
		lines(rhousars,margcorrusaSMArs[i,], lwd=1.0)}

	par(old.par)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

