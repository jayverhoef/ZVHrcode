sec_path = 'Rcode/Chapter7/Section 7.5/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)

################################################################################
#-------------------------------------------------------------------------------
#     Relationships between spatial weights and marginal correlations
#-------------------------------------------------------------------------------
################################################################################

# -------------------------- binary weights ------------------------------------

# First for a 3x3 square lattice without row standardization
xy = expand.grid(x = seq(1, 3), y = seq(1, 3))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 1.01)*1
diag(W) = 0
K <- diag(9)
# reciprocals of smallest and largest eigenvalues of this W are +-0.35355
1/min(eigen(W)$values)
1/max(eigen(W)$values)
# get first-order only indexes for pairwise correlations from 3 x 3 grid (9 locations)
nborpairs3 <- NULL
for(i in 1:8){
	for(j in (i + 1):9){
		if(W[i,j]!=0) nborpairs3 = rbind(nborpairs3, c(i,j))
	}
}
rho3 = rep(NA, times = 707)
margcorrCAR3 = matrix(NA, dim(nborpairs3)[1], 707)
margcorrSAR3 = matrix(NA, dim(nborpairs3)[1], 707)
for(i in 1:dim(nborpairs3)[1]){
	for(j in 1:707){
		rho3[j] <- (j-354)/1000
		SigmaCAR <- solve(diag(9)-rho3[j]*W)%*%K
		SigmaSAR <- solve(diag(9)-rho3[j]*W)%*%K%*%solve(diag(9)-rho3[j]*t(W))
		margcorrCAR3[i,j] = SigmaCAR[nborpairs3[i,1],nborpairs3[i,2]]/
			sqrt(SigmaCAR[nborpairs3[i,1],nborpairs3[i,1]]*
			SigmaCAR[nborpairs3[i,2],nborpairs3[i,2]])
		margcorrSAR3[i,j] = SigmaSAR[nborpairs3[i,1],nborpairs3[i,2]]/
			sqrt(SigmaSAR[nborpairs3[i,1],nborpairs3[i,1]]*
			SigmaSAR[nborpairs3[i,2],nborpairs3[i,2]])
	}
}

# Next for a 6x6 square lattice without row standardization
xy = expand.grid(x = seq(1, 6), y = seq(1, 6))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 1.01)*1
diag(W) = 0
K <- diag(36)
# reciprocals of smallest and largest eigenvalues of this W are +-0.27747
1/min(eigen(W)$values)
1/max(eigen(W)$values)
# get first-order only indexes for pairwise correlations from 6 x 6 grid (9 locations)
nborpairs6 <- NULL
for(i in 1:35){
	for(j in (i + 1):36){
		if(W[i,j]!=0) nborpairs6 = rbind(nborpairs6, c(i,j))
	}
}
rho6 = rep(NA, times = 545)
margcorrCAR6 = matrix(NA, dim(nborpairs6)[1], 545)
margcorrSAR6 = matrix(NA, dim(nborpairs6)[1], 545)
for(i in 1:dim(nborpairs6)[1]){
	for(j in 1:545){
		rho6[j] <- (j-278)/1000
		SigmaCAR <- solve(diag(36)-rho6[j]*W)%*%K
		SigmaSAR <- solve(diag(36)-rho6[j]*W)%*%K%*%solve(diag(36)-rho6[j]*t(W))
		margcorrCAR6[i,j] = SigmaCAR[nborpairs6[i,1],nborpairs6[i,2]]/
			sqrt(SigmaCAR[nborpairs6[i,1],nborpairs6[i,1]]*
			SigmaCAR[nborpairs6[i,2],nborpairs6[i,2]])
		margcorrSAR6[i,j] = SigmaSAR[nborpairs6[i,1],nborpairs6[i,2]]/
			sqrt(SigmaSAR[nborpairs6[i,1],nborpairs6[i,1]]*
			SigmaSAR[nborpairs6[i,2],nborpairs6[i,2]])
	}
}

file_name = "MargCorrBinaryWts"

pdf(paste0(file_name,'.pdf'), width = 11, height = 11)   

	layout(matrix(1:4, ncol = 2, byrow = TRUE))

	old.par = par(mar = c(5,5,5,1))
	plot(rho3, margcorrSAR3[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[SAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = '3 x 3 grid', cex.main = 2)
	for(i in 2:dim(nborpairs3)[1]){
		lines(rho3,margcorrSAR3[i,], lwd=1.0)}
	plot(rho3, margcorrCAR3[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[CAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = '3 x 3 grid', cex.main = 2)
	for(i in 2:dim(nborpairs3)[1]){
		lines(rho3,margcorrCAR3[i,], lwd=1.0)}
	plot(rho6, margcorrSAR6[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[SAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = '6 x 6 grid', cex.main = 2)
	for(i in 2:dim(nborpairs6)[1]){
		lines(rho6,margcorrSAR6[i,], lwd=1.0)}
	plot(rho6, margcorrCAR6[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[CAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = '6 x 6 grid', cex.main = 2)
	for(i in 2:dim(nborpairs6)[1]){
		lines(rho6,margcorrCAR6[i,], lwd=1.0)}
	par(old.par)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))


# --------------------- row-standardized weights -------------------------------

# First for a 3x3 square lattice WITH row standardization
xy = expand.grid(x = seq(1, 3), y = seq(1, 3))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 1.01)*1
diag(W) = 0
K = diag(1/apply(W, 1, sum))
W = W/apply(W,1,sum)

# reciprocals of smallest and largest eigenvalues of this W are +-1
1/min(eigen(W)$values)
1/max(eigen(W)$values)
# get first-order only indexes for pairwise correlations from 3 x 3 grid (9 locations)
nborpairs3 <- NULL
for(i in 1:8){
	for(j in (i + 1):9){
		if(W[i,j]!=0) nborpairs3 = rbind(nborpairs3, c(i,j))
	}
}
rho3 = rep(NA, times = 1999)
margcorrCAR3 = matrix(NA, dim(nborpairs3)[1], 1999)
margcorrSAR3 = matrix(NA, dim(nborpairs3)[1], 1999)
for(i in 1:dim(nborpairs3)[1]){
	for(j in 1:1999){
		rho3[j] <- (j-1000)/1000
		SigmaCAR <- solve(diag(9)-rho3[j]*W)%*%K
		SigmaSAR <- solve(diag(9)-rho3[j]*W)%*%K%*%solve(diag(9)-rho3[j]*t(W))
		margcorrCAR3[i,j] = SigmaCAR[nborpairs3[i,1],nborpairs3[i,2]]/
			sqrt(SigmaCAR[nborpairs3[i,1],nborpairs3[i,1]]*
			SigmaCAR[nborpairs3[i,2],nborpairs3[i,2]])
		margcorrSAR3[i,j] = SigmaSAR[nborpairs3[i,1],nborpairs3[i,2]]/
			sqrt(SigmaSAR[nborpairs3[i,1],nborpairs3[i,1]]*
			SigmaSAR[nborpairs3[i,2],nborpairs3[i,2]])
	}
}

# Next for a 6x6 square lattice without row standardization
xy = expand.grid(x = seq(1, 6), y = seq(1, 6))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 1.01)*1
diag(W) = 0
K = diag(1/apply(W, 1, sum))
W = W/apply(W,1,sum)
# get first-order only indexes for pairwise correlations from 6 x 6 grid (36 locations)
nborpairs6 <- NULL
for(i in 1:35){
	for(j in (i + 1):36){
		if(W[i,j]!=0) nborpairs6 = rbind(nborpairs6, c(i,j))
	}
}
rho6 = rep(NA, times = 1999)
margcorrCAR6 = matrix(NA, dim(nborpairs6)[1], 1999)
margcorrSAR6 = matrix(NA, dim(nborpairs6)[1], 1999)
for(i in 1:dim(nborpairs6)[1]){
	for(j in 1:1999){
		rho6[j] <- (j-1000)/1000
		SigmaCAR <- solve(diag(36)-rho6[j]*W)%*%K
		SigmaSAR <- solve(diag(36)-rho6[j]*W)%*%K%*%solve(diag(36)-rho6[j]*t(W))
		margcorrCAR6[i,j] = SigmaCAR[nborpairs6[i,1],nborpairs6[i,2]]/
			sqrt(SigmaCAR[nborpairs6[i,1],nborpairs6[i,1]]*
			SigmaCAR[nborpairs6[i,2],nborpairs6[i,2]])
		margcorrSAR6[i,j] = SigmaSAR[nborpairs6[i,1],nborpairs6[i,2]]/
			sqrt(SigmaSAR[nborpairs6[i,1],nborpairs6[i,1]]*
			SigmaSAR[nborpairs6[i,2],nborpairs6[i,2]])
	}
}

file_name = "MargCorrRowStandWts"

pdf(paste0(file_name,'.pdf'), width = 11, height = 11)   

	layout(matrix(1:4, ncol = 2, byrow = TRUE))

	old.par = par(mar = c(5,5,5,1))
	plot(rho3, margcorrSAR3[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[SAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = '3 x 3 grid', cex.main = 2)
	for(i in 2:dim(nborpairs3)[1]){
		lines(rho3,margcorrSAR3[i,], lwd=1.0)}
	plot(rho3, margcorrCAR3[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[CAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = '3 x 3 grid', cex.main = 2)
	for(i in 2:dim(nborpairs3)[1]){
		lines(rho3,margcorrCAR3[i,], lwd=1.0)}
	plot(rho6, margcorrSAR6[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[SAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = '6 x 6 grid', cex.main = 2)
	for(i in 2:dim(nborpairs6)[1]){
		lines(rho6,margcorrSAR6[i,], lwd=1.0)}
	plot(rho6, margcorrCAR6[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[CAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = '6 x 6 grid', cex.main = 2)
	for(i in 2:dim(nborpairs6)[1]){
		lines(rho6,margcorrCAR6[i,], lwd=1.0)}
	par(old.par)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))


# -------------------- Spatial weights for US States ---------------------------

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
margcorrusaCAR = matrix(NA, dim(nborpairsusa)[1], 53)
margcorrusaSAR = matrix(NA, dim(nborpairsusa)[1], 53)
for(i in 1:dim(nborpairsusa)[1]){
	for(j in 1:53){
		rhousa[j] <- -0.08+(j-27)/100
		SigmaCAR <- solve(diag(49)-rhousa[j]*W)%*%K
		SigmaSAR <- solve(diag(49)-rhousa[j]*W)%*%K%*%solve(diag(49)-rhousa[j]*t(W))
		margcorrusaCAR[i,j] <- SigmaCAR[nborpairsusa[i,1],nborpairsusa[i,2]]/
			sqrt(SigmaCAR[nborpairsusa[i,1],nborpairsusa[i,1]]*
			SigmaCAR[nborpairsusa[i,2],nborpairsusa[i,2]])
		margcorrusaSAR[i,j] <- SigmaSAR[nborpairsusa[i,1],nborpairsusa[i,2]]/
			sqrt(SigmaSAR[nborpairsusa[i,1],nborpairsusa[i,1]]*
			SigmaSAR[nborpairsusa[i,2],nborpairsusa[i,2]])
	}
}
# range of correlations for rho for rowstandardized weights
rhousars = rep(NA, times = 239)
margcorrusaCARrs = matrix(NA, dim(nborpairsusa)[1], 239)
margcorrusaSARrs = matrix(NA, dim(nborpairsusa)[1], 239)
for(i in 1:dim(nborpairsusa)[1]){
	for(j in 1:239){
		rhousars[j] <- -0.2+(j-120)/100
		SigmaCAR <- solve(diag(49) - rhousars[j]*Wrs) %*% Krs
		SigmaSAR <- solve(diag(49) - rhousars[j]*Wrs) %*% Krs %*%
			solve(diag(49) - rhousars[j]*t(Wrs))
		margcorrusaCARrs[i,j] <- SigmaCAR[nborpairsusa[i,1],nborpairsusa[i,2]]/
			sqrt(SigmaCAR[nborpairsusa[i,1],nborpairsusa[i,1]]*
			SigmaCAR[nborpairsusa[i,2],nborpairsusa[i,2]])
		margcorrusaSARrs[i,j] <- SigmaSAR[nborpairsusa[i,1],nborpairsusa[i,2]]/
			sqrt(SigmaSAR[nborpairsusa[i,1],nborpairsusa[i,1]]*
			SigmaSAR[nborpairsusa[i,2],nborpairsusa[i,2]])
	}
}

file_name = "MargCorrUSA"

pdf(paste0(file_name,'.pdf'), width = 11, height = 11)   

	layout(matrix(1:4, ncol = 2, byrow = TRUE))

	old.par = par(mar = c(5,5,5,1))
	plot(rhousa, margcorrusaSAR[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[SAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = 'Binary Weights', cex.main = 2)
	for(i in 2:dim(nborpairsusa)[1]){
		lines(rhousa,margcorrusaSAR[i,], lwd=1.0)}
	plot(rhousa, margcorrusaCAR[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[CAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = 'Binary Weights', cex.main = 2)
	for(i in 2:dim(nborpairsusa)[1]){
		lines(rhousa,margcorrusaCAR[i,], lwd=1.0)}
	plot(rhousars, margcorrusaSARrs[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[SAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = 'Row-standardized Weights', cex.main = 2)
	for(i in 2:dim(nborpairsusa)[1]){
		lines(rhousars,margcorrusaSARrs[i,], lwd=1.0)}
	lines(c(-1,-1),c(-1,1), lty = 2, lwd = 2)
	plot(rhousars, margcorrusaCARrs[1,], type = 'l', ylim=c(-1.0,1.0), lwd = 1.0, 
		xlab = expression(paste(rho[CAR])),
		ylab = 'Marginal Correlation', cex.axis = 1.5, cex.lab = 2, 
		main = 'Row-standardized Weights', cex.main = 2)
	for(i in 2:dim(nborpairsusa)[1]){
		lines(rhousars,margcorrusaCARrs[i,], lwd=1.0)}
	lines(c(-1,-1),c(-1,1), lty = 2, lwd = 2)
	par(old.par)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
