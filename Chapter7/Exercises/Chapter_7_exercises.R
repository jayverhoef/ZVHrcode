sec_path = 'Rcode/Chapter7/Section 7.7/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)

################################################################################
#-------------------------------------------------------------------------------
#                         Exercise 7.10a
#-------------------------------------------------------------------------------
################################################################################

# Investigate behavior of SAR model's marginal correlations for rho_SAR
# outside the interval within which I-rho_CAR*W is positive definite.

W0 <- matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,byrow=T)
id3 <- diag(3)
ze3 <- matrix(0,3,3)
W1 <- cbind(W0,id3,ze3)
W2 <- cbind(id3,W0,id3)
W3 <- cbind(ze3,id3,W0)
W <- rbind(W1,W2,W3)
K <- diag(9)
1/min(eigen(W)$values)
1/max(eigen(W)$values)
# reciprocals of smallest and largest eigenvalues of this W are +-0.35355
nborpairs <- NULL
for(i in 1:8){
	for(j in i:9){
		if(W[i,j]!=0) nborpairs = rbind(nborpairs, c(i,j))
	}
}
rho = rep(NA, times = 21)
margcorrSAR = matrix(NA, dim(nborpairs)[1], 21)
for(i in 1:dim(nborpairs)[1]){
	for(j in 1:21){
		rho[j] <- (j-11)/10
		SigmaSAR <- solve(diag(9) - rho[j]*W) %*% K %*%
			solve(diag(9) - rho[j]*t(W))
		margcorrSAR[i,j] <- SigmaSAR[nborpairs[i,1],nborpairs[i,2]]/
			sqrt(SigmaSAR[nborpairs[i,1],nborpairs[i,1]]*
			SigmaSAR[nborpairs[i,2],nborpairs[i,2]])
	}
}

old.par = par(mar = c(5,5,1,1))
plot(rho, margcorrSAR[1,], ylim=c(-1.0,1.0), type="l", lwd=0.5,
	ylab="marginal correlation", xlab=expression(paste(rho[SAR])),
	cex.lab = 2, cex.axis = 1.5)
for(i in 2:12)
	lines(rho,margcorrSAR[i,], lwd=0.5)
par(old.par)

# Then for a 6x6 square lattice
xy = expand.grid(x = seq(1, 6), y = seq(1, 6))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 1.01)*1
diag(W) = 0
K <- diag(36)
# reciprocals of smallest and largest eigenvalues of W 
1/min(eigen(W)$values)
1/max(eigen(W)$values)
# reciprocals of smallest and largest eigenvalues of this W are +-0.27748
nborpairs6 <- NULL
for(i in 1:35){
	for(j in (i + 1):36){
		if(W[i,j]!=0) nborpairs6 = rbind(nborpairs6, c(i,j))
	}
}
rho6 = rep(NA, times = 21)
margcorrSAR6 <- matrix(0,60,21)
for(i in 1:dim(nborpairs6)[1]){
	for(j in 1:21){
		rho6[j] <- (j-11)/10
		SigmaSAR <- solve(diag(36) - rho6[j]*W) %*% K %*%
			solve(diag(36) - rho6[j]*t(W))
		margcorrSAR6[i,j] <- SigmaSAR[nborpairs6[i,1],nborpairs6[i,2]]/
			sqrt(SigmaSAR[nborpairs6[i,1],nborpairs6[i,1]]*
			SigmaSAR[nborpairs6[i,2],nborpairs6[i,2]])
	}
}

old.par = par(mar = c(5,5,1,1))
plot(rho6,margcorrSAR6[1,], ylim=c(-1.0,1.0), type="l", lwd=0.5, 
	ylab="Marginal Correlation", xlab=expression(paste(rho[SAR])),
	cex.lab = 2, cex.axis = 1.5)
for(i in 2:60)
	lines(rho6,margcorrSAR6[i,], lwd=0.5)
par(old.par)

################################################################################
#-------------------------------------------------------------------------------
#                         Exercise 7.10b
#-------------------------------------------------------------------------------
################################################################################

# Now, for Exercise 7.10b, do something similar for the unconstrained CAR model:
# First for the 3x3 lattice
xy = expand.grid(x = seq(1, 3), y = seq(1, 3))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 1.01)*1
diag(W) = 0
margcorrCARpet <- matrix(0,12,401)
C <- matrix(0,9,9)
phi <- rep(NA, times = 401)
rowsum <- apply(W, 1, sum)
D <- diag(as.vector(rowsum))
for(i in 1:dim(nborpairs)[1]){
	for(k in 1:401){
		phi[k] <- (k - 201)/2
		Kinv <- diag(9) + abs(phi[k])*D
		K <- solve(Kinv)
		for(l in 1:9){
			for(m in 1:9){
				C[l,m] <- phi[k]*W[l,m]/(1 + abs(phi[k])*rowsum[l])
			}
		}
		SigmaCAR <- solve(diag(9) - C) %*% K
		margcorrCARpet[i,k] <- SigmaCAR[nborpairs[i,1],nborpairs[i,2]]/
			sqrt(SigmaCAR[nborpairs[i,1],nborpairs[i,1]]*SigmaCAR[nborpairs[i,2],nborpairs[i,2]])
	}
}

old.par = par(mar = c(5,5,1,1))
plot(phi, margcorrCARpet[1,], xlim=c(min(phi),max(phi)), ylim=c(-1,1), type="l", 
	lwd=0.5,ylab="Marginal Correlation",xlab="phi")
for(i in 2:12)
	lines(phi, margcorrCARpet[i,], lwd=0.5)
par(old.par)

# Then for 6x6 lattice
xy = expand.grid(x = seq(1, 6), y = seq(1, 6))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 1.01)*1
diag(W) = 0
margcorrCARpet6 <- matrix(0,dim(nborpairs6)[1],401)
C <- matrix(0,36,36)
phi6 <- rep(NA, times = 401)
rowsum <- apply(W, 1, sum)
D <- diag(as.vector(rowsum))
for(i in 1:dim(nborpairs6)[1]){
	for(k in 1:401){
		phi6[k] <- (k - 201)/2
		Kinv <- diag(36) + abs(phi6[k])*D
		K <- solve(Kinv)
		for(l in 1:36){
			for(m in 1:36){
				C[l,m] <- phi6[k]*W[l,m]/(1 + abs(phi6[k])*rowsum[l])
			}
		}
		SigmaCAR <- solve(diag(36) - C) %*% K
		margcorrCARpet6[i,k] <- SigmaCAR[nborpairs6[i,1],nborpairs6[i,2]]/
			sqrt(SigmaCAR[nborpairs6[i,1],nborpairs6[i,1]]*SigmaCAR[nborpairs6[i,2],nborpairs6[i,2]])
	}
}

old.par = par(mar = c(5,5,1,1))
plot(phi6, margcorrCARpet6[1,], xlim=c(min(phi),max(phi)), ylim=c(-1,1), type="l", 
	lwd=0.5,ylab="Marginal Correlation",xlab="phi")
for(i in 2:dim(nborpairs6)[1])
	lines(phi6, margcorrCARpet6[i,], lwd=0.5)
par(old.par)


################################################################################
#-------------------------------------------------------------------------------
#                         Exercise 7.10c
#-------------------------------------------------------------------------------
################################################################################

# Now, for Exercise 7.10c, do something similar for the "autocorrelated CAR model":
# First for 3x3 lattice
xy = expand.grid(x = seq(1, 3), y = seq(1, 3))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 1.01)*1
diag(W) = 0
rowsum <- apply(W,1,sum)
K <- diag(1/apply(W,1,sum))
Wac <- matrix(0,9,9)
for(i in 1:9){
	for(j in 1:9){
		Wac[i,j] <- sqrt(rowsum[j])/sqrt(rowsum[i])
	}
}
Wac <- Wac*W
# reciprocals of smallest and largest eigenvalues of W 
1/min(eigen(Wac)$values)
1/max(eigen(Wac)$values)

rhoac = rep(NA, times = 707)
margcorrCARac <- matrix(NA, dim(nborpairs)[1], 707)
for(i in 1:dim(nborpairs)[1]){
 for(j in 1:707){
  rhoac[j] <- (j-354)/1000
  SigmaCARac <- solve(diag(9) - rhoac[j]*Wac) %*% K
  margcorrCARac[i,j] <- SigmaCARac[nborpairs[i,1],nborpairs[i,2]]/
		sqrt(SigmaCARac[nborpairs[i,1],nborpairs[i,1]]*
		SigmaCARac[nborpairs[i,2],nborpairs[i,2]])
	}
}
old.par = par(mar = c(5,5,1,1))
plot(rhoac,margcorrCARac[1,], ylim=c(-1.0,1.0), type="l", lwd=0.5, 
	ylab="Marginal Correlation", xlab=expression(paste(rho[CAR])),
	cex.lab = 2, cex.axis = 1.5)
for(i in 2:12)
	lines(rhoac,margcorrCARac[i,], lwd=0.5)
par(old.par)

# Then for a 6x6 lattice:
xy = expand.grid(x = seq(1, 6), y = seq(1, 6))
# create distance matrix among spatial locations
Distmat = as.matrix(dist(xy))
# create first-order neighbor matrix (rook's move) from distances
W = (Distmat < 1.01)*1
diag(W) = 0
rowsum <- apply(W,1,sum)
K <- diag(1/apply(W,1,sum))
Wac <- matrix(0,36,36)
for(i in 1:36){
	for(j in 1:36){
		Wac[i,j] <- sqrt(rowsum[j])/sqrt(rowsum[i])
	}
}
Wac <- Wac*W
# reciprocals of smallest and largest eigenvalues of W 
1/min(Re(eigen(Wac)$values))
1/max(Re(eigen(Wac)$values))

rhoac6 = rep(NA, times = 555)
margcorrCARac6 <- matrix(0,60,555)
for(i in 1:dim(nborpairs6)[1]){
	for(j in 1:555){
		rhoac6[j] = (j - 278)/1000
		SigmaCAR <- solve(diag(36) - rhoac6[j]*Wac) %*% K
		margcorrCARac6[i,j] <- SigmaCAR[nborpairs6[i,1],nborpairs6[i,2]]/
			sqrt(SigmaCAR[nborpairs6[i,1],nborpairs6[i,1]]*
			SigmaCAR[nborpairs6[i,2],nborpairs6[i,2]])
	}
}

old.par = par(mar = c(5,5,1,1))
plot(rhoac6,margcorrCARac6[1,], ylim=c(-1.0,1.0), type="l", lwd=0.5, 
	ylab="Marginal Correlation", xlab=expression(paste(rho[CAR])),
	cex.lab = 2, cex.axis = 1.5)
for(i in 2:12)
	lines(rhoac6,margcorrCARac6[i,], lwd=0.5)
par(old.par)

