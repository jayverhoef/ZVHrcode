sec_path = 'Rcode/Chapter9/Section 9.3/'
setwd(paste0(SLEDbook_path,sec_path))

library(spmodel)
library(xtable)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#          Table 9.3
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Code for universal kriging example (with planar mean surface), spherical model in Table 9.3
x = c(1, 1, 2, 2, 5)
y = c(4, 3, 3, 2, 4)
addx = c(2, 4, 4)
addy = c(4, 2, 4)
X <- cbind(rep(1, times = 5), x, y)
x01 = c(1, addx[1], addy[1])
x02 = c(1, addx[2], addy[2])
x03 = c(1, addx[3], addy[3])

distMat = as.matrix(dist(cbind(c(x,addx), c(y, addy))))
Rall = 1 - (3*distMat/8) + (distMat^3/128)
R = Rall[1:5, 1:5]
r01 = Rall[1:5, 6]
r02 = Rall[1:5, 7]
r03 = Rall[1:5, 8]

# function to compute universal kriging weights
ukwts = function(X,R,r0,x0){
	t(x0) %*% solve(t(X) %*% solve(R) %*%X) %*% t(X) %*% solve(R) +
	t(r0) %*% solve(R) %*% (diag(dim(X)[1]) - X %*% 
	solve(t(X) %*% solve(R) %*% X) %*% t(X) %*% solve(R))
}
	
# function to compute prediction variance for universal kriging
ukpev = function(X,R,r0,x0,r00){
	r00 - t(r0) %*% solve(R) %*% r0 + (t(x0) - t(r0) %*%
	solve(R) %*% X) %*% solve(t(X) %*% solve(R) %*% X) %*% 
	(x0 - t(X) %*% solve(R) %*% r0)
}

table_ukwts_pev = cbind(
	c(ukwts(X, R, r01, x01), ukpev(X, R, r01, x01, 1)),
	c(ukwts(X, R, r02, x02), ukpev(X, R, r02, x02, 1)),
	c(ukwts(X, R, r03, x03), ukpev(X, R, r03, x03, 1))
)


print(
    xtable(table_ukwts_pev, 
      align = c('l',rep('l', times = length(table_ukwts_pev[1,]))),
      digits = c(0, rep(3, times = 3))
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
#          Figure 9.3
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

xy <- data.frame( x = c(0.1, 0.15, 0.30, 0.32, 0.35, 0.40, 0.50, 
	0.60, 0.75, 0.90), y = c(rep(0, 10)))
new_xy <- data.frame(x = seq(0, 1, by = 0.001), y = rep(0, 1001))

# simulate data
set.seed(17)
spcov_params_val <- spcov_params("exponential", de = 1, ie = 0.0001, 
	range = 0.05)
simdata_exp1 = sprnorm(spcov_params_val, data = xy, xcoord = x, ycoord = y)

D1 = data.frame(ysim = simdata_exp1, x = xy$x, y = xy$y)
# create spmodel with fixed covariance parameters
splm_exp = splm(ysim ~ 1, data = D1, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
	de = 1, ie = 0.00001, range = 0.05, known = 'given'))
predexp = predict(splm_exp, new_xy, se.fit = TRUE, interval = 'prediction')

# universal kriging model with x coordinate as explanatory variable
splm_expUK = splm(ysim ~ x, data = D1, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
	de = 1, ie = 0.00001, range = 0.05, known = 'given'))
predexpUK = predict(splm_exp, new_xy, se.fit = TRUE, interval = 'prediction')

# universal kriging equations (will work for ordinary kriging too)
# X is design matrix
# XRiX_i = solve(t(X) %*% solve(Covmat) %*% X) [no need to solve for every 
#      prediction
# Ri = solve(Covmat) [no need to solve for every prediction]
# x0 = design vector for prediction site
# r0 = covariance between observed data and prediction site
# y = observed data
ukpred = function(X, XRiX_i, Ri, x0, r0, y){
	t(x0) %*% XRiX_i %*% t(X) %*% Ri %*% y +
	t(r0) %*% Ri %*% (diag(length(y)) - X %*% 
	XRiX_i %*% t(X) %*% Ri) %*% y
}
	
# universal kriging prediction variances
# X is design matrix
# XRiX_i = solve(t(X) %*% solve(Covmat) %*% X) [no need to solve for every 
#      prediction
# Ri = solve(Covmat) [no need to solve for every prediction]
# x0 = design vector for prediction site
# r0 = covariance between observed data and prediction site
# r00 = variance at prediction location
ukpev = function(X, XRiX_i, Ri, x0, r0, r00){
	r00 - t(r0) %*% Ri %*% r0 + (t(x0) - t(r0) %*%
	Ri %*% X) %*% XRiX_i %*% 
	(x0 - t(X) %*% Ri %*% r0)
}

# function to compute prediction variance for ordinary kriging
okpev = function(X,Ri,r0,r00){
	r00-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
}

# combine observed and prediction locations
xy_all = rbind(xy, new_xy)
# get distance matrix
distMat = as.matrix(dist(xy_all))
# create exponential covariance matrix
CovMat = exp(-distMat/0.05)
# create design matrix for observed locations
XOK = rep(1, times = 10)
# create design matrix for prediction locations
XpOK = rep(1, times = dim(new_xy)[1])
# subset the covariance matrix for only observed locations
Robs = CovMat[1:10, 1:10]
# inverse of covariance matrix for only observed locations
Ri = solve(Robs)
# solve(t(X) %*% solve(Covmat) %*% X) [no need to solve for every 
#      prediction
XRiX_i = solve(t(XOK) %*% Ri %*%XOK)
#create empty vectors to hold results
okpreds = rep(NA, times = dim(new_xy)[1])
okpevs = okpreds
for(i in 1:dim(new_xy)[1]) {
	okpreds[i] = ukpred(XOK, XRiX_i, Ri, XpOK[i], CovMat[1:10,10 + i], D1$ysim)
	okpevs[i] = ukpev(XOK, XRiX_i, Ri, XpOK[i], CovMat[1:10,10 + i], 
		CovMat[10 + i, 10 + i])
}

CovMat[1:10,111]

# create design matrix for observed locations
XUK = cbind(rep(1, times = 10), xy[,1])
# create design matrix for prediction locations
XpUK = cbind(rep(1, times = dim(new_xy)[1]), new_xy[,1])
# solve(t(X) %*% solve(Covmat) %*% X) [no need to solve for every 
#      prediction
XRiX_i = solve(t(XUK) %*% Ri %*%XUK)
#create empty vectors to hold results
ukpreds = rep(NA, times = dim(new_xy)[1])
ukpevs = ukpreds
i = 2
for(i in 1:dim(new_xy)[1]) {
	ukpreds[i] = ukpred(XUK, XRiX_i, Ri, XpUK[i,], CovMat[1:10,10 + i], D1$ysim)
	ukpevs[i] = ukpev(XUK, XRiX_i, Ri, XpUK[i,], CovMat[1:10,10 + i], 
		CovMat[10 + i, 10 + i])
}

# create exponential covariance matrix
CovMat = 0.5*exp(-distMat/0.05) + 0.5*diag(dim(distMat)[1])
# make sure that covariance is exactly 1 when prediction site is same
# as observed site
# which prediction sites are same as observed sites?
ind = which(round(new_xy[,1],3) %in% xy[,1])
# make sure that covariance is exactly 1 when prediction site is same
# as observed site
for(i in 1:10)
CovMat[i, 10 + ind[i]] = 1
# create design matrix for observed locations
XOK = rep(1, times = 10)
# create design matrix for prediction locations
XpOK = rep(1, times = dim(new_xy)[1])
# subset the covariance matrix for only observed locations
Robs = CovMat[1:10, 1:10]
# inverse of covariance matrix for only observed locations
Ri = solve(Robs)
# solve(t(X) %*% solve(Covmat) %*% X) [no need to solve for every 
#      prediction
XRiX_i = solve(t(XOK) %*% Ri %*%XOK)
#create empty vectors to hold results
okpreds50 = rep(NA, times = dim(new_xy)[1])
okpevs50 = okpreds
for(i in 1:dim(new_xy)[1]) {
	okpreds50[i] = ukpred(XOK, XRiX_i, Ri, XpOK[i], CovMat[1:10,10 + i], D1$ysim)
	okpevs50[i] = ukpev(XOK, XRiX_i, Ri, XpOK[i], CovMat[1:10,10 + i], 
		CovMat[10 + i, 10 + i])
}

# create exponential covariance matrix
CovMat = 0.5*exp(-distMat/0.05) + 0.5*diag(dim(distMat)[1])
# create design matrix for observed locations
XOK = rep(1, times = 10)
# create design matrix for prediction locations
XpOK = rep(1, times = dim(new_xy)[1])
# subset the covariance matrix for only observed locations
Robs = CovMat[1:10, 1:10]
# inverse of covariance matrix for only observed locations
Ri = solve(Robs)
# solve(t(X) %*% solve(Covmat) %*% X) [no need to solve for every 
#      prediction
XRiX_i = solve(t(XOK) %*% Ri %*%XOK)
#create empty vectors to hold results
okpredsNoiseless = rep(NA, times = dim(new_xy)[1])
okpevsNoiseless = okpredsNoiseless
for(i in 1:dim(new_xy)[1]) {
	okpredsNoiseless[i] = ukpred(XOK, XRiX_i, Ri, XpOK[i], 
		CovMat[1:10,10 + i], D1$ysim)
	okpevsNoiseless[i] = ukpev(XOK, XRiX_i, Ri, XpOK[i], 
		CovMat[1:10,10 + i], CovMat[10 + i, 10 + i])
}

file_name = 'figures/UK1DNoiseless'
pdf(paste0(file_name,'.pdf'), width = 12, height = 8)

	padj = -.4
	adj = -.2
	ylim = c(-1.8, 4.1)
	axat = -1:4
	leg_cex = 1.7
	cex.right = 1.8

	layout(matrix(1:4, nrow = 2, byrow = TRUE))

	par(mar = c(5,5,4,6))
	plot(new_xy[,1], okpreds, type = 'l', ylim = ylim,
		xlab = 'Location Coordinate', ylab = 'Response Value', 
		cex.lab = 2, cex.axis = 1.5, lwd = 2)
	axis(4, at = axat, cex.lab = 2, cex.axis = 1.5)
	mtext(side = 4, line = 3, "Variance", cex = 2)
	points(D1$x, D1$ysim, pch = 19, cex = 2)
	lines(c(0,1),c(0,0), lty = 2, lwd = 2)
	lines(new_xy[,1], okpevs, lty = 2, col = 'red', lwd = 2)
	legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
				 lty = c(1, 2), lwd = c(2,2), cex = leg_cex)
	mtext('A', adj = adj, cex = 3, padj = padj)

	plot(new_xy[,1], ukpreds, type = 'l', ylim = ylim,
		xlab = 'Location Coordinate', ylab = 'Response Value', 
		cex.lab = 2, cex.axis = 1.5, lwd = 2)
	axis(4, at = axat, cex.lab = 2, cex.axis = 1.5)
	mtext(side = 4, line = 3, "Variance", cex = cex.right)
	points(D1$x, D1$ysim, pch = 19, cex = 2)
	lines(c(0,1),c(0,0), lty = 2, lwd = 2)
	lines(new_xy[,1], ukpevs, lty = 2, col = 'red', lwd = 2)
	legend("topright", legend = c("UK", "UK var"), col = c("black", "red"),
				 lty = c(1, 2), lwd = c(2,2), cex = leg_cex)
	mtext('B', adj = adj, cex = 3, padj = padj)

	plot(new_xy[,1], okpreds50, type = 'l', ylim = ylim,
		xlab = 'Location Coordinate', ylab = 'Response Value', 
		cex.lab = 2, cex.axis = 1.5, lwd = 2)
	axis(4, at = axat, cex.lab = 2, cex.axis = 1.5)
	mtext(side = 4, line = 3, "Variance", cex = cex.right)
	points(D1$x, D1$ysim, pch = 19, cex = 2)
	lines(c(0,1),c(0,0), lty = 2, lwd = 2)
	lines(new_xy[,1], okpevs50, lty = 2, col = 'red', lwd = 2)
	legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
				 lty = c(1, 2), lwd = c(2,2), cex = leg_cex)
	mtext('C', adj = adj, cex = 3, padj = padj)

	plot(new_xy[,1], okpredsNoiseless, type = 'l', ylim = ylim,
		xlab = 'Location Coordinate', ylab = 'Response Value', 
		cex.lab = 2, cex.axis = 1.5, lwd = 2)
	axis(4, at = axat, cex.lab = 2, cex.axis = 1.5)
	mtext(side = 4, line = 3, "Variance", cex = cex.right)
	points(D1$x, D1$ysim, pch = 19, cex = 2)
	lines(c(0,1),c(0,0), lty = 2, lwd = 2)
	lines(new_xy[,1], okpevsNoiseless, lty = 2, col = 'red', lwd = 2)
	legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
				 lty = c(1, 2), lwd = c(2,2), cex = leg_cex)
	mtext('D', adj = adj, cex = 3, padj = padj)

	layout(1)

dev.off()
		
system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

