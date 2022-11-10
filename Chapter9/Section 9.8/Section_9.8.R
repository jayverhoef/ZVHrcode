sec_path = 'Rcode/Chapter9/Section 9.8/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Figure to Illustrate Screening Effect
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Code for plotting Figure 9.4
fixedsite <- matrix(c(1.1, 0), nrow = 1, ncol = 2, byrow = T)
predsite <- matrix(c(0, 0), nrow = 1, ncol = 2, byrow = T)
transient <- matrix(c(-1, 0, -sqrt(2)/2, sqrt(2)/2, 0, 1, sqrt(2)/2,
	sqrt(2)/2, cos(pi/8), sin(pi/8), cos(pi/16), sin(pi/16), cos(pi/32),
	sin(pi/32), 1, 0), nrow=8, ncol=2, byrow = T)

file_name = 'figures/sites4screening_effect'

pdf(paste0(file_name,'.pdf'), width = 7.2, height = 4.5)

	par(mar = c(5,5,1,1))
	plot(transient[,1], transient[,2], xlim=c(-1.2,1.2), ylim = c(-.2, 1.2),
		xlab = "", ylab = "", pch = 19, cex = 2, cex.axis = 1.4)
	text(transient[,1]-0.08, transient[,2], labels=c(1:8), cex = 1.5)
	points(fixedsite, pch = 1, cex = 2, lwd = 2)
	points(predsite, pch = 4, cex = 2, lwd = 3)

dev.off()
		
system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Table for Screening Effect
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Some useful functions


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

# universal kriging equations (will work for ordinary kriging too)
# X is design matrix
# XRiX_i = solve(t(X) %*% solve(Covmat) %*% X) [no need to solve for every 
#      prediction
# Ri = solve(Covmat) [no need to solve for every prediction]
# x0 = design vector for prediction site
# r0 = covariance between observed data and prediction site
ukpredwts = function(X, XRiX_i, Ri, x0, r0){
	t(x0) %*% XRiX_i %*% t(X) %*% Ri +
	t(r0) %*% Ri %*% (diag(dim(R)[1]) - X %*% 
	XRiX_i %*% t(X) %*% Ri)
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


# Code for illustrating the screening effect, Gaussian model

# get all distances, transients first, then fixed observed site, then
# prediction site
locs = rbind(transient, c(1.1,0), c(0, 0))
distMat = as.matrix(dist(locs))
X <- rep(1, times = 2)

# spherical model
okwts_sph = NULL
okpev_sph = NULL
for(i in 1:8) {
	Rall = 1 - 1.5*distMat[c(i, 9, 10),c(i, 9, 10)]/4 + 
		0.5*(distMat[c(i, 9, 10),c(i, 9, 10)]/4)^3
	# observed locations
	R = Rall[1:2, 1:2]
	# prediction location
	r0 = Rall[1:2, 3]
	# kriging weights
	okwts_sph = c(okwts_sph, 
		ukpredwts(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0)[2])
	# prediction variance
	okpev_sph = c(okpev_sph,
		ukpev(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0, 1))
}

# exponential model
okwts_exp = NULL
okpev_exp = NULL
for(i in 1:8) {
	Rall = exp(-3*distMat[c(i, 9, 10),c(i, 9, 10)]/4)
	# observed locations
	R = Rall[1:2, 1:2]
	# prediction location
	r0 = Rall[1:2, 3]
	# kriging weights
	okwts_exp = c(okwts_exp, 
		ukpredwts(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0)[2])
	# prediction variance
	okpev_exp = c(okpev_exp,
		ukpev(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0, 1))
}

# Gaussian model
okwts_gau = NULL
okpev_gau = NULL
for(i in 1:8) {
	Rall = exp(-3*distMat[c(i, 9, 10),c(i, 9, 10)]^2/16)
	# observed locations
	R = Rall[1:2, 1:2]
	# prediction location
	r0 = Rall[1:2, 3]
	# kriging weights
	okwts_gau = c(okwts_gau, 
		ukpredwts(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0)[2])
	# prediction variance
	okpev_gau = c(okpev_gau,
		ukpev(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0, 1))
}

# create the table
screening_wts = cbind(
	okwts_sph, okpev_sph,
	okwts_exp, okpev_exp,
	okwts_gau, okpev_gau
)

print(
    xtable(screening_wts, 
      align = c('l',rep('l', times = length(screening_wts[1,]))),
      digits = c(0, rep(3, times = 6))
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
#           Illustrating the Screening Effect (shadow)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

tran_locs = data.frame(x = kronecker(seq(1, 3, by = 0.01), rep(1, times = 301)),
	y = kronecker(rep(1, times = 201), seq(-1.5, 1.5, by = 0.01)))
ntran = dim(tran_locs)[1]
X <- rep(1, times = 2)


gau_wts = matrix(NA, nrow = ntran, ncol = 2)
i = 1
for(i in 1:ntran) {
	distMat = as.matrix(dist(rbind(tran_locs[i,], c(1,0),c(0,0))))
	Rall = exp(-distMat^2/3) + diag(rep(1e-8, times = 3))
	# observed locations
	R = Rall[1:2, 1:2]
	# prediction location
	r0 = Rall[1:2, 3]
	# kriging weights
	gau_wts[i,] =  ukpredwts(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0)
}

file_name = 'figures/screening_shadow'

pdf(paste0(file_name,'.pdf'), width = 9, height = 9)

	par(mar = c(5,5,1,1))
	plot(0, 0, pch = 19, type = 'n', xlim = c(0,3), ylim = c(-1.5, 1.5),
		xlab = '', ylab = '', cex.axis = 1.5)
	grid(NULL, NULL, col = "gray70", lty = "dotted", lwd = 2)
	points(tran_locs[gau_wts[,1] < 0,], col = 'grey70', pch = 19, cex = .6)
	points(0, 0, pch = 1, cex = 2.5, lwd = 3)
	points(1, 0, pch = 19, cex = 2.5)
	text(0.11,-.12, label = expression(bold(s[0])), cex = 2.5)
	text(0.91,-.12, label = expression(bold(s[1])), cex = 2.5)

dev.off()
		
system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

png(paste0(file_name,'.png'), width = 640, height = 640)

	par(mar = c(5,5,1,1))
	plot(0, 0, pch = 19, type = 'n', xlim = c(0,3), ylim = c(-1.5, 1.5),
		xlab = '', ylab = '', cex.axis = 1.5)
	grid(NULL, NULL, col = "gray70", lty = "dotted", lwd = 2)
	points(tran_locs[gau_wts[,1] < 0,], col = 'grey70', pch = 19, cex = .6)
	points(0, 0, pch = 1, cex = 2.5, lwd = 3)
	points(1, 0, pch = 19, cex = 2.5)
	text(0.11,-.12, label = expression(bold(s[0])), cex = 2.5)
	text(0.91,-.12, label = expression(bold(s[1])), cex = 2.5)

dev.off()

		
