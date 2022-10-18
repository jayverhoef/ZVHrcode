sec_path = 'Rcode/Chapter9/Section 9.2/'
setwd(paste0(SLEDbook_path,sec_path))

library(spmodel)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#          Figure Spatial Configuration of Toy Example
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
file_name = 'figures/figure9pt1'
pdf(paste0(file_name,'.pdf'), width = 6, height = 6)

	# Code for creating Figure 9.1
	x = c(1,1,2,2,5)
	y = c(4,3,3,2,4)
	addx = c(2,4,4)
	addy = c(4,2,4)
	par(mar = c(5,5,1,1))
	plot(x, y, pch=16, xlim=c(1,5), ylim=c(1,5), cex=3,
		xlab = 'x-coordinate', ylab = 'y-coordinate', cex.lab = 2, 
		cex.axis = 1.5)
	points(addx, addy, cex=3, lwd = 3)
	text(x, y, labels = c(expression(bold(s)[1]), expression(bold(s)[2]),
		expression(bold(s)[3]),	expression(bold(s)[4]), expression(bold(s)[5])),
		cex = 2, pos = 3, offset = 1)
	text(addx, addy, labels = c(expression(bold(s)['01']), 
		expression(bold(s)['02']), expression(bold(s)['03'])), 
		cex = 2, pos = 3, offset = 1)

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
#          Matrices
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

sph_cor = function(distmat) {
	(1 - 3*distmat/8 + distmat^3/128)*(distmat < 4)
}

distMat = as.matrix(dist(cbind(c(x,addx), c(y, addy))))

covMat_all = sph_cor(distMat)
covMat_y = covMat_all[1:5,1:5]
blank = matrix(NA, nrow = 5, ncol = 5)
blank[upper.tri(covMat_y, diag = TRUE)] = 
	covMat_y[upper.tri(covMat_y, diag = TRUE)]
covMat_y_print = blank

library(xtable)

# correlation among all observed locations
print(
    xtable(covMat_y_print, 
      align = c('l',rep('l', times = length(covMat_y_print[1,]))),
      digits = c(0, rep(3, times = 5))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# correlation between observed locations and 1st prediction location
print(
    xtable(matrix(covMat_all[1:5,6], ncol = 1), 
      align = c('l','l'),
      digits = c(0, 3)
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# correlation between observed locations and 1st prediction location
print(
    xtable(matrix(covMat_all[1:5,7], ncol = 1), 
      align = c('l','l'),
      digits = c(0, 3)
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# correlation between observed locations and 1st prediction location
print(
    xtable(matrix(covMat_all[1:5,8], ncol = 1), 
      align = c('l','l'),
      digits = c(0, 3)
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
#          Table 9.1
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# function to compute ordinary kriging weights
okwts = function(X,R,r0){
	solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(dim(X)[1])-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
}
	
# function to compute prediction variance for ordinary kriging
okpev = function(X,R,r0,r00){
	r00-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
}

X <- matrix(1,5,1)

sph_wtsvar = cbind(
	c(
		okwts(X, covMat_y, covMat_all[1:5,6]),
		okpev(X, covMat_y, covMat_all[1:5,6], 1)),
	c(
		okwts(X, covMat_y, covMat_all[1:5,7]),
		okpev(X, covMat_y, covMat_all[1:5,7], 1)),
	c(
		okwts(X, covMat_y, covMat_all[1:5,8]),
		okpev(X, covMat_y, covMat_all[1:5,8], 1))
)


exp_cor = function(distmat) {
	exp(-3*distmat/4)
}
expMat_all = exp_cor(distMat)
expMat_y = expMat_all[1:5,1:5]

exp_wtsvar = cbind(
	c(
		okwts(X, expMat_y, expMat_all[1:5,6]),
		okpev(X, expMat_y, expMat_all[1:5,6], 1)),
	c(
		okwts(X, expMat_y, expMat_all[1:5,7]),
		okpev(X, expMat_y, expMat_all[1:5,7], 1)),
	c(
		okwts(X, expMat_y, expMat_all[1:5,8]),
		okpev(X, expMat_y, expMat_all[1:5,8], 1))
)

gau_cor = function(distmat) {
	exp(-3*distmat^2/16)
}
gauMat_all = gau_cor(distMat)
gauMat_y = gauMat_all[1:5,1:5]

gau_wtsvar = cbind(
	c(
		okwts(X, gauMat_y, gauMat_all[1:5,6]),
		okpev(X, gauMat_y, gauMat_all[1:5,6], 1)),
	c(
		okwts(X, gauMat_y, gauMat_all[1:5,7]),
		okpev(X, gauMat_y, gauMat_all[1:5,7], 1)),
	c(
		okwts(X, gauMat_y, gauMat_all[1:5,8]),
		okpev(X, gauMat_y, gauMat_all[1:5,8], 1))
)


table_wts_pev = cbind(sph_wtsvar, exp_wtsvar, gau_wtsvar)

print(
    xtable(table_wts_pev, 
      align = c('l',rep('l', times = length(table_wts_pev[1,]))),
      digits = c(0, rep(3, times = 9))
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
#          Figure 9.2
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

xy <- data.frame( x = c(0.1, 0.15, 0.30, 0.32, 0.35, 0.40, 0.50, 
	0.60, 0.75, 0.90), y = c(rep(0, 10)))
new_xy <- data.frame(x = seq(0, 1, by = 0.001), y = rep(0, 1001))

# simulate data
set.seed(13)
spcov_params_val <- spcov_params("exponential", de = 1, ie = 0.0001, 
	range = 0.05)
simdata_exp1 = sprnorm(spcov_params_val, data = xy, xcoord = x, ycoord = y)

D1 = data.frame(ysim = simdata_exp1, x = xy$x, y = xy$y)
# create spmodel with fixed covariance parameters
splm_exp = splm(ysim ~ 1, data = D1, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
	de = 1, ie = 0.00001, range = 0.05, known = 'given'))
splm_sph = splm(ysim ~ 1, data = D1, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'spherical', 
	de = 1, ie = 0.00001, range = 0.15, known = 'given'))
splm_gau = splm(ysim ~ 1, data = D1, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'gaussian', 
	de = 1, ie = 0.00001, range = 0.15/sqrt(3), known = 'given'))
predexp = predict(splm_exp, new_xy, se.fit = TRUE, interval = 'prediction')
predsph = predict(splm_sph, new_xy, se.fit = TRUE, interval = 'prediction')
predgau = predict(splm_gau, new_xy, se.fit = TRUE, interval = 'prediction')

file_name = 'figures/OK1DwithVar'
pdf(paste0(file_name,'.pdf'), width = 12, height = 8)

	padj = -.4
	adj = -.2
	ylim = c(-1,4)
	axat = -1:4
	leg_cex = 1.7
	layout(matrix(1:4, nrow = 2, byrow = TRUE))
	par(mar = c(5,5,4,6))
	plot(D1$x, D1$ysim, pch = 19, ylim = ylim, cex = 2,
		xlab = 'Location Coordinate', ylab = 'Response Value', 
		cex.lab = 2, cex.axis = 1.5)
	lines(c(0,1),c(0,0), lty = 2, lwd = 2)
	mtext('A', adj = adj, cex = 3, padj = padj)

	plot(new_xy[,1], predexp$fit[,1], type = 'l', ylim = ylim,
		xlab = 'Location Coordinate', ylab = 'Response Value', 
		cex.lab = 2, cex.axis = 1.5, lwd = 2)
	axis(4, at = axat, cex.lab = 2, cex.axis = 1.5)
	mtext(side = 4, line = 3, "Variance", cex = 2)
	points(D1$x, D1$ysim, pch = 19, cex = 2)
	lines(c(0,1),c(0,0), lty = 2, lwd = 2)
	lines(new_xy[,1], predexp$se.fit^2, lty = 2, col = 'red', lwd = 2)
	legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
				 lty = c(1, 2), lwd = c(2,2), cex = leg_cex)
	mtext('B', adj = adj, cex = 3, padj = padj)

	plot(new_xy[,1], predsph$fit[,1], type = 'l', ylim = ylim,
		xlab = 'Location Coordinate', ylab = 'Response Value', 
		cex.lab = 2, cex.axis = 1.5, lwd = 2)
	axis(4, at = axat, cex.lab = 2, cex.axis = 1.5)
	mtext(side = 4, line = 3, "Variance", cex = 2)
	points(D1$x, D1$ysim, pch = 19, cex = 2)
	lines(c(0,1),c(0,0), lty = 2, lwd = 2)
	lines(new_xy[,1], predsph$se.fit^2, lty = 2, col = 'red', lwd = 2)
	legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
				 lty = c(1, 2), lwd = c(2,2), cex = leg_cex)
	mtext('C', adj = adj, cex = 3, padj = padj)


	plot(new_xy[,1], predgau$fit[,1], type = 'l', ylim = ylim,
		xlab = 'Location Coordinate', ylab = 'Response Value', 
		cex.lab = 2, cex.axis = 1.5, lwd = 2)
	axis(4, at = axat, cex.lab = 2, cex.axis = 1.5)
	mtext(side = 4, line = 3, "Variance", cex = 2)
	points(D1$x, D1$ysim, pch = 19, cex = 2)
	lines(c(0,1),c(0,0), lty = 2, lwd = 2)
	lines(new_xy[,1], predgau$se.fit^2, lty = 2, col = 'red', lwd = 2)
	legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
				 lty = c(1, 2), lwd = c(2,2), cex = leg_cex)
	mtext('D', adj = adj, cex = 3, padj = padj)

	layout(1)

dev.off()

labels = c(expression(y[1]), expression(y[2]), expression(y[3]),
		expression(y[4]), expression(y[5]))
		
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
#          Table 9.2
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# Code for ordinary kriging example in the context of a spatial-weights linear model (Table 9.2)
xcoords = c(2, 1, 2, 3, 1, 2, 3)
ycoords = c(3, 2, 2, 2, 1, 1, 1)
W = (as.matrix(dist(cbind(xcoords, ycoords))) < 1.1)*1
diag(W) = 0

K <- diag(7)
SigmaSAR <- solve(diag(7)-0.3*W)%*%K%*%solve(diag(7)-0.3*t(W))
SigmaCAR <- solve(diag(7)-0.3*W)%*%K
corrSAR <- solve(diag(sqrt(diag(SigmaSAR))))%*%SigmaSAR%*%solve(diag(sqrt(diag(SigmaSAR))))
corrCAR <- solve(diag(sqrt(diag(SigmaCAR))))%*%SigmaCAR%*%solve(diag(sqrt(diag(SigmaCAR))))
R <- SigmaSAR[-1,-1]
r0 <- as.vector(SigmaSAR[2:7,1])
r00 <- SigmaSAR[1,1]
X <- matrix(1,nrow=6,ncol=1,byrow=T)
S1 = rep(NA, times = 8)
S1[2:7] = okwts(X, R, r0)
S1[8] = okpev(X, R, r0, r00)

R <- SigmaSAR[-3,-3]
r0 <- SigmaSAR[-3,3]
r00 <- SigmaSAR[3,3]
X <- matrix(1,nrow=6,ncol=1,byrow=T)
S2 = rep(NA, times = 8)
S2[c(1,2,4:7)] = okwts(X, R, r0)
S2[8] = okpev(X, R, r0, r00)

R <- SigmaSAR[-7,-7]
r0 <- SigmaSAR[-7,7]
r00 <- SigmaSAR[7,7]
X <- matrix(1,nrow=6,ncol=1,byrow=T)
S3 = rep(NA, times = 8)
S3[1:6] = okwts(X, R, r0)
S3[8] = okpev(X, R, r0, r00)

SARwtspev = cbind(S1, S2, S3)

print(
    xtable(SARwtspev, 
      align = c('l',rep('l', times = length(SARwtspev[1,]))),
      digits = c(0, rep(3, times = 3))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

