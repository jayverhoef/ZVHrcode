sec_path = 'Rcode/Chapter11/Section 11.5/'
setwd(paste0(SLEDbook_path,sec_path))

library(gtools)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Average variance of GLS estimates of treatment differences for the
#         9-observation example of Section 10.6
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# function for generalized inverse
mginv <-function(X, tol = sqrt(.Machine$double.eps)) { 
	dnx <- dimnames(X) 
	if(is.null(dnx)) dnx <- vector("list", 2) 
	s <- svd(X) 
	nz <- s$d > tol * s$d[1] 
	structure( 
		if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X, 
	dimnames = dnx[2:1]) 
}

# all 6^3 orders of 3 treatments within 3 blocks
base_order3 = permutations(3,3,v=1:3)
all_orders = matrix(NA, nrow = 6^3, ncol = 9)
mm = 0
for(i in 1:6) {
	for(j in 1:6) {
		for(k in 1:6) {
			mm = mm + 1
			all_orders[mm, 1:3] = base_order3[i,]
			all_orders[mm, 4:6] = base_order3[j,]
			all_orders[mm, 7:9] = base_order3[k,]
		}
	}
}

# rho = 0.5 -------------------------

rho = 0.5

# 3 x 3 correlation matrix within blocks
cov=matrix(c(
	1, rho, rho^2,
	rho, 1, rho,
	rho^2, rho, 1), 
nrow = 3, byrow = TRUE)

# full correlation matrix that includes across blocks
V = kronecker(diag(3),cov)

# a matrix where each row is one of the 3 contrasts
con = rbind(c(0, 0, 0, 1, -1, 0),
		c(0, 0, 0, 1, 0, -1),
		c(0, 0, 0, 0, 1, -1))

# vector to store results
var_store.5 = rep(NA, times = 6^3)
for(i in 1:6^3) {
	# create a data.frame with ith treatment within block order
	df = data.frame(blk = as.factor(kronecker(1:3, rep(1, times = 3))),
		trt = as.factor(all_orders[i,]))
	# over-parameterized design matrix
	X = cbind(model.matrix(~ - 1 + blk, data = df),
		model.matrix(~ - 1 + trt, data = df))

	# generalized inverse of the covariance matrix of fixed effects
	ginvvar = mginv(t(X) %*% solve(V,X)) 

	# the sum of the 3 contrast variances
	var_store.5[i] = sum(diag(con %*% ginvvar %*% t(con)))
}

#var.5 = unique(round(var_store.5,4))/(9/4)
#var.5
var.5 = round(unique(round(var_store.5,4))/3,4)
var.5

# worst to best design percentage improvement
(var.5[1] - var.5[3])/var.5[1]*100

# average of whole population
mean(var_store.5)/3
(var.5[1] - mean(var_store.5)/3)/(mean(var_store.5)/3)*100

# rho = 0.8 -------------------------

rho = 0.8
var_store.8 = rep(NA, times = 6^3)
for(i in 1:6^3) {
	df = data.frame(blk = as.factor(kronecker(1:3, rep(1, times = 3))),
		trt = as.factor(all_orders[i,]))
	# over-parameterized design matrix
	X = cbind(model.matrix(~ - 1 + blk, data = df),
		model.matrix(~ - 1 + trt, data = df))

	cov=matrix(c(
		1, rho, rho^2,
		rho, 1, rho,
		rho^2, rho, 1), 
	nrow = 3, byrow = TRUE)

	V = kronecker(diag(3),cov)

	con = rbind(c(0, 0, 0, 1, -1, 0),
		c(0, 0, 0, 1, 0, -1),
		c(0, 0, 0, 0, 1, -1))

	ginvvar = mginv(t(X) %*% solve(V,X))

	var_store.8[i] = sum(diag(con %*% ginvvar %*% t(con)))
}

# var.8 = unique(round(var_store.8,4))/(9/4)
# var.8
var.8 = round(unique(round(var_store.8,4))/3,4)
var.8

# worst to best design percentage improvement
(var.8[1] - var.8[3])/var.8[1]*100

# average of whole population
mean(var_store.8)/3
(var.8[1] - mean(var_store.8)/3)/(mean(var_store.8)/3)*100


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Average variance of estimated treatment difference for a 
#       5x5 quasi-complete Latin square design
#          under a first-order CAR model
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# create a data.frame with ith treatment within row and column order
df = data.frame(
	trt = as.factor(c(1,2,3,4,5,
	                  2,4,1,5,3,
	                  3,1,5,2,4,
	                  4,5,2,3,1,
	                  5,3,4,1,2))
  )

# over-parameterized design matrix
X = model.matrix(~ - 1 + trt, data = df)

# easy way to create binary weights matrix based on rook's move in rectangular
# grid
nrook = (as.matrix(dist(cbind(kronecker(1:5, rep(1, times = 5)),
	kronecker(rep(1, times = 5), 1:5)))) < 1.1)*1
diag(nrook) = 0

# inverse covariance matrix for CAR model
Vinv=diag(25) - 0.25*nrook
# covariance matrix of fixed effects
ginvvar=mginv(t(X) %*% Vinv %*% X);
	
conMat=rbind(c(1,-1,0,0,0),
	c(1,0,-1,0,0),
	c(1,0,0,-1,0),
	c(1,0,0,0,-1),
	c(0,1,-1,0,0),
	c(0,1,0,-1,0),
	c(0,1,0,0,-1),
	c(0,0,1,-1,0),
	c(0,0,1,0,-1),
	c(0,0,0,1,-1))
	
# average variance for all contrasts
var_qLatin_CARrook = mean(diag(conMat %*% ginvvar %*% t(conMat)))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Code to compute loss in efficiency using CRD rather 
#       Latin square for a 5x5 arrangement of
#          plots under a 1st-order CAR model
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# simulate CRD designs to get average variance of contrasts
nsim = 10000
var_CRD_CARrook = rep(NA, times = nsim)
for(i in 1:nsim) {
	# create a data.frame with a random ordering of treatments
	df = data.frame(
		trt = as.factor(sample(kronecker(1:5, rep(1, times = 5)), 25)) )
	# over-parameterized design matrix
	X = model.matrix(~ - 1 + trt, data = df)
	# covariance matrix of fixed effects
	ginvvar=mginv(t(X) %*% Vinv %*% X);

	# average variance for all contrasts
	var_CRD_CARrook[i] = mean(diag(conMat %*% ginvvar %*% t(conMat)))
}
round(100*(mean(var_CRD_CARrook)/var_qLatin_CARrook - 1),1)
round((1 - var_qLatin_CARrook/mean(var_CRD_CARrook))*100,1)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#             Graph of Block Designs
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = 'figures/NNdesigns'
pdf(paste0(file_name,'.pdf'), width = 15, height = 7)

layout(matrix(c(1,2,2,3,3,3), nrow = 1), widths = c(1.2,1.2,1))

cex_mtext = 4
# A
plot(0:3, 0:3, type = 'n', xlab = '', xaxt = 'n', ylab ='', yaxt = 'n',
	bty = 'n', asp = 1)
for(i in 1:2) {
	segments(i,0,i,3, lwd = 2, lty = 2)
	segments(0,i,3,i, lwd = 6)
}
	i = 0
	segments(i,0,i,3, lwd = 2)
	segments(0,i,3,i, lwd = 2)
	i = 3
	segments(i,0,i,3, lwd = 2)
	segments(0,i,3,i, lwd = 2)
	x = kronecker(rep(1, times = 3), 1:3) - .5
	y = kronecker(3:1, rep(1, times = 3)) - .5
	labs = as.character(c(1,2,3,3,1,2,1,3,2))
	text(x, y, labels = labs, cex = 3)
	mtext('A', 	padj = 3.9, adj = 0, cex = cex_mtext)

# B
par(mar = c(0,12,0,0))
plot(0:5, 0:5, type = 'n', xlab = '', xaxt = 'n', ylab ='', yaxt = 'n',
	bty = 'n', asp = 1)
for(i in 1:4) {
	segments(i,0,i,5, lwd = 2, lty = 2)
	segments(0,i,5,i, lwd = 6)
}
	i = 0
	segments(i,0,i,5, lwd = 2)
	segments(0,i,5,i, lwd = 2)
	i = 5
	segments(i,0,i,5, lwd = 2)
	segments(0,i,5,i, lwd = 2)
	x = kronecker(rep(1, times = 5), 1:5) - .5
	y = kronecker(5:1, rep(1, times = 5)) - .5
	labs = as.character(c(
		1,2,3,4,5,
		2,4,1,5,3,
		3,1,5,2,4,
		4,5,2,3,1,
		5,3,4,1,2))
	text(x, y, labels = labs, cex = 3)
	mtext('B', 	padj = 3.5, adj = 0, cex = cex_mtext)

# C
par(mar = c(0,0,4,0))
plot(c(0,5), c(0,10), type = 'n', xlab = '', xaxt = 'n', ylab ='', yaxt = 'n',
	bty = 'n', asp = 1)
for(i in 1:9)
	segments(0,i,5,i, lwd = 6)
for(i in 1:4) 
	segments(i,0,i,10, lwd = 2, lty = 2)
	i = 0
	segments(i,0,i,10, lwd = 2)
	segments(0,i,5,i, lwd = 2)
	i = 5
	segments(i,0,i,10, lwd = 2)
	i = 10
	segments(0,i,5,i, lwd = 2)
	x = kronecker(rep(1, times = 10), 1:5) - .5
	y = kronecker(10:1, rep(1, times = 5)) - .5
	labs = as.character(c(
		5,4,1,3,2,
		2,5,4,1,3,
		3,2,5,4,1,
		1,3,2,5,4,
		4,1,3,2,5,
		5,1,2,4,3,
		3,5,1,2,4,
		4,3,5,1,2,
		2,4,3,5,1,
		1,2,4,3,5))
	text(x, y, labels = labs, cex = 3)
	mtext('C', 	padj = 0.2, adj = .29, cex = cex_mtext)

layout(1)

dev.off()
		
system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))



