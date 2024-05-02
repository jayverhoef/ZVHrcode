#set a path as a working directory
sec_path = 'Rcode/Chapter8/Section 8.3/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                    graph explaining asymptotics 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# create graphic of two different types of asymptotics

temp1 <- seq(0, 1/3, length.out = 6)[2:5]
x1 <- rep(temp1, 4)
y1 <- rep(temp1, rep(4, 4))


temp2 <- seq(0, 1/2, length.out = 9)[2:8]
x2 <- rep(temp2, 7)
y2 <- rep(temp2, rep(7, 7))


temp3 <- seq(0, 1, length.out = 15)[2:14]
x3 <- rep(temp3, 13)
y3 <- rep(temp3, rep(13, 13))

temp4 <- seq(0+0.05, 1-0.05, length.out = 4)
x4 <- rep(temp4, 4)
y4 <- rep(temp4, rep(4, 4))


temp5 <- seq(0+0.05, 1-0.05, length.out = 7)
x5 <- rep(temp5, 7)
y5 <- rep(temp5, rep(7, 7))


temp6 <- seq(0, 1, length.out = 15)[2:14]
x6 <- rep(temp3, 13)
y6 <- rep(temp3, rep(13, 13))

file_name = "figures/asymptotics"
pdf(file = paste0(file_name,'.pdf'), width = 9, height = 6)

	adj = -.15
	padj = -.15
	layout(matrix(1:6, nrow = 2, byrow = TRUE)) 

		par(mar = c(1,4,4,1))
		plot(x1, y1, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
				 xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
		mtext('A', adj = adj, padj = padj, cex = 3)

		plot(x2, y2, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
				 xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
		mtext('B', adj = adj, padj = padj, cex = 3)

		plot(x3, y3, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
				 xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
		mtext('C', adj = adj, padj = padj, cex = 3)

		plot(x4, y4, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
				 xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
		mtext('D', adj = adj, padj = padj, cex = 3)

		plot(x5, y5, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
				 xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
		mtext('E', adj = adj, padj = padj, cex = 3)

		plot(x6, y6, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
				 xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
		mtext('F', adj = adj, padj = padj, cex = 3)

	layout(1)

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
#             Asymptotics Simulation, MLE, Table 8.3 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# function to create uniform distances
create_grid = function(type, n){
	if(type == 'infill') out = (1:n)/(n+1)
	if(type == 'outfill') out = (1:n)/51
	out
}

# profiled negative 2 times the loglikelihood to be minimized
# use analytical inverses for 1-D, it will be much faster than spmodel
m2LLprof = function(theta, z, del, n)
{
	# analytical inverse is available for 1-D exponential model
	# sqrt(D)%*%U matrix in Ver Hoef et al. (2009), Computing of Some Generalized 
	# Linear Mixed Pseudo-Models with Temporal Autocorrelation. Computational
	# Statistics 25(1): 39–55. DOI: 10.1007/s00180-009-0160-1, Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/theta)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/theta)/sqrt(1 - exp(-2*del/theta))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% z
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical determinant, Section 3.2 in Ver Hoef et al. (2009)
	SigLogDet = (n-1)*log(1 - exp(-2*del/theta))
	# this is equal to determinant(exp(-distmat/theta), logarithm = TRUE)$modulus
	# minus 2 times profiled likelihood, equation 8.4
	n*log(zQz) + SigLogDet + n*log(2*pi)
}

# theta is constant for all simulations
theta = 0.2

# set the random number seed so results are reproducible
set.seed(5001)

#-------------------------------------------------------------------------------
#                    Fixed-Domain, Sample Size 50
#-------------------------------------------------------------------------------

# pick sample size
n = 50
# create the grid
grd = create_grid('infill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_infill_50 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	thetahat
	# analytical inverse is available for 1-D exponential model
	# sqrt of inverse is sqrt(D)%*%U matrix 
	# in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/
		sqrt(1 - exp(-2*del/thetahat))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% zsim
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical solution for partial sill
	sigmahat = zQz/n
	
	#check it with spmodel
#	DF = data.frame(z = zsim, x = grd, y = rep(1, times = length(grd)))
#	library(spmodel)
#	splm(z ~ 1, xcoord = x, ycoord = y, data = DF,
#		spcov_initial = spcov_initial(spcov_type = 'exponential', 
#		ie = 0, known = 'ie'),
#		control = list(reltol = 1e-12), estmethod = 'ml')
	
	hold_infill_50[ithsim, 1] = thetahat
	hold_infill_50[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_infill_50[,1] - .2)^2)
mean((hold_infill_50[,2] - 1)^2)
mean((hold_infill_50[,2]/hold_infill_50[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    Fixed-Domain, Sample Size 250
#-------------------------------------------------------------------------------

# pick sample size
n = 250
# create the grid
grd = create_grid('infill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_infill_250 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	# analytical inverse is available for 1-D exponential model
	# sqrt(D)%*%U matrix in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/
		sqrt(1 - exp(-2*del/thetahat))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% zsim
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical solution for partial sill
	sigmahat = zQz/n
	hold_infill_250[ithsim, 1] = thetahat
	hold_infill_250[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_infill_250[,1] - .2)^2)
mean((hold_infill_250[,2] - 1)^2)
mean((hold_infill_250[,2]/hold_infill_250[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    Fixed-Domain, Sample Size 1000
#-------------------------------------------------------------------------------

# pick sample size
n = 1000
# create the grid
grd = create_grid('infill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_infill_1000 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	# analytical inverse is available for 1-D exponential model
	# sqrt(D)%*%U matrix in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/
		sqrt(1 - exp(-2*del/thetahat))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% zsim
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical solution for partial sill
	sigmahat = zQz/n
	hold_infill_1000[ithsim, 1] = thetahat
	hold_infill_1000[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_infill_50[,2]/hold_infill_50[,1] - 5)^2)
mean((hold_infill_250[,2]/hold_infill_250[,1] - 5)^2)
mean((hold_infill_1000[,2]/hold_infill_1000[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    Increasing-Domain, Sample Size 50
#-------------------------------------------------------------------------------

# pick sample size
n = 50
# create the grid
grd = create_grid('outfill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_outfill_50 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# analytical inverse is available for 1-D exponential model
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	# sqrt(D)%*%U matrix in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/sqrt(1 - exp(-2*del/thetahat))
	DDz = DD %*% zsim
	# analytical estimate of sigma^2
	sigmahat = sum(DDz^2)/n
	hold_outfill_50[ithsim, 1] = thetahat
	hold_outfill_50[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_outfill_50[,1] - .2)^2)
mean((hold_outfill_50[,2] - 1)^2)
mean((hold_outfill_50[,2]/hold_outfill_50[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    Increasing-Domain, Sample Size 250
#-------------------------------------------------------------------------------

# pick sample size
n = 250
# create the grid
grd = create_grid('outfill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_outfill_250 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	# sqrt(D)%*%U matrix in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/sqrt(1 - exp(-2*del/thetahat))
	DDz = DD %*% zsim
	# analytical estimate of sigma^2
	sigmahat = sum(DDz^2)/n
	hold_outfill_250[ithsim, 1] = thetahat
	hold_outfill_250[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_outfill_250[,1] - .2)^2)
mean((hold_outfill_250[,2] - 1)^2)
mean((hold_outfill_250[,2]/hold_outfill_250[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    Increasing-Domain, Sample Size 1000
#-------------------------------------------------------------------------------

# pick sample size
n = 1000
# create the grid
grd = create_grid('outfill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_outfill_1000 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	# sqrt(D)%*%U matrix in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/sqrt(1 - exp(-2*del/thetahat))
	DDz = DD %*% zsim
	# analytical estimate of sigma^2
	sigmahat = sum(DDz^2)/n
	hold_outfill_1000[ithsim, 1] = thetahat
	hold_outfill_1000[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_outfill_1000[,1] - .2)^2)
mean((hold_outfill_1000[,2] - 1)^2)
mean((hold_outfill_1000[,2]/hold_outfill_1000[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    MLE summary using MSE
#-------------------------------------------------------------------------------

#create output for latex, this is for MLE
hold_results = cbind(
	c(50,250,1000),
	c(
		mean((hold_infill_50[,2] - 1)^2),
		mean((hold_infill_250[,2] - 1)^2),
		mean((hold_infill_1000[,2] - 1)^2)
	),
	c(
		mean((hold_infill_50[,1] - .2)^2),
		mean((hold_infill_250[,1] - .2)^2),
		mean((hold_infill_1000[,1] - .2)^2)
	),
	c(
		mean((hold_infill_50[,2]/hold_infill_50[,1] - 5)^2),
		mean((hold_infill_250[,2]/hold_infill_250[,1] - 5)^2),
		mean((hold_infill_1000[,2]/hold_infill_1000[,1] - 5)^2)
	),
	c(
		mean((hold_outfill_50[,2] - 1)^2),
		mean((hold_outfill_250[,2] - 1)^2),
		mean((hold_outfill_1000[,2] - 1)^2)
	),
	c(
		mean((hold_outfill_50[,1] - .2)^2),
		mean((hold_outfill_250[,1] - .2)^2),
		mean((hold_outfill_1000[,1] - .2)^2)
	),
	c(
		mean((hold_outfill_50[,2]/hold_outfill_50[,1] - 5)^2),
		mean((hold_outfill_250[,2]/hold_outfill_250[,1] - 5)^2),
		mean((hold_outfill_1000[,2]/hold_outfill_1000[,1] - 5)^2)
	)
)

\boldsymbol{
#-------------------------------------------------------------------------------
#                    MLE summary using median rather than MSE
#-------------------------------------------------------------------------------

#create output for latex, this is for MLE
# use median rather than mean
hold_results = cbind(
	c(50,250,1000),
	c(
		median((hold_infill_50[,2] - 1)^2),
		median((hold_infill_250[,2] - 1)^2),
		median((hold_infill_1000[,2] - 1)^2)
	),
	c(
		median((hold_infill_50[,1] - .2)^2),
		median((hold_infill_250[,1] - .2)^2),
		median((hold_infill_1000[,1] - .2)^2)
	),
	c(
		median((hold_infill_50[,2]/hold_infill_50[,1] - 5)^2),
		median((hold_infill_250[,2]/hold_infill_250[,1] - 5)^2),
		median((hold_infill_1000[,2]/hold_infill_1000[,1] - 5)^2)
	),
	c(
		median((hold_outfill_50[,2] - 1)^2),
		median((hold_outfill_250[,2] - 1)^2),
		median((hold_outfill_1000[,2] - 1)^2)
	),
	c(
		median((hold_outfill_50[,1] - .2)^2),
		median((hold_outfill_250[,1] - .2)^2),
		median((hold_outfill_1000[,1] - .2)^2)
	),
	c(
		median((hold_outfill_50[,2]/hold_outfill_50[,1] - 5)^2),
		median((hold_outfill_250[,2]/hold_outfill_250[,1] - 5)^2),
		median((hold_outfill_1000[,2]/hold_outfill_1000[,1] - 5)^2)
	)
)

print(
    xtable(hold_results, 
      align = c('l',rep('l', times = length(hold_results[1,]))),
      digits = c(0,0,3,4,2,3,4,2),
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        Asymptotics Simulation, REMLE, not given as Table in book 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# function to create uniform distances
create_grid = function(type, n){
	if(type == 'infill') out = (1:n)/(n+1)
	if(type == 'outfill') out = (1:n)/51
	out
}

# profiled negative 2 times the loglikelihood to be minimized
# use analytical inverses for 1-D, it will be much faster than spmodel
m2LLprof = function(theta, z, del, n)
{
	# analytical inverse is available for 1-D exponential model
	# sqrt(D)%*%U matrix in Ver Hoef et al. (2009), Computing of Some Generalized 
	# Linear Mixed Pseudo-Models with Temporal Autocorrelation. Computational
	# Statistics 25(1): 39–55. DOI: 10.1007/s00180-009-0160-1, Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/theta)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/theta)/sqrt(1 - exp(-2*del/theta))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% z
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical determinant, Section 3.2 in Ver Hoef et al. (2009)
	SigLogDet = (n-1)*log(1 - exp(-2*del/theta))
	# this is equal to determinant(exp(-distmat/theta), logarithm = TRUE)$modulus
	# minus 2 times profiled likelihood, the equation after 8.5
	(n-1)*log(zQz) + SigLogDet + log(sum(DD1^2)) + (n-1)*log(2*pi) 
	# note that determinant(t(X) %*% solve(R) %*% X) is a scalar
}


# theta is constant for all simulations
theta = 0.2

# set the random number seed so results are reproducible
set.seed(5003)

#-------------------------------------------------------------------------------
#                    Fixed-Domain, Sample Size 50
#-------------------------------------------------------------------------------

# pick sample size
n = 50
# create the grid
grd = create_grid('infill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_infill_50 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	thetahat
	# analytical inverse is available for 1-D exponential model
	# sqrt of inverse is sqrt(D)%*%U matrix 
	# in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/
		sqrt(1 - exp(-2*del/thetahat))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% zsim
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical solution for partial sill
	sigmahat = zQz/(n-1)
	
	#check it with spmodel
#	DF = data.frame(z = zsim, x = grd, y = rep(1, times = length(grd)))
#	library(spmodel)
#	splm(z ~ 1, xcoord = x, ycoord = y, data = DF,
#		spcov_initial = spcov_initial(spcov_type = 'exponential', 
#		ie = 0, known = 'ie'),
#		control = list(reltol = 1e-12), estmethod = 'reml')
	
	hold_infill_50[ithsim, 1] = thetahat
	hold_infill_50[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_infill_50[,1] - .2)^2)
mean((hold_infill_50[,2] - 1)^2)
mean((hold_infill_50[,2]/hold_infill_50[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    Fixed-Domain, Sample Size 250
#-------------------------------------------------------------------------------

# pick sample size
n = 250
# create the grid
grd = create_grid('infill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_infill_250 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	thetahat
	# analytical inverse is available for 1-D exponential model
	# sqrt of inverse is sqrt(D)%*%U matrix 
	# in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/
		sqrt(1 - exp(-2*del/thetahat))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% zsim
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical solution for partial sill
	sigmahat = zQz/(n-1)
	hold_infill_250[ithsim, 1] = thetahat
	hold_infill_250[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_infill_250[,1] - .2)^2)
mean((hold_infill_250[,2] - 1)^2)
mean((hold_infill_250[,2]/hold_infill_250[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    Fixed-Domain, Sample Size 1000
#-------------------------------------------------------------------------------

# pick sample size
n = 1000
# create the grid
grd = create_grid('infill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_infill_1000 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	thetahat
	# analytical inverse is available for 1-D exponential model
	# sqrt of inverse is sqrt(D)%*%U matrix 
	# in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/
		sqrt(1 - exp(-2*del/thetahat))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% zsim
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical solution for partial sill
	sigmahat = zQz/(n-1)
	hold_infill_1000[ithsim, 1] = thetahat
	hold_infill_1000[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_infill_50[,2]/hold_infill_50[,1] - 5)^2)
mean((hold_infill_250[,2]/hold_infill_250[,1] - 5)^2)
mean((hold_infill_1000[,2]/hold_infill_1000[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    Increasing-Domain, Sample Size 50
#-------------------------------------------------------------------------------

# pick sample size
n = 50
# create the grid
grd = create_grid('outfill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_outfill_50 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	thetahat
	# analytical inverse is available for 1-D exponential model
	# sqrt of inverse is sqrt(D)%*%U matrix 
	# in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/
		sqrt(1 - exp(-2*del/thetahat))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% zsim
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical solution for partial sill
	sigmahat = zQz/(n-1)
	hold_outfill_50[ithsim, 1] = thetahat
	hold_outfill_50[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_outfill_50[,1] - .2)^2)
mean((hold_outfill_50[,2] - 1)^2)
mean((hold_outfill_50[,2]/hold_outfill_50[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    Increasing-Domain, Sample Size 250
#-------------------------------------------------------------------------------

# pick sample size
n = 250
# create the grid
grd = create_grid('outfill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_outfill_250 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	thetahat
	# analytical inverse is available for 1-D exponential model
	# sqrt of inverse is sqrt(D)%*%U matrix 
	# in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/
		sqrt(1 - exp(-2*del/thetahat))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% zsim
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical solution for partial sill
	sigmahat = zQz/(n-1)
	hold_outfill_250[ithsim, 1] = thetahat
	hold_outfill_250[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_outfill_250[,1] - .2)^2)
mean((hold_outfill_250[,2] - 1)^2)
mean((hold_outfill_250[,2]/hold_outfill_250[,1] - 5)^2)

#-------------------------------------------------------------------------------
#                    Increasing-Domain, Sample Size 1000
#-------------------------------------------------------------------------------

# pick sample size
n = 1000
# create the grid
grd = create_grid('outfill', n)
# get distance matrix
distmat = as.matrix(dist(grd))
# create covariance matrix
Sigma = exp(-distmat/theta)
# get Cholesky root of covariance matrix
SigChol = chol(Sigma)

# now start 1000 simulations
nsim = 1000
# create matrix to store results
hold_outfill_1000 = matrix(0, nrow = nsim, ncol = 2)
for(ithsim in 1:nsim) {
	cat("\r", "Simulation: ", ithsim)
	# simulate the data
	zsim = t(SigChol) %*% rnorm(n)
	# distance between successive grid points
	del = grd[2] - grd[1]
	# one dimensional optimization for theta
	thetahat = optimize(m2LLprof, lower = 0.00001, upper = 5, 
		z = zsim, del = del, n = n)$minimum
	thetahat
	# analytical inverse is available for 1-D exponential model
	# sqrt of inverse is sqrt(D)%*%U matrix 
	# in Ver Hoef et al. (2009), Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/thetahat)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/thetahat)/
		sqrt(1 - exp(-2*del/thetahat))
	# now t(DD) %*% DD = solve(exp(-distmat/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% zsim
	# sqrt of inverse of covariance matrix times vector of ones (row sums)
	DD1 = apply(DD,1,sum)
	# t(z) %*% Q(theta) %*% z where Q(theta) is scaled from equation (8.2)
	# because DDz and DD1 are vectors, we can take some computational shortcuts
	zQz = sum(DDz^2) - (sum(DD1*DDz))^2/sum(DD1^2)
	# analytical solution for partial sill
	sigmahat = zQz/(n-1)
	hold_outfill_1000[ithsim, 1] = thetahat
	hold_outfill_1000[ithsim, 2] = sigmahat
}
cat("\n")

mean((hold_outfill_1000[,1] - .2)^2)
mean((hold_outfill_1000[,2] - 1)^2)
mean((hold_outfill_1000[,2]/hold_outfill_1000[,1] - 5)^2)

#-------------------------------------------------------------------------------
#           REMLE summary using MSE
#-------------------------------------------------------------------------------

#create output for latex, this is for REMLE
hold_results = cbind(
	c(50,250,1000),
	c(
		mean((hold_infill_50[,2] - 1)^2),
		mean((hold_infill_250[,2] - 1)^2),
		mean((hold_infill_1000[,2] - 1)^2)
	),
	c(
		mean((hold_infill_50[,1] - .2)^2),
		mean((hold_infill_250[,1] - .2)^2),
		mean((hold_infill_1000[,1] - .2)^2)
	),
	c(
		mean((hold_infill_50[,2]/hold_infill_50[,1] - 5)^2),
		mean((hold_infill_250[,2]/hold_infill_250[,1] - 5)^2),
		mean((hold_infill_1000[,2]/hold_infill_1000[,1] - 5)^2)
	),
	c(
		mean((hold_outfill_50[,2] - 1)^2),
		mean((hold_outfill_250[,2] - 1)^2),
		mean((hold_outfill_1000[,2] - 1)^2)
	),
	c(
		mean((hold_outfill_50[,1] - .2)^2),
		mean((hold_outfill_250[,1] - .2)^2),
		mean((hold_outfill_1000[,1] - .2)^2)
	),
	c(
		mean((hold_outfill_50[,2]/hold_outfill_50[,1] - 5)^2),
		mean((hold_outfill_250[,2]/hold_outfill_250[,1] - 5)^2),
		mean((hold_outfill_1000[,2]/hold_outfill_1000[,1] - 5)^2)
	)
)

print(
    xtable(hold_results, 
      align = c('l',rep('l', times = length(hold_results[1,]))),
      digits = c(0,0,3,4,2,3,4,2),
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

#-------------------------------------------------------------------------------
#             REMLE summary using median rather than MSE
#-------------------------------------------------------------------------------

#create output for latex, this is for REMLE
# use median rather than mean
hold_results = cbind(
	c(50,250,1000),
	c(
		median((hold_infill_50[,2] - 1)^2),
		median((hold_infill_250[,2] - 1)^2),
		median((hold_infill_1000[,2] - 1)^2)
	),
	c(
		median((hold_infill_50[,1] - .2)^2),
		median((hold_infill_250[,1] - .2)^2),
		median((hold_infill_1000[,1] - .2)^2)
	),
	c(
		median((hold_infill_50[,2]/hold_infill_50[,1] - 5)^2),
		median((hold_infill_250[,2]/hold_infill_250[,1] - 5)^2),
		median((hold_infill_1000[,2]/hold_infill_1000[,1] - 5)^2)
	),
	c(
		median((hold_outfill_50[,2] - 1)^2),
		median((hold_outfill_250[,2] - 1)^2),
		median((hold_outfill_1000[,2] - 1)^2)
	),
	c(
		median((hold_outfill_50[,1] - .2)^2),
		median((hold_outfill_250[,1] - .2)^2),
		median((hold_outfill_1000[,1] - .2)^2)
	),
	c(
		median((hold_outfill_50[,2]/hold_outfill_50[,1] - 5)^2),
		median((hold_outfill_250[,2]/hold_outfill_250[,1] - 5)^2),
		median((hold_outfill_1000[,2]/hold_outfill_1000[,1] - 5)^2)
	)
)

print(
    xtable(hold_results, 
      align = c('l',rep('l', times = length(hold_results[1,]))),
      digits = c(0,0,3,4,2,3,4,2),
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

