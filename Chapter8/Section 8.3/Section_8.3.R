################################################################################
#-------------------------------------------------------------------------------
#              Asymptotics for Chapter 8
#-------------------------------------------------------------------------------
################################################################################

# function to create uniform distances on a line for Table 8.2
create_grid = function(type, n){
	if(type == 'infill') out = (1:n)/(n+1)
	if(type == 'outfill') out = (1:n)/51
	out
}

# profiled negative 2 times the loglikelihood to be minimized
m2LLprof = function(theta, z, del, n)
{
	# sqrt(D)%*%U matrix in Ver Hoef et al. (2009), Computing of Some Generalized 
	# Linear Mixed Pseudo-Models with Temporal Autocorrelation. Computational
	# Statistics 25(1): 39–55. DOI: 10.1007/s00180-009-0160-1, Section 3.2
	DD = diag(sqrt(rep(1/(1 - exp(-2*del/theta)), times = n)))
	DD[n,n] = 1
	DD[row(DD) - col(DD) == -1] = -exp(-del/theta)/sqrt(1 - exp(-2*del/theta))
	# sqrt of inverse of covariance matrix times data vector, Section 3.3 in
	# Ver Hoef et al. (2009)
	DDz = DD %*% z
	# analytical determinant, Section 3.2 in Ver Hoef et al. (2009)
	SigLogDet = (n-1)*log(1 - exp(-2*del/theta))
	# minus 2 times profiled likelihood, equation 8.4
	n*log(sum(DDz^2)) + SigLogDet
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

	m2LLprof(.2, zsim, del, n)

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

	m2LLprof(.2, zsim, del, n)

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

	m2LLprof(.2, zsim, del, n)

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

	m2LLprof(.2, zsim, del, n)

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

	m2LLprof(.2, zsim, del, n)

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

	m2LLprof(.2, zsim, del, n)

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

#create output for latex
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

library(xtable)
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
