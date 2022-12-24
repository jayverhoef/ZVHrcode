# inverse logit function
expit = function(x) {exp(x)/(1 + exp(x))}

# loglikelihood using Laplace approximation
# theta: covariance parameters in order 1) nugget, 2) partial sill, 3) range
#   these are on logit scale, and are converted to max*expit(theta[i])
#   internally where max is maxvar or maxrange
# y: vector of response variable
# X: design matrix of explanatory variables
# distmat: matrix of Euclidean distance among all points
# autocor_fun: geostatistical autocorrelation function that take distmat
#       and range parameter as it arguments
# maxvar: maximum variance for either partial sill or nugget effect
# maxrange: maximum for the range parameter
logLik_Laplace = function(theta, y, sampsize = NULL, X, distmat, autocor_fun,
	maxvar, maxrange, family)
{
	# transform theta
	gam_0 = maxvar*expit(theta[1]) + 1e-6
	gam_1 = maxvar*expit(theta[2])
	gam_2 = maxrange*expit(theta[3])

	# number of observed locations
	n = length(y)
	# for binomial models
	if(is.null(sampsize)) sampsize = rep(1, times = n)
	# starting values for w
	if(family == 'poisson') w = 0.5*log(y + 1)
	if(family == 'binomial') w = 0.5*(y - 0.5*sampsize)
	
	# covariance matrix
	CovMat = gam_1*autocor_fun(distmat, gam_2) + gam_0*diag(n)
	# inverse of covariance matrix
	CovMati = solve(CovMat)
	# inverse covariance matrix times X
	SigiX = CovMati %*% X
	# covariance matrix of fixed effects
	covbeta = solve(t(X) %*% SigiX)
	
	# compute spatial projection matrix
	mPtheta = CovMati - SigiX %*% covbeta %*% t(SigiX)
	
	wdiffmax = 1e+32
	niter = 0
#	browser()
#	for(i in 1:30) {
	while(wdiffmax > 1e-8 & niter < 100) {
		niter = niter + 1
		# the gradient vector
		# compute d based on family
		if(family == 'poisson') d = y - exp(w)
		if(family == 'binomial') d = y - sampsize*expit(w)
		g =  d - mPtheta %*% w
		# Next, compute H
		# compute D based on family
		if(family == 'poisson') D = -diag(as.vector(exp(w)))
		if(family == 'binomial') 
			D = -diag(as.vector(sampsize*expit(w)/(1 + exp(w))))
		H = D - mPtheta
		# compute new w
		solveHg = solve(H, g)
		wnew = w - solveHg
		# if g is not shrinking towards zero, decrease stepsize
		if(family == 'poisson') dnew = y - exp(wnew)
		if(family == 'binomial') dnew = y - sampsize*expit(wnew)
		gnew = dnew - mPtheta %*% wnew
		if(max(abs(gnew)) > max(abs(g))) 
			wnew = w - 0.1*solveHg
		# compute change in w for convergence
		wdiffmax = max(abs(wnew - w))
		# update w
		w = wnew
	}

	wts_beta = covbeta %*% t(SigiX)
	# estimation of fixed effects
	betahat =  wts_beta %*% w

	#now compute the -2*likelihood using Laplace approximation
		# data model part
	if(family == 'poisson')
		dm = -2*sum(dpois(y, lambda = exp(w), log = TRUE)) 
	if(family == 'binomial')
		dm = -2*sum(dbinom(y, size = sampsize, prob = expit(w), log = TRUE)) 
	Likelihood =
		# data model part
		dm + 		
		# multivariate normal part
		determinant(CovMat, logarithm = TRUE)$modulus + 
			t(w - X %*% betahat) %*% CovMati %*% (w - X %*% betahat) +
		# REML part
		determinant(t(X) %*% SigiX, logarithm = TRUE)$modulus +
		# Laplace part
	  determinant(-H, logarithm = TRUE)$modulus 
	return(as.numeric(Likelihood))
}
