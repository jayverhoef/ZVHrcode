#dinvgauss = function(y, mean, phi, log = TRUE)
#{
#0.5*log(phi*mean/(2*pi*y^3)) - phi*(y - mean)^2/(2*mean*y)
#}
dinvgauss = statmod:::dinvgauss

# inverse logit function
expit = function(x) {exp(x)/(1 + exp(x))}
#  logit function
logit = function(x) {log(x/(1 - x))}

k_0 = function(w, phi, y) {
	digamma(phi*expit(w)) - digamma(phi/(1 + exp(w))) + log(1/y - 1)
}
k_1 = function(w, phi, y) {
	phi*(trigamma(phi*expit(w)) + trigamma(phi/(1 + exp(w)))) -
	2*sinh(w)*(k_0(w, phi, y) + 2*atanh(1 - 2*y))
}

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
# family: family for the data model, one of either "poisson", "binomial", 
#       "negbinomial", "gamma", or "invgaussian"				
logLik_Laplace = function(theta, y, sampsize = NULL, X, distmat = NULL, 
	autocor_fun, maxvar, maxrange, maxphi = NULL, Z1 = NULL, Z2 = NULL, family, 
	mlmeth = 'reml', mask = NULL, W = NULL, M = NULL, use.nugget = TRUE)
{
	if(any(abs(theta) > 20)) return(1e+32)
	# transform theta
	gam_0 = + 1e-6
	gam_1 = maxvar*expit(theta[1])
	gam_2 = maxrange*expit(theta[2])
	nparsofar = 2
	if(use.nugget == TRUE) {
		nparsofar = nparsofar + 1
		gam_0 = maxvar*expit(theta[nparsofar]) + 1e-6
	}
	if(family %in% c('negbinomial','gamma', 'invgauss', 'beta')) {
		phi = maxphi*expit(theta[nparsofar + 1])
		nparsofar = nparsofar + 1
	}
	if(!is.null(Z1)) {
		 revar1 = maxvar*expit(theta[nparsofar + 1])
		 nparsofar = nparsofar + 1	
	}	 
	if(!is.null(Z2)) {
		 revar2 = maxvar*expit(theta[nparsofar + 1])
		 nparsofar = nparsofar + 1	
	}	 

	# number of observed locations
	n = length(y)
	# for binomial models
	if(is.null(sampsize)) sampsize = rep(1, times = n)
	# starting values for w
	if(family == 'poisson') w = 0.5*log(y + 1)
	if(family == 'binomial') w = 0.5*(y - 0.5*sampsize)
	if(family == 'negbinomial') w = 0.5*log(y + 1)
	if(family == 'gamma') w = 0.5*log(y + 1)
	if(family == 'invgauss') w = 0.5*log(y + .1)
	if(family == 'beta') w = 0.5*logit(y)
	
	# covariance matrix
	if(!is.null(W)){
		CovMat = gam_1*autocor_fun(W, M, gam_2) + gam_0*diag(n)
	} else {
		CovMat = gam_1*autocor_fun(distmat, gam_2) + gam_0*diag(n)
	}
	if(!is.null(Z1)) CovMat = CovMat + revar1*Z1 %*% t(Z1)
	if(!is.null(Z2)) CovMat = CovMat + revar2*Z2 %*% t(Z2)
	if(!is.null(mask)) CovMat = CovMat * mask
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
	while(wdiffmax > 1e-6 & niter < 100) {
		niter = niter + 1
		niter
		# the gradient vector
		# compute d based on family
		if(family == 'poisson') d = y - exp(w)
		if(family == 'binomial') d = y - sampsize*expit(w)
		if(family == 'negbinomial') d = phi*(y - sampsize*exp(w))/(phi + exp(w))
		if(family == 'gamma') d = -phi + y*phi*exp(-w)
		if(family == 'invgauss') d = phi*(y/(2*exp(w)) - exp(w)/(2*y)) + 0.5
		if(family == 'beta') d = -phi*exp(w)*k_0(w, phi, y)/(exp(w) + 1)^2
		g =  d - mPtheta %*% w
		# Next, compute H
		# compute D based on family
		if(family == 'poisson') D = -diag(as.vector(exp(w)))
		if(family == 'binomial') 
			D = -diag(as.vector(sampsize*expit(w)/(1 + exp(w))))
		if(family == 'negbinomial') 
			D = -diag(as.vector(phi*exp(w)*(phi + y)/(phi + exp(w))^2))
		if(family == 'gamma') 
			D = -diag(as.vector(y*phi*exp(-w)))
		if(family == 'invgauss') 
			D = diag(as.vector(-phi*(exp(2*w) + y^2)/(2*y*exp(w))))
		if(family == 'beta') 
			D = diag(as.vector(-phi*exp(2*w)*k_1(w, phi, y)/(exp(w) + 1)^4))
		H = D - mPtheta 
		# compute new w
		solveHg = solve(H, g)
		wnew = w - solveHg
		# if g is not shrinking towards zero, decrease stepsize
		if(family == 'poisson') dnew = y - exp(wnew)
		if(family == 'binomial') dnew = y - sampsize*expit(wnew)
		if(family == 'negbinomial') dnew = 
			phi*(y - sampsize*exp(wnew))/(phi + exp(wnew))
		if(family == 'gamma') dnew = -phi + y*phi*exp(-wnew)
		if(family == 'invgauss') dnew = phi*(y/(2*exp(wnew)) - 
			exp(wnew)/(2*y)) + 0.5
		if(family == 'beta') dnew = -phi*exp(wnew)*k_0(wnew, phi, y)/
			(exp(wnew) + 1)^2
		gnew = dnew - mPtheta %*% wnew
		if(max(abs(gnew)) > max(abs(g)))  
			wnew = w - 0.1*solveHg
		# compute change in w for convergence
#		wdiffmax = max(abs(wnew - w))
		# update w
		wdiffmax = max(abs(gnew))
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
	if(family == 'negbinomial')
		dm = -2*sum(dnbinom(y, mu = exp(w), size = phi, log = TRUE)) 
	if(family == 'gamma')
		dm = -2*sum(dgamma(y, shape = phi, scale = exp(w)/phi, log = TRUE)) 
	if(family == 'invgauss')
		dm = -2*sum(dinvgauss(y, mean = exp(w), shape = phi*exp(w), log = TRUE)) 
	if(family == 'beta')
		dm = -2*sum(dbeta(y, shape1 = phi*expit(w), shape2 = phi*(1- expit(w)), 
		log = TRUE)) 
	mlmethpart = 0
	if(mlmeth == 'reml')
		mlmethpart = determinant(t(X) %*% SigiX, logarithm = TRUE)$modulus
	Likelihood =
		# data model part
		dm + mlmethpart +  		
		# multivariate normal part
		determinant(CovMat, logarithm = TRUE)$modulus + 
		t(w - X %*% betahat) %*% CovMati %*% (w - X %*% betahat) +
		# Laplace part
	  determinant(-H, logarithm = TRUE)$modulus 
	return(as.numeric(Likelihood))
}

