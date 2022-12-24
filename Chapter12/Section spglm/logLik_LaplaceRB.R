logLik_LaplaceRB = function(theta, y, X, dist_ok, kernel_fun, stepsize = 1,
	maxiter = 50)
{
	if(max(theta) > 10) return(1e+30)
	gam_1 = exp(theta[1])
	gam_2 = exp(theta[2])
	# set the nugget effect to a small value and hold constant
	gam_0 = 0.0001
	n = length(y)
	# starting values for w
	w = rep(0, times = n) 
	
	# covariance matrix
	Z_ok = kernel_fun(dist_ok, gam_2)
	CovMat = gam_1*Z_ok %*% t(Z_ok) + gam_0*diag(n)
	# inverse of covariance matrix
	CovMati = solve(CovMat)
	# inverse covariance matrix times X
	SigiX = CovMati %*% X
	# covariance matrix of fixed effects
	covbeta = solve(t(X) %*% SigiX)
		
	CovMati = solve(CovMat)
	
	# compute spatial projection matrix
	mPtheta = CovMati - SigiX %*% covbeta %*% t(SigiX)
	
	wdiffmax = 1e+32
	niter = 0
#	browser()
#	for(i in 1:30) {
	while(wdiffmax > 1e-8 & niter < maxiter) {
		niter = niter + 1
		# the gradient vector
		g =  y - exp(w) - mPtheta %*% w
		# Next, compute H
		H = -diag(as.vector(exp(w))) - mPtheta
		# update w
		wnew = w - stepsize*solve(H, g)
		wdiffmax = max(abs(wnew - w))
		w = wnew
	}

	wts_beta = covbeta %*% t(SigiX)
	# estimation of fixed effects
	betahat =  wts_beta %*% w

	#now compute the -2*likelihood using Laplace approximation
	Likelihood =
		# Poisson part
		-2*sum(log(dpois(y, lambda = exp(w)))) +
		# multivariate normal part
		determinant(CovMat, logarithm = TRUE)$modulus + 
		t(w - X %*% betahat) %*% CovMati %*% (w - X %*% betahat) +
		# REML part
		determinant(t(X) %*% SigiX, logarithm = TRUE)$modulus +
		# Laplace part
	  determinant(-H, logarithm = TRUE)$modulus
	return(as.numeric(Likelihood))
}

