logLik_Laplace = function(theta, y, X, Hdist, autocor_fun, shrink = 1)
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
	CovMat = gam_1*autocor_fun(Hdist,gam_2) + gam_0*diag(n)
	# inverse of covariance matrix
	CovMati = solve(CovMat)
	# covariance matrix of fixed effects
	covbeta = solve(t(X) %*% CovMati %*% X)
	# inverse covariance matrix times X
	SigiX = CovMati %*% X
	
	# define constants used outside of Newton-Raphson looping
	Constant1 = covbeta  %*% t(SigiX)
	Constant2 = -CovMati + SigiX %*% covbeta %*% t(SigiX)
	# loop to get the maximum value, where the gradient will be flat
	# (g is all zeros)
	wdiffmax = 1e+32
	niter = 0
#	browser()
#	for(i in 1:30) {
	while(wdiffmax > 1e-10 & niter < 50) {
		niter = niter + 1
		betahat = Constant1 %*% w
		# compute the d vector
		d = -exp(w) + y
		# and then the gradient vector
		g = d - CovMati %*% w + CovMati %*% X %*% betahat
		# Next, compute H
		H = diag(as.vector(-exp(w))) + Constant2
		# update w
		wnew = w - shrink*solve(H, g)
		wdiffmax = max(abs(wnew - w))
		w = wnew
	}

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
