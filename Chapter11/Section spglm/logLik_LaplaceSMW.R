logLik_LaplaceSMW = function(theta, y, X, dist_ok, kernel_fun, stepsize = 1,
	w_start = NULL)
{
	if(max(theta) > 10) return(1e+30)
	gam_1 = exp(theta[1])
	gam_2 = exp(theta[2])
	# set the nugget effect to a small value and hold constant
	gam_0 = 0.0001
	n = length(y)
	nknots = dim(dist_ok)[2]
	# starting values for w
	if(is.null(w_start)) {
		w = rep(0, times = n) } else {
		w = w_start
	}
	
	# covariance matrix
	Z_ok = kernel_fun(dist_ok, gam_2)
	CovMat = gam_1*Z_ok %*% t(Z_ok) + gam_0*diag(n)
	# inverse covariance matrix times X
	SigiX = X/gam_0 - (Z_ok/gam_0) %*% solve(diag(nknots)/gam_1 + 
		t(Z_ok) %*% Z_ok/gam_0, t(Z_ok) %*% X/gam_0)
	# covariance matrix of fixed effects
	covbeta = solve(t(X) %*% SigiX)
			
	wdiffmax = 1e+32
	niter = 0
#	browser()
#	for(i in 1:30) {
	while(wdiffmax > 1e-8 & niter < 50) {
		niter = niter + 1
		# the gradient vector
		w_star = w/gam_0 - (Z_ok/gam_0) %*% solve(diag(nknots)/gam_1 + 
			t(Z_ok) %*% Z_ok/gam_0, t(Z_ok) %*% w/gam_0)
		g = y - exp(w) - w_star + SigiX %*% covbeta %*% t(SigiX) %*% w
		# g = y - exp(w) - mPtheta %*% w
		# Next, compute H
		k = 1/(-exp(w) - 1/gam_0)
		Z_star = as.vector(k)*Z_ok/gam_0
		g_starstar = k*g - Z_star %*% solve(diag(nknots)/gam_1 + 
			t(Z_ok) %*% Z_ok/gam_0 + t(Z_star) %*% Z_ok/gam_0, t(Z_star) %*% g)
		X_starstar = k*SigiX - Z_star %*% solve(diag(nknots)/gam_1 + 
			t(Z_ok) %*% Z_ok/gam_0 + t(Z_star) %*% Z_ok/gam_0, t(Z_star) %*% SigiX)
		Hinvg = g_starstar - X_starstar %*% solve(t(X) %*% SigiX + 
			t(SigiX) %*% X_starstar, t(X_starstar) %*% g)
		# update w
		wnew = w - stepsize*Hinvg
		wdiffmax = max(abs(wnew - w))
		w = wnew
	}

	wts_beta = covbeta %*% t(SigiX)
	# estimation of fixed effects
	betahat =  wts_beta %*% w
	r = w - X %*% betahat
	r_star = r/gam_0 - (Z_ok/gam_0) %*% solve(diag(nknots)/gam_1 + 
			t(Z_ok) %*% Z_ok/gam_0, t(Z_ok) %*% r/gam_0)


	#now compute the -2*likelihood using Laplace approximation
	Likelihood =
		# Poisson part
		-2*sum(log(dpois(y, lambda = exp(w)))) +
		# multivariate normal part
		determinant(CovMat, logarithm = TRUE)$modulus + 
		sum(r*r_star) +
		# REML part
		- determinant(covbeta, logarithm = TRUE)$modulus +
		# Laplace part
		determinant(t(X) %*% SigiX - t(SigiX) %*% X_starstar, 
			logarithm = TRUE)$modulus -
		determinant(t(X) %*% SigiX, logarithm = TRUE)$modulus +
		determinant(diag(nknots)/gam_1 + t(Z_ok) %*% Z_ok/gam_0 + 
			t(Z_star) %*% Z_ok/gam_0, logarithm = TRUE)$modulus -
		determinant(diag(nknots)/gam_1 + t(Z_ok) %*% Z_ok/gam_0, 
			logarithm = TRUE)$modulus	+
		sum(log(abs(1/k)))
	return(as.numeric(Likelihood))
}

