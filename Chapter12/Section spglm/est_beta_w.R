est_beta_w = function(theta, y, X, Hdist, autocor_fun, stepsize)
{
# after optimizing for covariance parameters, we need betahat and w's
# we just need one pass through the likelihood function to obtain them
# so just copy that same code
	gam_1 = exp(theta[1])
	gam_2 = exp(theta[2])
	# set the nugget effect to a small value
	gam_0 = 0.00001
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
		wnew = w - stepsize*solve(H, g)
		wdiffmax = max(abs(wnew - w))
		w = wnew
	}

	# Now we need the variance covariance matrix of estimated betas
	# The covbeta above does not account for the fact that w's are latent
	# and unobserved. We need to take into account the variance of w's, 
	# which we have from Newton Raphson as (-H)^(-1)
	mHi = solve(-H)
	# Then the variance of (X'Sig^{-1}X)^{-1}X'Sig^{-1}w is given below
	covbetaHM = covbeta %*% t(SigiX) %*% mHi %*% SigiX %*% covbeta
	
	# return a list of estimated w and beta, 
	# and covariance matrix of w and beta
	list(w = w, betahat = betahat, cov_w = mHi, covbetaHM = covbetaHM,
		covbeta = covbeta, covbeta2 = covbetaHM + covbeta)
}
