estpred = function(theta, y, X, dist, autocor_fun, maxvar, maxrange,
	Xp = NULL, dist_op = NULL, dist_pp = NULL)
{
# after optimizing for covariance parameters, we need betahat and w's
# we just need one pass through optimizing w's to obtain them
# so just copy that same code from likelihood function
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
	
	# compute spatial projection matrix
	mPtheta = CovMati - SigiX %*% covbeta %*% t(SigiX)
	
	wdiffmax = 1e+32
	niter = 0
#	browser()
#	for(i in 1:30) {
	while(wdiffmax > 1e-10 & niter < 50) {
		niter = niter + 1
		# the gradient vector
		g = y - exp(w) - mPtheta %*% w
		# Next, compute H
		H = -diag(as.vector(exp(w))) - mPtheta
		# update w
		wnew = w - stepsize*solve(H, g)
		wdiffmax = max(abs(wnew - w))
		w = wnew
	}

	mHi = solve(-H)

	wts_beta = covbeta %*% t(SigiX)
	# estimation of fixed effects
	betahat =  wts_beta %*% w
	# Then the variance of (X'Sig^{-1}X)^{-1}X'Sig^{-1}w is given below
	covbetaHM = wts_beta %*% mHi %*% t(wts_beta)
	
	w_pred = NULL
	w_se = NULL
	if(!is.null(Xp)) {
		npred = dim(Xp)[1]
		nobs = dim(X)[1]
		R_op = gam_1*autocor_fun(Dist_op,gam_2)
		R_pp = gam_1*autocor_fun(Dist_pp,gam_2) + + gam_0*diag(npred)

		wts_pred = (Xp %*% wts_beta + t(R_op) %*% CovMati %*% 
			(diag(nobs) - X %*% wts_beta))
		w_pred = wts_pred %*% w
		d1 = (Xp - t(R_op) %*% CovMati %*% X)
		w_se = sqrt(diag(d1 %*% covbeta %*% t(d1) -
			t(R_op) %*% CovMati %*% R_op + R_pp + 
			wts_pred %*% mHi %*% t(wts_pred)))
	}
	
	# return a list of estimated w and beta, 
	# and covariance matrix of w and beta, and predictions and 
	# prediction standard errors
	list(w = w, betahat = betahat, cov_w = mHi, covbetaHM = covbetaHM,
		covbeta = covbeta, covbeta2 = covbetaHM + covbeta,
		w_pred = w_pred, w_se = w_se)
}
