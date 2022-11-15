est_beta_w = function(theta, y, X, Hdist, autocor_fun, shrink)
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
	
	# loop to get the maximum value, where the gradient will be flat
	# (g is all zeros)
	for(i in 1:30) {
		betahat = Constant1 %*% w
		# compute the d vector
		d = -exp(w) + y
		# and then the gradient vector
		g = d - CovMati %*% w + CovMati %*% X %*% betahat
		# Next, compute H
		H = diag(as.vector(-exp(w))) + Constant2
		# update w
		w = w - shrink*solve(H, g)
	}

	# Now we need the variance covariance matrix of estimated betas
	# The covbeta above does not account for the fact that w's are latent
	# and unobserved. We need to take into account the variance of w's, 
	# which we have from Newton Raphson as (-H)^(-1)
	mHi = solve(-H)
	# kriging predictor using w and H
	xyobs = sim_data[sim_data$obspred == 'obs',c('xcoord','ycoord')]
	xypred = sim_data[sim_data$obspred == 'pred',c('xcoord','ycoord')]
	nobs = dim(xyobs)[1]
	npred = dim(xypred)[1]
	Dist_op = sqrt((outer(xyobs[,1], rep(1, times = npred)) - 
		outer(rep(1, times = nobs),xypred[,1]))^2 +
		(outer(xyobs[,2], rep(1, times = npred)) - 
		outer(rep(1, times = nobs),xypred[,2]))^2 )
	R_op = gam_1*autocor_fun(Dist_op,gam_2)
	Xp = model.matrix(~ 1, data = sim_data[sim_data$obspred == 'pred',])

	wtsMat = (Xp %*% Constant1 + t(R_op) %*% CovMati %*% 
		(diag(nobs) - X %*% Constant1))
	w_pred = wtsMat %*% w
	plot(sim_data[sim_data$obspred == 'pred','w_true'], w_pred, pch = 19)
	w_se = sqrt(diag(wtsMat %*% mHi %*% t(wtsMat)))
	w_pred - 1.645*w_se < sim_data[sim_data$obspred == 'pred','w_true'] &
		w_pred + 1.645*w_se > sim_data[sim_data$obspred == 'pred','w_true']
	# return a list of estimated w and beta, 
	# and covariance matrix of w and beta
	list(w = w, betahat = betahat, cov_w = mHi, covbetaHM = covbetaHM,
		covbeta = covbeta)
}
