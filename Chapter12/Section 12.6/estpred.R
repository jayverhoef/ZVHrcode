# prediction and estimation after estimating covariance parameters
# theta: covariance parameters in order 1) nugget, 2) partial sill, 3) range
#   these are on logit scale, and are converted to max*expit(theta[i])
#   internally where max is maxvar or maxrange
# y: vector of response variable
# X: design matrix of explanatory variables for observed locations
# distmat: matrix of Euclidean distance among all points
# autocor_fun: geostatistical autocorrelation function that takes distmat
#       and range parameter as it arguments
# maxvar: maximum variance for either partial sill or nugget effect
# maxrange: maximum for the range parameter
# family: distribution conditional on spatial random effects.  Currently,
#   only 'poisson' and 'binomial' are implemented
# Xp: design matrix of explanatory variables for prediction locations
# dist_op: Euclidean distance between observed and prediction locations
# dist_pp: Euclidean distance among prediction locations
estpred = function(theta, y, sampsize = NULL, X, distmat = NULL, autocor_fun, 
	maxvar, maxrange, maxphi = NULL, family, Z1 = NULL, Z2 = NULL, mask = NULL,
	Xp = NULL, Z1p = NULL, Z2p = NULL, dist_op = NULL, dist_pp = NULL,
	maskop = NULL, maskpp = NULL, W = NULL, M = NULL, use.nugget = TRUE)
{
# after optimizing for covariance parameters, we need betahat and w's
# we just need one pass through optimizing w's to obtain them
# so just copy that same code from likelihood function
	nobs = length(y)
	# transform theta
	gam_0 = NULL
	phi = NA
	revar1 = NA
	revar2 = NA
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
  theta_trans = c(gam_0, gam_1, gam_2, phi, revar1, revar2)
  names(theta_trans) = c('nugget', 'psill', 'range', 'phi', 'revar1', 'revar2')
  
	# number of observed locations
	n = length(y)
	# for binomial models
	if(is.null(sampsize)) sampsize = rep(1, times = n)
	# starting values for w
	if(family == 'poisson') w = 0.5*log(y + 1)
	if(family == 'binomial') w = 0.5*(y - 0.5*sampsize)
	if(family == 'negbinomial') w = 0.5*log(y + 1)
	if(family == 'gamma') w = 0.5*log(y + 1)
	if(family == 'invgauss') w = 0.5*log(y + 1)
	if(family == 'beta') w = 0.5*logit(y)
	
	# covariance matrix
	if(!is.null(W)){
		CovMat = gam_1*autocor_fun(W, M, gam_2) + gam_0*diag(n)
	} else {
		CovMat = gam_1*autocor_fun(distmat, gam_2) + gam_0*diag(n)
	}
	if(!is.null(Z1)) CovMat = CovMat + revar1*Z1 %*% t(Z1)
	if(!is.null(Z2)) CovMat = CovMat + revar1*Z2 %*% t(Z2)
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
		if(max(abs(gnew)) > max(abs(g))) wnew = w - 0.1*solveHg
		# compute change in w for convergence
#		wdiffmax = max(abs(wnew - w))
		# update w
		wdiffmax = max(abs(gnew))
		w = wnew
	}

	# estimated covariance matrix of the w's
	mHi = solve(-H)

	# weights for betahat
	wts_beta = covbeta %*% t(SigiX)
	# estimation of fixed effects
	betahat =  wts_beta %*% w
	# Then the variance of [wts_beta %*% w] is given below
	betawtsvarw = wts_beta %*% mHi %*% t(wts_beta)
	# adjusted covariance of betahats due to uncertainty in w
	covbeta_adj = covbeta + betawtsvarw
	
	# if we also want to make predictions
	w_pred = NULL
	wpred_se = NULL
	wpred_sew = NULL
	wpred_se_adj = NULL
	# do this part of Xp is not NULL
	if(!is.null(Xp)) {
		# number of prediction locations
		npred = dim(Xp)[1]
		# covariance matrix between observed and prediction locations
		R_op = gam_1*autocor_fun(dist_op,gam_2)
		# covariance matrix among prediction locations
		R_pp = gam_1*autocor_fun(dist_pp,gam_2) + gam_0*diag(npred)
		if(!is.na(revar1)) R_pp = R_pp + revar1*diag(dim(R_pp)[1])
		if(!is.null(Z1p)) {
			R_op = R_op + revar1*Z1 %*% t(Z1p)
			R_pp = R_pp + revar1*Z1p %*% t(Z1p)
		}
		if(!is.na(revar2)) R_pp = R_pp + revar2*diag(dim(R_pp)[1])
		if(!is.null(Z2p)) {
			R_op = R_op + revar2*Z2 %*% t(Z2p)
			R_pp = R_pp + revar2*Z2p %*% t(Z2p)
		}
		if(!is.null(maskop)) R_op = R_op * maskop
		if(!is.null(maskpp)) R_pp = R_pp * maskpp
		# prediction (kriging) weights
		wts_pred = Xp %*% wts_beta + t(R_op) %*% CovMati %*% 
			(diag(nobs) - X %*% wts_beta)
		# prediction
		w_pred = wts_pred %*% w
		# matrix of partial results to simplify algebra
		d1 = (Xp - t(R_op) %*% CovMati %*% X)
		# typical prediction (kriging) covariance matrix
		covpred = d1 %*% covbeta %*% t(d1) -
			t(R_op) %*% CovMati %*% R_op + R_pp
		# Then the variance of [wts_pred %*% w] is given below  
		predwtsvarw = wts_pred %*% mHi %*% t(wts_pred)
		# typical prediction standard errors
		wpred_se = sqrt(diag(covpred))
		# prediction standard errors due to variance of w
		wpred_sew = sqrt(diag(predwtsvarw))
		# adjust prediction standard errors due to uncertainty in w
		wpred_se_adj = sqrt(diag(covpred + predwtsvarw))
	}
	
	# return a list of estimated w and beta, 
	# and covariance matrix of w and beta, and predictions and 
	# prediction standard errors
	list(theta_trans = theta_trans, w = w, betahat = betahat, cov_w = mHi, 
		betawtsvarw = betawtsvarw,
		covbeta = covbeta, covbeta_adj = covbeta_adj,
		w_pred = w_pred, wpred_se = wpred_se, wpred_sew = wpred_sew,
		wpred_se_adj = wpred_se_adj)
}
