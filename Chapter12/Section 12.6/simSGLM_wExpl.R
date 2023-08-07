expit = function(x){exp(x)/(1 + exp(x))}
	
simSGLM_wExpl = function(n, betas, gammas, loc_type = 'random', pred = TRUE,
	autocorr_type = 'covariance', autocor_fun = rho_exp,
	nknots = NULL, sampsize = NULL, family = 'poisson')
{
	if(family == 'invgauss') require('statmod')
	
	if(loc_type == 'grid') {
		ndim = round(sqrt(n),0)
		xycoord = data.frame(x = kronecker(((1:ndim) - 0.5)/ndim, 
			rep(1, times = ndim)), y = kronecker(rep(1, times = ndim), 
			((1:ndim) - 0.5)/ndim))
		n = ndim^2
		nall = n
		if(pred == TRUE) {
		# create a prediction grid
			pgrid = data.frame(x = kronecker(((1:10) - 0.5)/10, 
				rep(1, times = 10)), y = kronecker(rep(1, times = 10), 
				((1:10) - 0.5)/10))
			xycoord = rbind(xycoord, pgrid)
			nall = n + dim(pgrid)[1]
		}
	}
	if(loc_type == 'random') {
		# generate coordinates randomly
		xycoord = data.frame(x = runif(n), y = runif(n))
		# get distances between all pairs of points
		Hdist = as.matrix(dist(xycoord))
		# check for minimum distance
		nall = n
		if(pred == TRUE) {
		# create a prediction grid
			pgrid = data.frame(x = kronecker(((1:10) - 0.5)/10, 
				rep(1, times = 10)), y = kronecker(rep(1, times = 10), 
				((1:10) - 0.5)/10))
			xycoord = rbind(xycoord, pgrid)
			nall = n + dim(pgrid)[1]
		}
	}
	# create a data.frame with simulated covariates
	# x_1 is continuous from normal, x_2 is Bernoulli (binary factor variable)
	data_dummy = data.frame(x_1 = rnorm(nall),
		x_2 = as.factor(rbinom(nall, 1, 0.5)))

	# create a design matrix that includes x_1, x_2, and their interaction
	X = model.matrix(~ x_1*x_2, data= data_dummy) 

	# choose some true betas for simulation
	# betas = c(-.5, .5, -.5, .5)
	# choose some true covariance parameters for simulation
	# gammas = c(1, 1, 0.0001)

	if(is.null(sampsize)) sampsize = rep(1, times = n)
	if(autocorr_type == 'covariance') {
		# create true covariance matrix
		Hdist = as.matrix(dist(xycoord))
		Sig_true = gammas[1]*autocor_fun(Hdist,gammas[2]) + 
				gammas[3]*diag(nall)
		# cholesky of true covariance matrix
		Sig_chol = t(chol(Sig_true))
		# create true w's for simulating data
		w_true = X %*% betas + Sig_chol %*% rnorm(nall)
		# simulate y conditional on w's
		if(family == 'poisson') y = rpois(nall, lambda = exp(w_true))
		if(family == 'binomial') y = rbinom(nall, size = sampsize, 
			p = expit(w_true))
		if(family == 'negbinomial') y = rnbinom(nall, mu = exp(w_true), 
			size = exp(gammas[4]))
		if(family == 'beta') y = rbeta(nall, shape1 = gammas[4]*expit(w_true), 
			shape2 = gammas[4]*(1- expit(w_true)))
		if(family == 'invgauss') y = rinvgauss(nall, mean = exp(w_true), 
			shape = exp(w_true)*gammas[4])
		if(family == 'gamma') y = rgamma(nall, shape = gammas[4], 
			scale = exp(w_true)/gammas[4])
	}
	if(autocorr_type == 'radialbasis') {
		knots = kmeans(xycoord, nknots)$centers
		# distance between observed locations and knots
		Z_dist = as.matrix(pdist(xycoord, knots))
		# create w values as linear combination of kernel distances
		w_true = gammas[1]*kern_epi(Z_dist, gammas[2]) %*% rnorm(nknots) + 0.0001
		# simulate y conditional on w's
		if(family == 'poisson') y = rpois(nall, lambda = exp(w_true))
		if(family == 'binomial') y = rbinom(nall, size = sampsize, 
			p = expit(w_true))
		if(family == 'negbinomial') y = rnbinom(nall, mu = exp(w_true), 
			size = exp(gammas[4]))
		if(family == 'beta') y = rbeta(nall, shape1 = gammas[4]*expit(w_true), 
			shape2 = gammas[4]*(1- expit(w_true)))
		if(family == 'invgauss') y = rinvgauss(nall, mean = exp(w_true), 
			shape = exp(w_true)*exp(gammas[4]))
	}
	obspred = rep('obs', times = nall)
	if(pred == TRUE) obspred = c(rep('obs', times = n), 
		rep('pred', times = nall - n))
	# create a data.frame
	sim_data = data.frame(y = y, x_1 = data_dummy$x_1, x_2 = data_dummy$x_2,
	 xcoord = xycoord$x, ycoord = xycoord$y, w_true = w_true, 
	 obspred = obspred)
	sim_data
}
