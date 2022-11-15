simSGLM_wExpl = function(n, autocor_fun, betas, gammas, pred = TRUE, type = 'random')
{
	# set sample size
	# n = 200
	if(type == 'grid') {
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
	if(type == 'random') {
		# set mindist to a small number
		mindist = 1e-32
		# simulate coordinates while minimum distance is too small
		while(mindist < 0.001) {
			# generate coordinates randomly
			xycoord = data.frame(x = runif(n), y = runif(n))
			# get distances between all pairs of points
			Hdist = as.matrix(dist(xycoord))
			# check for minimum distance
			# if too small, covariance matrix may be singular
			mindist = min(Hdist[upper.tri(Hdist)])
		}
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

	# create true covariance matrix
	Hdist = as.matrix(dist(xycoord))
	Sig_true = gammas[1]*autocor_fun(Hdist,gammas[2]) + 
			gammas[3]*diag(nall)
	# check the correlation matrix
	cor_true = (1/sqrt(diag(Sig_true)))*t((1/sqrt(diag(Sig_true)))*Sig_true)
	cor_true[1:5,1:5]

	# cholesky of true covariance matrix
	Sig_chol = t(chol(Sig_true))
	# create true w's for simulating data
	w_true = X %*% betas + Sig_chol %*% rnorm(nall)
	# simulate y conditional on w's
	y = rpois(nall, lambda = exp(w_true))

	obspred = rep('obs', times = nall)
	if(pred == TRUE) obspred = c(rep('obs', times = n), 
		rep('pred', times = nall - n))
	# create a data.frame
	sim_data = data.frame(y = y, x_1 = data_dummy$x_1, x_2 = data_dummy$x_2,
	 xcoord = xycoord$x, ycoord = xycoord$y, w_true = w_true, 
	 obspred = obspred)
	sim_data
}
