sec_path = 'Rcode/Chapter11/Section spglm/'
setwd(paste0(SLEDbook_path,sec_path))

source('simSGLM_wExpl.R')
source('autocorr_functions.R')
source('logLik_Laplace.R')
source('est_beta_w.R')

library(viridis)
library(classInt)

#set seed so reproducible
set.seed(1009)

# simulate data from spatial generalized linear model
# betas[1] is overall mean (on log scale)
# betas[2] is regression coefficient on simulated N(0,1) explanatory variables
# betas[3] is regression coefficient on simulated binary variables
# betas[4] is regression coefficient on normal*binary variable interaction
# gammas[1] is partial sill
# gammas[2] is range
# gammas[3] is nugget
betas = c(-.5, .5, -.5, .5)
betas = c(2, 0, 0, 0)
gammas = c(1, 1, 0.0001)
sim_data = simSGLM_wExpl(20^2, autocor_fun = rho_exp, betas = betas,
	gammas = gammas, type = 'grid', pred = FALSE)

cip = classIntervals(sim_data[,'y'], 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(sim_data$xcoord, sim_data$ycoord, col = cip_colors, pch = 15, cex = 4,
	cex.lab = 2, cex.axis = 1.5, xlab = 'x-coordinate',
	ylab = 'y-coordinate')

# create design matrix for observed data
X = model.matrix(~ x_1*x_2, data = sim_data[sim_data$obspred == 'obs',])
X = model.matrix(~ 1, data = sim_data[sim_data$obspred == 'obs',])
# get observed values
y = sim_data[sim_data$obspred == 'obs','y']
# get distances among observed data locations
Hdist = as.matrix(dist(sim_data[sim_data$obspred == 
	'obs',c('xcoord', 'ycoord')]))

#initial value for optim
theta = log(c(1,1))
# optimize for covariance parameters
# undebug(logLik_Laplace)
optout = optim(theta, logLik_Laplace, method = 'BFGS',
	y = y, X = X, Hdist = Hdist, autocor_fun = rho_exp, shrink = 0.1)
# covariance parameters
exp(optout$par)
# set theta as the optimized parameters on log scale
theta = optout$par

cip = classIntervals(sim_data$w_true, 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(sim_data$xcoord, sim_data$ycoord, col = cip_colors, pch = 15, cex = 4,
	cex.lab = 2, cex.axis = 1.5, xlab = 'x-coordinate',
	ylab = 'y-coordinate')

est_beta_w_out = est_beta_w(optout$par, y, X, Hdist, rho_exp, shrink = 0.1)

cip = classIntervals(est_beta_w_out$w, 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(sim_data$xcoord, sim_data$ycoord, col = cip_colors, pch = 15, cex = 4,
	cex.lab = 2, cex.axis = 1.5, xlab = 'x-coordinate',
	ylab = 'y-coordinate')



# make a graph for the likelihood surface for partial sill and range
# center on estimated values and scale by estimated value/10
theta1 = optout$par[1] + (-15:15)*optout$par[1]/5
theta2 = optout$par[2] + (-15:15)*optout$par[2]/5
# matrix to hold likelihood surface values
llgrid = matrix(NA, ncol = 3, nrow = length(theta1)*length(theta2))

# loop through various parameter values and compute -2*loglikelihood
iter = 0
for(i in 1:length(theta1)) {
	for(j in 1:length(theta2)) {
		iter = iter + 1
		llgrid[iter,1] = theta1[i]
		llgrid[iter,2] = theta2[j]
		llgrid[iter,3] = logLik_Laplace(theta = c(theta1[i],theta2[j]), 
			y = y, X = X, Hdist = Hdist, autocor_fun = rho_exp, shrink = .1)
	}
}

save(store,file = 'llgrid.rda')

# make a plot of the -2*loglikelihood surface
brks = quantile(llgrid[,3], probs = (0:9)/9)
cip = classIntervals(llgrid[,3], style = 'fixed', fixedBreaks= brks)
cip = classIntervals(llgrid[,3], style = 'fisher')
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(llgrid[,1:3], col = cip_colors, pch = 15, cex = 2,
	cex.lab = 2, cex.axis = 1.5, xlab = 'log(partial sill)',
	ylab = 'log(range)')
points(optout$par[1], optout$par[2], pch = 19, col = 'white', cex = 1.5)
points(0, 0, pch = 19, col = 'white', cex = 3)

exp(optout$par)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Simulation to test for bias and coverage when estimating fixed effects
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#set seed so reproducible
set.seed(1007)

niter = 2000
store = vector(mode='list', length = niter)
for(iter in 1:niter) {
	
	sim_data = simSGLM_wExpl(200, betas = c(-.5, .5, -.5, .5),
		gammas = c(1, 1, 0.0001))

	# create design matrix for observed data
	X = model.matrix(~ x_1*x_2, data = sim_data[sim_data$obspred == 'obs',])
	# get observed values
	y = sim_data[sim_data$obspred == 'obs','y']
	# get distances among observed data locations
	Hdist = as.matrix(dist(sim_data[sim_data$obspred == 
		'obs',c('xcoord', 'ycoord')]))

	#initial value for optim
	theta = log(c(1,1))
	# undebug(logLik_Laplace)
	# optimize for covariance parameters
	optout = optim(theta, logLik_Laplace, shrink = 0.1, # method = 'BFGS',
		y = y, X = X, Hdist = Hdist, autocor_fun = rho_exp)

	est_beta_w_out = est_beta_w(optout$par, y, X, Hdist, rho_exp, shrink = 0.1)
	est_beta_w_out$betahat

	# compare standard errors when using only covbeta versus the corrected one
	store[[iter]] = data.frame(True = betas, betahat = est_beta_w_out$betahat, 
		SE_corrected = sqrt(diag(est_beta_w_out$covbetaHM)), 
		CI90_corr = est_beta_w_out$betahat - 
		1.645*sqrt(diag(est_beta_w_out$covbetaHM)) < betas & 
			betas < est_beta_w_out$betahat + 
			1.645*sqrt(diag(est_beta_w_out$covbetaHM)),
		SE_uncorr = sqrt(diag(est_beta_w_out$covbeta)),
		CI90_uncorr = est_beta_w_out$betahat - 
			1.645*sqrt(diag(est_beta_w_out$covbeta)) < betas & 
			betas < est_beta_w_out$betahat + 
			1.645*sqrt(diag(est_beta_w_out$covbeta))
			)

}

save(store,file = 'store.rda')

i = 1
bias = rep(0, times = 4)
cover_cor = rep(0, times = 4)
cover_uncor = rep(0, times = 4)

for(i in 1:length(store)) {
	bias = bias + store[[i]]$betahat - store[[i]]$True
	cover_cor = cover_cor + store[[i]]$CI90_corr
	cover_uncor = cover_uncor + store[[i]]$CI90_uncorr
}
bias/niter
cover_cor/niter
cover_uncor/niter

# Run the code again with different (or random) seed.

# I haven't done it, but we need to do a simulation to make sure that
# the corrected standard errors have proper coverage.  I will make the 
# graduate student at Iowa State University do that.

# Note that I think that we will also need to modify the prediction equations
# to use (-H)^{-1} as the covariance matrix as well.  Again, I am going
# to have the ISU student work that out.

# make some graphs to compare simulated values to fits
# and true w's to estimated w's
layout(matrix(1:2, ncol = 1))
	par(mar = c(4,4,2,1))
	plot(y, exp(w), pch = 19)
	plot(w_true, w, pch = 19)
layout(1)


# create design matrix for observed data
X = model.matrix(~ x_1*x_2, data = sim_data[sim_data$obspred == 'obs',])
# get observed values
y = sim_data[sim_data$obspred == 'obs','y']
# get distances among observed data locations
Hdist = as.matrix(dist(sim_data[sim_data$obspred == 
	'obs',c('xcoord', 'ycoord')]))
