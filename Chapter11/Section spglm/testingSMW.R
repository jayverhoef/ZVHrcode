sec_path = 'Rcode/Chapter11/Section spglm/'
setwd(paste0(SLEDbook_path,sec_path))

source('simSGLM_wExpl.R')
source('autocorr_functions.R')
source('logLik_Laplace.R')
source('logLik_LaplaceSMW.R')
source('est_beta_w.R')
source('addBreakColorLegend.R')
source('mginv.R')

library(viridis)
library(classInt)
library(xtable)
library(pdist)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Simulate SGLM Poisson data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#set seed so reproducible
set.seed(1012)

# simulate data from spatial generalized linear model
# betas[1] is overall mean (on log scale)
# betas[2] is regression coefficient on simulated N(0,1) explanatory variables
# betas[3] is regression coefficient on simulated binary variables
# betas[4] is regression coefficient on normal*binary variable interaction
# gammas[1] is partial sill
# gammas[2] is range
# gammas[3] is nugget
#betas = c(-.5, .5, -.5, .5)
betas = c(2, 0, 0, 0)
gammas = c(1, 1, 0.0001)
sim_data = simSGLM_wExpl(20^2, autocor_fun = rho_exp, betas = betas,
	gammas = gammas, type = 'grid', pred = FALSE)
xyobs = sim_data[sim_data$obspred == 'obs',c('xcoord','ycoord')]
xypred = sim_data[sim_data$obspred == 'pred',c('xcoord','ycoord')]
npred = dim(xypred)[1]
nobs = dim(xyobs)[1]

epikern = function(dist, gam_2) 1 - (dist/gam_2)^2*(dist < gam_2)

set.seed(2021)
nknots = 80
knots = kmeans(xyobs, nknots)$centers
dist_ok = as.matrix(pdist(xyobs, knots))

# create design matrix for observed data
X = model.matrix(~ 1, data = sim_data[sim_data$obspred == 'obs',])

# get observed values
y = sim_data[sim_data$obspred == 'obs','y']

layout(matrix(1:4, nrow = 1), widths = rep(c(3,1), times = 2))

brks_cex = 1.5
mtext_cex = 3.5
leg_right = 0.6
adj = -0.21
padj = -0.3
cex_plot = 4
cip = classIntervals(sim_data[,'y'], 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(sim_data$xcoord, sim_data$ycoord, col = cip_colors, pch = 15, 
	cex = cex_plot, cex.lab = 2, cex.axis = 1.5, xlab = 'x-coordinate',
	ylab = 'y-coordinate')
mtext('A', cex = mtext_cex, adj = adj, padj = padj)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = "1.0")

cip = classIntervals(sim_data$w_true, 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(sim_data$xcoord, sim_data$ycoord, col = cip_colors, pch = 15, 
	cex = cex_plot, cex.lab = 2, cex.axis = 1.5, xlab = 'x-coordinate',
	ylab = 'y-coordinate')
points(knots, pch = 19, cex = 2, col = 'white')
mtext('B', cex = mtext_cex, adj = adj, padj = padj)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = "1.2")
layout(1)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Estimate covariance parameters using Laplace
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#initial value for optim
theta = log(c(1,1))
# optimize for covariance parameters
# undebug(logLik_Laplace)
optout = optim(theta, logLik_LaplaceSMW, # method = 'BFGS',
	y = y, X = X, dist_ok = dist_ok, kernel_fun = epikern, stepsize = 1)
# covariance parameters
exp(optout$par)
# set theta as the optimized parameters on log scale
theta = optout$par

logLik_LaplaceSMW(theta + .1, y = y, X = X, dist_ok = dist_ok, 
	kernel_fun = epikern, stepsize = 1)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Estimate fixed effects and spatial random effects after estimating
#       covariance parameters
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

pred_w_out = pred_wSMW(optout$par, y, X, dist_ok, epikern, stepsize = 1)
pred_w_out$betahat

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Create values on a grid for the likelihood surface
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

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
			y = y, X = X, Hdist = Hdist, autocor_fun = rho_exp, stepsize = 1)
	}
}

#  save(llgrid, file = 'llgrid.rda')
#  load('llgrid.rda', verbose = TRUE)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Figure of simulated data, spatial random effect and their prediction,
#                and the likelihood surface
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/SGLM_likelihood_estimation"

pdf(paste0(file_name,'.pdf'), width = 11, height = 9.2)

brks_cex = 1.5
mtext_cex = 3.5
leg_right = 0.6
adj = -0.21
padj = -0.3
cex_plot = 4

layout(matrix(1:8, nrow = 2, byrow = TRUE), widths = rep(c(3,1), times = 4))

cip = classIntervals(sim_data[,'y'], 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(sim_data$xcoord, sim_data$ycoord, col = cip_colors, pch = 15, 
	cex = cex_plot, cex.lab = 2, cex.axis = 1.5, xlab = 'x-coordinate',
	ylab = 'y-coordinate')
mtext('A', cex = mtext_cex, adj = adj, padj = padj)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = "1.0")

cip = classIntervals(sim_data$w_true, 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(sim_data$xcoord, sim_data$ycoord, col = cip_colors, pch = 15, 
	cex = cex_plot, cex.lab = 2, cex.axis = 1.5, xlab = 'x-coordinate',
	ylab = 'y-coordinate')
mtext('B', cex = mtext_cex, adj = adj, padj = padj)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = "1.2")

# make a plot of the -2*loglikelihood surface
brks = quantile(llgrid[,3], probs = (0:9)/9)
cip = classIntervals(llgrid[,3], style = 'fixed', fixedBreaks= brks)
cip = classIntervals(llgrid[,3], style = 'fisher')
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(llgrid[,1:3], col = cip_colors, pch = 15, cex = cex_plot,
	cex.lab = 2, cex.axis = 1.5, xlab = 'log(partial sill)',
	ylab = 'log(range)')
points(optout$par[1], optout$par[2], pch = 19, col = 'white', cex = 1.5)
points(0, 0, pch = 19, col = 'white', cex = 3)
mtext('C', cex = mtext_cex, adj = adj, padj = padj)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = "1.1")

cip = classIntervals(pred_w_out$w, 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(sim_data$xcoord, sim_data$ycoord, col = cip_colors, pch = 15, 
	cex = cex_plot, cex.lab = 2, cex.axis = 1.5, xlab = 'x-coordinate',
	ylab = 'y-coordinate')
mtext('D', cex = mtext_cex, adj = adj, padj = padj)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = "1.2")


layout(1)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Simulation to test for bias and coverage when estimating fixed effects
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#set seed so reproducible
set.seed(1007)

betas = c(-.5, .5, -.5, .5)
gammas = c(1, 1, 0.0001)
niter = 2000
store = vector(mode='list', length = niter)
start_time = Sys.time()
for(iter in 1:niter) {
	cat("\r", "iteration: ", iter)
	sim_data = simSGLM_wExpl(200, autocor_fun = rho_exp, 
		betas = betas, gammas = gammas, 
		type = 'random', pred = TRUE)
	xyobs = sim_data[sim_data$obspred == 'obs',c('xcoord','ycoord')]
	xypred = sim_data[sim_data$obspred == 'pred',c('xcoord','ycoord')]
	npred = dim(xypred)[1]
	nobs = dim(xyobs)[1]

	# create design matrix for observed data
	X = model.matrix(~ x_1*x_2, data = sim_data[sim_data$obspred == 'obs',])
	Xp = model.matrix(~ x_1*x_2, data = sim_data[sim_data$obspred == 'pred',])

	# get observed values
	y = sim_data[sim_data$obspred == 'obs','y']
	# get distances among observed data locations
	Hdist = as.matrix(dist(xyobs))
	Dist_op = sqrt((outer(xyobs[,1], rep(1, times = npred)) - 
		outer(rep(1, times = nobs),xypred[,1]))^2 +
		(outer(xyobs[,2], rep(1, times = npred)) - 
		outer(rep(1, times = nobs),xypred[,2]))^2 )
	Dist_pp = as.matrix(dist(xypred))

	# starting values for spatial random effects
	m1 = glm(y ~ x_1*x_2, data = sim_data[sim_data$obspred == 'obs',], 
		family="poisson")
	# use signed values of log of absolute values of residuals
	w_start = ((resid(m1) < 0)*-1 + (resid(m1) > 0)*1)*log(abs(resid(m1)))
	#initial covariance parameters values for optim
	theta = log(c(1,1))
	# undebug(logLik_Laplace)
	# optimize for covariance parameters
	optout = optim(theta, logLik_Laplace, stepsize = 1, # method = 'BFGS',
		y = y, X = X, Hdist = Hdist, autocor_fun = rho_exp)

  pred_out = pred_w(optout$par, y, X, Hdist, rho_exp, stepsize = 1, 
		Xp = Xp, dist_op = dist_op, dist_pp = dist_pp)
	est_beta_w_out = pred_out
  w_true = sim_data$w_true[sim_data$obspred == 'pred']
  bias_pred = mean(pred_out$w_pred - w_true)
  cover_pred = sum(pred_out$w_pred - 1.645*pred_out$w_se < w_true & 
		w_true < pred_out$w_pred + 1.645*pred_out$w_se)/100
	cbind(w_true, pred_out$w_pred, pred_out$w_se)
	plot(w_true, pred_out$w_pred)
	# compare standard errors when using only covbeta versus the corrected one
#	store[[iter]] = 
	data.frame(True = betas, betahat = est_beta_w_out$betahat, 
		SE_corrected = sqrt(diag(est_beta_w_out$covbetaHM)), 
		CI90_corr = est_beta_w_out$betahat - 
		1.645*sqrt(diag(est_beta_w_out$covbetaHM)) < betas & 
			betas < est_beta_w_out$betahat + 
			1.645*sqrt(diag(est_beta_w_out$covbetaHM)),
		SE_uncorr = sqrt(diag(est_beta_w_out$covbeta)),
		CI90_uncorr = est_beta_w_out$betahat - 
			1.645*sqrt(diag(est_beta_w_out$covbeta)) < betas & 
			betas < est_beta_w_out$betahat + 
			1.645*sqrt(diag(est_beta_w_out$covbeta)),
		SE_corr2 = sqrt(diag(est_beta_w_out$covbeta2)),
		CI90_corr2 = est_beta_w_out$betahat - 
			1.645*sqrt(diag(est_beta_w_out$covbeta2)) < betas & 
			betas < est_beta_w_out$betahat + 
			1.645*sqrt(diag(est_beta_w_out$covbeta2)),
		bias_pred = bias_pred,
		cover_pred = cover_pred
	)
exp(optout$par)

}
cat("\n")
end_time = Sys.time()
difftime(end_time, start_time)
#  save(store,file = 'store.r)da')
#  load('store.rda')

i = 1
bias = rep(0, times = 4)
cover_cor = rep(0, times = 4)
cover_cor2 = rep(0, times = 4)
cover_uncor = rep(0, times = 4)
avevar_cor = rep(0, times = 4)
pred_bias = 0
cover_pred = 0

for(i in 1:length(store)) {
	bias = bias + store[[i]]$betahat - store[[i]]$True
	pred_bias = pred_bias + store[[i]]$bias_pred[1]
	avevar_cor = avevar_cor + store[[i]]$SE_corrected^2
	cover_cor = cover_cor + store[[i]]$CI90_corr
	cover_cor2 = cover_cor2 + store[[i]]$CI90_corr2
	cover_uncor = cover_uncor + store[[i]]$CI90_uncorr
	cover_pred = cover_pred + store[[i]]$cover_pred[1]

}
sglm_fe = data.frame(bias =  bias/niter,
	cover_corrected = cover_cor2/niter,
#	cover_corHonly = cover_cor/niter,
	cover_uncorr = cover_uncor/niter
)
sglm_fe = rbind(sglm_fe,
	c(pred_bias/niter, cover_pred/niter, NA))

print(
    xtable(sglm_fe, 
      align = c('l',rep('l', times = length(sglm_fe[1,]))),
      digits = c(0, rep(3, times = 3))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)


# bias by effect
bvec = NULL
for(i in 1:length(store)) 
	bvec = cbind(bvec, store[[i]]$betahat - store[[i]]$True)
hist(bvec[3,])

# Run the code again with different (or random) seed.

# I haven't done it, but we need to do a simulation to make sure that
# the corrected standard errors have proper coverage.  I will make the 
# graduate student at Iowa State University do that.

# Note that I think that we will also need to modify the prediction equations
# to use (-H)^{-1} as the covariance matrix as well.  Again, I am going
# to have the ISU student work that out.

# make some graphs to compare simulated values to fits
# and true w's to estimated w's
layout(matrix(1:4, ncol = 2, byrow = TRUE))
	par(mar = c(4,4,2,1))
	plot(exp(sim_data$w_true), y, pch = 19)
	plot(exp(est_beta_w_out$w), y, pch = 19)
	plot(sim_data$w_true, est_beta_w_out$w, pch = 19)
layout(1)


est_beta_w_out$w

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Testing Estimation
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
	stepsize = 1
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
		w = w - stepsize*solve(H, g)
	}

	# Now we need the variance covariance matrix of estimated betas
	# The covbeta above does not account for the fact that w's are latent
	# and unobserved. We need to take into account the variance of w's, 
	# which we have from Newton Raphson as (-H)^(-1)
	mHi = solve(-H)
	# Then the variance of (X'Sig^{-1}X)^{-1}X'Sig^{-1}w is given below
	covbetaHM = covbeta %*% t(SigiX) %*% mHi %*% SigiX %*% covbeta

data.frame(True = betas, betahat = betahat, 
	SE_corrected = sqrt(diag(covbetaHM)), 
	CI90_corr = betahat - 
	1.645*sqrt(diag(covbetaHM)) < betas & 
		betas < betahat + 
		1.645*sqrt(diag(covbetaHM)),
	SE_uncorr = sqrt(diag(covbeta)),
	CI90_uncorr = betahat - 
		1.645*sqrt(diag(covbeta)) < betas & 
		betas < betahat + 
		1.645*sqrt(diag(covbeta)),
	SE_corr2 = sqrt(diag(covbeta + covbetaHM)),
	CI90_uncorr = betahat - 
		1.645*sqrt(diag(covbeta + covbetaHM)) < betas & 
		betas < betahat + 
		1.645*sqrt(diag(covbeta + covbetaHM))
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Testing Prediction
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


#set seed so reproducible
set.seed(1091)

shrink = .1
# simulate data from spatial generalized linear model
# betas[1] is overall mean (on log scale)
# betas[2] is regression coefficient on simulated N(0,1) explanatory variables
# betas[3] is regression coefficient on simulated binary variables
# betas[4] is regression coefficient on normal*binary variable interaction
# gammas[1] is partial sill
# gammas[2] is range
# gammas[3] is nugget
betas = c(-.5, .5, -.5, .5)
betas = c(1, 0, 0, 0)
gammas = c(1, 1, 0.0001)
sim_data = simSGLM_wExpl(200, autocor_fun = rho_exp, betas = betas,
	gammas = gammas, type = 'random', pred = TRUE)

sim_data$y

# create design matrix for observed data
X = model.matrix(~ x_1*x_2, data = sim_data[sim_data$obspred == 'obs',])
X = model.matrix(~ 1, data = sim_data[sim_data$obspred == 'obs',])
# get observed values
y = sim_data[sim_data$obspred == 'obs','y']
# get distances among observed data locations
Hdist = as.matrix(dist(sim_data[sim_data$obspred == 
	'obs',c('xcoord', 'ycoord')]))

#-------------------------------------------------------------------------------
#   Estimate covariance parameters using Laplace
#-------------------------------------------------------------------------------

#initial value for optim
theta = log(c(1,1))
# optimize for covariance parameters
# undebug(logLik_Laplace)
optout = optim(theta, logLik_Laplace, method = 'BFGS',
	y = y, X = X, Hdist = Hdist, autocor_fun = rho_exp, shrink = 1)
# covariance parameters
exp(optout$par)
# set theta as the optimized parameters on log scale
theta = optout$par

#-------------------------------------------------------------------------------
#   prediction
#-------------------------------------------------------------------------------

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
	mH = -H
	# kriging predictor using w and H
	xyobs = sim_data[sim_data$obspred == 'obs',c('xcoord','ycoord')]
	xypred = sim_data[sim_data$obspred == 'pred',c('xcoord','ycoord')]
	nobs = dim(xyobs)[1]
	npred = dim(xypred)[1]
	Dist_op = sqrt((outer(xyobs[,1], rep(1, times = npred)) - 
		outer(rep(1, times = nobs),xypred[,1]))^2 +
		(outer(xyobs[,2], rep(1, times = npred)) - 
		outer(rep(1, times = nobs),xypred[,2]))^2 )
	Dist_pp = as.matrix(dist(xypred))
	R_op = gam_1*autocor_fun(Dist_op,gam_2)
	R_pp = gam_1*autocor_fun(Dist_pp,gam_2) + + gam_0*diag(npred)
	Xp = model.matrix(~ 1, data = sim_data[sim_data$obspred == 'pred',])

	wtsMat = (Xp %*% Constant1 + t(R_op) %*% CovMati %*% 
		(diag(nobs) - X %*% Constant1))
	w_pred = wtsMat %*% w
	w_true = sim_data[sim_data$obspred == 'pred','w_true']
	
	plot(w_true, w_pred, pch = 19)
	d1 = (Xp - t(R_op) %*% CovMati %*% X)
	w_se = sqrt(diag(d1 %*% solve(t(X) %*% CovMati %*% X) %*% t(d1) -
		t(R_op) %*% CovMati %*% R_op + R_pp + 
		wtsMat %*% mHi %*% t(wtsMat)))
	pred_cover = w_pred - 1.645*w_se < w_true &
		w_pred + 1.645*w_se > w_true
	mean(pred_cover)
	exp(theta)

junk1 = wtsMat %*% CovMat %*% t(wtsMat) + R_pp
junk2 = wtsMat %*% R_op  + t(wtsMat %*% R_op)
cbind(diag(junk1), diag(junk2))
