setwd('/mnt/ExtraDrive1/Work/desktop_data/2022_papers/hglmm/Github/hglmm')
source('simSGLM_wExpl.R')
source('autocorr_functions.R')
source('logLik_Laplace.R')
source('estpred.R')
source('addBreakColorLegend.R')

library(viridis)
library(classInt)
library(pdist)
library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Simulation to test for bias and coverage when estimating fixed effects
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#set seed so reproducible
set.seed(17877)

betas = c(0.5, 0.5, -0.5, 0.5)
gammas = c(1, 1, 0.0001)
niter = 2000
store = vector(mode='list', length = niter)
start_time = Sys.time()
for(iter in 1:niter) {
	cat("\r", "iteration: ", iter)
	sim_data = simSGLM_wExpl(200, autocor_fun = rho_exp, 
		betas = betas, gammas = gammas, 
		loc_type = 'random', pred = TRUE)
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
	distmat = as.matrix(dist(xyobs))
	dist_op = as.matrix(pdist(xyobs, xypred))
	dist_pp = as.matrix(dist(xypred))

	#initial covariance parameters values for optim
	theta = rep(-2, times = 3)
	# undebug(logLik_Laplace)
	# optimize for covariance parameters
	maxvar = 10*var(log(y+1))
	maxrange = 10*max(distmat)
	# undebug(logLik_Laplace)
	optout = optim(theta, logLik_Laplace, y = y, X = X, distmat = distmat, 
		autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
		maxphi = maxphi, family = 'poisson', mlmeth = 'reml')
	# covariance parameters
	maxvar*expit(optout$par[1])
	maxrange*expit(optout$par[2])
	maxvar*expit(optout$par[3])
	# set theta as the optimized parameters on scaled logit transformation
	theta = optout$par

	pred_out = estpred(optout$par, y = y, X = X, distmat = distmat, 
		autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
		family = 'poisson', Xp = Xp, dist_op = dist_op, dist_pp = dist_pp )
  pred_out$betahat
  
  w_true = sim_data$w_true[sim_data$obspred == 'pred']
  bias_pred = mean(pred_out$w_pred - w_true)
  mse_pred = mean((pred_out$w_pred - w_true)^2)
  cover_pred = sum(pred_out$w_pred - 1.645*pred_out$wpred_se < w_true & 
		w_true < pred_out$w_pred + 1.645*pred_out$wpred_se)/100
  cover_pred_adj = sum(pred_out$w_pred - 1.645*pred_out$wpred_se_adj < w_true & 
		w_true < pred_out$w_pred + 1.645*pred_out$wpred_se_adj)/100
		
	plot(c(min(w_true,pred_out$w_pred), max(w_true,pred_out$w_pred)),
		c(min(w_true,pred_out$w_pred), max(w_true,pred_out$w_pred)), type = 'l',
		xlab = 'True w', ylab = 'Predicted w')
	points(w_true, pred_out$w_pred)

	# compare standard errors when using only covbeta versus the corrected one
	store[[iter]] = 
	data.frame(True = betas, betahat = pred_out$betahat, 
		SE_uncorr = sqrt(diag(pred_out$covbeta)),
		CI90_uncorr = pred_out$betahat - 
			1.645*sqrt(diag(pred_out$covbeta)) < betas & 
			betas < pred_out$betahat + 
			1.645*sqrt(diag(pred_out$covbeta)),
		SE_corr = sqrt(diag(pred_out$covbeta_adj)),
		CI90_corr = pred_out$betahat - 
			1.645*sqrt(diag(pred_out$covbeta_adj)) < betas & 
			betas < pred_out$betahat + 
			1.645*sqrt(diag(pred_out$covbeta_adj)),
		bias_pred = bias_pred,
		mse_pred = mse_pred,
		cover_pred = cover_pred,
		cover_pred_adj = cover_pred_adj
	)
}
cat("\n")
end_time = Sys.time()
difftime(end_time, start_time)

store_pois = store
save(store_pois,file = 'store_pois.rda')
#  load('store.rda')

i = 1
est_bias = rep(0, times = 4)
est_mse = rep(0, times = 4)
cover_est_uncor = rep(0, times = 4)
cover_est_cor = rep(0, times = 4)
pred_bias = 0
pred_mse = 0
cover_pred_uncor = 0
cover_pred_cor = 0
for(i in 1:niter) {
	est_bias = est_bias + store[[i]]$betahat - store[[i]]$True
	est_mse = est_mse + (store[[i]]$betahat - store[[i]]$True)^2
	cover_est_uncor = cover_est_uncor + store[[i]]$CI90_uncorr
	cover_est_cor = cover_est_cor + store[[i]]$CI90_corr
	pred_bias = pred_bias + store[[i]]$bias_pred[1]
	pred_mse = pred_mse + store[[i]]$mse_pred[1]
	cover_pred_uncor = cover_pred_uncor + store[[i]]$cover_pred[1]
	cover_pred_cor = cover_pred_cor + store[[i]]$cover_pred_adj[1]
}

sglm_fe = data.frame(est_bias =  est_bias/niter,
	est_mse = est_mse/niter,
	bias2_prop = (est_bias/niter)^2/(est_mse/niter),
	cover_est_uncor = cover_est_uncor/niter,
	cover_est_cor = cover_est_cor/niter
)

pred_mse/niter
sglm_simsumm = rbind(sglm_fe,
	c(pred_bias/niter, pred_mse/niter, (pred_bias/niter)^2/(pred_mse/niter),
	cover_pred_uncor/niter, cover_pred_cor/niter))

print(
    xtable(sglm_simsumm, 
      align = c('l',rep('l', times = length(sglm_fe[1,]))),
      digits = c(0, rep(3, times = 5))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
