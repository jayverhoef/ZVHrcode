sec_path = 'Rcode/Chapter12/Section 12.6/'
setwd(paste0(SLEDbook_path,sec_path))

source('simSGLM_wExpl.R')
source('autocorr_functions.R')
source('logLik_Laplace.R')
source('estpred.R')
source('addBreakColorLegend.R')

library(viridis)
library(classInt)
library(pdist)
library(spmodel)
library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Simulation to test for bias and coverage when estimating fixed effects
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# a function spatial prediction as if the w's were observed
# takes a fitted model from spmodel using the spglm function

predvar_naive = function(modfit, data)
{
	theta = coef(modfit, type = 'spcov')
	xyobs = sim_data[sim_data$obspred == 'obs',c('xcoord','ycoord')]
	xypred = sim_data[sim_data$obspred == 'pred',c('xcoord','ycoord')]
	distmat = as.matrix(dist(xyobs))
	dist_op = as.matrix(pdist(xyobs, xypred))
	dist_pp = as.matrix(dist(xypred))
	Sigma = theta['de']*exp(-distmat/theta['range']) + 
		theta['ie']*diag(length(modfit$y))
	# covariance matrix between observed and prediction locations
	R_op = theta['de']*exp(-dist_op/theta['range'])
	# covariance matrix among prediction locations
	R_pp = theta['de']*exp(-dist_pp/theta['range']) + 
		theta['ie']*diag(dim(xypred)[1])
	X = model.matrix(modfit)
	form = modfit$formula
	form[[2]] = NULL
	Xp = model.matrix(form, data[data$obspred == 'pred',])
	
	CovMati = solve(Sigma)
	covbeta = solve(t(X) %*% (CovMati %*% X))
	# matrix of partial results to simplify algebra
	d1 = (Xp - t(R_op) %*% (CovMati %*% X))
	# typical prediction (kriging) covariance matrix
	covpred = d1 %*% covbeta %*% t(d1) -
		t(R_op) %*% CovMati %*% R_op + R_pp
	sqrt(diag(covpred))
}

#set seed so reproducible
set.seed(11033)

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
	sim_data_predNA = sim_data
	sim_data_predNA[sim_data_predNA$obspred == 'pred', 'y'] = NA

	poismod <- spglm(y ~ x_1*x_2, family = "poisson", data = sim_data_predNA, 
		spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, 
		control = list(reltol = 1e-8))	

	xyobs = sim_data[sim_data$obspred == 'obs',c('xcoord','ycoord')]
	xypred = sim_data[sim_data$obspred == 'pred',c('xcoord','ycoord')]
	npred = dim(xypred)[1]
	nobs = dim(xyobs)[1]

  vcovnai = vcov(poismod, var_correct = FALSE)
  vcovadj = vcov(poismod) 
  predadj = predict(poismod, se.fit = TRUE)
  predsenai = predvar_naive(poismod, sim_data_predNA) 
  
  w_true = sim_data$w_true[sim_data$obspred == 'pred']
  bias_pred = mean(predadj$fit - w_true)
  MSPE = mean((predadj$fit - w_true)^2) 
  cover_pred_adj = sum(predadj$fit - 1.645*predadj$se.fit < w_true & 
		w_true < predadj$fit + 1.645*predadj$se.fit)/100
  cover_pred_nai = sum(predadj$fit - 1.645*predsenai < w_true & 
		w_true < predadj$fit + 1.645*predsenai)/100
	plot(c(min(w_true, predadj$fit),max(w_true, predadj$fit)),
		c(min(w_true, predadj$fit),max(w_true, predadj$fit)), type = 'l',
		xlab = 'true', ylab = 'predicted')
	points(w_true, predadj$fit)
	# compare standard errors when using only covbeta versus the corrected one
	store[[iter]] = 
	data.frame(True = betas, betahat = coef(poismod), 
		SE_uncorr = sqrt(diag(vcovnai)),
		CI90_uncorr = coef(poismod) - 
			1.645*sqrt(diag(vcovnai)) < betas & 
			betas < coef(poismod) + 
			1.645*sqrt(diag(vcovnai)),
		SE_corr = sqrt(diag(vcovadj)),
		CI90_corr = coef(poismod) - 
			1.645*sqrt(diag(vcovadj)) < betas & 
			betas < coef(poismod) + 
			1.645*sqrt(diag(vcovadj)),
		bias_pred = bias_pred,
		MSPE = MSPE,
		cover_pred = cover_pred_nai,
		cover_pred_adj = cover_pred_adj
	)
}
cat("\n")
end_time = Sys.time()
difftime(end_time, start_time)

#  save(store,file = 'store.rda')
#  load('store.rda')

i = 1
est_bias = rep(0, times = 4)
MSE = rep(0, times = 4)
cover_est_uncor = rep(0, times = 4)
cover_est_cor = rep(0, times = 4)
pred_bias = 0
cover_pred_uncor = 0
cover_pred_cor = 0
MSPE = 0

for(i in 1:niter) {
	est_bias = est_bias + store[[i]]$betahat - store[[i]]$True
	MSE = MSE + (store[[i]]$betahat - store[[i]]$True)^2
	cover_est_uncor = cover_est_uncor + store[[i]]$CI90_uncorr
	cover_est_cor = cover_est_cor + store[[i]]$CI90_corr
	pred_bias = pred_bias + store[[i]]$bias_pred[1]
	cover_pred_uncor = cover_pred_uncor + store[[i]]$cover_pred[1]
	cover_pred_cor = cover_pred_cor + store[[i]]$cover_pred_adj[1]
	MSPE = MSPE + store[[i]]$MSPE[1]
}
sglm_fe = data.frame(est_bias =  est_bias/niter,
	MSE = MSE/niter,
	cover_est_uncor = cover_est_uncor/niter,
	cover_est_cor = cover_est_cor/niter
)
(est_bias/niter)^2/(MSE/niter)
MSPE = MSPE/niter
sglm_simsumm = rbind(sglm_fe,
	c(pred_bias/niter, MSPE, cover_pred_uncor/niter, cover_pred_cor/niter))

print(
    xtable(sglm_simsumm, 
      align = c('l',rep('l', times = length(sglm_fe[1,]))),
      digits = c(0, rep(3, times = 4))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
