sec_path = 'Rcode/Chapter12/Section 12.6/'
setwd(paste0(SLEDbook_path,sec_path))

source('simSGLM_wExpl.R')
source('addBreakColorLegend.R')

library(viridis)
library(classInt)
library(spmodel)
library(pdist)
library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Simulate SGLM Poisson data and estimate parameters and make predictions
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
# exponential autocorrelation function
rho_exp = function(H,gamma_2)
{
	exp(-H/gamma_2)
}

betas = c(2, 0, 0, 0)
gammas = c(1, 1, 0.0001)
sim_data = simSGLM_wExpl(20^2, autocor_fun = rho_exp, betas = betas,
	gammas = gammas, loc_type = 'grid', pred = FALSE)
#save(sim_data, file = 'sim_data.rda')

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
par_orig = par(mar = c(5,5,5,1))
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

par(par_orig)
layout(1)

#-------------------------------------------------------------------------------
#   Estimate covariance parameters and fixed effects using Laplace
#   as implemented in spmodel
#-------------------------------------------------------------------------------

poismod <- spglm(y ~ 1, family = "poisson", data = sim_data, 
	spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, 
	control = list(reltol = 1e-8))	
summary(poismod)
	
#-------------------------------------------------------------------------------
#   Plot true values versus spatial random effects (w's) after estimating
#       covariance parameters
#-------------------------------------------------------------------------------

plot(sim_data$w_true, fitted(poismod, type = 'link'))

#-------------------------------------------------------------------------------
#   Create values on a grid for the likelihood surface
#-------------------------------------------------------------------------------

coef(poismod, type = "spcov")

# make a graph for the likelihood surface for partial sill and range
# centered on estimated values and scale by estimated value/5
theta1 = log(coef(poismod, type = "spcov")['de']) + (-15:15)/5
theta2 = log(coef(poismod, type = "spcov")['range']) + (-15:15)/5
# matrix to hold likelihood surface values
llgrid = matrix(NA, ncol = 3, nrow = length(theta1)*length(theta2))

# loop through various parameter values and compute -2*loglikelihood
# this takes a while; if you don't want to wait, use
# load('llgrid_pois.rda', verbose = TRUE)
iter = 0
for(i in 1:length(theta1)) {
	for(j in 1:length(theta2)) {
		cat("\r", "i: ", i, " out of ", length(theta1),"   j: ", 
			j, " out of ", length(theta1))
		iter = iter + 1
		llgrid[iter,1] = theta1[i]
		llgrid[iter,2] = theta2[j]
		llgrid[iter,3] = -2*logLik(spglm(y ~ 1, family = "poisson", data = sim_data, 
			spcov_initial = spcov_initial("exponential", de = exp(theta1[i]), 
				range = exp(theta2[j]), ie = coef(poismod, type = "spcov")['ie'],
				known = 'given'),
			xcoord = xcoord, ycoord = ycoord))	
	}
}

llgrid_pois = llgrid
#  save(llgrid_pois, file = 'llgrid_pois.rda')
#  load('llgrid_pois.rda', verbose = TRUE)
#  llgrid = llgrid_pois

#-------------------------------------------------------------------------------
#   Figure of simulated data, spatial random effect and their prediction,
#                and the likelihood surface
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
#cip = classIntervals(llgrid[,3], n = 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(llgrid[,1:2], col = cip_colors, pch = 15, cex = cex_plot,
	cex.lab = 2, cex.axis = 1.5, xlab = 'log(partial sill)',
	ylab = 'log(range)')
points(log(coef(poismod, type = "spcov")['de']), 
	log(coef(poismod, type = "spcov")['range']), 
	pch = 19, col = 'white', cex = 1.5)
#points(0, 0, pch = 19, col = 'white', cex = 3)
mtext('C', cex = mtext_cex, adj = adj, padj = padj)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = "1.1")

cip = classIntervals(fitted(poismod, type = 'link'), 9, style = 'fisher')
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
#        and making predictions
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
# The simulation takes a long time, so you can load the stored data instead
# load('store.rda')

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
#  restore the stored data if you want to skip the simulation
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
	MSE = sqrt(MSE/niter),
	cover_est_uncor = cover_est_uncor/niter,
	cover_est_cor = cover_est_cor/niter
)
(est_bias/niter)^2/(MSE/niter)
MSPE = MSPE/niter
sglm_simsumm = rbind(sglm_fe,
	c(pred_bias/niter, sqrt(MSPE), cover_pred_uncor/niter, cover_pred_cor/niter))

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
