sec_path = 'Rcode/Chapter11/Section spglm/'
setwd(paste0(SLEDbook_path,sec_path))

source('simSGLM_wExpl.R')
source('autocorr_functions.R')
source('logLik_Laplace.R')
source('estpred.R')
source('addBreakColorLegend.R')

library(viridis)
library(classInt)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Simulate SGLM binomial data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#set seed so reproducible
set.seed(2014)

# simulate data from spatial generalized linear model
# betas[1] is overall mean (on log scale)
# betas[2] is regression coefficient on simulated N(0,1) explanatory variables
# betas[3] is regression coefficient on simulated binary variables
# betas[4] is regression coefficient on normal*binary variable interaction
# gammas[1] is partial sill
# gammas[2] is range
# gammas[3] is nugget
#betas = c(-.5, .5, -.5, .5)
betas = c(0, 0, 0, 0)
gammas = c(5, .5, 0.0001)
#undebug(simSGLM_wExpl)
sim_data = simSGLM_wExpl(20^2, autocor_fun = rho_exp, betas = betas,
	gammas = gammas, loc_type = 'grid', pred = FALSE, family = 'binomial')

layout(matrix(1:4, nrow = 1), widths = rep(c(3,1), times = 2))

brks_cex = 1.5
mtext_cex = 3.5
leg_right = 0.6
adj = -0.21
padj = -0.3
cex_plot = 4

cip = classIntervals(sim_data[,'y'], 2, style = 'fisher')
palp = viridis(2)
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

# create design matrix for observed data
X = model.matrix(~ 1, data = sim_data[sim_data$obspred == 'obs',])
#X = model.matrix(~ 1, data = sim_data[sim_data$obspred == 'obs',])
# get observed values
y = sim_data[sim_data$obspred == 'obs','y']
# get distances among observed data locations
distmat = as.matrix(dist(sim_data[sim_data$obspred == 
	'obs',c('xcoord', 'ycoord')]))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Estimate covariance parameters using Laplace
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#initial value for optim
theta = c(-2,-2,-2)
# optimize for covariance parameters
# undebug(logLik_Laplace)
maxvar = 10
maxrange = 5*max(distmat)
# undebug(logLik_Laplace)
start_time = Sys.time()
optout = optim(theta, logLik_Laplace,
	y = y, X = X, distmat = distmat, autocor_fun = rho_exp,
	maxvar = maxvar, maxrange = maxrange, family = 'binomial')
end_time = Sys.time()
difftime(end_time, start_time, units = 'secs')
# covariance parameters
maxvar*expit(optout$par[1])
maxvar*expit(optout$par[2])
maxrange*expit(optout$par[3])
# set theta as the optimized parameters on log scale
theta = optout$par
	
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Estimate fixed effects and spatial random effects after estimating
#       covariance parameters
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

beta_out = estpred(optout$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	family = 'binomial' )
beta_out$betahat

plot(sim_data[,'w_true'], beta_out$w)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Create values on a grid for the likelihood surface
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# make a graph for the likelihood surface for partial sill and range
# center on estimated values and scale by estimated value/10
theta2 = optout$par[2] + (-15:15)/5
theta3 = optout$par[3] + (-15:15)/5
# matrix to hold likelihood surface values
llgrid = matrix(NA, ncol = 3, nrow = length(theta2)*length(theta3))

# loop through various parameter values and compute -2*loglikelihood
iter = 0
i = 1
j = 15
for(i in 1:length(theta2)) {
	for(j in 1:length(theta3)) {
		iter = iter + 1
		llgrid[iter,1] = theta2[i]
		llgrid[iter,2] = theta3[j]
		llgrid[iter,3] = logLik_Laplace(
			theta = c(optout$par[1], theta2[i],theta3[j]), 
			y = y, X = X, distmat = distmat, autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, family = 'binomial')
	}
}

llgrid_binom = llgrid
#  save(llgrid_binom, file = 'llgrid_binom.rda')
#  load('llgrid_binom.rda', verbose = TRUE)
#  llgrid = llgrid_binom

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Figure of simulated data, spatial random effect and their prediction,
#                and the likelihood surface
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/binom_likelihood_estimation"

pdf(paste0(file_name,'.pdf'), width = 11, height = 9.2)

brks_cex = 1.5
mtext_cex = 3.5
leg_right = 0.6
adj = -0.21
padj = -0.3
cex_plot = 4

layout(matrix(1:8, nrow = 2, byrow = TRUE), widths = rep(c(3,1), times = 4))

brks = c(-0.01, 0.01, 1.01)
cip = classIntervals(y, n = 2, style = 'fixed', fixedBreaks= brks)
palp = viridis(2)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(sim_data$xcoord, sim_data$ycoord, col = cip_colors, pch = 15, 
	cex = cex_plot, cex.lab = 2, cex.axis = 1.5, xlab = 'x-coordinate',
	ylab = 'y-coordinate')
mtext('A', cex = mtext_cex, adj = adj, padj = padj)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .35, xright = .2, ytop = .55,
  breaks = brks, colors = palp, cex = brks_cex, printFormat = "1.0")

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
	cex.lab = 2, cex.axis = 1.5, xlab = 'logit(partial sill)',
	ylab = 'logit(range)')
points(optout$par[2], optout$par[3], pch = 19, col = 'white', cex = 1.5)
#points(0, 0, pch = 19, col = 'white', cex = 3)
mtext('C', cex = mtext_cex, adj = adj, padj = padj)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = "1.1")

cip = classIntervals(beta_out$w, 9, style = 'fisher')
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
