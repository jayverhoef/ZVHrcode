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
#   Simulate SGLM Poisson data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#set seed so reproducible
set.seed(1013)

# simulate data from spatial generalized linear model
# betas[1] is overall mean (on log scale)
# betas[2] is regression coefficient on simulated N(0,1) explanatory variables
# betas[3] is regression coefficient on simulated binary variables
# betas[4] is regression coefficient on normal*binary variable interaction
# gammas[1] is partial sill
# gammas[2] is range
# gammas[3] is nugget
betas = c(-.5, .5, -.5, .5)
#betas = c(5, 0, 0, 0)
gammas = c(1, 1, 0.0001, .3)
sim_data = simSGLM_wExpl(20^2, autocor_fun = rho_exp, betas = betas,
	gammas = gammas, loc_type = 'grid', pred = FALSE, family = 'negbinomial')
sqrt(var(sim_data$y))

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

# create design matrix for observed data
X = model.matrix(~ 1, data = sim_data[sim_data$obspred == 'obs',])
X = model.matrix(~ x_1 + x_2 + x_1:x_2, 
	data = sim_data[sim_data$obspred == 'obs',])
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

#initial covariance parameters values for optim
theta = rep(-2, times = 3)
thetanb = c(theta, -7)
# undebug(logLik_Laplace)
# optimize for covariance parameters
maxvar = 2*sqrt(var(sim_data$y))
maxrange = 5*max(distmat)
maxphi = 1000
# undebug(logLik_Laplace)
# undebug(logLik_Laplace)
optout1 = optim(theta, logLik_Laplace, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	family = 'poisson', mlmeth = 'ml')
optout2 = optim(thetanb, logLik_Laplace, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, family = 'negbinomial', mlmeth = 'ml')
optout3 = optim(theta, logLik_Laplace, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	family = 'poisson', mlmeth = 'reml')
optout4 = optim(thetanb, logLik_Laplace, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	maxphi = maxphi, family = 'negbinomial', mlmeth = 'reml')
# covariance parameters
round(cbind(
c(
maxvar*expit(optout1$par[1]),
maxvar*expit(optout1$par[2]),
maxrange*expit(optout1$par[3]),
NA),
c(
maxvar*expit(optout2$par[1]),
maxvar*expit(optout2$par[2]),
maxrange*expit(optout2$par[3]),
maxphi*expit(optout2$par[4])),
c(
maxvar*expit(optout3$par[1]),
maxvar*expit(optout3$par[2]),
maxrange*expit(optout3$par[3]),
NA),
c(
maxvar*expit(optout4$par[1]),
maxvar*expit(optout4$par[2]),
maxrange*expit(optout4$par[3]),
maxphi*expit(optout4$par[4]))
),4)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Estimate fixed effects and spatial random effects after estimating
#       covariance parameters
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

beta_out1 = estpred(optout1$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	family = 'poisson')
beta_out2 = estpred(optout2$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, maxphi = 1000,
	family = 'negbinomial' )
beta_out3 = estpred(optout3$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, 
	family = 'poisson')
beta_out4 = estpred(optout4$par, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, maxphi = 1000,
	family = 'negbinomial' )
beta_out1$betahat
beta_out2$betahat
beta_out3$betahat
beta_out4$betahat

layout(matrix(1:4, nrow = 2, byrow = TRUE))
	plot(sim_data$w_true, beta_out1$w, pch = 19)
	plot(sim_data$w_true, beta_out2$w, pch = 19)
	plot(sim_data$w_true, beta_out3$w, pch = 19)
	plot(sim_data$w_true, beta_out4$w, pch = 19)
layout(1)

cbind(beta_out1$w, beta_out2$w, beta_out3$w, beta_out4$w)

layout(matrix(1:8, nrow = 2, byrow = TRUE), widths = rep(c(3,1), times = 4))

cip = classIntervals(beta_out1$w, 9, style = 'fisher')
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
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = "1.2")

cip = classIntervals(beta_out2$w, 9, style = 'fisher')
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

cip = classIntervals(beta_out3$w, 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(sim_data$xcoord, sim_data$ycoord, col = cip_colors, pch = 15, 
	cex = cex_plot, cex.lab = 2, cex.axis = 1.5, xlab = 'x-coordinate',
	ylab = 'y-coordinate')
mtext('C', cex = mtext_cex, adj = adj, padj = padj)
par(mar = c(0,0,0,0))
plot(c(0,leg_right),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
  xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
  breaks = cip$brks, colors = palp, cex = brks_cex, printFormat = "1.2")

cip = classIntervals(beta_out4$w, 9, style = 'fisher')
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
for(i in 1:length(theta2)) {
	for(j in 1:length(theta3)) {
		iter = iter + 1
		llgrid[iter,1] = theta2[i]
		llgrid[iter,2] = theta3[j]
		llgrid[iter,3] = logLik_Laplace(
			theta = c(optout$par[1], theta2[i],theta3[j]), 
			y = y, X = X, distmat = distmat, autocor_fun = rho_exp, maxvar = maxvar, maxrange = maxrange, family = 'poisson')
	}
}

llgrid_pois = llgrid
#  save(llgrid_pois, file = 'llgrid_pois.rda')
#  load('llgrid_pois.rda', verbose = TRUE)
#  llgrid = llgrid_pois

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
cip = classIntervals(llgrid[,3], n = 9, style = 'fisher')
palp = viridis(9)
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(llgrid[,1:2], col = cip_colors, pch = 15, cex = cex_plot,
	cex.lab = 2, cex.axis = 1.5, xlab = 'log(partial sill)',
	ylab = 'log(range)')
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


