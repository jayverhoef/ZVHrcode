sec_path = 'Rcode/Chapter12/Section 12.6/'
setwd(paste0(SLEDbook_path,sec_path))

source('simSGLM_wExpl.R')
source('addBreakColorLegend.R')

library(viridis)
library(classInt)
library(spmodel)

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
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Estimate covariance parameters and fixed effects using Laplace
#   as implemented in spmodel
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

poismod <- spglm(y ~ 1, family = "poisson", data = sim_data, 
	spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, 
	control = list(reltol = 1e-8))	
summary(poismod)
	
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Plot true values versus spatial random effects (w's) after estimating
#       covariance parameters
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

plot(sim_data$w_true, fitted(poismod, type = 'link'))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Create values on a grid for the likelihood surface
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

coef(poismod, type = "spcov")

# make a graph for the likelihood surface for partial sill and range
# centered on estimated values and scale by estimated value/5
theta1 = log(coef(poismod, type = "spcov")['de']) + (-15:15)/5
theta2 = log(coef(poismod, type = "spcov")['range']) + (-15:15)/5
# matrix to hold likelihood surface values
llgrid = matrix(NA, ncol = 3, nrow = length(theta1)*length(theta2))

# loop through various parameter values and compute -2*loglikelihood
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

