sec_path = 'Rcode/Chapter12/Section 12.6/'
setwd(paste0(SLEDbook_path,sec_path))

source('simSGLM_wExpl.R')
source('estpred_4Mike.R')
library(pdist)
library(spmodel)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Define a few functions first
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# exponential autocorrelation function
rho_exp = function(H,gamma_2)
{
	exp(-H/gamma_2)
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#   Testing prediction standard errors in spmodel
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#set seed so reproducible
set.seed(11033)

betas = c(0.5, 0.5, -0.5, 0.5)
gammas = c(1, 1, 0.0001)

sim_data = simSGLM_wExpl(200, autocor_fun = rho_exp, 
	betas = betas, gammas = gammas, 
	loc_type = 'random', pred = TRUE)
sim_data_predNA = sim_data
sim_data_predNA[sim_data_predNA$obspred == 'pred', 'y'] = NA

poismod <- spglm(y ~ x_1*x_2, family = "poisson", data = sim_data_predNA, 
	spcov_type = "exponential", xcoord = xcoord, ycoord = ycoord, 
	control = list(reltol = 1e-8))	
predadj = predict(poismod, se.fit = TRUE)

# ------ now use my code with covariance parameter estimates from spglm --------

# coordinate data and sample sizes
xyobs = sim_data[sim_data$obspred == 'obs',c('xcoord','ycoord')]
xypred = sim_data[sim_data$obspred == 'pred',c('xcoord','ycoord')]
npred = dim(xypred)[1]
nobs = dim(xyobs)[1]

# get observed values
y = sim_data[sim_data$obspred == 'obs','y']
# create design matrix for observed data
X = model.matrix(~ x_1*x_2, data = sim_data[sim_data$obspred == 'obs',])
Xp = model.matrix(~ x_1*x_2, data = sim_data[sim_data$obspred == 'pred',])

# get distances among observed data and prediction locations
distmat = as.matrix(dist(xyobs))
dist_op = as.matrix(pdist(xyobs, xypred))
dist_pp = as.matrix(dist(xypred))

# create covariance parameter vector for my function
theta4jay = coef(poismod, type = 'spcov')[c('de', 'range', 'ie')]

# simplified version of my function for doing estimation and prediction
predjay = estpred_4Mike(theta4jay, y = y, X = X, distmat = distmat, 
	autocor_fun = rho_exp, family = 'poisson', Xp = Xp, 
	dist_op = dist_op, dist_pp = dist_pp )

# are our w's the same?
cbind(predjay$w, fitted(poismod, type = 'link'))
# are our standard errors of fixed effects the same?
cbind(sqrt(diag(predjay$covbeta_adj)), sqrt(diag(vcov(poismod))))
# are our prediction standard errors the same?
cbind(predjay$wpred_se_adj, predadj$se.fit)

