sec_path = 'Rcode/Chapter11/Section 11.4/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Get the Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# attach data library
library(ZVHdata)
library(sp)
library(viridis)
library(classInt)
library(colorspace)
library(spdep)
library(spmodel)
library(xtable)


# load data for graphics and analysis
data(caribouDF)


################################################################################
#-------------------------------------------------------------------------------
#                     Section 11.3.1  Simulated Data
#-------------------------------------------------------------------------------
################################################################################

# get parameters from caribou data set using geostatistical exponential model
geo_exp = splm(z ~ water*tarp, data = caribouDF,  
	xcoord = 'x', ycoord = 'y', 
	spcov_type = 'exponential', estmethod = 'reml')
summary(geo_exp)

set.seed(102)
sim1D = sprnorm(spcov_params = coef(geo_exp, type = "spcov"),
	data = data.frame(x = 1:1000, y = rep(1, times = 1000)),
	xcoord = x, ycoord = y)

# randomized treatments
set.seed(2002)
trt_ran = sample(rep(1:6, times = 5))

d3 = data.frame(x = 1:30, y = rep(1, times = 30), 
	trt = as.factor(trt_ran), z = sim1D[901:930])
sim_expran = splm(z ~ trt, data = d3, xcoord = x, ycoord = y,
	spcov_type = 'exponential')
summary(sim_expran)
sim_expr50 = splm(z ~ trt, data = d3, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = 200, ie = 1e-5, de = 75,
		known = c('range')) )
summary(sim_expr50)

################################################################################
#-------------------------------------------------------------------------------
#      11.4  BLUPTE for the simulated data and the simulation study
#-------------------------------------------------------------------------------
################################################################################

ell = c(1, 0, 0, 0, 0, 0)
BLUPTE4394 = ell %*% coef(sim_expran) + mean(resid(sim_expran))
seBLUPTE4394 = sqrt(t(ell - apply(model.matrix(z ~ trt, data = d3),2,mean)) %*% 
	vcov(sim_expran) %*% (ell - apply(model.matrix(z ~ trt, data = d3),2,mean)))
BLUPTE4394
seBLUPTE4394
BLUPTE4394 - 1.96*seBLUPTE4394
BLUPTE4394 + 1.96*seBLUPTE4394

BLUPTE50 = ell %*% coef(sim_expr50) + mean(resid(sim_expr50))
seBLUPTE50 = sqrt(t(ell - apply(model.matrix(z ~ trt, data = d3),2,mean)) %*% 
	vcov(sim_expr50) %*% (ell - apply(model.matrix(z ~ trt, data = d3),2,mean)))
BLUPTE50
seBLUPTE50
BLUPTE50 - 1.96*seBLUPTE50
BLUPTE50 + 1.96*seBLUPTE50

# same Simulation Experiment as in Section 11.3.3

set.seed(4001)
trt_ran = sample(rep(1:6, times = 5))
i  = 1
ind_mspe_true0 = rep(NA, times = 971)
exp_mspe_true0 = rep(NA, times = 971)
ind_T1err_true0 = rep(NA, times = 971)
exp_T1err_true0 = rep(NA, times = 971)
ind_mspe_rlzdmean = rep(NA, times = 971)
exp_mspe_rlzdmean = rep(NA, times = 971)
ind_T1err_rlzdmean = rep(NA, times = 971)
exp_T1err_rlzdmean = rep(NA, times = 971)
ind_mspe_cntrst = rep(NA, times = 971)
exp_mspe_cntrst = rep(NA, times = 971)
ind_T1err_cntrst = rep(NA, times = 971)
exp_T1err_cntrst = rep(NA, times = 971)
ell = c(1, 0, 0, 0, 0, 0)
BLUPTE_mspe_rlzdmean = rep(NA, times = 971)
BLUPTE_T1err_rlzdmean = rep(NA, times = 971)
for(i in 1:971) {
	dtemp = data.frame(x = 1:30, y = rep(1, times = 30), 
		trt = as.factor(sample(rep(1:6, times = 5))), 
		z = sim1D[i:(i + 29)])
	ind_temp = splm(z ~ trt, data = dtemp, xcoord = x, ycoord = y,
		spcov_type = 'none')
	exp_temp = splm(z ~ trt, data = dtemp, xcoord = x, ycoord = y,
		spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = coef(geo_exp, type = 'spcov')['range'],
		de = coef(geo_exp, type = 'spcov')['de'],
		ie = coef(geo_exp, type = 'spcov')['ie'],
		known = c('range', 'de', 'ie')))
	ind_mspe_true0[i] = coef(ind_temp)[1]^2
	exp_mspe_true0[i] = coef(exp_temp)[1]^2
	ind_T1err_true0[i] = (abs(coef(ind_temp)[1])/
		summary(ind_temp)$coefficients$fixed$Std_Error[1]) > 2.063
	exp_T1err_true0[i] =  (abs(coef(exp_temp)[1])/
		summary(exp_temp)$coefficients$fixed$Std_Error[1]) > 1.96
	rlzdmean = mean(sim1D[i:(i + 29)])
	ind_mspe_rlzdmean[i] = (coef(ind_temp)[1] - rlzdmean)^2
	exp_mspe_rlzdmean[i] = (coef(exp_temp)[1] - rlzdmean)^2
	ind_T1err_rlzdmean[i] = (abs(coef(ind_temp)[1] - rlzdmean)/
		summary(ind_temp)$coefficients$fixed$Std_Error[1]) > 2.063
	exp_T1err_rlzdmean[i] = (abs(coef(exp_temp)[1] - rlzdmean)/
		summary(exp_temp)$coefficients$fixed$Std_Error[1]) > 1.96
	ind_mspe_cntrst[i] = coef(ind_temp)[2]^2
	exp_mspe_cntrst[i] = coef(exp_temp)[2]^2
	ind_T1err_cntrst[i] = (abs(coef(ind_temp)[2])/
		summary(ind_temp)$coefficients$fixed$Std_Error[2]) > 2.063
	exp_T1err_cntrst[i] = (abs(coef(exp_temp)[2])/
		summary(exp_temp)$coefficients$fixed$Std_Error[2]) > 1.96  
  BLUPTE_mspe_rlzdmean[i] = (ell %*% coef(exp_temp) + mean(resid(exp_temp))
		- rlzdmean)^2
  BLUPTE_T1err_rlzdmean[i] = (abs(ell %*% coef(exp_temp) + 
		mean(resid(exp_temp)) - rlzdmean)/
		sqrt(t(ell - apply(model.matrix(z ~ trt, 
		data = dtemp),2,mean)) %*% vcov(exp_temp) %*% 
		(ell - apply(model.matrix(z ~ trt, data = dtemp),2,mean)))) > 1.96
}

sqrt(mean(BLUPTE_mspe_rlzdmean))
mean(BLUPTE_T1err_rlzdmean)


#BLUPTE values to be compared to those values in Table 11.3
sqrt(mean(ind_mspe_rlzdmean))
sqrt(mean(exp_mspe_rlzdmean))

# Values for Table 11.3
sqrt(mean(ind_mspe_true0))
sqrt(mean(exp_mspe_true0))

sqrt(mean(ind_mspe_cntrst))
sqrt(mean(exp_mspe_cntrst))

mean(ind_T1err_rlzdmean)
mean(exp_T1err_rlzdmean)

mean(ind_T1err_true0)
mean(exp_T1err_true0)

mean(ind_T1err_cntrst)
mean(exp_T1err_cntrst)
