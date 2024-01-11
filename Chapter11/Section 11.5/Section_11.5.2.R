sec_path = 'Rcode/Chapter11/Section 11.5/'
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
#                     Section 11.5  
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


################################################################################
#-------------------------------------------------------------------------------
#       11.5.2 Optimal Spatial Design of Experiments, Simulation Experiment
#-------------------------------------------------------------------------------
################################################################################

# Simulation Experiment

set.seed(4010)
ell = c(1, 0, 0, 0, 0, 0)
results_table = matrix(rep(0, times = 24), nrow = 3, ncol = 8)
for(i in 1:971) {
	dtemp1 = data.frame(x = 1:30, y = rep(1, times = 30), 
		trt = as.factor(rep(1:6, times = 5)), z = sim1D[i:(i + 29)])
	dtemp1$trt = relevel(dtemp1$trt, ref = '3')
	dtemp2 = data.frame(x = 1:30, y = rep(1, times = 30), 
		trt = as.factor(kronecker(1:6,rep(1,times = 5))), z = sim1D[i:(i + 29)])
	dtemp2$trt = relevel(dtemp2$trt, ref = '3')
	dtemp3 = data.frame(x = 1:30, y = rep(1, times = 30), 
		trt = as.factor(sample(rep(1:6, times = 5))), 
		z = sim1D[i:(i + 29)])
	dtemp3$trt = relevel(dtemp3$trt, ref = '3')
	rlzdmean = mean(sim1D[i:(i + 29)])
	
	ind_temp = splm(z ~ trt, data = dtemp3, xcoord = x, ycoord = y,
		spcov_type = 'none')
	exp_temp1 = splm(z ~ trt, data = dtemp1, xcoord = x, ycoord = y,
		spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = coef(geo_exp, type = 'spcov')['range'],
		de = coef(geo_exp, type = 'spcov')['de'],
		ie = coef(geo_exp, type = 'spcov')['ie'],
		known = c('range', 'de', 'ie')))
	exp_temp2 = splm(z ~ trt, data = dtemp2, xcoord = x, ycoord = y,
		spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = coef(geo_exp, type = 'spcov')['range'],
		de = coef(geo_exp, type = 'spcov')['de'],
		ie = coef(geo_exp, type = 'spcov')['ie'],
		known = c('range', 'de', 'ie')))
	exp_temp3 = splm(z ~ trt, data = dtemp3, xcoord = x, ycoord = y,
		spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = coef(geo_exp, type = 'spcov')['range'],
		de = coef(geo_exp, type = 'spcov')['de'],
		ie = coef(geo_exp, type = 'spcov')['ie'],
		known = c('range', 'de', 'ie')))
		
	results_table[1,1] = results_table[1,1] + (coef(ind_temp)[1] - rlzdmean)^2
	results_table[1,2] = results_table[1,2] + (ell %*% coef(exp_temp3) + 
		mean(resid(exp_temp3)) - rlzdmean)^2
	results_table[1,3] = results_table[1,3] + (ell %*% coef(exp_temp1) + 
		mean(resid(exp_temp1)) - rlzdmean)^2
	results_table[1,4] = results_table[1,4] + (ell %*% coef(exp_temp2) + 
		mean(resid(exp_temp2)) - rlzdmean)^2
		
	results_table[1,5] = results_table[1,5] + ((abs(coef(ind_temp)[1]- rlzdmean)/
		summary(ind_temp)$coefficients$fixed$Std_Error[1]) > 2.063)*1
	results_table[1,6] = results_table[1,6] + ((abs(ell %*% coef(exp_temp3) + 
		mean(resid(exp_temp3)) - rlzdmean)/
		sqrt(t(ell - apply(model.matrix(z ~ trt, 
		data = dtemp3),2,mean)) %*% vcov(exp_temp3) %*% 
		(ell - apply(model.matrix(z ~ trt, data = dtemp3),2,mean)))) > 1.96)*1
	results_table[1,7] = results_table[1,7] + ((abs(ell %*% coef(exp_temp1) + 
		mean(resid(exp_temp1)) - rlzdmean)/
		sqrt(t(ell - apply(model.matrix(z ~ trt, 
		data = dtemp1),2,mean)) %*% vcov(exp_temp1) %*% 
		(ell - apply(model.matrix(z ~ trt, data = dtemp1),2,mean)))) > 1.96)*1
	results_table[1,8] = results_table[1,8] + ((abs(ell %*% coef(exp_temp2) + 
		mean(resid(exp_temp2)) - rlzdmean)/
		sqrt(t(ell - apply(model.matrix(z ~ trt, 
		data = dtemp2),2,mean)) %*% vcov(exp_temp2) %*% 
		(ell - apply(model.matrix(z ~ trt, data = dtemp2),2,mean)))) > 1.96)*1

	results_table[2,1] = results_table[2,1] + coef(ind_temp)[1]^2
	results_table[2,2] = results_table[2,2] + coef(exp_temp3)[1]^2
	results_table[2,3] = results_table[2,3] + coef(exp_temp1)[1]^2
	results_table[2,4] = results_table[2,4] + coef(exp_temp2)[1]^2
			
	results_table[2,5] = results_table[2,5] +  ((abs(coef(ind_temp)[1])/
		summary(ind_temp)$coefficients$fixed$Std_Error[1]) > 2.063)*1
	results_table[2,6] = results_table[2,6] +  ((abs(coef(exp_temp3)[1])/
		summary(exp_temp3)$coefficients$fixed$Std_Error[1]) > 1.96)*1
	results_table[2,7] = results_table[2,7] +  ((abs(coef(exp_temp1)[1])/
		summary(exp_temp1)$coefficients$fixed$Std_Error[1]) > 1.96)*1
	results_table[2,8] = results_table[2,8] +  ((abs(coef(exp_temp2)[1])/
		summary(exp_temp2)$coefficients$fixed$Std_Error[1]) > 1.96)*1
		
	results_table[3,1] = results_table[3,1] + coef(ind_temp)[3]^2
	results_table[3,2] = results_table[3,2] + coef(exp_temp3)[3]^2
	results_table[3,3] = results_table[3,3] + coef(exp_temp1)[3]^2
	results_table[3,4] = results_table[3,4] + coef(exp_temp2)[3]^2
			
	results_table[3,5] = results_table[3,5] +  ((abs(coef(ind_temp)[3])/
		summary(ind_temp)$coefficients$fixed$Std_Error[3]) > 2.063)*1
	results_table[3,6] = results_table[3,6] +  ((abs(coef(exp_temp3)[3])/
		summary(exp_temp3)$coefficients$fixed$Std_Error[3]) > 1.96)*1
	results_table[3,7] = results_table[3,7] +  ((abs(coef(exp_temp1)[3])/
		summary(exp_temp1)$coefficients$fixed$Std_Error[3]) > 1.96)*1
	results_table[3,8] = results_table[3,8] +  ((abs(coef(exp_temp2)[3])/
		summary(exp_temp2)$coefficients$fixed$Std_Error[3]) > 1.96)*1

}

results_table = results_table/971
results_table[1:3,1:4] = sqrt(results_table[1:3,1:4])

################################################################################
#             Table 11.4
################################################################################

print(
    xtable(results_table, 
      align = c('l',rep('l', times = length(results_table[1,]))),
      digits = c(0,rep(3, times = 4),rep(3, times = 4)),
      caption = 'Estimates and standard errors',
      label = 'tab:SealsFixEff'
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)


