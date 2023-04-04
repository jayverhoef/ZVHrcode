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
#       11.5.1 Optimal Spatial Design of Experiments, Simulated Data
#-------------------------------------------------------------------------------
################################################################################

#  Optimal Experimental Design

d1 = data.frame(x = 1:30, y = rep(1, times = 30), 
	trt = as.factor(rep(1:6, times = 5)), z = sim1D[901:930])
d2 = data.frame(x = 1:30, y = rep(1, times = 30), 
	trt = as.factor(kronecker(1:6,rep(1,times = 5))), z = sim1D[901:930])
# fit uncorrelated models to all 3 designs
sim_indran = splm(z ~ trt, data = d3, xcoord = x, ycoord = y,
	spcov_type = 'none' )
sim_exprand = splm(z ~ trt, data = d3, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential',
	de = coef(geo_exp, type = "spcov")['de'],
	ie = coef(geo_exp, type = "spcov")['ie'],
	range = coef(geo_exp, type = "spcov")['range'], 
	known = c('ie','de','range')) )
sim_exp1 = splm(z ~ trt, data = d1, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential',
	de = coef(geo_exp, type = "spcov")['de'],
	ie = coef(geo_exp, type = "spcov")['ie'],
	range = coef(geo_exp, type = "spcov")['range'], 
	known = c('ie','de','range')) )
sim_exp2 = splm(z ~ trt, data = d2, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential',
	de = coef(geo_exp, type = "spcov")['de'],
	ie = coef(geo_exp, type = "spcov")['ie'],
	range = coef(geo_exp, type = "spcov")['range'], 
	known = c('ie','de','range')) )
summary(sim_indran)
summary(sim_exprand)
summary(sim_exp1)
summary(sim_exp2)
ell = c(1, 0, 0, 0, 0, 0)
BLUPTE_exprand = ell %*% coef(sim_exprand) + mean(resid(sim_exprand))
seBLUPTE_exprand = sqrt(t(ell - apply(model.matrix(z ~ trt, data = d3),2,mean)) %*% 
	vcov(sim_exprand) %*% (ell - apply(model.matrix(z ~ trt, data = d3),2,mean)))
BLUPTE_exp1 = ell %*% coef(sim_exp1) + mean(resid(sim_exp1))
seBLUPTE_exp1 = sqrt(t(ell - apply(model.matrix(z ~ trt, data = d1),2,mean)) %*% 
	vcov(sim_exp1) %*% (ell - apply(model.matrix(z ~ trt, data = d1),2,mean)))
BLUPTE_exp2 = ell %*% coef(sim_exp2) + mean(resid(sim_exp2))
seBLUPTE_exp2 = sqrt(t(ell - apply(model.matrix(z ~ trt, data = d2),2,mean)) %*% 
	vcov(sim_exp2) %*% (ell - apply(model.matrix(z ~ trt, data = d2),2,mean)))

sim_design = rbind(
	c(coef(sim_indran)[1], BLUPTE_exprand, BLUPTE_exp1, BLUPTE_exp2, 
		summary(sim_indran)$coefficients$fixed$Std_Error[1], 
		seBLUPTE_exprand, seBLUPTE_exp1, seBLUPTE_exp2),
	c(coef(sim_indran)[1], coef(sim_exprand)[1], coef(sim_exp1)[1], coef(sim_exp2)[1], 
		summary(sim_indran)$coefficients$fixed$Std_Error[1],
		summary(sim_exprand)$coefficients$fixed$Std_Error[1],
		summary(sim_exp1)$coefficients$fixed$Std_Error[1], 
		summary(sim_exp2)$coefficients$fixed$Std_Error[1]),
	c(coef(sim_indran)[2], coef(sim_exprand)[2], coef(sim_exp1)[2], coef(sim_exp2)[2], 
		summary(sim_indran)$coefficients$fixed$Std_Error[2],
		summary(sim_exprand)$coefficients$fixed$Std_Error[2],
		summary(sim_exp1)$coefficients$fixed$Std_Error[2], 
		summary(sim_exp2)$coefficients$fixed$Std_Error[2])
)

print(
    xtable(sim_design, 
      align = c('l',rep('l', times = length(sim_design[1,]))),
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

1.063/0.496


sim_exp2_fitted = splm(z ~ trt, data = d2, xcoord = x, ycoord = y,
	spcov_type = 'exponential')
summary(sim_exp2_fitted)

ell %*% coef(sim_exp2_fitted) + mean(resid(sim_exp2_fitted))
sqrt(t(ell - apply(model.matrix(z ~ trt, data = d2),2,mean)) %*% 
	vcov(sim_exp2_fitted) %*% (ell - apply(model.matrix(z ~ trt, data = d2),2,mean)))

(ell %*% coef(sim_exp2_fitted) + mean(resid(sim_exp2_fitted)))/
sqrt(t(ell - apply(model.matrix(z ~ trt, data = d2),2,mean)) %*% 
	vcov(sim_exp2_fitted) %*% (ell - apply(model.matrix(z ~ trt, data = d2),2,mean)))

