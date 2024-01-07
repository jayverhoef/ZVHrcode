#set a path as a working directory
sec_path = 'Rcode/Chapter8/Section 8.2/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)
library(spmodel)

# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 12
locxy <- cbind(rep(1:K,K), rep(1:K, each = K))
n <- K^2
theta <- 0.5

# Define a covariance structure on the grid: an exponential covariance function with correlation theta at
# unit distance
d = as.matrix(dist(cbind(locxy)))
Gtrue = theta^d
tcholGtrue = t(chol(Gtrue))

# ----------------------------- Intercept Only ---------------------------------

mle1 <- matrix(0, 1000, 4)
remle1 <- matrix(0, 1000, 4)
set.seed(506)
for(k in 1:1000){
	cat("\r", "iteration: ", k)
	# simulate data, all true regression coefficients are zero
	y1 <- 0 + tcholGtrue %*% rnorm(n)
	# create data.frame for spmodel
	DF = data.frame(y = y1, xcoord = locxy[,1], ycoord = locxy[,2])
	# ------- ML ESTIMATION
	# fit model using spmodel, holding nugget effect at zero
	fitout = splm(y ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'ml')
	# estimated range parameter, reparameterized from spmodel to rho^distance
	# as given in the book
	mle1[k,1] = exp(-1/coef(fitout, type = 'spcov')['range'])
	# estimated variance parameter
	mle1[k,2] = coef(fitout, type = 'spcov')['de']
	# average estimated 90% confidebce interval width
	mle1[k,3] = mean(summary(fitout)$coefficients$fixed[,'Std_Error']*2*1.645)
	# 90% confidebce interval coverage (of true zero value)
	mle1[k,4] = mean(summary(fitout)$coefficients$fixed[,'p'] > 0.1)
	# ------- REML ESTIMATION
	fitout = splm(y ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'reml')
	remle1[k,1] = exp(-1/coef(fitout, type = 'spcov')['range'])
	remle1[k,2] = coef(fitout, type = 'spcov')['de']
	remle1[k,3] = mean(summary(fitout)$coefficients$fixed[,'Std_Error']*2*1.645)
	remle1[k,4] = mean(summary(fitout)$coefficients$fixed[,'p'] > 0.1)
}

round(c(mean(mle1[,2] - 1.0),
mean(remle1[,2] - 1.0),
mean(mle1[,1] - 0.5),
mean(remle1[,1] - 0.5),
mean((mle1[,2] - 1.0)^2),
mean((remle1[,2] - 1.0)^2),
mean((mle1[,1] - 0.5)^2),
mean((remle1[,1] - 0.5)^2),
mean(mle1[,3]),
mean(remle1[,3]),
mean(mle1[,4]),
mean(remle1[,4])
),3)

# ----------------------------- Row Effects ------------------------------------

mle2 <- matrix(0, 1000, 4)
remle2 <- matrix(0, 1000, 4)
set.seed(607)
for(k in 1:1000){
	cat("\r", "iteration: ", k)
	# simulate data, all true regression coefficients are zero
	y1 <- 0 + tcholGtrue %*% rnorm(n)
	# create data.frame for spmodel
	DF = data.frame(y = y1, xcoord = locxy[,1], ycoord = locxy[,2])
	# ------- ML ESTIMATION
	# fit model using spmodel, holding nugget effect at zero
	fitout = splm(y ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'ml')
	# estimated range parameter
	mle2[k,1] = exp(-1/coef(fitout, type = 'spcov')['range'])
	# estimated variance parameter
	mle2[k,2] = coef(fitout, type = 'spcov')['de']
	# average estimated 90% confidebce interval width
	mle2[k,3] = mean(summary(fitout)$coefficients$fixed[,'Std_Error']*2*1.645)
	# 90% confidebce interval coverage (of true zero value)
	mle2[k,4] = mean(summary(fitout)$coefficients$fixed[,'p'] > 0.1)
	# ------- REML ESTIMATION
	fitout = splm(y ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'reml')
	remle2[k,1] = exp(-1/coef(fitout, type = 'spcov')['range'])
	remle2[k,2] = coef(fitout, type = 'spcov')['de']
	remle2[k,3] = mean(summary(fitout)$coefficients$fixed[,'Std_Error']*2*1.645)
	remle2[k,4] = mean(summary(fitout)$coefficients$fixed[,'p'] > 0.1)
}

round(c(mean(mle2[,2] - 1.0),
mean(remle2[,2] - 1.0),
mean(mle2[,1] - 0.5),
mean(remle2[,1] - 0.5),
mean((mle2[,2] - 1.0)^2),
mean((remle2[,2] - 1.0)^2),
mean((mle2[,1] - 0.5)^2),
mean((remle2[,1] - 0.5)^2),
mean(mle2[,3]),
mean(remle2[,3]),
mean(mle2[,4]),
mean(remle2[,4])
),3)

# ---------------------- Row and Column Effects --------------------------------

mle3 <- matrix(0, 1000, 4)
remle3 <- matrix(0, 1000, 4)
set.seed(508)
for(k in 1:1000){
	cat("\r", "iteration: ", k)
	# simulate data, all true regression coefficients are zero
	y1 <- 0 + tcholGtrue %*% rnorm(n)
	# create data.frame for spmodel
	DF = data.frame(y = y1, xcoord = locxy[,1], ycoord = locxy[,2])
	# ------- ML ESTIMATION
	# fit model using spmodel, holding nugget effect at zero
	fitout = splm(y ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'ml')
	# estimated range parameter
	mle3[k,1] = exp(-1/coef(fitout, type = 'spcov')['range'])
	# estimated variance parameter
	mle3[k,2] = coef(fitout, type = 'spcov')['de']
	# average estimated 90% confidebce interval width
	mle3[k,3] = mean(summary(fitout)$coefficients$fixed[,'Std_Error']*2*1.645)
	# 90% confidebce interval coverage (of true zero value)
	mle3[k,4] = mean(summary(fitout)$coefficients$fixed[,'p'] > 0.1)
	# ------- REML ESTIMATION
	fitout = splm(y ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'reml')
	remle3[k,1] = exp(-1/coef(fitout, type = 'spcov')['range'])
	remle3[k,2] = coef(fitout, type = 'spcov')['de']
	remle3[k,3] = mean(summary(fitout)$coefficients$fixed[,'Std_Error']*2*1.645)
	remle3[k,4] = mean(summary(fitout)$coefficients$fixed[,'p'] > 0.1)
}

round(c(mean(mle3[,2] - 1.0),
mean(remle3[,2] - 1.0),
mean(mle3[,1] - 0.5),
mean(remle3[,1] - 0.5),
mean((mle3[,2] - 1.0)^2),
mean((remle3[,2] - 1.0)^2),
mean((mle3[,1] - 0.5)^2),
mean((remle3[,1] - 0.5)^2),
mean(mle3[,3]),
mean(remle3[,3]),
mean(mle3[,4]),
mean(remle3[,4])
),3)


simDF = 
rbind(
	c(mean(mle1[,2] - 1.0),
		mean(remle1[,2] - 1.0),
		mean(mle1[,1] - 0.5),
		mean(remle1[,1] - 0.5),
		mean((mle1[,2] - 1.0)^2),
		mean((remle1[,2] - 1.0)^2),
		mean((mle1[,1] - 0.5)^2),
		mean((remle1[,1] - 0.5)^2),
		mean(mle1[,3]),
		mean(remle1[,3]),
		mean(mle1[,4]),
		mean(remle1[,4])
	),
	c(mean(mle2[,2] - 1.0),
		mean(remle2[,2] - 1.0),
		mean(mle2[,1] - 0.5),
		mean(remle2[,1] - 0.5),
		mean((mle2[,2] - 1.0)^2),
		mean((remle2[,2] - 1.0)^2),
		mean((mle2[,1] - 0.5)^2),
		mean((remle2[,1] - 0.5)^2),
		mean(mle2[,3]),
		mean(remle2[,3]),
		mean(mle2[,4]),
		mean(remle2[,4])
	),
	c(mean(mle3[,2] - 1.0),
		mean(remle3[,2] - 1.0),
		mean(mle3[,1] - 0.5),
		mean(remle3[,1] - 0.5),
		mean((mle3[,2] - 1.0)^2),
		mean((remle3[,2] - 1.0)^2),
		mean((mle3[,1] - 0.5)^2),
		mean((remle3[,1] - 0.5)^2),
		mean(mle3[,3]),
		mean(remle3[,3]),
		mean(mle3[,4]),
		mean(remle3[,4])
	)
)

# Create Tables 8.1 and 8.2
simDF1 = simDF[,1:8]
print(
    xtable(simDF1, 
      align = c('l',rep('l', times = length(simDF1[1,]))),
      digits = c(0,rep(3, times = 8)),
      caption = 'Covariates used for model-fitting',
      label = 'tab:covList'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

simDF2 = simDF[,9:12]
print(
    xtable(simDF2, 
      align = c('l',rep('l', times = length(simDF2[1,]))),
      digits = c(0,rep(3, times = 4)),
      caption = 'Covariates used for model-fitting',
      label = 'tab:covList'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

