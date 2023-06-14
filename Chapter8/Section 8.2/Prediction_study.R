#set a path as a working directory
sec_path = 'Rcode/Chapter8/Section 8.2/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)
library(spmodel)

# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 12
locxy <- cbind(rep(1:K,K), rep(1:K, each = K))
n <- K^2
theta <- 0.99

# Define a covariance structure on the grid: an exponential covariance function with correlation theta at
# unit distance
d = as.matrix(dist(cbind(locxy)))
Gtrue = theta^d
tcholGtrue = t(chol(Gtrue))

# ----------------------------- Intercept Only ---------------------------------

nsim = 1000
mle1 <- matrix(0, nsim, 4)
remle1 <- matrix(0, nsim, 4)
set.seed(508)
for(k in 1:nsim){
	cat("\r", "iteration: ", k)
	# simulate data, all true regression coefficients are zero
	y1 <- 0 + tcholGtrue %*% rnorm(n)
	# create data.frame for spmodel
	DF = data.frame(y = y1, xcoord = locxy[,1], ycoord = locxy[,2])
	indNA = sort(sample(1:n, 10))
	yp = DF[indNA,'y']
	DF[indNA,'y'] = NA
	# ------- ML ESTIMATION
	# fit model using spmodel, holding nugget effect at zero
	fitout = splm(y ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'ml')
	summary(fitout)
	predout = predict(fitout, interval = 'prediction', level = 0.90)
	#bias
	mle1[k,1] = mean(predout[,'fit'] - yp)
	# rmspe
	mle1[k,2] = sqrt(mean((predout[,'fit'] - yp)^2))
	# prediction interval length
	mle1[k,3] = mean(predout[,'upr'] - predout[,'lwr'])
	# coverage
	mle1[k,4]	= mean(predout[,'lwr'] < yp & yp < predout[,'upr'])
	# ------- REML ESTIMATION
	fitout = splm(y ~ 1, data = DF, xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'reml')
	predout = predict(fitout, interval = 'prediction', level = 0.90)
	#bias
	remle1[k,1] = mean(predout[,'fit'] - yp)
	# rmspe
	remle1[k,2] = sqrt(mean((predout[,'fit'] - yp)^2))
	# prediction interval length
	remle1[k,3] = mean(predout[,'upr'] - predout[,'lwr'])
	# coverage
	remle1[k,4]	= mean(predout[,'lwr'] < yp & yp < predout[,'upr'])
}

round(c(mean(mle1[,1]),
mean(remle1[,1]),
mean(mle1[,2]),
mean(remle1[,2]),
mean(mle1[,3]),
mean(remle1[,3]),
mean(mle1[,4]),
mean(remle1[,4])
),4)

# ----------------------------- Row Effects ------------------------------------


nsim = 1000
mle2 <- matrix(0, nsim, 4)
remle2 <- matrix(0, nsim, 4)
set.seed(608)
for(k in 1:nsim){
	cat("\r", "iteration: ", k)
	# simulate data, all true regression coefficients are zero
	y1 <- 0 + tcholGtrue %*% rnorm(n)
	# create data.frame for spmodel
	DF = data.frame(y = y1, xcoord = locxy[,1], ycoord = locxy[,2])
	indNA = sort(sample(1:n, 10))
	yp = DF[indNA,'y']
	DF[indNA,'y'] = NA
	# ------- ML ESTIMATION
	# fit model using spmodel, holding nugget effect at zero
	fitout = splm(y ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'ml')
	summary(fitout)
	predout = predict(fitout, interval = 'prediction', level = 0.90)
	#bias
	mle2[k,1] = mean(predout[,'fit'] - yp)
	# rmspe
	mle2[k,2] = sqrt(mean((predout[,'fit'] - yp)^2))
	# prediction interval length
	mle2[k,3] = mean(predout[,'upr'] - predout[,'lwr'])
	# coverage
	mle2[k,4]	= mean(predout[,'lwr'] < yp & yp < predout[,'upr'])
	# ------- REML ESTIMATION
	fitout = splm(y ~ I(as.factor(xcoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'reml')
	predout = predict(fitout, interval = 'prediction', level = 0.90)
	#bias
	remle2[k,1] = mean(predout[,'fit'] - yp)
	# rmspe
	remle2[k,2] = sqrt(mean((predout[,'fit'] - yp)^2))
	# prediction interval length
	remle2[k,3] = mean(predout[,'upr'] - predout[,'lwr'])
	# coverage
	remle2[k,4]	= mean(predout[,'lwr'] < yp & yp < predout[,'upr'])
}

round(c(mean(mle2[,1]),
mean(remle2[,1]),
mean(mle2[,2]),
mean(remle2[,2]),
mean(mle2[,3]),
mean(remle2[,3]),
mean(mle2[,4]),
mean(remle2[,4])
),4)

# ---------------------- Row and Column Effects --------------------------------


nsim = 1000
mle3 <- matrix(0, nsim, 4)
remle3 <- matrix(0, nsim, 4)
set.seed(708)
for(k in 1:nsim){
	cat("\r", "iteration: ", k)
	# simulate data, all true regression coefficients are zero
	y1 <- 0 + tcholGtrue %*% rnorm(n)
	# create data.frame for spmodel
	DF = data.frame(y = y1, xcoord = locxy[,1], ycoord = locxy[,2])
	indNA = sort(sample(1:n, 10))
	yp = DF[indNA,'y']
	DF[indNA,'y'] = NA
	# ------- ML ESTIMATION
	# fit model using spmodel, holding nugget effect at zero
	fitout = splm(y ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'ml')
	summary(fitout)
	predout = predict(fitout, interval = 'prediction', level = 0.90)
	#bias
	mle3[k,1] = mean(predout[,'fit'] - yp)
	# rmspe
	mle3[k,2] = sqrt(mean((predout[,'fit'] - yp)^2))
	# prediction interval length
	mle3[k,3] = mean(predout[,'upr'] - predout[,'lwr'])
	# coverage
	mle3[k,4]	= mean(predout[,'lwr'] < yp & yp < predout[,'upr'])
	# ------- REML ESTIMATION
	fitout = splm(y ~ I(as.factor(xcoord)) + I(as.factor(ycoord)), data = DF, 
		xcoord = xcoord, ycoord = ycoord,
		spcov_initial = spcov_initial('exponential', ie = 0, known = c('ie')),
		estmethod = 'reml')
	predout = predict(fitout, interval = 'prediction', level = 0.90)
	#bias
	remle3[k,1] = mean(predout[,'fit'] - yp)
	# rmspe
	remle3[k,2] = sqrt(mean((predout[,'fit'] - yp)^2))
	# prediction interval length
	remle3[k,3] = mean(predout[,'upr'] - predout[,'lwr'])
	# coverage
	remle3[k,4]	= mean(predout[,'lwr'] < yp & yp < predout[,'upr'])
}

round(c(mean(mle3[,1]),
mean(remle3[,1]),
mean(mle3[,2]),
mean(remle3[,2]),
mean(mle3[,3]),
mean(remle3[,3]),
mean(mle3[,4]),
mean(remle3[,4])
),4)








simDF = 
rbind(
	c(mean(mle1[,1]),
		mean(remle1[,1]),
		mean(mle1[,2]),
		mean(remle1[,2]),
		mean(mle1[,3]),
		mean(remle1[,3]),
		mean(mle1[,4]),
		mean(remle1[,4])
	),
	c(mean(mle2[,1]),
		mean(remle2[,1]),
		mean(mle2[,2]),
		mean(remle2[,2]),
		mean(mle2[,3]),
		mean(remle2[,3]),
		mean(mle2[,4]),
		mean(remle2[,4])
	),
	c(mean(mle3[,1]),
		mean(remle3[,1]),
		mean(mle3[,2]),
		mean(remle3[,2]),
		mean(mle3[,3]),
		mean(remle3[,3]),
		mean(mle3[,4]),
		mean(remle3[,4])
	)
)

print(
    xtable(simDF, 
      align = c('l',rep('l', times = length(simDF[1,]))),
      digits = c(0,rep(4, times = 8)),
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

