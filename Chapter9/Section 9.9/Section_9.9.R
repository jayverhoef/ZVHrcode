sec_path = 'Rcode/Chapter9/Section 9.9/'
setwd(paste0(SLEDbook_path,sec_path))

library(spmodel)
library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Kriging in Practice: ML vs REML Predictions
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# Code for comparing quality of predictions based on ML and REML estimation (Table 9.6)
# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 8
locxy <- data.frame(x = c(matrix(rep(1:K,K),ncol=1,byrow=T),4.5, 4.5),
	y = c(matrix(rep(1:K,each=K),ncol=1,byrow=T), 4.5, 1))
n <- K^2

#theta^d, where d is distance, is equal to an exponential model with more
# typical form e^(-d/alpha), where alpha = -1/log(theta).  So, 
# theta = 0.3 implies alpha = 0.8305835
spcov_params_val <- spcov_params("exponential", de = 1, ie = 0, 
	range = -1/log(0.3))

#number of simulations
nsim = 1000
# store ml performance at theta = 0.3
ml03_width = matrix(NA, nrow = nsim, ncol = 2)
ml03_cover = matrix(NA, nrow = nsim, ncol = 2)
reml03_width = matrix(NA, nrow = nsim, ncol = 2)
reml03_cover = matrix(NA, nrow = nsim, ncol = 2)
for(kk in 1:nsim) {
	# simulate errors
	errsim = sprnorm(spcov_params_val, data = locxy, xcoord = x, ycoord = y)
	# constant but unknown mean = 0, simulated values at observed and prediction
	# locations
	z = errsim[1:64]
	zp = errsim[65:66]

	# create a data.frame for observed and prediction locations
	obsDF = data.frame(x = locxy$x[1:64], y = locxy$y[1:64], z = z)
	predDF = data.frame(x = locxy$x[65:66], y = locxy$y[65:66])

	# ml prediction at both locations
	ml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial("exponential", ie = 1e-8, known = c("ie")),
		estmethod = 'ml')
	ml_pred = predict(ml_out, predDF, se.fit = TRUE)
	# reml prediction at both locations
	reml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial("exponential", ie = 1e-8, known = c("ie")),
		estmethod = 'reml')
	reml_pred = predict(reml_out, predDF, se.fit = TRUE)
	# width of 90% interval using t-distribution
	ml03_width[kk,] = 2*ml_pred$se.fit*qt(.90, n - p)
	# coverage
	ml03_cover[kk,] = ml_pred$fit - ml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < ml_pred$fit + ml_pred$se.fit*qt(.90,n - p)
	# width of 90% interval using t-distribution
	reml03_width[kk,] = 2*reml_pred$se.fit*qt(.90, n - p)
	# coverage
	reml03_cover[kk,] = reml_pred$fit - reml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < reml_pred$fit + reml_pred$se.fit*qt(.90,n - p)
}
apply(ml03_width,2,mean)
apply(reml03_width,2,mean)
apply(ml03_cover,2,mean)
apply(reml03_cover,2,mean)



#theta^d, where d is distance, is equal to an exponential model with more
# typical form e^(-d/alpha), where alpha = -1/log(theta).  So, 
# theta = 0.7 implies alpha = 2.803673
spcov_params_val <- spcov_params("exponential", de = 1, ie = 0, 
	range = -1/log(0.7))

#number of simulations
nsim = 1000
# store ml performance at theta = 0.3
ml07_width = matrix(NA, nrow = nsim, ncol = 2)
ml07_cover = matrix(NA, nrow = nsim, ncol = 2)
reml07_width = matrix(NA, nrow = nsim, ncol = 2)
reml07_cover = matrix(NA, nrow = nsim, ncol = 2)
for(kk in 1:nsim) {
	# simulate errors
	errsim = sprnorm(spcov_params_val, data = locxy, xcoord = x, ycoord = y)
	# constant but unknown mean = 0, simulated values at observed and prediction
	# locations
	z = errsim[1:64]
	zp = errsim[65:66]

	# create a data.frame for observed and prediction locations
	obsDF = data.frame(x = locxy$x[1:64], y = locxy$y[1:64], z = z)
	predDF = data.frame(x = locxy$x[65:66], y = locxy$y[65:66])

	# ml prediction at both locations
	ml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial("exponential", ie = 1e-8, known = c("ie")),
		estmethod = 'ml')
	ml_pred = predict(ml_out, predDF, se.fit = TRUE)
	# reml prediction at both locations
	reml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_initial = spcov_initial("exponential", ie = 1e-8, known = c("ie")),
		estmethod = 'reml')
	reml_pred = predict(reml_out, predDF, se.fit = TRUE)
	# width of 90% interval using t-distribution
	ml07_width[kk,] = 2*ml_pred$se.fit*qt(.90, n - p)
	# coverage
	ml07_cover[kk,] = ml_pred$fit - ml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < ml_pred$fit + ml_pred$se.fit*qt(.90,n - p)
	# width of 90% interval using t-distribution
	reml07_width[kk,] = 2*reml_pred$se.fit*qt(.90, n - p)
	# coverage
	reml07_cover[kk,] = reml_pred$fit - reml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < reml_pred$fit + reml_pred$se.fit*qt(.90,n - p)
}
apply(ml07_width,2,mean)
apply(reml07_width,2,mean)
apply(ml07_cover,2,mean)
apply(reml07_cover,2,mean)


# create the table
ml_vs_repl_pred = rbind(
	c(apply(ml03_width,2,mean), apply(ml07_width,2,mean)),
  c(apply(reml03_width,2,mean), apply(reml07_width,2,mean)),
	c(apply(ml03_cover,2,mean), apply(ml07_cover,2,mean)),
	c(apply(reml03_cover,2,mean), apply(reml07_cover,2,mean))
)

print(
    xtable(ml_vs_repl_pred, 
      align = c('l',rep('l', times = length(ml_vs_repl_pred[1,]))),
      digits = c(0, rep(2, times = 4))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)














#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Kriging in Practice: ML vs REML Predictions
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# Code for comparing quality of predictions based on ML and REML estimation (Table 9.6)
# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 13
locxy <- data.frame(x = c(matrix(rep(1:K,K),ncol=1,byrow=T), 7, 7),
	y = c(matrix(rep(1:K,each=K),ncol=1,byrow=T), 7, 1))
n <- K^2

# theta^d, where d is distance, is equal to an exponential model with more
# typical form: e^(-d/alpha), where alpha = -1/log(theta).  So, 
# theta = 0.3 implies alpha = 0.8305835
spcov_params_val <- spcov_params("exponential", de = 0.9, ie = 0.1, 
	range = -1/log(0.3))

#number of simulations
nsim = 500
# store ml performance at theta = 0.3
ml03_width = matrix(NA, nrow = nsim, ncol = 2)
ml03_cover = matrix(NA, nrow = nsim, ncol = 2)
reml03_width = matrix(NA, nrow = nsim, ncol = 2)
reml03_cover = matrix(NA, nrow = nsim, ncol = 2)
for(kk in 1:nsim) {
	# simulate errors
	errsim = sprnorm(spcov_params_val, data = locxy, xcoord = x, ycoord = y)
	# constant but unknown mean = 0, simulated values at observed and prediction
	# locations
	z = errsim[1:n]
	zp = errsim[(n + 1):(n + 2)]

	# create a data.frame for observed and prediction locations
	obsDF = data.frame(x = locxy$x[1:n], y = locxy$y[1:n], z = z)
	predDF = data.frame(x = locxy$x[(n + 1):(n + 2)], 
		y = locxy$y[(n + 1):(n + 2)])

	# ml prediction at both locations
	ml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_type = "exponential", 
		estmethod = 'ml')
	ml_pred = predict(ml_out, predDF, se.fit = TRUE)
	# reml prediction at both locations
	reml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_type = "exponential", 
		estmethod = 'reml')
	reml_pred = predict(reml_out, predDF, se.fit = TRUE)
	# width of 90% interval using t-distribution
	ml03_width[kk,] = 2*ml_pred$se.fit*qt(.90, n - p)
	# coverage
	ml03_cover[kk,] = ml_pred$fit - ml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < ml_pred$fit + ml_pred$se.fit*qt(.90, n - p)
	# width of 90% interval using t-distribution
	reml03_width[kk,] = 2*reml_pred$se.fit*qt(.90, n - p)
	# coverage
	reml03_cover[kk,] = reml_pred$fit - reml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < reml_pred$fit + reml_pred$se.fit*qt(.90,n - p)
}
apply(ml03_width,2,mean)
apply(reml03_width,2,mean)
apply(ml03_cover,2,mean)
apply(reml03_cover,2,mean)



#theta^d, where d is distance, is equal to an exponential model with more
# typical form e^(-d/alpha), where alpha = -1/log(theta).  So, 
# theta = 0.7 implies alpha = 2.803673
spcov_params_val <- spcov_params("exponential", de = 0.9, ie = 0.1, 
	range = -1/log(0.7))

#number of simulations
nsim = 500
# store ml performance at theta = 0.3
ml07_width = matrix(NA, nrow = nsim, ncol = 2)
ml07_cover = matrix(NA, nrow = nsim, ncol = 2)
reml07_width = matrix(NA, nrow = nsim, ncol = 2)
reml07_cover = matrix(NA, nrow = nsim, ncol = 2)
for(kk in 1:nsim) {
	# simulate errors
	errsim = sprnorm(spcov_params_val, data = locxy, xcoord = x, ycoord = y)
	# constant but unknown mean = 0, simulated values at observed and prediction
	# locations
	z = errsim[1:n]
	zp = errsim[(n + 1):(n + 2)]

	# create a data.frame for observed and prediction locations
	obsDF = data.frame(x = locxy$x[1:n], y = locxy$y[1:n], z = z)
	predDF = data.frame(x = locxy$x[(n + 1):(n + 2)], 
		y = locxy$y[(n + 1):(n + 2)])

	# ml prediction at both locations
	ml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_type = "exponential", 
		estmethod = 'ml')
	ml_pred = predict(ml_out, predDF, se.fit = TRUE)
	# reml prediction at both locations
	reml_out = splm(z ~ 1, data = obsDF, xcoord = 'x', ycoord = 'y',
		spcov_type = "exponential", 
		estmethod = 'reml')
	reml_pred = predict(reml_out, predDF, se.fit = TRUE)
	# width of 90% interval using t-distribution
	ml07_width[kk,] = 2*ml_pred$se.fit*qt(.90, n - p)
	# coverage
	ml07_cover[kk,] = ml_pred$fit - ml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < ml_pred$fit + ml_pred$se.fit*qt(.90,n - p)
	# width of 90% interval using t-distribution
	reml07_width[kk,] = 2*reml_pred$se.fit*qt(.90, n - p)
	# coverage
	reml07_cover[kk,] = reml_pred$fit - reml_pred$se.fit*qt(.90,n - p) < zp & 
		zp < reml_pred$fit + reml_pred$se.fit*qt(.90,n - p)
}
apply(ml07_width,2,mean)
apply(reml07_width,2,mean)
apply(ml07_cover,2,mean)
apply(reml07_cover,2,mean)


# create the table
ml_vs_repl_pred = rbind(
	c(apply(ml03_width,2,mean), apply(ml07_width,2,mean)),
  c(apply(reml03_width,2,mean), apply(reml07_width,2,mean)),
	c(apply(ml03_cover,2,mean), apply(ml07_cover,2,mean)),
	c(apply(reml03_cover,2,mean), apply(reml07_cover,2,mean))
)

print(
    xtable(ml_vs_repl_pred, 
      align = c('l',rep('l', times = length(ml_vs_repl_pred[1,]))),
      digits = c(0, rep(2, times = 4))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

