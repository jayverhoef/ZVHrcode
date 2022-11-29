sec_path = 'Rcode/Chapter5/Section 5.5/'
setwd(paste0(SLEDbook_path,sec_path))

library(sp)
library(viridis)
library(classInt)
library(spmodel)

source('pointSimSyst.R')
source('corModels.R')
source('distGeoAni.R')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#     Investigate spatial patterning of eigenvectors for systematic points 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Figure of Eigenvalues

#create systematic grid of points
xy = data.frame(
	x = (kronecker(1:40, rep(1, times = 40)) - 0.5)/40, 
	y = (kronecker(rep(1, times = 40), 1:40) - 0.5)/40 )
	
#create distance matrix 
dismat = as.matrix(dist(xy))
#create covariance matrix
covmat = exp(-dismat/0.5) + diag(rep(0.000001, times = 40*40))
#eigen decomposition
eig_covmat = eigen(covmat)

#make a plot of eigenvalue, and map eigenvectors spatially
# make coordinates for image plot
x = ((1:40) - .5)/40
y = x
# pick a set of 25 eigenvectors to plot
evecset = c(1, 2, 3, 4, 5,
	6, 8, 10, 12, 14, 
	20, 30, 40, 50, 60, 
	100, 200, 300, 400, 500,
	1596, 1597, 1598, 1599, 1600 )
# number of color classes
nbrks = 20

file_name = "figures/Eigval_bg"

pdf(file = paste0(file_name,'.pdf'), width = 10, height = 10)
  oldpar = par(mar = c(5,5,1,1))
  plot(sqrt(eig_covmat$values), pch = 19, cex = 1, cex.axis = 1.5, cex.lab = 2,
  xlab = 'index', ylab = 'Square Root of Eigenvalue')
  par(oldpar)
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

# Figure of Eigenvectors

file_name = "figures/Eigvec25_bg"
nbrks = 20

pdf(file = paste0(file_name,'.pdf'), width = 8, height = 8)
  oldpar = par(mar = c(2,2,3,1))
  layout(matrix(1:25, nrow = 5, byrow = TRUE))
  for(i in 1:length(evecset)) {
    z = eig_covmat$vector[,evecset[i]]*sqrt(eig_covmat$values[evecset[i]])
    z = matrix(z, nrow = 40)
    brks = quantile(z, probs = (0:nbrks)/nbrks)
    cramp = viridis(nbrks)
    image(x, y, z, breaks = brks, col = cramp, 
      main = paste0("Eigenvector ",evecset[i]),
      cex.main = 1.2, cex.axis = .8, xlab = '', ylab = '')
  #	points(xy[,1],xy[,2], pch = 19, cex = 2)
  }
  oldpar = par(mar = c(2,2,3,1))
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#      Simulations to see if we get bias or lack of proper coverage, 
#                  and check out RMSE for GLS vs. OLS
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# pick number of points to simulate
npts = 120
# Number of simulations
nsim = 1000
# set of eigenvalues
evi_range = c(1:10, 20, 30, 40, 50, 60, 70, 100, 115, 119, 120)
# choose regression coefficients
beta0 = 0  # overall mean
beta1 = 1  # coefficient for chosen eigenvector
beta2 = 1  # coefficient for X2, generated from IID normal(0,1)

#set the random number seed so results are reproducible
set.seed(121)

# empty vectors to store results
bias_OLS = NULL
bias_GLS = NULL
RMSE_OLS = NULL
RMSE_GLS = NULL
CI90_OLS = NULL
CI90_GLS = NULL

startTime = Sys.time()
for(evi in evi_range) {

  # pick out an eigenvector to copy over to design matrix
  eveci1 = evi

  # empty vectors to store results of fits and their standard errors
  geobeta1est = rep(NA, times = nsim)
  geobeta1se = rep(NA, times = nsim)
  geobeta2est = rep(NA, times = nsim)
  geobeta2se = rep(NA, times = nsim)
  indbeta1est = rep(NA, times = nsim)
  indbeta1se = rep(NA, times = nsim)
  indbeta2est = rep(NA, times = nsim)
  indbeta2se = rep(NA, times = nsim)

  cat("\n")
  for(kk in 1:nsim) {
    
    cat("\r", "Eigenvector: ", evi, "      Iter Number: ", kk)
    
    # simulate points as spatially random
    xy = data.frame(x = runif(npts), y = runif(npts))
    #create distance matrix 
    dismat = as.matrix(dist(xy))
    #create covariance matrix. Different options can be found in corModels.R file
    covmat = exp(-dismat) + diag(rep(0.0001, times = npts))
    #eigen decomposition
    eig_covmat = eigen(covmat)

    # scale the eigenvectors by the square root of the corresponding eigenvalue
    X1 = eig_covmat$vector[,eveci1]*sqrt(eig_covmat$values[eveci1])
    # add two completely independent covariates
    X2 = rnorm(npts, 0, 1)

    # get Cholesky of covariance matrix
    L = chol(covmat)
    # simulate spatially autocorrelated errors with prescribed covariance matrix
    errsim = t(L) %*% rnorm(npts, 0, 1)

    # create the response variable from the linear model
    # the vectors X1 and X2 are in fixed effects, and eigenvectors of covariance
    y = beta0 + beta1*X1 + beta2*X2 + errsim

    # create a data.frame of the data
    sim_reg_DF = data.frame(y = y, X1 = X1, X2 = X2, 
			xcoord = xy$x, ycoord = xy$y)

    # turn it into a SpatialPointsDataFrame class from sp package
    sp_sim_reg_DF = SpatialPointsDataFrame(xy,sim_reg_DF)

    # the function splmm fits spatial linear models to SpatialPointsDataFrame  
    # for geostatistical data using REML
#    fit_geo = splmm(y ~ X1 + X2, sp_sim_reg_DF, varComps = "exponential")	 
    fit_spmodel = splm(y ~ X1 + X2, data = sim_reg_DF, 
			xcoord = xcoord, ycoord = ycoord, spcov_type = 'exponential')
  # store results of fits and their standard errors
    geobeta1est[kk] = coef(fit_spmodel)['X1']
    geobeta2est[kk] = coef(fit_spmodel)['X2']
    geobeta1se[kk] = sqrt(diag(vcov(fit_spmodel)))['X1']
    geobeta2se[kk] = sqrt(diag(vcov(fit_spmodel)))['X2']
		sqrt(diag(vcov(fit_spmodel)))['X1']
   
  # Now, which assumes independence (so no confounding)
    fit_ind = summary(lm(y ~ X1 + X2, sim_reg_DF))
    
    indbeta1est[kk] = coef(fit_ind)['X1','Estimate']
    indbeta2est[kk] = coef(fit_ind)['X2','Estimate']
    indbeta1se[kk] = coef(fit_ind)['X1','Std. Error']
    indbeta2se[kk] = coef(fit_ind)['X2','Std. Error']

  }

  # GLS bias, where first two are confounded covariates, second two are not
  bias_GLS = rbind(bias_GLS,
    c(mean(geobeta1est) - 1, mean(geobeta2est) - 1)
  )
  # OLS bias, where first two are confounded covariates, second two are not
  bias_OLS = rbind(bias_OLS,
    c(mean(indbeta1est) - 1, mean(indbeta2est) - 1)
  )

  # GLS root mean-squared error
  RMSE_GLS = rbind(RMSE_GLS,
    sqrt(c(mean((geobeta1est - 1)^2), mean((geobeta2est - 1)^2)))
  )
  # OLS root mean-squared error
  RMSE_OLS = rbind(RMSE_OLS,
    sqrt(c(mean((indbeta1est - 1)^2), mean((indbeta2est - 1)^2)))
  )

  # GLS 90% confidence interval coverage
  CI90_GLS = rbind(CI90_GLS,
    c(
      sum(geobeta1est - 1.645*geobeta1se < 1 & 
				1 < geobeta1est + 1.645*geobeta1se)/nsim,
      sum(geobeta2est - 1.645*geobeta2se < 1 & 
				1 < geobeta2est + 1.645*geobeta2se)/nsim
    )
  )
  # OLS 90% confidence interval coverage
  CI90_OLS = rbind(CI90_OLS,
    c(
      sum(indbeta1est - 1.645*indbeta1se < 1 & 
				1 < indbeta1est + 1.645*indbeta1se)/nsim,
      sum(indbeta2est - 1.645*indbeta2se < 1 & 
				1 < indbeta2est + 1.645*indbeta2se)/nsim
    )
  )

}
  
endTime = Sys.time()
difftime(endTime, startTime, units = 'mins')

VI_ev = data.frame(EV = evi_range, bias_OLS.1 = bias_OLS[,1], bias_OLS.2 = bias_OLS[,2],
  bias_GLS.1 = bias_GLS[,1], bias_GLS.2 = bias_GLS[,2],
  RMSE_OLS.1 = RMSE_OLS[,1], RMSE_OLS.2 = RMSE_OLS[,2],
  RMSE_GLS.1 = RMSE_GLS[,1], RMSE_GLS.2 = RMSE_GLS[,2],
  CI90_OLS.1 = CI90_OLS[,1], CI90_OLS.2 = CI90_OLS[,2],
  CI90_GLS.1 = CI90_GLS[,1], CI90_GLS.2 = CI90_GLS[,2])
 
#save the VI_ev object so you don't have to do simulations again.
save(VI_ev, file = 'VI_ev.Rdata')
# If you want to start here (for making graphics, etc.)
#load('VI_ev.Rdata')   
  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#       RMSE Figure
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/zinX_exact_RMSE"
pdf(file = paste0(file_name,'.pdf'), width = 10, height = 8)
  old.par = par(mar = c(5,6,1,1), mgp=c(3.7,1,0))
  plot(1:dim(VI_ev)[1],VI_ev$RMSE_OLS.1, pch = 1, ylim = c(0,max(VI_ev[,6:9])),
    ylab = 'RMSE', xlab = expression(italic(i)*th~Simulation~Set~with~Ordered~Eigenvector~bold("q")[italic(i)]^symbol("\052")==bold("x")[2]), cex = 3.5, 
    lwd = 2, cex.lab = 2, cex.axis = 1.5, xaxt = 'n')
  axis(1, at = 1:dim(VI_ev)[1], labels = c(evi_range))
  points(1:dim(VI_ev)[1],VI_ev$RMSE_GLS.1, pch = 19, cex = 2.5)
  points(1:dim(VI_ev)[1],VI_ev$RMSE_OLS.2, pch = 5, cex = 3)
  points(1:dim(VI_ev)[1],VI_ev$RMSE_GLS.2, pch = 18, cex = 3)
  legend(6, 2.8, legend = c(expression(OLS~estimation~of~beta~"for"~bold("x")[2]==bold("q")[italic(i)]^symbol("\052")),
    expression(GLS~estimation~of~beta~"for"~bold("x")[2]==bold("q")[italic(i)]^symbol("\052")),
    expression(OLS~estimation~of~beta~"for"~bold("x")[3]~(Indep.)),
    expression(GLS~estimation~of~beta~"for"~bold("x")[3]~(Indep.))),
    pch = c(1, 19, 5, 18), cex = 2.05)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#       90% Confidence Interval Coverage Figure
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/zinX_exact_CI90"
pdf(file = paste0(file_name,'.pdf'), width = 10, height = 8)
  old.par = par(mar = c(5,6,1,1), mgp=c(3.7,1,0))
  plot(1:dim(VI_ev)[1],VI_ev$CI90_OLS.1, pch = 1, ylim = c(0,1),
    ylab = '90 % Confidence Interval Coverage', xlab = expression(italic(i)*th~Simulation~Set~with~Ordered~Eigenvector~bold("q")[italic(i)]^symbol("\052")==bold("x")[2]), cex = 3, lwd = 2,
    cex.lab = 2, cex.axis = 1.5, xaxt = 'n')
  axis(1, at = 1:dim(VI_ev)[1], labels = c(evi_range))
  lines(c(1,20), c(.9,.9), lty = 2, lwd = 3)
  points(1:dim(VI_ev)[1],VI_ev$CI90_GLS.1, pch = 19, cex = 2.5)
  points(1:dim(VI_ev)[1],VI_ev$CI90_OLS.2, pch = 5, cex = 3)
  points(1:dim(VI_ev)[1],VI_ev$CI90_GLS.2, pch = 18, cex = 3)
  legend(6, .34, legend = c(expression(OLS~estimation~of~beta~"for"~bold("x")[2]==bold("q")[italic(i)]^symbol("\052")),
    expression(GLS~estimation~of~beta~"for"~bold("x")[2]==bold("q")[italic(i)]^symbol("\052")),
    expression(OLS~estimation~of~beta~"for"~bold("x")[3]~(Indep.)),
    expression(GLS~estimation~of~beta~"for"~bold("x")[3]~(Indep.))),
    pch = c(1, 19, 5, 18), cex = 2.05)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#       Type I error rates associated with the red shift
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Code for obtaining the Type I error rates associated with the red shift, as presented in Table 5.3.
# Set up a small example of 25 sites on a 5x5 square grid with unit spacing:
locx <- matrix(rep(1:5, 5), ncol = 1, byrow = T)
locy <- matrix(rep(1:5, each = 5), ncol = 1, byrow = T)

# Define a covariance structure on the grid: an exponential function with correlation 0.5 at one unit of distance
G1 <- matrix(0, 25, 25)
G2 <- G1
for(i in 1:25){
 for(j in 1:25){
  G1[i,j] <- 0.5^(sqrt((locx[i, 1] - locx[j, 1])^2 + 
		(locy[i, 1] - locy[j, 1])^2))
  G2[i,j] <- 0.25^(sqrt((locx[i, 1] - locx[j, 1])^2 + 
		(locy[i, 1] - locy[j, 1])^2))
}}

G1i = solve(G1)
G2i = solve(G2)
cholG1 = chol(G1)
cholG2 = chol(G2)

# Create intercept
int1 <- matrix(1, 25, 1)

# Simulations: simulate response and regressors, compute ols and gls estimators of intercept and slopes, and count the number of times (out of 10000) that the nominal 95% t-based confidence intervals for the slopes contains the true value (0.0).  Default setting is for y to be autocorrelated and for some of the regressors to be autocorrelated; relevant lines of codes should be replaced by lines with percent signs for the cases of non-autocorrelated y and non-autocorrelated regressors.
store_ols_no_no = rep(0, times = 3)
store_ols_yes_no = rep(0, times = 3)
store_ols_no_yes = rep(0, times = 3)
store_ols_yes_yes = rep(0, times = 3)
store_gls_yes_no = rep(0, times = 3)
store_gls_yes_yes = rep(0, times = 3)
set.seed(1069)
niter = 100000
cat("\r")
for(i in 1:niter){
	cat("\r", "Iter Number: ", i)
	y_no <- int1 + rnorm(25)
	y_yes <- int1+t(cholG1)%*%rnorm(25)
	x1_no <- rnorm(25)
	x2_no <- rnorm(25)
	x3 <- rnorm(25)
	x1_yes <- t(cholG1)%*%rnorm(25)
	x2_yes <- t(cholG2)%*%rnorm(25)
	Xreg <- cbind(x1_no, x2_no, x3, x1_yes, x2_yes)
	# rescaled by post-multiplying by the inverse square root 
	# of the sample covariance matrix,
	covXreg <- cov(Xreg)
	covXreghalf <- eigen(covXreg)$vectors %*% 
		diag(sqrt(eigen(covXreg)$values)) %*% t(eigen(covXreg)$vectors)
	XregcovXreghalfi <- Xreg %*% solve(covXreghalf)
	# create design matrices
	X_no <- cbind(int1, XregcovXreghalfi[,1:3])
	X_yes <- cbind(int1, XregcovXreghalfi[,c(4,5,3)])
	# fixed effects correlation matrices
	covols_no <- solve(t(X_no) %*% X_no)
	covols_yes <- solve(t(X_yes) %*% X_yes)
	covgls_no <- solve(t(X_no) %*% G1i %*% X_no)
	covgls_yes <- solve(t(X_no) %*% G1i %*% X_no)
	# fixed effects estimates
	ols_no_no <- covols_no %*% t(X_no) %*% y_no
	ols_yes_no <- covols_no %*% t(X_no) %*% y_yes
	ols_no_yes <- covols_yes %*% t(X_yes) %*% y_no
	ols_yes_yes <- covols_yes %*% t(X_yes) %*% y_yes
	gls_yes_no <- covgls_no %*% t(X_no) %*% G1i %*% y_yes
	gls_yes_yes <- covgls_yes %*%	t(X_yes) %*% G1i %*% y_yes
	# estimated standard errors from sigma^2 times the correlation matrix
	olsse_no_no <- sqrt(diag(as.numeric(t(y_no - X_no %*% ols_no_no) %*% 
		(y_no - X_no %*% ols_no_no)/21) * covols_no))
	olsse_yes_no <- sqrt(diag(as.numeric(t(y_yes - X_no %*% ols_yes_no) %*% 
		(y_yes - X_no %*% ols_yes_no)/21) * covols_no))
	olsse_no_yes <- sqrt(diag(as.numeric(t(y_no - X_yes %*% ols_no_yes) %*% 
		(y_no - X_yes %*% ols_no_yes)/21) * covols_yes))
	olsse_yes_yes <- sqrt(diag(as.numeric(t(y_yes - X_yes %*% ols_yes_yes) %*% 
		(y_yes - X_yes %*% ols_yes_yes)/21) * covols_yes))
	glsse_yes_no <- sqrt(diag(as.numeric(t(y_yes - X_no %*% gls_yes_no) %*% 
		solve(G1) %*% (y_yes - X_no %*% gls_yes_no)/21) * covgls_no))
	glsse_yes_yes <- sqrt(diag(as.numeric(t(y_yes - X_yes %*% gls_yes_yes) %*% 
		solve(G1) %*% (y_yes - X_yes %*% gls_yes_yes)/21) * covgls_yes))
	# indicator of confidence interval coverage at 95% using t with 21 df.
	store_ols_no_no = store_ols_no_no + 
		((ols_no_no - 2.079614*olsse_no_no < 0 & 
		ols_no_no + 2.079614*olsse_no_no > 0)[2:4])*1
	store_ols_yes_no = store_ols_yes_no + 
		((ols_yes_no - 2.079614*olsse_yes_no < 0 & 
		ols_yes_no + 2.079614*olsse_yes_no > 0)[2:4])*1
	store_ols_no_yes = store_ols_no_yes + 
		((ols_no_yes - 2.079614*olsse_no_yes < 0 & 
		ols_no_yes + 2.079614*olsse_no_yes > 0)[2:4])*1
	store_ols_yes_yes = store_ols_yes_yes + 
		((ols_yes_yes - 2.079614*olsse_yes_yes < 0 & 
		ols_yes_yes + 2.079614*olsse_yes_yes > 0)[2:4])*1
	store_gls_yes_no = store_gls_yes_no + 
		((gls_yes_no - 2.079614*glsse_yes_no < 0 & 
		gls_yes_no + 2.079614*glsse_yes_no > 0)[2:4])*1
	store_gls_yes_yes = store_gls_yes_yes + 
		((gls_yes_yes - 2.079614*glsse_yes_yes < 0 & 
		gls_yes_yes + 2.079614*glsse_yes_yes > 0)[2:4])*1
}

store_results = rbind(
	1 - store_ols_no_no/niter,
	1 - store_ols_yes_no/niter,
	1 - store_gls_yes_no/niter,
	1 - store_ols_no_yes/niter,
	1 - store_ols_yes_yes/niter,
	1 - store_gls_yes_yes/niter)

store_df = data.frame(y_autoc = c('no', 'yes', 'yes', 'no', 'yes', 'yes'),
	X_autoc = c('no', 'no', 'no', 'yes', 'yes', 'yes'),
	OLS_GLS = c('OLS', 'OLS', 'GLS', 'OLS', 'OLS', 'GLS'),
	x1 = store_results[,1],
	x2 = store_results[,2],
	x3 = store_results[,3]
)

print(
    xtable(store_df, 
      align = c('l',rep('l', times = length(store_df[1,]))),
      digits = c(0, rep(0, times = 3), rep(4, times = 3))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
