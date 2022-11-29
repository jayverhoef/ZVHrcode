sec_path = 'Rcode/Chapter5/Section 5.6/'
setwd(paste0(SLEDbook_path,sec_path))

library(sp)
library(viridis)
library(classInt)

source('pointSimSyst.R')
source('pointSimCSR.R')
source('geostatSim.R')
source('distGeoAni.R')
source('corModels.R')
source('splmm.R')
source('covParmIni.R')
source('m2LLi.R')
source('m2LL.R')
source('makeCovMat.R')

#-------------------------------------------------------------------------------
# Simulations to see if we get bias or lack of proper coverage, 
# and check out RMSE for GLS vs. OLS
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
    xy = pointSimCSR(npts, lower.x.lim = 0, upper.x.lim = 1, 
        lower.y.lim = 0, upper.y.lim = 1) 
    #create distance matrix 
    dismat = distGeoAni(xy$x, xy$y, xy$x, xy$y, rotate = 90, range = 1.0, minorp = 1)
    #create covariance matrix. Different options can be found in corModels.R file
    covmat = corModelExponential(dismat) + diag(rep(0.0001, times = npts))
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
    sim_reg_DF = data.frame(y = y, X1 = X1, X2 = X2)

    # turn it into a SpatialPointsDataFrame class from sp package
    sp_sim_reg_DF = SpatialPointsDataFrame(xy,sim_reg_DF)

    # the function splmm fits spatial linear models to SpatialPointsDataFrame  
    # for geostatistical data using REML
    fit_geo = splmm(y ~ X1 + X2, sp_sim_reg_DF, varComps = "exponential")	 
     
  # store results of fits and their standard errors
    geobeta1est[kk] = fit_geo$betaHat[2]
    geobeta2est[kk] = fit_geo$betaHat[3]
    geobeta1se[kk] = sqrt(diag(fit_geo$covBetaHat))[2]
    geobeta2se[kk] = sqrt(diag(fit_geo$covBetaHat))[3]
   
  # Now, which assumes independence (so no confounding)
    fit_ind = summary(lm(y ~ X1 + X2, sim_reg_DF))
    
    indbeta1est[kk] = fit_ind$coefficients[2,1]
    indbeta2est[kk] = fit_ind$coefficients[3,1]
    indbeta1se[kk] = fit_ind$coefficients[2,2]
    indbeta2se[kk] = fit_ind$coefficients[3,2]

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
      sum(geobeta1est - 1.645*geobeta1se < 1 & 1 < geobeta1est + 1.645*geobeta1se)/nsim,
      sum(geobeta2est - 1.645*geobeta2se < 1 & 1 < geobeta2est + 1.645*geobeta2se)/nsim
    )
  )
  # OLS 90% confidence interval coverage
  CI90_OLS = rbind(CI90_OLS,
    c(
      sum(indbeta1est - 1.645*indbeta1se < 1 & 1 < indbeta1est + 1.645*indbeta1se)/nsim,
      sum(indbeta2est - 1.645*indbeta2se < 1 & 1 < indbeta2est + 1.645*indbeta2se)/nsim
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
  
file_name = "zinX_exact_RMSE"
pdf(file = paste0(file_name,'.pdf'), width = 10, height = 8)
  old.par = par(mar = c(5,6,1,1), mgp=c(3.7,1,0))
  plot(1:dim(VI_ev)[1],VI_ev$RMSE_OLS.1, pch = 1, ylim = c(0,max(VI_ev[,6:9])),
    ylab = 'RMSE', xlab = expression(italic(i)*th~Simulation~Set~with~Ordered~Eigenvector~bold("q")[italic(i)]^symbol("\052")==bold("x")[2]), cex = 2, cex.lab = 2, cex.axis = 1.5, xaxt = 'n')
  axis(1, at = 1:dim(VI_ev)[1], labels = c(evi_range))
  points(1:dim(VI_ev)[1],VI_ev$RMSE_GLS.1, pch = 19, cex = 1.5)
  points(1:dim(VI_ev)[1],VI_ev$RMSE_OLS.2, pch = 5, cex = 2)
  points(1:dim(VI_ev)[1],VI_ev$RMSE_GLS.2, pch = 18, cex = 2)
  legend(10, 2.8, legend = c(expression(OLS~estimation~of~beta~"for"~bold("x")[2]==bold("q")[italic(i)]^symbol("\052")),
    expression(GLS~estimation~of~beta~"for"~bold("x")[2]==bold("q")[italic(i)]^symbol("\052")),
    expression(OLS~estimation~of~beta~"for"~bold("x")[3]~(Indep.)),
    expression(GLS~estimation~of~beta~"for"~bold("x")[3]~(Indep.))),
    pch = c(1, 19, 5, 18), cex = 1.5)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

file_name = "zinX_exact_CI90.pdf"
pdf(file = paste0(file_name,'.pdf'), width = 10, height = 8)
  old.par = par(mar = c(5,6,1,1), mgp=c(3.7,1,0))
  plot(1:dim(VI_ev)[1],VI_ev$CI90_OLS.1, pch = 1, ylim = c(0,1),
    ylab = '90 % Confidence Interval Coverage', xlab = expression(italic(i)*th~Simulation~Set~with~Ordered~Eigenvector~bold("q")[italic(i)]^symbol("\052")==bold("x")[2]), cex = 2,
    cex.lab = 2, cex.axis = 1.5, xaxt = 'n')
  axis(1, at = 1:dim(VI_ev)[1], labels = c(evi_range))
  lines(c(1,20), c(.9,.9), lty = 2, lwd = 3)
  points(1:dim(VI_ev)[1],VI_ev$CI90_GLS.1, pch = 19, cex = 1.5)
  points(1:dim(VI_ev)[1],VI_ev$CI90_OLS.2, pch = 5, cex = 2)
  points(1:dim(VI_ev)[1],VI_ev$CI90_GLS.2, pch = 18, cex = 2)
  legend(10, .5, legend = c(expression(OLS~estimation~of~beta~"for"~bold("x")[2]==bold("q")[italic(i)]^symbol("\052")),
    expression(GLS~estimation~of~beta~"for"~bold("x")[2]==bold("q")[italic(i)]^symbol("\052")),
    expression(OLS~estimation~of~beta~"for"~bold("x")[3]~(Indep.)),
    expression(GLS~estimation~of~beta~"for"~bold("x")[3]~(Indep.))),
    pch = c(1, 19, 5, 18), cex = 1.5)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

