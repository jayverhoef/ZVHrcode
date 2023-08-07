sec_path = 'Rcode/Chapter9/Section 9.12/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Get the Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# attach data library
library(ZVHdata)
library(sf)
library(viridis)
library(classInt)
library(colorspace)
library(gstat)
library(spmodel)
library(stringr)
library(xtable)
source('addBreakColorLegend.R')


# load data for graphics and analysis
data(MOSSobs)
data(MOSSpreds)
data(CAKRboundary)


# transform some of the variables
DF = MOSSobs
DF$easting = st_coordinates(DF)[,1]/1e+3
DF$northing = st_coordinates(DF)[,2]/1e+3
DF = st_drop_geometry(DF)
DF$year = as.factor(DF$year)
DF$field_dup = as.factor(DF$field_dup)
DF$dist2road = log(DF$dist2road)
DF$Pb = log(DF$Pb)



#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Conditional Simulation
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# fit the model
PbOut = splm(Pb ~ year + dist2road + sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing', spcov_type = 'exponential',
	partition_factor = ~ year,
	random = ~ (1 | sample) + (1 |sample:field_dup),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(PbOut)

# get the estimated covariance parameters for spatial model
theta_sp = coef(PbOut, type = 'spcov')[1:3]
# get the estimated variances of the random effectds
theta_re = coef(PbOut, type = 'randcov')
# total sample size
n = dim(DF)[1]

# create a mask with zeros and ones to force covariance matrix
# between years to be zero
part_mask = outer(DF$year == '2001', DF$year == '2001') + 
	outer(DF$year == '2006', DF$year == '2006')
# get the spatial distance between all records
D = as.matrix(dist(DF[,c('easting','northing')]))

# random effects design matrices
Z1 = model.matrix(~ -1 + sample, data = DF)
Z2 = model.matrix(~ -1 + sample:field_dup, data = DF)
Z2 = Z2[,apply(Z2,2,sum) > 0]

# fixed effects design matrix
X = model.matrix(~ year + dist2road + sideroad:dist2road, data = DF)
# get fitted covariance matrix
Sigma = covmatrix(PbOut)
Sigmai = solve(Sigma)
# Projection matrix
P = Sigmai - Sigmai %*% X %*% solve(t(X) %*% Sigmai %*% X, t(X)) %*% Sigmai

# derivative of the covariance matrix with respect to partial sill parameter
dpsill = theta_sp['de']*exp(-D/theta_sp['range'])*part_mask
# derivative of the covariance matrix with respect to nugget parameter
dnugg = diag(n)
# derivative of the covariance matrix with respect to range parameter
drange = D*exp(-D/theta_sp['range'])*part_mask/theta_sp['range']^2
# derivative of the covariance matrix with respect to random effects
dre1 = Z1 %*% t(Z1)
dre2 = Z2 %*% t(Z2)
# observed Fisher Information matrix
I_R = matrix(NA, nrow = 5, ncol = 5)
I_R[1,1] = I_R[1,1] = 0.5*stats:::Tr(P %*% dpsill %*% P %*% dpsill)
I_R[1,2] = I_R[2,1] = 0.5*stats:::Tr(P %*% dpsill %*% P %*% dnugg)
I_R[1,3] = I_R[3,1] = 0.5*stats:::Tr(P %*% dpsill %*% P %*% drange)
I_R[1,4] = I_R[4,1] = 0.5*stats:::Tr(P %*% dpsill %*% P %*% dre1)
I_R[1,5] = I_R[5,1] = 0.5*stats:::Tr(P %*% dpsill %*% P %*% dre2)
I_R[2,2] = I_R[2,2] = 0.5*stats:::Tr(P %*% dnugg %*% P %*% dnugg)
I_R[2,3] = I_R[3,2] = 0.5*stats:::Tr(P %*% dnugg %*% P %*% drange)
I_R[2,4] = I_R[4,2] = 0.5*stats:::Tr(P %*% dnugg %*% P %*% dre1)
I_R[2,5] = I_R[5,2] = 0.5*stats:::Tr(P %*% dnugg %*% P %*% dre2)
I_R[3,3] = I_R[3,3] = 0.5*stats:::Tr(P %*% drange %*% P %*% drange)
I_R[3,4] = I_R[4,3] = 0.5*stats:::Tr(P %*% drange %*% P %*% dre1)
I_R[3,5] = I_R[5,3] = 0.5*stats:::Tr(P %*% drange %*% P %*% dre2)
I_R[4,4] = I_R[4,4] = 0.5*stats:::Tr(P %*% dre1 %*% P %*% dre1)
I_R[4,5] = I_R[5,4] = 0.5*stats:::Tr(P %*% dre1 %*% P %*% dre2)
I_R[5,5] = I_R[5,5] = 0.5*stats:::Tr(P %*% dre1 %*% P %*% dre2)

# create one covariance parameter vector
theta = c(theta_sp, theta_re)
# estimated covariance matrix of covariance parameters
theta_cov = solve(I_R)

# cholesky decomposition of estimated covariance matrix of 
# covariance parameters to be used to simulate values
tcov_chol = chol(theta_cov)

# prepare fixed quantities for predictions in strata
DFp = MOSSpreds
DFp$easting = st_coordinates(DFp)[,1]/1e+3
DFp$northing = st_coordinates(DFp)[,2]/1e+3
DFp = st_drop_geometry(DFp)
DFp$dist2road = log(DFp$dist2road)
DFp01 = DFp
DFp06 = DFp
DFp01[,'year'] = 2001
DFp06[,'year'] = 2006
DFpy = rbind(DFp01,DFp06)
DFpy$year = as.factor(DFpy$year)
	
#strata1
DFp1 = DFpy[DFpy$strat == 1,]
DF1 = rbind(DF[,c('easting','northing','year')],
	DFp1[,c('easting','northing','year')])
D1 = as.matrix(dist(DF1[,c('easting','northing')]))
part_mask1 = outer(DF1$year == '2001', DF1$year == '2001') + 
	outer(DF1$year == '2006', DF1$year == '2006')
n1 = dim(D1)[1]

#strata2
DFp2 = DFpy[DFpy$strat == 2,]
DF2 = rbind(DF[,c('easting','northing','year')],
	DFp2[,c('easting','northing','year')])
D2 = as.matrix(dist(DF2[,c('easting','northing')]))
part_mask2 = outer(DF2$year == '2001', DF2$year == '2001') + 
	outer(DF2$year == '2006', DF2$year == '2006')
n2 = dim(D2)[1]

#strata3
DFp3 = DFpy[DFpy$strat == 3,]
DF3 = rbind(DF[,c('easting','northing','year')],
	DFp3[,c('easting','northing','year')])
D3 = as.matrix(dist(DF3[,c('easting','northing')]))
part_mask3 = outer(DF3$year == '2001', DF3$year == '2001') + 
	outer(DF3$year == '2006', DF3$year == '2006')
n3 = dim(D3)[1]

#strata4
DFp4 = DFpy[DFpy$strat == 4,]
DF4 = rbind(DF[,c('easting','northing','year')],
	DFp4[,c('easting','northing','year')])
D4 = as.matrix(dist(DF4[,c('easting','northing')]))
part_mask4 = outer(DF4$year == '2001', DF4$year == '2001') + 
	outer(DF4$year == '2006', DF4$year == '2006')
n4 = dim(D4)[1]

#strata5
DFp5 = DFpy[DFpy$strat == 5,]
DF5 = rbind(DF[,c('easting','northing','year')],
	DFp5[,c('easting','northing','year')])
D5 = as.matrix(dist(DF5[,c('easting','northing')]))
part_mask5 = outer(DF5$year == '2001', DF5$year == '2001') + 
	outer(DF5$year == '2006', DF5$year == '2006')
n5 = dim(D5)[1]

#iterations start here
set.seed(9001)
nsim = 500
strat1_all = matrix(NA, nrow = dim(DFp1)[1], ncol = nsim)
strat2_all = matrix(NA, nrow = dim(DFp2)[1], ncol = nsim)
strat3_all = matrix(NA, nrow = dim(DFp3)[1], ncol = nsim)
strat4_all = matrix(NA, nrow = dim(DFp4)[1], ncol = nsim)
strat5_all = matrix(NA, nrow = dim(DFp5)[1], ncol = nsim)
k = 1
for(k in 1:nsim) {
	cat("\r", "Simulation: ", k)

  # draw covariance parameter using information matrix as covariance matrix
	get_pos = 0
	while(get_pos == 0) {
		theta_k = theta + t(tcov_chol) %*% rnorm(5)
		if(all(theta_k > 0)) get_pos = 1
	}

	# form covariance matrix from random draw of covariance parameters
	Sigma_k =  theta_k[1]*exp(-D/theta_k[3])*part_mask + 
		theta_k[2]*diag(n) + theta_k[4]*Z1%*%t(Z1) + theta_k[5]*Z2%*%t(Z2)
	Sigma_ki = solve(Sigma_k)

	#draw fixed effects parameters given covariance matrix
	covb = solve(t(X) %*% Sigma_ki %*% X)
	beta = covb %*% t(X) %*% Sigma_ki %*% DF$Pb
	betacov_chol = chol(covb)
	beta_k = beta + t(betacov_chol) %*% rnorm(length(beta))
  
	#strat 1
	# covariance matrix for data and prediction locations without nugget,
	# years are independent by using part_mask
	Sigma_k1 = theta_k[1]*exp(-D1/theta_k[3])*part_mask1
	# subset to prediction locations, add nugget effect
	Sigma_k1pp = Sigma_k1[(n+1):n1,(n+1):n1] + sum(theta_k[c(2,4,5)])*diag(n1 - n)
	# subset to covariance between observed and prediction locations
	Sigma_k1op = Sigma_k1[1:n,(n+1):n1]
	# design matrix for prediction locations
	Xp1 = model.matrix(~ year + dist2road + sideroad:dist2road, data = DFp1)
	# conditional mean given the sampled beta and covariance matrix
	mup_k1 = Xp1 %*% beta_k + t(Sigma_k1op) %*% Sigma_ki %*% 
		(DF$Pb - X %*% beta_k)
	# conditional variance given the sampled beta and covariance parameters
	# equation 9.3, broken into several parts
	Sigp_k1 = Sigma_k1pp - t(Sigma_k1op) %*% Sigma_ki %*% Sigma_k1op
	XD = Xp1 - t(Sigma_k1op) %*% Sigma_ki %*% X
	Sigp_k1 = Sigp_k1 + XD %*% covb %*% t(XD)
	# simulate a realization
	simstrat1 = mup_k1 + t(chol(Sigp_k1)) %*% rnorm(dim(Sigp_k1)[1])
	strat1_all[,k] = simstrat1

	
	#strat2
	Sigma_k2 = theta_k[1]*exp(-D2/theta_k[3])*part_mask2
	Sigma_k2pp = Sigma_k2[(n+1):n2,(n+1):n2] + sum(theta_k[c(2,4,5)])*diag(n2 - n)
	Sigma_k2op = Sigma_k2[1:n,(n+1):n2]
	Xp2 = model.matrix(~ year + dist2road + sideroad:dist2road, data = DFp2)
	mup_k2 = Xp2 %*% beta_k + t(Sigma_k2op) %*% Sigma_ki %*% (DF$Pb - X %*% beta_k)
	Sigp_k2 = Sigma_k2pp - t(Sigma_k2op) %*% Sigma_ki %*% Sigma_k2op
	XD = Xp2 - t(Sigma_k2op) %*% Sigma_ki %*% X
	Sigp_k2 = Sigp_k2 + XD %*% covb %*% t(XD)
	simstrat2 = mup_k2 + t(chol(Sigp_k2)) %*% rnorm(dim(Sigp_k2)[1])
	strat2_all[,k] = simstrat2
		
	#strat3
	Sigma_k3 = theta_k[1]*exp(-D3/theta_k[3])*part_mask3
	Sigma_k3pp = Sigma_k3[(n+1):n3,(n+1):n3] + sum(theta_k[c(2,4,5)])*diag(n3 - n)
	Sigma_k3op = Sigma_k3[1:n,(n+1):n3]
	Xp3 = model.matrix(~ year + dist2road + sideroad:dist2road, data = DFp3)
	mup_k3 = Xp3 %*% beta_k + t(Sigma_k3op) %*% Sigma_ki %*% (DF$Pb - X %*% beta_k)
	Sigp_k3 = Sigma_k3pp - t(Sigma_k3op) %*% Sigma_ki %*% Sigma_k3op
	XD = Xp3 - t(Sigma_k3op) %*% Sigma_ki %*% X
	Sigp_k3 = Sigp_k3 + XD %*% covb %*% t(XD)
	simstrat3 = mup_k3 + t(chol(Sigp_k3)) %*% rnorm(dim(Sigp_k3)[1])
	strat3_all[,k] = simstrat3

	#strat4
	Sigma_k4 = theta_k[1]*exp(-D4/theta_k[3])*part_mask4
	Sigma_k4pp = Sigma_k4[(n+1):n4,(n+1):n4] + sum(theta_k[c(2,4,5)])*diag(n4 - n)
	Sigma_k4op = Sigma_k4[1:n,(n+1):n4]
	Xp4 = model.matrix(~ year + dist2road + sideroad:dist2road, data = DFp4)
	mup_k4 = Xp4 %*% beta_k + t(Sigma_k4op) %*% Sigma_ki %*% (DF$Pb - X %*% beta_k)
	Sigp_k4 = Sigma_k4pp - t(Sigma_k4op) %*% Sigma_ki %*% Sigma_k4op
	XD = Xp4 - t(Sigma_k4op) %*% Sigma_ki %*% X
	Sigp_k4 = Sigp_k4 + XD %*% covb %*% t(XD)
	simstrat4 = mup_k4 + t(chol(Sigp_k4)) %*% rnorm(dim(Sigp_k4)[1])
	strat4_all[,k] = simstrat4

	#strat5
	Sigma_k5 = theta_k[1]*exp(-D5/theta_k[3])*part_mask5
	Sigma_k5pp = Sigma_k5[(n+1):n5,(n+1):n5] + sum(theta_k[c(2,4,5)])*diag(n5 - n)
	Sigma_k5op = Sigma_k5[1:n,(n+1):n5]
	Xp5 = model.matrix(~ year + dist2road + sideroad:dist2road, data = DFp5)
	mup_k5 = Xp5 %*% beta_k + t(Sigma_k5op) %*% Sigma_ki %*% (DF$Pb - X %*% beta_k)
	Sigp_k5 = Sigma_k5pp - t(Sigma_k5op) %*% Sigma_ki %*% Sigma_k5op
	XD = Xp5 - t(Sigma_k5op) %*% Sigma_ki %*% X
	Sigp_k5 = Sigp_k5 + XD %*% covb %*% t(XD)
	simstrat5 = mup_k5 + t(chol(Sigp_k5)) %*% rnorm(dim(Sigp_k5)[1])
	strat5_all[,k] = simstrat5

}

Sigp_k1[1:5,1:5]
plot(strat4_all[6,], type = 'l')

save(strat1_all,file = 'strat1_all.rda')
save(strat2_all,file = 'strat2_all.rda')
save(strat3_all,file = 'strat3_all.rda')
save(strat4_all,file = 'strat4_all.rda')
save(strat5_all,file = 'strat5_all.rda')
# load('strat1_all.rda')
# load('strat2_all.rda')
# load('strat3_all.rda')
# load('strat4_all.rda')
# load('strat5_all.rda')

np1 = dim(strat1_all)[1]/2
np2 = dim(strat2_all)[1]/2
np3 = dim(strat3_all)[1]/2
np4 = dim(strat4_all)[1]/2
np5 = dim(strat5_all)[1]/2

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        Histogram of yearly differences at prediction locations
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

strat1_yeardiff = apply(exp(strat1_all[1:np1,]) - 
	exp(strat1_all[(np1 + 1):(2*np1),]),2,mean)
strat2_yeardiff = apply(exp(strat2_all[1:np2,]) - 
	exp(strat2_all[(np2 + 1):(2*np2),]),2,mean)
strat3_yeardiff = apply(exp(strat3_all[1:np3,]) - 
	exp(strat3_all[(np3 + 1):(2*np3),]),2,mean)
strat4_yeardiff = apply(exp(strat4_all[1:np4,]) - 
	exp(strat4_all[(np4 + 1):(2*np4),]),2,mean)
strat5_yeardiff = apply(exp(strat5_all[1:np5,]) - 
	exp(strat5_all[(np5 + 1):(2*np5),]),2,mean)

file_name = 'figures/Moss_diffhist'
pdf(paste0(file_name,'.pdf'), width = 6, height = 10)

layout(matrix(1:5, ncol = 1))
par(mar = c(5,5,5,1))
hist(strat1_yeardiff, main = 'Stratum 1', cex.lab = 2,
	cex.axis = 1.5, xlab = '2001 Pb conc. - 2006 Pb conc.',
	cex.main = 2)
hist(strat2_yeardiff, main = 'Stratum 2', cex.lab = 2,
	cex.axis = 1.5, xlab = '2001 Pb conc. - 2006 Pb conc.',
	cex.main = 2)
hist(strat3_yeardiff, main = 'Stratum 3', cex.lab = 2,
	cex.axis = 1.5, xlab = '2001 Pb conc. - 2006 Pb conc.',
	cex.main = 2)
hist(strat4_yeardiff, main = 'Stratum 4', cex.lab = 2,
	cex.axis = 1.5, xlab = '2001 Pb conc. - 2006 Pb conc.',
	cex.main = 2)
hist(strat5_yeardiff, main = 'Stratum 5', cex.lab = 2,
	cex.axis = 1.5, xlab = '2001 Pb conc. - 2006 Pb conc.',
	cex.main = 2)
layout(1)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

2*(1 - pnorm(mean(strat1_yeardiff)/sqrt(var(strat1_yeardiff))))
2*(1 - pnorm(mean(strat2_yeardiff)/sqrt(var(strat2_yeardiff))))
2*(1 - pnorm(mean(strat3_yeardiff)/sqrt(var(strat3_yeardiff))))
2*(1 - pnorm(mean(strat4_yeardiff)/sqrt(var(strat4_yeardiff))))

sum(strat5_yeardiff > 0)
2*(1 - sum(strat5_yeardiff > 0)/500)
2*(1 - pnorm(mean(strat5_yeardiff)/sqrt(var(strat5_yeardiff))))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Map of Percent Change from 2001 to 2006
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = 'figures/Moss_diffmaps'
pdf(paste0(file_name,'.pdf'), width = 8, height = 11)

layout(matrix(c(1,1,2,3), nrow = 2, byrow = TRUE), heights = c(1.8, 1))
#stratum 4
preds01 = exp(strat4_all[1:np4,])
preds06 = exp(strat4_all[(np4 + 1):(2*np4),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
quantile(cbind(preddiffmn, preddiffup, preddifflo), probs = (0:9)/9)
brks_pred = c(-(6:0)/6.5,0.3, 0.7, 2)
cip = classIntervals(preddiffmn, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
plot(DFp4[,c('easting','northing')], col = cip_colors, pch = 19, cex = 2,
	cex.lab = 2, cex.axis = 1.5, main = 'Mean Proportional Change from 2001 to 2006',
	cex.main = 2)

#stratum 3
preds01 = exp(strat3_all[1:np3,])
preds06 = exp(strat3_all[(np3 + 1):(2*np3),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
quantile(cbind(preddiffmn, preddiffup, preddifflo), probs = (0:9)/9)
brks_pred = c(-(6:0)/6,0.3, 0.7, 2)
cip = classIntervals(preddiffmn, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
points(DFp3[,c('easting','northing')], col = cip_colors, pch = 19, cex = 1.5)

#stratum 2
preds01 = exp(strat2_all[1:np2,])
preds06 = exp(strat2_all[(np2 + 1):(2*np2),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
quantile(cbind(preddiffmn, preddiffup, preddifflo), probs = (0:9)/9)
brks_pred = c(-(6:0)/6,0.3, 0.7, 2)
cip = classIntervals(preddiffmn, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
points(DFp2[,c('easting','northing')], col = cip_colors, pch = 19, cex = 1)

#stratum 1
preds01 = exp(strat1_all[1:np1,])
preds06 = exp(strat1_all[(np1 + 1):(2*np1),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
quantile(cbind(preddiffmn, preddiffup, preddifflo), probs = (0:9)/9)
brks_pred = c(-(6:0)/6,0.3, 0.7, 2)
cip = classIntervals(preddiffmn, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
points(DFp1[,c('easting','northing')], col = cip_colors, pch = 19, cex = .3)

addBreakColorLegend(xleft = -402, ybottom = 1986, xright = -400, ytop = 1999,
			breaks = brks_pred, colors = palp, cex = 1.5, printFormat = "4.3")

######## B ########

#stratum 4
preds01 = exp(strat4_all[1:np4,])
preds06 = exp(strat4_all[(np4 + 1):(2*np4),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
brks_pred = c(-(6:0)/6.5,0.3, 0.7, 2)
cip = classIntervals(preddifflo, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
par(mar = c(0,1,5,1))
plot(DFp4[,c('easting','northing')], type = 'n', bty = 'n', xlab = '',
	ylab = '', xaxt = 'n', yaxt = 'n', main = 'Lower 90% Bound',
	cex.main = 2)
points(DFp4[,c('easting','northing')], col = cip_colors, pch = 19, cex = 1.0)

#stratum 3
preds01 = exp(strat3_all[1:np3,])
preds06 = exp(strat3_all[(np3 + 1):(2*np3),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
brks_pred = c(-(6:0)/6.5,0.3, 0.7, 2)
cip = classIntervals(preddifflo, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
points(DFp3[,c('easting','northing')], col = cip_colors, pch = 19, cex = 0.7)

#stratum 2
preds01 = exp(strat2_all[1:np2,])
preds06 = exp(strat2_all[(np2 + 1):(2*np2),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
quantile(cbind(preddiffmn, preddiffup, preddifflo), probs = (0:9)/9)
brks_pred = c(-(6:0)/6,0.3, 0.7, 2)
cip = classIntervals(preddifflo, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
points(DFp2[,c('easting','northing')], col = cip_colors, pch = 19, cex = .5)

#stratum 1
preds01 = exp(strat1_all[1:np1,])
preds06 = exp(strat1_all[(np1 + 1):(2*np1),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
quantile(cbind(preddiffmn, preddiffup, preddifflo), probs = (0:9)/9)
brks_pred = c(-(6:0)/6,0.3, 0.7, 2)
cip = classIntervals(preddifflo, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
points(DFp1[,c('easting','northing')], col = cip_colors, pch = 19, cex = .1)


######## C ########

#stratum 4
preds01 = exp(strat4_all[1:np4,])
preds06 = exp(strat4_all[(np4 + 1):(2*np4),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
brks_pred = c(-(6:0)/6.5,0.3, 0.7, 2)
cip = classIntervals(preddiffup, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
par(mar = c(0,1,5,1))
plot(DFp4[,c('easting','northing')], type = 'n', bty = 'n', xlab = '',
	ylab = '', xaxt = 'n', yaxt = 'n', main = 'Upper 90% Bound',
	cex.main = 2)
points(DFp4[,c('easting','northing')], col = cip_colors, pch = 19, cex = 1.0)

#stratum 3
preds01 = exp(strat3_all[1:np3,])
preds06 = exp(strat3_all[(np3 + 1):(2*np3),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
brks_pred = c(-(6:0)/6.5,0.3, 0.7, 2)
cip = classIntervals(preddiffup, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
par(mar = c(5,5,5,1))
points(DFp3[,c('easting','northing')], col = cip_colors, pch = 19, cex = 0.7)

#stratum 2
preds01 = exp(strat2_all[1:np2,])
preds06 = exp(strat2_all[(np2 + 1):(2*np2),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
quantile(cbind(preddiffmn, preddiffup, preddifflo), probs = (0:9)/9)
brks_pred = c(-(6:0)/6,0.3, 0.7, 2)
cip = classIntervals(preddiffup, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
points(DFp2[,c('easting','northing')], col = cip_colors, pch = 19, cex = .5)

#stratum 1
preds01 = exp(strat1_all[1:np1,])
preds06 = exp(strat1_all[(np1 + 1):(2*np1),])
preddiff = (preds06 - preds01)/preds01
preddiffmn = apply(preddiff,1,mean)
preddiffup = apply(preddiff,1,quantile, prob = 0.95)
preddifflo = apply(preddiff,1,quantile, prob = 0.05)
quantile(cbind(preddiffmn, preddiffup, preddifflo), probs = (0:9)/9)
brks_pred = c(-(6:0)/6,0.3, 0.7, 2)
cip = classIntervals(preddiffup, style = 'fixed', 
		fixedBreaks= brks_pred)
palp = viridis(12)[c(1:6,10:12)]
cip_colors = findColours(cip, palp)
points(DFp1[,c('easting','northing')], col = cip_colors, pch = 19, cex = .1)

layout(1)

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
#           Proportion Above a Threshold
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

p_above_thresh = rbind(
c(
quantile(apply((exp(strat1_all[1:np1,]) > 55)*1, 2, mean), prob = c(0.025, 0.5, 0.975)),
quantile(apply((exp(strat1_all[(np1 + 1):(2*np1),]) > 55)*1, 2, mean), prob = c(0.025, 0.5, 0.975))),
c(
quantile(apply((exp(strat2_all[1:np1,]) > 55)*1, 2, mean), prob = c(0.025,0.5, 0.975)),
quantile(apply((exp(strat2_all[(np2 + 1):(2*np2),]) > 55)*1, 2, mean), prob = c(0.025, 0.5, 0.975))),
c(
quantile(apply((exp(strat3_all[1:np1,]) > 55)*1, 2, mean), prob = c(0.025, 0.5, 0.975)),
quantile(apply((exp(strat3_all[(np3 + 1):(2*np3),]) > 55)*1, 2, mean), prob = c(0.025, 0.5, 0.975)))
)

print(
    xtable(p_above_thresh, 
      align = c('l',rep('l', times = length(p_above_thresh[1,]))),
      digits = c(0, rep(3, times = 6)),
      caption = 'Metrics for model selection',
      label = 'tab:ModSelMetrics'
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
