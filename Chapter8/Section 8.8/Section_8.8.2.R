sec_path = 'Rcode/Chapter8/Section 8.8/'
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
library(Matrix)
library(xtable)
library(spmodel)
library(vioplot)
library(numDeriv)
source('addBreakColorLegend.R')

# load data for graphics and analysis
data(sealPolys)
seals_sf = st_as_sf(sealPolys)
seals_sf$stockname = as.factor(as.character(seals_sf$stockname))

################################################################################
#-------------------------------------------------------------------------------
#                    Create Neighborhood Matrices
#-------------------------------------------------------------------------------
################################################################################

nTot = length(sealPolys)
nObs = sum(!is.na(sealPolys$Estimate))
nMiss = nTot - nObs

# a function to create matrix from lists and numbers of neighbors
Neighmat <- function(adj, num, n)
{
	N.mat <- matrix( 0, nrow = n, ncol = n )
	k <- 1
	for (i in 1:n){
		if(num[i] > 0) {
				N.mat[i,adj[[i]]] <- 1
			}
		}
	N.mat
}

# use spdep to find touching neighors, and then add rest manually
Nlist = poly2nb(sealPolys, snap = 2000)
Nlist[[79]] = as.integer(c(211, 463))
Nlist[[211]] = as.integer(c(Nlist[[211]]), 79)
Nlist[[463]] = as.integer(c(Nlist[[463]]), 79)
Nlist[[130]] = as.integer(302)
Nlist[[302]] = as.integer(c(Nlist[[302]],130))
Nlist[[325]] = as.integer(c(326, 353))
Nlist[[326]] = as.integer(c(Nlist[[326]],325))
Nlist[[353]] = as.integer(c(Nlist[[353]],325))
Nlist[[435]] = as.integer(437)
Nlist[[437]] = as.integer(c(Nlist[[437]],435))
Nlist[[436]] = as.integer(c(86, 88))
Nlist[[86]] = as.integer(c(Nlist[[86]],436))
Nlist[[88]] = as.integer(c(Nlist[[88]],436))
Nlist[[437]] = as.integer(87)
Nlist[[87]] = as.integer(c(Nlist[[87]],437))
Nlist[[438]] = as.integer(436)
Nlist[[436]] = as.integer(c(Nlist[[436]],438))
Nlist[[439]] = as.integer(346)
Nlist[[346]] = as.integer(c(Nlist[[346]],439))
Nlist[[443]] = as.integer(281)
Nlist[[281]] = as.integer(c(Nlist[[281]],443))
Nlist[[463]] = as.integer(79)
attr(Nlist,'polyid') = as.factor(as.character(sealPolys@data$polyid))
attr(Nlist,'stockid') = as.factor(as.character(sealPolys@data$stockid))
num = lapply(Nlist, function(x) length(x))
num = unlist(num)
Nmat = Neighmat(Nlist, num, length(num))
Nmat1 = pmax(Nmat,t(Nmat))
Nmat2 = (((Nmat1 %*% Nmat1 > 0)*1 + Nmat1) > 0)*1 - diag(dim(Nmat1)[1])
Nmat3 = (((Nmat2 %*% Nmat1 > 0)*1 + Nmat2) > 0)*1 - diag(dim(Nmat1)[1])
Nmat4 = (((Nmat3 %*% Nmat1 > 0)*1 + Nmat3) > 0)*1 - diag(dim(Nmat1)[1])
Nmat5 = (((Nmat4 %*% Nmat1 > 0)*1 + Nmat4) > 0)*1 - diag(dim(Nmat1)[1])
Nmat6 = (((Nmat5 %*% Nmat1 > 0)*1 + Nmat5) > 0)*1 - diag(dim(Nmat1)[1])

polycentroids = st_coordinates(st_centroid(seals_sf$geometry))
rownames(polycentroids) = rownames(seals_sf)
distMat = as.matrix(dist(polycentroids))/1000


################################################################################
#-------------------------------------------------------------------------------
#                  CAR and SAR Likelihoods
#-------------------------------------------------------------------------------
################################################################################

# reproduce parts of Figure 5 in Ver Hoef et al. 2018
# investigate CAR versus SAR, and row-standardized versus binary weights

formulas = c(Estimate ~ 1, Estimate ~ stockname)
Ws = list(Nmat1, Nmat2, Nmat4)
formlabs = c('m','X')
Wlabs = c('1','2','4')
mods = c('car','sar')
rowstand = c(TRUE, FALSE)
store_results = matrix(nrow = 24, ncol = 7)
ii = 0
for(i in 1:2) {
	for(j in 1:3) {
		for(k in 1:2) {
			for(m in 1:2) {
				spfit = spautor( data = seals_sf, estmethod = 'ml', 
						control = list(reltol = 1e-7),
						formulas[[i]], 
						W = Ws[[j]], 
						spcov_type = mods[k],
						row_st = rowstand[m] )
				ii = ii + 1
				store_results[ii,1] = 3*(i - 1) + j
				store_results[ii,2] = i
				store_results[ii,3] = j
				store_results[ii,4] = k
				store_results[ii,5] = m
				store_results[ii,6] = -2*logLik(spfit)
				store_results[ii,7] = AIC(spfit)
			}
		}
	}
}
# fix AIC results
store_results[1:12,7] = store_results[1:12,6] + 6
store_results[13:24,7] = store_results[13:24,6] + 14
store_results[1:12,1] = store_results[1:12,1] + 1
store_results[13:24,1] = store_results[13:24,1] + 2
lmout_m = splm(Estimate ~ 1, data = seals_sf, spcov_type = 'none', 
	estmeth = 'ml')
lmout_X = splm(Estimate ~ stockname, data = seals_sf, spcov_type = 'none', 
	estmeth = 'ml')
store_jitter = store_results
store_jitter[,1] = store_jitter[,1] + ((1:4) - 2.5)/10

file_name = "figures/Seals_m2LL_AIC"
pdf(paste0(file_name,'.pdf'), width = 12, height = 6)

	layout(matrix(1:2, nrow = 1))
	labs = c('m0', 'm1','m2', 'm4', 'X0','X1', 'X2', 'X4')
	sar_col = '#4daf4a'
	car_col = '#e41a1c'
	padj = 0
	adj = -.2
	cex_mtext = 3.2
	cex_all = 1.8
	par(mar = c(5,6,4,1))
	plot(store_jitter[,c(1,6)], xlim = c(.9, max(store_jitter[,1]) + .1),
		ylim = c(min(-2*logLik(lmout_m), -2*logLik(lmout_X), store_jitter[,6]),
			max(-2*logLik(lmout_m), -2*logLik(lmout_X), store_jitter[,6]) + 10),
		type = 'n', xlab = '', xaxt = 'n',
		ylab = expression("-2"*italic(L)(tilde(bold(theta))~";"~bold(y))), cex.lab = 2, cex.axis = 1.5)
	points(1, -2*logLik(lmout_m), pch = 19, cex = cex_all)
	points(5, -2*logLik(lmout_X), pch = 19, cex = cex_all)
	ind = store_jitter[,4] == 1 & store_jitter[,5] == 1
	points(store_jitter[ind,1],store_jitter[ind,6], col = car_col, 
		pch = 15, cex = cex_all)
	ind = store_jitter[,4] == 1 & store_jitter[,5] == 2
	points(store_jitter[ind,1],store_jitter[ind,6], col = car_col, 
		pch = 19, cex = cex_all)
	ind = store_jitter[,4] == 2 & store_jitter[,5] == 1
	points(store_jitter[ind,1],store_jitter[ind,6], col = sar_col, 
		pch = 15, cex = cex_all)
	ind = store_jitter[,4] == 2 & store_jitter[,5] == 2
	points(store_jitter[ind,1],store_jitter[ind,6], col = sar_col, 
		pch = 19, cex = cex_all)
	axis(1, at = 1:8, labels = labs, las = 1, cex.axis = 1.5)
	legend(4.3,-338, legend = 
		c('Independence','CAR unstandardized','CAR row-standard',
			'SAR unstandardized','SAR row-standard'),
		pch = c(19, 19, 15, 19, 15), col = c('black',car_col, car_col,
			sar_col, sar_col), cex = 1.3)
	mtext('A', adj = adj, cex = cex_mtext, padj = padj)

	par(mar = c(5,5,4,1))
	plot(store_jitter[,c(1,6)], xlim = c(1, max(store_jitter[,1])),
		ylim = c(min(-2*logLik(lmout_m), -2*logLik(lmout_X), store_jitter[,6]),
			max(-2*logLik(lmout_m), -2*logLik(lmout_X), store_jitter[,6]) + 10),
		type = 'n', xlab = '', xaxt = 'n',
		ylab = 'AIC', cex.lab = 2, cex.axis = 1.5)
	points(1, AIC(lmout_m), pch = 19, cex = cex_all)
	points(5, AIC(lmout_X), pch = 19, cex = cex_all)
	ind = store_jitter[,4] == 1 & store_jitter[,5] == 1
	points(store_jitter[ind,1],store_jitter[ind,7], col = car_col, 
		pch = 15, cex = cex_all)
	ind = store_jitter[,4] == 1 & store_jitter[,5] == 2
	points(store_jitter[ind,1],store_jitter[ind,7], col = car_col, 
		pch = 19, cex = cex_all)
	ind = store_jitter[,4] == 2 & store_jitter[,5] == 1
	points(store_jitter[ind,1],store_jitter[ind,7], col = sar_col, 
		pch = 15, cex = cex_all)
	ind = store_jitter[,4] == 2 & store_jitter[,5] == 2
	points(store_jitter[ind,1],store_jitter[ind,7], col = sar_col, 
		pch = 19, cex = cex_all)
	axis(1, at = 1:8, labels = labs, las = 1, cex.axis = 1.5)
	mtext('B', adj = adj, cex = cex_mtext, padj = padj)
	
	layout(1)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

################################################################################
#-------------------------------------------------------------------------------
#          Should Outlier be Removed?
#-------------------------------------------------------------------------------
################################################################################

summary(spautor(Estimate ~ stockname, data = seals_sf, estmethod = 'ml',  
	W = Nmat1, spcov_type = 'car', row_st = TRUE))

which(seals_sf$Estimate > .7)

summary(spautor(Estimate ~ stockname, 
	data = seals_sf[!(1:dim(seals_sf)[1] == 135), ],
	estmethod = 'ml', W = Nmat1, spcov_type = 'car', row_st = TRUE))
						
################################################################################
#-------------------------------------------------------------------------------
#          Likelihood Surface and Profile Likelihood
#-------------------------------------------------------------------------------
################################################################################

spfit = spautor(Estimate ~ stockname, data = seals_sf, estmethod = 'reml', 
	control = list(reltol = 1e-7), W = Nmat4, spcov_type = 'car', row_st = TRUE )
summary(spfit)

#set grids sequences and dimensions
seqs = seq(-1,1,by = 0.05)
grid_dim = length(seqs)

# create a sequence of values to try for range and partial sill
de_seq = log(coef(spfit, type = "spcov")['de']) + .2*seqs
de_seq = exp(de_seq)
range_seq = 0.99*(seqs + 1)/2

zA = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		spautorout = spautor(
			Estimate ~ stockname, data = seals_sf, estmethod = 'reml',
			control = list(reltol = 1e-7), W = Nmat4, row_st = TRUE,
			spcov_initial = spcov_initial(spcov_type = 'car', 
				de = de_seq[i], range = range_seq[j],
				known = c('de','range'))
		)	
		zA[i,j] = logLik(spautorout)
	}
}

store_range = rep(NA, times = grid_dim)
for(i in 1:grid_dim) {
		spautorout = spautor(
			Estimate ~ stockname, data = seals_sf, estmethod = 'reml',
			control = list(reltol = 1e-7), W = Nmat4, row_st = TRUE,
			spcov_initial = spcov_initial(spcov_type = 'car', 
				range = range_seq[i],
				known = c('range'))
		)	
	store_range[i] = logLik(spautorout)
}

store_de = rep(NA, times = grid_dim)
for(i in 1:grid_dim) {
		spautorout = spautor(
			Estimate ~ stockname, data = seals_sf, estmethod = 'reml',
			control = list(reltol = 1e-7), W = Nmat4, row_st = TRUE,
			spcov_initial = spcov_initial(spcov_type = 'car', 
				de = de_seq[i],
				known = c('de'))
		)	
	store_de[i] = logLik(spautorout)
}

# use profile likelihood to get confidence interval on autocorrelation parameter	
minrhoindx = min(which(2*store_range > 
	(2*logLik(spfit) - qchisq(0.95, df = 1))))
maxrhoindx = max(which(2*store_range > 
	(2*logLik(spfit) - qchisq(0.95, df = 1))))

# linear interpolation for bounds of confidence interval
lo_range = mean(range_seq[minrhoindx],range_seq[minrhoindx-1])
up_range = mean(range_seq[maxrhoindx],range_seq[maxrhoindx+1])
lo_range
up_range

# use profile likelihood to get confidence interval on variance parameter	
minsigindx = min(which(2*store_de > 
	(2*logLik(spfit) - qchisq(0.95, df = 1))))
maxsigindx = max(which(2*store_de > 
	(2*logLik(spfit) - qchisq(0.95, df = 1))))

# linear interpolation for bounds of confidence interval
lo_sig = mean(de_seq[minsigindx],de_seq[minsigindx-1])
up_sig = mean(de_seq[maxsigindx],de_seq[maxsigindx+1])
lo_sig
up_sig

file_name = 'figures/Seals_logLik'
#pdf(paste0(file_name,'.pdf'), width = 12, height = 6)
tiff(paste0(file_name,'.tiff'), width = 960, height = 480)

	padj = -.5
	adj = -.17
	cex_mtext = 3.3
	cex_lab = 3.5
	cex_axis = 2
	
	layout(matrix(c(1,1,2,2,
									1,1,3,3), nrow = 2, byrow = TRUE))
	nbrks = 20
	brks = quantile(zA, probs = (0:nbrks)/nbrks)
	cramp = viridis(nbrks)
	par(mar = c(7,7,5,2), mgp=c(4, 1.3, 0))
	image(de_seq, range_seq, zA, breaks = brks, col = cramp,
		cex.main = 2, xlab = '', ylab = expression(rho),
		cex.axis = cex_axis, cex.lab = cex_lab)
	title(xlab = expression(sigma^2), line = 5, cex.lab = cex_lab)
	points(coef(spfit, type = "spcov")['de'], 
		coef(spfit, type = "spcov")['range'], 
		pch = 19, cex = 2.5, col = 'black')
	mtext('A', adj = adj, cex = cex_mtext, padj = padj)

	cex_lab = 2.7
	adj = -.13

	par(mar = c(6,9,5,1), mgp=c(4, 1.3, 0))
	plot(de_seq,2*store_de, type = 'l', lwd = 3, 
		xlab = '', 
		ylab = expression("2"*italic(L)[italic("-i,R")](sigma^2~";"~tilde(rho),bold(y))), 
		cex.lab = cex_lab, cex.axis = cex_axis)
	points(coef(spfit, type = "spcov")['de'], 2*logLik(spfit), pch = 19,
		cex = 2.5)
	lines(c(min(de_seq), max(de_seq)), 
		rep(2*logLik(spfit) - qchisq(.95, df = 1), times = 2), 
		lty = 2, lwd = 3)
	lines(c(lo_sig, up_sig), 
		rep(2*logLik(spfit) - qchisq(.95, df = 1), times = 2), 
		lwd = 7)
	title(xlab = expression(sigma^2), line = 5, cex.lab = 3)
	mtext('B', adj = adj, cex = cex_mtext, padj = padj)


	plot(range_seq,2*store_range, type = 'l', lwd = 3, 
		xlab = '', 
		ylab = expression("2"*italic(L)[italic("-i,R")](rho~";"~tilde(sigma)^2,bold(y))), 
		cex.lab = cex_lab, cex.axis = cex_axis)
	points(coef(spfit, type = "spcov")['range'], 2*logLik(spfit), pch = 19,
		cex = 2.5)
	lines(c(-1, 1), rep(2*logLik(spfit) - qchisq(.95, df = 1), times = 2), 
		lty = 2, lwd = 3)
	lines(c(lo_range, up_range), 
		rep(2*logLik(spfit) - qchisq(.95, df = 1), times = 2), 
		lwd = 7)
	title(xlab = expression(rho), line = 5, cex.lab = 3)
	mtext('C', adj = adj, cex = cex_mtext, padj = padj)
		
	layout(1)

dev.off()
	
system(paste0('tiff2pdf -o','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.tiff','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'.tiff','\''))

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))


################################################################################
#-------------------------------------------------------------------------------
# Hessian, Fisher Information, and Covariance Matrix for Covariance Parameters
#-------------------------------------------------------------------------------
################################################################################

spfit = spautor(Estimate ~ stockname, data = seals_sf, estmethod = 'reml', 
	control = list(reltol = 1e-7), W = Nmat4, spcov_type = 'car', row_st = TRUE )
theta = coef(spfit, type = 'spcov')
summary(spfit)

# ---------------------- Analytical Approach ----------------------------

# create the covariance matrix and take its inverse
D = diag(apply(Nmat4,1,sum))
DW = D - theta['range']*Nmat4
DWi = solve(DW)
Vi = DW/theta['de']
V = theta['de']*DWi
A = -theta['de']*DWi %*% Nmat4 %*% DWi
X = model.matrix(~ stockname, data = seals_sf)
P = Vi - Vi %*% X %*% solve(t(X) %*% Vi %*% X, t(X)) %*% Vi
		
# compute Fisher Information element by element using trace formula
FI = matrix(rep(NA, times = 4), nrow = 2)
FI[1,1] = sum(diag(P %*% DWi %*% P %*% DWi))
FI[1,2] = FI[2,1] = sum(diag(P %*% DWi %*% P %*% A))
FI[2,2] = sum(diag(P %*% A %*% P %*% A))
FI = 0.5*FI
asycov = solve(FI)
sqrt(diag(asycov))
theta[c('de','range')] - 1.96*sqrt(diag(asycov))
theta[c('de','range')] + 1.96*sqrt(diag(asycov))

# asymptotic correlation matrix
diag(1/sqrt(diag(asycov))) %*% asycov %*% diag(1/sqrt(diag(asycov)))

# ---------------------- Numerical Hessian Approach ----------------------------

# a function to return log-likelihood for specified theta
fixed_parms = function(theta) {
	#initialize covariance parameters and hold them constant
	spini = spcov_initial = spcov_initial(spcov_type = 'car', de = theta[1], 
		range = theta[2], known = c('de','range'))
	ModelFit = spautor(Estimate ~ stockname, data = seals_sf, estmethod = 'reml', 
		control = list(reltol = 1e-7), W = Nmat4, row_st = TRUE, 
		spcov_initial = spini)
	logLik(ModelFit)
}
fixed_parms(theta)
# find numerical Hessian
FishInf = -hessian(fixed_parms, theta[c(1,3)])
# asymptotic covariance matrix
asycov_comp = solve(FishInf)

sqrt(diag(asycov_comp))

theta[c('de','range')] - 1.96*sqrt(diag(asycov_comp))
theta[c('de','range')] + 1.96*sqrt(diag(asycov_comp))

# asymptotic correlation matrix
diag(1/sqrt(diag(asycov_comp))) %*% asycov_comp %*% 
	diag(1/sqrt(diag(asycov_comp)))

################################################################################
#-------------------------------------------------------------------------------
#          Nonstationarity
#-------------------------------------------------------------------------------
################################################################################

spfit4_rs = spautor(Estimate ~ stockname, data = seals_sf, estmethod = 'reml', 
	control = list(reltol = 1e-7), W = Nmat4, spcov_type = 'car', row_st = TRUE)
rhohat_rs = coef(spfit4_rs, type = 'spcov')['range']
sig2hat_rs = coef(spfit4_rs, type = 'spcov')['de']

spfit4_un = spautor(Estimate ~ stockname, data = seals_sf, estmethod = 'ml', 
	control = list(reltol = 1e-7), W = Nmat4, spcov_type = 'car', row_st = FALSE)
rhohat_un = coef(spfit4_un, type = 'spcov')['range']
sig2hat_un = coef(spfit4_un, type = 'spcov')['de']

# create fitted CAR covariance matrix for 4th order neighbors, both
# row-standardized and unstandardized
n = dim(Nmat4)[1]
Sigma_rs = sig2hat_rs*
    solve(diag(apply(Nmat4,1,sum)) - 
		rhohat_rs*Nmat4)
margvar_rs = diag(Sigma_rs)
plot(apply(Nmat4,1,sum), margvar_rs)

Sigma_un = sig2hat_un*
    solve(diag(dim(Nmat4)[1]) - 
		rhohat_un*Nmat4)
margvar_un = diag(Sigma_un)
plot(apply(Nmat4,1,sum), margvar_un)

# correlation matrix for Sigma_rs
Cormat_rs = diag(1/sqrt(diag(Sigma_rs))) %*% Sigma_rs %*%
	diag(1/sqrt(diag(Sigma_rs)))
cN1 = Cormat_rs[Nmat1 == 1]
cN2 = Cormat_rs[Nmat2 - Nmat1 == 1]
cN3 = Cormat_rs[Nmat3 - Nmat2 == 1]
cN4 = Cormat_rs[Nmat4 - Nmat3 == 1]
cN5 = Cormat_rs[Nmat5 - Nmat4 == 1]
cN6 = Cormat_rs[Nmat6 - Nmat5 == 1]

corNei = rbind(cbind(cN1,1), cbind(cN2,2), cbind(cN3,3), cbind(cN4,4),
	cbind(cN5,5), cbind(cN6,6)) 
corNei = data.frame(cor = corNei[,1], ordNei = corNei[,2])
vioplot(cor ~ ordNei, data = corNei)

dist_cor = data.frame(dist = distMat[upper.tri(distMat)],
	cor = Cormat_rs[upper.tri(Cormat_rs)])
plot(dist_cor, pch = 19, cex = .5, col = rgb(0,0,0,.1), xlim = c(0,200) )

file_name = 'figures/Seals_nonstationary'
#pdf(paste0(file_name,'.pdf'), width = 6, height = 6)
png(paste0(file_name,'.png'), width = 960, height = 960)

	padj = -.25
	adj = -.23
	cex_mtext = 3.5
	cex_lab = 3.5
	cex_axis = 2
	layout(matrix(1:4, nrow = 2, byrow = TRUE))
	# A
	par(mar = c(7,7,5,2), mgp=c(4, 1.3, 0))
	plot(apply(Nmat4,1,sum), margvar_rs, ylab = 'Marginal Variance',
		xlab = 'Number of Neighbors', cex.lab = cex_lab, 
		cex.axis = cex_axis, cex = 1.5)
	mtext('A', adj = adj, cex = cex_mtext, padj = padj)
	# B
	plot(apply(Nmat4,1,sum), margvar_un, ylab = 'Marginal Variance',
		xlab = 'Number of Neighbors', cex.lab = cex_lab, 
		cex.axis = cex_axis, cex = 1.5)
	mtext('B', adj = adj, cex = cex_mtext, padj = padj)
	# C
	vioplot(cor ~ ordNei, data = corNei, xlab = '',
		ylab = '', cex.lab = cex_lab+2, cex.axis = cex_axis)
	title(ylab="Autocorrelation", xlab="Neighbor Order", cex.lab = cex_lab)
	axis(1, at = 1:4, labels = as.character(1:4), las = 1, cex.axis = cex_axis,
		cex.lab = cex_lab, xlab = 'Number of')
	mtext('C', adj = adj, cex = cex_mtext, padj = padj)
	# D
	plot(dist_cor, pch = 19, cex = .5, col = rgb(0,0,0,.1), xlim = c(0,200),
		xlab = 'Distance (km)', ylab = 'Autocorrelation', cex.lab = cex_lab, 
		cex.axis = cex_axis)
	mtext('D', adj = adj, cex = cex_mtext, padj = padj)

	layout(1)

dev.off()
	
system(paste0('img2pdf ','\'',SLEDbook_path,
	sec_path,file_name,'.png','\'', ' -o ', '\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'.png','\''))

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

################################################################################
#-------------------------------------------------------------------------------
#                         Other Models
#-------------------------------------------------------------------------------
################################################################################

# model that allows islands
# let spmodel determine neighbors by those that share any border
spautor_3parms = spautor(Estimate ~ -1 + stockname, data = seals_sf, 
	estmethod = 'ml', spcov_type = 'car', row_st = TRUE,
	control = list(reltol = 1e-8))
summary(spautor_3parms)
-2*logLik(spautor_3parms)
-2*logLik(spautor_3parms) + 2*5 + 2*3

# try a geostatistical model based on centroids
splm_ml = splm(Estimate ~ stockname, data = seals_sf, 
	spcov_type = "circular",
	xcoord = polycentroids$X, ycoord = polycentroids$Y,
	estmethod = 'ml')
summary(splm_ml)
-2*logLik(splm_ml)
-2*logLik(splm_ml) + 2*5 + 2*3

# try the Tieseldorf weights
Tdorfsum = apply(Nmat4,1,sum)
W_Tdorf = (1/sqrt(Tdorfsum))*Nmat4
K_Tdorf_vec = 1/sqrt(Tdorfsum)
1/max(eigen(W_Tdorf)$values)
1/min(eigen(W_Tdorf)$values)

spautor_Tdorf = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'ml', spcov_type = 'car',
	W = W_Tdorf, M = K_Tdorf_vec, row_st = FALSE)
summary(spautor_Tdorf)
-2*logLik(spautor_Tdorf)
-2*logLik(spautor_Tdorf) + 2*5 + 2*2

################################################################################
#-------------------------------------------------------------------------------
#          Confidence Intervals and ANOVA on Fixed Effects
#-------------------------------------------------------------------------------
################################################################################

smry_spautor_3parms = summary(spautor_3parms)
smry_spautor_3parms
smry_spautor_3parms$coefficients$fixed

# confidence interval for Dixon/Cape Decision
c(smry_spautor_3parms$coefficients$fixed[2,1] - 
		1.96*smry_spautor_3parms$coefficients$fixed[2,2],
	smry_spautor_3parms$coefficients$fixed[2,1] + 
		1.96*smry_spautor_3parms$coefficients$fixed[2,2])
		
# Estimate the difference between the mean of the southern
# three stocks minus the mean of the northern two stocks
 
ell = c(1/3,1/3,-1/2,-1/2,1/3)
DCest = t(ell) %*% summary(spautor_3parms)$coefficient$fixed$estimates
DCest
#standard error
DCse = sqrt(t(ell) %*% vcov(spautor_3parms) %*% ell)
DCse
# confidence interval
c(DCest - qnorm(.975)*DCse, DCest + qnorm(.975)*DCse)


# set-first-treatment-to-zero parameterization (default for R)
spautor_3parms_set1 = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'ml', spcov_type = 'car', row_st = TRUE,
	control = list(reltol = 1e-8))
anova(spautor_3parms_set1)

# check it manually
# intercept is contrasted to all other treatments, so test for only intercept
L = matrix(c(1,0,0,0,0), nrow = 1)
L
Chi2 = t(L %*% coef(spautor_3parms)) %*% 
	solve(L %*% vcov(spautor_3parms) %*% t(L)) %*% 
	(L %*% coef(spautor_3parms))
Chi2
1 - pchisq(Chi2,df = 1)

#now test for stock, where each level is already contrasted with intercept
L = cbind(rep(0, times = 4),  diag(4))
L
Chi2 = t(L %*% coef(spautor_3parms_set1)) %*% 
	solve(L %*% vcov(spautor_3parms_set1) %*% t(L)) %*% 
	(L %*% coef(spautor_3parms_set1))
Chi2
1 - pchisq(Chi2,df = 4)

# ------------------------------------------------------------------------------
#        Use F-test with stock effects, no intercept
# ------------------------------------------------------------------------------

spautor_4Ftest = spautor(Estimate ~ -1 + stockname, data = seals_sf, 
	estmethod = 'ml', spcov_type = 'car', row_st = TRUE)
summary(spautor_4Ftest)

#stock
L = cbind(diag(4), rep(0, times = 4)) - cbind(rep(0, times = 4),  diag(4))
L
# see equation (5.6)
Fval = t(L %*% coef(spautor_4Ftest)) %*% 
	solve(L %*% vcov(spautor_4Ftest) %*% t(L)) %*% 
	(L %*% coef(spautor_4Ftest))/
	(4*max(coef(spautor_4Ftest, type = 'spcov')[c('de','extra')]))
# note that using max(coef(spautor_4Ftest, type = 'spcov')[c('de','extra')])
# is equivalent to factoring out an overall variance parameter, and then
# leaving the ratio, and 1-ration of the two variance parameters 
# for multiplying R and I.
1 - pf(Fval, df1 = sum(!is.na(seals_sf$Estimate)) - 5, df2 = 4)

# ------------------------------------------------------------------------------
#        likelihood ratio test
# ------------------------------------------------------------------------------

# likelihood ratio test
spautor_3parms_meanonly = spautor(Estimate ~ 1, data = seals_sf, 
	estmethod = 'ml', spcov_type = 'car', row_st = TRUE)
summary(spautor_3parms_meanonly)
-2*logLik(spautor_3parms_meanonly)
anova(spautor_3parms_meanonly, spautor_3parms)
# check it manually
chi2 = 2*logLik(spautor_3parms) - 2*logLik(spautor_3parms_meanonly)
chi2
1-pchisq(chi2, df = 4)


################################################################################
#-------------------------------------------------------------------------------
#                  Fixed Effects Table
#-------------------------------------------------------------------------------
################################################################################


summary(spautor_3parms)
library(xtable)
FE_table = summary(spautor_3parms)$coefficient$fixed
FE_table
print(
    xtable(FE_table, 
      align = c('l',rep('l', times = length(FE_table[1,]))),
      digits = c(0,4,4,3,5),
      caption = 'Fitted fixed effects',
      label = 'tab:SealsFixEff'
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = TRUE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
