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

AllPolyCentroids = data.frame(x = coordinates(sealPolys)[,1], 
    y = coordinates(sealPolys)[,2], 
    stockid = as.factor(as.character(sealPolys@data$stockid)),
    polyid = as.factor(as.character(sealPolys@data$polyid)))
polycentroids = st_coordinates(st_centroid(seals_sf$geometry))
rownames(polycentroids) = rownames(seals_sf)
distMat = as.matrix(dist(AllPolyCentroids[,c('x','y')]))/1000
distMat = as.matrix(dist(polycentroids))/1000

# some useful transformations
logit = function(x) {log(x/(1 - x))}
expit = function(x) {exp(x)/(1 + exp(x))}


coords = coordinates(sealPolys)

Nlist2 = apply(Nmat2, 1, function(x) as.integer(which(x > 0)))
  class(Nlist2) = 'nb'
  attr(Nlist2,'type') = 'queen'
  attr(Nlist2,'sym') = TRUE
  attr(Nlist2,'polyid') = attr(Nlist,'polyid')
  attr(Nlist2,'stockid') = attr(Nlist,'polyid')

Nlist4 = apply(Nmat4, 1, function(x) as.integer(which(x > 0)))
  class(Nlist4) = 'nb'
  attr(Nlist4,'type') = 'queen'
  attr(Nlist4,'sym') = TRUE
  attr(Nlist4,'polyid') = attr(Nlist,'polyid')
  attr(Nlist4,'stockid') = attr(Nlist,'polyid')

xleft = 1050000
ybottom = 1030000
xright = xleft + 122000
ytop = ybottom + 62000
xexp = (xright - xleft)*.06
yexp = (ytop - ybottom)*.35

file_name = "Seals_neighbors"
#pdf(paste0(file_name,'.pdf'), width = 9, height = 9)
tiff(paste0(file_name,'.tiff'), width = 720, height = 720)

	layout(matrix(1:4, ncol = 2, byrow = TRUE))
	par(mar = c(.5,.5,.5,.5))

	text_cex = 3.0
	plot(sealPolys)
	rect(xleft - xexp, ybottom - yexp, xright + xexp, ytop + yexp, 
		col = 'grey80', border = NA)
	plot(sealPolys, add = TRUE)
	text(914000, 1242000, 'A', cex = text_cex)

	plot(sealPolys, xlim = c(xleft, xright), ylim = c(ybottom, ytop))
	rect(xleft - xexp, ybottom - 2*yexp, xright + xexp, ytop + 2*yexp, 
		col = 'grey80', border = NA)
	plot(sealPolys, xlim = c(xleft, xright), ylim = c(ybottom, ytop),
		add = TRUE)
	plot(Nlist, coords, add = TRUE, lwd = 2)
	text(1050000, 1120000, 'B', cex = text_cex)

	plot(sealPolys, xlim = c(xleft, xright), ylim = c(ybottom, ytop))
	rect(xleft - xexp, ybottom - 2*yexp, xright + xexp, ytop + 2*yexp, 
		col = 'grey80', border = NA)
	plot(sealPolys, xlim = c(xleft, xright), ylim = c(ybottom, ytop),
		add = TRUE)
	plot(Nlist2, coords, add = TRUE, lwd = 1)
	text(1050000, 1120000, 'C', cex = text_cex)

	plot(sealPolys, xlim = c(xleft, xright), ylim = c(ybottom, ytop))
	rect(xleft - xexp, ybottom - 2*yexp, xright + xexp, ytop + 2*yexp, 
		col = 'grey80', border = NA)
	plot(sealPolys, xlim = c(xleft, xright), ylim = c(ybottom, ytop),
		add = TRUE)
	plot(Nlist4, coords, add = TRUE, lwd = 0.5)
	text(1050000, 1120000, 'D', cex = text_cex)

	layout(1)

dev.off()

system(paste0('tiff2pdf -o','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\' ','\'',SLEDbook_path,
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
#                  CAR and SAR Likelihoods
#-------------------------------------------------------------------------------
################################################################################

#reproduce parts of Figure 5 in Ver Hoef et al. 2018
#investigate CAR versus SAR, and Row-standardized versus unstandardized

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

file_name = "Seals_m2LL_AIC"
pdf(paste0(file_name,'.pdf'), width = 12, height = 6)

	layout(matrix(1:2, nrow = 1))
	labs = c('m0', 'm1','m2', 'm4', 'X0','X1', 'X2', 'X4')
	sar_col = '#4daf4a'
	car_col = '#e41a1c'
	padj = 0
	adj = -.2
	cex_mtext = 3.2
	cex_all = 1.8
	par(mar = c(5,5,4,1))
	plot(store_jitter[,c(1,6)], xlim = c(.9, max(store_jitter[,1]) + .1),
		ylim = c(min(-2*logLik(lmout_m), -2*logLik(lmout_X), store_jitter[,6]),
			max(-2*logLik(lmout_m), -2*logLik(lmout_X), store_jitter[,6]) + 10),
		type = 'n', xlab = '', xaxt = 'n',
		ylab = expression("-2"*italic(L)(bold(theta)~";"~bold(y))), cex.lab = 2, cex.axis = 1.5)
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

file_name = 'Seals_logLik'
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
		ylab = expression("2"*italic(L)[italic("-i,R")](sigma^2~";"~hat(rho),bold(y))), 
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
		ylab = expression("2"*italic(L)[italic("-i,R")](rho~";"~hat(sigma)^2,bold(y))), 
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

# ---------------------- Analytical Approach ----------------------------

# create the covariance matrix and take its inverse
D = diag(apply(Nmat4,1,sum))
DW = D - theta['range']*Nmat4
DWi = solve(DW)
Vi = DW/theta['de']
V = theta['de']*DWi
A = -theta['de']*DWi %*% Nmat4 %*% DWi
		
# compute Fisher Information element by element using trace formula
FI = matrix(rep(NA, times = 4), nrow = 2)
FI[1,1] = sum(diag(Vi %*% DWi %*% Vi %*% DWi))
FI[1,2] = FI[2,1] = sum(diag(Vi %*% DWi %*% Vi %*% A))
FI[2,2] = sum(diag(Vi %*% A %*% Vi %*% A))
FI = 0.5*FI
asycov = solve(FI)
theta[c('de','range')] - 1.96*sqrt(diag(asycov))
theta[c('de','range')] + 1.96*sqrt(diag(asycov))

# asymptotic correlation matrix
diag(1/sqrt(diag(asycov))) %*% asycov %*% diag(1/sqrt(diag(asycov)))

# ---------------------- Numerical Hessian Approach ----------------------------

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

# some useful transformations
logit = function(x) {log(x/(1 - x))}
expit = function(x) {exp(x)/(1 + exp(x))}

# rhohat_rs = rhohat
# sig2hat_rs = sig2hat
# rhohat_un = rhohat
# sig2hat_un = sig2hat

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

file_name = 'Seals_nonstationary'
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
		xlab = 'Distance', ylab = 'Autocorrelation', cex.lab = cex_lab, 
		cex.axis = cex_axis)
	mtext('D', adj = adj, cex = cex_mtext, padj = padj)

	layout(1)

dev.off()
	
system(paste0('img2pdf ','\'',SLEDbook_path,
	sec_path,file_name,'.png','\'', ' -o ', '\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))

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
spautor_3parms = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'ml', spcov_type = 'car', row_st = TRUE)
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
#          ANOVA on Fixed Effects
#-------------------------------------------------------------------------------
################################################################################

smry_spautor_3parms = summary(spautor_3parms)
smry_spautor_3parms
smry_spautor_3parms$coefficients$fixed

# asymptotic marginal test
anova(spautor_3parms)

# check it manually
# intercept
L = matrix(c(1,0,0,0,0), nrow = 1)
L
Chi2 = t(L %*% coef(spautor_3parms)) %*% 
	solve(L %*% vcov(spautor_3parms) %*% t(L)) %*% 
	(L %*% coef(spautor_3parms))
Chi2
1 - pchisq(Chi2,df = 1)
#stock
L = cbind(rep(0, times = 4),  diag(4))
L
Chi2 = t(L %*% coef(spautor_3parms)) %*% 
	solve(L %*% vcov(spautor_3parms) %*% t(L)) %*% 
	(L %*% coef(spautor_3parms))
Chi2
1 - pchisq(Chi2,df = 4)

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
#                  Estimating Fixed Effects
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

ell = c(1,1,0,0,0)
# Dixon/Cape Decision estimate
DCest = t(ell) %*% summary(spautor_3parms)$coefficient$fixed$estimates
DCest
# Dixon/Cape Decision standard error
DCse = sqrt(t(ell) %*% vcov(spautor_3parms) %*% ell)
DCse
# Dixon/Cape Decision confidence interval
DCest - qnorm(.975)*DCse
DCest + qnorm(.975)*DCse

# Estimate the difference between the mean of the southern
# three stocks minus the mean of the northern two stocks 

ell = c(0,1/3,-1/2,-1/2,1/3)
# Dixon/Cape Decision estimate
DCest = t(ell) %*% summary(spautor_3parms)$coefficient$fixed$estimates
DCest
# Dixon/Cape Decision standard error
DCse = sqrt(t(ell) %*% vcov(spautor_3parms) %*% ell)
DCse
# Dixon/Cape Decision confidence interval
DCest - qnorm(.975)*DCse
DCest + qnorm(.975)*DCse



################################################################################
#-------------------------------------------------------------------------------
#                  Prediction
#-------------------------------------------------------------------------------
################################################################################


lmfit_X = splm(Estimate ~ stockname, data = seals_sf, estmethod = 'reml', 
	control = list(reltol = 1e-7), spcov_type = 'none')
summary(lmfit_X)
spfit_m = spautor(Estimate ~ 1, data = seals_sf, 
	estmethod = 'reml', spcov_type = 'car', row_st = TRUE)
summary(spfit_m)
spfit_X = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'reml', spcov_type = 'car', row_st = TRUE)
summary(spfit_m)

seals_pred = seals_sf
seals_pred[,'Estimate_se'] = NA
seals_pred_m = seals_sf
seals_pred_m[,'Estimate_se'] = NA
predout = predict(spfit_X, se.fit = TRUE)
predout_m = predict(spfit_m, se.fit = TRUE)
ind = is.na(seals_pred$Estimate)
seals_pred[ind,'Estimate'] = predout$fit
seals_pred[ind,'Estimate_se'] = predout$se.fit
seals_pred_m[ind,'Estimate'] = predout_m$fit
seals_pred_m[ind,'Estimate_se'] = predout_m$se.fit
loocvout = loocv(spfit_X, cv_predict = TRUE, se.fit = TRUE)
loocvout_m = loocv(spfit_m, cv_predict = TRUE, se.fit = TRUE)
seals_smooth = seals_pred
seals_smooth[!ind,'Estimate'] = loocvout$cv_predict
seals_smooth[!ind,'Estimate_se'] = loocvout$se.fit
seals_smooth_m = seals_pred_m
seals_smooth_m[!ind,'Estimate'] = loocvout_m$cv_predict
seals_smooth_m[!ind,'Estimate_se'] = loocvout_m$se.fit

file_name = 'seal_predhist'
pdf(paste0(file_name,'.pdf'), width = 12, height = 6)

	layout(matrix(1:2, nrow = 1))
	padj = -.5
	adj = -.19
	cex.mtext = 3
	par(mar = c(5,5,4,1))
	hist(seals_sf$Estimate,  breaks = seq(-0.6, 1, by = .04),
		xlim = c(-.42, .5), col = rgb(.8, .8, .8, .3), ylim = c(0,150),
		main = '', cex.lab = 2, cex.axis = 1.5, xlab = 'log(Trend)')
	hist(seals_smooth$Estimate,  breaks = seq(-0.6, 1, by = .02),
		xlim = c(-.5, .5), add = TRUE, col = rgb(.1, .1, .1, .4))
	legend(x = .05, y = 145, legend =c('Raw Values','LOOCV Predictions'),
		pch = 15, col = c(rgb(.7, .7, .7, .5),rgb(.1, .1, .1, .6)),
		cex = 1.2)
	mtext('A', adj = adj, cex = cex.mtext, padj = padj)


	h1out = hist(predict(spfit_m), breaks = seq(-0.15, 0.25, by = .01), plot = FALSE)
	h2out = hist(predict(spfit_X), breaks = seq(-0.15, 0.25, by = .01), plot = FALSE)
	h3out = hist(predict(lmfit_X), breaks = seq(-0.15, 0.25, by = .01), plot = FALSE)

	par(mar = c(5,5,4,1), bg = 'white')
	lwd_bars = 3
	offset = .003
	padj = -.5
	adj = -.15
	cex.mtext = 3
	plot(h1out$mids, h1out$counts, type = 'n', bty = 'n',
		cex.lab = 2, cex.axis = 1.5, ylab = 'Frequency', xlab = 'log(Trend)',
		xaxt = 'n')
	for(i in 1:40) {
		if(h1out$counts[i] > 0)
			lines(rep(h1out$mids[i],2), c(0,h1out$counts[i]), lwd = lwd_bars,
				col = 'gray50', lend = 2)
	}
	for(i in 1:40) {
		if(h2out$counts[i] > 0)
			lines(rep(h2out$mids[i]-offset, 2), c(0,h2out$counts[i]), lwd = lwd_bars,
				col = 'gray80', lend = 2)
	}
	for(i in 1:40) {
		if(h3out$counts[i] > 0)
			lines(rep(h3out$mids[i]+offset, 2), c(0,h3out$counts[i]), lwd = lwd_bars,
				col = 'black', lend = 2)
	}
	lines(c(-.145,.245), c(-.1,-.1))
	axis(1, at = c(-.1,-.05, 0,.05,.1,.15, .2), cex.axis = 1.5)
	d1out = density(predict(spfit_m), bw = .02)
	lines(d1out$x, d1out$y, col = 'gray50', lwd = 4)
	d2out = density(predict(spfit_X), bw = .02)
	lines(d2out$x, d2out$y, col = 'gray80', lwd = 4)
	legend(x = .05, y = 42, legend =c('Indep. Mean','Autocorr. Mean', 'Autocorr. Stock'),
		lwd = 3, col = c('black','gray50','gray80'),
		cex = 1.2)
	mtext('B', adj = adj, cex = cex.mtext, padj = padj)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))




spfit_X_norowst = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'reml', spcov_type = 'car', row_st = FALSE)
spfit_X_Nmat1 = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'reml', W = Nmat1, spcov_type = 'car', row_st = TRUE)
spfit_X_Nmat4 = spautor(Estimate ~ stockname, data = seals_sf, 
	estmethod = 'reml', W = Nmat4, spcov_type = 'car', row_st = TRUE)
loocv(spfit_X)
loocv(spfit_X_norowst)
loocv(spfit_X_Nmat1)
loocv(spfit_X_Nmat4)
loocv(lmfit_X)
loocv(spfit_m)




file_name = 'seal_maps'
	pdf(paste0(file_name,'.pdf'), width = 8, height = 12)

	layout(matrix(1:6, ncol = 2, byrow = TRUE))
	par(mar = c(0,0,0,0), bg = 'gray80')
	brks_pred = quantile(seals_pred$Estimate, probs = seq(0, 1, .125))
	brks_se = quantile(seals_pred[ind,]$Estimate_se, probs = seq(0, 1, .125))
	brks_smoo = quantile(seals_smooth$Estimate, probs = seq(0, 1, .125))
	brks_smoo_se = quantile(seals_smooth$Estimate_se, probs = seq(0, 1, .125))
	brks_smoo_m = quantile(seals_smooth_m$Estimate, probs = seq(0, 1, .125))
	brks_smoo_m_se = quantile(seals_smooth_m$Estimate_se, probs = seq(0, 1, .125))

	# make a color map of predictions

	# A

	cip = classIntervals(seals_sf$Estimate, style = 'fixed', 
		fixedBreaks= brks_pred)
	palp = viridis(8)
	cip_colors = findColours(cip, palp)
	plot(st_geometry(seals_sf), col = cip_colors)
	polycent = st_centroid(seals_pred)
	ind = is.na(seals_sf$Estimate)
	cip = classIntervals(seals_pred[ind,]$Estimate, style = 'fixed', 
		fixedBreaks = brks_pred)
	cip_colors = findColours(cip, palp)
	plot(st_geometry(polycent[ind,]), add = TRUE, pch = 19, cex = 1.3, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_pred, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'A', cex = 4)

	# B

	cip = classIntervals(seals_pred[ind,]$Estimate_se, style = 'fixed', 
		fixedBreaks = brks_se)
	palp = plasma(8)
	cip_colors = findColours(cip, palp)
	plot(st_geometry(seals_sf))
	plot(st_geometry(polycent[ind,]), add = TRUE, pch = 19, cex = 1.3, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_se, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'B', cex = 4)

	# C
	cip = classIntervals(seals_smooth$Estimate, style = 'fixed', 
		fixedBreaks= brks_smoo)
	palp = viridis(8)
	cip_colors = findColours(cip, palp)
	# make a color map of prediction standard errors
	plot(st_geometry(seals_smooth), col = cip_colors)
	polycent = st_centroid(seals_smooth)
	plot(st_geometry(polycent), add = TRUE, pch = 19, cex = 1, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_smoo, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'C', cex = 4)

	# D
	cip = classIntervals(seals_smooth$Estimate_se, style = 'fixed', 
		fixedBreaks= brks_smoo_se)
	palp = plasma(8)
	cip_colors = findColours(cip, palp)
	# make a color map of prediction standard errors
	plot(st_geometry(seals_smooth), col = cip_colors)
	polycent = st_centroid(seals_smooth)
	plot(st_geometry(polycent), add = TRUE, pch = 19, cex = 1, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_smoo_se, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'D', cex = 4)

	# E
	cip = classIntervals(seals_smooth_m$Estimate, style = 'fixed', 
		fixedBreaks= brks_smoo_m)
	palp = viridis(8)
	cip_colors = findColours(cip, palp)
	# make a color map of prediction standard errors
	plot(st_geometry(seals_smooth_m), col = cip_colors)
	polycent = st_centroid(seals_smooth_m)
	plot(st_geometry(polycent), add = TRUE, pch = 19, cex = 1, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_smoo_m, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'E', cex = 4)

	# F
	cip = classIntervals(seals_smooth_m$Estimate_se, style = 'fixed', 
		fixedBreaks= brks_smoo_m_se)
	palp = plasma(8)
	cip_colors = findColours(cip, palp)
	# make a color map of prediction standard errors
	plot(st_geometry(seals_smooth), col = cip_colors)
	polycent = st_centroid(seals_smooth)
	plot(st_geometry(polycent), add = TRUE, pch = 19, cex = 1, col = cip_colors)
	addBreakColorLegend(xleft = 1310000, ybottom = 986649, xright = 1340000, ytop = 1201343,
			breaks = brks_smoo_m_se, colors = palp, cex = 1.5, printFormat = "4.3")
	text(930000, 1180000, 'F', cex = 4)

layout(1)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))


