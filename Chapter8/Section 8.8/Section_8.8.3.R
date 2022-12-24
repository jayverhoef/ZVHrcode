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
library(spmodel)
library(xtable)


# load data for graphics and analysis
data(caribouDF)

#-------------------------------------------------------------------------------
#                    Create Neighborhood Matrices
#-------------------------------------------------------------------------------

nTot = length(caribouDF)

# get Euclidean distance between centroids of plots
Distmat = as.matrix(dist(caribouDF[,c('x','y')]))
# create first-order neighbor matrix (rook's move) from distances
Nmat1 = (Distmat < 1.01)*1
diag(Nmat1) = 0
Nmat1rst = Nmat1/apply(Nmat1,1,sum)
# create second-order neighbor matrix (neighbors of neighbors) from first-order
Nmat2 = (Nmat1 %*% Nmat1 > 0 | Nmat1 > 0)*1
diag(Nmat2) = 0

################################################################################
#-------------------------------------------------------------------------------
#                     Compare independence model to
#         another independence model with row, column fixed effects
#                exponential, spherical (geostat models) to
#                    CAR and SAR (autoregressive models)
#         for full covariate model (2 main effects and interaction)
#-------------------------------------------------------------------------------
################################################################################

# independence model
lm_all = splm(z ~ water*tarp, data = caribouDF, spcov_type = 'none',
	estmeth = 'ml')
summary(lm_all)
anova(lm_all)
-2*logLik(lm_all)
-2*logLik(lm_all) + 2*6 + 2*1
-2*logLik(lm_all) + 7*log(30)
loocv(lm_all)

# independence model with row and column effects
lm_rcol = splm(z ~ water*tarp + I(as.factor(x)) + I(as.factor(y)), 
	data = caribouDF, spcov_type = 'none',
	estmeth = 'ml')
summary(lm_rcol)
anova(lm_rcol)
-2*logLik(lm_rcol)
-2*logLik(lm_rcol) + 2*15 + 2*1
-2*logLik(lm_rcol) + 16*log(30)
loocv(lm_rcol)

# CAR 1st-order binary weighting
car_1bin = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	estmethod = 'ml', row_st = FALSE, spcov_type = 'car')
summary(car_1bin)
anova(car_1bin)
-2*logLik(car_1bin)
-2*logLik(car_1bin) + 2*6 + 2*2
-2*logLik(car_1bin) + 8*log(30)
loocv(car_1bin)

# CAR 1st-order row-standardized weighting
car_1rst = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	spcov_type = 'car', estmethod = 'ml', row_st = TRUE)
summary(car_1rst)
anova(car_1rst)
-2*logLik(car_1rst)
-2*logLik(car_1rst) + 2*6 + 2*2
-2*logLik(car_1rst) + 8*log(30)
loocv(car_1rst)

# CAR 2nd-order binary weighting
car_2bin = spautor(z ~ water*tarp, data = caribouDF, W = Nmat2, 
	estmethod = 'ml', row_st = FALSE, spcov_type = 'car')
summary(car_2bin)
anova(car_2bin)
-2*logLik(car_2bin)
-2*logLik(car_2bin) + 2*6 + 2*2
-2*logLik(car_2bin) + 8*log(30)
loocv(car_2bin)

# CAR 2nd-order row-standardized weighting
car_2rst = spautor(z ~ water*tarp, data = caribouDF, W = Nmat2, 
	spcov_type = 'car', estmethod = 'ml', row_st = TRUE)
summary(car_2rst)
anova(car_2rst)
-2*logLik(car_2rst)
-2*logLik(car_2rst) + 2*6 + 2*2
-2*logLik(car_2rst) + 8*log(30)
loocv(car_2rst)

# SAR 1st-order binary weighting
sar_1bin = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	estmethod = 'ml', row_st = FALSE, spcov_type = 'sar')
summary(sar_1bin)
anova(sar_1bin)
-2*logLik(sar_1bin)
-2*logLik(sar_1bin) + 2*6 + 2*2
-2*logLik(sar_1bin) + 8*log(30)
loocv(sar_1bin)

# SAR 1st-order row-standardized weighting
sar_1rst = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	spcov_type = 'sar', estmethod = 'ml', row_st = TRUE)
summary(sar_1rst)
anova(sar_1rst)
-2*logLik(sar_1rst)
-2*logLik(sar_1rst) + 2*6 + 2*2
-2*logLik(sar_1rst) + 8*log(30)
loocv(sar_1rst)

# SAR 2nd-order binary weighting
sar_2bin = spautor(z ~ water*tarp, data = caribouDF, W = Nmat2, 
	estmethod = 'ml', row_st = FALSE, spcov_type = 'sar')
summary(sar_2bin)
anova(sar_2bin)
-2*logLik(sar_2bin)
-2*logLik(sar_2bin) + 2*6 + 2*2
-2*logLik(sar_2bin) + 8*log(30)
loocv(sar_2bin)

# SAR 2nd-order row-standardized weighting
sar_2rst = spautor(z ~ water*tarp, data = caribouDF, W = Nmat2, 
	spcov_type = 'sar', estmethod = 'ml', row_st = TRUE)
summary(sar_2rst)
anova(sar_2rst)
-2*logLik(sar_2rst)
-2*logLik(sar_2rst) + 2*6 + 2*2
-2*logLik(sar_2rst) + 8*log(30)
loocv(sar_2rst)

# spherical geostatistical model 
geo_sph = splm(z ~ water*tarp, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical',
	control = list(reltol = 1e-6), estmeth = 'ml')
summary(geo_sph)
anova(geo_sph)
-2*logLik(geo_sph)
-2*logLik(geo_sph) + 2*6 + 2*3
-2*logLik(geo_sph) + 9*log(30)
loocv(geo_sph)

# exponential geostatistical model 
geo_exp = splm(z ~ water*tarp, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'exponential',
	control = list(reltol = 1e-6), estmethod = 'ml')
summary(geo_exp)
anova(geo_exp)
-2*logLik(geo_exp)
-2*logLik(geo_exp) + 2*6 + 2*3
-2*logLik(geo_exp) + 9*log(30)
loocv(geo_exp)

Model_compare = rbind(
data.frame(model = 'lm_all', m2LL = -2*logLik(lm_all),
	AIC = -2*logLik(lm_all) + 2*6 + 2*1,
	BIC = -2*logLik(lm_all) + 7*log(30),
	LOOCV = loocv(lm_all)),
data.frame(model = 'lm_rcol', m2LL = -2*logLik(lm_rcol),
	AIC = -2*logLik(lm_rcol) + 2*15 + 2*1,
	BIC = -2*logLik(lm_rcol) + 16*log(30),
	LOOCV = loocv(lm_rcol)),
data.frame(model = 'car_1bin', m2LL = -2*logLik(car_1bin),
	AIC = -2*logLik(car_1bin) + 2*6 + 2*2,
	BIC = -2*logLik(car_1bin) + 8*log(30),
	LOOCV = loocv(car_1bin)),
data.frame(model = 'car_1rst', m2LL = -2*logLik(car_1rst),
	AIC = -2*logLik(car_1rst) + 2*6 + 2*2,
	BIC = -2*logLik(car_1rst) + 8*log(30),
	LOOCV = loocv(car_1rst)),
data.frame(model = 'car_2bin', m2LL = -2*logLik(car_2bin),
	AIC = -2*logLik(car_2bin) + 2*6 + 2*2,
	BIC = -2*logLik(car_2bin) + 8*log(30),
	LOOCV = loocv(car_1bin)),
data.frame(model = 'car_2rst', m2LL = -2*logLik(car_2rst),
	AIC = -2*logLik(car_2rst) + 2*6 + 2*2,
	BIC = -2*logLik(car_2rst) + 8*log(30),
	LOOCV = loocv(car_2rst)),
data.frame(model = 'sar_1bin', m2LL = -2*logLik(sar_1bin),
	AIC = -2*logLik(sar_1bin) + 2*6 + 2*2,
	BIC = -2*logLik(sar_1bin) + 8*log(30),
	LOOCV = loocv(sar_1bin)),
data.frame(model = 'sar_1rst', m2LL = -2*logLik(sar_1rst),
	AIC = -2*logLik(sar_1rst) + 2*6 + 2*2,
	BIC = -2*logLik(sar_1rst) + 8*log(30),
	LOOCV = loocv(sar_1rst)),
data.frame(model = 'sar_2bin', m2LL = -2*logLik(sar_2bin),
	AIC = -2*logLik(sar_2bin) + 2*6 + 2*2,
	BIC = -2*logLik(sar_2bin) + 8*log(30),
	LOOCV = loocv(sar_1bin)),
data.frame(model = 'sar_2rst', m2LL = -2*logLik(sar_2rst),
	AIC = -2*logLik(sar_2rst) + 2*6 + 2*2,
	BIC = -2*logLik(sar_2rst) + 8*log(30),
	LOOCV = loocv(sar_2rst)),
data.frame(model = 'geo_sph', m2LL = -2*logLik(geo_sph),
	AIC = -2*logLik(geo_sph) + 2*6 + 2*3,
	BIC = -2*logLik(geo_sph) + 9*log(30),
	LOOCV = loocv(geo_sph)),
data.frame(model = 'geo_exp', m2LL = -2*logLik(geo_exp),
	AIC = -2*logLik(geo_exp) + 2*6 + 2*3,
	BIC = -2*logLik(geo_exp) + 9*log(30),
	LOOCV = loocv(geo_exp))
)
Model_compare

print(
    xtable(Model_compare, 
      align = c('l',rep('l', times = length(Model_compare[1,]))),
      digits = c(0,0,2,2,2,4),
      caption = 'Fitted fixed effects',
      label = 'tab:SealsFixEff'
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

################################################################################
#-------------------------------------------------------------------------------
#          Likelihood Surface and Profile Likelihood
#-------------------------------------------------------------------------------
################################################################################

spfit = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	spcov_type = 'sar', estmethod = 'reml', row_st = TRUE)

summary(spfit)

#set grids sequences and dimensions
seqs = seq(-1,1,by = 0.05)
grid_dim = length(seqs)

# create a sequence of values to try for range and partial sill
de_seq = log(coef(spfit, type = "spcov")['de']) + .9*seqs
de_seq = exp(de_seq)
range_seq = .9*(seqs*.75 + .25)

zA = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		spautorout = spautor(
			z ~ water*tarp, data = caribouDF, estmethod = 'reml',
			control = list(reltol = 1e-7), W = Nmat1, row_st = TRUE,
			spcov_initial = spcov_initial(spcov_type = 'sar', 
				de = de_seq[i], range = range_seq[j],
				known = c('de','range'))
		)	
		zA[i,j] = logLik(spautorout)
	}
}

store_range = rep(NA, times = grid_dim)
for(i in 1:grid_dim) {
		spautorout = spautor(
			z ~ water*tarp, data = caribouDF, estmethod = 'reml',
			control = list(reltol = 1e-7), W = Nmat1, row_st = TRUE,
			spcov_initial = spcov_initial(spcov_type = 'sar', 
				range = range_seq[i],
				known = c('range'))
		)	
	store_range[i] = logLik(spautorout)
}

store_de = rep(NA, times = grid_dim)
for(i in 1:grid_dim) {
		spautorout = spautor(
			z ~ water*tarp, data = caribouDF, estmethod = 'reml',
			control = list(reltol = 1e-7), W = Nmat1, row_st = TRUE,
			spcov_initial = spcov_initial(spcov_type = 'sar', 
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

file_name = 'figures/Caribou_proflike'
pdf(paste0(file_name,'.pdf'), width = 12, height = 12)
#tiff(paste0(file_name,'.tiff'), width = 960, height = 960)

	padj = -.5
	adj = -.17
	cex_mtext = 4
	cex_lab = 2.7
	cex_axis = 2
	
#	layout(matrix(c(1,1,2,2,
#									1,1,3,3), nrow = 2, byrow = TRUE))
layout(matrix(1:2, nrow = 2))
#	nbrks = 20
#	brks = quantile(zA, probs = (0:nbrks)/nbrks)
#	cramp = viridis(nbrks)
#	par(mar = c(7,7,5,2), mgp=c(4, 1.3, 0))
#	image(de_seq, range_seq, zA, breaks = brks, col = cramp,
#		cex.main = 2, xlab = '', ylab = expression(rho),
#		cex.axis = cex_axis, cex.lab = cex_lab)
#	title(xlab = expression(sigma^2), line = 5, cex.lab = cex_lab)
#	points(coef(spfit, type = "spcov")['de'], 
#		coef(spfit, type = "spcov")['range'], 
#		pch = 19, cex = 2.5, col = 'black')
#	mtext('A', adj = adj, cex = cex_mtext, padj = padj)


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
	mtext('A', adj = adj, cex = cex_mtext, padj = padj)


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
	mtext('B', adj = adj, cex = cex_mtext, padj = padj)
		
	layout(1)

dev.off()
	
#system(paste0('tiff2pdf -o','\'',SLEDbook_path,
#	sec_path,file_name,'.pdf','\' ','\'',SLEDbook_path,
#	sec_path,file_name,'.tiff','\''))
#system(paste0('rm ','\'',SLEDbook_path,
#		sec_path,file_name,'.tiff','\''))

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))


################################################################################
#-------------------------------------------------------------------------------
#            Spatial weights compared to Geostatistical Models
#-------------------------------------------------------------------------------
################################################################################

sar_1rst = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	spcov_type = 'sar', estmethod = 'reml', row_st = TRUE)
summary(sar_1rst)
Sigma_sar1rst = coef(sar_1rst, type = 'spcov')['de']*solve(diag(30) - 
	coef(sar_1rst, type = 'spcov')['range']*Nmat1rst) %*% solve(diag(30) - 
	coef(sar_1rst, type = 'spcov')['range']*t(Nmat1rst))
var_sar1rst = diag(Sigma_sar1rst)
Cor_sar1rst = t(1/sqrt(var_sar1rst) * t(1/sqrt(var_sar1rst) * Sigma_sar1rst))

car_1rst = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	spcov_type = 'car', estmethod = 'reml', row_st = TRUE)
summary(car_1rst)
Sigma_car1rst = coef(car_1rst, type = 'spcov')['de']*solve(diag(30) - 
	coef(car_1rst, type = 'spcov')['range']*Nmat1rst) %*% solve(diag(30) - 
	coef(car_1rst, type = 'spcov')['range']*t(Nmat1rst))
var_car1rst = diag(Sigma_car1rst)
Cor_car1rst = t(1/sqrt(var_car1rst) * t(1/sqrt(var_car1rst) * Sigma_car1rst))

geo_exp = splm(z ~ water*tarp, data = caribouDF,  
	xcoord = 'x', ycoord = 'y', 
	spcov_type = 'exponential', estmethod = 'reml')
summary(geo_exp)
Sigma_geoexp = coef(geo_exp, type = 'spcov')['de']*
	exp(-Distmat/coef(geo_exp, type = 'spcov')['range'] ) + 
	coef(geo_exp, type = 'spcov')['ie']*diag(30)
var_geoexp = diag(Sigma_geoexp)
Cor_geoexp = t(1/sqrt(var_geoexp) * t(1/sqrt(var_geoexp) * Sigma_geoexp))
	
geo_exp5 = splm(z ~ water*tarp, data = caribouDF,  
	xcoord = 'x', ycoord = 'y', 
	estmethod = 'reml',
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = 20, known = 'range'))
summary(geo_exp5)
Sigma_geoexp5 = coef(geo_exp5, type = 'spcov')['de']*
	exp(-Distmat/coef(geo_exp5, type = 'spcov')['range'] ) + 
	coef(geo_exp5, type = 'spcov')['ie']*diag(30)
var_geoexp5 = diag(Sigma_geoexp5)
Cor_geoexp5 = t(1/sqrt(var_geoexp5) * t(1/sqrt(var_geoexp5) * Sigma_geoexp5))

#set grids sequences and dimensions
seqs = seq(-1,1,by = 0.05)
grid_dim = length(seqs)

# A.   SAR

# create a sequence of values to try for range and partial sill
de_seqSAR = log(coef(sar_1rst, type = "spcov")['de']) + .9*seqs
de_seqSAR = exp(de_seqSAR)
range_seqSAR = .9*(seqs*.75 + .25)

zA = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		spautorout = spautor(
			z ~ water*tarp, data = caribouDF, estmethod = 'reml',
			W = Nmat1, row_st = TRUE,
			spcov_initial = spcov_initial(spcov_type = 'sar', 
				de = de_seqSAR[i], range = range_seqSAR[j],
				known = c('de','range'))
		)	
		zA[i,j] = logLik(spautorout)
	}
}


# B.   CAR

# create a sequence of values to try for range and partial sill
de_seqCAR = log(coef(car_1rst, type = "spcov")['de']) + .9*seqs
de_seqCAR = exp(de_seqCAR)
range_seqCAR = .9*(seqs*.75 + .25)

zB = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		spautorout = spautor(
			z ~ water*tarp, data = caribouDF, estmethod = 'reml',
			W = Nmat1, row_st = TRUE,
			spcov_initial = spcov_initial(spcov_type = 'car', 
				de = de_seqCAR[i], range = range_seqCAR[j],
				known = c('de','range'))
		)	
		zB[i,j] = logLik(spautorout)
	}
}


# C.   Exponential

de_seqEXP = log(coef(geo_exp, type = "spcov")['de']) + 1*seqs
de_seqEXP = exp(de_seqEXP)
range_seqEXP = coef(geo_exp, type = "spcov")['range'] + 49*seqs
# compute minus 2 times loglikelihood on a grid around REMLE estimates
# empty matrix to store results
zC = matrix(NA, nrow = grid_dim, ncol = grid_dim)
for(i in 1:grid_dim) {
	for(j in 1:grid_dim) {
		splmout = splm(
			z ~ water*tarp, data = caribouDF, 
				xcoord = 'x', ycoord = 'y', estmethod = 'reml',
				spcov_initial = spcov_initial(spcov_type = 'exponential', 
				de = de_seqEXP[i], range = range_seqEXP[j], 
				ie = coef(geo_exp, type = "spcov")['ie'],
				known = c('de','range','ie'))
		)	
		zC[i,j] = logLik(splmout)
	}
}

file_name = 'figures/Caribou_compmod'
#pdf(paste0(file_name,'.pdf'), width = 12, height = 6)
tiff(paste0(file_name,'.tiff'), width = 960, height = 960)

	layout(matrix(1:4, nrow = 2, byrow = TRUE))
		padj = -.5
		adj = -.17
		cex_mtext = 3.3
		cex_lab = 3.5
		cex_axis = 2

		nbrks = 20
		brks = quantile(zA, probs = (0:nbrks)/nbrks)
		cramp = viridis(nbrks)
		par(mar = c(7,7,5,2), mgp=c(4, 1.3, 0))
		image(de_seqSAR, range_seqSAR, zA, breaks = brks, col = cramp,
			cex.main = 2, xlab = '', ylab = expression(rho),
			cex.axis = cex_axis, cex.lab = cex_lab)
		title(xlab = expression(sigma^2), line = 5, cex.lab = cex_lab)
		points(coef(sar_1rst, type = "spcov")['de'], 
			coef(sar_1rst, type = "spcov")['range'], 
			pch = 1, cex = 2.5, col = 'black', lwd = 4)
		mtext('A', adj = adj, cex = cex_mtext, padj = padj)

		nbrks = 20
		brks = quantile(zB, probs = (0:nbrks)/nbrks)
		cramp = viridis(nbrks)
		par(mar = c(7,7,5,2), mgp=c(4, 1.3, 0))
		image(de_seqCAR, range_seqCAR, zB, breaks = brks, col = cramp,
			cex.main = 2, xlab = '', ylab = expression(rho),
			cex.axis = cex_axis, cex.lab = cex_lab)
		title(xlab = expression(sigma^2), line = 5, cex.lab = cex_lab)
		points(coef(car_1rst, type = "spcov")['de'], 
			coef(car_1rst, type = "spcov")['range'], 
			pch = 3, cex = 2.5, col = 'black', lwd = 4)
		mtext('B', adj = adj, cex = cex_mtext, padj = padj)

		nbrks = 20
		brks = quantile(zC, probs = (0:nbrks)/nbrks)
		cramp = viridis(nbrks)
		par(mar = c(7,7,5,2), mgp=c(4, 1.3, 0))
		image(de_seqEXP, range_seqEXP, zC, breaks = brks, col = cramp,
			cex.main = 2, xlab = '', ylab = expression(rho),
			cex.axis = cex_axis, cex.lab = cex_lab)
		title(xlab = expression(sigma^2), line = 5, cex.lab = cex_lab)
		points(coef(geo_exp, type = "spcov")['de'], 
			coef(geo_exp, type = "spcov")['range'], 
			pch = 19, cex = 2.5, col = 'black')
		points(coef(geo_exp5, type = "spcov")['de'], 
			coef(geo_exp5, type = "spcov")['range'], 
			pch = 15, cex = 2.5, col = 'black', lwd = 3)
		mtext('C', adj = adj, cex = cex_mtext, padj = padj)

		cex_all = 2
		dist_cor = data.frame(dist = Distmat[upper.tri(Distmat)],
			cor = Cor_sar1rst[upper.tri(Cor_sar1rst)])
		par(mar = c(7,7,5,2), mgp=c(4, 1.3, 0))
		plot(dist_cor, pch = 1, cex = cex_all, ylim = c(0,1),
			xlab = 'Distance', ylab = 'Autocorrelation',
			cex.lab = cex_lab, cex.axis = cex_axis)
		dist_cor = data.frame(dist = Distmat[upper.tri(Distmat)],
			cor = Cor_car1rst[upper.tri(Cor_sar1rst)])
		points(dist_cor, pch = 3, cex = cex_all)
		dist_cor = data.frame(dist = Distmat[upper.tri(Distmat)],
			cor = Cor_geoexp[upper.tri(Cor_geoexp)])
		points(dist_cor, pch = 19, cex = cex_all)
		dist_cor = data.frame(dist = Distmat[upper.tri(Distmat)],
			cor = Cor_geoexp5[upper.tri(Cor_geoexp5)])
		points(dist_cor, pch = 15, cex = cex_all)
		mtext('D', adj = adj, cex = cex_mtext, padj = padj)

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


