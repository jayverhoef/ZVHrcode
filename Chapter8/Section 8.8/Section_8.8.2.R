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

file_name = 'Caribou_proflike'
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

file_name = 'Caribou_compmod'
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


est_se = cbind(summary(lm_all)$coefficients$fixed[,1],
	summary(sar_1rst)$coefficients$fixed[,1],
	summary(car_1rst)$coefficients$fixed[,1],
	summary(geo_exp)$coefficients$fixed[,1],
	summary(geo_exp5)$coefficients$fixed[,1],
	summary(lm_all)$coefficients$fixed[,2],
	summary(sar_1rst)$coefficients$fixed[,2],
	summary(car_1rst)$coefficients$fixed[,2],
	summary(geo_exp)$coefficients$fixed[,2],
	summary(geo_exp5)$coefficients$fixed[,2])

print(
    xtable(est_se, 
      align = c('l',rep('l', times = length(est_se[1,]))),
      digits = c(0,rep(3, times = 5),rep(3, times = 5)),
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

sar_1rst_alt = spautor(z ~ -1 + trt, data = caribouDF, W = Nmat1, 
	spcov_type = 'sar', estmethod = 'reml', row_st = TRUE)
summary(sar_1rst_alt)
summary(sar_1rst)

covb_sar1rstalt = vcov(sar_1rst_alt)
ell = c(1, 0, -1, 1, 0, -1)
t(ell) %*% coef(sar_1rst_alt)
sqrt(t(ell) %*% covb_sar1rstalt %*% ell)
t(ell) %*% coef(sar_1rst_alt)/
sqrt(t(ell) %*% covb_sar1rstalt %*% ell)

################################################################################
#-------------------------------------------------------------------------------
#                             Simulation
#-------------------------------------------------------------------------------
################################################################################

set.seed(102)
sim1D = sprnorm(spcov_params = coef(geo_exp, type = "spcov"),
	data = data.frame(x = 1:1000, y = rep(1, times = 1000)),
	xcoord = x, ycoord = y)

# randomized treatments
set.seed(2001)
trt_ran = sample(rep(1:6, times = 5))

file_name = 'Caribou_sim1000'
pdf(paste0(file_name,'.pdf'), width = 12, height = 9)

	layout(matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))
	
		padj = -.5
		adj = -.28
		cex_mtext = 3.3
		cex_lab = 2.4
		cex_axis = 1.5

		par(mar = c(7,7,4,2), mgp=c(4, 1.3, 0))
		plot(1:1000, sim1D, pch = 19, cex = .5, xlab = 'coordinate',
			ylab = 'Simulated Value', cex.lab = cex_lab, cex.axis = cex_axis)
		points(321:350, sim1D[321:350], pch = 1, cex = 1.5)
		points(901:930, sim1D[901:930], pch = 1, cex = 1.5)
		mtext('A', adj = -.115, cex = cex_mtext, padj = padj)
		
		plot(321:350, sim1D[321:350], cex = 2, pch = 19, cex.lab = 2.5, cex.axis = 1.5,
			xlab = 'Coordinate', ylab = 'Simulated Value')
		mtext('B', adj = adj, cex = cex_mtext, padj = padj)

		plot(901:930, sim1D[901:930], cex = 2, pch = 19, cex.lab = 2.5, cex.axis = 1.5,
			xlab = 'Coordinate', ylab = 'Simulated Value', ylim = c(.2, 1.65))
		text(901:930, rep(1.3, times = 30), 
			labels = as.character(trt_ran), cex = 1.7)
		text(901:930, rep(1.45, times = 30), 
			labels = as.character(kronecker(1:6,rep(1,times = 5))), cex = 1.7)
	  text(901:930, rep(1.6, times = 30), 
			labels = as.character(rep(1:6, times = 5)), cex = 1.7)
		mtext('C', adj = adj, cex = cex_mtext, padj = padj)

	layout(1)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))


d1 = data.frame(x = 1:30, y = rep(1, times = 30), 
	trt = as.factor(rep(1:6, times = 5)), z = sim1D[901:930])
d2 = data.frame(x = 1:30, y = rep(1, times = 30), 
	trt = as.factor(kronecker(1:6,rep(1,times = 5))), z = sim1D[901:930])
d3 = data.frame(x = 1:30, y = rep(1, times = 30), 
	trt = as.factor(trt_ran), z = sim1D[901:930])
# fit uncorrelated models to all 3 designs
sim_ind1 = splm(z ~ -1 + trt, data = d1, xcoord = x, ycoord = y,
	spcov_type = 'none')
sim_ind2 = splm(z ~ -1 + trt, data = d2, xcoord = x, ycoord = y,
	spcov_type = 'none')
sim_ind3 = splm(z ~ -1 + trt, data = d3, xcoord = x, ycoord = y,
	spcov_type = 'none')
# fit autocorrelated models to all 3 designs with known covariance
sim_exp1 = splm(z ~ -1 + trt, data = d1, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = coef(geo_exp, type = 'spcov')['range'],
		de = coef(geo_exp, type = 'spcov')['de'],
		ie = coef(geo_exp, type = 'spcov')['ie'],
		known = c('range', 'de', 'ie')) )
sim_exp2 = splm(z ~ -1 + trt, data = d2, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = coef(geo_exp, type = 'spcov')['range'],
		de = coef(geo_exp, type = 'spcov')['de'],
		ie = coef(geo_exp, type = 'spcov')['ie'],
		known = c('range', 'de', 'ie')) )
sim_exp3 = splm(z ~ -1 + trt, data = d3, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = coef(geo_exp, type = 'spcov')['range'],
		de = coef(geo_exp, type = 'spcov')['de'],
		ie = coef(geo_exp, type = 'spcov')['ie'],
		known = c('range', 'de', 'ie')) )
covb_ind1 = vcov(sim_ind1)
covb_ind2 = vcov(sim_ind2)
covb_ind3 = vcov(sim_ind3)
covb_exp1 = vcov(sim_exp1)
covb_exp2 = vcov(sim_exp2)
covb_exp3 = vcov(sim_exp3)
ell1 = c(1, -1, 0, 0, 0, 0)
ell2 = c(1, 0, 0, -1, 0, 0)
ell3 = c(1, 0, 0, 0, 0, -1)

sim_cont = rbind(
	cbind(summary(sim_ind3)$coefficients$fixed$estimates,
			summary(sim_exp1)$coefficients$fixed$estimates,
			summary(sim_exp2)$coefficients$fixed$estimates,
			summary(sim_exp3)$coefficients$fixed$estimates,
		summary(sim_ind3)$coefficients$fixed$Std_Error,
			summary(sim_exp1)$coefficients$fixed$Std_Error,
			summary(sim_exp2)$coefficients$fixed$Std_Error,
			summary(sim_exp3)$coefficients$fixed$Std_Error
	),
	cbind(
		# contrasts estimates
		cbind(
			c(ell1 %*% coef(sim_ind3),
				ell2 %*% coef(sim_ind3),
				ell3 %*% coef(sim_ind3)),
			c(ell1 %*% coef(sim_exp1),
				ell2 %*% coef(sim_exp1),
				ell3 %*% coef(sim_exp1)),
			c(ell1 %*% coef(sim_exp2),
				ell2 %*% coef(sim_exp2),
				ell3 %*% coef(sim_exp2)),
			c(ell1 %*% coef(sim_exp3),
				ell2 %*% coef(sim_exp3),
				ell3 %*% coef(sim_exp3))
		),
		# standard errors
		sqrt(cbind(
			c(t(ell1) %*% covb_ind3 %*% ell1,
				t(ell2) %*% covb_ind3 %*% ell2,
				t(ell3) %*% covb_ind3 %*% ell3),
			c(t(ell1) %*% covb_exp1 %*% ell1,
				t(ell2) %*% covb_exp1 %*% ell2,
				t(ell3) %*% covb_exp1 %*% ell3),
			c(t(ell1) %*% covb_exp2 %*% ell1,
				t(ell2) %*% covb_exp2 %*% ell2,
				t(ell3) %*% covb_exp2 %*% ell3),
			c(t(ell1) %*% covb_exp3 %*% ell1,
				t(ell2) %*% covb_exp3 %*% ell2,
				t(ell3) %*% covb_exp3 %*% ell3)
		))
	)
)

print(
    xtable(sim_cont, 
      align = c('l',rep('l', times = length(sim_cont[1,]))),
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

mean(coef(sim_ind3)^2)
mean(coef(sim_exp3)^2)

i  = 1
ind_mspe_true0 = rep(NA, times = 971)
exp_mspe_true0 = rep(NA, times = 971)
ind_T1err_true0 = rep(NA, times = 971)
exp_T1err_true0 = rep(NA, times = 971)
ind_mspe_rlzdmean = rep(NA, times = 971)
exp_mspe_rlzdmean = rep(NA, times = 971)
ind_T1err_rlzdmean = rep(NA, times = 971)
exp_T1err_rlzdmean = rep(NA, times = 971)
for(i in 1:971) {
	dtemp = data.frame(x = 1:30, y = rep(1, times = 30), 
		trt = as.factor(trt_ran), z = sim1D[i:(i + 29)])
	ind_temp = splm(z ~ -1 + trt, data = dtemp, xcoord = x, ycoord = y,
		spcov_type = 'none')
	exp_temp = splm(z ~ -1 + trt, data = dtemp, xcoord = x, ycoord = y,
		spcov_type = 'exponential')
	ind_mspe_true0[i] = mean(coef(ind_temp)^2)
	exp_mspe_true0[i] = mean(coef(exp_temp)^2)
	ind_T1err_true0[i] = mean(abs(coef(ind_temp))/
		summary(ind_temp)$coefficients$fixed$Std_Error > 2.06)
	exp_T1err_true0[i] =  mean(abs(coef(exp_temp))/
		summary(exp_temp)$coefficients$fixed$Std_Error > 2.06)
	rlzdmean = mean(sim1D[i:(i + 29)])
	ind_mspe_rlzdmean[i] = mean((coef(ind_temp) - rlzdmean)^2)
	exp_mspe_rlzdmean[i] = mean((coef(exp_temp) - rlzdmean)^2)
	ind_T1err_rlzdmean[i] = mean(abs(coef(ind_temp) - rlzdmean)/
		summary(ind_temp)$coefficients$fixed$Std_Error > 2.063)
	exp_T1err_rlzdmean[i] =  mean(abs(coef(exp_temp)- rlzdmean)/
		summary(exp_temp)$coefficients$fixed$Std_Error > 2.063)
}

plot(ind_mspe_true0, type = 'l', lty = 2)
lines(exp_mspe_true0)
mean(ind_mspe_true0)
mean(exp_mspe_true0)

mean(ind_T1err_true0)
mean(exp_T1err_true0)

mean(ind_mspe_rlzdmean)
mean(exp_mspe_rlzdmean)

mean(ind_T1err_rlzdmean)
mean(exp_T1err_rlzdmean)

mean(sim1D[901:930])

mean((coef(sim_ind3) - mean(sim1D[901:930]))^2)
mean((coef(sim_exp3) - mean(sim1D[901:930]))^2)

0.452/0.106
0.836/0.467

-0.323/0.149

d4 = data.frame(x = 1:30, y = rep(1, times = 30), 
	trt = as.factor(c(rep(c(1,2), times = 5), rep(c(3,4), times = 5), 
		rep(c(5,6), times = 5))), 
	z = sim1D[901:930])
sim_exp4 = splm(z ~ -1 + trt, data = d4, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = coef(geo_exp, type = 'spcov')['range'],
		de = coef(geo_exp, type = 'spcov')['de'],
		ie = coef(geo_exp, type = 'spcov')['ie'],
		known = c('range', 'de', 'ie')) )
ell4 = c(0, 0, 1, -1, 0, 0)
ell5 = c(0, 0, 0, 0, 1, -1)
covb_exp4 = vcov(sim_exp4)
sqrt(t(ell1) %*% covb_exp1 %*% ell1)
sqrt(t(ell4) %*% covb_exp1 %*% ell4)
sqrt(t(ell5) %*% covb_exp1 %*% ell5)
sqrt(t(ell1) %*% covb_exp4 %*% ell1)
sqrt(t(ell4) %*% covb_exp4 %*% ell4)
sqrt(t(ell5) %*% covb_exp4 %*% ell5)

################################################################################
#-------------------------------------------------------------------------------
#                     Model selection on covariates
#-------------------------------------------------------------------------------
################################################################################

# look at raw means with labels
spfit = spautor(z ~ -1 + trt, data = caribouDF, W = Nmat1, 
	spcov_type = 'sar', estmethod = 'reml', row_st = TRUE)
spfit_sumry = summary(spfit)

file_name = "Caribou_means_interact"
pdf(paste0(file_name,'.pdf'), width = 6, height = 8)

	offset = .05
	par(mar = c(5,5,1,1))
	plot(c(1 - offset, 3 + offset), c(1.75, 2.45), type = 'n',
		 xaxt = 'n', xlab = '', ylab = '% Nitrogen', cex.lab = 2, cex.axis = 1.5)
	axis(1, at = c(1:3), labels = c('shade', 'clear', 'none'), 
		cex.axis = 2)
	trt = 1
	trp = 1
	points(trp - offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 19)
	lines(c(trp - offset,trp - offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	trt = 2
	trp = 2
	points(trp - offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 19)
	lines(c(trp - offset,trp - offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	trt = 3
	trp = 3
	points(trp - offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 19)
	lines(c(trp - offset,trp - offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	trt = 4
	trp = 1
	points(trp + offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 1)
	lines(c(trp + offset,trp + offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	trt = 5
	trp = 2
	points(trp + offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 1)
	lines(c(trp + offset,trp + offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	trt = 6
	trp = 3
	points(trp + offset,
		spfit_sumry$coefficients$fixed$estimates[trt], cex = 3, pch = 1)
	lines(c(trp + offset,trp + offset), 
		c(spfit_sumry$coefficients$fixed$estimates[trt] - 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt], spfit_sumry$coefficients$fixed$estimates[trt] + 1.96*spfit_sumry$coefficients$fixed$Std_Error[trt]), lwd = 5)
	lines(c((1:3) - offset), coef(spfit)[1:3], lty = 2, lwd = 3)
	lines(c((1:3) + offset), coef(spfit)[4:6], lty = 2, lwd = 3)
	legend(2,2.45, legend = c('Watered','Control'), pch = c(19,1), cex = 2)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))





lm_mains = lm(z ~ water + tarp, data = caribouDF)

car_mains = spautor(z ~ water + tarp, data = caribouDF, W = Nmat1, 'sar',
  estmethod = 'ml')

geo_mains = splm(z ~ water + tarp, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical',
	control = list(reltol = 1e-6), estmethod = 'ml')

lm_tarp = lm(z ~ tarp, data = caribouDF)

car_tarp = spautor(z ~ tarp, data = caribouDF, W = Nmat1, 'car',
  estmethod = 'ml')

geo_tarp = splm(z ~ tarp, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical',
	control = list(reltol = 1e-6), estmethod = 'ml')

lm_mean = lm(z ~ 1, data = caribouDF)

car_mean = spautor(z ~ 1, data = caribouDF, W = Nmat1, 'car',
  estmethod = 'ml')

geo_mean = splm(z ~ 1, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical',
	control = list(reltol = 1e-6), estmethod = 'ml')

#-------------------------------------------------------------------------------
#                    Testing anova() function
#-------------------------------------------------------------------------------

#behavior for single objects
anova(lm_all)
anova(car_tarp)
summary(car_tarp)

#behavior for 2 objects
anova(lm_all, lm_mains)

anova(car_all, car_mains)
car_all$optim$value
car_mains$optim$value

anova(geo_all, geo_mains)
geo_all$optim$value
geo_mains$optim$value

#behavior for multiple objects
anova(lm_all, lm_mains, lm_tarp, lm_mean)
#effect of order
anova(lm_all, lm_mean, lm_tarp, lm_mains)

anova(car_all, car_mains, car_tarp)
anova(geo_all, geo_mains, geo_tarp)

#effect of defaults
geo_all_def = splm(z ~ water*tarp, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical')
geo_mains_def = splm(z ~ water + tarp, data = caribouDF,
	xcoord = 'x', ycoord = 'y', spcov_type = 'spherical')
anova(geo_all_def, geo_mains_def)
geo_all_def$optim$value
geo_mains_def$optim$value

# contrasts for geo model
colnames(model.matrix(geo_tarp))
tidy(anova(geo_tarp, L = list(c(0, 1, -1))))
#by hand
L = matrix(c(0, 1, -1),ncol = 1)
t(L) %*% geo_tarp$coefficients$fixed
sqrt(t(L) %*% geo_tarp$vcov$fixed %*% L)
tval = abs(as.numeric(t(L) %*% geo_tarp$coefficients$fixed/
	sqrt(t(L) %*% geo_tarp$vcov$fixed %*% L)))
tval
Fval = tval^2
Fval
# use t-distribution
2*(1 - pt(2.80963, 24))
# use standard normal
2*(1 - pnorm(2.80963))


# contrasts for car model
colnames(model.matrix(car_tarp))
tidy(anova(car_tarp, L = list(c(0, 1, -1))))
#by hand
L = matrix(c(0, 1, -1),ncol = 1)
t(L) %*% car_tarp$coefficients$fixed
sqrt(t(L) %*% car_tarp$vcov$fixed %*% L)
tval = abs(as.numeric(t(L) %*% car_tarp$coefficients$fixed/
	sqrt(t(L) %*% car_tarp$vcov$fixed %*% L)))
tval
Fval = tval^2
Fval
# use t-distribution
2*(1 - pt(2.80963, 24))
# use standard normal
2*(1 - pnorm(2.80963))


#-------------------------------------------------------------------------------
#
#           makeCovMat
#
#-------------------------------------------------------------------------------

#' make a CAR/SAR covariance matrix for modeling
#'
#' make a CAR/SAR covariance matrix for modeling
#'
#' @param theta covariance parameters, with overall variance parameter profiled out.
#' @param indComp an additive independent component to the model.  Default is TRUE.  
#' @param Nmat neighborhood matrix
#' @param distMat distance matrix
#' @param indSamp indicator vector for wheter location was sampled. Zero, or FALSE, indicates missing value
#' @param model either 'CAR' or 'SAR' to be the model associated with the Nmat argument
#' @param logical value on whether Nmat should be row-standardized
#' @param rhobound a vector of two elements containing bounds for rho.  This should be determined from the eigenvalues of Nmat prior to running function
#'
#' @return two times the negative log-likelihood
#'
#' @author Jay Ver Hoef jverhoef
#' @rdname m2LL
#' @export m2LL 

makeCovMat = function(theta, nN, indComp = TRUE, Nmat = NULL, 
  model = 'CAR', rowStand = TRUE, rhoBound = c(-1,1))
{
  nN = dim(distMat)[1]
  V = matrix(1, nrow = nN, ncol = nN )
  diag(V) = 0
  itheta = 0
  if(!is.null(Nmat)) {
    V = as(Nmat, 'sparseMatrix')
  }
  if(!is.null(model)) {
    itheta = itheta + 1
    rho = rhoBound[1] + .00005 + exp(theta[itheta])/
      (1 + exp(theta[itheta]))*.9999*(rhoBound[2] - rhoBound[1])
    attr(theta,'names')[itheta] = 'logitRho'
    rs = rep(1, times = nN)
    if(rowStand) rs = apply(V,1,sum)
    if(model == 'CAR')  V = diag(rs) - rho*V
    if(model == 'SAR') V = (diag(nN) - rho*(1/rs)*V) %*%
      (diag(nN) - rho*t((1/rs)*V))
  }
  if(indComp & is.null(Nmat)) {
    itheta = itheta + 1
    relEps = exp(theta[itheta])
    attr(theta,'names')[itheta] = 'relEps'
    V = V -  V %*% solve(V + diag(rep(relEps,times = nN)),V)
  }
  if(indComp & is.null(Nmat)) {
    V = diag(nN)
  }
  V
}

#-------------------------------------------------------------------------------
#
#           m2LL
#
#-------------------------------------------------------------------------------

#' two times the negative log-likelihood
#'
#' two times the negative log-likelihood
#'
#' @param theta covariance parameters, with overall variance parameter profiled out.
#' @param X design matrix for fixed effects
#' @param y vector of data for response variable
#' @param indComp an additive independent component to the model.  Default is TRUE.  
#' @param Nmat neighborhood matrix
#' @param distMat distance matrix
#' @param indSamp indicator vector for wheter location was sampled. Zero, or FALSE, indicates missing value
#' @param model either 'CAR' or 'SAR' to be the model associated with the Nmat argument
#' @param logical value on whether Nmat should be row-standardized
#' @param rhobound a vector of two elements containing bounds for rho.  This should be determined from the eigenvalues of Nmat prior to running function
#'
#' @return two times the negative log-likelihood
#'
#' @author Jay Ver Hoef jverhoef
#' @rdname m2LL
#' @export m2LL 

m2LL = function(theta, X, y, indComp = TRUE, Nmat = NULL, 
  model = 'CAR', rowStand = TRUE, rhoBound = c(-1,1), MLmeth = 'REMLE')
{
  if(any(abs(theta) > 10.1)) return(1e+32)

  ntheta = 0
  nN = length(y)
	Vi.oo = makeCovMat(theta = theta, nN = nN, indComp = indComp, Nmat = Nmat, 
		model = model, rowStand = rowStand, rhoBound = rhoBound)
  XVi = t(X) %*% Vi.oo
  covbi = XVi %*% X
  covb = solve(covbi)
  bHat = covb %*% XVi %*% y
  r = y - X %*% bHat
  n = length(y)
  p = length(X[1,])
  if(MLmeth == 'MLE') {
  m2LL = n*log(t(r) %*% Vi.oo %*% r) - 
    as.numeric(determinant(Vi.oo, logarithm = TRUE)$modulus) +
    n*(log(2*pi) + 1 - log(n))
	} else if(MLmeth == 'REMLE') {
  m2LL = (n-p)*log(t(r) %*% Vi.oo %*% r) - 
    as.numeric(determinant(Vi.oo, logarithm = TRUE)$modulus) +
    as.numeric(determinant(XVi %*% X, logarithm = TRUE)$modulus) +
    (n - p)*(log(2*pi) + 1 - log((n - p)))
	} else {return('MLmeth argument must be either MLE or REMLE')}
 
	attr(m2LL,'covParms') = theta
	as.numeric(m2LL)
}


X0 = as.matrix(model.matrix(z ~ 1, data = caribouDF))
X1 = as.matrix(model.matrix(z ~ water + tarp + water:tarp, data = caribouDF))
X2 = as.matrix(model.matrix(z ~ water + tarp, data = caribouDF))
X3 = as.matrix(model.matrix(z ~ water + tarp, data = caribouDF))
y = caribouDF$z

ntheta = 2

if(ntheta == 1) {
	# undebug(m2LL)
	# undebug(makeCovMat)
  optOut = optimize(m2LL, interval = c(-10,10), X = X1, y = y, 
		Nmat = Nmat4, indComp = FALSE)
  theta = optOut$minimum
  m2LLargmin = optOut$objective
  } else {
	# undebug(m2LL)
	# undebug(makeCovMat)
  optOut = optim(rep(0, times = ntheta), m2LL, X = X1, y = y, 
		Nmat = Nmat1, indComp = TRUE, model = 'SAR')
  theta = optOut$par
  m2LLargmin = optOut$value
}
X = X1
Vi.oo = makeCovMat(theta, nN = length(y), Nmat = Nmat1, indComp = TRUE,
	model = 'SAR')
XVi = t(X) %*% Vi.oo
covbi = XVi %*% X
covb = solve(covbi)
bHat = covb %*% XVi %*% y
bHat
r = as.matrix(y - X %*% bHat)
n = length(y)
p = length(X[1,])
sigma = as.numeric((t(r) %*% Vi.oo %*% r)/n)
bHat_se = sqrt(sigma*diag(covb))
bHat_se
bHat/bHat_se

optOut = optim(rep(0, times = 3), m2LL, X = X1, y = y, 
	indSamp = indSamp, Nmat = Nmat4, distMat = distMat4, indComp = FALSE,
	MLmeth = 'MLE')
theta = optOut$par
m2LLargmin = optOut$value

X = X1
WMi = makeCovMat(theta, indSamp = indSamp, Nmat = Nmat4, distMat = distMat4,
	indComp = FALSE)
WMi.oo = WMi[indSamp,indSamp] 
WMi.uu = WMi[!indSamp,!indSamp]
WMi.uo = WMi[!indSamp,indSamp]
WMi.ou = WMi[indSamp,!indSamp]
Vi.oo = WMi.oo - WMi.ou %*% solve(WMi.uu, WMi.uo)
XVi = t(X) %*% Vi.oo
covbi = XVi %*% X
covb = solve(covbi)
bHat = covb %*% XVi %*% y
bHat
r = as.matrix(y - X %*% bHat)
n = length(y)
p = length(X[1,])
sigma = as.numeric((t(r) %*% Vi.oo %*% r)/n)
bHat_se = sqrt(sigma*diag(covb))
bHat_se
