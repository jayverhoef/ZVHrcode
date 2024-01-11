sec_path = 'Rcode/Chapter11/Section 11.6/'
setwd(paste0(SLEDbook_path,sec_path))

# attach data library
library(ZVHdata)
library(sp)
library(viridis)
library(classInt)
library(colorspace)
library(spdep)
library(spmodel)
library(xtable)

################################################################################
#-------------------------------------------------------------------------------
#                     Caribou Example
#-------------------------------------------------------------------------------
################################################################################

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

# Table 11.5

# fit 3 different models
spfitIND = splm(z ~ water*tarp, data = caribouDF, xcoord = 'x', ycoord = 'y', 
	spcov_type = 'none', estmethod = 'reml')
smryIND = summary(spfitIND)
spfitSAR = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	spcov_type = 'sar', estmethod = 'reml', row_st = TRUE)
smrySAR = summary(spfitSAR)
spfitEXP = splm(z ~ water*tarp, data = caribouDF, xcoord = 'x', ycoord = 'y', 
	spcov_type = 'exponential', estmethod = 'reml')
smryEXP = summary(spfitEXP)

# compute BLUPTE
ell = c(1, 0, 0, 0, 0, 0)
drvlmSARest = (ell %*% coef(spfitSAR) + mean(resid(spfitSAR)))
drvlmSARse = sqrt(t(ell - apply(model.matrix(z ~ water*tarp, 
	data = caribouDF),2,mean)) %*% vcov(spfitSAR) %*% 
	(ell - apply(model.matrix(z ~ water*tarp, data = caribouDF),2,mean)))

drvlmEXPest = (ell %*% coef(spfitEXP) + mean(resid(spfitEXP)))
drvlmEXPse = sqrt(t(ell - apply(model.matrix(z ~ water*tarp, 
	data = caribouDF),2,mean)) %*% vcov(spfitEXP) %*% 
	(ell - apply(model.matrix(z ~ water*tarp, data = caribouDF),2,mean)))

caribouFE = 
	rbind(
		c(smryIND$coefficients$fixed[1,1], drvlmSARest, drvlmEXPest,
			smryIND$coefficients$fixed[1,2], drvlmSARse, drvlmEXPse),
		cbind(
			smryIND$coefficients$fixed[,1],
			smrySAR$coefficients$fixed[,1],
			smryEXP$coefficients$fixed[,1],
			smryIND$coefficients$fixed[,2],
			smrySAR$coefficients$fixed[,2],
			smryEXP$coefficients$fixed[,2])
	)	

################################################################################
#              Table 11.5
################################################################################

print(
    xtable(caribouFE, 
      align = c('l',rep('l', times = length(caribouFE[1,]))),
      digits = c(0,rep(3, times = 3),rep(3, times = 3)),
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

