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


# load data for graphics and analysis
data(caribouDF)

#-------------------------------------------------------------------------------
#                    Create Neighborhood Matrices
#-------------------------------------------------------------------------------

nTot = dim(caribouDF)[1]

# get Euclidean distance between centroids of plots
Distmat = as.matrix(dist(caribouDF[,c('x','y')]))
# create first-order neighbor matrix (rook's move) from distances
Nmat1 = (Distmat < 1.01)*1
diag(Nmat1) = 0
Nmat2 = (Nmat1 %*% Nmat1 > 0 | Nmat1 > 0)*1
diag(Nmat2) = 0
M1 = apply(Nmat1, 1, sum)
M1i = 1/apply(Nmat1, 1, sum)
M11 = rep(1, times = length(M1))

# do row-standardization manually
Nmat1st = Nmat1*M1i

#-------------------------------------------------------------------------------
#                    Testing spautor CAR models
#-------------------------------------------------------------------------------

lm_all = lm(z ~ water*tarp, data = caribouDF)
summary(lm_all)

# simple weighting, M all 1's           ***OK***
car_all = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, M = M11, 
	spcov_type = 'car', estmethod = 'ml', row_st = FALSE)
summary(car_all)
  
# simple weighting, leave M unspecified       ***OK, same as above***
car_all_1 = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	spcov_type = 'car', estmethod = 'ml', row_st = FALSE)
summary(car_all_1)
  
# manual row-standardization, M = inverse sum of rows before standardizing          
car_all_RSm = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1st, M = M1i, 
	spcov_type = 'car', estmethod = 'ml', row_st = FALSE)
summary(car_all_RSm)

# automatic row-standardization, M unspecified          
car_all_RSa = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1,
	spcov_type = 'car', estmethod = 'ml', row_st = TRUE)
summary(car_all_RSa)

# simple weighting, leave M unspecified       ***OK, same as above***
sar_all_1 = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	spcov_type = 'sar', estmethod = 'ml', row_st = FALSE)
summary(sar_all_1)

# manual row-standardization, leave M unspecified     
sar_all_RSm = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1st, 
	spcov_type = 'sar', estmethod = 'ml', row_st = FALSE)
summary(sar_all_RSm)

# automatic row-standardization, leave M unspecified  ***OK, same as above***  
sar_all_RSa = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, 
	spcov_type = 'sar', estmethod = 'ml', row_st = TRUE)
summary(sar_all_RSa)

