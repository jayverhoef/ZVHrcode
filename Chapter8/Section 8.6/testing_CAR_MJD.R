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

# do row-standardization manually
M_rs = 1/apply(Nmat1, 1, sum)
Minv_rs = apply(Nmat1, 1, sum)
Nmat1st = Nmat1*M_rs
# precision matrix
diag(Minv_rs) %*% (diag(nTot) - .9*Nmat1st)
# marginal variances for CAR
diag(solve(diag(Minv_rs) %*% (diag(nTot) - .9*Nmat1st)))/
max(diag(solve(diag(Minv_rs) %*% (diag(nTot) - .9*Nmat1st))))
# marginal variances for SAR
diag(solve((diag(nTot) - .9*Nmat1st) %*% (diag(nTot) - .9*Nmat1st)))/
max(diag(solve((diag(nTot) - .9*Nmat1st) %*% (diag(nTot) - .9*Nmat1st))))

# create Tieseldorf weights from binary weights
M_Tie = sqrt(apply(Nmat1,1,sum))
Minv_Tie = 1/M_Tie
Nmat1Tie = Nmat1*M_Tie
# precision matrix
diag(Minv_Tie) %*% (diag(nTot) - .9*Nmat1Tie)
# marginal variances for CAR
diag(solve(diag(Minv_Tie) %*% (diag(nTot) - .9*Nmat1Tie)))/
max(diag(solve(diag(Minv_Tie) %*% (diag(nTot) - .9*Nmat1Tie))))
# marginal variances for SAR
diag(solve((diag(nTot) - .9*Nmat1Tie) %*% (diag(nTot) - .9*Nmat1Tie)))/
max(diag(solve((diag(nTot) - .9*Nmat1Tie) %*% (diag(nTot) - .9*Nmat1Tie))))

# fit using Tieseldorf weights
car_Tie = spautor(z ~ water*tarp, data = caribouDF, W = NmatTie, M = M_Tie,
	spcov_type = 'car', estmethod = 'ml', row_st = FALSE)
summary(car_Tie)

#-------------------------------------------------------------------------------
#                    Testing spautor CAR models
#-------------------------------------------------------------------------------

lm_all = lm(z ~ water*tarp, data = caribouDF)
summary(lm_all)

# simple weighting, M all 1's           ***OK***
car_all = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, M = M11,
	spcov_type = 'car', estmethod = 'ml', row_st = FALSE)
summary(car_all)
# MD: ***OK*** now

# simple weighting, leave M unspecified       ***OK, same as above***
car_all = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1,
	spcov_type = 'car', estmethod = 'ml', row_st = FALSE)
summary(car_all)
# MD: ***OK*** now

# manual row-standardization, M = sum of rows before standardizing
car_all = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1st, M = M1,
	spcov_type = 'car', estmethod = 'ml', row_st = FALSE)
# Error: (I - rho * W) * 1 / M must satisfy the CAR symmetry condition
Covmati = solve(diag(M1i)) %*% (diag(nTot) - 0.9*Nmat1st) 
all(Covmati == t(Covmati))
# MD: This should error given current parameterization of M
Covmati[1,2] ==
Covmati[2,1]

# manual row-standardization, M = inverse sum of rows before standardizing
car_all = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1st, M = M1i,
	spcov_type = 'car', estmethod = 'ml', row_st = FALSE)
# Error: (I - rho * W) * 1 / M must satisfy the CAR symmetry condition
Covmat = (diag(nTot) - 0.9*Nmat1st) %*% diag(M1i)
all(Covmat == t(Covmat))
eigen(Covmat)$values  # this specification should not give error
# MD: ***OK*** now

# manual row-standardization, M = all 1's  (sum of rows after standardizing)
car_all = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1st, M = M11,
	spcov_type = 'car', estmethod = 'ml', row_st = FALSE)
# Error: (I - rho * W) * 1 / M must satisfy the CAR symmetry condition
Covmat = (diag(nTot) - 0.9*Nmat1st) %*% diag(M11)
all(Covmat == t(Covmat))
# MD: This should error given current parameterization of M

# automatic row-standardization, M unspecified
car_all = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1,
	spcov_type = 'car', estmethod = 'ml', row_st = TRUE)
# Error: (I - rho * W) * 1 / M must satisfy the CAR symmetry condition
# MD: ***OK*** now

# automatic row-standardization, M = sum of rows before standardizing
car_all = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, M = M1,
	spcov_type = 'car', estmethod = 'ml', row_st = TRUE)
# Error: (I - rho * W) * 1 / M must satisfy the CAR symmetry condition
# MD: ***OK*** now -- warning that we override M when row_st = TRUE given
# that Mii = tau2 / Wi,+ (i.e., it is proportional to 1 / rowsums(W))

# automatic row-standardization, M = sum of rows after standardizing
car_all = spautor(z ~ water*tarp, data = caribouDF, W = Nmat1, M = M11,
	spcov_type = 'car', estmethod = 'ml', row_st = TRUE)
# Error: (I - rho * W) * 1 / M must satisfy the CAR symmetry condition
# MD: ***OK*** now -- warning that we override M when row_st = TRUE given
# that Mii = tau2 / Wi,+ (i.e., it is proportional to 1 / rowsums(W))
