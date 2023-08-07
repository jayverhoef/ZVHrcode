sec_path = 'Rcode/Chapter10/Section 12.4/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        Cokriging weights and variances versus kriging weights and variance
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#
# cokriging weights and variance
#

# spatial locations along a transect, including prediction location last
s = c(0, 2, 5, 3)
# distance between locations
dis = as.matrix(dist(s))
# covariance matrices
C11 = C22 = exp(-dis/2)
# i is row and for variable 1, j is column for variable 2, 
# and the differences in i - j for s = c(0, 2, 5, 3)
ijdiff = matrix(rep(s, times = 4), nrow = 4) - 
	t(matrix(rep(s, times = 4), nrow = 4))
ijdiff
# now abs(ijdiff) is the standard distance, but we want
# distances computed with the shift a_1 = -1 and a_2 = 1
# and then take the absolute value
dis1 = abs(ijdiff - 1)
dis2 = abs(ijdiff + 1)
# note that "distance" between i = 4 (prediction location, variable 1)
# and j = 2 (second observed location, variable 2) is now zero for dis1
dis1
# and here are the "distances" when a_2 = 1, which is the transpose of dis1
dis2
dis1 == t(dis2)
# cross-covariance matrices using asymmetric distances
C12 = exp(-dis1/2)
C21 = exp(-dis2/2)
# note that unscaled cross-correlation between variable 1 (rows) and
# variable 2 (columns) is 1 at C[4,2]
C12
# joint covariance matrix, scale the cross-covariances by 0.7
Sigall = rbind(cbind(C11,0.7*C12),cbind(0.7*C21,C22))
# covariance matrix subsetted to observed data
Sig = Sigall[c(1:3,5:7),c(1:3,5:7)]
# covariance between prediction location (variable 1, 4th index)
# and observed data (columns 1:3, skip 4, and then 5:7)
# from full joint covariance matrix
c12 = Sigall[4, c(1:3,5:7)]
## inverse of the covariance matrix for observed locations
Sigi = solve(Sig)
# design matrix for observed locations
X = rep(1, times = 6)
# weights for estimating the mean
bhat_wts = (Sigi %*% X)/sum( Sigi)
# prediction weights
round(t(bhat_wts) + c12 %*% Sigi %*% (diag(6) - X %*% t(bhat_wts)), 3)
# prediction variance (using short-cuts because X is all 1's)
1 - c12 %*% Sigi %*% c12 + (1 - sum(Sigi %*% c12))^2/sum(Sigi)

#
# kriging weights and variance
#

# covariance matrix subsetted to observed data for variable 1
Sig = Sigall[1:3,1:3]
# inverse of the covariance matrix for observed locations
Sigi = solve(Sig)
# covariance between prediction location and observed data
c12 = Sigall[4, 1:3]
# design matrix for observed locations
X = rep(1, times = 3)
# weights for estimating the mean
bhat_wts = (Sigi %*% X)/as.numeric((t(X) %*% Sigi %*% X))
# prediction weights
t(bhat_wts) + c12 %*% Sigi %*% (diag(3) - X %*% t(bhat_wts))
# prediction variance (using short-cuts because X is all 1's)
1 - c12 %*% Sigi %*% c12 + (1 - sum(Sigi %*% c12))^2/sum(Sigi)
