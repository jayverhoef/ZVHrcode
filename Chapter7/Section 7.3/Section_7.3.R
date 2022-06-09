sec_path = 'Rcode/Chapter7/Section 7.3/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)

################################################################################
#-------------------------------------------------------------------------------
#              Standardization of Spatial Weights
#-------------------------------------------------------------------------------
################################################################################

# Computation of marginal variances and correlations for the SAR and CAR models 
# on the 7-site example of Section 7.2

# create W-matrix
W0 <- matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,byrow=T)
id3 <- diag(3)
ze3 <- matrix(0,3,3)
W1 <- cbind(W0,id3,ze3)
W2 <- cbind(id3,W0,id3)
W3 <- cbind(ze3,id3,W0)
W <- rbind(W1,W2,W3)
W <- W[-1,-1]
W <- W[-2,-2]
# row-standardized weights
Wbar = W/apply(W,1,sum)
# row-standardized K matrix
Kbar = diag(1/apply(W,1,sum))	
# SAR covariance and correlation matrices for row-standardized weights
SigmaWbarSAR <- solve(diag(7) - 0.8*Wbar)%*%Kbar%*%solve(diag(7) - 0.8*t(Wbar))
corrWbarSAR = diag(sqrt(1/diag(SigmaWbarSAR))) %*% 
	SigmaWbarSAR%*%diag(sqrt(1/diag(SigmaWbarSAR)))
# create matrix with variances on diagonal and correlations on off-diagonal
SigmaWbarSAR[lower.tri(SigmaWbarSAR)] = NA
SigmaWbarSAR[upper.tri(SigmaWbarSAR)] = corrWbarSAR[upper.tri(SigmaWbarSAR)]
# CAR covariance and correlation matrices for row-standardized weights
SigmaWbarCAR <- solve(diag(7) - 0.8*Wbar) %*% Kbar
corrWbarCAR = diag(sqrt(1/diag(SigmaWbarCAR))) %*% 
	SigmaWbarCAR%*%diag(sqrt(1/diag(SigmaWbarCAR)))
# create matrix with variances on diagonal and correlations on off-diagonal
SigmaWbarCAR[lower.tri(SigmaWbarCAR)] = NA
SigmaWbarCAR[upper.tri(SigmaWbarCAR)] = corrWbarCAR[upper.tri(corrWbarCAR)]

# print matrices
print(
    xtable(SigmaWbarSAR, 
      align = c('c',rep('c', times = length(W[1,]))),
      digits = c(0,rep(2, times = 7))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

print(
    xtable(SigmaWbarCAR, 
      align = c('c',rep('c', times = length(W[1,]))),
      digits = c(0,rep(2, times = 7))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

################################################################################
#-------------------------------------------------------------------------------
#              Variance Stabilizing Spatial Weights
#-------------------------------------------------------------------------------
################################################################################

# variance stabilized weights
Wvstab = W/sqrt(apply(W,1,sum))
# variance stabilized K matrix
Kvstab = diag(1/sqrt(apply(W,1,sum)))	

# SAR covariance and correlation matrices for row-standardized weights
SigmaWvstabSAR <- solve(diag(7) - 0.5*Wvstab)%*% Kvstab %*% 
	solve(diag(7) - 0.5*t(Wvstab))
corrWvstabSAR = diag(sqrt(1/diag(SigmaWvstabSAR))) %*% 
	SigmaWvstabSAR%*%diag(sqrt(1/diag(SigmaWvstabSAR)))
# create matrix with variances on diagonal and correlations on off-diagonal
SigmaWvstabSAR[lower.tri(SigmaWvstabSAR)] = NA
SigmaWvstabSAR[upper.tri(SigmaWvstabSAR)] = corrWvstabSAR[upper.tri(SigmaWvstabSAR)]
# CAR covariance and correlation matrices for row-standardized weights
SigmaWvstabCAR <- solve(diag(7) - 0.5*Wvstab) %*% Kvstab
corrWvstabCAR = diag(sqrt(1/diag(SigmaWvstabCAR))) %*% 
	SigmaWvstabCAR%*%diag(sqrt(1/diag(SigmaWvstabCAR)))
# create matrix with variances on diagonal and correlations on off-diagonal
SigmaWvstabCAR[lower.tri(SigmaWvstabCAR)] = NA
SigmaWvstabCAR[upper.tri(SigmaWvstabCAR)] = corrWvstabCAR[upper.tri(corrWvstabCAR)]

# print matrices
print(
    xtable(SigmaWvstabSAR, 
      align = c('c',rep('c', times = length(W[1,]))),
      digits = c(0,rep(2, times = 7))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

print(
    xtable(SigmaWvstabCAR, 
      align = c('c',rep('c', times = length(W[1,]))),
      digits = c(0,rep(2, times = 7))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)






