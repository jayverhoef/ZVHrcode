sec_path = 'Rcode/Chapter9/Section 9.5/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        Table 9.4 Weights for ordinary kriging using exponential
#  model of various strengths, nuggets, anisotropies.  IDW weights and
#  variances are also included
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Some useful functions


# universal kriging equations (will work for ordinary kriging too)
# X is design matrix
# XRiX_i = solve(t(X) %*% solve(Covmat) %*% X) [no need to solve for every 
#      prediction
# Ri = solve(Covmat) [no need to solve for every prediction]
# x0 = design vector for prediction site
# r0 = covariance between observed data and prediction site
# y = observed data
ukpred = function(X, XRiX_i, Ri, x0, r0, y){
	t(x0) %*% XRiX_i %*% t(X) %*% Ri %*% y +
	t(r0) %*% Ri %*% (diag(length(y)) - X %*% 
	XRiX_i %*% t(X) %*% Ri) %*% y
}

# universal kriging equations (will work for ordinary kriging too)
# X is design matrix
# XRiX_i = solve(t(X) %*% solve(Covmat) %*% X) [no need to solve for every 
#      prediction
# Ri = solve(Covmat) [no need to solve for every prediction]
# x0 = design vector for prediction site
# r0 = covariance between observed data and prediction site
ukpredwts = function(X, XRiX_i, Ri, x0, r0){
	t(x0) %*% XRiX_i %*% t(X) %*% Ri +
	t(r0) %*% Ri %*% (diag(length(y)) - X %*% 
	XRiX_i %*% t(X) %*% Ri)
}

# universal kriging prediction variances
# X is design matrix
# XRiX_i = solve(t(X) %*% solve(Covmat) %*% X) [no need to solve for every 
#      prediction
# Ri = solve(Covmat) [no need to solve for every prediction]
# x0 = design vector for prediction site
# r0 = covariance between observed data and prediction site
# r00 = variance at prediction location
ukpev = function(X, XRiX_i, Ri, x0, r0, r00){
	r00 - t(r0) %*% Ri %*% r0 + (t(x0) - t(r0) %*%
	Ri %*% X) %*% XRiX_i %*% 
	(x0 - t(X) %*% Ri %*% r0)
}


# basic data set up
x = c(1, 1, 2, 2, 5)
y = c(4, 3, 3, 2, 4)
addx = 2
addy = 4
X <- rep(1, times = 5)
distMat = as.matrix(dist(cbind(c(x,addx), c(y, addy))))


# Model 1 with practical range = 2

Rall = exp(-3*distMat/2)
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
# kriging weights
okwts_mod1 = ukpredwts(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0)
# prediction variance
okpev_mod1 = ukpev(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0, 1)


# Model 2 (effective range twice as large = 4)

Rall = exp(-3*distMat/4)
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
# kriging weights
okwts_mod2 = ukpredwts(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0)
# prediction variance
okpev_mod2 = ukpev(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0, 1)


# Model 3 (effective range 4 times as large = 8)

Rall = exp(-3*distMat/8)
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
# kriging weights
okwts_mod3 = ukpredwts(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0)
# prediction variance
okpev_mod3 = ukpev(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0, 1)


# Model 4 (25% nugget, but same sill)

Rall = 0.75*exp(-3*distMat/4) + diag(rep(0.25, times = 6))
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
# kriging weights
okwts_mod4 = ukpredwts(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0)
# prediction variance
okpev_mod4 = ukpev(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0, 1)

# Model 5 (50% nugget, but same sill)

Rall = 0.50*exp(-3*distMat/4) + diag(rep(0.50, times = 6))
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
# kriging weights
okwts_mod5 = ukpredwts(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0)
# prediction variance
okpev_mod5 = ukpev(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0, 1)



# Model 6 (geometrically anisotropic with major axis in E-W direction)
# basic data set up
distMatGA = as.matrix(dist(cbind(c(x,addx)/2, c(y, addy))))
Rall = exp(-3*distMatGA/4)
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
# kriging weights
okwts_mod6 = ukpredwts(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0)
# prediction variance
okpev_mod6 = ukpev(X, solve(t(X) %*% solve(R, X)), solve(R), 1, r0, 1)


#-------------------------------------------------------------------------------

# IDW and IDSW for example of Section 9.6

dist0 = distMat[1:5,6]
# inverse distance weights
IDWwts <- (1/sum(1/dist0))*(1/dist0)
# inverse distance squared weights
IDSWwts <- (1/sum(1/dist0^2))*(1/dist0^2)

# IDW and IDSW variance under Model 1
Rall = exp(-3*distMat/2)
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
#variances 
IDWpev_mod1 <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev_mod1 <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1

# IDW and IDSW variance under Model 2
Rall = exp(-3*distMat/4)
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
#variances 
IDWpev_mod2 <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev_mod2 <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1

# IDW and IDSW variance under Model 3
Rall = exp(-3*distMat/8)
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
#variances 
IDWpev_mod3 <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev_mod3 <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1

# IDW and IDSW variance under Model 4
Rall = 0.75*exp(-3*distMat/4) + diag(rep(0.25, times = 6))
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
#variances 
IDWpev_mod4 <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev_mod4 <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1

# IDW and IDSW variance under Model 5
Rall = 0.50*exp(-3*distMat/4) + diag(rep(0.50, times = 6))
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
#variances 
IDWpev_mod5 <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev_mod5 <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1

# IDW and IDSW variance under Model 6
Rall = exp(-3*distMatGA/4)
# observed locations
R = Rall[1:5, 1:5]
# prediction location
r0 = Rall[1:5, 6]
#variances 
IDWpev_mod6 <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev_mod6 <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1

# create the table
OKandIDW_wts = cbind(
	c(okwts_mod1, okpev_mod1, IDWpev_mod1, IDSWpev_mod1),
	c(okwts_mod2, okpev_mod2, IDWpev_mod2, IDSWpev_mod2),
	c(okwts_mod3, okpev_mod3, IDWpev_mod3, IDSWpev_mod3),
	c(okwts_mod4, okpev_mod4, IDWpev_mod4, IDSWpev_mod4),
	c(okwts_mod5, okpev_mod5, IDWpev_mod5, IDSWpev_mod5),
	c(okwts_mod6, okpev_mod6, IDWpev_mod6, IDSWpev_mod6),
	c(IDWwts, NA, NA, NA),
	c(IDSWwts, NA, NA, NA)
)

# 
print(
    xtable(OKandIDW_wts, 
      align = c('l',rep('l', times = length(OKandIDW_wts[1,]))),
      digits = c(0, rep(3, times = 8))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

