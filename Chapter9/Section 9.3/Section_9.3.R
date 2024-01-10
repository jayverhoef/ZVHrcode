sec_path = 'Rcode/Chapter9/Section 9.3/'
setwd(paste0(SLEDbook_path,sec_path))

library(spmodel)
library(xtable)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#          Table 9.3
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Code for universal kriging example (with planar mean surface), spherical model in Table 9.3
x = c(1, 1, 2, 2, 5)
y = c(4, 3, 3, 2, 4)
addx = c(2, 4, 4)
addy = c(4, 2, 4)
X <- cbind(rep(1, times = 5), x, y)
x01 = c(1, addx[1], addy[1])
x02 = c(1, addx[2], addy[2])
x03 = c(1, addx[3], addy[3])

distMat = as.matrix(dist(cbind(c(x,addx), c(y, addy))))
Rall = 1 - (3*distMat/8) + (distMat^3/128)
R = Rall[1:5, 1:5]
r01 = Rall[1:5, 6]
r02 = Rall[1:5, 7]
r03 = Rall[1:5, 8]

# function to compute universal kriging weights
ukwts = function(X, R, r0, x0){
	t(x0) %*% solve(t(X) %*% solve(R) %*%X) %*% t(X) %*% solve(R) +
		t(r0) %*% solve(R) %*% (diag(dim(X)[1]) - X %*% 
		solve(t(X) %*% solve(R) %*% X) %*% t(X) %*% solve(R))
}
	
# function to compute prediction variance for universal kriging
ukpev = function(X, R, r0, x0, r00){
	r00 - t(r0) %*% solve(R) %*% r0 + (t(x0) - t(r0) %*%
		solve(R) %*% X) %*% solve(t(X) %*% solve(R) %*% X) %*% 
		(x0 - t(X) %*% solve(R) %*% r0)
}

table_ukwts_pev = cbind(
	c(ukwts(X, R, r01, x01), ukpev(X, R, r01, x01, 1)),
	c(ukwts(X, R, r02, x02), ukpev(X, R, r02, x02, 1)),
	c(ukwts(X, R, r03, x03), ukpev(X, R, r03, x03, 1))
)


print(
    xtable(table_ukwts_pev, 
      align = c('l',rep('l', times = length(table_ukwts_pev[1,]))),
      digits = c(0, rep(3, times = 3))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

