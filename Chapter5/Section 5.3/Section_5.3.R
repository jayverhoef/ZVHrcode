sec_path = 'Rcode/Chapter5/Section 5.3/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#              Table 5.2 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Compute some results in Section 4.2
xpts <- matrix(c(1, 1, 2, 2, 5), 5, 1)
ypts <- matrix(c(4, 3, 3, 2, 4), 5, 1)
Dmat = as.matrix(dist(cbind(xpts, ypts)))
Rmat = (1 - 3*Dmat/8 + Dmat^3/128)*(Dmat < 4)
Rmat
Rhalf <- eigen(Rmat)$vectors %*% diag(sqrt(eigen(Rmat)$values)) %*% 
	t(eigen(Rmat)$vectors)
y <- c(1,0,3,1,5)
X <- matrix(1, 5, 1)

RhalfiX <- solve(Rhalf)%*%X
RiX <- solve(Rmat)%*%X
P <- RiX %*% solve(t(RiX) %*% X) %*% t(RiX)
QQ <- solve(Rmat) %*% (y - X %*% mutilde)


#Obtain ordinary and generalized leverages as listed on page 104
Xplanar <- cbind(X, xpts, ypts)
planarreg <- lm(y ~ Xplanar)
RhalfiXplanar <- solve(Rhalf) %*% Xplanar
hsharpplanar <- diag(RhalfiXplanar %*% 
	solve(t(RhalfiXplanar) %*% RhalfiXplanar) %*% t(RhalfiXplanar))
genlev = data.frame(
	hii = hatvalues(planarreg),
	hsharp = hsharpplanar
)

print(
    xtable(genlev, 
      align = c('l',rep('l', times = length(genlev[1,]))),
      digits = c(0, rep(3, times = 2))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

