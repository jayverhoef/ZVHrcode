#-------------------------------------------------------------------------------
#
#          m2LL
#
#-------------------------------------------------------------------------------

#' minus 2 times loglikelihood for case 12
#'
#' fits a geostatistical linear model
#'
#' @param theta vector of estimated covariance parameters
#' @param z vector of data
#' @param X design matrix for fixed effects
#' @param Z list of design matrices for each random effect 
#' @param xcoords vector with the x-coordinates
#' @param ycoords vector with the y-coordinates
#' @param estMeth estimation method.  Default is "REML" for restricted maximum likelihood.  Other options are "ML" for maximum likelihood
#'
#' @return minus 2 times the loglikelihood
#'
#' @author Jay Ver Hoef
m2LL <- function(theta, z, X, Z, xcoords, ycoords, 
	estMeth, tranf, scale, type, label)
{
		
	covMat <-	makeCovMat(theta = theta, Z1 = Z, Z2 = Z, xcoords1 = xcoords, 
		ycoords1 = ycoords, xcoords2 = xcoords, ycoords2 = ycoords, 
		tranf = tranf, scale = scale, type = type, label = label)

	n <- length(z)
	p <- sum(svd(X)$d>1e-10)
	nX <- length(X[1,])

	logDetV <- as.numeric(determinant(covMat, logarithm = TRUE)$modulus)
	Vi <- solve(covMat)
	XViX <- t(X)%*%Vi%*%X
	logDetXViX <- as.numeric(determinant(XViX, logarithm = TRUE)$modulus)
	covBetaHat <- solve(XViX)
	betaHat <- covBetaHat %*% t(X)%*%Vi%*%z
	r <- z - X %*% betaHat
	rVir <- t(r) %*% Vi %*% r
	minus2LL <- logDetV + rVir + n*log(2*pi)
	if(estMeth == "REML") minus2LL <- minus2LL + logDetXViX - p*log(2*pi)
	return( as.numeric(minus2LL) )
}

