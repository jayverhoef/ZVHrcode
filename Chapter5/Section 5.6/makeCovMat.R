#-------------------------------------------------------------------------------
#
#          makeCovMat
#
#-------------------------------------------------------------------------------

#' creates an appropriate spatial mixed model covariance matrix from parameter vector
#'
#' creates an appropriate spatial mixed model covariance matrix from parameter vector
#'
#' @param theta vector of estimated covariance parameters
#' @param Z1 list of design matrices for each random effect for first data set 
#' @param Z2 list of design matrices for each random effect for second data set 
#' @param xcoords1 vector with the x-coordinates for first data set 
#' @param ycoords1 vector with the y-coordinates for first data set 
#' @param xcoords2 vector with the x-coordinates for second data set 
#' @param ycoords2 vector with the y-coordinates for second data set 
#' @param tranf name of transformation function. Input theta is transformed by
#'				 tranf(theta)*scale
#' @param scale vector of scaling values. Input theta is transformed by
#'				 tranf(theta)*scale
#' @param type type of covariance parameter.  Can be one of "parsil", "range",
#'				"rotate", "minorp", "extrap", "vc", or "nugget"
#' @param label covariance model associated with parameters.
#'
#' @return a covariance matrix
#'
#' @author Jay Ver Hoef
#' @export
makeCovMat <- function(theta, Z1, Z2, xcoords1, ycoords1, xcoords2, ycoords2,
	tranf, scale, type, label, pred = FALSE)
{
	#expand theta by transformation and scaling factors
	theta[tranf == "exp"] <- exp(theta[tranf == 
		"exp"])*scale[tranf == "exp"]
	theta[tranf == "exp/(1+exp)"] <- exp(theta[tranf == "exp/(1+exp)"])/(1 + 
		exp(theta[tranf == "exp/(1+exp)"]))*scale[tranf == "exp/(1+exp)"]
	n <- length(xcoords1)
	covMat <- matrix(0, nrow = length(xcoords1), ncol = length(xcoords2))
	if(sum(type == "range") > 0) {
		spModels <- label[type == "range"]
		for(i in 1:length(spModels)) {
			rotatei <- 0
			minorpi <- 1
			rangei <- theta[type == "range" & 
				label == spModels[i]]
			parsili <- theta[type == "parsil" & 
				label == spModels[i]]
			if(any(type == "rotate" & 
				label == spModels[i])) {
					rotatei <- theta[type == "rotate" & 
						label == spModels[i]]
					minorpi <- theta[type == "minorp" & 
						label == spModels[i]]
			}
			if(any(type == "extrap" & 
				label == spModels[i])) {
					extrapi <- theta[type == "extrap" & 
						label == spModels[i]]
			}
			dismat <- distGeoAni(xcoords1, ycoords1, xcoords2, ycoords2, rotate = rotatei, 
				range = rangei, minorp = minorpi) 
			# compute correlation matrix for scaled distance matrix
			if(spModels[i] == "exponential") CorMat <- corModelExponential(dismat)
			if(spModels[i] == "expRadon2") CorMat <- corModelExpRadon2(dismat)
			if(spModels[i] == "expRadon4") CorMat <- corModelExpRadon4(dismat)
			if(spModels[i] == "gaussian") CorMat <- corModelGaussian(dismat)
			if(spModels[i] == "stable") CorMat <- corModelStable(dismat, extrapi)
			if(spModels[i] == "rationalQuad") CorMat <- corModelRationalQuad(dismat)
			if(spModels[i] == "cauchyGrav") CorMat <- corModelCauchyGrav(dismat)
			if(spModels[i] == "cauchyMag") CorMat <- corModelCauchyMag(dismat)
			if(spModels[i] == "cauchy") CorMat <- corModelCauchy(dismat, extrapi)
			if(spModels[i] == "circular") CorMat <- corModelCircular(dismat)
			if(spModels[i] == "spherical") CorMat <- corModelSpherical(dismat)
			if(spModels[i] == "cubic") CorMat <- corModelCubic(dismat)
			if(spModels[i] == "penta") CorMat <- corModelPenta(dismat)
			if(spModels[i] == "cardinalSine") CorMat <- corModelCardinalSine(dismat)
			if(spModels[i] == "besselK") CorMat <- corModelBesselK(dismat, extrapi)
			if(spModels[i] == "besselJ") CorMat <- corModelBesselJ(dismat, extrapi)
			# create the full covariance matrix with random effects for location and nugget
			covMat <- covMat + parsili*CorMat 
		}
	}
	if(sum(type == "vc") > 0) {
		reModels <- label[type == "vc"]
		for(i in 1:length(reModels)) {
			vc <- theta[type == "vc" & 
				label == reModels[i]] 
			Z1i <- getElement(Z1,reModels[i])
			Z2i <- getElement(Z2,reModels[i])
			covMat <- covMat + vc * Z1i %*% t(Z2i)
		}
	}
	if(!pred) {
		nugget <- theta[type == "nugget"]
		nugplus <- 0
		if(any(type == "parsil")) 
			nugplus <- nugplus + .000001*sum(theta[type == "parsil"]) 
		if(any(type == "vc"))
			nugplus <- nugplus + .000001*sum(theta[type == "vc"]) 
		covMat <- covMat + nugget*diag(n) + nugplus*diag(n)
	}
	return(covMat)
}

