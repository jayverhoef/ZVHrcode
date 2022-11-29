#-------------------------------------------------------------------------------
#
#          m2LLi
#
#-------------------------------------------------------------------------------

#' minus 2 times loglikelihood for selected covariance parameter
#'
#' minus 2 times loglikelihood for selected covariance parameter. Used for one-dimensional optimization
#'
#' @param x real-valued variable
#' @param i ith component of vector of covariance parameters
#' @param theta vector of estimated covariance parameters
#' @param z vector of data
#' @param X design matrix for fixed effects
#' @param Z list of design matrices for each random effect 
#' @param xcoords vector with the x-coordinates
#' @param ycoords vector with the y-coordinates
#' @param estMeth estimation method.  Default is "REML" for restricted maximum likelihood.  Other options are "ML" for maximum likelihood
#' @param tranf name of transformation function. Input theta is transformed by
#'				 tranf(theta)*scale
#' @param scale vector of scaling values. Input theta is transformed by
#'				 tranf(theta)*scale
#' @param type type of covariance parameter.  Can be one of "parsil", "range",
#'				"rotate", "minorp", "extrap", "vc", or "nugget"
#' @param label covariance model associated with parameters.
#'
#' @return minus 2 times the loglikelihood
#'
#' @author Jay Ver Hoef
	m2LLi <- function(x, i, theta, z, X, Z, xcoords, ycoords, estMeth,
		tranf, scale, type, label){
		theta[i] <- x
		m2LL(theta, z = z, X = X, Z = Z,
			xcoords = xcoords, ycoords = ycoords, 
			estMeth = estMeth, tranf = tranf, 
			scale = scale, type = type, 
			label = label)
	}

