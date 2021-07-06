#-------------------------------------------------------------------------------
#
#           covParmIni
#
#-------------------------------------------------------------------------------

#' Initializes the covariance parameter vector
#'
#' Initializes the covariance parameter vector
#'
#' @param varComps variance components passed from slmGeo
#' @param useAnisotropy logical vector indicating use of anisotropy for
#'   each component in varComps
#' @param z a vecotr with the observed response variable
#' @param xcoords the vector of the column in the data data.frame with the x-coordinates, default is "x"
#' @param ycol the name of the column in the data data.frame with the y-coordinates, default is "y"
#' @param EstMeth estimation method.  Default is "REML" for restricted maximum likelihood.  Other options are "ML" for maximum likelihood
#' @param CorModel spatial autocorrelation model for errors.  Default is "Exponential".
#' @param group the name of the column for a grouping variable; all groups are considered independent when fitting the model
#' @param use.spatial Fit a spatial model?  Default is "TRUE"
#' @param use.nugget include a nugget effect in parameter estimation?  Default is "TRUE"
#' @param use.anisotropy include anistropy in parameter estimation?  Default is "FALSE"
#' @param use.locrep include variance parameter for replicate measurements per location in parameter estimation?  Default is "FALSE"
#'
#' @return a list with fitted model objects: ...$theta (fitted covariance parameters)  ...$m2LL (minus 2 time the loglikelilhood),	...$Asycov (asymptotic covariance matrix of covariance parameters, from negative inverse Hessian), ...$V (fitted covariance matrix), ...$Vi (inverse of fitted covariance matrix)
#'
#' @author Jay Ver Hoef
covParmIni <- function(varComps, useAnisotropy, z, X, xcoords, ycoords)
{
	# ------------- set up parameter vector with attributes
	# a list of current models that are implemented
	modelList <- c("exponential","expRadon2","expRadon4",
		"gaussian","stable","rationalQuad","cauchyGrav","cauchyMag",
		"cauchy","circular","spherical","cubic","penta","cardinalSine",
		"besselK","besselJ")
	extraList <- c("stable","cauchy","besselK","besselJ")
	# character vector of names of spatial models
	spMods <- NULL
	if(sum(varComps %in% modelList) > 0)
		spMods <- modelList[modelList %in% varComps]
	# column names of any random effects
	reff <- NULL
	if(sum(!(varComps %in% modelList)) > 0)
		reff <- varComps[!(varComps %in% modelList)]
	vCList <- list(spatial = spMods, reff = reff)
	# create a character vector of named covariance components
	#   along with type of parameter
	covLab <- "nugget"
	covType <- "nugget"
	if(!is.null(vCList$spatial)) {
		for(i in 1:length(vCList$spatial)){
			# add partial sill and range parameters common to all models
			covLab <- c(covLab,vCList$spatial[i], vCList$spatial[i])
			covType <- c(covType,"parsil","range")
			# if anisotropy for this model, add rotate and minorp parameters
			if(useAnisotropy[i]){
				covLab <- c(covLab,vCList$spatial[i], vCList$spatial[i])
				covType <- c(covType,"minorp","rotate")
			}
			if(any(vCList$spatial[i] %in%extraList)){
				covLab <- c(covLab,vCList$spatial[i])
				covType <- c(covType,"extrap")			
			}
		}
	}
	if(!is.null(vCList$reff)) {
		for(i in 1:length(vCList$reff)){
			# add variance component parameter for random effect
			covLab <- c(covLab,vCList$reff[i])
			covType <- c(covType,"vc")
		}
	}
	covParms <- matrix(rep(0, times = length(covLab)), ncol  = 1)
	attr(covParms,"type") <- covType
	attr(covParms,"label") <- covLab
	rownames(covParms) <- 1:length(covLab)
	colnames(covParms) <- "parms"
	scaleFactors <- rep(1, times = length(covParms))
	vcs <- attr(covParms,"type") == "vc" | attr(covParms,"type") == 
		"parsil" | attr(covParms,"type") == "nugget"
	scaleFactors[vcs] <- var(z)/sum(vcs)
	scaleFactors[attr(covParms,"type") == "range"]  <- 4*
		sqrt((var(xcoords) + var(ycoords))/2)
	scaleFactors[attr(covParms,"type") == "rotate"]  <- 180
	scaleFactors[attr(covParms,"type") == "extrap"]  <- 2
	attr(covParms,"scale") <- scaleFactors
	transfType <- rep("exp", times = length(covParms))
	transfType[attr(covParms,"type") == "rotate"] <- "exp/(1+exp)"
	transfType[attr(covParms,"type") == "minorp"]  <- "exp/(1+exp)"
	transfType[attr(covParms,"type") == "extrap"] <- "exp/(1+exp)"
	attr(covParms,"tranf") <- transfType
	covParms
}

