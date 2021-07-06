#-------------------------------------------------------------------------------
#
#           splmm
#
#-------------------------------------------------------------------------------

#' Fits geostatistical linear model
#'
#' fits a geostatistical linear model
#'
#' @param formula an R linear model formula
#' @param spdata an sp SpatialPointsDataFrame
#' @param estMeth estimation method.  Default is "REML" for restricted maximum likelihood.  Other options are "ML" for maximum likelihood
#' @param varComps a list of variance components, including spatial autocorrelation models 
#' for errors and traditional random effects.  The list of spatial autocorrelation 
#' models is "exponential","expRadon2","expRadon4","gaussian","stable",
#' "rationalQuad","cauchyGrav","cauchyMag","cauchy","circular","spherical",
#' "cubic","penta","cardinalSine","besselK","besselJ"  Default is "exponential".
#' Any names in the list not given above will be searched among the columns in the data set and
#' used as a factor variable for levels of a traditional random effect.
#' @param useAnisotropy include anistropy in parameter estimation?  Default is "FALSE"
#'
#' @return a list of class "splmm".  The functions "summary" and "print" are used to obtain and print a summary. "anova" returns just the analysis of variance table...
#'
#' @author Jay Ver Hoef
#' @export
splmm <- function(formula, spdata, estMeth = "REML", varComps = "exponential", 
	useAnisotropy = FALSE) 
{

	# ----------------------------------------------------------------------------
	# prepare data
	# ----------------------------------------------------------------------------
	#get name of response variable in data set
	if(class(spdata) != "SpatialPointsDataFrame")
		return("Error: spdata must be have class SpatialPointsDataFrame")		
	data <- spdata@data
	trms <- terms(formula, data = data)
	respCol <- as.character(as.list(attr(trms,"variables")[-1]))[1]
	covList <- attr(trms,"term.labels")
	#make sure the response variable is numeric
	if(is.factor(data[,respCol]) | is.character(data[,respCol]))
		return("Error: response variable is character or factor")
	#design matrix
	X <- model.matrix(formula, data)
	modelList <- c("exponential","expRadon2","expRadon4",
		"gaussian","stable","rationalQuad","cauchyGrav","cauchyMag",
		"cauchy","circular","spherical","cubic","penta","cardinalSine",
		"besselK","besselJ")
	spModels <- NULL
	spModels <- varComps[varComps %in% modelList]
	reModels <- NULL
	reModels <- varComps[varComps %in% names(data)]
	if(length(reModels) > 0L) {
		Z <- rep(list(),length(reModels))
		for(i in 1L:length(reModels)) {
			data[,reModels[i]] <- as.factor(as.character(data[,reModels[i]]))
			Z[[i]] <- model.matrix(as.formula(paste("~ -1 +",reModels[i])),data)
		}
	}
	if(length(reModels) > 0) names(Z) <- reModels
	#subset to all records without missing covariates or random effects
	allRowNames <- rownames(data)
	ind <- rep(TRUE, times = length(allRowNames))
	ind <- allRowNames %in% rownames(X) & ind
	if(length(reModels) > 0) {
		for(i in 1L:length(reModels)) 
			ind <- allRowNames %in% rownames(Z[[i]]) & ind
	}
	X <- X[allRowNames[ind], , drop = F]
	if(length(reModels) > 0) {
		for(i in 1L:length(reModels)) 
			Z[[i]] <- Z[[i]][allRowNames[ind],]
	}
	# attach logical attribute to data indicating missing data
	#   due to either missing response or covariates
 	attr(data,"NAs") <- !ind
	#vector of response variable
	z <- as.matrix(data[allRowNames[ind], respCol, drop = FALSE])
	if(any(is.infinite(z))) return("Response variable has infinite values")
	#vectors of spatial coordinates
	xcoords <- spdata@coords[ind,1]
	ycoords <- spdata@coords[ind,2]

	n <- length(z)
	p <- sum(svd(X)$d>1e-10)
	nX <- length(X[1,])

	# ----------------------------------------------------------------------------
	# initialize vector of covariance parameters
	# ----------------------------------------------------------------------------

	thetaIni <- covParmIni(varComps = varComps, useAnisotropy = useAnisotropy, 
		z = z, X = X, xcoords = xcoords, ycoords = ycoords)

	if(length(thetaIni) == 1L) { 
		parmest <- optimize(m2LLi, interval = c(-5,5), i = 1, 
			theta = thetaIni, z = z, X = X, 
			Z = Z, xcoords = xcoords, ycoords = ycoords, 
			estMeth = estMeth, tranf = attr(thetaIni,"tranf"), 
			scale = attr(thetaIni,"scale"), type = attr(thetaIni,"type"), 
			label = attr(thetaIni,"label"))
		theta <- thetaIni
		theta[1] <- parmest$minimum
		thetap <- parmest$minimum
		m2LL <- parmest$objective
	} else {
	#do some 1-D optimizations first
		for(j in 1:3L){
			for(i in 1L:length(thetaIni)) {
				thetaIni[i] <- as.numeric(optimize(m2LLi, interval = c(-5,5), i = i, 
					theta = thetaIni, z = z, X = X, 
					Z = Z, xcoords = xcoords, ycoords = ycoords, 
					estMeth = estMeth, tranf = attr(thetaIni,"tranf"), 
					scale = attr(thetaIni,"scale"), type = attr(thetaIni,"type"), 
					label = attr(thetaIni,"label"))$minimum)
			}
		}
		parmest <- optim(thetaIni, m2LL, z = z, X = X, Z = Z,
			xcoords = xcoords, ycoords = ycoords, 
			estMeth = estMeth, tranf = attr(thetaIni,"tranf"), 
			scale = attr(thetaIni,"scale"), type = attr(thetaIni,"type"), 
			label = attr(thetaIni,"label"), hessian = TRUE)
		m2LL <- parmest$value
		thetap <- parmest$par
		theta <- parmest$par
	}		
		theta[attr(theta,"tranf") == "exp"] <- exp(theta[attr(theta,"tranf") == 
			"exp"])*attr(thetaIni,"scale")[attr(theta,"tranf") == "exp"]
		theta[attr(theta,"tranf") == "exp/(1+exp)"] <- exp(theta[attr(theta,
			"tranf") == "exp/(1+exp)"])/(1 + exp(theta[attr(theta,"tranf") == 
			"exp/(1+exp)"]))*attr(theta,"scale")[attr(theta,"tranf") == 
			"exp/(1+exp)"]
		nugplus <- 0
		if(any(attr(theta,"type")=="parsil")) 
			nugplus <- nugplus + .000001*sum(theta[attr(theta,"type")=="parsil"]) 
		if(any(attr(theta,"type")=="vc"))
			nugplus <- nugplus + .000001*sum(theta[attr(theta,"type")=="vc"])
		theta[attr(theta,"type")=="nugget"] <- theta[attr(theta,"type")=="nugget"] +
			nugplus

	covMat <-	makeCovMat(theta = thetap, Z1 = Z, Z2 = Z, xcoords1 = xcoords, 
		ycoords1 = ycoords, xcoords2 = xcoords, ycoords2 = ycoords, 
		tranf = attr(theta,"tranf"), scale = attr(theta,"scale"), 
		type = attr(theta,"type"), label = attr(theta,"label"))
	Vi <- solve(covMat)
	XViX <- t(X)%*%Vi%*%X
	covBetaHat <- solve(XViX)
	betaHat <- covBetaHat %*% t(X)%*%Vi%*%z

	outpt <- list(
		call = match.call(),
		terms = trms,
		formula = formula,
		data = data,
		xcoords = spdata@coords[,1],
		ycoords = spdata@coords[,2],
		data.sample.size = length(data[,1]),
		obs.sample.size = length(rownames(X)),
		rank = p,
		varComps = varComps,
		estMeth = estMeth,
		useAnisotropy = useAnisotropy,
		theta = theta,
		covBetaHat = covBetaHat,
		betaHat = betaHat,
		z = z,
		V = covMat,
		Vi = Vi,
		XViX = XViX,
		X = X,
		m2LL = m2LL
	)
	class(outpt) <- "splmm"
	outpt
}

