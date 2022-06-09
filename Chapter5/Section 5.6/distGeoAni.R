#-------------------------------------------------------------------------------
#
#           distGeoAni
#
#-------------------------------------------------------------------------------

#' Compute anistropy corrected distance between two sets of data
#'
#' computes anistropy corrected distance between two sets of data
#'
#' @param xrow vector with x-coordinates that will form rows of distance matrix 
#' @param yrow vector with y-coordinates that will form rows of distance matrix, must be of same length as xrow
#' @param xcol vector with x-coordinates that will form columns of distance matrix 
#' @param ycol vector with y-coordinates that will form columns of distance matrix, must be of same length as xcol
#' @param rotate rotation of anisotropic axes, default = 0
#' @param range range of autocorrelation model, default = 1
#' @param minorp proportion of range in x direction to that of y direction for unrotated anisotropic model, default = 1
#'
#' @return matrix of distances
#'
#' @author Jay Ver Hoef
distGeoAni <- function(xrow, yrow, xcol, ycol, rotate = 0, range = 1, minorp = 1)
{
	# expand all x-coordinates
	  sxr = matrix(xrow, nrow = length(xrow), ncol = length(xcol))
		sxc <- matrix(xcol, nrow = length(xrow), ncol = length(xcol), 
			byrow = TRUE)
		syr = matrix(yrow, nrow = length(yrow), ncol = length(ycol))
		syc <- matrix(ycol, nrow = length(yrow), ncol = length(ycol), 
			byrow = TRUE)
	# find difference in coordinates between all pairwise locations
		sxdif <- sxr - sxc
		sydif <- syr - syc
	# rotate coordinates
		newx <- cos(rotate*.0174533)*sxdif - sin(rotate*.0174533)*sydif
		newy <- sin(rotate*.0174533)*sxdif + cos(rotate*.0174533)*sydif
	# scale coordinates by minor and major axes */
		newx <- newx/(range*minorp)
		newy <- newy/range
	# compute distance for the scaled and rotated coordinates */
		sqrt(newx^2 + newy^2)
}
