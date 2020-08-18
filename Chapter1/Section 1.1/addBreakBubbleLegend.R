#-------------------------------------------------------------------------------
#
#           addBreakColorLegend
#
#-------------------------------------------------------------------------------

#' Adds a Legend to Plots with Bubble Sizes and Breaks
#'
#' Adds a Legend to Plots with Bubble Sizes and Breaks
#'
#' @param xleft a scalar of the left x position.
#' @param ybottom a scalar of bottom y position.
#' @param xright a scalar of the right x position.
#' @param ytop a scalar of top y position.
#' @param brks the breaks between color classes.  Needs a lower bound, and upper bound, so there should be one more break than number of classes.
#' @param cexes a vector of character expansions (cex), should have one less element than the number of breaks.
#' @param printFormat: a character variable in '4.2' format where the 4 control the number of digits before the decimal, and 2 controls the number of digits after the decimal.
#'
#' @seealso \code{\link{plotPointsRGB}}, \code{\link{rect}}, 

#' @return add a color ramp legend as a rectangle to the currently active plot
#'
#' @author Jay Ver Hoef
#' @rdname addBreakColorLegend
#' @export addBreakColorLegend 


addBreakBubbleLegend <- function(ybottom, ytop, bubble_labels, 
	cexes, printFormat = "4.2", ...) 
{
  old.par = par(mar = c(0,0,0,0))
  ncexes <- length(cexes)
  plot(c(0,1),c(0,1), type ='n', xlab = '', ylab ='',xaxt = 'n', yaxt = 'n',
    bty = 'n')
  for(i in 1:ncexes)
    points(0.3, ybottom + (i-1)/ncexes*(ytop - ybottom), cex = cexes[i])
  for(i in 1:ncexes)
    text(0.5, ybottom + (i-1)/ncexes*(ytop - ybottom), bubble_labels[i])
  par(old.par)
  return(invisible())	
}
