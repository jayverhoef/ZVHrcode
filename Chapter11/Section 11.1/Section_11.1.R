sec_path = 'Rcode/Chapter11/Section 11.1/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  3-D sphere showing arc and chord
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

library("rgl")
options(rgl.printRglwidget = TRUE)

normalize <- function(v) v/sqrt(sum(v^2))

normalize_hold3 = function(x) {
	c(x[1]/sqrt((x[1] - x[3])^2 + (x[2]-x[3])^2), 
		 x[2]/sqrt((x[1] - x[3])^2 + (x[2]-x[3])^2),
		 x[3]) 
}

# These vectors all have the same length

from <- normalize(c(1,0.2,0))
to <- normalize(c(1,3,2))
center <- c(0, 0, 0)
top <- c(0,0,1)
down <- c(0,0,-1)
to_xy <- c(-to[1],-to[2], to[3])
from_xy <- c(-from[1], -from[2], from[3])
#seq of points for straight line
x <- seq(from[1], to[1], length.out = 50)
y <- seq(from[2], to[2], length.out = 50)
z <- seq(from[3], to[3], length.out = 50)
#new 3d plot
open3d()


rgl.clear( type = "shapes" )
rgl.sphere(0, 0, 0, radius = 1, col = "white", alpha = 0, axes = FALSE,
	xlab = '', ylab = '', zlab = '')

alpha = 10
zlist = c(-rev(((1:30)/alpha)^2), 0, ((1:30)/alpha)^2)
latcol = 'gray80'
for(i in 1:length(zlist)) {
	lat = c(1,1,zlist[i])
	lat = normalize(lat)

	arc3d(lat, c(-lat[1], lat[2], lat[3]), c(0,0,lat[3]), 
		col = latcol)
	arc3d(c(-lat[1], lat[2], lat[3]), c(-lat[1], -lat[2], lat[3]), 
		c(0,0,lat[3]), col = latcol)
	arc3d(c(-lat[1], -lat[2], lat[3]), c(lat[1], -lat[2], lat[3]), 
		c(0,0,lat[3]), col = latcol)
	arc3d(c(lat[1], -lat[2], lat[3]), c(lat[1], lat[2], lat[3]), 
		c(0,0,lat[3]), col = latcol)
}
#theta = (1:32)/16
#for(i in 1:length(theta)) {
#	lon = c(sin(pi*theta[i]),cos(pi*theta[i]),0)
#	lon = normalize(lon)
#	arc3d(top, lon, center, col = latcol)
#	arc3d(lon, down, center, col = latcol)
#}
latbands_col = 'green4'
lwd_curves = 5
for(i in 1:1) {

	arc3d(top, to, center, col = 'blue', lwd = lwd_curves)
	arc3d(to, down, center, col = 'blue', lwd = lwd_curves)
	arc3d(top, from, center, col = 'blue', lwd = lwd_curves)
	arc3d(from, down, center, col = 'blue', lwd = lwd_curves)

	arc3d(to, c(-to[1], to[2], to[3]), 
		c(0,0,to[3]), col = latbands_col, lwd = lwd_curves )
	arc3d(c(to[1], -to[2], to[3]), c(to[1], to[2], to[3]), 
		c(0,0,to[3]), col = latbands_col, lwd = lwd_curves)
	arc3d(from, normalize_hold3(c(-to[1], to[2], from[3])), 
		c(0,0,from[3]), col = latbands_col, lwd = lwd_curves, lty = 1)
	arc3d(normalize_hold3(c(to[1], -to[2], from[3])), from, 
		c(0,0,from[3]), col = latbands_col, lwd = lwd_curves)

	#curve between two points
	arc3d(from, to, center, col = "red", lwd = lwd_curves)
	#straight line
	rgl.lines(x,y,z, col = "black", lwd = lwd_curves)
}
#add text
text3d(rbind(1.05*from, 1.2*to), 
       text = c(expression(s[1]), expression(s[2])),
       cex = 3, adj = 1)


