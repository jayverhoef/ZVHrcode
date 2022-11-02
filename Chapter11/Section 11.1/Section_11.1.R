sec_path = 'Rcode/Chapter9/Section 9.2/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  3-D sphere showing arc and chord
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

library("rgl")
options(rgl.printRglwidget = FALSE)

rgl.sphere <- function (x, y=NULL, z=NULL, ng=50, radius = 1, color="white", add=F, ...) {
  lat <- matrix(seq(90, -90, len = ng)*pi/180, ng, ng, byrow = TRUE)
  long <- matrix(seq(-180, 180, len = ng)*pi/180, ng, ng)

  vertex  <- rgl:::rgl.vertex(x, y, z)
  nvertex <- rgl:::rgl.nvertex(vertex)
  radius  <- rbind(vertex, rgl:::rgl.attr(radius, nvertex))[4,]
  color  <- rbind(vertex, rgl:::rgl.attr(color, nvertex))[4,]

  for(i in 1:nvertex) {
    add2 <- if(!add) i>1 else T
    x <- vertex[1,i] + radius[i]*cos(lat)*cos(long)
    y <- vertex[2,i] + radius[i]*cos(lat)*sin(long)
    z <- vertex[3,i] + radius[i]*sin(lat)
    persp3d(x, y, z, specular="white", add=add2, color=color[i], ...)
  }
}


normalize <- function(v) v/sqrt(sum(v^2))

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
#add new sphere
#spheres3d(center, radius = 1, col = "white", alpha = 0.1, 
# fastTransparency = FALSE)
rgl.sphere(0, 0, 0, radius = 1, col = "white", alpha = 0, axes = FALSE,
	xlab = '', ylab = '', zlab = '')

#curve between two points
arc3d(from, to, center, col = "red", lwd = 2)

arc3d(top, to, center, col = 'blue')
arc3d(to, down, center, col = 'blue')
##arc3d(top, to_xy, center, col = 'blue')
##arc3d(to_xy, down, center, col = 'blue')
arc3d(top, from, center, col = 'blue')
arc3d(from, down, center, col = 'blue')
##arc3d(top, from_xy, center, col = 'blue')
##arc3d(from_xy, down, center, col = 'blue')

#straight line
rgl.lines(x,y,z, col = "black", lwd = 2)

#lines3d(c(from[1], to[1]), c(from[2], to[2]),c(from[3], to[3]), col = "black")
arc3d(to, c(-to[1], to[2], to[3]), c(0,0,to[3]), col = 'green')
arc3d(c(-to[1], to[2], to[3]), c(-to[1], -to[2], to[3]), c(0,0,to[3]), col = 'green')
arc3d(c(-to[1], -to[2], to[3]), c(to[1], -to[2], to[3]), c(0,0,to[3]), col = 'green')
arc3d(c(to[1], -to[2], to[3]), c(to[1], to[2], to[3]), c(0,0,to[3]), col = 'green')

arc3d(from, c(-from[1], from[2], from[3]), c(0,0,from[3]), col = 'green')
arc3d(c(-from[1], from[2], from[3]), c(-from[1], -from[2], from[3]), c(0,0,from[3]), col = 'green')
arc3d(c(-from[1], -from[2], from[3]), c(from[1], -from[2], from[3]), c(0,0,from[3]), col = 'green')
arc3d(c(from[1], -from[2], from[3]), c(from[1], from[2], from[3]), c(0,0,from[3]), col = 'green')

# turn on html viewer widget
options(rgl.printRglwidget = TRUE)

#add text
text3d(rbind(from, to), 
       text = c(expression(s[1]), expression(s[2])),
       depth_mask = FALSE, depth_test = "always", cex = 2)

# print image from web browser and use krop to crop it


