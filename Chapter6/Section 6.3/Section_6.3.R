sec_path = 'Rcode/Chapter6/Section 6.3/figures/'
setwd(paste0(SLEDbook_path,sec_path))

library(gstat)
library(viridis)
library(classInt)
source('addBreakColorLegend.R')

################################################################################
#-------------------------------------------------------------------------------
#              Simulation Anistropic Plots
#-------------------------------------------------------------------------------
################################################################################

file_name = "SimAniso"

pdf(paste0(file_name,'.pdf'), width = 14, height = 8)

cex_main = 2.0    
source('addBreakColorLegend.R')
layout(matrix(1:12, ncol = 6, byrow = TRUE), widths = c(4,1,4,1,4,1))

# create a grid of points
xy <- expand.grid(seq(0.1, 10, length.out = 100), seq(0.1, 10, length.out = 100))
names(xy) <- c("x", "y")


#-------------------------------------------------------------------------------
## Geometric anisotropic model c(45, 1/4)
#-------------------------------------------------------------------------------

# gstat uses geographical coordinates, so anisotropy of 45 degrees is
# pi/4 in polar coordinates, the other anisotropy parameter is the ratio
# of the y-axis (long axis) to x-axis (short axis) in gstat parameterization
# but these are reversed in polar coordinates, as, prior to rotation,
# the y-axis becomes the short axis by the ratio 1/4 to the longer x-axis.
# These same considerations hold for all geometric anisotropy plots in this
# figure
g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(1, "Exp", 2/3, anis = c(45, 1/4)), nmax = 20)
set.seed(171)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Exponential(", alpha, " = 2/3, ", 
		phi,"=", pi,"/4, ", tau, "=1/4)"))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Geometric anisotropic model c(60, 1/4)
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(1, "Exp", 2/3, anis = c(60, 1/4)), nmax = 20)
set.seed(172)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Exponential(", alpha, " = 2/3, ", 
		phi,"=", pi,"/6, ", tau, "=1/4)"))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Exponential Separable model theta_2 = 1, theta_3 = 2)
#-------------------------------------------------------------------------------

# a couple of functions for separable models
# based on the special pattern of the separable covariance matrix. 
my_E1 <- function(column1, theta1 = 1, theta3){
    mat <- matrix(NA, nrow = dim(column1)[1], ncol = dim(column1)[1])
    for (i in 1 : dim(column1)[1]){
        for (j in 1 : dim(column1)[1]){
            diff <- abs(column1[i, 2] - column1[j, 2])
            mat[i, j] <- theta1 * exp(- theta3 * diff)
        }
    }
    mat
} 

my_E2 <- function(row1, theta2){
    mat <- matrix(NA, nrow = dim(column1)[1], ncol = dim(column1)[1])
    for (i in 1 : dim(column1)[1]){
        for (j in 1 : dim(column1)[1]){
            diff <- abs(row1[i, 1] - row1[j, 1])
            mat[i, j] <- exp(- theta2 * diff)
        }
    }
    mat
} 

column1 <- cbind(rep(0.1, 100), seq(0.1, 10, length.out = 100))
row1 <- cbind(seq(0.1, 10, length.out = 100), rep(0.1, 100))

E1 <- my_E1(column1, 1, 2)
E1_sqrt <- chol(E1)

E2 <- my_E2(row1, 1)
E2_sqrt <- chol(E2)

V_sqrt <- kronecker(E1_sqrt, E2_sqrt)

set.seed(173)
temp <- rnorm(10000, 0, 1)
sim1 <- V_sqrt %*% temp
simdata1 <- cbind(xy, sim1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Separable (", alpha[u], " = 1, ", 
		alpha[v], " = 2)"))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Geometric anisotropic model c(45, 1/2)
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(1, "Exp", 2/3, anis = c(45, 1/2)), nmax = 20)
set.seed(179)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Exponential(", alpha, " = 2/3, ", 
		phi,"=", pi,"/4, ", tau, "=1/2)"))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Geometric anisotropic model c(60, 1/4)
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(1, "Exp", 2/3, anis = c(60, 1/2)), nmax = 20)
set.seed(175)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Exponential(", alpha, " = 2/3, ", 
		phi,"=", pi,"/6, ", tau, "=1/2)"))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Exponential Separable model theta_2 = 2, theta_3 = 1)
#-------------------------------------------------------------------------------


E1 <- my_E1(column1, 2, 1)
E1_sqrt <- chol(E1)

E2 <- my_E2(row1, 2)
E2_sqrt <- chol(E2)

V_sqrt <- kronecker(E1_sqrt, E2_sqrt)

set.seed(176)
temp <- rnorm(10000, 0, 1)
sim1 <- V_sqrt %*% temp
simdata1 <- cbind(xy, sim1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Separable (", alpha[u], " = 2, ", 
		alpha[v], " = 1)"))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

par(old.par)
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

