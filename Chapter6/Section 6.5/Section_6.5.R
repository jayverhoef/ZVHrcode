sec_path = 'Rcode/Chapter6/Section 6.5/'
setwd(paste0(SLEDbook_path,sec_path))

library(gstat)

################################################################################
#-------------------------------------------------------------------------------
#              Semivariogram Figure
#-------------------------------------------------------------------------------
################################################################################

file_name = "figures/Semivariograms"

pdf(paste0(file_name,'.pdf'), width = 11, height = 8.5)

###### Semivariograms (Figure 6.6) ###### 

layout(matrix(1:4, ncol = 2, byrow = TRUE))

#-------------------------------------------------------------------------------
## No-correlation 
#-------------------------------------------------------------------------------

old.par = par(mar = c(5,5,4,1))
plot(variogramLine(vgm(psill = 0, "Nug",  0, nugget = 1), covariance = FALSE, maxdist = 4), 
     ylim = c(0, 1), xlab = "r", ylab = "Semivariogram", main = "No-correlation", type = "l", cex.main=2, lwd = 3, cex.lab = 2, cex.axis = 1.5)
points(0, 0, pch = 19, cex = 2)
legend("bottomright", cex = 1.5, lwd = 2,
  legend=c(expression(paste(sigma^2, " = 1", ""))),lty=1)

#-------------------------------------------------------------------------------
## Exponential
#-------------------------------------------------------------------------------

plot(variogramLine(vgm(psill = 1, "Exp",  range = 1/3, nugget = 0), 
  covariance = FALSE, maxdist = 4), ylim = c(0, 1), xlab = "r", lwd = 3,
  ylab = "Semivariogram", main = "Exponential", type = "l", cex.main=2, cex.lab =2,
  cex.axis = 1.5,)
lines(variogramLine(vgm(psill = 1, "Exp",  range = 2/3, nugget = 0), 
  covariance = FALSE, maxdist = 4), type = "l", lty = 2, lwd = 3)
lines(variogramLine(vgm(psill = 1, "Exp",  range = 1, nugget = 0), 
  covariance = FALSE, maxdist = 4), type = "l", lty = 3, lwd = 3)
legend("bottomright", cex = 1.5, lwd = 2,
  legend=c(expression(paste(sigma^2, " = 1, ", alpha, " = 1/3")), 
  expression(paste(sigma^2, " = 1, ", alpha, " = 2/3")), 
  expression(paste(sigma^2, " = 1, ", alpha, " = 1"))),lty=1:3)

#-------------------------------------------------------------------------------
## Linear
#-------------------------------------------------------------------------------

# linear variogram
vgmline = variogramLine(vgm(psill = 1, "Pow", 1, nugget = 0), 
  covariance = FALSE, maxdist = 4)
# scale variogram by 0.5 
 
plot(vgmline[,1], 0.5*vgmline[,2], ylim = c(0, 6), xlab = "r", 
  ylab = "Semivariogram", main = "Linear", type = "l", cex.main=2, 
  lwd = 3, cex.lab = 2, cex.axis = 1.5)
lines(vgmline[,1], vgmline[,2], type = "l", lty = 2, lwd = 3)
lines(vgmline[,1], 1.5*vgmline[,2], type = "l", lty = 3, lwd = 3)
legend("topleft", cex = 1.5, lwd = 2,
  legend=c(expression(paste(theta[1], " = 0.5", "")), 
  expression(paste(theta[1], " = 1", "")), 
  expression(paste(theta[1], " = 1.5", ""))),lty=1:3)

#-------------------------------------------------------------------------------
## Power
#-------------------------------------------------------------------------------

plot(variogramLine(vgm(psill = 1, "Pow",  range = 0.5, nugget = 0), 
  covariance = FALSE, maxdist = 4), ylim = c(0, 8), xlab = "r", lwd = 3,
  ylab = "Semivariogram", main = "Power", type = "l", cex.main=2, cex.lab =2,
  cex.axis = 1.5,)
lines(variogramLine(vgm(psill = 1, "Pow",  range = 1, nugget = 0), 
  covariance = FALSE, maxdist = 4), type = "l", lty = 2, lwd = 3)
lines(variogramLine(vgm(psill = 1, "Pow",  range = 1.5, nugget = 0), 
  covariance = FALSE, maxdist = 4), type = "l", lty = 3, lwd = 3)
legend("topleft", cex = 1.5, lwd = 2,
  legend=c(expression(paste(theta[1], " = 1, ", theta[2], " = 0.5")), 
  expression(paste(theta[1], " = 1, ", theta[2], " = 1")), 
  expression(paste(theta[1], " = 1, ", theta[2], " = 1.5"))),lty=1:3)

par(old.par)
dev.off()


system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

