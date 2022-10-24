sec_path = 'Rcode/Chapter6/Section 6.2/figures/'
setwd(paste0(SLEDbook_path,sec_path))

library(gstat)
library(viridis)
library(classInt)
source('addBreakColorLegend.R')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                Single Parameter Covariogram Figure
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "Covariogram1"

pdf(paste0(file_name,'.pdf'), width = 8.5, height = 11)

###### Covariogram1 (Figure 6.1) ###### 

layout(matrix(1:8, ncol = 2, byrow = TRUE))

#-------------------------------------------------------------------------------
## No-correlation #-------------------------------------------------------------------------------

old.par = par(mar = c(5,5,4,1))
plot(variogramLine(vgm(psill = 1, "Nug",  0, nugget = 0), covariance = TRUE, maxdist = 4), 
     ylim = c(0, 1), xlab = "r", ylab = "Covariance", main = "No-correlation", type = "l", cex.main=2, lwd = 3, cex.lab = 2, cex.axis = 1.5)
points(0, 1, pch = 19, cex = 2)

#-------------------------------------------------------------------------------
## Triangular (Linear) #-------------------------------------------------------------------------------

plot(variogramLine(vgm(psill = 1, "Lin",  range = 1, nugget = 0), 
  covariance = TRUE, maxdist = 4), ylim = c(0, 1), xlab = "r", lwd = 3,
  ylab = "Covariance", main = "Triangular", type = "l", cex.main=2, cex.lab =2,
  cex.axis = 1.5,)
lines(variogramLine(vgm(psill = 1, "Lin",  range = 2, nugget = 0), 
  covariance = TRUE, maxdist = 4), type = "l", lty = 2, lwd = 3)
lines(variogramLine(vgm(psill = 1, "Lin",  range = 3, nugget = 0), 
  covariance = TRUE, maxdist = 4), type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2,
  legend=c(expression(paste(alpha, " = 1", "")), 
  expression(paste(alpha, " = 2", "")), 
  expression(paste(alpha, " = 3", ""))),lty=1:3)

#-------------------------------------------------------------------------------
## Circular 
#-------------------------------------------------------------------------------

plot(variogramLine(vgm(psill = 1, "Cir",  range = 1, nugget = 0),
  covariance = TRUE, maxdist = 4), ylim = c(0, 1), xlab = "r", lwd = 3, 
  ylab = "Covariance", main = "Circular", type = "l", cex.main=2, cex.lab = 2,
  cex.axis = 1.5)
lines(variogramLine(vgm(psill = 1, "Cir",  range = 2, nugget = 0), 
  covariance = TRUE, maxdist = 4), lwd = 3,
  type = "l", lty = 2)
lines(variogramLine(vgm(psill = 1, "Cir",  range = 3, nugget = 0), 
  covariance = TRUE, maxdist = 4), lwd = 3,
  type = "l", lty = 3)
legend("topright", cex = 1.5, lwd = 2,
  legend=c(expression(paste(alpha, " = 1", "")), 
  expression(paste(alpha, " = 2", "")), 
  expression(paste(alpha, " = 3", ""))),lty=1:3)

#-------------------------------------------------------------------------------
## Spherical
#-------------------------------------------------------------------------------

plot(variogramLine(vgm(psill = 1, "Sph",  range = 1, nugget = 0),
  covariance = TRUE, maxdist = 4), ylim = c(0, 1), xlab = "r", lwd = 3, 
  ylab = "Covariance", main = "Spherical", type = "l", cex.main=2, cex.lab = 2,
  cex.axis = 1.5)
lines(variogramLine(vgm(psill = 1, "Sph",  range = 2, nugget = 0), 
  covariance = TRUE, maxdist = 4), lwd = 3,
  type = "l", lty = 2)
lines(variogramLine(vgm(psill = 1, "Sph",  range = 3, nugget = 0), 
  covariance = TRUE, maxdist = 4), lwd = 3,
  type = "l", lty = 3)
legend("topright", cex = 1.5, lwd = 2,
  legend=c(expression(paste(alpha, " = 1", "")), 
  expression(paste(alpha, " = 2", "")), 
  expression(paste(alpha, " = 3", ""))),lty=1:3)

#-------------------------------------------------------------------------------
## Cosine
#-------------------------------------------------------------------------------

my_cos_fnc <- function(r, theta2){
    cos(r / theta2)
}
x <- seq(0, 4, 0.01)
plot(x, my_cos_fnc(x, theta2 = 0.2), xlab = "r", ylab = "Covariance", 
  main = "Cosine", type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3)
lines(x, my_cos_fnc(x, theta2 = 0.5),type = "l", lty = 2, lwd = 3)
lines(x, my_cos_fnc(x, theta2 = 1),type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2,bg = 'white',
  legend=c(expression(paste(alpha, " = 0.2", "")), 
  expression(paste(alpha, " = 0.5", "")), 
  expression(paste(alpha, " = 1", ""))),lty=1:3)

#-------------------------------------------------------------------------------
## Wave
#-------------------------------------------------------------------------------

my_wave_fnc <- function(r, theta2){
    theta2 * sin(r / theta2) / r
}
x <- seq(0, 4, 0.01)
plot(x, my_wave_fnc(x, theta2 = 0.1), xlab = "r", ylab = "Covariance", 
  main = "Wave", type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3)
lines(x, my_wave_fnc(x, theta2 = 0.2),type = "l", lty = 2, lwd = 3)
lines(x, my_wave_fnc(x, theta2 = 0.5),type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2,bg = 'white',
  legend=c(expression(paste(alpha, " = 0.1", "")), 
  expression(paste(alpha, " = 0.2", "")), 
  expression(paste(alpha, " = 0.5", ""))),lty=1:3)

#-------------------------------------------------------------------------------
## Exponential
#-------------------------------------------------------------------------------

plot(variogramLine(vgm(psill = 1, "Exp",  range = 1/3, nugget = 0),
  covariance = TRUE, maxdist = 4), ylim = c(0, 1), xlab = "r", lwd = 3, 
  ylab = "Covariance", main = "Exponential", type = "l", cex.main=2, cex.lab = 2,
  cex.axis = 1.5)
lines(variogramLine(vgm(psill = 1, "Exp",  range = 2/3, nugget = 0), 
  covariance = TRUE, maxdist = 4), lwd = 3,
  type = "l", lty = 2)
lines(variogramLine(vgm(psill = 1, "Exp",  range = 1, nugget = 0), 
  covariance = TRUE, maxdist = 4), lwd = 3,
  type = "l", lty = 3)
legend("topright", cex = 1.5, lwd = 2,
  legend=c(expression(paste(alpha, " = 1/3", "")), 
  expression(paste(alpha, " = 2/3", "")), 
  expression(paste(alpha, " = 1", ""))),lty=1:3)

#-------------------------------------------------------------------------------
## Gaussian
#-------------------------------------------------------------------------------

plot(variogramLine(vgm(psill = 1, "Gau",  range = 1/sqrt(3), nugget = 0),
  covariance = TRUE, maxdist = 4), ylim = c(0, 1), xlab = "r", lwd = 3, 
  ylab = "Covariance", main = "Gaussian", type = "l", cex.main=2, cex.lab = 2,
  cex.axis = 1.5)
lines(variogramLine(vgm(psill = 1, "Gau",  range = 2/sqrt(3), nugget = 0), 
  covariance = TRUE, maxdist = 4), lwd = 3,
  type = "l", lty = 2)
lines(variogramLine(vgm(psill = 1, "Gau",  range = 3/sqrt(3), nugget = 0), 
  covariance = TRUE, maxdist = 4), lwd = 3,
  type = "l", lty = 3)
legend("topright", cex = 1.5, lwd = 2,
  legend=c(expression(paste(alpha, " = ", 1/sqrt(3), "")), 
  expression(paste(alpha, " = ", 2/sqrt(3), "")), 
  expression(paste(alpha, " = ", 3/sqrt(3), ""))),lty=1:3)

par(old.par)
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#      2-Parameter Covariogram Figures Cachy, Powed Exponential, Matern
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "Covariogram2"

pdf(paste0(file_name,'.pdf'), width = 8.5, height = 11)

layout(matrix(1:6, ncol = 2, byrow = TRUE))
old.par = par(mar = c(5,5,4,1))

#-------------------------------------------------------------------------------
## Cauchy
#-------------------------------------------------------------------------------

my_cau_fnc <- function(r, theta2, theta3){
    (1 + (r / theta2)^2)^{-theta3}
}
x <- seq(0, 4, 0.01)

plot(x, my_cau_fnc(x, theta2 = 1/2, theta3 = 1), xlab = "r", ylab = "Covariance", 
  main = "Cauchy", type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3)
lines(x, my_cau_fnc(x, theta2 = 1, theta3 = 1),type = "l", lty = 2, lwd = 3)
lines(x, my_cau_fnc(x, theta2 = 2, theta3 = 1),type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2, bg = 'white',
  legend=c(expression(paste(alpha, " = 0.5, ", phi, " = 1.0", "")), 
          expression(paste(alpha, " = 1.0, ", phi, " = 1.0", "")), 
          expression(paste(alpha, " = 2.0, ", phi, " = 1.0", ""))),
          lty=1:3)
          
plot(x, my_cau_fnc(x, theta2 = 1, theta3 = 0.5), xlab = "r", ylab = "Covariance", 
  main = "Cauchy", type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1))
lines(x, my_cau_fnc(x, theta2 = 1, theta3 = 1),type = "l", lty = 2, lwd = 3)
lines(x, my_cau_fnc(x, theta2 = 1, theta3 = 1.5),type = "l", lty = 3, lwd = 3)
lines(x, my_cau_fnc(x, theta2 = 1, theta3 = 2),type = "l", lty = 4, lwd = 3)
legend("topright", cex = 1.5, lwd = 1.5, bg = 'white',
  legend=c(expression(paste(alpha, " = 1.0, ", phi, " = 0.5", "")), 
          expression(paste(alpha, " = 1.0, ", phi, " = 1.0", "")), 
          expression(paste(alpha, " = 1.0, ", phi, " = 1.5", "")),
          expression(paste(alpha, " = 1.0, ", phi, " = 2.0", ""))),
          lty=1:4)

#-------------------------------------------------------------------------------
## Powered Exponential
#-------------------------------------------------------------------------------

my_pe_fnc <- function(r, theta2, theta3){
    exp(- (r^theta3) / theta2)
}

plot(x, my_pe_fnc(x, theta2 = 1, theta3 = 1.5), xlab = "r", ylab = "Covariance", 
  main = "Powered Exponential", type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1))
lines(x, my_pe_fnc(x, theta2 = 2, theta3 = 1.5),type = "l", lty = 2, lwd = 3)
lines(x, my_pe_fnc(x, theta2 = 3, theta3 = 1.5),type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2, bg = 'white',
  legend=c(expression(paste(alpha, " = 1.0, ", theta, " = 1.5", "")), 
          expression(paste(alpha, " = 2.0, ", theta, " = 1.5", "")), 
          expression(paste(alpha, " = 3.0, ", theta, " = 1.5", ""))),
          lty=1:3)

plot(x, my_pe_fnc(x, theta2 = 2, theta3 = 0.5), xlab = "r", ylab = "Covariance", 
  main = "Powered Exponential", type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1))
lines(x, my_pe_fnc(x, theta2 = 2, theta3 = 1.0),type = "l", lty = 2, lwd = 3)
lines(x, my_pe_fnc(x, theta2 = 2, theta3 = 1.5),type = "l", lty = 3, lwd = 3)
lines(x, my_pe_fnc(x, theta2 = 2, theta3 = 2.0),type = "l", lty = 4, lwd = 3)
legend("topright", cex = 1.5, lwd = 1.5, bg = 'white',
  legend=c(expression(paste(alpha, " = 2.0, ", theta, " = 0.5", "")), 
          expression(paste(alpha, " = 2.0, ", theta, " = 1.0", "")), 
          expression(paste(alpha, " = 2.0, ", theta, " = 1.5", "")),
          expression(paste(alpha, " = 2.0, ", theta, " = 2.0", ""))),
          lty=1:4)

#-------------------------------------------------------------------------------
## Matern
#-------------------------------------------------------------------------------

my_mat_fnc_2 <- function(r, theta2, theta3){
    ifelse(r > 0, 2^{1 - theta3} / gamma(theta3) * ((sqrt(2 * theta3) * r/ theta2) ^ {theta3}) * besselK(x = sqrt(2 * theta3) * r / theta2, theta3), 1)
}


plot(x, my_mat_fnc_2(x, theta2 = 1, theta3 = 1.5), xlab = "r", 
	ylab = "Covariance", 
  main = "Mat\uE9rn", type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1))
lines(x, my_mat_fnc_2(x, theta2 = 2, theta3 = 1.5),type = "l", lty = 2, lwd = 3)
lines(x, my_mat_fnc_2(x, theta2 = 3, theta3 = 1.5),type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2, bg = 'white',
  legend=c(expression(paste(alpha, " = 1.0, ", nu, " = 1.5", "")), 
          expression(paste(alpha, " = 2.0, ", nu, " = 1.5", "")), 
          expression(paste(alpha, " = 3.0, ", nu, " = 1.5", ""))),
          lty=1:3)

plot(x, my_mat_fnc_2(x, theta2 = 2, theta3 = 0.5), xlab = "r", ylab = "Covariance", 
  main = "Mat\uE9rn", type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1))
lines(x, my_mat_fnc_2(x, theta2 = 2, theta3 = 1.5),type = "l", lty = 2, lwd = 3)
lines(x, my_mat_fnc_2(x, theta2 = 2, theta3 = 2.5),type = "l", lty = 3, lwd = 3)
lines(variogramLine(vgm(psill = 1, "Gau",  range = 2 * sqrt(2), nugget = 0), 
  covariance = TRUE, maxdist = 4), type = "l", lty = 4, lwd = 3)
legend("topright", cex = 1.5, lwd = 1.5, bg = 'white',
  legend=c(expression(paste(alpha, " = 2.0, ", nu, " = 0.5", "")), 
          expression(paste(alpha, " = 2.0, ", nu, " = 1.5", "")), 
          expression(paste(alpha, " = 2.0, ", nu, " = 2.5", "")),
          expression(paste(alpha, " = 2.0, ", nu, " = ", infinity, ""))),
          lty=1:4)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#      2-Parameter Covariogram Figures Askey and Wendland
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "Covariogram3"

pdf(paste0(file_name,'.pdf'), width = 8.5, height = 11)

layout(matrix(1:6, ncol = 2, byrow = TRUE))
old.par = par(mar = c(5,5,4,1))

x <- seq(0, 4, 0.01)

#-------------------------------------------------------------------------------
## Askey
#-------------------------------------------------------------------------------

my_ask_fnc <- function(r, alpha, mu){
	if(r/alpha < 1) {(1 - r/alpha)^mu} else { 0 }
}

fit = x*0
for(i in 1:length(x)) fit[i] = my_ask_fnc(x[i], alpha = 1, mu = 1.5)
plot(x, fit , xlab = "r", 
	ylab = "Covariance", 
  main = "Askey", type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1))
fit = x*0
for(i in 1:length(x)) fit[i] = my_ask_fnc(x[i], alpha = 2, mu = 1.5)
lines(x, fit, type = "l", lty = 2, lwd = 3)
fit = x*0
for(i in 1:length(x)) fit[i] = my_ask_fnc(x[i], alpha = 3, mu = 1.5)
lines(x, fit, type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2, bg = 'white',
  legend=c(expression(paste(alpha, " = 1.0, ", mu, " = 1.5", "")), 
          expression(paste(alpha, " = 2.0, ", mu, " = 1.5", "")), 
          expression(paste(alpha, " = 3.0, ", mu, " = 1.5", ""))),
          lty=1:3)

fit = x*0
for(i in 1:length(x)) fit[i] = my_ask_fnc(x[i], alpha = 3, mu = 1.5)
plot(x, fit , xlab = "r", 
	ylab = "Covariance", 
  main = "Askey", type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1))
fit = x*0
for(i in 1:length(x)) fit[i] = my_ask_fnc(x[i], alpha = 3, mu = 2.5)
lines(x, fit, type = "l", lty = 2, lwd = 3)
fit = x*0
for(i in 1:length(x)) fit[i] = my_ask_fnc(x[i], alpha = 3, mu = 3.5)
lines(x, fit, type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2, bg = 'white',
  legend=c(expression(paste(alpha, " = 3.0, ", mu, " = 1.5", "")), 
          expression(paste(alpha, " = 3.0, ", mu, " = 2.5", "")), 
          expression(paste(alpha, " = 3.0, ", mu, " = 3.5", ""))),
          lty=1:3)

#-------------------------------------------------------------------------------
## C2-Wendland
#-------------------------------------------------------------------------------

my_C2W_fnc <- function(r, alpha, mu){
	if(r/alpha < 1) {(1 + (mu + 1)*r/alpha)*(1 - r/alpha)^(mu + 1)} else { 0 }
}

fit = x*0
for(i in 1:length(x)) fit[i] = my_C2W_fnc(x[i], alpha = 1, mu = 2.5)
plot(x, fit , xlab = "r", 
	ylab = "Covariance", 
  main = expression(bold(paste(C^2,"-Wendland"))),
  type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1.0))
fit = x*0
for(i in 1:length(x)) fit[i] = my_C2W_fnc(x[i], alpha = 2, mu = 2.5)
lines(x, fit, type = "l", lty = 2, lwd = 3)
fit = x*0
for(i in 1:length(x)) fit[i] = my_C2W_fnc(x[i], alpha = 3, mu = 1.5)
lines(x, fit, type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2, bg = 'white',
  legend=c(expression(paste(alpha, " = 1.0, ", mu, " = 2.5", "")), 
          expression(paste(alpha, " = 2.0, ", mu, " = 2.5", "")), 
          expression(paste(alpha, " = 3.0, ", mu, " = 2.5", ""))),
          lty=1:3)

fit = x*0
for(i in 1:length(x)) fit[i] = my_C2W_fnc(x[i], alpha = 3, mu = 2.5)
plot(x, fit , xlab = "r", 
	ylab = "Covariance", 
  main = expression(bold(paste(C^2,"-Wendland"))),
  type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1.0))
fit = x*0
for(i in 1:length(x)) fit[i] = my_C2W_fnc(x[i], alpha = 3, mu = 5)
lines(x, fit, type = "l", lty = 2, lwd = 3)
fit = x*0
for(i in 1:length(x)) fit[i] = my_C2W_fnc(x[i], alpha = 3, mu = 10)
lines(x, fit, type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2, bg = 'white',
  legend=c(expression(paste(alpha, " = 3.0, ", mu, " = 2.5", "")), 
          expression(paste(alpha, " = 3.0, ", mu, " = 5", "")), 
          expression(paste(alpha, " = 3.0, ", mu, " = 10", ""))),
          lty=1:3)
          
#-------------------------------------------------------------------------------
## C4-Wendland
#-------------------------------------------------------------------------------

my_C4W_fnc <- function(r, alpha, mu){
	if(r/alpha < 1) {((1 + (mu + 2)*r/alpha) +
		((mu + 2)^2 - 1)/3*(r/alpha)^2)*(1 - r/alpha)^(mu + 2)} else { 0 }
}

fit = x*0
for(i in 1:length(x)) fit[i] = my_C4W_fnc(x[i], alpha = 1, mu = 3.5)
plot(x, fit , xlab = "r", 
	ylab = "Covariance", 
  main = expression(bold(paste(C^4,"-Wendland"))),
  type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1.0))
fit = x*0
for(i in 1:length(x)) fit[i] = my_C4W_fnc(x[i], alpha = 2, mu = 3.5)
lines(x, fit, type = "l", lty = 2, lwd = 3)
fit = x*0
for(i in 1:length(x)) fit[i] = my_C4W_fnc(x[i], alpha = 3, mu = 1.5)
lines(x, fit, type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2, bg = 'white',
  legend=c(expression(paste(alpha, " = 1.0, ", mu, " = 3.5", "")), 
          expression(paste(alpha, " = 2.0, ", mu, " = 3.5", "")), 
          expression(paste(alpha, " = 3.0, ", mu, " = 3.5", ""))),
          lty=1:3)

fit = x*0
for(i in 1:length(x)) fit[i] = my_C2W_fnc(x[i], alpha = 3, mu = 3.5)
plot(x, fit , xlab = "r", 
	ylab = "Covariance", 
  main = expression(bold(paste(C^4,"-Wendland"))),
  type = "l", cex.main=2, cex.lab = 2, cex.axis = 1.5, 
  lwd = 3, ylim = c(0,1.0))
fit = x*0
for(i in 1:length(x)) fit[i] = my_C2W_fnc(x[i], alpha = 3, mu = 5)
lines(x, fit, type = "l", lty = 2, lwd = 3)
fit = x*0
for(i in 1:length(x)) fit[i] = my_C2W_fnc(x[i], alpha = 3, mu = 10)
lines(x, fit, type = "l", lty = 3, lwd = 3)
legend("topright", cex = 1.5, lwd = 2, bg = 'white',
  legend=c(expression(paste(alpha, " = 3.0, ", mu, " = 3.5", "")), 
          expression(paste(alpha, " = 3.0, ", mu, " = 5", "")), 
          expression(paste(alpha, " = 3.0, ", mu, " = 10", ""))),
          lty=1:3)

layout(1)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))


################################################################################
#-------------------------------------------------------------------------------
#              Simulation Plot 1 
#-------------------------------------------------------------------------------
################################################################################

file_name = "SimPlot1"

pdf(paste0(file_name,'.pdf'), width = 14, height = 8)
   
cex_main = 2.5
 
source('addBreakColorLegend.R')

layout(matrix(1:12, ncol = 6, byrow = TRUE), widths = c(4,1,4,1,4,1))

# create a grid of points
xy <- expand.grid(seq(0.1, 10, length.out = 100), seq(0.1, 10, length.out = 100))
names(xy) <- c("x", "y")

#-------------------------------------------------------------------------------
## No correlation
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(1, "Sph", range=0.01), nmax = 20)
set.seed(150)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("No-correlation"))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Spherical 1
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(1, "Sph", range=1.0), nmax = 20)
set.seed(151)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Spherical (", alpha, " = 1)", ""))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Wave 1
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(1, "Wav", 0.1 * pi), nmax = 20)
set.seed(152)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Wave (", alpha, " = 0.1)", ""))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Spherical plus nugget
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(psill=0.75, "Sph", range=2.0, nugget=0.25), nmax = 20)
set.seed(153)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Spherical (", alpha, " = 2, ", 
		sigma[0]^2, " = 0.25)"))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2.0)

#-------------------------------------------------------------------------------
## Spherical 2
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(1, "Sph", range=2.0), nmax = 20)
set.seed(154)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Spherical (", alpha, " = 2)", ""))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Wave 2
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(1, "Wav", 0.5 * pi), nmax = 20)
set.seed(155)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Wave (", alpha, " = 0.5)", ""))),
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






################################################################################
#-------------------------------------------------------------------------------
#              Simulation Plot 2 
#-------------------------------------------------------------------------------
################################################################################

file_name = "SimPlot2"

pdf(paste0(file_name,'.pdf'), width = 14, height = 8)
    
cex_main = 2.5

source('addBreakColorLegend.R')

layout(matrix(1:12, ncol = 6, byrow = TRUE), widths = c(4,1,4,1,4,1))

# create a grid of points
xy <- expand.grid(seq(0.1, 10, length.out = 100), seq(0.1, 10, length.out = 100))
names(xy) <- c("x", "y")

#-------------------------------------------------------------------------------
## Exponential 1
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(psill = 1, "Exp",  range = 1/3, nugget = 0), nmax = 20)
set.seed(161)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Exponential (", alpha, " = 1/3)", ""))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Gaussian 1
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(psill = 1, "Gau",  range = 1/sqrt(3), nugget = 0), nmax = 20)
set.seed(162)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Gaussian (", alpha, " = 1/",sqrt(3),")"))), 
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Matern 1
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
  model = vgm(psill = 1, "Mat",  range = 1, nugget = 0, kappa = 1.5), 
  nmax = 20)
set.seed(163)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Mat\uE9rn (", alpha, " = 1, ", 
		nu, " = 1.5)", ""))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Exponential 2
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(psill = 1, "Exp",  range = 1, nugget = 0), nmax = 20)
set.seed(164)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Exponential (", alpha, " = 1)", ""))),
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Gaussian 2
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                   model = vgm(psill = 1, "Gau",  range = 3/sqrt(3), nugget = 0), nmax = 20)
set.seed(165)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Gaussian (", alpha, " = 3/",sqrt(3),")"))), 
  breaks = cip$brks, cex.main = cex_main,
  xlim = c(0.05,10.05), ylim = c(0.05,10.05))
par(mar = c(0,0,0,0))
plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
addBreakColorLegend(xleft = 0.02, ybottom = .05, xright = .3, ytop = .85,
    breaks = cip$brks, colors = palp, cex = 2)

#-------------------------------------------------------------------------------
## Matern 2
#-------------------------------------------------------------------------------

g.dummy.1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
  model = vgm(psill = 1, "Mat",  range = 1, nugget = 0, kappa = 2.5), 
  nmax = 20)
set.seed(166)
simdata1 <- predict(g.dummy.1, newdata = xy, nsim = 1)
old.par = par(mar = c(1,2,3,0))
cip = classIntervals(simdata1$sim1, n = 12, style = 'fisher')
palp = viridis(12)
image(simdata1,col = palp, bty = 'n', xaxt = 'n', yaxt = 'n', 
  main = expression(bold(paste("Mat\uE9rn (", alpha, " = 1, ", 
		nu, " = 2.5)", ""))),
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

