library(gstat)
library(sp)
library(splancs) #Dependency of {geoR}
library(RandomFieldsUtils)
library(RandomFields) #Dependency of {geoR}, no longer exists on CRAN
library(geoR) #No longer exists on CRAN
library(spmodel)

setwd("C:/Users/rainb/Desktop/20220701")
RNGkind(sample.kind = "Rounding")
xy <- cbind(c(0.1, 0.15, 0.30, 0.32, 0.35, 0.40, 0.50, 0.60, 0.75, 0.90), c(rep(0, 10)))
xy <- as.data.frame(xy)
names(xy) <- c("x", "y")
new.xy <- as.data.frame(cbind(c(seq(0, 1, by = 0.001)), c(rep(0, 1001))))
names(new.xy) <- c("x", "y")

########Figure 9.2########----------------------------------------------------------------------------
g.dummy.exp1 <- gstat(formula = z ~ 1, locations = ~x+y, dummy = T, beta = 0,
                      model = vgm(1, "Exp", 0.05), nmax = 20)
set.seed(10)
simdata.exp1 <- predict(g.dummy.exp1, newdata = xy, nsim = 4)
# simulate data
set.seed(13)
spcov_params_val <- spcov_params("exponential", de = 1, ie = 0.0001, 
	range = 0.05)
simdata_exp1 = sprnorm(spcov_params_val, data = xy, xcoord = x, ycoord = y)

D1 = data.frame(ysim = simdata_exp1, x = xy$x, y = xy$y)
# create spmodel with fixed covariance parameters
splm_exp = splm(ysim ~ 1, data = D1, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
	de = 1, ie = 0.00001, range = 0.05, known = 'given'))
splm_sph = splm(ysim ~ 1, data = D1, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'spherical', 
	de = 1, ie = 0.00001, range = 0.15, known = 'given'))
splm_gau = splm(ysim ~ 1, data = D1, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'gaussian', 
	de = 1, ie = 0.00001, range = 0.15/sqrt(3), known = 'given'))
predexp = predict(splm_exp, new.xy, se.fit = TRUE, interval = 'prediction')
predsph = predict(splm_sph, new.xy, se.fit = TRUE, interval = 'prediction')
predgau = predict(splm_gau, new.xy, se.fit = TRUE, interval = 'prediction')

layout(matrix(1:4, nrow = 2, byrow = TRUE))
par(mar = c(5,5,1,1))
plot(D1$x, D1$ysim, pch = 19, ylim = c(-3, 3), cex = 2,
	xlab = 'Location Coordinate', ylab = 'Response Value', 
  cex.lab = 2, cex.axis = 1.5)

plot(new.xy[,1], predexp$fit[,1], type = 'l', ylim = c(-3,4),
  xlab = 'Location Coordinate', ylab = 'Response Value', 
  cex.lab = 2, cex.axis = 1.5, lwd = 2)
points(D1$x, D1$ysim, pch = 19, cex = 2)
lines(new.xy[,1], predexp$se.fit^2, lty = 2, col = 'red', lwd = 2)
legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
       lty = c(1, 2), cex = 1.5)

plot(new.xy[,1], predsph$fit[,1], type = 'l', ylim = c(-3,4),
  xlab = 'Location Coordinate', ylab = 'Response Value', 
  cex.lab = 2, cex.axis = 1.5, lwd = 2)
points(D1$x, D1$ysim, pch = 19, cex = 2)
lines(new.xy[,1], predsph$se.fit^2, lty = 2, col = 'red', lwd = 2)
legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
       lty = c(1, 2), cex = 1.5)

plot(new.xy[,1], predgau$fit[,1], type = 'l', ylim = c(-3,4),
  xlab = 'Location Coordinate', ylab = 'Response Value', 
  cex.lab = 2, cex.axis = 1.5, lwd = 2)
points(D1$x, D1$ysim, pch = 19, cex = 2)
lines(new.xy[,1], predgau$se.fit^2, lty = 2, col = 'red', lwd = 2)
legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
       lty = c(1, 2), cex = 1.5)

layout(1)

## prediction
KC_exp1 <- krige.control(type.krige = "ok", cov.model = "exponential", 
	cov.pars = c(1, 0.05))
pred.exp1.4 <- krige.conv(coords = xy, data = simdata.exp1$sim4, 
	locations = new.xy, borders = NULL, krige = KC_exp1)
var.exp1.4 <- krige.conv(coords = xy, data = simdata.exp1$sim4, 
	locations = new.xy, borders = NULL, krige = KC_exp1)$krige.var


## first panel
png("plot1_1.png", width = 600, height = 300)
plot(xy$x, simdata.exp1$sim4, pch = 4, xlab = "", ylab = "", ylim = c(-2.5, 2.5))
dev.off()

## second panel
png("plot1_2.png", width = 600, height = 300)
plot(new.xy[,1], pred.exp1.4$predict, type = "l", xlab = "", ylab = "", ylim = c(-2.5, 2.5))
points(xy$x, simdata.exp1$sim4, pch = 4)
lines(new.xy[,1], var.exp1.4, lty = 2, col = "red")
legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
       lty = c(1, 2))
dev.off()

## third panel
KC_sph1 <- krige.control(type.krige = "ok", cov.model = "spherical", cov.pars = c(1, 0.15))
pred.sph1.4 <- krige.conv(coords = xy, data = simdata.exp1$sim4, locations = new.xy, borders = NULL,
                          krige = KC_sph1)
var.sph1.4 <- krige.conv(coords = xy, data = simdata.exp1$sim4, locations = new.xy, borders = NULL,
                         krige = KC_sph1)$krige.var

png("plot1_3.png", width = 600, height = 300)
plot(new.xy[,1], pred.sph1.4$predict, type = "l", xlab = "", ylab = "", ylim = c(-2.5, 2.5))
points(xy$x, simdata.exp1$sim4, pch = 4)
lines(new.xy[,1], var.sph1.4, lty = 2, col = "red")
legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
       lty = c(1, 2))
dev.off()


## fourth panel
KC_gau1 <- krige.control(type.krige = "ok", cov.model = "gaussian", cov.pars = c(1, 0.15 / sqrt(3)))
pred.gau1.4 <- krige.conv(coords = xy, data = simdata.exp1$sim4, locations = new.xy, borders = NULL,
                          krige = KC_gau1)
var.gau1.4 <- krige.conv(coords = xy, data = simdata.exp1$sim4, locations = new.xy, borders = NULL,
                         krige = KC_gau1)$krige.var

png("plot1_4.png", width = 600, height = 300)
plot(new.xy[,1], pred.gau1.4$predict, type = "l", xlab = "", ylab = "", ylim = c(-2.5, 2.5))
points(xy$x, simdata.exp1$sim4, pch = 4)
lines(new.xy[,1], var.gau1.4, lty = 2, col = "red")
legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
       lty = c(1, 2))
dev.off()

########Figure 9.3########-----------------------------------------------------------------------------
##second panel (universal kriging)
KC_exp_uk <- krige.control(type.krige = "ok",  trend.d = ~1 + xy$x, trend.l = ~ 1 + new.xy$x,
                           cov.model = "exponential", cov.pars = c(1, 0.05))
pred.exp.uk<- krige.conv(coords = xy, data = simdata.exp1$sim4, locations = new.xy, borders = NULL,
                          krige = KC_exp_uk)
var.exp.uk <- krige.conv(coords = xy, data = simdata.exp1$sim4, locations = new.xy, borders = NULL,
                         krige = KC_exp_uk)$krige.var

png("plot2_2.png", width = 600, height = 300)
plot(new.xy$x, pred.exp.uk$predict, type = "l", xlab = "", ylab = "", ylim = c(-2.5, 2.5))
points(xy$x, simdata.exp1$sim4, pch = 4)
lines(new.xy$x, var.exp.uk, lty = 2, col = "red")
legend("topright", legend = c("UK", "UK var"), col = c("black", "red"),
       lty = c(1, 2))
dev.off()

##third panel
KC_exp4 <- krige.control(type.krige = "ok", cov.model = "exponential", cov.pars = c(0.5, 0.05), nugget = 0.5)
pred.exp4.4 <- krige.conv(coords = xy, data = simdata.exp1$sim4, locations = new.xy, borders = NULL,
                          krige = KC_exp4)
var.exp4.4 <- krige.conv(coords = xy, data = simdata.exp1$sim4, locations = new.xy, borders = NULL,
                         krige = KC_exp4)$krige.var
png("plot2_3.png", width = 600, height = 300)
plot(new.xy$x, pred.exp4.4$predict, type = "l", xlab = "", ylab = "", ylim = c(-2.5, 2.5))
points(xy$x, simdata.exp1$sim4, pch = 4)
lines(new.xy$x, var.exp4.4, lty = 2, col = "red")
legend("topright", legend = c("OK", "OK var"), col = c("black", "red"),
       lty = c(1, 2))
dev.off()

##fourth panel
Dist1  <- matrix(c(0, 0.05, 0.20, 0.22, 0.25, 0.30, 0.40, 0.50, 0.65, 0.80, 0,
                   0,    0, 0.15, 0.17, 0.20, 0.25, 0.35, 0.45, 0.60, 0.75, 0,
                   0,    0,    0, 0.02, 0.05, 0.10, 0.20, 0.30, 0.45, 0.60, 0,
                   0,    0,    0,    0, 0.03, 0.08, 0.18, 0.28, 0.43, 0.58, 0,
                   0,    0,    0,    0,    0, 0.05, 0.15, 0.25, 0.40, 0.55, 0,
                   0,    0,    0,    0,    0,    0, 0.10, 0.20, 0.35, 0.50, 0,
                   0,    0,    0,    0,    0,    0,    0, 0.10, 0.25, 0.40, 0,
                   0,    0,    0,    0,    0,    0,    0,    0, 0.15, 0.30, 0,
                   0,    0,    0,    0,    0,    0,    0,    0,    0, 0.15, 0,
                   0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 0,
                   0,    0,    0,    0,    0,    0,    0,    0,    0,    0, 0),
                 ncol = 11, byrow = TRUE)
Dist0 <- t(Dist1) + Dist1

nl.out <- function(ind, model, sigma.ep2, obs){
  Gamma0 <- apply(Dist0, c(1, 2), model)
  Gamma0[1 : 10, 11] <- rep(1, 10)
  Gamma0[11, 1 : 10] <- rep(1, 10)
  diag(Gamma0) <- rep(0, 11)
  dist0 <- abs(xy[, 1] - xy[ind, 1])
  gamma0 <- c(model(dist0), 1)
  gamma0[ind] <- gamma0[ind] + sigma.ep2
  lambda <- solve(Gamma0) %*% gamma0
  pred <- obs %*% lambda[1 : 10]
  mspe <- t(lambda) %*% gamma0 - sigma.ep2
  c(pred = pred, mspe = mspe)
}

exp1 <- function(r){
  1 - exp(- r / 0.05)
}

# nl.out(1, exp1, 0, simdata.exp1$sim4)

new.xy <- new.xy[-c(101, 151, 301, 321, 351, 401, 501, 601, 751, 901), ]

exp.out.nl <- matrix(0, ncol = 2, nrow = 10)
for (i in 1 : 10){
  exp.out.nl[i, ] <- nl.out(i, exp1, 0.5, simdata.exp1$sim4)
}

pred.exp.nl <- krige.conv(coords = xy, data = simdata.exp1$sim4, locations = new.xy, borders = NULL,
                          krige = KC_exp4)
var.exp.nl <- krige.conv(coords = xy, data = simdata.exp1$sim4, locations = new.xy, borders = NULL,
                         krige = KC_exp4)$krige.var - 0.5

png("plot2_4.png", width = 600, height = 300)
plot(new.xy[,1], pred.exp.nl$predict, type = "l", xlab = "x", ylab = "y",
     ylim = c(-2.5, 2.5))
points(xy$x, exp.out.nl[,1], cex = 0.5)
points(xy$x, simdata.exp1$sim4, pch = 4)

lines(new.xy[,1], var.exp.nl, lty = 2, col = "red")
points(xy$x, exp.out.nl[,2], cex = 0.5,  col = "red")
legend("topright", legend = c("prediction", "MSPE"), col = c("black", "red"),
       lty = c(1, 2))
dev.off()

