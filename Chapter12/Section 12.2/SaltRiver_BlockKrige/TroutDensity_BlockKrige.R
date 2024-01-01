
## Load SSN package into R
library("SSN")

## Set working directory to location of .ssn directory
setwd("C:/...")

## Import the data from the .ssn directory and create a SpatialStreamNetwork
## object with basic set of prediction points for all reach midpoints
SaltWQ <- importSSN("LSN_TroutDensity_BlockKrige.ssn", predpts = "preds")

## Import prediction points spaced at 100m intervals for block-kriging estimates of individual streams
SaltWQ <- importPredpts(SaltWQ, "cottonwood", "ssn")
SaltWQ <- importPredpts(SaltWQ, "crow", "ssn")
SaltWQ <- importPredpts(SaltWQ, "dry", "ssn")
SaltWQ <- importPredpts(SaltWQ, "jackknife", "ssn")
SaltWQ <- importPredpts(SaltWQ, "salt", "ssn")
SaltWQ <- importPredpts(SaltWQ, "spring", "ssn")
SaltWQ <- importPredpts(SaltWQ, "strawberry", "ssn")
SaltWQ <- importPredpts(SaltWQ, "stump", "ssn")
SaltWQ <- importPredpts(SaltWQ, "swift", "ssn")
SaltWQ <- importPredpts(SaltWQ, "tincup", "ssn")
SaltWQ <- importPredpts(SaltWQ, "willow", "ssn")
SaltWQ <- importPredpts(SaltWQ, "saltriver", "ssn")

## Import prediction points spaced at 100m intervals for block-kriging estimate of full network
SaltWQ <- importPredpts(SaltWQ, "network", "ssn")

## Create distance matrices among stream prediction points - these have been provided with the
## example dataset.
## createDistMat(SaltWQ, predpts = "preds", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Cottonwood", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Crow", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Dry", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Jackknife", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Salt", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Spring", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Strawberry", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Stump", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Swift", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Tincup", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "Willow", o.write = TRUE, amongpreds = TRUE)
## createDistMat(SaltWQ, predpts = "SaltRiver", o.write = TRUE, amongpreds = TRUE)

## Create distance matrix among network prediction points (calculations may require
## a few minutes)
createDistMat(SaltWQ, predpts = "network", o.write = TRUE, amongpreds = TRUE)

## Describe the names of the variables in the point.data data.frame for each
## observed and prediction data set
names(SaltWQ)

## Plot Salt River network and locations of 108 trout density observations
plot(SaltWQ, lwdLineCol = "afvArea", lwdLineEx = 5, lineCol = "blue",
             pch = 19, xlab = "x-coordinate (m)",
             ylab = "y-coordinate (m)", asp = 1)

## Plot values of 108 trout density observations (Figure 1)
brks <- plot(SaltWQ, "trout_100m", lwdLineCol = "afvArea",
             lwdLineEx = 5, lineCol = "black", xlab = "x-coordinate" ,
             ylab = "y-coordinate", asp=1 )

##plot Torgegram based on 108 trout density observations (Figure 2)
SaltWQ.Torg <- Torgegram(SaltWQ, "trout_100m", nlag = 15,
             nlagcutoff = 1, maxlag = 50000)
plot(SaltWQ.Torg)

## Fit nonspatial multiple linear regression (MLR) in Table 2
SaltWQ.glmssn0 <- glmssn(trout_100m ~ SLOPE + S1_93_11 + CANOPY, SaltWQ,
             CorModels = NULL, use.nugget = TRUE, EstMeth = "REML")
summary(SaltWQ.glmssn0)

## Fit SSN1 in Table 2.
SaltWQ.glmssn1 <- glmssn(trout_100m ~ SLOPE + S1_93_11 + CANOPY, SaltWQ,
             CorModels = c("Exponential.tailup", "Exponential.taildown"),
             addfunccol = "afvArea", EstMeth = "REML")
summary(SaltWQ.glmssn1)

## Fit SSN2 in Table 2.
SaltWQ.glmssn2 <- glmssn(trout_100m ~ SLOPE + S1_93_11 + CANOPY, SaltWQ,
             CorModels = c("Exponential.tailup", "Exponential.taildown",
             "Exponential.Euclid"),addfunccol = "afvArea", EstMeth = "REML")
summary(SaltWQ.glmssn2)

## Fit SSN3 in Table 2.
SaltWQ.glmssn3 <- glmssn(trout_100m ~ S1_93_11, SaltWQ,
             CorModels = c("Exponential.tailup", "Exponential.taildown"),
             addfunccol = "afvArea", EstMeth = "REML")
summary(SaltWQ.glmssn3)

## Fit SSN4 in Table 2.
SaltWQ.glmssn4 <- glmssn(trout_100m ~ 1, SaltWQ,
             CorModels = c("Exponential.tailup", "Exponential.taildown"),
             addfunccol = "afvArea", EstMeth = "REML")
summary(SaltWQ.glmssn4)

## Fit SSN5 in Table 2.
SaltWQ.glmssn5 <- glmssn(trout_100m ~ 1, SaltWQ,
             CorModels = c("Exponential.tailup", "Exponential.taildown",
             "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "REML")
summary(SaltWQ.glmssn5)

## Report AIC values (Use ML instead of REML in above model fits to estimate
## correct AIC values)
AIC(SaltWQ.glmssn0)
AIC(SaltWQ.glmssn1)
AIC(SaltWQ.glmssn2)
AIC(SaltWQ.glmssn3)
AIC(SaltWQ.glmssn4)
AIC(SaltWQ.glmssn5)

## Report cross-validation statistics and confidence intervals
CrossValidationStatsSSN(SaltWQ.glmssn0)
CrossValidationStatsSSN(SaltWQ.glmssn1)
CrossValidationStatsSSN(SaltWQ.glmssn2)
CrossValidationStatsSSN(SaltWQ.glmssn3)
CrossValidationStatsSSN(SaltWQ.glmssn4)
CrossValidationStatsSSN(SaltWQ.glmssn5)

## Report variance composition among covariate effects and autocovariance
## functions
varcomp(SaltWQ.glmssn0)
varcomp(SaltWQ.glmssn1)
varcomp(SaltWQ.glmssn2)
varcomp(SaltWQ.glmssn3)
varcomp(SaltWQ.glmssn4)
varcomp(SaltWQ.glmssn5)

## Plot graphs of leave-one-out cross-validation (LOOCV) predictions & SEs
cv.out <- CrossValidationSSN(SaltWQ.glmssn2)
par(mfrow = c(1, 1))
plot(SaltWQ.glmssn2$sampinfo$z,cv.out[, "cv.pred"], pch = 19,
     xlab = "Observed Data", ylab = "LOOCV Prediction")

## Save LOOCV predictions & SEs to working directory file
write.csv(cv.out, "cv_out_trout100m.csv", row.names = FALSE)

## Calculate & plot model residuals & influence measures
resids <- residuals(SaltWQ.glmssn2)
class(resids)
resids.df <- getSSNdata.frame(resids)
names(resids.df)
plot(resids)
hist(resids, xlab = "Raw Residuals")
qqnorm(resids)

## Save residuals & influence measures to working directory file
write.csv(resids.df, "resids_trout100m.csv", row.names = FALSE)

## Plot 108 observation sites as large circles prior to block-kriging
## prediction points
plot(SaltWQ, "trout_100m", pch = 1, cex = 3, xlab = "x-coordinate (m)",
     ylab = "y-coordinate (m)",
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

## Plot 100m prediction points for individual stream blocks
SaltWQ.glmssn2.cottonwood <- predict(SaltWQ.glmssn2, "cottonwood")
plot(SaltWQ.glmssn2.cottonwood, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.crow <- predict(SaltWQ.glmssn2, "crow")
plot(SaltWQ.glmssn2.crow, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.dry <- predict(SaltWQ.glmssn2, "dry")
plot(SaltWQ.glmssn2.dry, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.jackknife <- predict(SaltWQ.glmssn2, "jackknife")
plot(SaltWQ.glmssn2.jackknife, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.salt <- predict(SaltWQ.glmssn2, "salt")
plot(SaltWQ.glmssn2.salt, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.spring <- predict(SaltWQ.glmssn2, "spring")
plot(SaltWQ.glmssn2.spring, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.strawberry <- predict(SaltWQ.glmssn2, "strawberry")
plot(SaltWQ.glmssn2.strawberry, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.stump <- predict(SaltWQ.glmssn2, "stump")
plot(SaltWQ.glmssn2.stump, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.swift <- predict(SaltWQ.glmssn2, "swift")
plot(SaltWQ.glmssn2.swift, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.tincup <- predict(SaltWQ.glmssn2, "tincup")
plot(SaltWQ.glmssn2.tincup, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.willow <- predict(SaltWQ.glmssn2, "willow")
plot(SaltWQ.glmssn2.willow, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

SaltWQ.glmssn2.saltriver <- predict(SaltWQ.glmssn2, "saltriver")
plot(SaltWQ.glmssn2.saltriver, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

## Obtain block-kriging estimates of mean trout density & SEs
## for individual stream blocks
SaltWQ.glmssn2.cottonwood <- BlockPredict(SaltWQ.glmssn2, "cottonwood")
SaltWQ.glmssn2.cottonwood

SaltWQ.glmssn2.crow <- BlockPredict(SaltWQ.glmssn2, "crow")
SaltWQ.glmssn2.crow

SaltWQ.glmssn2.dry <- BlockPredict(SaltWQ.glmssn2, "dry")
SaltWQ.glmssn2.dry

SaltWQ.glmssn2.jackknife <- BlockPredict(SaltWQ.glmssn2, "jackknife")
SaltWQ.glmssn2.jackknife

SaltWQ.glmssn2.salt <- BlockPredict(SaltWQ.glmssn2, "salt")
SaltWQ.glmssn2.salt

SaltWQ.glmssn2.spring <- BlockPredict(SaltWQ.glmssn2, "spring")
SaltWQ.glmssn2.spring

SaltWQ.glmssn2.strawberry <- BlockPredict(SaltWQ.glmssn2, "strawberry")
SaltWQ.glmssn2.strawberry

SaltWQ.glmssn2.stump <- BlockPredict(SaltWQ.glmssn2, "stump")
SaltWQ.glmssn2.stump

SaltWQ.glmssn2.swift <- BlockPredict(SaltWQ.glmssn2, "swift")
SaltWQ.glmssn2.swift

SaltWQ.glmssn2.tincup <- BlockPredict(SaltWQ.glmssn2, "tincup")
SaltWQ.glmssn2.tincup

SaltWQ.glmssn2.willow <- BlockPredict(SaltWQ.glmssn2, "willow")
SaltWQ.glmssn2.willow

SaltWQ.glmssn2.SaltRiver <- BlockPredict(SaltWQ.glmssn2, "saltriver")
SaltWQ.glmssn2.SaltRiver

## Save values of block predictions & SEs at 100m prediction points to working
## directory file
SaltWQ.cottonwood <- predict(SaltWQ.glmssn2, "cottonwood")
pred1df <- getSSNdata.frame(SaltWQ.cottonwood, "cottonwood")
write.csv(pred1df, "SaltWQ_trout100m_SSN2_cottonwood_BlockPredictions.csv",
     row.names = FALSE)

## Plot 100m prediction points for full network
SaltWQ.glmssn2.network <- predict(SaltWQ.glmssn2, "network")
plot(SaltWQ.glmssn2.network, "trout_100m", add = TRUE,
     xlim = c(1715000,1770000), ylim = c(1370000,1450000))

## Obtain block-kriging estimates of mean trout density & SEs for full network
## (calculation requires several minutes)
SaltWQ.glmssn2.network <- BlockPredict(SaltWQ.glmssn2, "network")
SaltWQ.glmssn2.network

## Plot predictions for full network at reach midpoints with symbol
## size inverse to SEs
SaltWQ.preds <- predict(SaltWQ.glmssn2, "preds")
plot(SaltWQ.preds, SEcex.max = 1.4, SEcex.min = .7/3*2,
     breaktype = "user", brks = brks)

## Is this cool or what?


