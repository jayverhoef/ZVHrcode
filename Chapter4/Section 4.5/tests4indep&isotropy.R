> library(spmodel)
> library(sf)
> library(sm)
Linking to GEOS 3.10.2, GDAL 3.4.3, PROJ 8.2.1; sf_use_s2() is TRUE
> sul_coord <- st_coordinates(sulfate$geometry)
> sul.raw <- data.frame(sulfate = sulfate$sulfate, x = sul_coord[,1], y = sul_coord[,2])

# Remove 3 spatial outliers
> sul <- sul.raw[!(sul.raw$sulfate == 44.023 |
> sul.raw$sulfate == 4.875 | sul.raw$sulfate == 5.225),]

# Take square root of responses
> sul$sulfate <- sqrt(sul$sulfate)

# Estimate 2nd-order polynomial surface by OLS
> spm2 <- splm(sul$sulfate ~ polym(x,y,degree=2),data=sul,xcoord=x,ycoord=y,spcov_type="none")

> summary(spm2)

Call:
splm(formula = sul$sulfate ~ polym(x, y, degree = 2), data = sul, 
    spcov_type = "none", xcoord = x, ycoord = y)

Residuals:
     Min       1Q   Median       3Q      Max 
-2.10114 -0.51571 -0.01522  0.50395  2.05308 

Coefficients (fixed):
                            Estimate Std. Error z value Pr(>|z|)    
(Intercept)                  3.13725    0.05511  56.926  < 2e-16 ***
polym(x, y, degree = 2)1.0  17.25645    0.82098  21.019  < 2e-16 ***
polym(x, y, degree = 2)2.0  -0.62026    0.82119  -0.755   0.4501    
polym(x, y, degree = 2)0.1  -1.28083    0.82249  -1.557   0.1194    
polym(x, y, degree = 2)1.1 -37.25754   11.73575  -3.175   0.0015 ** 
polym(x, y, degree = 2)0.2  -4.53251    0.80795  -5.610 2.02e-08 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Pseudo R-squared: 0.7329

Coefficients ( spatial covariance):
    ie 
0.5777

# Compute the square root of the residual mean square
> sqrt(.5777)
[1] 0.7600658

# Perform Diblasi-Bowman tests for independence, stationarity, and isotropy of the residuals from the fitted 2nd-order surface
> sph.m2 <- cbind(x=sul$x,y=sul$y,res=summary(spm2)[["residuals"]]$raw)
> sm.variogram(sph.m2[,1:2],sph.m2[,3],model="independent")
Test of spatial independence: p =  0.006
> sm.variogram(sph.m2[,1:2],sph.m2[,3],model="isotropic",display="none")
Test of isotropy: p =  0.314 
> sm.variogram(sph.m2[,1:2],sph.m2[,3],model="stationary")
Test of stationarity: p =  0.242 

# Estimate 3rd-order polynomial surface by OLS
> spm3 <- splm(sul$sulfate ~ polym(x,y,degree=3),data=sul,xcoord=x,ycoord=y,
+ spcov_type="none")
> summary(spm3)

Call:
splm(formula = sul$sulfate ~ polym(x, y, degree = 3), data = sul, 
    spcov_type = "none", xcoord = x, ycoord = y)

Residuals:
     Min       1Q   Median       3Q      Max 
-1.72920 -0.30947  0.07489  0.32183  1.29428 

Coefficients (fixed):
                            Estimate Std. Error z value Pr(>|z|)    
(Intercept)                  3.07306    0.04266  72.032  < 2e-16 ***
polym(x, y, degree = 3)1.0  16.10199    0.59030  27.278  < 2e-16 ***
polym(x, y, degree = 3)2.0  -2.01378    0.65129  -3.092  0.00199 ** 
polym(x, y, degree = 3)3.0  -5.81609    0.63000  -9.232  < 2e-16 ***
polym(x, y, degree = 3)0.1  -1.37730    0.65446  -2.104  0.03534 *  
polym(x, y, degree = 3)1.1  17.62682   10.03098   1.757  0.07888 .  
polym(x, y, degree = 3)2.1  51.94175    9.26901   5.604 2.10e-08 ***
polym(x, y, degree = 3)0.2  -4.25673    0.62845  -6.773 1.26e-11 ***
polym(x, y, degree = 3)1.2 -56.61714    9.13202  -6.200 5.65e-10 ***
polym(x, y, degree = 3)0.3  -1.10205    0.59886  -1.840  0.06573 .  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Pseudo R-squared: 0.873

Coefficients ( spatial covariance):
    ie 
0.2806

# Compute square root of residual mean square
> sqrt(.2806)
[1] 0.5297169

# Perform Diblasi-Bowman tests for independence and stationarity of the residuals from the fitted 3rd-order surface
> sph.m3 <- cbind(x=sul$x,y=sul$y,res=summary(spm3)[["residuals"]]$raw)
> sm.variogram(sph.m3[,1:2],sph.m3[,3],model="independent")
Test of spatial independence: p =  0.22
> sm.variogram(sph.m3[,1:2],sph.m3[,3],model="isotropic",display="none")
Test of isotropy: p =  0.028
> sm.variogram(sph.m3[,1:2],sph.m3[,3],model="stationary")
Test of stationarity: p =  0.463 

# Prepare data for Guan-Sherman-Calvin test and Maity-Sherman test.  Unlike the previous tests, we need to scale 
# things in a certain way.  Specifically, we scale the study area to the rectangle [0, 26.4] × [0, 16]. This is
# because we need to keep the ratio of width to height equal to the ratio of the difference in longitude between
# the westernmost and easternmost data location (for the numerator) to the difference in latitude between the
# northernmost and southernmost data location (for the denominator). Since the width-to-height ratio is
# close to 1.65 (up to 2 decimals), we set the width to 16 × 1.65 = 26.4.
# For the sampling region, I consider the first example provided in the illustration document of these functions
# in the spTest package. The tricky thing here is that the sampling region (defined as the xlim times ylim)
# should be a multiple for the window sizes. Otherwise, these functions would not work. To make the sampling
# region has a similar scale of the scaled study region, I set the xlim = c(0, 27) and ylim = c(0, 16) with
# corresponding window size is window.dims = c(3, 4).)

> sul_d <- st_coordinates(sulfate$geometry)
> (max(sul_d[,1]) - min(sul_d[,1])) / (max(sul_d[,2]) - min(sul_d[,2]))
[1] 1.651473

> sul.raw <- data.frame(sulfate = sulfate$sulfate,
x = 26.4*(sul_d[,1] - min(sul_d[,1]) ) / (max(sul_d[,1]) - min(sul_d[,1])),
y = 16*(sul_d[,2] - min(sul_d[,2]) ) / (max(sul_d[,2]) - min(sul_d[,2])))

> # remove spatial outliers
sul <- sul.raw[!(sul.raw$sulfate == 44.023 | sul.raw$sulfate == 4.875 | sul.raw$sulfate == 5.225),]

> # square root on responses
sul$sulfate <- sqrt(sul$sulfate)

> m2 <- splm(sulfate ~ polym(x, y, degree = 2), data = sul,
xcoord = x, ycoord = y, spcov_type = "none", estmethod = "reml")

> sph.m2 <- cbind(x = sul$x, y = sul$y, res = summary(m2)[["residuals"]]$raw)

> library("spTest")
# Now perform Guan-Sherman-Calvin Test
mylags = rbind(c(1,0), c(0, 1), c(1, 1), c(-1,1))
myA = rbind(c(1, -1, 0 , 0), c(0, 0, 1, -1))
my.grid = c(1,1)
my.windims = c(3, 4)
myh = 0.7 #0.7
myh.sb = 0.7 #0.8
my.xlims = c(0, 27)
my.ylims = c(0, 16)
trace(GuanTestUnif, edit = TRUE)
## [1] "GuanTestUnif"
## in line 6
## changing "class(spdata)" to "class(spdata)[1]"
### one may refer to https://youtu.be/-taaGUHxOng
### or https://stackoverflow.com/questions/24331690/modify-package-function
# order 2
set.seed(6350)
Tracing function "GuanTestUnif" in package "spTest"
[1] "GuanTestUnif"
> GuanTestUnif(sph.m2, mylags, myA, df = 2, myh, "norm", 1.5,
my.xlims, my.ylims, my.grid, my.windims, myh.sb)

	Test of isotropy from Guan et. al. (2004) for uniformly distributed
	sampling locations using the sample semivariogram.

data:  sph.m2
Chi-sq = 1.6632, df = 2, p-value = 0.4353
p-value (finite adj.) = 0.1598, number of subblocks: 244
alternative hypothesis: true difference in directional semivariograms is not equal to 0

sample estimates: (lag value)
    (1,0)     (0,1)     (1,1)    (-1,1) 
0.1940539 0.1854232 0.2206749 0.2585273 

estimated asymp. variance-covariance matrix:
           [,1]       [,2]       [,3]       [,4]
[1,] 0.11619425 0.07893779 0.09577720 0.10870999
[2,] 0.07893779 0.12083585 0.09916332 0.13144003
[3,] 0.09577720 0.09916332 0.16157435 0.08387689
[4,] 0.10870999 0.13144003 0.08387689 0.24180342

Warning messages:
1: In if (spdata.class == "geodata") { :
  the condition has length > 1 and only the first element will be used
2: In if (spdata.class == "SpatialPointsDataFrame") { :
  the condition has length > 1 and only the first element will be used
3: In est_subblock_irreg(spdata, lagmat, xlims, ylims, window.dims,  :
  72  subblocks discarded due to inadequate data in the block (n < 4); consider increasing block size

# Now perform Maity-Sherman Test
> trace(MaityTest, edit = TRUE)
## [1] "MaityTest"
Tracing function "MaityTest" in package "spTest"
[1] "MaityTest"
> ## in line 6
## changing "class(spdata)" to "class(spdata)[1]"
# order 2
set.seed(6350)
MaityTest(sph.m2, mylags, myA, 2,
xlims = my.xlims, ylims = my.ylims, block.dims = my.windims)

	Test of isotropy from Maity and Sherman (2012) for sampling locations
	with general sampling design using the sample covariogram.

data:  sph.m2
Chi-sq = 2.0812, df = 2, p-value = 0.3532

alternative hypothesis: true difference in directional covariograms is not equal to 0

sample estimates: (lag value)
    (1,0)     (0,1)     (1,1)    (-1,1) 
0.2973276 0.3475199 0.3041686 0.2646470 

estimated asymp. variance-covariance matrix:
            [,1]        [,2]        [,3]        [,4]
[1,] 0.003864858 0.003219634 0.002194664 0.002537900
[2,] 0.003219634 0.004539704 0.002895372 0.003226149
[3,] 0.002194664 0.002895372 0.002590042 0.001929606
[4,] 0.002537900 0.003226149 0.001929606 0.003191972

