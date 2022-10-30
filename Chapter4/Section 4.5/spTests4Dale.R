# set working directory to directory with data and spTest files
# setwd(....)
sec_path = 'Rcode/Chapter4/Section 4.5/'
setwd(paste0(SLEDbook_path,sec_path))

source('GuanTestGrid.R')
source('GuanTestUnif.R')
source('GuanTestGrid_help.R')
source('GuanTestUnif_help.R')
source('MaityTest.R')
source('MaityTest_help.R')
source('htestIso_class.R')

# read in the CSV file
# the column z is the square root of SO4 for cleaned data set
# easting and northing are coordinates in meters from CONUS Albers projection
DF = read.csv('DF.csv')

# fit the linear model on the coordinates
fit_2c <- lm(z ~ poly(easting, northing, degree = 2), data = DF)
#residuals from quadratic fit
r2 = residuals(fit_2c)

# create a data.frame with first two columns the spatial coordinates, and
# 3rd column the response, as required by GuanTestGrid help file,
# which is included as the header to GuanTestGrid.R
# we use the residuals as the response
# transform the coordinates by scaling them to be much smaller, 
# and translate them so they have lowest value at zero
dt_res<- cbind(DF[,c('easting','northing')]/100000, r2)
dt_res[,1] = dt_res[,1] - min(dt_res[,1])
dt_res[,2] = dt_res[,2] - min(dt_res[,2])

#check limits
min(dt_res[,1])
max(dt_res[,1])
min(dt_res[,2])
max(dt_res[,2])


# run the tests
# but I have no idea what the options mean
GuanTestUnif(as.matrix(dt_res), 
  lagmat = rbind(c(1, 0), c(0, 1), c(1,1),c(-1, 1)), 
  A = rbind(c(1, -1, 0, 0), c(0, 0, 1, -1)), df = 2, 
  xlims = c(0, 45), 
  ylims = c(0, 30), 
  grid.spacing = c(1,1),
  window.dims = c(15,15))

set.seed(10001)
MaityTest(as.matrix(dt_res), 
  lagmat = rbind(c(1, 0), c(0, 1), c(1, 1),c(-1, 1)), 
  A = rbind(c(1, -1, 0, 0), c(0, 0, 1, -1)), df = 2, 
  xlims = c(0, 45), 
  ylims = c(0, 30), 
  block.dims = c(15,10))
