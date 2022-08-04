sec_path = 'Rcode/Chapter8/Section 8.6/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Get the Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# attach data library
library(ZVHdata)
library(sp)
library(spmodel)

# load data for graphics and analysis
data(MOSSobs)
data(MOSSpreds)

DF = data.frame(MOSSobs@data, easting = MOSSobs@coords[,1]/1e+3,
	northing = MOSSobs@coords[,2]/1e+3)
DF$year = as.factor(DF$year)
DF$field_dup = as.factor(DF$field_dup)
DF$Zn = log(DF$Zn)
DF$Pb = log(DF$Pb)
DF$dist2road = log(DF$dist2road)

#undebug(splm)
#undebug(spmodel:::cov_estimate_gloglik_splm)
#undebug(spmodel:::get_data_object_splm)

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing', spcov_type = 'spherical',
#	partition_factor = ~ year,
#	spcov_initial = spcov_initial('spherical', range = 4, known = c('range')),
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(PbOut)
PbOut$optim$value

# change coordinate names from 'easting' and 'northing' to 'xcoords' and 'ycoords'
DF$xcoords = DF$easting
DF$ycoords = DF$northing
PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'xcoords', ycoord = 'ycoords', spcov_type = 'spherical',
#	partition_factor = ~ year,
#	spcov_initial = spcov_initial('spherical', range = 4, known = c('range')),
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-12), estmethod = 'reml')
summary(PbOut)
PbOut$optim$value

