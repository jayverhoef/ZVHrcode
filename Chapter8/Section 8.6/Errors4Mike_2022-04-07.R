# attach data library
library(ZVHdata)
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

#------------------------- test all of them ------------------------------------

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'exponential',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'spherical',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'gaussian',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'circular',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'cubic',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'penta',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'wave',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'jbessel',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'gravity',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'rquad',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'magnetic',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'matern',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')
	
PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'cauchy',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')

PbOut = splm(Pb ~ year + dist2road + sideroad + sideroad:dist2road + 
	year:dist2road + year:sideroad + year:sideroad:dist2road, data = DF,
	xcoord = 'easting', ycoord = 'northing',   spcov_type = 'pexponential',
	partition_factor = ~ year,
	random = ~ -1 + sample + I(as.factor(paste(sample,field_dup))),
	control = list(reltol = 1e-7), estmethod = 'reml')
