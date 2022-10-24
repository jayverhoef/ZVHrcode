sec_path = 'Rcode/Chapter3/Section 3.7/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(spdep)
library(viridis)
library(classInt)
source('addBreakColorLegend.R')
data(sealPolys)


# ------------------------------------------------------------------------------
#     We want to do a correlogram on observed data only
#-------------------------------------------------------------------------------

# get a list of all polygons that have observed values
ind_samp = which(!is.na(sealPolys$Estimate))
# subset the SpatialPolygonsDataFrame to those with observed values
sealPolys_obs = sealPolys[ind_samp,]
# use poly2nb to create neighborhood list for each polygon
Nlist_obs = poly2nb(sealPolys_obs, snap = 2000)
# some polygons have no neighbors, so remove them as well
sealPolys_obs = sealPolys_obs[which(unlist(lapply(
	Nlist_obs,function(x){all(x!=0)}))),]
# again use poly2nb to create neighborhood list for each polygon
Nlist_obs = poly2nb(sealPolys_obs, snap = 2000)
attr(Nlist_obs,'polyid') = as.factor(as.character(sealPolys_obs$polyid))
attr(Nlist_obs,'stockid') = as.factor(as.character(sealPolys_obs$stockid))
coords_obs = st_drop_geometry(st_centroid(sealPolys_obs)[,c('x','y')])


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         plot_LISA
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = 'figures/plot_LISA'

pdf(file = paste0(file_name,'.pdf'), width = 8, height = 8)

  old_par = par(mar = c(5,5,3,1))
  infl_out = moran.plot(sealPolys_obs$Estimate, 
    listw = nb2listw(Nlist_obs, style = 'W'), zero.policy=TRUE, cex =2,
    xlab = 'Centered Estimate', ylab = 'Sum of Weighted Neighbors', cex.lab = 2,
    cex.axis = 1.5, labels = FALSE, pch = 19, 
    xlim = c(-.6, .85), ylim = c(-.6, .85))
  par(old_par)


dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))


# reproduce plot_LISA from scratch

# a contains two lists.  The "neighbor" list gives indices of neighbors
# "weigths" list item gives weights for neighbors.  We use 'W' option on
# style which is row-standardized weights
a = nb2listw(Nlist_obs, style = 'W')
# create a vector of observed values
b = sealPolys_obs$Estimate

s = NULL
for(i in 1:length(b))
	# for each record, get weighted average of neighbors
  s = c(s, sum(b[a$neighbours[[i]]]*a$weights[[i]]))
# plot each response minus the overall mean by the weighted average of neighbors
plot((b - mean(b)),s)
# get the slope from the regression of weighted neighbors on response deviation
lm(s ~ I(b - mean(b)))

# Is the slope of the plot_LISA equal to Moran's I, as computed by spdep? 
moran(sealPolys_obs$Estimate, 
  nb2listw(Nlist_obs, style = 'W'), 
  length(Nlist_obs),
  Szero(nb2listw(Nlist_obs, style = 'W')))
# Yes!

# get record of possible outlier with centered estimate > 0.6
which(b - mean(b) > 0.6)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Cloropleth_LISA
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#  Compute local Moran's I

lmoran_out = localmoran(sealPolys_obs$Estimate, 
  listw = nb2listw(Nlist_obs, style = 'W'), zero.policy=TRUE,
  alternative = 'two.sided')

# everything here can be computed from two columns "Ii" and "Var.Ii"
# Z.Ii is a standard normal, so mean/std.dev
cbind(lmoran_out[,'Z.Ii'],
	lmoran_out[,'Ii']/sqrt(lmoran_out[,'Var.Ii']))

#assuming Z.Ii is standard normal, computed probabilities that it 
# is not equal to 0
cbind(2*(1-pnorm(abs(lmoran_out[,'Z.Ii']))),
	lmoran_out[,'Pr(z != E(Ii))'])


# compute some Bonferonni corrections
# at alpha 0.05 (for two-tailed test)
boncor05 = 0.025/length(sealPolys_obs$Estimate)
# at alpha 0.10 (for two-tailed test)
boncor10 = 0.05/length(sealPolys_obs$Estimate)
# at alpha 0.20 (for two-tailed test)
boncor20 = 0.10/length(sealPolys_obs$Estimate)
# get the indices of possible outliers at each alpha level for local-LISA
ind05 = which(abs(lmoran_out[,'Z.Ii']) > qnorm(1-boncor05))
ind10 = which(abs(lmoran_out[,'Z.Ii']) > qnorm(1-boncor10))
ind20 = which(abs(lmoran_out[,'Z.Ii']) > qnorm(1-boncor20))
ind05
ind10
ind20

# Get max and min local LISA's for creating fixed break points
max(lmoran_out[,'Z.Ii'])
min(lmoran_out[,'Z.Ii'])
# use style 'fisher' to find initial breaks
# f10 = classIntervals(lmoran_out[,'Z.Ii'], n = 10, style = 'fisher')
# modify breaks
fixedBreaks = c(-2.9, -1.8, -0.7, 0, 0.9, 0.3, 0.75, 1.3, 2, 2.75, 3.6)
f10 = classIntervals(lmoran_out[,'Z.Ii'], style = 'fixed', 
			fixedBreaks = fixedBreaks)
pal = viridis(length(fixedBreaks)-1)
f10Colours = findColours(f10, pal)

# plot Cloropleth_LISA
file_name = 'figures/Cloropleth_LISA'

pdf(file = paste0(file_name,'.pdf'), width = 16, height = 8)

	layout(matrix(1:2, nrow = 1))
	par(mar = c(0,0,0,0))
	plot(st_geometry(sealPolys_obs), col = f10Colours)
		addBreakColorLegend(1320000, 960000, 1350000, 1190000, f10$brks, 
			colors = pal, cex = 1.8, printFormat = "1.2")
	text(x = 990000, y = 1200000, labels = 'A', cex = 4)
	# plot LISA values
	plot(st_geometry(sealPolys_obs))
	# plot alpha = 0.01 per site
	points(coords_obs[lmoran_out[,'Pr(z != E(Ii))'] < 0.01 & 
			lmoran_out[,'Ii'] > 0,], pch = 19, col = 'red3', cex = 2)
	# plot Bonferonni-corrected alpha = 0.20
	# nothing changes for Bonferonni-corrected alpha = 0.10 
	points(coords_obs[ind20,], pch = 19, col = 'blue3', cex = 2.0)
	# there are no LISA outliers at Bonferonni-corrected alpha = 0.05
	# Plot the site with a global LISA value greater than 0.6
	points(coords_obs[infl_out$x > 0.6, ], pch = 19, col = 'green3',
			cex = 2.0)
	text(x = 990000, y = 1200000, labels = 'B', cex = 4)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
