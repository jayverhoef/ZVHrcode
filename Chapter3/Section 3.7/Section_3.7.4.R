library(ZVHdata)
library(spdep)
library(viridis)
library(classInt)
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
sealPolys_obs = sealPolys_obs[which(unlist(lapply(Nlist_obs,function(x){all(x!=0)}))),]
# again use poly2nb to create neighborhood list for each polygon
Nlist_obs = poly2nb(sealPolys_obs, snap = 2000)
attr(Nlist_obs,'polyid') = as.factor(as.character(sealPolys_obs@data$polyid))
attr(Nlist_obs,'stockid') = as.factor(as.character(sealPolys_obs@data$stockid))
coords_obs = coordinates(sealPolys_obs)

# ------------------------------------------------------------------------------
#     Compute local Moran's I
#-------------------------------------------------------------------------------

lmoran_out = localmoran(sealPolys_obs$Estimate, 
  listw = nb2listw(Nlist_obs, style = 'W'), zero.policy=TRUE,
  p.adjust.method = 'bonferroni', alternative = 'two.sided')

source('/media/jay/data/desktop_data/ActiveBooks/GeostatDale/Rcode/Chapter3/Section 3.7.4/addBreakColorLegend.R')
# use style 'fisher' to find initial breaks
f10 = classIntervals(lmoran_out[,'Z.Ii'], n = 10, style = 'fisher')
# modify breaks
fixedBreaks = c(-3.4, -1.8, -0.7, 0, 0.9, 0.3, 0.75, 1.3, 2, 2.75, 3.6)
f10 = classIntervals(lmoran_out[,'Z.Ii'], style = 'fixed', 
			fixedBreaks = fixedBreaks)
pal = viridis(length(fixedBreaks)-1)
f10Colours = findColours(f10, pal)
setwd('/media/jay/data/desktop_data/ActiveBooks/GeostatDale/Rcode/Chapter3/Section 3.7.4')
pdf(file = paste0('Cloropleth_LISA.pdf'), width = 10, height = 10)
  old_par = par(mar = c(5,5,3,1))
  plot(sealPolys_obs, col = f10Colours)
  points(coords_obs[lmoran_out[,'Pr(z != 0)'] < 0.05 & lmoran_out[,'Ii'] < 0,], 
    pch = 4, cex = 2, lwd = 4, col = 'maroon') 
  points(coords_obs[lmoran_out[,'Pr(z != 0)'] < 0.05 & lmoran_out[,'Ii'] > 0,], 
    pch = 3, cex = 2, lwd = 4, col = 'maroon') 
  addBreakColorLegend(1320000, 960000, 1350000, 1190000, f10$brks, 
    colors = pal, cex = 1.8, printFormat = "1.2")
dev.off()
system("pdfcrop '/media/jay/data/desktop_data/ActiveBooks/GeostatDale/Rcode/Chapter3/Section 3.7.4/Cloropleth_LISA.pdf'")

setwd('/mnt/ExtraDrive1/Work/desktop_data/ActiveBooks/GeostatDale/Rcode/Chapter3/Section 3.7')
pdf(file = paste0('plot_LISA.pdf'), width = 8, height = 8)
  old_par = par(mar = c(5,5,3,1))
  infl_out = moran.plot(sealPolys_obs$Estimate, 
    listw = nb2listw(Nlist_obs, style = 'W'), zero.policy=TRUE, cex =2,
    xlab = 'Centered Estimate', ylab = 'Sum of Weighted Neighbors', cex.lab = 2,
    cex.axis = 1.5, labels = FALSE, pch = 19, 
    xlim = c(-.6, .85), ylim = c(-.6, .85))
  par(old_par)
dev.off()
system("pdfcrop '/mnt/ExtraDrive1/Work/desktop_data/ActiveBooks/GeostatDale/Rcode/Chapter3/Section 3.7/plot_LISA.pdf'")


lmoran_out = localmoran(sealPolys_obs$Estimate, 
  listw = nb2listw(Nlist_obs, style = 'W'), zero.policy=TRUE)
lm(lmoran_out[,'Ii']~sealPolys_obs$Estimate)

lag.listw(nb2listw(Nlist_obs, style = 'B'), sealPolys_obs$Estimate, zero.policy = TRUE)

a = nb2listw(Nlist_obs, style = 'W')
b = sealPolys_obs$Estimate

s = NULL
for(i in 1:length(b))
s = c(s, sum(b[a$neighbours[[i]]]*a$weights[[i]]))
plot((b - mean(b)),s)
lm(s ~ I(b - mean(b))) 
moran(sealPolys_obs$Estimate, 
  nb2listw(Nlist_obs, style = 'W'), 
  length(Nlist_obs),
  Szero(nb2listw(Nlist_obs, style = 'W')))
