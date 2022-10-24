sec_path = 'Rcode/Chapter3/Section 3.7/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(xtable)
library(spdep)
data(sealPolys)


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Neighbor_lags
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#     We want to do a correlogram on observed data only

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
attr(Nlist_obs,'polyid') = as.factor(as.character(sealPolys_obs$polyid))
attr(Nlist_obs,'stockid') = as.factor(as.character(sealPolys_obs$stockid))
coords_obs = st_drop_geometry(st_centroid(sealPolys_obs)[,c('x','y')])

#     First, make sure we understand neighbor and order defintions

# to help vizualize, plot 1st, 2nd, and 3rd order lags
# set zoom-in coordinates to examine more closely
xleft = 1268000
xright = 1390000
ybottom = 738000
ytop = 800000
xexp = (xright - xleft)*.06
yexp = (ytop - ybottom)*.35

# use nblag to create lists of neighbors by order
seal.lags <- nblag(Nlist_obs, 3)

file_name = 'figures/Neighbor_lags'

pdf(file = paste0(file_name,'.pdf'), width = 10, height = 10)

  # create the zoomed-in plot of all polygons
  plot(st_geometry(sealPolys), xlim = c(xleft, xright), 
    ylim = c(ybottom, ytop), )
  # add the zoomed-in plot of polygons with values
  plot(st_geometry(sealPolys_obs), add = TRUE,
    xlim = c(xleft, xright), ylim = c(ybottom, ytop), col = 'grey60')
  # plot the polygon centroids and 1st order connections
  plot(Nlist_obs, coords_obs, add = TRUE, lwd = 7, col = '#8dd3c7', cex = 3,
    pch = 19)
  # add the 2nd order connections
  plot(seal.lags[[2]], coords_obs, add=TRUE, col='#ffffb3', lwd = 4)
  # add the 3rd order connections
  plot(seal.lags[[3]], coords_obs, add=TRUE, col='#bebada', lwd = 3)

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
#                   Table: Number of neighbors by order
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# make a table of number of neighbors by order
# for correlogram, response variable is sealPolys_obs$Estimate
# compute the correlogram using Moran's I, which also computes number of nodes by
#    number of edges and lag
cgram_I = sp.correlogram(Nlist_obs, sealPolys_obs$Estimate,  style = "B",
  order=15, method="I", zero.policy=TRUE)
# maximum number of neighbors across all lags
nblagmax = max(unlist(lapply(cgram_I$cardnos, length)))
# maximum lag
nlags = length(cgram_I$cardnos)
# make a blank table for number of neighbors by lag
nb_lag = matrix(0, nrow = nblagmax, ncol = nlags)
# fill in table with numbers of neighbors by lag
for(i in 1:nlags)
  nb_lag[as.integer(row.names(cgram_I$cardnos[[i]])) + 1,i] = 
    cgram_I$cardnos[[i]]
    
print(
    xtable(nb_lag, 
      align = c('r',rep('r', times = length(nb_lag[1,]))),
      digits = rep(0, times = 16),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                   Correlogram using Moran's I
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Default plot of Moran's I
plot(cgram_I)
 
file_name = 'figures/Correlogram_I'

pdf(file = paste0(file_name,'.pdf'), width = 8, height = 8)

  old_par = par(mar = c(5,5,3,1))
  plot_vals = c(print(cgram_I)[,'expectation'] - 
      1.96*sqrt(print(cgram_I)[,'variance']),
    print(cgram_I)[,'expectation'] + 
      1.96*sqrt(print(cgram_I)[,'variance']),
    print(cgram_I)[,'estimate'])
  plot(c(.5,15.5), c(min(plot_vals),max(plot_vals)), type = 'n', xlab = 'Lag',
    ylab = "Moran's I", cex.lab = 2, cex.axis = 1.5)
  points(1:15, print(cgram_I)[,'expectation'], pch = 1, cex = 2)
  lines(1:15, print(cgram_I)[,'expectation'], lty = 2, lwd = 2)
  points(1:15, print(cgram_I)[,'expectation'] - 
    1.96*sqrt(print(cgram_I)[,'variance']), pch = 8, cex = 2)
  lines(1:15, print(cgram_I)[,'expectation'] - 
    1.96*sqrt(print(cgram_I)[,'variance']), lty = 2, lwd = 2)
  points(1:15, print(cgram_I)[,'expectation'] + 
    1.96*sqrt(print(cgram_I)[,'variance']), pch = 8, cex = 2)
  lines(1:15, print(cgram_I)[,'expectation'] + 
    1.96*sqrt(print(cgram_I)[,'variance']), lty = 2, lwd = 2)
  points(1:15, print(cgram_I)[,'estimate'], pch = 19, cex = 4)
  lines(1:15, print(cgram_I)[,'estimate'], lwd = 4)
  par(old_par)
  
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
