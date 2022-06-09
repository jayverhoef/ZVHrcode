sec_path = 'Rcode/Chapter3/Section 3.7/'
setwd(paste0(SLEDbook_path,sec_path))

################################################################################
#-------------------------------------------------------------------------------
#                            sim-MoranGeary
#-------------------------------------------------------------------------------
################################################################################
# nearest neighbor (NN) scatterplot
file_name = 'sim-MoranGeary'
library(spdep)
Site <- c(1:9)
set.seed(1007)
Response <- 1+Site+rnorm(9)
pdf(paste0(file_name,'.pdf'))
  old.par = par(mar = c(5,5,1,1))
  plot(Site, Response, xlim=c(0,10), ylim=c(0,10), xlab="Site", ylab="Response",
    pch = 19, cex = 2, cex.lab = 2, cex.axis = 1.5)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

#creat neighborhood list by hand
Nlist_sim = list(as.integer(2), as.integer(c(1,3)), as.integer(c(2,4)), 
  as.integer(c(3,5)), as.integer(c(4,6)), as.integer(c(5,7)), 
  as.integer(c(6,8)), as.integer(c(7,9)), as.integer(c(8)))
class(Nlist_sim) = 'nb'
attr(Nlist_sim,'type') = 'rook'
attr(Nlist_sim,'sym') = TRUE
str(Nlist)
# Moran's I
moran(Response, 
  nb2listw(Nlist_sim, style = 'B'), 
  length(Nlist_sim),
  Szero(nb2listw(Nlist_sim, style = 'B')))
# Normality test for Moran's I  
moran.test(Response, 
  nb2listw(Nlist_sim, style = 'B'),
  randomisation=FALSE)
# Randomization test for Moran's I
moran.test(Response, 
  nb2listw(Nlist_sim, style = 'B'),
  randomisation=TRUE)
# Geary's C
geary(Response, 
  nb2listw(Nlist_sim, style = 'B'), 
  length(Nlist_sim), length(Nlist_sim) - 1,
  Szero(nb2listw(Nlist_sim, style = 'B')))
# Normality test for Geary's C
geary.test(Response, 
  nb2listw(Nlist_sim, style = 'B'),
  randomisation=FALSE)
# Randomization test for Geary's C  
geary.test(Response, 
  nb2listw(Nlist_sim, style = 'B'),
  randomisation=TRUE)



library(ZVHdata)
library(xtable)
data(sealPolys)
# ------------------------------------------------------------------------------
#     Code from Introduction, a Reminder about the Structure of the Data
#-------------------------------------------------------------------------------
#find neighboring polygons using poly2nb function from spdep package
Nlist = poly2nb(sealPolys, row.names = sealPolys$polyid, snap = 2000)
# Nlist is a list object.  However, some polygons have no neighbors. We want
# all polygons to have at least one neighbor, so we added the nearest neighbor(s)
# by visual inspection
Nlist[[79]] = as.integer(c(211, 463))
Nlist[[211]] = as.integer(c(Nlist[[211]]), 79)
Nlist[[463]] = as.integer(c(Nlist[[463]]), 79)
Nlist[[130]] = as.integer(302)
Nlist[[302]] = as.integer(c(Nlist[[302]],130))
Nlist[[325]] = as.integer(c(326, 353))
Nlist[[326]] = as.integer(c(Nlist[[326]],325))
Nlist[[353]] = as.integer(c(Nlist[[353]],325))
Nlist[[435]] = as.integer(437)
Nlist[[437]] = as.integer(c(Nlist[[437]],435))
Nlist[[436]] = as.integer(c(86, 88))
Nlist[[86]] = as.integer(c(Nlist[[86]],436))
Nlist[[88]] = as.integer(c(Nlist[[88]],436))
Nlist[[437]] = as.integer(87)
Nlist[[87]] = as.integer(c(Nlist[[87]],437))
Nlist[[438]] = as.integer(436)
Nlist[[436]] = as.integer(c(Nlist[[436]],438))
Nlist[[439]] = as.integer(346)
Nlist[[346]] = as.integer(c(Nlist[[346]],439))
Nlist[[443]] = as.integer(281)
Nlist[[281]] = as.integer(c(Nlist[[281]],443))
Nlist[[463]] = as.integer(79)
# is Nlist in the same order as sealPolys@data?
all(attr(Nlist,'region.id') == sealPolys@data$polyid)
# attributes for region.id has extra factor levels, so condense
attr(Nlist,'region.id') = as.factor(as.character(sealPolys@data$polyid))
# add stockid as an attribute
attr(Nlist,'stockid') = as.factor(as.character(sealPolys@data$stockid))

# get samples sizes for total and observed and missing data
nTot = length(sealPolys)
nObs = sum(!is.na(sealPolys$Estimate))
nMiss = nTot - nObs

# number of neighbors for each polygon
num = card(Nlist)
# a function to create a neighbor matrix from a list
# adj contains the indexes of neighbors
# num contains the number of neighbors
# n is the total number of polygons
Neighmat = function (adj, num, n) 
{
    N.mat <- matrix(0, nrow = n, ncol = n)
    k <- 1
    for (i in 1:n) {
        if (num[i] > 0) {
            N.mat[i, adj[[i]]] <- 1
        }
    }
    row.names(N.mat) = attr(adj,'region.id')
    colnames(N.mat) = attr(adj,'region.id')
    N.mat
}

Nmat = Neighmat(Nlist, num, nTot)
# make it symmetric
Nmat1 = pmax(Nmat,t(Nmat))
# find neighbors of neighbors 
Nmat2 = (Nmat1 %*% Nmat1 > 0 | Nmat1 > 0)*1
diag(Nmat2) = 0
# find neighbors of neighbors of neighbors
Nmat4 = (Nmat2 %*% Nmat2 > 0 | Nmat2 > 0)*1
diag(Nmat4) = 0

# visualize the neighbors
xleft = 1268000
xright = 1390000
ybottom = 738000
ytop = 800000
xexp = (xright - xleft)*.06
yexp = (ytop - ybottom)*.35
  
# make a plot of neighbors
coords = coordinates(sealPolys)
layout(matrix(c(1,1,2,1,1,1,3,1,1), nrow = 3, byrow = TRUE))
par(mar = c(0,0,0,0))
plot(sealPolys)
rect(xleft - xexp, ybottom - yexp, xright + xexp, ytop + yexp, col = 'grey80', 
border = NA)
plot(sealPolys, add = TRUE)
plot(Nlist, coords, add = TRUE, lwd = 2)
text(920000, 1190000, labels = 'A', cex = 6)
Nlist2 = apply(Nmat2, 1, function(x) as.integer(which(x > 0)))
class(Nlist2) = 'nb'
attr(Nlist2,'type') = 'queen'
attr(Nlist2,'sym') = TRUE
attr(Nlist2,'polyid') = attr(Nlist,'polyid')
attr(Nlist2,'stockid') = attr(Nlist,'polyid')
plot(sealPolys, xlim = c(xleft, xright), ylim = c(ybottom, ytop))
plot(Nlist2, coords, add = TRUE, lwd = 1)
mtext('B', adj = -.15, padj = 1.1, cex = 4)
Nlist4 = apply(Nmat4, 1, function(x) as.integer(which(x > 0)))
class(Nlist4) = 'nb'
attr(Nlist4,'type') = 'queen'
attr(Nlist4,'sym') = TRUE
attr(Nlist4,'polyid') = attr(Nlist,'polyid')
attr(Nlist4,'stockid') = attr(Nlist,'polyid')
plot(sealPolys, xlim = c(xleft, xright), ylim = c(ybottom, ytop))
plot(Nlist4, coords, add = TRUE, lwd = .5)
mtext('C', adj = 0, padj = -.4, cex = 4)

# ------------------------------------------------------------------------------
#     We want to do Moran's I and Geary's c on observed data only
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
# plot the neighbors for observed polygons only
plot(sealPolys_obs)
coords_obs = coordinates(sealPolys_obs)
plot(Nlist_obs, coords_obs, add = TRUE, lwd = 2)

# number of neighbors
num_obs = lapply(Nlist_obs, function(x) length(x))
num_obs = unlist(num_obs)
# turn list into a matrix
Nmat_obs = Neighmat(Nlist_obs, num_obs, length(num_obs))
# make it symmetric
Nmat1_obs = pmax(Nmat_obs,t(Nmat_obs))
# find neighbors of neighbors 
Nmat2_obs = (Nmat1_obs %*% Nmat1_obs > 0 | Nmat1_obs > 0)*1
diag(Nmat2_obs) = 0
# find neighbors of neighbors of neighbors
Nmat4_obs = (Nmat2_obs %*% Nmat2_obs > 0 | Nmat2_obs > 0)*1
diag(Nmat4_obs) = 0
# now turn matrices back into neighborhood lists
Nlist2_obs = apply(Nmat2_obs, 1, function(x) as.integer(which(x > 0)))
class(Nlist2_obs) = 'nb'
attr(Nlist2_obs,'type') = 'queen'
attr(Nlist2_obs,'sym') = TRUE
attr(Nlist2_obs,'polyid') = attr(Nlist,'polyid')
attr(Nlist2_obs,'stockid') = attr(Nlist,'polyid')
Nlist4_obs = apply(Nmat4_obs, 1, function(x) as.integer(which(x > 0)))
class(Nlist4_obs) = 'nb'
attr(Nlist4_obs,'type') = 'queen'
attr(Nlist4_obs,'sym') = TRUE
attr(Nlist4_obs,'polyid') = attr(Nlist,'polyid')
attr(Nlist4_obs,'stockid') = attr(Nlist,'polyid')
# plot the 2nd-order neighborhood map
plot(sealPolys_obs)
plot(Nlist2_obs, coords_obs, add = TRUE, lwd = 1)
# plot the 4th-order neighborhood map
plot(sealPolys_obs)
plot(Nlist4_obs, coords_obs, add = TRUE, lwd = 1)

# now use spdep functions on new dataset
# Here we use the binary weights
moran1_I = moran(sealPolys_obs$Estimate, 
  nb2listw(Nlist_obs, style = 'B'), 
  length(Nlist_obs),
  Szero(nb2listw(Nlist_obs, style = 'B')))
moran1_Itest_norm = moran.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist_obs, style = 'B'),
  randomisation=FALSE)
moran1_Itest_rand = moran.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist_obs, style = 'B'),
  randomisation=TRUE)
moran1_Itest_MC = moran.mc(sealPolys_obs$Estimate, 
  listw = nb2listw(Nlist_obs, style = 'B'), nsim = 10000)

geary1_C = geary(sealPolys_obs$Estimate, 
  nb2listw(Nlist_obs, style = 'B'), 
  length(Nlist_obs), length(Nlist_obs) - 1,
  Szero(nb2listw(Nlist_obs, style = 'B')))
geary1_Ctest_norm = geary.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist_obs, style = 'B'),
  randomisation=FALSE)
geary1_Ctest_rand = geary.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist_obs, style = 'B'),
  randomisation=TRUE)
geary1_Ctest_MC = geary.mc(sealPolys_obs$Estimate, 
  listw = nb2listw(Nlist_obs, style = 'B'), nsim = 10000)

# compute by adding 2nd order neighbors
moran2_I = moran(sealPolys_obs$Estimate, 
  nb2listw(Nlist2_obs, style = 'B'), 
  length(Nlist2_obs),
  Szero(nb2listw(Nlist2_obs, style = 'B')))
moran2_Itest_norm = moran.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist2_obs, style = 'B'),
  randomisation=FALSE)
moran2_Itest_rand = moran.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist2_obs, style = 'B'),
  randomisation=TRUE)
moran2_Itest_MC = moran.mc(sealPolys_obs$Estimate, 
  listw = nb2listw(Nlist2_obs, style = 'B'), nsim = 10000)

geary2_C = geary(sealPolys_obs$Estimate, 
  nb2listw(Nlist2_obs, style = 'B'), 
  length(Nlist2_obs), length(Nlist_obs) - 1,
  Szero(nb2listw(Nlist2_obs, style = 'B')))
geary2_Ctest_norm = geary.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist2_obs, style = 'B'),
  randomisation=FALSE)
geary2_Ctest_rand = geary.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist2_obs, style = 'B'),
  randomisation=TRUE)
geary2_Ctest_MC = geary.mc(sealPolys_obs$Estimate, 
  listw = nb2listw(Nlist2_obs, style = 'B'), nsim = 10000)

# compute by adding 4th order neighbors
moran4_I = moran(sealPolys_obs$Estimate, 
  nb2listw(Nlist4_obs, style = 'B'), 
  length(Nlist4_obs),
  Szero(nb2listw(Nlist4_obs, style = 'B')))
moran4_Itest_norm = moran.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist4_obs, style = 'B'),
  randomisation=FALSE)
moran4_Itest_rand = moran.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist4_obs, style = 'B'),
  randomisation=TRUE)
moran4_Itest_MC = moran.mc(sealPolys_obs$Estimate, 
  listw = nb2listw(Nlist4_obs, style = 'B'), nsim = 10000)

geary4_C = geary(sealPolys_obs$Estimate, 
  nb2listw(Nlist4_obs, style = 'B'), 
  length(Nlist4_obs), length(Nlist_obs) - 1,
  Szero(nb2listw(Nlist4_obs, style = 'B')))
geary4_Ctest_norm = geary.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist4_obs, style = 'B'),
  randomisation=FALSE)
geary4_Ctest_rand = geary.test(sealPolys_obs$Estimate, 
  nb2listw(Nlist4_obs, style = 'B'),
  randomisation=TRUE)
geary4_Ctest_MC = geary.mc(sealPolys_obs$Estimate, 
  listw = nb2listw(Nlist4_obs, style = 'B'), nsim = 10000)

#create a data.frame to use in a latex table
scheme = c('Order 1', 'Order 1 and 2', 'Order 1, 2, and 4')
moran_I = c(moran1_I$I, moran2_I$I, moran4_I$I)
moran_I_Pvals = c(
    paste0(formatC(moran1_Itest_norm$p.value,format = 'e', digits = 2),',', 
    formatC(moran1_Itest_MC$p.value,format = 'e', digits = 2)),
    paste0(formatC(moran2_Itest_norm$p.value,format = 'e', digits = 2),',', 
    formatC(moran2_Itest_MC$p.value,format = 'e', digits = 2)),
    paste0(formatC(moran4_Itest_norm$p.value,format = 'e', digits = 2),',', 
    formatC(moran4_Itest_MC$p.value,format = 'e', digits = 2))
  )
geary_C = c(geary1_C$C, geary2_C$C, geary4_C$C)
geary_C_Pvals = c(
    paste0(formatC(geary1_Ctest_norm$p.value,format = 'e', digits = 2),',', 
    formatC(geary1_Ctest_MC$p.value,format = 'e', digits = 2)),
    paste0(formatC(geary2_Ctest_norm$p.value,format = 'e', digits = 2),',', 
    formatC(geary2_Ctest_MC$p.value,format = 'e', digits = 2)),
    paste0(formatC(geary4_Ctest_norm$p.value,format = 'e', digits = 2),',', 
    formatC(geary4_Ctest_MC$p.value,format = 'e', digits = 2))
  )
moran_geary_all = data.frame(scheme = scheme, MoranI = moran_I, MIpvals = moran_I_Pvals,
  GearyC = geary_C, GCpvals = geary_C_Pvals)
moran_geary_all

# create a table to paste into latex
print(
    xtable(moran_geary_all, 
      align = c('l',rep('l', times = length(moran_geary_all[1,]))),
      digits = c(0, 0, 3, 0, 3, 0),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

# as an aside, here is how to compute Moran's I and Geary's C using the 
# neighborhood (or any weights) matrix, as a matrix operation

# compute deviations from mean
y_dev = sealPolys_obs$Estimate - mean(sealPolys_obs$Estimate)
#Moran's I computed using matrix operations
length(y_dev)*t(y_dev) %*% Nmat1_obs %*% y_dev/
  (sum(y_dev^2)*sum(Nmat1_obs))
#Moran's I from spdep function
moran(sealPolys_obs$Estimate, nb2listw(Nlist_obs, style = 'B'), length(Nlist_obs),
  Szero(nb2listw(Nlist_obs, style = 'B')))$I

# let v be our vector of values
v = sealPolys_obs$Estimate
# compute all pairwise differences as a matrix
dif2 = (outer(v, rep(1, times = length(v))) -
  outer(rep(1, times = length(v)), v))^2
# now compute Geary's C
(length(v)- 1)*sum(Nmat1_obs*dif2)/
  (2*sum(y_dev^2)*sum(Nmat1_obs))
#Geary's C from spdep function
geary(sealPolys_obs$Estimate, 
  nb2listw(Nlist_obs, style = 'B'), 
  length(Nlist_obs), length(Nlist_obs) - 1,
  Szero(nb2listw(Nlist_obs, style = 'B')))$C
