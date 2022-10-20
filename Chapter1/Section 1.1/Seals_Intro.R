sec_path = 'Rcode/Chapter1/Section 1.1/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(shape)
library(maps)
library(spdep)
library(sf)
source('addBreakColorLegend.R')
data(sealPolys)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Seal Stocks in Southeast Alaska
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/seal_Stocks"
pdf(file = paste0(file_name,'.pdf'), width = 11, height = 11)
  colPalStock = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3')
  plot(st_geometry(sealPolys), col = colPalStock[sealPolys$stockid - 7])
  text(980000,1085000,'8', cex = 6)
  text(1123764,1130000,'9', cex = 6)
  text(1134476,948464,'10', cex = 6)
  text(1230000,810000,'11', cex = 6)
  text(1380000,930000,'12', cex = 6)
  par(usr = c(-268,-128,19.4,73))
  map("world", c("USA:Alaska"), add = TRUE, lwd = 2)
  lines(c(-138,-129),c(59.5,59.5), lwd = 3, col = 'red')
  lines(c(-138,-129),c(54.5,54.5), lwd = 3, col = 'red')
  lines(c(-138,-138),c(54.5,59.5), lwd = 3, col = 'red')
  lines(c(-129,-129),c(54.5,59.5), lwd = 3, col = 'red')
  Arrows(x0 = -245, y0 = 28, x1 = -245, y1 = 38, lwd = 6)
  text(-245, 42, labels = 'N', cex = 4)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

#### as a png file

png(file = paste0(file_name,'.png'), width = 960, height = 960)
  colPalStock = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3')
  plot(st_geometry(sealPolys), col = colPalStock[sealPolys$stockid - 7])
  text(980000,1085000,'8', cex = 6)
  text(1123764,1130000,'9', cex = 6)
  text(1134476,948464,'10', cex = 6)
  text(1230000,810000,'11', cex = 6)
  text(1380000,930000,'12', cex = 6)
  par(usr = c(-268,-128,19.4,73))
  map("world", c("USA:Alaska"), add = TRUE, lwd = 2)
  lines(c(-138,-129),c(59.5,59.5), lwd = 3, col = 'red')
  lines(c(-138,-129),c(54.5,54.5), lwd = 3, col = 'red')
  lines(c(-138,-138),c(54.5,59.5), lwd = 3, col = 'red')
  lines(c(-129,-129),c(54.5,59.5), lwd = 3, col = 'red')
  Arrows(x0 = -245, y0 = 28, x1 = -245, y1 = 38, lwd = 6)
  text(-245, 42, labels = 'N', cex = 4)
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Trends in Seal Stocks
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/seal_Trends"
pdf(file = paste0(file_name,'.pdf'), width = 11, height = 11)
  col.pal7 = c("#440154FF", "#443A83FF", "#31688EFF", "#21908CFF", "#35B779FF", 
    "#8FD744FF", "#FDE725FF")
  missingCol = 'grey80'
  coords = st_drop_geometry(st_centroid(sealPolys)[,c('x','y')])
  brks7 = c(-5, -.15, -.07, -.02, .02, .07, .15, 5)
  #  par(bg = 'grey70')
  par(mar = c(0,0,0,0))
  plot(st_geometry(sealPolys), col = missingCol)
  points(coords, pch = 19, cex = 2, col = missingCol)
    #  xlim = c(896054, 1498427),ylim = c(712289, 1216913))
  plot(st_geometry(sealPolys), add = TRUE, col = col.pal7[as.integer(
		cut(sealPolys$Estimate, breaks = brks7))])
  points(coords, pch = 19, cex = 2, 
    col = col.pal7[as.integer(cut(sealPolys$Estimate, breaks = brks7))]) 
  brks7[1] = min(sealPolys$Estimate, na.rm = TRUE) - .00001
  brks7[8] = max(sealPolys$Estimate, na.rm = TRUE)
  addBreakColorLegend(1304916, 981392, 1357619, 1198651, brks7, 
    colors = col.pal7, cex = 1.5, printFormat = "4.3")
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
  
#### as a png file

png(file = paste0(file_name,'.png'), width = 960, height = 960)
  col.pal7 = c("#440154FF", "#443A83FF", "#31688EFF", "#21908CFF", "#35B779FF", 
    "#8FD744FF", "#FDE725FF")
  missingCol = 'grey80'
  coords = st_drop_geometry(st_centroid(sealPolys)[,c('x','y')])
  brks7 = c(-5, -.15, -.07, -.02, .02, .07, .15, 5)
  #  par(bg = 'grey70')
  par(mar = c(0,0,0,0))
  plot(st_geometry(sealPolys), col = missingCol)
  points(coords, pch = 19, cex = 2, col = missingCol)
    #  xlim = c(896054, 1498427),ylim = c(712289, 1216913))
  plot(st_geometry(sealPolys), add = TRUE, col = col.pal7[as.integer(
		cut(sealPolys$Estimate, breaks = brks7))])
  points(coords, pch = 19, cex = 2, 
    col = col.pal7[as.integer(cut(sealPolys$Estimate, breaks = brks7))]) 
  brks7[1] = min(sealPolys$Estimate, na.rm = TRUE) - .00001
  brks7[8] = max(sealPolys$Estimate, na.rm = TRUE)
  addBreakColorLegend(1304916, 981392, 1357619, 1198651, brks7, 
    colors = col.pal7, cex = 1.5, printFormat = "4.3")
dev.off()

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Neighbors in Seal Stocks
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

nTot = length(sealPolys)
nObs = sum(!is.na(sealPolys$Estimate))
nMiss = nTot - nObs

source('Neighmat.R')
Nlist = poly2nb(sealPolys, snap = 2000)
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
attr(Nlist,'polyid') = as.factor(as.character(sealPolys@data$polyid))
attr(Nlist,'stockid') = as.factor(as.character(sealPolys@data$stockid))
num = lapply(Nlist, function(x) length(x))
num = unlist(num)
Nmat = Neighmat(Nlist, num, length(num))
Nmat1 = pmax(Nmat,t(Nmat))
Nmat2 = (Nmat1 %*% Nmat1 > 0 | Nmat1 > 0)*1
diag(Nmat2) = 0
Nmat4 = (Nmat2 %*% Nmat2 > 0 | Nmat2 > 0)*1
diag(Nmat4) = 0
xleft = 1268000
xright = 1390000
ybottom = 738000
ytop = 800000
xexp = (xright - xleft)*.06
yexp = (ytop - ybottom)*.35
  
file_name = "figures/seal_Neighbors"
pdf(file = paste0(file_name,'.pdf'), width = 11, height = 11)
  layout(matrix(c(1,1,2,1,1,1,3,1,1), nrow = 3, byrow = TRUE))
  old.par = par(mar = c(0,0,0,0))
  plot(st_geometry(sealPolys))
  rect(xleft - xexp, ybottom - yexp, xright + xexp, ytop + yexp, col = 'grey80', 
  border = NA)
  plot(st_geometry(sealPolys), add = TRUE)
  plot(Nlist, coords, add = TRUE, lwd = 2)
  text(920000, 1190000, labels = 'A', cex = 6)
  Nlist2 = apply(Nmat2, 1, function(x) as.integer(which(x > 0)))
  class(Nlist2) = 'nb'
  attr(Nlist2,'type') = 'queen'
  attr(Nlist2,'sym') = TRUE
  attr(Nlist2,'polyid') = attr(Nlist,'polyid')
  attr(Nlist2,'stockid') = attr(Nlist,'polyid')
  plot(st_geometry(sealPolys), xlim = c(xleft, xright), 
		ylim = c(ybottom, ytop))
  plot(Nlist2, coords, add = TRUE, lwd = 1)
  mtext('B', adj = -.15, padj = 1.1, cex = 4)
  Nlist4 = apply(Nmat4, 1, function(x) as.integer(which(x > 0)))
  class(Nlist4) = 'nb'
  attr(Nlist4,'type') = 'queen'
  attr(Nlist4,'sym') = TRUE
  attr(Nlist4,'polyid') = attr(Nlist,'polyid')
  attr(Nlist4,'stockid') = attr(Nlist,'polyid')
  plot(st_geometry(sealPolys), xlim = c(xleft, xright), ylim = c(ybottom, ytop))
  plot(Nlist4, coords, add = TRUE, lwd = .5)
  mtext('C', adj = 0, padj = -.4, cex = 4)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
  
#### as a png file

png(file = paste0(file_name,'.png'), width = 960, height = 960)
  layout(matrix(c(1,1,2,1,1,1,3,1,1), nrow = 3, byrow = TRUE))
  old.par = par(mar = c(0,0,0,0))
  plot(st_geometry(sealPolys))
  rect(xleft - xexp, ybottom - yexp, xright + xexp, ytop + yexp, col = 'grey80', 
  border = NA)
  plot(st_geometry(sealPolys), add = TRUE)
  plot(Nlist, coords, add = TRUE, lwd = 2)
  text(920000, 1190000, labels = 'A', cex = 6)
  Nlist2 = apply(Nmat2, 1, function(x) as.integer(which(x > 0)))
  class(Nlist2) = 'nb'
  attr(Nlist2,'type') = 'queen'
  attr(Nlist2,'sym') = TRUE
  attr(Nlist2,'polyid') = attr(Nlist,'polyid')
  attr(Nlist2,'stockid') = attr(Nlist,'polyid')
  plot(st_geometry(sealPolys), xlim = c(xleft, xright), 
		ylim = c(ybottom, ytop))
  plot(Nlist2, coords, add = TRUE, lwd = 1)
  mtext('B', adj = -.15, padj = 1.1, cex = 4)
  Nlist4 = apply(Nmat4, 1, function(x) as.integer(which(x > 0)))
  class(Nlist4) = 'nb'
  attr(Nlist4,'type') = 'queen'
  attr(Nlist4,'sym') = TRUE
  attr(Nlist4,'polyid') = attr(Nlist,'polyid')
  attr(Nlist4,'stockid') = attr(Nlist,'polyid')
  plot(st_geometry(sealPolys), xlim = c(xleft, xright), ylim = c(ybottom, ytop))
  plot(Nlist4, coords, add = TRUE, lwd = .5)
  mtext('C', adj = 0, padj = -.4, cex = 4)
  par(old.par)
dev.off()

