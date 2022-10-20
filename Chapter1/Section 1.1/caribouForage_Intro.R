sec_path = 'Rcode/Chapter1/Section 1.1/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(classInt)
library(viridis)
source('addBreakColorLegend.R')
data(caribouDF)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Experimental Design for Caribou Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

cip = classIntervals(caribouDF$z, n = 4, style = 'fisher')
palp = viridis(5)[1:4]
cip_colors = findColours(cip, palp)

file_name = "figures/caribou_Design"
pdf(file = paste0(file_name,'.pdf'), width = 10, height = 8)
  layout(matrix(1:2, nrow = 1), widths = c(4,1))
  old.par = par(mar = c(5,5,1,1))
  plot(caribouDF$x,caribouDF$y, pch =15, cex = 12, 
    xlim = c(0.7,5.3), ylim = c(0.7,6.3),
    col = cip_colors, cex.lab = 2.5, cex.axis = 2, xlab = 'x-coordinate',
    ylab = 'y-coordinate')
  text(caribouDF$x,caribouDF$y,labels = caribouDF$tarp, pos = 3, cex = 2, col = 'white')
  text(caribouDF$x,caribouDF$y,labels = caribouDF$water, pos = 1, cex = 2, col = 'white')
  par(mar = c(0,0,0,0))
  plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
  addBreakColorLegend(xleft = 0.02, ybottom = .2, xright = .3, ytop = .8,
    breaks = cip$brks, colors = palp, cex = 2.5)
  par(old.par)
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
#           Boxplots of Caribou Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/caribou_boxplots"
pdf(file = paste0(file_name,'.pdf'), width = 10, height = 7)
  old.par = par(mar = c(2,6,1,1))
  boxplot(z ~ tarp + water, data = caribouDF, cex.lab = 2.3, cex.axis = 1.5,
    ylab = 'Percent Nitrogen in Prostrate Willows', pch = 19, cex = 2, lwd = 2)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))


