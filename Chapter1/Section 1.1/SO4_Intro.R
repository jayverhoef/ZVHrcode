sec_path = 'Rcode/Chapter1/Section 1.1/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(classInt)
library(viridis)
library(colorspace)

# load data from package
data(SO4obs)
data(USboundary)

# Make traditional plots of predictions and standard errors
cip = classIntervals(SO4obs@data$SO4, n = 6, style = 'fisher')
palp = viridis(6)
cip_colors = findColours(cip, palp)

file_name = "SO4_Intro"
source('addBreakColorLegend.R')
pdf(file = paste0(file_name,'.pdf'), width = 11, height = 8.5)
  old.par = par(mar = c(0,0,5,0))
  layout(matrix(1:2, nrow = 1, byrow = TRUE), widths = c(3,1))
  plot(SO4obs, col = cip_colors, pch = 19, cex = 1.5)
  plot(USboundary, add = TRUE, border = 'black')
  par(mar = c(0,0,0,0))
  plot(c(0,1),c(0,1), type = 'n', xaxt = 'n', yaxt = 'n',
    xlab = '', ylab = '', bty = 'n')
  addBreakColorLegend(xleft = 0, ybottom = .2, xright = .2, ytop = .7,
    breaks = cip$brks, colors = palp, cex = 1.5)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
