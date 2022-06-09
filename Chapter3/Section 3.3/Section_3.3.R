sec_path = 'Rcode/Chapter3/Section 3.3/'
setwd(paste0(SLEDbook_path,sec_path))

################################################################################
#-------------------------------------------------------------------------------
#                            so4-3dscatterplot
#-------------------------------------------------------------------------------
################################################################################
# Create a 3-D scatterplot of the wet sulfate deposition data (Figure 3.1)
file_name = "so4-3dscatterplot"
library(ZVHdata)
library(scatterplot3d)
data(SO4obs)
pdf(paste0(file_name,'.pdf'))
scatterplot3d(SO4obs@coords[,1]/1000, SO4obs@coords[,2]/1000, SO4obs@data$SO4,
  pch=16, type="h", angle=120, grid=FALSE, box=FALSE, 
  xlab="Albers CONUS North (km)", ylab="Albers CONUS West (km)", zlab="SO4", 
  cex.lab = 1.5, cex.axis = 1.1)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

################################################################################
#-------------------------------------------------------------------------------
#                            so4-bubbleplot
#-------------------------------------------------------------------------------
################################################################################

# Create bubble of the wet sulfate deposition data (Figure 3.2)
file_name = "so4-bubbleplot"
source('addBreakBubbleLegend.R')
library(sp)
data(USboundary)
# create categories from continuous data
cip = cut(SO4obs@data$SO4, breaks = c(0, 10, 20, 30, 40, 50))
# scale the bubbles according to categories
cex = as.integer(cip)/2 + .5
# make the plot
pdf(paste0(file_name,'.pdf'), width = 12, height = 12)
  layout(matrix(c(1:2), nrow = 1), widths = c(3,1))
  plot(USboundary)
  points(SO4obs, cex = cex)
  addBreakBubbleLegend(ybottom = .4, ytop = .6, 
    cexes = sort(unique(cex)), bubble_labels = levels(cip), cex = 1.7,
    header = 'SO4 (kg/ha)', cex.header = 2, x.header = 0)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
