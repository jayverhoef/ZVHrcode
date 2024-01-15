sec_path = 'Rcode/Chapter3/Section 3.3/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(scatterplot3d)
library(sf)
source('addBreakBubbleLegend.R')

data(SO4obs)
SO4data = st_drop_geometry(SO4obs)
SO4data$xcoord = st_coordinates(SO4obs)[,1]
SO4data$ycoord = st_coordinates(SO4obs)[,2]

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           so4-3dscatterplot
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Create a 3-D scatterplot of the wet sulfate deposition data (Figure 3.1)
file_name = "figures/so4-3dscatterplot"
pdf(paste0(file_name,'.pdf'), width =  9, height = 9)

	scatterplot3d(SO4data$xcoord/1000, SO4data$ycoord/1000, 
		SO4data$SO4, pch=16, type="h", angle=120, grid=FALSE, box=FALSE, 
		xlab="Albers CONUS east (km)", ylab="", # mar = c(9,9,9,9),
		zlab="SO4 (kg/ha)", cex.lab = 1.5, cex.axis = 1.1)
	text(x = -5.1, y = 2.0, "Albers CONUS North (km)", srt = 120, cex = 1.4)
		
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
#           so4-bubbleplot
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Create bubble of the wet sulfate deposition data (Figure 3.2)
file_name = "figures/so4-bubbleplot"
data(USboundary)
# create categories from continuous data
cip = cut(SO4data$SO4, breaks = c(0, 10, 20, 30, 40, 50))
# scale the bubbles according to categories
cex = as.integer(cip)/2 + .5
# make the plot
pdf(paste0(file_name,'.pdf'), width = 12, height = 12)

  layout(matrix(c(1:2), nrow = 1), widths = c(3,1))
  plot(st_geometry(USboundary))
  points(SO4data[,c('xcoord', 'ycoord')], cex = cex)
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
