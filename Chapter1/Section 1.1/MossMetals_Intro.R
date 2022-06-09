sec_path = 'Rcode/Chapter1/Section 1.1/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(viridis)

# load data for graphics and analysis
# these are all sp objects
# outline of Alaska
data(AKboundary)
# outline of Cape Krusenstern National Preserve
data(CAKRboundary)
data(MOSSobs)
data(MOSSpreds)


# first plot Alaska, with Cape Krusenstern as red polygon, with big red
# box around it
file_name = "RedDog_Alaska"
pdf(file = paste0(file_name,'.pdf'), width = 11, height = 11)
  old.par = par(mar = c(0,0,0,0))
  plot(AKboundary, lwd = 2)
  plot(CAKRboundary, add = TRUE, col = 'red')
  lines(c(-450000, -450000, -370000, -370000, -450000),
    c(1900000, 2030000, 2030000, 1900000, 1900000), type = 'l',
    lwd = 5, col = 'red')
  #text(-1520000, 2200000, 'A', cex = 4)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

# make maps showing study area, data, and prediction locations
file_name = "RedDog_sampling"
pdf(file = paste0(file_name,'.pdf'), width = 15, height = 10)
  old.par = par(mar = c(0,1,1,0))
  layout(matrix(c(1,2), nrow = 1, ncol = 2, byrow = TRUE))
  # show the sampling locations in 
  plot(CAKRboundary, lwd = 2)
  plot(MOSSobs[MOSSobs$year == 2001,], add = TRUE, pch = 19, cex = .8, 
    col = 'darkorchid2')
  plot(MOSSobs[MOSSobs$year == 2006,], add = TRUE, pch = 19, cex = .8, 
    col = 'darkolivegreen3')
  text(-437345.1, 2012454, 'A', cex = 4)
  legend(-390000, 1990000, legend = c('2001', '2006'),
    pch = c(19, 19), col = c('darkorchid2','darkolivegreen3'),
    cex = 2)
  pal = viridis(7)
  par(mar = c(0,1,1,0))
  plot(CAKRboundary, lwd = 2)
  plot(MOSSpreds[MOSSpreds$strat==1,], add = TRUE, pch = 19, cex = .2, col=pal[1])
  plot(MOSSpreds[MOSSpreds$strat==2,], add = TRUE, pch = 19, cex = .2, col=pal[2])
  plot(MOSSpreds[MOSSpreds$strat==3,], add = TRUE, pch = 19, cex = .2, col=pal[3])
  plot(MOSSpreds[MOSSpreds$strat==4,], add = TRUE, pch = 19, cex = .5, col=pal[6])
  plot(MOSSpreds[MOSSpreds$strat==5,], add = TRUE, pch = 19, cex = 1, col=pal[7])
  text(-437345.1, 2012454, 'B', cex = 4)
  legend(-390000, 1990000, 
    legend = c('strat 1', 'strat 2', 'strat 3', 'strat 4', 'strat 5'),
    pch = c(19, 19), col = c(pal[1:3],pal[6:7]),
    cex = 2)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

MOSS_data = MOSSobs@data
MOSS_data01 = MOSS_data[MOSS_data$year == 2001,]
MOSS_data06 = MOSS_data[MOSS_data$year == 2006,]

# Show plot with concentration in relation to distance from road
file_name = "RedDog_DistRoad"
pdf(file = paste0(file_name,'.pdf'), width = 10, height = 10)
  old.par = par(mar = c(5,5,1,1))
  plot(
    c(-log(MOSS_data[MOSS_data$year == 2001 & 
        MOSS_data$sideroad == 'N','dist2road']),
      log(MOSS_data[MOSS_data$year == 2001 & 
        MOSS_data$sideroad == 'S','dist2road'])),
    c(log(MOSS_data[MOSS_data$year == 2001 & MOSS_data$sideroad == 'N','Zn']),
      log(MOSS_data[MOSS_data$year == 2001 & MOSS_data$sideroad == 'S','Zn'])),
    xlab = 'Log Distance From Road Towards South', 
    ylab = 'Log Zinc Concentration', 
    pch = 19, cex = 1.2, cex.lab = 2, cex.axis = 1.5
  )
  points(
    c(-log(MOSS_data[MOSS_data$year == 2006 & 
        MOSS_data$sideroad == 'N','dist2road']),
      log(MOSS_data[MOSS_data$year == 2006 & 
        MOSS_data$sideroad == 'S','dist2road'])),
    c(log(MOSS_data[MOSS_data$year == 2006 & MOSS_data$sideroad == 'N','Zn']),
      log(MOSS_data[MOSS_data$year == 2006 & MOSS_data$sideroad == 'S','Zn'])),
    pch = 3, cex = 2, lwd = 1.2
  )
  legend(6, 8.1, legend = c('2001', '2006'),
    pch = c(19, 3), cex = 2.5)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
