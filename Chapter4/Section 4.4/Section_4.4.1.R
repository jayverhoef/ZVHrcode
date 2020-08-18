sec_path = 'Rcode/Chapter4/Section 4.4/'
setwd(paste0(SLEDbook_path,sec_path))

file_name = "corr-not-trend"
# Create plot of simulated data on transect with strong positive correlation and no trend (Figure 4.3)
pdf(paste0(file_name,'.pdf'), height = 7, width = 9)
  y <- matrix(0,20,1)
  set.seed(1008)
  for(i in 2:20){
    y[i,1] <- 0.95*y[i-1,1]+rnorm(1,0,0.1)
  }
  i <- c(1:20)
  old.par = par(mar = c(4,4.5,.5,.5))
  plot(i, y, pch = 19, cex = 2.5, cex.lab = 2, cex.axis = 1.5)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

