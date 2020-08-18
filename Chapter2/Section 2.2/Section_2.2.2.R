sec_path = 'Rcode/Chapter2/Section 2.2/'
setwd(paste0(SLEDbook_path,sec_path))

file_name = "spatialconfigs"
# Create plots of loci of equal correlation for second-order stationary covariance functions (Figure 2.1)
library(plotrix)
pdf(paste0(file_name,'.pdf'), height = 7, width = 10)
  padj = -.5
  adj = -.2
  cex = 2.5
  old.par = par(mfrow = c(2,3),pty = "s", mar = c(3,1,5,1))
  plot(c(3,-3,3,-3), c(2,2,-2,-2), pch = 16, xlim = c(-4,4), ylim = c(-4,4), 
    xlab = "", ylab = "", cex.axis = 2, cex = cex)
  points(0, 0, pch = 4, cex = cex, lwd = 4)
  mtext('A', adj = adj, cex = 3, padj = padj)
  plot(c(-3,-2,2,3), c(-2,-3,3,2), pch = 16, xlim = c(-4,4), ylim = c(-4,4),
    xlab = "", ylab = "", cex.axis = 2, cex = cex)
  points(0, 0, pch = 4, cex = cex, lwd = 4)
  mtext('B', adj = adj, cex = 3, padj = padj)
  plot(c(-3,-3,-2,-2,2,2,3,3), c(-2,2,-3,3,-3,3,-2,2), pch = 16, 
    xlim = c(-4,4), ylim = c(-4,4), xlab = "", ylab = "", cex.axis = 2, cex = cex)
  points(0, 0, pch = 4, cex = cex, lwd = 4)
  mtext('C', adj = adj, cex = 3, padj = padj)
  plot(0, 0, pch = 4, xlim = c(-4,4), ylim = c(-4,4), xlab="", ylab="",
    cex.axis = 2, cex = cex, lwd = 4)
  draw.circle(0, 0, 3, lwd = 4)
  mtext('D', adj = adj, cex = 3, padj = padj)
  plot(0, 0, pch = 4, xlim = c(-4,4), ylim = c(-4,4), xlab="", ylab="",
    cex.axis = 2, cex = cex, lwd = 4)
  draw.ellipse( 0, 0, sqrt(13), 2, 45, lwd = 4)
  mtext('E', adj = adj, cex = 3, padj = padj)
  plot(0, 0, pch = 4, xlim = c(-4,4), ylim = c(-4,4), xlab="", ylab="",
    cex.axis = 2, cex = cex, lwd = 4)
  polygon(c(3,0,-3,0), c(0,3,0,-3), lwd = 4)
  mtext('F', adj = adj, cex = 3, padj = padj)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
