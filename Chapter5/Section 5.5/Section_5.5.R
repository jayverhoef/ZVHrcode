sec_path = 'Rcode/Chapter5/Section 5.5/'
setwd(paste0(SLEDbook_path,sec_path))

file_name = "spatialconfigs"
# Create plot of possibly non-orthogonal rectangular grids (Figure 5.1)
x1 <- c(1,5,12,14,3,7,14,16,4,8,15,17,6,10,17,19,9,13,20,22)
y1 <- c(1,1,1,1,4,4,4,4,6,6,6,6,9,9,9,9,14,14,14,14)
x2 <- c(1,4,7,10,2,5,8,11,3,6,9,12,4,7,10,13,5,8,11,14)
y2 <- c(1,1,1,1,4,4,4,4,7,7,7,7,10,10,10,10,13,13,13,13)
x3 <- c(1,5,8,10,1,5,8,10,1,5,8,10,1,5,8,10,1,5,8,10)
y3 <- c(1,1,1,1,4,4,4,4,5,5,5,5,8,8,8,8,13,13,13,13)
pdf(paste0(file_name,'.pdf'))
  old.par = par(mfrow=c(1,3),bty="n",pty="s", mar = c(.5,3,5,3))
  cex = 2
  lwd = 2
  padj = -.1
  adj = -.2
  plot(x1,y1,pch=4,xlab="",ylab="",xaxt="n",yaxt="n", cex = cex, lwd = lwd)
  mtext('A', adj = adj, cex = 3, padj = padj)
  plot(x2,y2,pch=4,xlab="",ylab="",xaxt="n",yaxt="n", cex = cex, lwd = lwd)
  mtext('B', adj = adj, cex = 3, padj = padj)
  plot(x3,y3,pch=4,xlab="",ylab="",xaxt="n",yaxt="n", cex = cex, lwd = lwd)
  mtext('C', adj = adj, cex = 3, padj = padj)
  par(old.par)
dev.off()
system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
