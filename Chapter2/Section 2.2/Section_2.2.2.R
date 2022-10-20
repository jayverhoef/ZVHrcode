sec_path = 'Rcode/Chapter2/Section 2.2/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#    Create plots to illustrate separability and additivity (Figure 2.1)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/sep_and_add"
pdf(paste0(file_name,'.pdf'), height = 7, width = 11)

	layout(matrix(1:2, nrow = 1))

		padj = -.5
		adj = -.2

		x <- c(1,1,3,3,5,5)
		y <- c(1.5,4.5,1.5,4.5,1.5,4.5)
		par(mar = c(5,5,5,1))
		plot(x,y,pch=16,xlim=c(0,6),ylim=c(0,6),xlab="",ylab="",cex=2.0,
			cex.axis = 2)
		segments(x0=1,y0=1.5,x1=3,y1=4.5,lty=1,lwd=3.0)
		segments(x0=1,y0=1.5,x1=5,y1=4.5,lty=5,lwd=3.0)
		segments(x0=3,y0=1.5,x1=5,y1=4.5,lty=3,lwd=5.0)
		segments(x0=1,y0=4.5,x1=3,y1=1.5,lty=1,lwd=3.0)
		segments(x0=1,y0=4.5,x1=5,y1=1.5,lty=5,lwd=3.0)
		segments(x0=3,y0=4.5,x1=5,y1=1.5,lty=3,lwd=5.0)
		points(x,y,pch=16,cex=2.0)
		x <- c(2,4,4)
		y <- c(1.5,1.5,4.5)
		mtext('A', adj = adj, cex = 3, padj = padj)

		plot(x,y,pch=16,xlim=c(0,6),ylim=c(0,6),xlab="",ylab="",cex=2.0,
			cex.axis = 2)
		segments(x0=2,y0=1.5,x1=4,y1=1.5,lty="dashed",lwd=3.0)
		segments(x0=4,y0=1.5,x1=4,y1=4.5,lty="dashed",lwd=3.0)
		segments(x0=2,y0=1.5,x1=4,y1=4.5,lty="solid",lwd=3.0)
		mtext('B', adj = adj, cex = 3, padj = padj)

	layout(1)

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
#    Create plots of loci of equal correlation for second-order stationary 
#                  covariance functions (Figure 2.2)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/spatialconfigs"
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


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#    Create plots of loci of equal correlation for second-order stationary
#                       covariance functions (Figure 2.3)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

x = c(1,1,2,2,5)
y = c(4,3,3,2,4)
labs = c(1,2,3,4,5)

file_name = "figures/locs4geocovfuns"
pdf(paste0(file_name,'.pdf'), height = 7, width = 7)
# Create plots of loci of equal correlation for second-order stationary covariance functions (Figure 2.3)

	par(mar = c(5,5,1,1))
	plot(x,y, xlim = c(1,5), ylim = c(1,5), pch = 19, cex = 2,
		cex.lab = 2.5, cex.axis = 1.8)
	text(x,y, labels = labs, pos = 3, cex = 2)

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
#                         Covariance Functions
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

C_1 = function(x,y){
	distmat = as.matrix(dist(cbind(x,y)))
	(1 - (3/8)*distmat + (1/128)*distmat^3)*(distmat <= 4)
}
C1mat = C_1(x,y)
C1mat[lower.tri(C1mat)] = NA
C1mat
print(
    xtable(C1mat, 
      align = c('c',rep('c', times = length(C1mat[1,]))),
      digits = c(0,rep(3, times = 5))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

C_2 = function(x,y){
	xdist = abs(outer(x, rep(1, times = length(x))) - 
		outer(rep(1, times = length(x)), x))
	ydist = abs(outer(y, rep(1, times = length(y))) - 
		outer(rep(1, times = length(y)), y))
	(1 - (1/5)*xdist)*(1 - (1/2)*ydist)*(xdist <= 5 & ydist <= 2)
}
C2mat = C_2(x,y)
C2mat[lower.tri(C2mat)] = NA
C2mat
print(
    xtable(C2mat, 
      align = c('c',rep('c', times = length(C2mat[1,]))),
      digits = c(0,rep(3, times = 5))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

C_3 = function(x,y){
	distmat = as.matrix(dist(cbind(x,y)))
  exp(-distmat/2) + 1*diag(length(x))
}
C3mat = C_3(x,y)
C3mat[lower.tri(C3mat)] = NA
C3mat
print(
    xtable(C3mat, 
      align = c('c',rep('c', times = length(C3mat[1,]))),
      digits = c(0,rep(3, times = 5))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
