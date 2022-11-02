sec_path = 'Rcode/Chapter10/Section 10.3/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#          Probability-based sampling designs for Figure 10.3
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Probability-based sampling designs

file_name = 'figures/Prob-based_designs'
pdf(paste0(file_name,'.pdf'), width = 10, height = 5)

	layout(matrix(1:3, nrow = 1, byrow = T))

	padj = 0
	adj = -.1
	cex_mtext = 3
	cex_plot = 1.5

	par(pty = "s", xaxt = "n", yaxt = "n", bty = "n", mar = c(2,2,2,2))

	# A
	set.seed(1001)
	x <- runif(25)
	y <- runif(25)
	plot(x, y, xlab = "", ylab = "", xlim = c(0, 1), 
		ylim = c(0, 1), pch = 19, cex = cex_plot)
	segments(c(0, 0, 1, 1), c(0, 1, 1, 0), c(0, 1, 1, 0), c(1, 1, 0, 0))
	mtext('A', adj = adj, cex = cex_mtext, padj = padj)
	
	#B
	x0 <- runif(1, min = 0, max = 0.2)
	y0 <- runif(1, min = 0, max = 0.2)
	x00 <- c(x0, x0 + 0.2, x0 + 0.4, x0 + 0.6, x0 + 0.8)
	y00 <- c(y0, y0+ 0.2, y0 + 0.4, y0 + 0.6, y0 + 0.8)
	x <- rep(x00, 5)
	y <- rep(y00, each = 5)
	plot(x, y, xlab = "", ylab = "", xlim = c(0, 1), 
		ylim = c(0, 1), pch = 19, cex = cex_plot)
	segments(c(0, 0, 1, 1), c(0, 1, 1, 0), c(0, 1, 1, 0), c(1, 1, 0, 0))
	mtext('B', adj = adj, cex = cex_mtext, padj = padj)
	
	#C
	x <- matrix(0, 25, 1)
	y <- matrix(0, 25, 1)
	for(i in 1:5){
	 for(j in 1:5){
		x[(i - 1)*5 + j, 1] <- runif(1, min = 0.2*(i - 1), max = 0.2*i)
		y[(i - 1)*5 + j, 1] <- runif(1, min = 0.2*(j - 1), max = 0.2*j)
	}}
	plot(x, y, xlab = "", ylab = "", xlim = c(0, 1), 
		ylim = c(0, 1), pch = 19, cex = cex_plot)
	segments(c(0, 0, 1, 1), c(0, 1, 1, 0), c(0, 1, 1, 0), c(1, 1, 0, 0))
	segments(c(0, 0, 0, 0), c(0.2, 0.4, 0.6, 0.8), c(1, 1, 1, 1), 
		c(0.2, 0.4, 0.6, 0.8))
	segments(c(0.2, 0.4, 0.6, 0.8), c(0, 0, 0, 0), c(0.2, 0.4, 0.6, 0.8), 
		c(1, 1, 1, 1))
	mtext('C', adj = adj, cex = cex_mtext, padj = padj)

	layout(1)

dev.off()
		
system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

