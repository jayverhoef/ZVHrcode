sec_path = 'Rcode/Chapter4/Section 4.4/figures/'
setwd(paste0(SLEDbook_path,sec_path))

file_name = "rowscolumns"

# Create plot of sites on an incomplete rectangular lattice (Figure 4.4)
pdf(paste0(file_name,'.pdf'), height = 5, width = 11)

	layout(matrix(1:3, nrow = 1))

		padj = -.1
		adj = -.2
		cex_pt = 2

		xpts <- kronecker(c(1,2,5,7,9), rep(1, times = 6))
		ypts <- kronecker(rep(1, times = 5), c(1,2,5,8,9,11))
		par(mar = c(1,5,5,1))
		plot(xpts, ypts, pch=16, bty="n", xaxt="n", yaxt="n", xlab="", ylab="",
			cex = cex_pt)
		segments(1,0,1,12)
		segments(2,0,2,12)
		segments(5,0,5,12)
		segments(7,0,7,12)
		segments(9,0,9,12)
		segments(0,1,10,1)
		segments(0,2,10,2)
		segments(0,5,10,5)
		segments(0,8,10,8)
		segments(0,9,10,9)
		segments(0,11,10,11)
		mtext('A', adj = adj, cex = 3, padj = padj)
		
		# Create plot of sites on an incomplete rectangular lattice (Figure 4.4)
		xpts <- c(1,1,2,2,2,5,5,7,9,9)
		ypts <- c(1,8,2,5,11,1,9,9,1,5)

		plot(xpts, ypts, pch=16, bty="n", xaxt="n", yaxt="n", xlab="", ylab="",
			cex = cex_pt)
		segments(1,0,1,12)
		segments(2,0,2,12)
		segments(5,0,5,12)
		segments(7,0,7,12)
		segments(9,0,9,12)
		segments(0,1,10,1)
		segments(0,2,10,2)
		segments(0,5,10,5)
		segments(0,8,10,8)
		segments(0,9,10,9)
		segments(0,11,10,11)
		mtext('B', adj = adj, cex = 3, padj = padj)

	# Create plot of sites on an incomplete rectangular lattice (Figure 4.4)
		xpts <- c(1,2,2,5,5,5,7,7,9,9)
		ypts <- c(8,2,1,11,8,5,9,2,11,5)
		plot(xpts, ypts, pch=16, bty="n", xaxt="n", yaxt="n", xlab="", ylab="",
			cex = cex_pt)
		segments(1,0,1,12)
		segments(2,0,2,12)
		segments(5,0,5,12)
		segments(7,0,7,12)
		segments(9,0,9,12)
		segments(0,1,10,1)
		segments(0,2,10,2)
		segments(0,5,10,5)
		segments(0,8,10,8)
		segments(0,9,10,9)
		segments(0,11,10,11)
		mtext('c', adj = adj, cex = 3, padj = padj)
		
	layout(1)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
