sec_path = 'Rcode/Chapter10/Section 10.4/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Augmented space-filling designs for Figure 10.5
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = 'figures/augment_designs'
pdf(paste0(file_name,'.pdf'), width = 10, height = 5)


	layout(matrix(1:3, nrow = 1, byrow = T))

	padj = 0
	adj = -.1
	cex_mtext = 3
	cex_plot = 1.0

	par(pty = "s", xaxt = "n", yaxt = "n", bty="n")
	
# A
	set.seed(3141)
	x0 <- c(1/14, 3/14, 5/14, 7/14, 9/14, 11/14, 13/14)
	y0 <- x0
	x <- rep(x0, 7)
	y <- rep(y0, each = 7)
	xsample <- sample(x, 15)
	ysample <- sample(y, 15)
	xdell <- matrix(0, 30, 1)
	ydell <- xdell
	for(i in 1:30){
		xdel <- runif(1, -1/14, 1/14)
		ydel <- runif(1, -1/14, 1/14)
		if(sqrt(xdel^2 + ydel^2) < 1/14){
			xdell[i,1] <- xdel
			ydell[i,1] <- ydel
		}
	}
	xnew <- xsample + xdell[1:15,1]
	ynew <- ysample + ydell[1:15,1]
	xall <- cbind(x, xnew)
	yall <- cbind(y, ynew)
	plot(xall, yall, xlab = "", ylab = "", xlim = c(0,1), 
		ylim = c(0,1), pch = 19, cex = cex_plot)
	segments(c(0, 0, 1, 1), c(0, 1, 1, 0), c(0, 1, 1, 0), c(1, 1, 0, 0))
	mtext('A', adj = adj, cex = cex_mtext, padj = padj)

# B
	set.seed(3143)
	x00 <- c(1/14, 3/14, 5/14, 7/14, 9/14, 11/14)
	y00 <- x00
	x6 <- rep(x00, 6)
	y6 <- rep(y00, each = 6)
	xsample2 <- sample(x6, 2)
	ysample2 <- sample(y6, 2)
	plot(x, y, xlab = "", ylab = "", xlim = c(0, 1), 
		ylim = c(0, 1), pch = 16, cex = cex_plot)
	segments(c(0, 0, 1, 1), c(0, 1, 1, 0), c(0, 1, 1, 0), c(1, 1, 0, 0))
	x000 <- c(0,(1/3)*(1/7),(2/3)*(1/7),1/7)
	y000 <- c(0,(1/3)*(1/7),(2/3)*(1/7),1/7)
	xdelta <- rep(x000, 4)
	ydelta <- rep(y000, each = 4)
	points(xsample2[1] + xdelta, ysample2[1] + ydelta, pch = 16)
	points(xsample2[2] + xdelta, ysample2[2] + ydelta, pch = 16)
	mtext('B', adj = adj, cex = cex_mtext, padj = padj)

# C
# library(spatstat)
	set.seed(2141)
	xyssipp <- rSSI(1/9, 49)
	xssisample <- sample(xyssipp$x, 15)
	yssisample <- sample(xyssipp$y, 15)
	for(i in 1:30){
		xssidel <- runif(1, -1/14, 1/14)
		yssidel <- runif(1, -1/14, 1/14)
		xssidell <- matrix(0, 30, 1)
		yssidell <- xssidell
		if(sqrt(xssidel^2 + yssidel^2) < 1/14){
			xssidell[i,1] <- xssidel
			yssidell[i,1] <- yssidel
		}
	}
	xssinew <- xssisample + xssidell[1:15,1]
	yssinew <- yssisample + yssidell[1:15,1]
	xssiall <- c(xyssipp$x, xssinew)
	yssiall <- c(xyssipp$y, yssinew)
	plot(xssiall, yssiall, xlab = "",ylab = "", xlim = c(0, 1), 
		ylim = c(0, 1), pch = 16, cex = cex_plot)
	segments(c(0, 0, 1, 1), c(0, 1, 1, 0), c(0, 1, 1, 0), c(1, 1, 0, 0))
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
