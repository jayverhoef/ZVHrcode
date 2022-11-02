sec_path = 'Rcode/Chapter10/Section 10.2/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        R code for displaying model-based sampling designs
#                  (toy examples) in Figure 10.1
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# R code for displaying model-based sampling designs 
x = rep(seq(0.1, 0.9, 0.2), 5)
y = rep(seq(0.9, 0.1, -0.2), each = 5)
z = cbind(x, y)

file_name = 'figures/LocOptToy5'
pdf(paste0(file_name,'.pdf'), width = 7, height = 5)

	layout(matrix(1:6, nrow = 2, byrow = T))

	par(pty="s",bty="n",xaxt="n",yaxt="n", mar = c(3,4,5,1))
	padj = -.4
	adj = -.2
	cex_mtext = 2.5
	cex_plot = 2.1

	# A
	black1 <- t(matrix(c(z[1,], z[5,], z[13,], z[21,], z[25,]), nrow = 2))
	plot(z, axes=F, xlab = "", ylab = "", pch = 1, cex = cex_plot)
	matpoints(black1[,1], black1[,2], cex = cex_plot, pch = 16)
	mtext('A', adj = adj, cex = cex_mtext, padj = padj)

	# B
	black2 <- t(matrix(c(z[1,],z[3,],z[5,],z[21,],z[25,]),nrow=2))
	plot(z, axes=F, xlab = "", ylab = "", pch = 1, cex = cex_plot)
	matpoints(black2[,1], black2[,2], cex = cex_plot, pch = 16)
	mtext('B', adj = adj, cex = cex_mtext, padj = padj)

	# C
	black3 <- t(matrix(c(z[2,],z[10,],z[13,],z[16,],z[24,]),nrow=2))
	plot(z, axes=F, xlab = "", ylab = "", pch = 1, cex = cex_plot)
	matpoints(black3[,1], black3[,2], cex = cex_plot, pch = 16)
	mtext('C', adj = adj, cex = cex_mtext, padj = padj)

	# D
	black4 <- t(matrix(c(z[3,],z[8,],z[13,],z[18,],z[23,]),nrow=2))
	plot(z, axes=F, xlab = "", ylab = "", pch = 1, cex = cex_plot)
	matpoints(black4[,1], black4[,2], cex = cex_plot, pch = 16)
	mtext('D', adj = adj, cex = cex_mtext, padj = padj)

	# E
	black5 <- t(matrix(c(z[1,],z[2,],z[3,],z[24,],z[25,]),nrow=2))
	plot(z, axes=F, xlab = "", ylab = "", pch = 1, cex = cex_plot)
	matpoints(black5[,1], black5[,2], cex = cex_plot, pch = 16)
	mtext('E', adj = adj, cex = cex_mtext, padj = padj)

	# F
	black6 <- t(matrix(c(z[1,],z[2,],z[14,],z[19,],z[24,]),nrow=2))
	plot(z, axes=F, xlab = "", ylab = "", pch = 1, cex = cex_plot)
	matpoints(black6[,1], black6[,2], cex = cex_plot, pch = 16)
	mtext('F', adj = adj, cex = cex_mtext, padj = padj)

	layout(1)

dev.off()
		
system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

