sec_path = 'Rcode/Chapter10/Section 10.2/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        R code for plots of the two "big example" sampling designs
#         (CPE-optimal and EK-optimal when rho=0.5 and nugget=0.5)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
	
# for Figure 10.2:

# CPElocations for the "big example" sampling design (output of B_ML.sas with rho=0.5 and nugget=.50):
CPElocations = read.table("CPEdesignlocations.txt")
colnames(CPElocations) = c('x', 'y')

# EKlocations for the "big example" sampling design (output of VPE_ML.sas for rho=0.5 and nugget=0.50):
EKlocations = read.table("trueEKdesignlocations.txt")
colnames(EKlocations) = c('x', 'y')

file_name = 'figures/CPE_EK_BigEx'
pdf(paste0(file_name,'.pdf'), width = 10, height = 5)

	layout(matrix(1:2, nrow = 1, byrow = T))

	padj = 0
	adj = -.05
	cex_mtext = 3
	cex_plot = .6

	par(mar = c(2, 2, 2, 2))
	plot(CPElocations, pch = 16, xaxt = "n", yaxt = "n", xlab = "", 
		ylab = "", bty = "n", cex = cex_plot)
	lines(c(0,0),c(0,1))
	lines(c(0,1),c(1,1))
	lines(c(1,1),c(1,0))
	lines(c(1,0),c(0,0))
	mtext('A', adj = adj, cex = cex_mtext, padj = padj)

	plot(EKlocations, pch = 16, xaxt = "n", yaxt = "n", xlab = "", 
		ylab = "", bty = "n", cex = cex_plot)
	lines(c(0,0),c(0,1))
	lines(c(0,1),c(1,1))
	lines(c(1,1),c(1,0))
	lines(c(1,0),c(0,0))
	mtext('B', adj = adj, cex = cex_mtext, padj = padj)

	layout(1)

dev.off()
		
system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))
