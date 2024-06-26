sec_path = 'Rcode/Chapter12/Section 12.2/'
setwd(paste0(SLEDbook_path,sec_path))

library(sf)
library(viridis)
library(classInt)
library(prettymapr)
source('addBreakColorLegend.R')

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Stream Network Distances
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = 'figures/Stream_Dist'
pdf(paste0(file_name,'.pdf'), width = 12, height = 8)

	layout(matrix(1:2, nrow = 1))
	seglwd = 2
	ptcex = 2
	padj = -.0
	adj = -.0
	botlabcex = 1.9

	par(mar = c(0,1,4,1))
	plot(c(-1.5,1.5),c(0,4), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n',
		xlab = '', ylab = '')
	segments(x0 = 0, y0 = .5, x1 = 0, y1 = 2, lwd = 1.5*seglwd )
	segments(x0 = 0, y0 = 2, x1 = -1, y1 = 3, lwd = seglwd)
	segments(x0 = 0, y0 = 2, x1 = 1, y1 = 3, lwd = seglwd)
	segments(x0 = -1, y0 = 3, x1 = -1.5, y1 = 4, lwd = 0.5*seglwd )
	segments(x0 = -1, y0 = 3, x1 = -0.5, y1 = 4, lwd = 0.5*seglwd)
	segments(x0 = 1, y0 = 3, x1 = 0.5, y1 = 4, lwd = 0.5*seglwd)
	segments(x0 = 1, y0 = 3, x1 = 1.5, y1 = 4, lwd = 0.5*seglwd)
	arrows(x0 = .5, y0 = 1.75, x1 = .5, y1 = 0.75, angle = 15, lwd = 7)
	text(.6, 1.25, label = 'Flow', pos = 4, cex = 2)
	points(-0.5, 2.5, pch = 19, cex = ptcex)
	points(-0.75, 3.5, pch = 19, cex = ptcex)
	text(-0.48, 2.5, label = expression(italic(t)[j]), pos = 4, cex = 3)
	text(-0.71, 3.45, label = expression(italic(s)[i]), pos = 4, cex = 3)
	text(0, 0.1, label = expression(paste("Stream distance = ", 
		italic(s[i] - t[j]))), cex = botlabcex)
	mtext('A', adj = adj, cex = 4, padj = padj)


	plot(c(-1.5,1.5),c(0,4), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n',
		xlab = '', ylab = '')
	segments(x0 = 0, y0 = .5, x1 = 0, y1 = 2, lwd = 1.5*seglwd )
	segments(x0 = 0, y0 = 2, x1 = -1, y1 = 3, lwd = seglwd)
	segments(x0 = 0, y0 = 2, x1 = 1, y1 = 3, lwd = seglwd)
	segments(x0 = -1, y0 = 3, x1 = -1.5, y1 = 4, lwd = 0.5*seglwd )
	segments(x0 = -1, y0 = 3, x1 = -0.5, y1 = 4, lwd = 0.5*seglwd)
	segments(x0 = 1, y0 = 3, x1 = 0.5, y1 = 4, lwd = 0.5*seglwd)
	segments(x0 = 1, y0 = 3, x1 = 1.5, y1 = 4, lwd = 0.5*seglwd)
	arrows(x0 = .5, y0 = 1.75, x1 = .5, y1 = 0.75, angle = 15, lwd = 7)
	text(.6, 1.25, label = 'Flow', pos = 4, cex = 2)
	points(0.5, 2.5, pch = 19, cex = ptcex)
	points(-0.75, 3.5, pch = 19, cex = ptcex)
	points(0, 2, pch = 19, cex = ptcex)
	text(0.18, 2.6, label = expression(italic(t)[j]), pos = 4, cex = 3)
	text(-0.71, 3.45, label = expression(italic(s)[i]), pos = 4, cex = 3)
	text(-0.45, 1.89, label = expression(italic(q)[ij]), pos = 4, cex = 3)
	text(0, 0.1, label = expression(paste("Stream distance = ", 
		italic((s[i] - q[ij]) + (t[j] -q[ij])))), cex = botlabcex)
	mtext('B', adj = adj, cex = 4, padj = padj)

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
#                  Salt River Map
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

setwd(paste0(SLEDbook_path,sec_path,'SaltRiver_BlockKrige/',
	'TroutDensity_BlockKrige.ssn'))
strms <- st_read("edges.shp")
sites <- st_read("sites.shp")
setwd(paste0(SLEDbook_path,sec_path))

file_name = 'figures/SaltRiver_Map'
pdf(paste0(file_name,'.pdf'), width = 10, height = 10)

par(mar = c(0,0,0,0))
plot(st_geometry(strms), lwd = 1.5 + strms$afvArea*5, col = 'lightblue2')
brks = c(0, 6.000001, 11.000001, 19.000001, 28.000001, 39.000001,
	64.000001, 132.000001)
cip = classIntervals(sites$trout_100m, style = 'fixed', fixedBreaks = brks)
palp = viridis(7)
cip_colors = findColours(cip, palp)
plot(st_geometry(sites), add = TRUE, pch = 19, col = cip_colors, cex = 1.5)
addBreakColorLegend(xleft = 1752000, ybottom = 1424484, xright = 1757000, ytop = 1445000, breaks = brks, colors = palp, cex = 1.5, printFormat = "4.0")
text(1745000, 1448000, 'Trout per 100m', pos = 4, cex = 2)
addnortharrow(pos = 'bottomleft', padin = c(1.9, 3.5))
addscalebar(pos = 'bottomleft', widthhint = 0.125, padin = c(1.7, 3.2),
	label.cex = 2)

dev.off()
		
system(paste0('pdfcrop ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
	sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
	sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
		sec_path,file_name,'-crop.pdf','\''))

