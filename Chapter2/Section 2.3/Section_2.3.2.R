sec_path = 'Rcode/Chapter2/Section 2.3/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#    Toy Example Box Polygon Layout (Figure 2.4)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Toy Example Box Polygon Layout (Figure 2.4)
file_name = "figures/toy5box"
pdf(paste0(file_name,'.pdf'), height = 7, width = 7)

		x <- c(2.5, 1.5, 2.5, 3.5, 3.5)
		y <- c(3.5, 2.5, 2.5, 2.5, 1.5)
		labs = c(1,2,3,4,5)
		par(mar = c(5,5,5,1))
		plot(c(1,4),c(1,4),type = 'n',xlab="",ylab="",xaxt = 'n', yaxt = 'n',
			bty = 'n')
		segments(x0=1,y0=2,x1=4,y1=2,lwd=3.0)
		segments(x0=1,y0=3,x1=4,y1=3,lwd=3.0)
		segments(x0=1,y0=2,x1=1,y1=3,lwd=3.0)
		segments(x0=2,y0=2,x1=2,y1=4,lwd=3.0)
		segments(x0=3,y0=1,x1=3,y1=4,lwd=3.0)
		segments(x0=4,y0=1,x1=4,y1=3,lwd=3.0)
		segments(x0=3,y0=1,x1=4,y1=1,lwd=3.0)
		segments(x0=2,y0=4,x1=3,y1=4,lwd=3.0)
		text(x,y,labels = labs, adj = .5, cex = 4)

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
#    Covariance Matrices Used in Text
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

BCM = matrix(0, nrow = 5, ncol = 5)
BCM[1,3] = BCM[3,1] = BCM[4,5] = BCM[5,4] = .5
BCM[2,3] = BCM[3,2] = BCM[3,4] = BCM[4,3] = .25
BCM
print(
    xtable(BCM, 
      align = c('c',rep('c', times = length(BCM[1,]))),
      digits = c(0,rep(2, times = 5))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

eigen(diag(5) - BCM)$values
eigen(BCM + diag(5))$values

Sig_CAR = solve(diag(5) - BCM)
Sig_CAR
Sig_CAR[lower.tri(Sig_CAR)] = NA
print(
    xtable(Sig_CAR, 
      align = c('c',rep('c', times = length(Sig_CAR[1,]))),
      digits = c(0,rep(2, times = 5))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

Sig_SAR = solve(diag(5) - BCM) %*% t(solve(diag(5) - BCM))
Sig_SAR
Sig_SAR[lower.tri(Sig_SAR)] = NA
print(
    xtable(Sig_SAR, 
      align = c('c',rep('c', times = length(Sig_SAR[1,]))),
      digits = c(0,rep(2, times = 5))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

Sig_SMA = (BCM + diag(5)) %*% (t(BCM) + diag(5))
Sig_SMA
Sig_SMA[lower.tri(Sig_SMA)] = NA
print(
    xtable(Sig_SMA, 
      align = c('c',rep('c', times = length(Sig_SMA[1,]))),
      digits = c(0,rep(2, times = 5))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

BCMall = rbind(cbind(BCM,c(0,0,0.5,0,0.25)),
	c(0,0,0.5,0,0.25,0))

print(
    xtable(BCMall, 
      align = c('c',rep('c', times = length(BCMall[1,]))),
      digits = c(0,rep(2, times = 6))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

Sig_CARall = solve(diag(6) - BCMall)
Sig_CARstar = Sig_CARall[1:5,1:5]
Sig_CARstar
Sig_CARstar[lower.tri(Sig_CARstar)] = NA
print(
    xtable(Sig_CARstar, 
      align = c('c',rep('c', times = length(Sig_CARstar[1,]))),
      digits = c(0,rep(2, times = 5))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

Sig_SARall = solve(diag(6) - BCMall) %*% t(solve(diag(6) - BCMall))
Sig_SARstar = Sig_SARall[1:5,1:5]
Sig_SARstar
Sig_SARstar[lower.tri(Sig_SARstar)] = NA
print(
    xtable(Sig_SARstar, 
      align = c('c',rep('c', times = length(Sig_SARstar[1,]))),
      digits = c(0,rep(2, times = 5))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

Sig_SMAall = (BCMall + diag(6)) %*% (t(BCMall) + diag(6))
Sig_SMAstar = Sig_SMAall[1:5,1:5]
Sig_SMAstar
Sig_SMAstar[lower.tri(Sig_SMAstar)] = NA
print(
    xtable(Sig_SMAstar, 
      align = c('c',rep('c', times = length(Sig_SMAstar[1,]))),
      digits = c(0,rep(2, times = 5))
    ),
    hline.after = NULL,
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
