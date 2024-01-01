sec_path = 'Rcode/Chapter4/Section 4.2/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Toy Figure for OLS and GLS 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

file_name = "figures/toyolsgls"

# Create plot of sites for the toy example of Section 4.2 (Figure 4.1)
x <- c(1,1,2,2,5)
y <- c(4,3,3,2,4)
labs = as.character(1:5)

pdf(paste0(file_name,'.pdf'))

	par(mar = c(5,5,1,1))
	plot(x,y, xlim = c(1,5), ylim = c(1,5), pch = 19, cex = 2,
		cex.lab = 2.5, cex.axis = 1.8, xlab = '', ylab = '')
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
#              Compute some results in Section 4.2 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Compute some results in Section 4.2
xpts <- matrix(c(1,1,2,2,5),5,1)
ypts <- matrix(c(4,3,3,2,4),5,1)
Dmat = as.matrix(dist(cbind(xpts, ypts)))
Rmat = (1 - 3*Dmat/8 + Dmat^3/128)*(Dmat < 4)
Rmat

int1 <- matrix(1,5,1)
varofolse <- t(int1) %*% Rmat %*% int1/25
varofolse

1/(t(int1) %*% int1)

(sum(diag(Rmat)) - t(int1) %*% Rmat %*% int1/5)/4

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#              Three histograms for toy example 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Create plot of three histograms for the toy example of Section 4.2 (Figure 4.2)
set.seed(1005)
hold <- matrix(0,10000,3)
for(i in 1:10000){
  y <- t(chol(Rmat))%*%rnorm(5)
  hold[i,1] <- mean(y)
  hold[i,2] <- t(y-mean(y))%*%(y-mean(y))
  hold[i,3] <- mean(y)/sqrt((hold[i,2]/4)/5)
}

normmuhat <- function(x){
  y=(1/sqrt(2*pi*.2))*exp(-x^2/(2*.2))
  return(y)
}

chisighat2 <- function(x){
  y=(x*exp(-x/2))/4
  return(y)
}

tdist <- function(x){
  y=(3/(4*sqrt(5)))*(1+x^2/4)^(-5/2)
  return(y)
}

file_name = "figures/threehists"

pdf(paste0(file_name,'.pdf'), width = 9, height = 9)

  old.par = par(mfrow=c(3,1), mar = c(5,5,5,1), lwd = 2)
  hist(hold[,1], breaks=24, xlab="", xlim=c(-3,3), ylim=c(0,.9),
    freq=FALSE, cex.lab = 2, cex.axis = 1.5, cex.main = 2,
    main="Relative frequency histogram of OLS estimator of mu")
  curve(normmuhat, from=-3, to=3, add=TRUE)
  mtext('A', adj = -.08, padj = -.3, cex = 3)
  hist(hold[,2], breaks=24, xlab="", freq=FALSE, 
    cex.lab = 2, cex.axis = 1.5, cex.main = 2,
    main="Relative frequency histogram of residual sum of squares")
  curve(chisighat2,from=0,to=20,add=TRUE)
  mtext('B', adj = -.08, padj = -.3, cex = 3)
  hist(hold[,3], breaks=96, xlab="", xlim=c(-10,10), ylim=c(0,0.35),
    freq=FALSE, cex.lab = 2, cex.axis = 1.5, cex.main = 2,
    main="Relative frequency histogram of t test statistic")
  curve(tdist, from=-10, to=10, add=TRUE)
  mtext('C', adj = -.08, padj = -.3, cex = 3)
  par(old.par)
  
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

