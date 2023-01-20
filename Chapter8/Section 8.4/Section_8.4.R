#set a path as a working directory
sec_path = 'Rcode/Chapter8/Section 8.4/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#             Code for Table on Information Content
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Code for Table 8.4 and the accompanying figure

# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 12
n <- K^2
locxy <- cbind(rep(1:K,K), rep(1:K, each = K))

# Set parameter values for Case I; modify straightforwardly for Cases II and III
psill <- 0.75
nugget <- 0.25
range <- 2.0

# create a function to extract information matrix parts
sp_info_exp = function(locxy, psill, nugget, range)
{
	# Obtain the covariance matrix for an exponential covariance function with partial sill, nugget, and range
	d = as.matrix(dist(locxy))
	Sigma = psill*exp(-d/range)+nugget*diag(144)

	# Obtain first-order derivatives of Sigma with respect to each parameter
	deriv1 <- exp(-d/range)
	deriv2 <- diag(144)
	deriv3 <- psill*(1/range^2)*d*deriv1

	# Obtain elements of information matrix
	Sigmainv <- solve(Sigma)
	info <- diag(3)
	info[1,1] <- 0.5*sum(diag(Sigmainv%*%deriv1%*%Sigmainv%*%deriv1))
	info[1,2] <- 0.5*sum(diag(Sigmainv%*%deriv1%*%Sigmainv%*%deriv2))
	info[1,3] <- 0.5*sum(diag(Sigmainv%*%deriv1%*%Sigmainv%*%deriv3))
	info[2,1] <- info[1,2]
	info[2,2] <- 0.5*sum(diag(Sigmainv%*%deriv2%*%Sigmainv%*%deriv2))
	info[2,3] <- 0.5*sum(diag(Sigmainv%*%deriv2%*%Sigmainv%*%deriv3))
	info[3,1] <- info[1,3]
	info[3,2] <- info[2,3]
	info[3,3] <- 0.5*sum(diag(Sigmainv%*%deriv3%*%Sigmainv%*%deriv3))

	# Invert information matrix and reciprocate diagonal elements 
	# to get information content
	infoinv <- solve(info)
	infocontent <- 1/diag(infoinv)

	# Form asymptotic correlation matrix from information matrix
	infoinvsds <- diag(sqrt(diag(infoinv)))
	asympcorrs <- solve(infoinvsds)%*%infoinv%*%solve(infoinvsds)
	list(infocontent = infocontent, asympcorrs = asympcorrs)
}

# get information matrix parts for various cases as given in text
case1 = sp_info_exp(locxy = locxy, psill = 0.75, nugget = 0.25, range = 2.0)
case2 = sp_info_exp(locxy = locxy, psill = 0.25, nugget = 0.75, range = 2.0)
case3 = sp_info_exp(locxy = locxy, psill = 0.75, nugget = 0.25, range = 3.0)
cluster <- c(0.25, 0.75, 2.25, 2.75, 4.25, 4.75, 6.25, 
	6.75, 8.25, 8.75, 10.25, 10.75)
loccl <- cbind(rep(cluster, K), rep(cluster, each = K))
case4 = sp_info_exp(locxy = loccl, psill = 0.75, nugget = 0.25, range = 2.0)

# create a matrix of results for table
Info_table = rbind(
	c(case1$infocontent[1:2], case1$asympcorrs[1,2:3], case1$asympcorrs[2,3]),
	c(case2$infocontent[1:2], case2$asympcorrs[1,2:3], case2$asympcorrs[2,3]),
	c(case3$infocontent[1:2], case3$asympcorrs[1,2:3], case3$asympcorrs[2,3]),
	c(case4$infocontent[1:2], case4$asympcorrs[1,2:3], case4$asympcorrs[2,3])
)

# print table in latex format
print(
    xtable(Info_table, 
      align = c('l',rep('l', times = length(Info_table[1,]))),
      digits = c(0,1,1,2,2,2),
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#             Figure with clustered locations 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Replace the assignment of locxy with the next two lines for Case IV
cluster <- c(0.25, 0.75, 2.25, 2.75, 4.25, 4.75, 6.25, 
	6.75, 8.25, 8.75, 10.25, 10.75)
loccl <- cbind(rep(cluster, K), rep(cluster, each = K))

file_name = "figures/spatial-config-info"
pdf(file = paste0(file_name,'.pdf'))
	par(pty="s")
	plot(loccl[,1], loccl[,2], xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
		pch = 16, cex = 1.5)
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
