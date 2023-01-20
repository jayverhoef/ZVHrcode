sec_path = 'Rcode/Chapter6/Section 6.1/figures/'
setwd(paste0(SLEDbook_path,sec_path))

### Double exponential MA function
double_exp <- function(x, sigma) exp(-abs(x)/sigma)
double_exp_plot <- function(x, sigma) exp(-abs(x-0.5)/sigma)

### Gaussian process generation
n1 <- 200
WN_x <- seq(0, 1, length.out = n1 + 1) 

file_name = "GP_generation"

pdf(paste0(file_name,'.pdf'), width = 11, height = 7)

	range <- c(0.008)
	for (k in 1:length(range)) {
		set.seed(200)
		## generate white noise (WN)
		WN_y <- rnorm(n1 + 1)
		
		n2 <- 1000
		## generate Gaussian process (GP)
		GP_x <- seq(0, 1, length.out = n2 + 1)
		GP_y <- rep(NA, length.out = n2 + 1)
		
		for (j in 1:(n2 + 1)) {
			add_term <- rep(NA, (n1 + 1))
			for (i in 1:(n1 + 1)) {
				X <- WN_y[i]
				num <- double_exp(GP_x[j] - WN_x[i], range[k])
				#dem <- 1
				#dem <- sum(double_exp(GP_x - WN_x[i], range[k]))
				dem <- sum(double_exp(GP_x[j] - WN_x, range[k]))
				add_term[i] <- X * num / dem
			}
			GP_y[j] <-sum(add_term)   
		}
		
		plot(WN_x, WN_y, type = 'h', xlab = '', ylab = '', ylim = c(-2.1, 4), main = '', axes = F)
		lines(GP_x, GP_y, type = 'l', lwd = 3)
	}

	## Moving-average function (double exponential)
	MA_x <- seq(0.4, 0.6, length.out = 1e5)
	MA_y <- double_exp_plot(MA_x, 0.008) + 3.2

	lines(MA_x, MA_y, type = 'l', lwd = 3, xlab = '', ylab = '')
	lines(MA_x, rep(3.17, 1e5), type = 'l')

	text(0.02,2.0, 'B', cex = 4)
	text(0.38,3.9, 'A', cex = 4)

dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))

# set.seed(200)
# WN_y <- rnorm(n1 + 1)
# 
# n2 <- 1000
# GP_x <- seq(0, 1, length.out = n2 + 1)
# GP_y <- rep(NA, length.out = n2 + 1)
# 
# for (j in 1:(n2 + 1)) {
#   add_term <- rep(NA, (n1 + 1))
#   for (i in 1:(n1 + 1)) {
#     X <- WN_y[i]
#     num <- double_exp(GP_x[j] - WN_x[i], 0.2)
#     #dem <- sum(double_exp(GP_x - WN_x[i], 0.2))
#     dem <- 1
#     add_term[i] <- X * num / dem
#   }
#   GP_y[j] <-sum(add_term)   
# }
# 
# plot(WN_x, WN_y, type = 'h', xlab = 'x', ylab = 'y', ylim = c(-5, 5))
# lines(GP_x, GP_y, type = 'l', lwd = 3)


