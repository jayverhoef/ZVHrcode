sec_path = 'Rcode/Chapter11/Section 11.3/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                         Get the Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# attach data library
library(ZVHdata)
library(sp)
library(viridis)
library(classInt)
library(colorspace)
library(spdep)
library(spmodel)
library(xtable)


# load data for graphics and analysis
data(caribouDF)

#-------------------------------------------------------------------------------
#                    Create Neighborhood Matrices
#-------------------------------------------------------------------------------

nTot = length(caribouDF)

# get Euclidean distance between centroids of plots
Distmat = as.matrix(dist(caribouDF[,c('x','y')]))
# create first-order neighbor matrix (rook's move) from distances
Nmat1 = (Distmat < 1.01)*1
diag(Nmat1) = 0
Nmat1rst = Nmat1/apply(Nmat1,1,sum)
# create second-order neighbor matrix (neighbors of neighbors) from first-order
Nmat2 = (Nmat1 %*% Nmat1 > 0 | Nmat1 > 0)*1
diag(Nmat2) = 0


################################################################################
#-------------------------------------------------------------------------------
#                     Section 11.3.1  Simulated Data
#-------------------------------------------------------------------------------
################################################################################

# get parameters from caribou data set using geostatistical exponential model
geo_exp = splm(z ~ water*tarp, data = caribouDF,  
	xcoord = 'x', ycoord = 'y', 
	spcov_type = 'exponential', estmethod = 'reml')
summary(geo_exp)

set.seed(102)
sim1D = sprnorm(spcov_params = coef(geo_exp, type = "spcov"),
	data = data.frame(x = 1:1000, y = rep(1, times = 1000)),
	xcoord = x, ycoord = y)

# randomized treatments
set.seed(2002)
trt_ran = sample(rep(1:6, times = 5))

################################################################################
# Figure 11.1
################################################################################

file_name = 'figures/Caribou_sim1000'
pdf(paste0(file_name,'.pdf'), width = 12, height = 9)

	layout(matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))
	
		padj = -.5
		adj = -.28
		cex_mtext = 3.3
		cex_lab = 2.4
		cex_axis = 1.5

		par(mar = c(7,7,4,2), mgp=c(4, 1.3, 0))
		plot(1:1000, sim1D, pch = 19, cex = .5, xlab = 'Coordinate',
			ylab = 'Simulated Value', cex.lab = cex_lab, cex.axis = cex_axis)
		points(321:350, sim1D[321:350], pch = 1, cex = 1.5)
		points(901:930, sim1D[901:930], pch = 1, cex = 1.5)
		mtext('A', adj = -.115, cex = cex_mtext, padj = padj)
		
		plot(321:350, sim1D[321:350], cex = 2, pch = 19, cex.lab = 2.5, 
			cex.axis = 1.5, xlab = 'Coordinate', ylab = 'Simulated Value')
		mtext('B', adj = adj, cex = cex_mtext, padj = padj)

		plot(901:930, sim1D[901:930], cex = 2, pch = 19, cex.lab = 2.5, 
			cex.axis = 1.5, xlab = 'Coordinate', ylab = 'Simulated Value', 
			ylim = c(.2, 1.65))
		text(901:930, rep(1.6, times = 30), 
			labels = as.character(trt_ran), cex = 1.7)
		text(901:930, rep(1.3, times = 30), 
			labels = as.character(kronecker(1:6,rep(1,times = 5))), cex = 1.7)
	  text(901:930, rep(1.45, times = 30), 
			labels = as.character(rep(1:6, times = 5)), cex = 1.7)
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

################################################################################
# Table 11.1
################################################################################

d3 = data.frame(x = 1:30, y = rep(1, times = 30), 
	trt = as.factor(trt_ran), z = sim1D[901:930])
# fit uncorrelated models to all 3 designs
sim_indran = splm(z ~ trt, data = d3, xcoord = x, ycoord = y,
	spcov_type = 'none')
summary(sim_indran)
sim_expran = splm(z ~ trt, data = d3, xcoord = x, ycoord = y,
	spcov_type = 'exponential')
summary(sim_expran)
sim_expr50 = splm(z ~ trt, data = d3, xcoord = x, ycoord = y,
	spcov_initial = spcov_initial(spcov_type = 'exponential', 
		range = 200, ie = 1e-5, de = 75,
		known = c('range')) )
summary(sim_expr50)
sim_cont = cbind(coef(sim_indran), coef(sim_expran), coef(sim_expr50),
	summary(sim_indran)$coefficients$fixed$Std_Error,
	summary(sim_expran)$coefficients$fixed$Std_Error,
	summary(sim_expr50)$coefficients$fixed$Std_Error)

print(
    xtable(sim_cont, 
      align = c('l',rep('l', times = length(sim_cont[1,]))),
      digits = c(0,rep(3, times = 3),rep(3, times = 3)),
      caption = 'Estimates and standard errors',
      label = 'tab:SealsFixEff'
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

mean(sim1D[901:930])
0.872 - 1.96*0.124
0.872 + 1.96*0.124
0.787 - 1.96*0.767
0.787 + 1.96*0.767

