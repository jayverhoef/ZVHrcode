sec_path = 'Rcode/Chapter3/Section 3.5/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(xtable)
library(sf)
data(SO4obs)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Table so4-regionalsd
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

new_dt = data.frame(so4x = st_coordinates(SO4obs)[,1], 
	so4y = st_coordinates(SO4obs)[,2],
  so4 = SO4obs$SO4)
  
# Split the study region into 8 subregions. 
upper_half <- new_dt[new_dt$so4y >= median(new_dt$so4y), ]
lower_half <- new_dt[new_dt$so4y < median(new_dt$so4y), ]

sumup = summary(upper_half$so4x)
sumup['1st Qu.']

upper_sub1 <- upper_half[upper_half$so4x <= sumup['1st Qu.'], ]
upper_sub2 <- upper_half[upper_half$so4x > sumup['1st Qu.'] & 
  upper_half$so4x <= sumup['Median'], ]
upper_sub3 <- upper_half[upper_half$so4x > sumup['Median'] & 
  upper_half$so4x <= sumup['3rd Qu.'], ]
upper_sub4 <- upper_half[upper_half$so4x > sumup['3rd Qu.'], ]

sumlo = summary(lower_half$so4x)

lower_sub1 <- lower_half[lower_half$so4x <= sumlo['1st Qu.'], ]
lower_sub2 <- lower_half[lower_half$so4x > sumlo['1st Qu.'] & 
  lower_half$so4x <= sumlo['Median'], ]
lower_sub3 <- lower_half[lower_half$so4x > sumlo['Median'] & 
  lower_half$so4x <= sumlo['3rd Qu.'], ]
lower_sub4 <- lower_half[lower_half$so4x > sumlo['3rd Qu.'], ]

# create data for table with std. dev. and sample size

so4_regionalsd = rbind(
  unlist(lapply(list(upper_sub1, upper_sub2, upper_sub3, upper_sub4),
    function(dt){sd(dt$so4)})),
  unlist(lapply(list(lower_sub1, lower_sub2, lower_sub3, lower_sub4),
    function(dt){sd(dt$so4)}))
)
so4_regionalsd = data.frame(NS = c('North','South'), 
  FarW = so4_regionalsd[,1], NearW = so4_regionalsd[,2],
  NearE = so4_regionalsd[,3], FarE = so4_regionalsd[,4])


so4_regionaln = rbind(
  unlist(lapply(list(upper_sub1, upper_sub2, upper_sub3, upper_sub4),
    function(dt){nrow(dt)})),
  unlist(lapply(list(lower_sub1, lower_sub2, lower_sub3, lower_sub4),
    function(dt){nrow(dt)}))
)
so4_regionaln = data.frame(NS = c(NA,NA), 
  FarW = so4_regionaln[,1], NearW = so4_regionaln[,2],
  NearE = so4_regionaln[,3], FarE = so4_regionaln[,4])

print(
    xtable(so4_regionalsd, 
      digits = c(0, 0, 2, 2, 2, 2),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

print(
    xtable(so4_regionaln, 
      digits = c(0, 0, 0, 0, 0, 0),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

#-------------------copy results into Latex text
\begin{table}
\centering
\begin{tabular}{c|cccc}
& \multicolumn{4}{c}{Longitude class} \\
Latitude class & Far west & Near west & Near east & Far east \\ \hline
North & 1.69 & 2.71 & 8.67 & 9.21 \\
& 25 & 25 & 24 & 25 \\ 
\hline
South & 1.69 & 5.03 & 5.03 & 6.66 \\
& 25 & 24 & 24 & 25 \\ 
\hline
\end{tabular}
\caption{Standard deviations of the wet sulfate deposition data within eight latitude-longitude rectangular regions.  Sample sizes are given in parentheses.
\label{tab:so4-regionalsd}}
\end{table}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#                  Figure so4-sd-versus-mean 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Scatterplot of standard deviation vs. mean
file_name = 'figures/so4-sd-versus-mean'

mean_sd = data.frame(
  mean = unlist(lapply(list(upper_sub1, upper_sub2, upper_sub3, upper_sub4,
    lower_sub1, lower_sub2, lower_sub3, lower_sub4),
    function(dt){mean(dt$so4)})),
  sd = unlist(lapply(list(upper_sub1, upper_sub2, upper_sub3, upper_sub4,
    lower_sub1, lower_sub2, lower_sub3, lower_sub4),
    function(dt){sd(dt$so4)}))
)

pdf(paste0(file_name,'.pdf'))

  old.par = par(mar = c(5,5,1,1))
  plot(mean_sd[,1], mean_sd[,2], xlab = "SO4 Mean (kg/ha)", 
  ylab = "SO4 Standard Deviation (kg/ha)", 
  pch = 19, cex = 3, cex.lab = 2, cex.axis = 1.5)
  par(old.par)
  
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
#                  Barlett's Homogeneity Test  
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

sd_n = data.frame(
  sd = unlist(lapply(list(upper_sub1, upper_sub2, upper_sub3, upper_sub4,
    lower_sub1, lower_sub2, lower_sub3, lower_sub4),
    function(dt){sd(dt$so4)})),
  n = unlist(lapply(list(upper_sub1, upper_sub2, upper_sub3, upper_sub4,
    lower_sub1, lower_sub2, lower_sub3, lower_sub4),
    function(dt){nrow(dt)}))
)
(sum(sd_n$n) - nrow(sd_n))*log(mean(sd_n$sd^2)) -
  sum((sd_n$n - 1)*log(sd_n$sd^2))
qchisq(.99,7)

