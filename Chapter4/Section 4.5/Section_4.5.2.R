sec_path = 'Rcode/Chapter4/Section 4.5/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(xtable)
library(spdep)

################################################################################
#-------------------------------------------------------------------------------
#              Caribou Forage Coefficients Table
#-------------------------------------------------------------------------------
################################################################################

# load data for graphics and analysis
data(caribouDF)
str(caribouDF)
caribouDF$i = as.factor(caribouDF$i)
caribouDF$j = as.factor(caribouDF$j)

#-------------------------------------------------------------------------------
# Table 4.2
#-------------------------------------------------------------------------------
# check out the R^2
summary(lm(z ~  water + tarp + water*tarp, data = caribouDF))
# create ANOVA table
out = summary(aov(z ~  water + tarp + water*tarp, data = caribouDF))
out
# turn into data.frame
outDF = as.data.frame(out[[1]])
outDF
# multiply Sum sq and Mean Sq by 1000 for readability of table
outDF[,2:3] = 1000*outDF[,2:3]
# add df and Sum sq as a row to the table
outDF = rbind(outDF, c(sum(outDF[,1]), sum(outDF[,2]), NA, NA, NA))
# add row names manually
outDF = cbind(c('Water', 'Tarp', 'Water $\\times$ Tarp', 
  'Error', 'Corrected Total'), outDF)
# add column names manually
colnames(outDF) = c('Source', 'df', 'Sum of squares', 'Mean square', 
  'F', '$P$-value')
# create table to paste into latex
  print(
    xtable(outDF, 
      align = c('l',rep('l', times = length(outDF[1,]))),
      digits = c(0, 0, 0, 2, 2, 2, 3),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.colnames.function = identity,
    only.contents = TRUE,
    include.colnames = TRUE,
    hline.after = c(-1, 0, nrow(outDF) - 1, nrow(outDF))
  )

#-------------------------------------------------------------------------------
# Moran's I and Geary's c on residuals
#-------------------------------------------------------------------------------

# get residuals
zr = residuals(lm(z ~  water + tarp + water*tarp, data = caribouDF))
# use distance to find rook's neighbors
TFmat = as.matrix(dist(caribouDF[,c('x','y')])) < 1.01 &
  as.matrix(dist(caribouDF[,c('x','y')])) > 0
# create an empty list
nei = vector('list',30)
# get neighbor indices
for(i in 1:30) nei[[i]] = as.integer(which(TFmat[i,]))
attr(nei,'class') = 'nb'
attr(nei,'type') = 'rook'
attr(nei,'row') = as.integer(caribouDF$i)
attr(nei,'col') = as.integer(caribouDF$j)
str(nei)
# Randomization test for Moran's I
moran.test(zr, 
  nb2listw(nei, style = 'W'),
  randomisation=TRUE)
# Randomization test for Geary's C  
geary.test(zr, 
  nb2listw(nei, style = 'W'),
  randomisation=TRUE)


#-------------------------------------------------------------------------------
# Table 4.3
#-------------------------------------------------------------------------------
# check out the R^2 and coefficient values for sum-to-zero constraints on rows
# and columns
summary(lm(z ~  i + j + water + tarp + water*tarp, data = caribouDF))
# make an ANOVA table where each factor is adjusted for all others first
# set row effect as last
out = summary(aov(terms(z ~  j + water + tarp + water*tarp + i, 
  keep.order = TRUE), data = caribouDF))
out
outDF = as.data.frame(out[[1]])[5,]
# set column effect as last
out = summary(aov(terms(z ~  i + water + tarp + water*tarp + 
  j, keep.order = TRUE), data = caribouDF))
out
outDF = rbind(outDF, as.data.frame(out[[1]])[5,])
# set water, tarp and tarp:water last, plus Error
out = summary(aov(terms(z ~  i + j + water + tarp + 
  water*tarp, keep.order = TRUE), data = caribouDF))
out
outDF = rbind(outDF, as.data.frame(out[[1]])[3:6,])
# add row of corrected Totals
outDF = rbind(outDF, c(sum(outDF[,1]), 
  var(caribouDF$z)*(length(caribouDF$z) - 1), NA, NA, NA))
# multiply Sum sq and Mean Sq by 1000 for readability of table
outDF[,2:3] = 1000*outDF[,2:3]
# add row names manually
outDF = cbind(c('Row', 'Column', 'Water', 'Tarp', 'Water $\\times$ Tarp', 
  'Error', 'Corrected Total'), outDF)
colnames(outDF) = c('Source', 'df', 'Sum of squares', 'Mean square', 
  'F', '$P$-value')

# create table to paste into latex
  print(
    xtable(outDF, 
      align = c('l',rep('l', times = length(outDF[1,]))),
      digits = c(0, 0, 0, 2, 2, 2, 3),
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.colnames.function = identity,
    only.contents = TRUE,
    include.colnames = TRUE,
    hline.after = c(-1, 0, nrow(outDF) - 1, nrow(outDF))
  )

# Drop the interaction term, and use sum-to-zero constraints on all factors
# other than tarp so intercept is at an "average" value, then get
# mean levels of tarp
contrasts(caribouDF$i) = contr.sum(6)
contrasts(caribouDF$j) = contr.sum(5)
contrasts(caribouDF$water) = contr.sum(2)
summary(lm(z ~  i + j + water + tarp, data = caribouDF))
# Notice that now column effects are significant
summary(aov(z ~  i + j + water + tarp, data = caribouDF))
# get tarp effects
# clear tarp
summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[1,1]
# no tarp
summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[1,1] +
  summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[12,1]
# shade tarp
summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[1,1] +
  summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[13,1]

#-------------------------------------------------------------------------------
# Moran's I and Geary's c on residuals
#-------------------------------------------------------------------------------

# get residuals
zr = residuals(lm(z ~  i + j + water + tarp + water*tarp, data = caribouDF))

# Randomization test for Moran's I
moran.test(zr, 
  nb2listw(nei, style = 'W'),
  randomisation=TRUE)
# Randomization test for Geary's C  
geary.test(zr, 
  nb2listw(nei, style = 'W'),
  randomisation=TRUE)

# get residuals without interaction term
zr = residuals(lm(z ~  i + j + water + tarp, data = caribouDF))

# Randomization test for Moran's I
moran.test(zr, 
  nb2listw(nei, style = 'W'),
  randomisation=TRUE)
# Randomization test for Geary's C  
geary.test(zr, 
  nb2listw(nei, style = 'W'),
  randomisation=TRUE)
