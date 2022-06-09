sec_path = 'Rcode/Chapter4/Section 4.5/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(xtable)
library(spdep)
library(car)

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
    hline.after = c(-1, 0, nrow(outDF) - 1)
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
# use sum-to-zero constraints on all factors
options(contrasts = c("contr.sum", "contr.poly"))
summary(lm(z ~  i + j + tarp + water + water*tarp, data = caribouDF))
# make an SAS type III ANOVA table 
summary(lm(z ~  i + j + water + tarp + water*tarp, data = caribouDF))
model = lm(z ~  i + j + water + tarp + water*tarp, data = caribouDF)
out = Anova(model, type = 'III')
out = out[2:7,]
out[,1] = 1000*out[,1]
out = cbind(out, out[,1]/out[,2])
out = rbind(out,
  c(1000*var(caribouDF$z)*sum(out[,2]), sum(out[,2]), NA, NA, NA))
outDF = data.frame(
  Source = c('Row','Column','Water','Tarp','Water $\\times$ Tarp','Error','Corrected Total'),
  df = out[,2], SS = out[,1], MS = out[,5], Fval = out[,3], 
  Pval = out[,4])
colnames(outDF) = c('Source','df','Sum of squares','Mean square','F','$P$-value')
outDF

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
    hline.after = c(-1, 0, nrow(outDF) - 1)
  )

# get tarp effects
# clear tarp
summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[1,1] +
  summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[12,1]
# no tarp
summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[1,1] +
  summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[13,1]
# shade tarp
summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[1,1] -
  summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[12,1] -
  summary(lm(z ~  i + j + water + tarp + water:tarp, data = caribouDF))$coef[13,1]

#drop nonsignificant interaction term
summary(lm(z ~  i + j + water + tarp, data = caribouDF))
model = lm(z ~  i + j + water + tarp, data = caribouDF)
out = Anova(model, type = 'III')
out = out[2:6,]
out[,1] = 1000*out[,1]
out = cbind(out, out[,1]/out[,2])
out = rbind(out,
  c(1000*var(caribouDF$z)*sum(out[,2]), sum(out[,2]), NA, NA, NA))
outDF = data.frame(
  Source = c('Row','Column','Water','Tarp','Error','Corrected Total'),
  df = out[,2], SS = out[,1], MS = out[,5], Fval = out[,3], 
  Pval = out[,4])
colnames(outDF) = c('Source','df','Sum of squares','Mean square','F','$P$-value')
outDF

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
    hline.after = c(-1, 0, nrow(outDF) - 1)
  )

# get tarp effects
# clear tarp
summary(lm(z ~  i + j + water + tarp, data = caribouDF))$coef[1,1] +
  summary(lm(z ~  i + j + water + tarp, data = caribouDF))$coef[12,1]
# no tarp
summary(lm(z ~  i + j + water + tarp, data = caribouDF))$coef[1,1] +
  summary(lm(z ~  i + j + water + tarp, data = caribouDF))$coef[13,1]
# shade tarp
summary(lm(z ~  i + j + water + tarp, data = caribouDF))$coef[1,1] -
  summary(lm(z ~  i + j + water + tarp, data = caribouDF))$coef[12,1] -
  summary(lm(z ~  i + j + water + tarp, data = caribouDF))$coef[13,1]

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
