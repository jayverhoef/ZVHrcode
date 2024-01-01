sec_path = 'Rcode/Chapter4/Section 4.5/'
setwd(paste0(SLEDbook_path,sec_path))

library(ZVHdata)
library(ZVHdata)
library(xtable)
library(sf)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Seal Trends Coefficients Table
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# load data for graphics and analysis
data(sealPolys)
DF = st_drop_geometry(sealPolys)
str(DF)
DF = DF[!is.na(DF$Estimate),]
table(DF$stockid)
out = summary(lm(Estimate ~ -1 + I(as.factor(stockid)), data = DF))
out
out = out$coefficients[,1:2]

rownames(out) = c('$\\alpha_8$', '$\\alpha_9$', '$\\alpha_10$', 
  '$\\alpha_11$', '$\\alpha_12$')
# create table to paste into latex
  print(
    xtable(out, 
      align = c('l',rep('l', times = length(out[1,]))),
      digits = c(0, 4, 4),
    ),
    size = 'footnotesize',
    include.rownames = TRUE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
  )

