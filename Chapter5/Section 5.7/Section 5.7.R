sec_path = 'Rcode/Chapter5/Section 5.7/'
setwd(paste0(SLEDbook_path,sec_path))

library(spmodel)
library(xtable)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#              Table 5.1 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Compute some results in Section 4.2
xpts <- matrix(c(1, 1, 2, 2, 5), 5, 1)
ypts <- matrix(c(4, 3, 3, 2, 4), 5, 1)
Dmat = as.matrix(dist(cbind(xpts, ypts)))
Rhat = (1 - 3*Dmat/6 + Dmat^3/54)*(Dmat < 3)
Rhat
Rlower = Rhat
Rlower[lower.tri(Rmat)] = NA

print(
    xtable(Rlower, 
      align = c('l',rep('l', times = length(Rlower[1,]))),
      digits = c(0, rep(3, times = 5))
    ),
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
