sec_path = 'Rcode/Chapter10/Section 11.5/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#        CPE-optimal Designs for 25 x 25 Grid with Mat(0.5) model
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

gdim = 4
nsamp = 4
# grid of spatial locations
xy = data.frame(x = kronecker(1:gdim, rep(1, times = gdim)),
	y = kronecker(rep(1, times = gdim), 1:gdim))
# distance among all locations
distmat = as.matrix(dist(xy))

# compute distance with an offset of (a,a) for second set
a = 1
distmatoff = sqrt((outer(xy$x, rep(1, times = gdim^2)) - 
	outer(rep(1, times = gdim^2), xy$x + a))^2 + 
	(outer(xy$y, rep(1, times = gdim^2)) - 
	outer(rep(1, times = gdim^2), xy$y + a))^2)


# autocorrelation parameter
rho = 0.8
# crosscorrelation parameter
rho_c = 0.8
# base correlation matrix for each variable
baseC = rho^distmat
# base crosscorrelation with offset
baseCoff = rho^distmatoff
# partial of correlation matrix with respect to rho
baseC_drho = distmat*rho^(distmat - 1)
# partial of correlation matrix with respect to rho with offset
baseCoff_drho = distmatoff*rho^(distmatoff - 1)
# partial of crosscorrelation matrix with respect to rho_c
baseC_drho_c = baseC
# partial of crosscorrelation matrix with respect to rho_c when offset used
baseCoff_drho_c = baseCoff
# matrix of all zeros
zeros = matrix(0, nrow = gdim^2, ncol = gdim^2)
# joint correlation and crosscorrelation matrix for both variables
bigC = rbind(cbind(baseC, rho_c*baseC),
			cbind(rho_c*baseC, baseC))
# joint correlation and crosscorrelation matrix for both variables with offset
bigCoff = rbind(cbind(baseC, rho_c*baseCoff),
			cbind(rho_c*t(baseCoff), baseC))
# partial of joint correlation and crosscorrelation matrix for both
# variables with respect to rho
bigC_drho = rbind(cbind(baseC_drho, rho_c*baseC_drho),
									cbind(rho_c*baseC_drho, baseC_drho))
# partial of joint correlation and crosscorrelation matrix for both
# variables with respect to rho with offset
bigCoff_drho = rbind(cbind(baseC_drho, rho_c*baseCoff_drho),
									cbind(rho_c*t(baseCoff_drho), baseC_drho))
# partial of joint correlation and crosscorrelation matrix for both
# variables with respect to rho_c
bigC_drho_c = rbind(cbind(zeros, baseC),
									cbind(baseC, zeros))
# partial of joint correlation and crosscorrelation matrix for both
# variables with respect to rho_c when using offset
bigCoff_drho_c = rbind(cbind(zeros, baseCoff),
									cbind(baseCoff, zeros))

# function to compute Fisher Information matrix
# samp: a list where first item is named v1 containing
#		indexes of sample locations for variable 1
#   and second item is named v2 containing
#		indexes of sample locations for variable 2
# bigC: autocorrelation and crosscorrelaton matrix for both variables
# bigC_drho: partial of autocorrelation/crosscorrelation matrix for rho
# bigC_drho_c: partial of autocorrelation/crosscorrelation matrix for rho_c
detFI = function(samp, bigC, bigC_drho, bigC_drho_c)
{
	# subset joint correlation matrix for sample locations
	sampC = bigC[c(samp$v1, samp$v2), c(samp$v1, samp$v2)]
	# inverse
	sampCi = solve(sampC)
	# subset joint partial matrices
	sampC_drho = bigC_drho[c(samp$v1, samp$v2), c(samp$v1, samp$v2)]
	sampC_drho_c = bigC_drho_c[c(samp$v1, samp$v2), c(samp$v1, samp$v2)]
	# compute sampCi %*% 1, where 1 is a column of 1's
	sampCi_rowsum = apply(sampCi, 1, sum)
	# projection matrix: Ci - Ci %*% X (t(X) %*% Ci %*% X)^{-1} %*% t(X) %*% Ci
	# for the simplified case where X is a column of 1's
	# e.g., (t(X) %*% Ci %*% X)^{-1} is simply 1/sum(Ci) and
	# Ci %*% X is a vector of row sums.
	Pmat = sampCi - outer(sampCi_rowsum, sampCi_rowsum)/sum(sampCi)

	# Fisher Information
	FI = matrix(NA, nrow = 2, ncol = 2)
	#1/2 trace of each element
	FI[1,1] = sum(diag(sampCi %*% sampC_drho %*% 
		sampCi %*% sampC_drho))/2
	FI[1,2] = FI[2,1] = sum(diag(sampCi %*% sampC_drho %*% 
		sampCi %*% sampC_drho_c))/2
	FI[2,2] = sum(diag(sampCi %*% sampC_drho_c %*% 
		sampCi %*% sampC_drho_c))/2
	# determinant of Fisher Information
	FI[1,1]*FI[2,2] - FI[1,2]^2
}

# initial samples, sample 1 is from integers from 1 to gdim^2
# second sample is from integers gdim^2 + 1 to 2*gdim^2
sample1 = sample(1:gdim^2, nsamp)
sample2 = sample(gdim^2 + 1:gdim^2, nsamp)
samp = list(v1 = sample1, v2 = sample2)
# initial samples, sample 1 is from integers from 1 to gdim^2
# second sample is from integers gdim^2 + 1 to 2*gdim^2
#sample1 = ((0:(gdim - 1))*gdim + 1)[1:nsamp] + 0:(nsamp - 1)
#sample2 = gdim^2 + ((0:(gdim - 1))*gdim + 1)[1:nsamp] + 0:(nsamp - 1)

# keep track of overall minimum during simulated annealing
minall = 1e+32
# set a counter to zero
j = 0
# annealing multiplier, can change if run below in batches. As anmult gets
# bigger, the algorithm gets greedier, and more prone to get stuck locally.
anmult = .01
big
# can run this in batches to keep track of progress
for(i in 1:100000) {
	# increment the counter
	j = j + 1
	# chose a new location for variable 1 from those that are not currently
	# chosen
	new1 = sample(which(!(1:gdim^2 %in% sample1)),1)
	# copy the current sample
	sample1_try = sample1
	# create the new sample with a single new location
	sample1_try[sample(1:nsamp,1)] = new1
	# keep them colocated
	sample2_col = sample1_try + gdim^2
	# chose a new location for variable 2 from those that are not currently
	# chosen
	new2 = sample((gdim^2 + 1:gdim^2)
		[which(!((gdim^2 + 1:gdim^2) %in% sample2))],1)
	# copy the current sample
	sample2_try = sample2
	# create the new sample with a single new location
	sample2_try[sample(1:nsamp,1)] = new2
	# keep them colocated
	sample1_col = sample2_try - gdim^2
	# 1/|Fisher-Information| for current sample
	old = 1/detFI(sample1, sample2, bigC, bigC_drho, 
		bigC_drho_c)
	# 1/|Fisher-Information| for trying new location for variable 1
	try1 = 1/detFI(sample1_try, sample2, bigC, bigC_drho, 
		bigC_drho_c)
	# 1/|Fisher-Information| for trying new location for variable 1 and colocated
	try1_col = 1/detFI(sample1_try, sample2_col, bigC, bigC_drho,
		bigC_drho_c)
	# 1/|Fisher-Information| for trying new location for variable 2
	try2 = 1/detFI(sample1, sample2_try, bigC, bigC_drho, 
		bigC_drho_c)
	# 1/|Fisher-Information| for trying new location for variable 2 and colocated
	try2_col = 1/detFI(sample1_col, sample2_try, bigC, bigC_drho,
		bigC_drho_c)
	# 1/|Fisher-Information| for trying new location for both variables
	try12 = 1/detFI(sample1_try, sample2_try, bigC, bigC_drho,
		bigC_drho_c)
	# which one has the minimum |Fisher-Information|?
	minDet = min(old, try1, try1_col, try2, try2_col, try12)
	# move to a new samples if it has minimum 1/|Fisher-Information|
	if(minDet == try1) 
		sample1 = sample1_try
	if(minDet == try1_col) {
		sample1 = sample1_try
		sample2 = sample2_col
	}
	if(minDet == try2) 
		sample2 = sample2_try
	if(minDet == try2_col) {
		sample1 = sample1_col
		sample2 = sample2_try
	}
	if(minDet == try12) { 
		sample1 = sample1_try
		sample2 = sample2_try
	}
	# we may move off the best sample due to the randomness of simulated
	# annealing, so keep track of best one ever tried
	if(minDet < minall) {
		minall = minDet
		bestsamp1 = sample1
		bestsamp2 = sample2
	}
	# randomly choose a sample with decaying annealing schedule where
	# the probability depends on the counter: 
	# (j*.01*log(j))/((j*.01*log(j))+1) < runif(1)
	# So it is chosen randomly less and less as j increases
	if((j*anmult*log(j))/((j*anmult*log(j))+1) < runif(1)) {
		whichone = sample(1:5,1)
		if(whichone == 1) sample1 = sample1_try
		if(whichone == 2) {
			sample1 = sample1_try
			sample2 = sample2_col
		}
		if(whichone == 3) sample2 = sample2_try
		if(whichone == 4) {
			sample1 = sample1_col
			sample2 = sample2_try
		}
		if(whichone == 5) {
			sample1 = sample1_try
			sample2 = sample2_try
		}
	}
}
(j*anmult*log(j))/((j*anmult*log(j))+1)
minDet
minall
# plot the currently chosen design
plot(xy[sample1,], cex = 2, xlim = c(1,gdim), ylim = c(1,gdim))
points(xy[sample2 - gdim^2,], pch = 19)

# plot the best one ever tested
plot(xy[bestsamp1,], cex = 2, xlim = c(1,gdim), ylim = c(1,gdim))
points(xy[bestsamp2 - gdim^2,], pch = 19)


# Figure 4 from Li and Zimmerman
LiZ_x = c(7,7,7,8,8,8,10,11,11,12,12,17,17,18,18)
LiZ_y = c(18,19,20,18,19,20,20,19,20,19,20,5,6,5,6)
plot(LiZ_x,LiZ_y, pch = 19, xlim = c(1,25), ylim = c(1,25))

LiZsamples = NULL
for(i in 1:15) LiZsamples = 
	c(LiZsamples, which(LiZ_x[i] == xy$x & LiZ_y[i] == xy$y))

1/detFI(LiZsamples, 25^2 + LiZsamples, bigC, bigC_drho, bigC_drho_c)

# LiZ simulated annealing cooling schedule
# colocted.  This is a completely new index drawn so no repeats of any current
# indexes.
LiZ_ProbNew = function(samp, samp_new, tau_0, Sig, Sig_drho, 
		Sig_drho_c)
{
	gam_samp_old = 1/detFI(samp, Sig, Sig_drho,
		Sig_drho_c)
	gam_samp_new = 1/detFI(samp_new, Sig, Sig_drho,
		Sig_drho_c)
	if(gam_samp_new <= gam_samp_old) prob_new = 1
	if(gam_samp_new > gam_samp_old) prob_new = 
		exp((gam_samp_old - gam_samp_new)/tau_0)
	prob_new
}

# My way to sample a new location for variable 1 or two, or keep them
# colocted.  This is a completely new index drawn so no repeats of any current
# indexes.
sample_jay = function(samp, whichvar, colocate = FALSE) {
	if(whichvar == 1) {
		# sample an index from those not already in sample
		new1 = sample(which(!(1:gdim^2 %in% samp$v1)),1)
		# replace one randomly chosen index in the sample with the new index
		samp$v1[sample(1:nsamp,1)] = new1
		if(colocate == TRUE) samp$v2 = gdim^2 + samp$v1
	}
	if(whichvar == 2) {
		# sample an index from those not already in sample
		new1 = gdim^2 + sample(which(!((gdim^2 + 1:gdim^2) %in% samp$v2)),1)
		# replace one randomly chosen index in the sample with the new index
		samp$v2[sample(1:nsamp,1)] = new1
		if(colocate == TRUE) samp$v1 = samp$v2 - gdim^2
	}
	samp
}


	1/detFI(samp, bigC, bigC_drho, bigC_drho_c)
	1/detFI(samp_new, bigC, bigC_drho, bigC_drho_c)

	samp_new = sample_jay(samp, 1, colocate = FALSE)

	LiZ_ProbNew(samp = samp, samp_new = samp_new, tau_0 = 0.0001, Sig = bigC, 
		Sig_drho = bigC_drho, Sig_drho_c = bigC_drho_c)
		

sample_LiZ = function(samp, whichvar, hmax, xy) {
	if(whichvar == 1) {
		# chose one of the indexes to change
		whichindex = sample(1:nsamp,1)
		# test condition for while loop
		testcond = TRUE
		while(testcond == TRUE) {
			# convert from sample indexes to coordinates
			xy_whichindex = xy[samp$v1[whichindex],]
			# choose one of the coordinates at random
			whichcoord = sample(1:2,1)
			# move the coordinate value according to hmax
			newcoord_try = xy_whichindex[,whichcoord] + 
				sample(c(1:hmax,-(1:hmax)),1)
			while(newcoord_try < 1 | newcoord_try > gdim)
				newcoord_try = xy_whichindex[,whichcoord] + 
				sample(c(1:hmax,-(1:hmax)),1)
			xy_new = xy_whichindex
			xy_new[,whichcoord] = newcoord_try
			new_index = which(xy[,1] == xy_new[,1] & xy[,2] == xy_new[,2])
			testcond = any(new_index == samp$v1)
		}
		samp$v1[whichindex] = new_index
	}
	if(whichvar == 2) {
		# chose one of the indexes to change
		whichindex = sample(1:nsamp,1)
		# test condition for while loop
		testcond = TRUE
		while(testcond == TRUE) {
			# convert from sample indexes to coordinates
			xy_whichindex = xy[samp$v2[whichindex] - gdim^2,]
			# choose one of the coordinates at random
			whichcoord = sample(1:2,1)
			# move the coordinate value according to hmax
			newcoord_try = xy_whichindex[,whichcoord] + 
				sample(c(1:hmax,-(1:hmax)),1)
			while(newcoord_try < 1 | newcoord_try > gdim)
				newcoord_try = xy_whichindex[,whichcoord] + 
				sample(c(1:hmax,-(1:hmax)),1)
			xy_new = xy_whichindex
			xy_new[,whichcoord] = newcoord_try
			new_index = gdim^2 + which(xy[,1] == xy_new[,1] & xy[,2] == xy_new[,2])
			testcond = any(new_index == samp$v2)
		}
		samp$v2[whichindex] = new_index
	}
	samp
}

sample1 = sample(1:gdim^2, nsamp)
sample2 = sample(gdim^2 + 1:gdim^2, nsamp)
samp = list(v1 = sample1, v2 = sample2)
j = 0

for(i in 1:1000) {
	j = j + 1
	samp_new =  sample_LiZ(samp, 1, 1, xy)
	U = LiZ_ProbNew(samp = samp, samp_new = samp_new, tau_0 = 1/j, 
		Sig = bigC, Sig_drho = bigC_drho, Sig_drho_c = bigC_drho_c)
	if(runif(1) < U) samp = samp_new	
	samp_new =  sample_LiZ(samp, 2, 1, xy)
	U = LiZ_ProbNew(samp = samp, samp_new = samp_new, tau_0 = 1/j, 
		Sig = bigC, Sig_drho = bigC_drho, Sig_drho_c = bigC_drho_c)
	if(runif(1) < U) samp = samp_new	
}
# plot the currently chosen design
plot(xy[samp$v1,], cex = 2, xlim = c(1,gdim), ylim = c(1,gdim))
points(xy[samp$v2 - gdim^2,], pch = 19)


sample1 = sample(1:gdim^2, nsamp)
sample2 = sample(gdim^2 + 1:gdim^2, nsamp)
samp = list(v1 = sample1, v2 = sample2)
j = 0
minall = 1e+32

for(i in 1:100000) {
	j = j + 1
	samp_new =  sample_jay(samp, 1, TRUE)
	U = LiZ_ProbNew(samp = samp, samp_new = samp_new, tau_0 = 10/j, 
		Sig = bigC, Sig_drho = bigC_drho, Sig_drho_c = bigC_drho_c)
	if(runif(1) < U) samp = samp_new
	minDet = 1/detFI(samp, bigC, bigC_drho,
		bigC_drho_c)
	if(minDet < minall) {
		minall = minDet
		bestsamp = samp
	}
	samp_new =  sample_jay(samp, 2, TRUE)
	U = LiZ_ProbNew(samp = samp, samp_new = samp_new, tau_0 = 10/j, 
		Sig = bigC, Sig_drho = bigC_drho, Sig_drho_c = bigC_drho_c)
	minDet = 1/detFI(samp, bigC, bigC_drho,
		bigC_drho_c)
	if(minDet < minall) {
		minall = minDet
		bestsamp = samp
	}
	if(runif(1) < U) samp = samp_new	
	samp_new =  sample_jay(samp, 1, FALSE)
	U = LiZ_ProbNew(samp = samp, samp_new = samp_new, tau_0 = 10/j, 
		Sig = bigC, Sig_drho = bigC_drho, Sig_drho_c = bigC_drho_c)
	minDet = 1/detFI(samp, bigC, bigC_drho,
		bigC_drho_c)
	if(minDet < minall) {
		minall = minDet
		bestsamp = samp
	}
	if(runif(1) < U) samp = samp_new	
	samp_new =  sample_jay(samp, 2, FALSE)
	U = LiZ_ProbNew(samp = samp, samp_new = samp_new, tau_0 = 10/j, 
		Sig = bigC, Sig_drho = bigC_drho, Sig_drho_c = bigC_drho_c)
	minDet = 1/detFI(samp, bigC, bigC_drho,
		bigC_drho_c)
	if(minDet < minall) {
		minall = minDet
		bestsamp = samp
	}
	if(runif(1) < U) samp = samp_new	
}
minall*1e+5
1/detFI(samp, bigC, bigC_drho,
		bigC_drho_c)*1e+5
# plot the currently chosen design
plot(xy[samp$v1,], cex = 2, xlim = c(1,gdim), ylim = c(1,gdim))
points(xy[samp$v2 - gdim^2,], pch = 19)


plot(xy[bestsamp$v1,], cex = 2, xlim = c(1,gdim), ylim = c(1,gdim))
points(xy[bestsamp$v2 - gdim^2,], pch = 19)
