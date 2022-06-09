#set a path as a working directory
sec_path = 'Rcode/Chapter8/Section 8.2/'
setwd(paste0(SLEDbook_path,sec_path))

library(xtable)

################################################################################
#-------------------------------------------------------------------------------
#              ML versus REML for Chapter 8
#-------------------------------------------------------------------------------
################################################################################

# define spherical covariance matrix as function of unscaled distance matrix
spherical <- function(distance.matrix)
{
	CovMat = (1 - 1.5*distance.matrix + 0.5*distance.matrix^3)
	CovMat[distance.matrix > 1] = 0
	CovMat
}

# define exponential covariance matrix as function of unscaled distance matrix
exponential <- function(distance.matrix)
{
	exp(-3*distance.matrix) 
}

# define Gaussian covariance matrix as function of unscaled distance matrix
Gaussian <- function(distance.matrix)
{
	exp(-sqrt(3)*distance.matrix^2) + 
		diag(rep(1e-11, times = dim(distance.matrix)[1])) 
}

# function to set up square grid with unit spacing
create_grid = function(K)
{
  locx <- matrix(rep(1:K,K),ncol=1,byrow=T)
  locy <- matrix(rep(1:K,each=K),ncol=1,byrow=T)
  cbind(locx, locy)
}

# function to simulate spatially uniform random locations
create_spatialUnif = function(n, coordmax = 1)
{
  locx <- coordmax*runif(n)
  locy <- coordmax*runif(n)
  cbind(locx, locy)
}

#useful functions
expit = function(x){exp(x)/(1 + exp(x))}
logit = function(x){log(x/(1 - x))}

# minus 2 times the log-likelihood (which has been profiled) to be optimized
LL_mle = function(theta, covf, distmat, XX, y1, maxrange = 1, estmeth = 'mle')
{
			rangetry = maxrange*expit(theta[1])
			ratiotry = expit(theta[2])
			n = length(y1)
			p = dim(XX)[2]
			G = ((1-ratiotry)*covf(distmat/rangetry) + ratiotry*diag(n))
			Ginv = solve(G)
			covB = solve(t(XX)%*%Ginv%*%XX)
			covBXGi = covB%*%t(XX)%*%Ginv
			Q = Ginv-Ginv%*%XX%*%covBXGi
			if(estmeth == 'mle') LL = log(det(G)) + n*log(t(y1)%*%Q%*%y1)
			if(estmeth =='remle') LL = log(det(G)) + log(det(t(XX)%*%Ginv%*%XX)) + 
				(n-p)*log(t(y1)%*%Q%*%y1)
			LL
}   

#set simulation seed for repeatability
set.seed(1111)

#set covariance parameters for the simulation
covf = exponential
sigmap = 1
#set fixed effect parameters
beta0 = 3
beta1 = 1
beta2 = 2
# set simulation parameters
no = 100
nsim = 1000


store_mle <- matrix(0,nsim,15)
store_remle <- matrix(0,nsim,15)
maxrange = 4*sqrt(2)
scenario = list(
	ratiop = c(.99, 0.5, 0.5, 0.01, 0.01),
	rangep = c(.01, 0.5, 1.5, 0.5, 1.5)
)	
hold_results = matrix(NA, nrow = 27, ncol = 10)
start.time = Sys.time()
for(m in 1:5) {

	ratiop = scenario$ratiop[m]
	rangep = scenario$rangep[m]

	#start a simulation
	k = 1
	for(k in 1:nsim){
		cat("\r", "Scenario: ", m, "   Simulation: ", k)
		# create spatial locations
		xy = create_spatialUnif(no)
		xy = rbind(xy,c(0.5,0.5),c(1,1))
		n <- dim(xy)[1]
		nu = n - no
		# simulate data
		distmat = as.matrix(dist(xy))
		Gtrue = sigmap*((1-ratiop)*covf(distmat/rangep) + ratiop*diag(n))
		beta = c(beta0, beta1, beta2)
		X <- cbind(rep(1, times = n), rnorm(n), xy[,1])
		p <- 3
		y1 = X %*% beta + t(chol(Gtrue))%*%rnorm(n)
		# partition into observed and unobserved parts
		yo = y1[1:no]
		yu = y1[(no+1):n]
		Xo = X[1:no,]
		Xu = X[(no+1):n,]
		distmato = distmat[1:no,1:no]
		distmatou = distmat[1:no,(no+1):n]
		distmatu = distmat[(no+1):n,(no+1):n]
		rec = 0
		gridvals = c((1:5)/1000, (1:5)/100, (1:9)/10, (95:99)/100, (995:999)/1000)
		hold = matrix(0,length(gridvals)^2,4)
		# grid search for both mle and remle
		for(i in gridvals){
			for(j in gridvals) {
				rec = rec + 1
				rangetry = i*maxrange
				ratiotry = j
				G = ((1-ratiotry)*covf(distmato/rangetry) + ratiotry*diag(no))
				Ginv = solve(G)
				covB = solve(t(Xo)%*%Ginv%*%Xo)
				covBXGi = covB%*%t(Xo)%*%Ginv
				Q = Ginv-Ginv%*%Xo%*%covBXGi
				L = -0.5*log(det(G))-(n/2)*log(t(yo)%*%Q%*%yo)
				LR = -0.5*log(det(G))- 0.5*log(det(t(Xo)%*%Ginv%*%Xo))-
					((n-p)/2)*log(t(yo)%*%Q%*%yo)
				hold[rec,] = c(rangetry, ratiotry, L, LR)
			}
		}
		# find optimum values from grid search
		mleRec = max(which(hold[,3] == max(hold[,3])))
		theta_mle = hold[mleRec,1:2]
		remleRec = max(which(hold[,4] == max(hold[,4])))
		theta_remle = hold[remleRec,1:2]
		# use optim to optimize around the best grid search values
		# transform values so that optimation is unconstrained
		theta_mle[1] = logit(theta_mle[1]/maxrange)
		theta_mle[2] = logit(theta_mle[2])
		theta_mle = optim(theta_mle, LL_mle, covf = covf, distmat = distmato, XX = Xo, 
			y1 = yo, maxrange = 2*sqrt(2), estmeth = 'mle')$par
		theta_mle[1] = expit(theta_mle[1])*maxrange
		theta_mle[2] = expit(theta_mle[2])
		theta_remle[1] = logit(theta_remle[1]/maxrange)
		theta_remle[2] = logit(theta_remle[2])
		theta_remle = optim(theta_remle, LL_mle, covf = covf, distmat = distmato, XX = Xo, 
			y1 = yo, maxrange = 2*sqrt(2), estmeth = 'remle')$par
		theta_remle[1] = expit(theta_remle[1])*maxrange
		theta_remle[2] = expit(theta_remle[2])
		# store mle quantities
		store_mle[k,1:2] = theta_mle
		G = (1-theta_mle[2])*covf(distmato/theta_mle[1]) + theta_mle[2]*diag(no)
		Ginv = solve(G)
		covB = solve(t(Xo)%*%Ginv%*%Xo)
		covBXGi = covB%*%t(Xo)%*%Ginv
		Q = Ginv-Ginv%*%Xo%*%covBXGi
		sig2 = as.numeric(t(yo)%*%Q%*%yo/no)
		Ginv = Ginv/sig2
		covB = sig2*covB
		store_mle[k,3] = sig2
		store_mle[k,4:6] = covBXGi%*%yo
		store_mle[k,7:9] = diag(covB)
		cc = sig2*(1-theta_mle[2])*covf(distmatou/theta_mle[1])
		preds = Xu %*% covBXGi %*% yo + t(cc) %*% 
			Ginv %*% (yo - Xo %*% covBXGi%*%yo)
		store_mle[k,10:11] = yu
		store_mle[k,12:13] = preds
		Gu = sig2*((1-theta_mle[2])*covf(distmatu/theta_mle[1]) + 
			theta_mle[2]*diag(nu))
		dd = Xu - t(cc) %*% Ginv %*% Xo
		store_mle[k,14:15] = diag(Gu - t(cc) %*% Ginv %*% cc + 
			dd %*% covB %*% t(dd))
		# store remle quantities
		store_remle[k,1:2] = theta_remle
		G = ((1-theta_remle[2])*covf(distmato/theta_remle[1]) + 
			theta_remle[2]*diag(no))
		Ginv = solve(G)
		covB = solve(t(Xo)%*%Ginv%*%Xo)
		covBXGi = covB%*%t(Xo)%*%Ginv
		Q = Ginv-Ginv%*%Xo%*%covBXGi
		sig2 = as.numeric(t(yo)%*%Q%*%yo/(no-p))
		Ginv = Ginv/sig2
		covB = sig2*covB
		store_remle[k,3] = sig2
		store_remle[k,4:6] = covBXGi%*%yo
		store_remle[k,7:9] = diag(covB)
		cc = sig2*(1-theta_remle[2])*covf(distmatou/theta_remle[1])
		preds = Xu %*% covBXGi %*% yo + t(cc) %*% 
			Ginv %*% (yo - Xo %*% covBXGi%*%yo)
		store_remle[k,10:11] = yu
		store_remle[k,12:13] = preds
		Gu = sig2*((1-theta_remle[2])*covf(distmatu/theta_remle[1]) + 
			theta_remle[2]*diag(nu))
		dd = Xu - t(cc) %*% Ginv %*% Xo
		store_remle[k,14:15] = diag(Gu - t(cc) %*% Ginv %*% cc + 
			dd %*% covB %*% t(dd))
	}
	cat("\n")

	#summarize simulation into useful results

	#bias for covariance paramters

	hold_results[1,2*(m-1)+1] = mean(store_mle[,1]) - rangep
	hold_results[1,2*(m-1)+2] = mean(store_remle[,1]) - rangep

	hold_results[2,2*(m-1)+1] = mean(store_mle[,2]) - ratiop
	hold_results[2,2*(m-1)+2] = mean(store_remle[,2]) - ratiop

	hold_results[3,2*(m-1)+1] = mean(store_mle[,3]) - sigmap
	hold_results[3,2*(m-1)+2] = mean(store_remle[,3]) - sigmap

	#rmse for covariance parameters

	hold_results[4,2*(m-1)+1] = sqrt(mean((store_mle[,1] - rangep)^2))
	hold_results[4,2*(m-1)+2] = sqrt(mean((store_remle[,1] - rangep)^2))

	hold_results[5,2*(m-1)+1] = sqrt(mean((store_mle[,2] - ratiop)^2))
	hold_results[5,2*(m-1)+2] = sqrt(mean((store_remle[,2] - ratiop)^2))

	hold_results[6,2*(m-1)+1] = sqrt(mean((store_mle[,3] - sigmap)^2))
	hold_results[6,2*(m-1)+2] = sqrt(mean((store_remle[,3] - sigmap)^2))

	#bias for autocovariances at specified distances

	cov01 = sigmap*((1-ratiop)*covf(0.01/rangep))
	hold_results[7,2*(m-1)+1] =
		mean(store_mle[,3]*((1- store_mle[,2])*covf(.01/store_mle[,1])) - cov01)
	hold_results[7,2*(m-1)+2] =
		mean(store_remle[,3]*((1- store_remle[,2])*covf(.01/store_remle[,1]))- cov01)

	cov5 = sigmap*((1-ratiop)*covf(0.5/rangep))
	hold_results[8,2*(m-1)+1] =
		mean(store_mle[,3]*((1- store_mle[,2])*covf(.5/store_mle[,1])) - cov5)
	hold_results[8,2*(m-1)+2] =
		mean(store_remle[,3]*((1- store_remle[,2])*covf(.5/store_remle[,1]))- cov5)

	cov14 = sigmap*((1-ratiop)*covf(sqrt(2)/rangep))
	hold_results[9,2*(m-1)+1] =
		mean(store_mle[,3]*((1- store_mle[,2])*covf(sqrt(2)/store_mle[,1])) - cov5)
	hold_results[9,2*(m-1)+2] =
		mean(store_remle[,3]*((1- store_remle[,2])*covf(sqrt(2)/store_remle[,1]))- cov5)

	#rmse for autocovariances at specified distances

	hold_results[10,2*(m-1)+1] =
		sqrt(mean((store_mle[,3]*(1- store_mle[,2])*covf(.01/store_mle[,1]) - cov01)^2))
	hold_results[10,2*(m-1)+2] =
		sqrt(mean((store_remle[,3]*(1- store_remle[,2])*covf(.01/store_remle[,1])- cov01)^2))

	hold_results[11,2*(m-1)+1] =
		sqrt(mean((store_mle[,3]*(1- store_mle[,2])*covf(.5/store_mle[,1]) - cov5)^2))
	hold_results[11,2*(m-1)+2] =
		sqrt(mean((store_remle[,3]*(1- store_remle[,2])*covf(.5/store_remle[,1])- cov5)^2))

	hold_results[12,2*(m-1)+1] =
		sqrt(mean((store_mle[,3]*(1- store_mle[,2])*covf(sqrt(2)/store_mle[,1]) - cov14)^2))
	hold_results[12,2*(m-1)+2] =
		sqrt(mean((store_remle[,3]*(1- store_remle[,2])*covf(sqrt(2)/store_remle[,1])- cov14)^2))

	#bias for fixed effects

	hold_results[13,2*(m-1)+1] = mean(store_mle[,4]) - beta0
	hold_results[13,2*(m-1)+2] = mean(store_remle[,4]) - beta0

	hold_results[14,2*(m-1)+1] = mean(store_mle[,5]) - beta1
	hold_results[14,2*(m-1)+2] = mean(store_remle[,5]) - beta1

	hold_results[15,2*(m-1)+1] = mean(store_mle[,6]) - beta2
	hold_results[15,2*(m-1)+2] = mean(store_remle[,6]) - beta2

	#rmse for fixed effects

	hold_results[16,2*(m-1)+1] = sqrt(mean((store_mle[,4] - beta0)^2))
	hold_results[16,2*(m-1)+2] = sqrt(mean((store_remle[,4] - beta0)^2))
	 
	hold_results[17,2*(m-1)+1] = sqrt(mean((store_mle[,5] - beta1)^2))
	hold_results[17,2*(m-1)+2] = sqrt(mean((store_remle[,5] - beta1)^2))

	hold_results[18,2*(m-1)+1] = sqrt(mean((store_mle[,6] - beta2)^2))
	hold_results[18,2*(m-1)+2] = sqrt(mean((store_remle[,6] - beta2)^2))

	#90% coverage for fixed effects
		
	hold_results[19,2*(m-1)+1] = 
		mean(store_mle[,4] - 1.645*sqrt(store_mle[,7]) < beta0 & 
		beta0 < store_mle[,4] + 1.645*sqrt(store_mle[,7]))
	hold_results[19,2*(m-1)+2] =
		mean(store_remle[,4] - 1.645*sqrt(store_remle[,7]) < beta0 & 
		beta0 < store_remle[,4] + 1.645*sqrt(store_remle[,7]))

	hold_results[20,2*(m-1)+1] = 
		mean(store_mle[,5] - 1.645*sqrt(store_mle[,8]) < beta1 & 
		beta1 < store_mle[,5] + 1.645*sqrt(store_mle[,8]))
	hold_results[20,2*(m-1)+2] =
		mean(store_remle[,5] - 1.645*sqrt(store_remle[,8]) < beta1 & 
		beta1 < store_remle[,5] + 1.645*sqrt(store_remle[,8]))

	hold_results[21,2*(m-1)+1] = 
		mean(store_mle[,6] - 1.645*sqrt(store_mle[,9]) < beta2 & 
		beta2 < store_mle[,6] + 1.645*sqrt(store_mle[,9]))
	hold_results[21,2*(m-1)+2] = 
		mean(store_remle[,6] - 1.645*sqrt(store_remle[,9]) < beta2 & 
		beta2 < store_remle[,6] + 1.645*sqrt(store_remle[,9]))

	#prediction bias

	hold_results[22,2*(m-1)+1] = mean(store_mle[,10] - store_mle[,12])
	hold_results[22,2*(m-1)+2] = mean(store_remle[,10] - store_remle[,12])

	hold_results[23,2*(m-1)+1] = mean(store_mle[,11] - store_mle[,13])
	hold_results[23,2*(m-1)+2] = mean(store_remle[,11] - store_remle[,13])

	#prediction rmspe

	hold_results[24,2*(m-1)+1] = sqrt(mean((store_mle[,10] - store_mle[,12])^2))
	hold_results[24,2*(m-1)+2] = sqrt(mean((store_remle[,10] - store_remle[,12])^2))

	hold_results[25,2*(m-1)+1] = sqrt(mean((store_mle[,11] - store_mle[,13])^2))
	hold_results[25,2*(m-1)+2] = sqrt(mean((store_remle[,11] - store_remle[,13])^2))

	#prediction 90% coverage

	hold_results[26,2*(m-1)+1] = 
		mean(store_mle[,10] - 1.645*sqrt(store_mle[,14]) < store_mle[,12] & 
		store_mle[,12] < store_mle[,10] + 1.645*sqrt(store_mle[,14]))
	hold_results[26,2*(m-1)+2] = 
		mean(store_remle[,10] - 1.645*sqrt(store_remle[,14]) < store_remle[,12] & 
		store_remle[,12] < store_remle[,10] + 1.645*sqrt(store_remle[,14]))

	hold_results[27,2*(m-1)+1] = 
		mean(store_mle[,11] - 1.645*sqrt(store_mle[,15]) < store_mle[,13] & 
		store_mle[,13] < store_mle[,11] + 1.645*sqrt(store_mle[,15]))
	hold_results[27,2*(m-1)+2] = 
		mean(store_remle[,11] - 1.645*sqrt(store_remle[,15]) < store_remle[,13] & 
		store_remle[,13] < store_remle[,11] + 1.645*sqrt(store_remle[,15]))

}
end.time = Sys.time()
difftime(end.time, start.time, units = 'min')
#save results to disk to be loaded later if needed
save(hold_results, file = 'hold_results_4maxrange.Rdata')
#load('hold_results_4maxrange.Rdata')

covparm_results = hold_results[1:12,]
covparm_results = data.frame(Measure = c('bias$^a$','bias$^b$','bias$^c$',
	'rmse$^a$','rmse$^b$','rmse$^c$','bias$^d$','bias$^e$','bias$^f$','rmse$^d$','rmse$^e$','rmse$^f$'), 
	covparm_results)

print(
    xtable(covparm_results, 
      align = c('l',rep('l', times = length(covparm_results[1,]))),
      digits = c(0,0,rep(3, times = 10)),
      caption = 'Covariates used for model-fitting',
      label = 'tab:covList'
    ),
    size = 'footnotesize',
    sanitize.text.function = identity,
    include.rownames = FALSE,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
 
fixeff_results = hold_results[13:21,]
fixeff_results = data.frame(Measure = c('bias$^a$','bias$^b$','bias$^c$',
	'rmse$^a$','rmse$^b$','rmse$^c$','cov90$^a$','cov90$^b$','cov90$^c$'), 
	fixeff_results)

print(
    xtable(fixeff_results, 
      align = c('l',rep('l', times = length(covparm_results[1,]))),
      digits = c(0,0,rep(3, times = 10)),
      caption = 'Covariates used for model-fitting',
      label = 'tab:covList'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)
 
pred_results = hold_results[22:27,]
pred_results = data.frame(Measure = c('bias$^a$','bias$^b$',
	'rmspe$^a$','rmspe$^b$','cov90$^a$','cov90$^b$'), 
	pred_results)

print(
    xtable(pred_results, 
      align = c('l',rep('l', times = length(covparm_results[1,]))),
      digits = c(0,0,rep(3, times = 10)),
      caption = 'Covariates used for model-fitting',
      label = 'tab:covList'
    ),
    size = 'footnotesize',
    include.rownames = FALSE,
    sanitize.text.function = identity,
    sanitize.rownames.function = identity,
    only.contents = TRUE,
    include.colnames = FALSE
)

