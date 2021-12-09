-------------------------- Multimodality investigations for Figure 8.1 ----------------------------

# (a) spherical case, used seeds 404, 405, and 406
# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 8
locx <- matrix(rep(1:K,K),ncol=1,byrow=T)
locy <- matrix(rep(1:K,each=K),ncol=1,byrow=T)
n <- K^2
theta <- 3

# Define a covariance structure on the grid: a spherical covariance function with range equal to 3
G <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  if(d<theta){G[i,j] <- 1-(3*d/(2*theta))+(d^3/(2*theta^3))}
}}

# Set up the model matrix: constant but unknown mean
int1 <- matrix(1,n,1)
set.seed(406)
X <- int1

# Simulate response from an SLM, with true intercept equal to 1.
y1 <- int1+t(chol(G))%*%rnorm(n)
y1

# Create the orthogonal projection matrix onto the column space of X
P <- X%*%solve(t(X)%*%X)%*%t(X)

# Evaluate the profile log-liklihood function at values of theta from 0.01 to 8.00 (by 0.01).
y1 <- int1+t(chol(G))%*%rnorm(n)
hold <- matrix(0,800,2)
for(m in 1:800){
theta <- m/100
G <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  if(d<theta){G[i,j] <- 1-(3*d/(2*theta))+(d^3/(2*theta^3))}
}}
Ginv <- solve(G)
Q <- Ginv-Ginv%*%X%*%solve(t(X)%*%Ginv%*%X)%*%t(X)%*%Ginv
L <- -0.5*log(det(G))-(n/2)*log(t(y1)%*%Q%*%y1)
hold[m,1] <- theta
hold[m,2] <- L
}
plot(hold[,1],hold[,2],xlab="Range",ylab="Profile log-likelihood")


# (b) exponential case, used seeds 401, 402, and 403
# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 8
locx <- matrix(rep(1:K,K),ncol=1,byrow=T)
locy <- matrix(rep(1:K,each=K),ncol=1,byrow=T)
n <- K^2
theta <- 3

# Define a covariance structure on the grid: an exponential function with effective range equal to 3
G <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  G[i,j] <- exp(-d/(theta/3))
}}

# Set up the model matrix: constant but unknown mean
int1 <- matrix(1,n,1)
set.seed(403)
X <- int1

# Simulate response from an SLM, with true intercept equal to 1.
y1 <- int1+t(chol(G))%*%rnorm(n)
y1

# Create the orthogonal projection matrix onto the column space of X
P <- X%*%solve(t(X)%*%X)%*%t(X)

# Evaluate the log-liklihood function at values of theta from 0.01 to 8.00 (by 0.01).
y1 <- int1+t(chol(G))%*%rnorm(n)
hold <- matrix(0,800,2)
for(m in 1:800){
theta <- m/100
G <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  G[i,j] <- exp(-d/(theta/3))
}}
Ginv <- solve(G)
Q <- Ginv-Ginv%*%X%*%solve(t(X)%*%Ginv%*%X)%*%t(X)%*%Ginv
L <- -0.5*log(det(G))-(n/2)*log(t(y1)%*%Q%*%y1)
hold[m,1] <- theta
hold[m,2] <- L
}
plot(hold[,1],hold[,2],xlab="Effective range",ylab="Profile log-likelihood")

# (c) Gaussian case, used seeds 407, 408, and 409
# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 8
locx <- matrix(rep(1:K,K),ncol=1,byrow=T)
locy <- matrix(rep(1:K,each=K),ncol=1,byrow=T)
n <- K^2
theta <- 3

# Define a covariance structure on the grid: a Gaussian covariance function with effective range equal to 3
G <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  G[i,j] <- exp(-d^2/theta)
}}

# Set up the model matrix: constant but unknown mean
int1 <- matrix(1,n,1)
set.seed(409)
X <- int1

# Simulate response from an SLM, with true intercept equal to 1.
y1 <- int1+t(chol(G))%*%rnorm(n)
y1

# Create the orthogonal projection matrix onto the column space of X
P <- X%*%solve(t(X)%*%X)%*%t(X)

# Evaluate the log-liklihood function at values of theta from 0.01 to 8.00 (by 0.01).
y1 <- int1+t(chol(G))%*%rnorm(n)
hold <- matrix(0,800,2)
for(m in 1:800){
theta <- m/100
G <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  G[i,j] <- exp(-d^2/theta)
}}
Ginv <- solve(G)
Q <- Ginv-Ginv%*%X%*%solve(t(X)%*%Ginv%*%X)%*%t(X)%*%Ginv
L <- -0.5*log(det(G))-(n/2)*log(t(y1)%*%Q%*%y1)
hold[m,1] <- theta
hold[m,2] <- L
}
plot(hold[,1],hold[,2],xlab="Effective range",ylab="Profile log-likelihood")



----------------------- ML versus REML for Table 8.1 -----------------------------------------------------------

# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 8
locx <- matrix(rep(1:K,K),ncol=1,byrow=T)
locy <- matrix(rep(1:K,each=K),ncol=1,byrow=T)
n <- K^2
theta <- 0.5

# Define a covariance structure on the grid: an exponential covariance function with correlation theta at
# unit distance
Gtrue <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  Gtrue[i,j] <- theta^d
}}

# Set up the model matrix: I, constant but unknown mean; II, row effects; III, row and column effects
int1 <- matrix(1,n,1)
set.seed(406)
X <- int1
p <- 1
# set.seed(407)
# X <- diag(K)%x%matrix(1,K,1)
# p <- 8
# set.seed(408)
# X1 <- diag(K)%x%matrix(1,K,1)
# X2 <- matrix(1,K,1)%x%diag(K)
# X2 <- X2[,1:K-1]
# X <- cbind(X1,X2)
# p <- 15

# Create the orthogonal projection matrix onto the column space of X
P <- X%*%solve(t(X)%*%X)%*%t(X)

# Simulate response from an SLM, with true intercept equal to 1, and evaluate the log-liklihood function at 
# values of theta from 0.001 to 0.999 (by 0.001).
mle <- matrix(0,100,2)
remle <- matrix(0,100,2)
count <- 0
for(k in 1:100){
y1 <- int1+t(chol(Gtrue))%*%rnorm(n)
hold <- matrix(0,999,4)
for(m in 1:999){
theta <- m/1000
G <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  G[i,j] <- theta^d
}}
Ginv <- solve(G)
Q <- Ginv-Ginv%*%X%*%solve(t(X)%*%Ginv%*%X)%*%t(X)%*%Ginv
L <- -0.5*log(det(G))-(n/2)*log(t(y1)%*%Q%*%y1)
LR <- -0.5*log(det(G))-0.5*log(det(t(X)%*%Ginv%*%X))-((n-p)/2)*log(t(y1)%*%Q%*%y1)
hold[m,1] <- theta
hold[m,2] <- L
hold[m,3] <- LR
hold[m,4] <- t(y1)%*%Q%*%y1
}
maxL <- max(hold[,2])
for(i in 1:999){
 if(hold[i,2]==maxL){mle[k,1] <- hold[i,1]}
 if(hold[i,2]==maxL){mle[k,2] <- hold[i,4]/n}
}
maxLR <- max(hold[,3])
for(i in 1:999){
 if(hold[i,3]==maxLR){remle[k,1] <- hold[i,1]}
 if(hold[i,3]==maxLR){remle[k,2] <- hold[i,4]/(n-p)}
}
if(remle[k,1]==0.999){remle[k,1] <- mle[k,1]; remle[k,2] <- mle[k,2]; count <- count+1}
}

theta=0.5, seed=406, K=8, overall constant mean structure:
> summary(mle[,1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0700  0.3362  0.4535  0.4386  0.5312  0.7670 
> summary(mle[,2])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.3942  0.7304  0.8779  0.9450  1.0731  2.8544 
> summary(remle[,1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0920  0.3772  0.5075  0.4970  0.5895  0.9990 
> summary(remle[,2])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  0.4106   0.8062   0.9748   4.5176   1.2294 333.4272 
# Empirical bias results
> mean(mle[,1]-0.5)
[1] -0.06137
> mean(mle[,2]-1.0)
[1] -0.05503914
> mean(remle[,1]-0.5)
[1] -0.00296
> mean(remle[,2]-1.0)
[1] 3.517614
# Empirical MSE results
> mean((mle[,1]-0.5)^2)
[1] 0.02398119
> mean((mle[,2]-1.0)^2)
[1] 0.1187557
> mean((remle[,1]-0.5)^2)
[1] 0.02710262
> mean((remle[,2]-1.0)^2)
[1] 1106.725


theta=0.5, seed=407, K=8, row effects mean structure:
> summary(mle[,1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0410  0.1975  0.2970  0.3170  0.4310  0.7270 
> summary(mle[,2])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.3308  0.5563  0.6200  0.6889  0.8026  1.6074 
> summary(remle[,1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.1290  0.3070  0.4155  0.4541  0.5913  0.9990 
> summary(remle[,2])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  0.4387   0.7150   0.8563   6.3434   1.1715 528.1425 
# Empirical bias results
> mean(mle[,1]-0.5)
[1] -0.18302
> mean(mle[,2]-1.0)
[1] -0.3110759
> mean(remle[,1]-0.5)
[1] -0.04589
> mean(remle[,2]-1.0)
[1] 5.343425rem
# Empirical MSE results
> mean((mle[,1]-0.5)^2)
[1] 0.05703284
> mean((mle[,2]-1.0)^2)
[1] 0.1423976
> mean((remle[,1]-0.5)^2)
[1] 0.03910949
> mean((remle[,2]-1.0)^2)
[1] 2779.187


theta=0.5, seed=408, K=8, row and column effects mean structure:
> summary(mle[,1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00100 0.07425 0.20300 0.19873 0.31075 0.55900 
> summary(mle[,2])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.2560  0.4021  0.4871  0.4901  0.5622  0.8596 
> summary(remle[,1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0590  0.3048  0.4650  0.4790  0.6335  0.9990 
> summary(remle[,2])
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
  0.3442   0.6540   0.9178  38.0879   1.2798 617.9802 
# Empirical bias results
> mean(mle[,1]-0.5)
[1] -0.30127
> mean(mle[,2]-1.0)
[1] -0.5098685
> mean(remle[,1]-0.5)
[1] -0.02096
> mean(remle[,2]-1.0)
[1] 37.0879
# Empirical MSE results
> mean((mle[,1]-0.5)^2)
[1] 0.1105446
> mean((mle[,2]-1.0)^2)
[1] 0.2752327
> mean((remle[,1]-0.5)^2)
[1] 0.0537609
> mean((remle[,2]-1.0)^2)
[1] 19920.13


> summary(mle[,1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
0.00100 0.07425 0.20300 0.19873 0.31075 0.55900 
> summary(mle[,2])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.2560  0.4021  0.4871  0.4901  0.5622  0.8596 
> summary(remle[,1])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.0590  0.3048  0.4355  0.4408  0.5845  0.8270 
> summary(remle[,2])
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.3442  0.6427  0.8362  0.9650  1.1716  3.2908 
> mean((remle[,1]-0.5)^2)
[1] 0.0368706
> mean((remle[,2]-1.0)^2)
[1] 0.2192062


--------- Consistency under increasing-domain versus fixed-domain asymptotic regimes for Table 8.2 ---------

# Code written by Hongbeng Lim (a student in my spatial stats class)

library(gstat)
library(geoR)
library(sp)
library(parallel)
library(tictoc)

#creates exponential covariance matrix
#n is number of points in interval of length 1
#lenint is length of spatial interval
makeexponcov <- function(sill=1,range=0.1,n=50,lenint=1){
 totaln <- n*lenint
 sill <- max(0,sill)
 if(range <= 0) return(sill*diag(totaln))
 row1 <- sill*exp(0:(-(totaln-1))/(range*(n+1)))
 toeplitz(row1)
}

#exponential covariance model likelihood for multivariate normal used in MLE procedure
#setting usemean=TRUE uses GLS mean in estimating covariance: else mean assumed 0
likexpon <- function(par=c(1,0.1),dat,n,lenint,usemean=FALSE){
 totaln <- n*lenint
 if(length(dat)!=totaln){
  stop("data does not match length of n and interval")
 }
 if(par[1]<=0){
  return(Inf)
 }
 exponcov <- makeexponov(sill=par[1],range=par[2],n=n,lenint=lenint)
 invmat <- solve(exponcov)
 if(usemean){
  glsest <- sum(invmat%*%dat)/sum(invmat)
  y <- dat-glsest
 } else {
  y <- dat
 }
 determinant(exponcov)$modulus[1]+y %*% (invmat %*% y)
}

#calculates MLE for sill and range
mlenorm <- function(par=c(1,0.1),dat,n=50,lenint=1,usemean=TRUE){
 likfun <- function(par) likexpon(par=par,dat=dat,n=n,lenint=lenint,usemean=usemean)
 optim(par=par,fn=likfun)$par
}

#simulation for MLE in asymptotic domain
simasymp2 <- function(M=10,lenint=1,n=50,currentrange=0.1,usemean=FALSE){
 totaln <- lenint*n
 x <- 1:(totaln)/(n+1)
 xy <- cbind(x,y=0)
 simdata <- grf(n=totaln,grid=xy,xlims=c(0,lenint),ylims=c(0,1),
  messages=FALSE,nsim=M,cov.model="exponential",cov.pars=c(1,currentrange))
 estcov <- matrix(nrow=M,ncol=2)
 if(M>1){
  for(i in 1:M){
   tic()
   estcov[i,]<- mlenorm(data=simdata$data[,i],n=n,lenint=lenint,par=c(1,0.1),usemean=usemean)
   toc()
  }
} else{
  tic()
  estcov[1,] <- mlenorm(dat=simdata$data[,i]n=n,lenint=lenint,par=c(1,0.1),usemean=usemean)
  toc()
}
estcov
}

#parallel wrapper for above fn
par.simasymp2 <- function(M=20,lenint=1,n=50,currentrange=0.1){
 RR <- distribute(M,nnodes)
 dat <- parLapply(cl=cl,RR,fun=simasymp2,lenint=lenint,n=n,currentrange=currentrange)
 do.call(rbind,dat)
}

#sets up parallel framework
nnodes <- detectCores()
cl <- makeCluster(nnodes)
clusterEvalQ(cl, {
 library(gstat)
 library(geoR)
 library(sp)
 library(tictoc)
})
clusterExport(cl=cl,varlist=c("likexpon","makeexponcov","mlenorm"))

$fixed-domain asympotic simulation
myn <- c(50,250, 1000)
myranges <- c(0.1,0.2)
listn <- list()
for(i in 1:3){
 for(j in 1:2){
  tic()
  listn[[2*(i-1)+j]] <- par.simasymp2(M=1000,lenint=1,n=myn[i],currentrange=myranges[j])
  toc()
  }
}

#increasing-domain simulation
myl <- c(1,5,20)
listl <- list()
for(i in 1:3){
 for(j in 1:2){
  tic()
  listl[[2*(i-1)+j]] <- par.simasymp2(M=1000,lenint=myl[i],n=50,currentrange=myranges[j])
  toc()
 }
}
stopCluster(cl)

#calculates bias, variance, and MSE for MLE
getstats2 <- function(mat,range=0.1){
 newmat <- cbind(mat,mat[,1]/mat[,2])
 varvec <- apply(newmat,2,var)
 means <- c(1,range,1/range)
 demeanedmat <- sweep(newmat,2,means)
 biasvec <- colMeans(demeanedmat)
 msevec <- colMeans(demeanedmat^2)
 rbind(biasvec,varvec,msevec)
}



