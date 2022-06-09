sec_path = 'Rcode/Chapter8/Section 8.1/'
setwd(paste0(SLEDbook_path,sec_path))

################################################################################
#-------------------------------------------------------------------------------
#              Multimodality Investigations
#-------------------------------------------------------------------------------
################################################################################

# (a) spherical case, used seeds 404, 405, and 406
# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K = 8
locx = matrix(rep(1:K,K),ncol=1,byrow=T)
locy = matrix(rep(1:K,each=K),ncol=1,byrow=T)
distmat = as.matrix(dist(cbind(rep(1:K,K),rep(1:K,each=K))))
n = K^2

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
	exp(-sqrt(3)*distance.matrix^2) + diag(rep(1e-11, times = dim(distance.matrix)[1])) 
}

# Set up the model matrix: constant but unknown mean
int1 = matrix(1, n, 1)
X = int1

# Create the orthogonal projection matrix onto the column space of X
P <- X%*%solve(t(X)%*%X,t(X))

sph.seed = 1330
exp.seed = 1350
gau.seed = 1370
theta.sph = 3
theta.exp = 3
theta.gau = 3
lwd_all = 4
# Simulate response from an SLM, with true intercept equal to 1 and spatial 
# range from spherical model.
G = spherical(distmat/theta.sph)
set.seed(sph.seed)
y1 = int1+t(chol(G))%*%rnorm(n)
# Evaluate the profile log-liklihood function at values of theta from 0.01 to 
# 8.00 (by 0.01).
hold1 = matrix(0,800,2)
for(i in 1:800){
  theta <- i/100
  G = spherical(distmat/theta)
  Ginv = solve(G)
  b = solve(t(X)%*%Ginv%*%X)%*%t(X)%*%Ginv%*%y1
  r = y1 - X%*%b 
  L = -(n/2)*log(t(r)%*%Ginv%*%r) - determinant(G, logarithm = TRUE)$modulus/2 
  hold1[i,1] = theta
  hold1[i,2] = L
}

hold = matrix(0,800,2)
for(i in 1:800){
  theta <- i/100
  G = spherical(distmat/theta)
Ginv <- solve(G)
Q <- Ginv-Ginv%*%X%*%solve(t(X)%*%Ginv%*%X)%*%t(X)%*%Ginv
L <- -0.5*log(det(G))-(n/2)*log(t(y1)%*%Q%*%y1)  
  hold[i,1] = theta
  hold[i,2] = L
}

# Simulate response from an SLM, with true intercept equal to 1 and spatial 
# range  from spherical model.
G = spherical(distmat/theta.sph)
set.seed(sph.seed + 1)
y1 = int1+t(chol(G))%*%rnorm(n)
# Evaluate the profile log-liklihood function at values of theta from 0.01 to 
# 8.00 (by 0.01).
hold2 = matrix(0,800,2)
for(i in 1:800){
  theta <- i/100
  G = spherical(distmat/theta)
  Ginv = solve(G)
  Q = Ginv - Ginv%*%X%*%solve(t(X)%*%Ginv%*%X, t(X))%*%Ginv
  L = -0.5*log(det(G)) - (n/2)*log(t(y1)%*%Q%*%y1)
  hold2[i,1] = theta
  hold2[i,2] = L
}

# Simulate response from an SLM, with true intercept equal to 1 and spatial 
# range from spherical model.
G = spherical(distmat/theta.sph)
set.seed(sph.seed + 2)
y1 = int1+t(chol(G))%*%rnorm(n)
# Evaluate the profile log-liklihood function at values of theta from 0.01 to 
# 8.00 (by 0.01).
hold3 = matrix(0,800,2)
for(i in 1:800){
  theta <- i/100
  G = spherical(distmat/theta)
  Ginv = solve(G)
  Q = Ginv - Ginv%*%X%*%solve(t(X)%*%Ginv%*%X, t(X))%*%Ginv
  L = -0.5*log(det(G)) - (n/2)*log(t(y1)%*%Q%*%y1)
  hold3[i,1] = theta
  hold3[i,2] = L
}

# Simulate response from an SLM, with true intercept equal to 1 and spatial 
# range from exponential model.

G = exponential(distmat/theta.exp)
set.seed(exp.seed)
y1 = int1+t(chol(G))%*%rnorm(n)
# Evaluate the profile log-liklihood function at values of theta from 0.01 to 
# 8.00 (by 0.01).
hold4 = matrix(0,800,2)
for(i in 1:800){
  theta <- i/100
  G = exponential(distmat/theta)
  Ginv = solve(G)
  Q = Ginv - Ginv%*%X%*%solve(t(X)%*%Ginv%*%X, t(X))%*%Ginv
  L = -0.5*log(det(G)) - (n/2)*log(t(y1)%*%Q%*%y1)
  hold4[i,1] = theta
  hold4[i,2] = L
}

# Simulate response from an SLM, with true intercept equal to 1 and spatial 
# range from exponential model.
G = exponential(distmat/theta.exp)
set.seed(exp.seed + 1)
y1 = int1+t(chol(G))%*%rnorm(n)
# Evaluate the profile log-liklihood function at values of theta from 0.01 to 
# 8.00 (by 0.01).
hold5 = matrix(0,800,2)
for(i in 1:800){
  theta <- i/100
  G = exponential(distmat/theta)
  Ginv = solve(G)
  Q = Ginv - Ginv%*%X%*%solve(t(X)%*%Ginv%*%X, t(X))%*%Ginv
  L = -0.5*log(det(G)) - (n/2)*log(t(y1)%*%Q%*%y1)
  hold5[i,1] = theta
  hold5[i,2] = L
}

# Simulate response from an SLM, with true intercept equal to 1 and spatial 
# range from exponential model.
G = exponential(distmat/theta.exp)
set.seed(exp.seed + 2)
y1 = int1+t(chol(G))%*%rnorm(n)
# Evaluate the profile log-liklihood function at values of theta from 0.01 to 
# 8.00 (by 0.01).
hold6 = matrix(0,800,2)
for(i in 1:800){
  theta <- i/100
  G = exponential(distmat/theta)
  Ginv = solve(G)
  Q = Ginv - Ginv%*%X%*%solve(t(X)%*%Ginv%*%X, t(X))%*%Ginv
  L = -0.5*log(det(G)) - (n/2)*log(t(y1)%*%Q%*%y1)
  hold6[i,1] = theta
  hold6[i,2] = L
}


# Simulate response from an SLM, with true intercept equal to 1 and spatial 
# range from Gaussian model.

G = Gaussian(distmat/theta.gau)
set.seed(gau.seed)
y1 = int1+t(chol(G))%*%rnorm(n)
# Evaluate the profile log-liklihood function at values of theta from 0.01 to 
# 8.00 (by 0.01).
hold7 = matrix(0,600,2)
for(i in 1:600){
  theta <- i/100
  G = Gaussian(distmat/theta)
  Ginv = solve(G)
  Q = Ginv - Ginv%*%X%*%solve(t(X)%*%Ginv%*%X, t(X))%*%Ginv
  L = -0.5*log(det(G)) - (n/2)*log(t(y1)%*%Q%*%y1)
  hold7[i,1] = theta
  hold7[i,2] = L
}
plot(hold7)

# Simulate response from an SLM, with true intercept equal to 1 and spatial 
# range from Gaussian model.
G = Gaussian(distmat/theta.gau)
set.seed(gau.seed + 1)
y1 = int1+t(chol(G))%*%rnorm(n)
# Evaluate the profile log-liklihood function at values of theta from 0.01 to 
# 8.00 (by 0.01).
hold8 = matrix(0,600,2)
for(i in 1:600){
  theta <- i/100
  G = Gaussian(distmat/theta)
  Ginv = solve(G)
  Q = Ginv - Ginv%*%X%*%solve(t(X)%*%Ginv%*%X, t(X))%*%Ginv
  L = -0.5*log(det(G)) - (n/2)*log(t(y1)%*%Q%*%y1)
  hold8[i,1] = theta
  hold8[i,2] = L
}

# Simulate response from an SLM, with true intercept equal to 1 and spatial 
# range from Gaussaian model.
G = Gaussian(distmat/theta.gau)
set.seed(gau.seed + 2)
y1 = int1+t(chol(G))%*%rnorm(n)
# Evaluate the profile log-liklihood function at values of theta from 0.01 to 
# 8.00 (by 0.01).
hold9 = matrix(0,600,2)
for(i in 1:600){
  theta <- i/100
  G = Gaussian(distmat/theta)
  Ginv = solve(G)
  Q = Ginv - Ginv%*%X%*%solve(t(X)%*%Ginv%*%X, t(X))%*%Ginv
  L = -0.5*log(det(G)) - (n/2)*log(t(y1)%*%Q%*%y1)
  hold9[i,1] = theta
  hold9[i,2] = L
}

file_name = "profileloglikelihoods"

pdf(paste0(file_name,'.pdf'), width = 8.5, height = 8.5)

layout(matrix(1:9, ncol = 3, byrow = TRUE), 
  widths = c(1.18,1,1),
  heights = c(1,1,1.11))
old.par = par(mar=c(2,6,5,2))
plot(hold1, type = 'l', xlab="", lwd = lwd_all ,
  ylab="Profile log-likelihood", cex.lab = 2, cex.axis = 1.5)
par(mar=c(2,1,5,2))
plot(hold2, type = 'l', xlab="", lwd = lwd_all ,
  ylab="", cex.lab = 2, cex.axis = 1.5,
  main = "Spherical", cex.main = 3)
par(mar=c(2,1,5,2))
plot(hold3, type = 'l', xlab="", lwd = lwd_all ,
  ylab="", cex.lab = 2, cex.axis = 1.5)

par(mar=c(2,6,5,2))
plot(hold4, type = 'l', xlab="", lwd = lwd_all ,
  ylab="Profile log-likelihood", cex.lab = 2, cex.axis = 1.5)
par(mar=c(2,1,5,2))
plot(hold5, type = 'l', xlab="", lwd = lwd_all ,
  ylab="", main = "Exponential", cex.main = 3, cex.axis = 1.5)
par(mar=c(2,1,5,2))
plot(hold6, type = 'l', xlab="", lwd = lwd_all ,
  ylab="", cex.axis = 1.5)

par(mar=c(5,6,5,2))
plot(hold7, type = 'l', xlab="(Effective) Range", lwd = lwd_all , 
  cex.lab = 2, cex.axis = 1.5, ylab="Profile log-likelihood")
par(mar=c(5,1,5,2))
plot(hold8, type = 'l', xlab="(Effective) Range", lwd = lwd_all ,
  ylab="", main = "Gaussian", cex.main = 3, cex.lab = 2, cex.axis = 1.5)
par(mar=c(5,1,5,2))
plot(hold9, type = 'l', xlab="(Effective) Range", lwd = lwd_all , 
  cex.lab = 2, cex.axis = 1.5, ylab="")

par(old.par)
dev.off()

system(paste0('pdfcrop ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('cp ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\' ','\'',SLEDbook_path,
  sec_path,file_name,'.pdf','\''))
system(paste0('rm ','\'',SLEDbook_path,
  sec_path,file_name,'-crop.pdf','\''))
