sec_path = 'Rcode/Chapter9/Section 9.2/'
setwd(paste0(SLEDbook_path,sec_path))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#           Spatial Configuration of Toy Example
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
file_name = 'figures/Moss_diffhist'
pdf(paste0(file_name,'.pdf'), width = 6, height = 10)

# Code for creating Figure 9.1
x = c(1,1,2,2,5)
y = c(3,4,3,2,4)
addx <- c(2,4,4)
addy <- c(4,2,4)
pdf("figure9pt1.pdf")
	plot(x,y,pch=16,xlim=c(1,5),ylim=c(1,5),xlab="",ylab="",cex=1.5)
	points(addx,addy,cex=1.5)
dev.off()

-----------------------------------------------------------------------------------------------------------

# Code for ordinary kriging, spherical model in Table 9.1
y <- matrix(c(1,0,3,1,5),5,1)
X <- matrix(1,5,1)
u <- matrix(c(1,1,2,2,5),5,1)
v <- matrix(c(4,3,3,2,4),5,1)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  if(dist<4){
   R[i,j] <- 1-(3*dist/8)+(dist^3/128)}
}}
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 if(dist<4){
  r0[i,1] <- 1-(3*dist/8)+(dist^3/128)}
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u02 <- 4
v02 <- 2
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u02)^2+(v[i,1]-v02)^2)
 if(dist<4){
  r0[i,1] <- 1-(3*dist/8)+(dist^3/128)}
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u03 <- 4
v03 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u03)^2+(v[i,1]-v03)^2)
 if(dist<4){
  r0[i,1] <- 1-(3*dist/8)+(dist^3/128)}
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#

# Code for ordinary kriging, exponential model in Table 9.1
y <- matrix(c(1,0,3,1,5),5,1)
X <- matrix(1,5,1)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- exp(-3*dist/4)
}}
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u02 <- 4
v02 <- 2
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u02)^2+(v[i,1]-v02)^2)
 r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u03 <- 4
v03 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u03)^2+(v[i,1]-v03)^2)
 r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#

# Code for ordinary kriging, Gaussian model in Table 9.1
y <- matrix(c(1,0,3,1,5),5,1)
X <- matrix(1,5,1)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- exp(-3*dist^2/16)
}}
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u02 <- 4
v02 <- 2
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u02)^2+(v[i,1]-v02)^2)
 r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u03 <- 4
v03 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u03)^2+(v[i,1]-v03)^2)
 r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#

-----------------------------------------------------------------------------------------------------------------------

# Code for ordinary kriging example in the context of a spatial-weights linear model (Table 9.2)
W0 <- matrix(c(0,1,0,1,0,1,0,1,0),nrow=3,byrow=T)
id3 <- diag(3)
ze3 <- matrix(0,3,3)
W1 <- cbind(W0,id3,ze3)
W2 <- cbind(id3,W0,id3)
W3 <- cbind(ze3,id3,W0)
W <- rbind(W1,W2,W3)
W <- W[-1,-1]
W <- W[-2,-2]
K <- diag(7)
SigmaSAR <- solve(diag(7)-0.3*W)%*%K%*%solve(diag(7)-0.3*t(W))
SigmaCAR <- solve(diag(7)-0.3*W)%*%K
corrSAR <- solve(diag(sqrt(diag(SigmaSAR))))%*%SigmaSAR%*%solve(diag(sqrt(diag(SigmaSAR))))
corrCAR <- solve(diag(sqrt(diag(SigmaCAR))))%*%SigmaCAR%*%solve(diag(sqrt(diag(SigmaCAR))))
R <- SigmaSAR[-1,-1]
r0 <- SigmaSAR[2:7,1]
r00 <- SigmaSAR[1,1]
X <- matrix(1,nrow=6,ncol=1,byrow=T)
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(6)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- r00-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
R <- SigmaSAR[-3,-3]
r0 <- SigmaSAR[-3,3]
r00 <- SigmaSAR[3,3]
X <- matrix(1,nrow=6,ncol=1,byrow=T)
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(6)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- r00-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
R <- SigmaSAR[-7,-7]
r0 <- SigmaSAR[-7,7]
r00 <- SigmaSAR[7,7]
X <- matrix(1,nrow=6,ncol=1,byrow=T)
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(6)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- r00-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#

---------------------------------------------------------------------------------------------------------------------

# Code for universal kriging example (with planar mean surface), spherical model in Table 9.3
y <- matrix(c(1,0,3,1,5),5,1)
X <- matrix(c(1,1,4,1,1,3,1,2,3,1,2,2,1,5,4),nrow=5,ncol=3,byrow=T)
u <- matrix(c(1,1,2,2,5),5,1)
v <- matrix(c(4,3,3,2,4),5,1)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  if(dist<4){
   R[i,j] <- 1-(3*dist/8)+(dist^3/128)}
}}
u01 <- 2
v01 <- 4
x0 <- matrix(c(1,u01,v01),nrow=3,ncol=1,byrow=T)
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 if(dist<4){
  r0[i,1] <- 1-(3*dist/8)+(dist^3/128)}
}
ukwts <- t(x0)%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
ukpev <- 1-t(r0)%*%solve(R)%*%r0+(t(x0)-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(x0-t(X)%*%solve(R)%*%r0)
#
u02 <- 4
v02 <- 2
x0 <- matrix(c(1,u02,v02),nrow=3,ncol=1,byrow=T)
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u02)^2+(v[i,1]-v02)^2)
 if(dist<4){
  r0[i,1] <- 1-(3*dist/8)+(dist^3/128)}
}
ukwts <- t(x0)%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
ukpev <- 1-t(r0)%*%solve(R)%*%r0+(t(x0)-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(x0-t(X)%*%solve(R)%*%r0)
#
u03 <- 4
v03 <- 4
x0 <- matrix(c(1,u03,v03),nrow=3,ncol=1,byrow=T)
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u03)^2+(v[i,1]-v03)^2)
 if(dist<4){
  r0[i,1] <- 1-(3*dist/8)+(dist^3/128)}
}
ukwts <- t(x0)%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
ukpev <- 1-t(r0)%*%solve(R)%*%r0+(t(x0)-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(x0-t(X)%*%solve(R)%*%r0)
#

--------------------------------------------------------------------------------------------------------------------------

# Code for Table 9.4 (ordinary kriging, exponential model of various strengths, nuggets, anisotropies)
y <- matrix(c(1,0,3,1,5),5,1)
X <- matrix(1,5,1)
# Model 1
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- exp(-3*dist/4)
}}
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
# Model 2 (effective range twice as large)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- exp(-3*dist/8)
}}
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- exp(-3*dist/8)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
# Model 3 (effective range half as large)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- exp(-3*dist/2)
}}
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- exp(-3*dist/2)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
# Model 4 (25% nugget, but same sill)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- 0.75*exp(-3*dist/4)
}}
R <- R+.25*diag(5)
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- 0.75*exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
# Model 5 (50% nugget, but same sill)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- 0.5*exp(-3*dist/4)
}}
R <- R+.5*diag(5)
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- 0.5*exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
# Model 6 (geometrically anisotropic with major axis in E-W direction)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2/4+(v[i,1]-v[j,1])^2)
  R[i,j] <- exp(-3*dist/4)
}}
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2/4+(v[i,1]-v01)^2)
 r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(5)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#

---------------------------------------------------------------------------------------------------------------------------

# IDW and IDSW for example of Section 9.6
y <- matrix(c(1,0,3,1,5),5,1)
X <- matrix(1,5,1)
# Model 1
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- exp(-3*dist/4)
}}
u01 <- 2
v01 <- 4
dist0 <- matrix(0,5,1)
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist0[i,1] <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- exp(-3*dist0[i,1]/4)
}
IDWwts <- (1/sum(1/dist0))*(1/dist0)
IDSWwts <- (1/sum(1/dist0^2))*(1/dist0^2)
IDWpev <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1
#

# Model 2 (effective range twice as large)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- exp(-3*dist/8)
}}
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- exp(-3*dist/8)
}
IDWpev <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1

# Model 3 (effective range half as large)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- exp(-3*dist/2)
}}
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- exp(-3*dist/2)
}
IDWpev <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1

# Model 4 (25% nugget, but same sill)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- 0.75*exp(-3*dist/4)
}}
R <- R+.25*diag(5)
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- 0.75*exp(-3*dist/4)
}
IDWpev <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1

# Model 5 (50% nugget, but same sill)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
  R[i,j] <- 0.5*exp(-3*dist/4)
}}
R <- R+.5*diag(5)
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
 r0[i,1] <- 0.5*exp(-3*dist/4)
}
IDWpev <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1

# Model 6 (geometrically anisotropic with major axis in E-W direction)
R <- matrix(0,5,5)
for(i in 1:5){
 for(j in 1:5){
  dist <- sqrt((u[i,1]-u[j,1])^2/4+(v[i,1]-v[j,1])^2)
  R[i,j] <- exp(-3*dist/4)
}}
u01 <- 2
v01 <- 4
r0 <- matrix(0,5,1)
for(i in 1:5){
 dist <- sqrt((u[i,1]-u01)^2/4+(v[i,1]-v01)^2)
 r0[i,1] <- exp(-3*dist/4)
}
IDWpev <- t(IDWwts)%*%R%*%IDWwts - 2*t(r0)%*%IDWwts + 1
IDSWpev <- t(IDSWwts)%*%R%*%IDSWwts - 2*t(r0)%*%IDSWwts + 1
#

------------------------------------------------------------------------------------------------------------------

# Code to generate Figure 9.5 (the negative-weights rainshadow)
library(gstat)
library(sp)
library(geoR)
KC_gau1 <- krige.control(type.krige = "ok", cov.model = "gaussian", cov.pars = c(1, sqrt(3)))
hor <- seq(1, 5, length.out = 1000)
ver <- seq(-2, 2, length.out = 1000)
gd <- expand.grid(x = hor, y = ver)
loc <- list()
for (i in 1 : 10 ^ 6){
    loc[[i]] <- matrix(c(1, 0, gd[i, 1], gd[i, 2]), ncol = 2, byrow = TRUE)
}
pred.loc <- c(0, 0)
pred.loc
gau.wts <- matrix(unlist(lapply(loc, function(x)krweights(x, pred.loc, KC_gau1))), ncol = 2,
                  byrow = TRUE)
pdf("Figure9pt5.pdf")
plot(NA, xlim = c(-1, 3), ylim = c(-2, 2), asp = 1, xlab = "", ylab = "")
grid(NULL, NULL, col = "lightgray", lty = "dotted")
points(gd[gau.wts[,2] < 0, 1], gd[gau.wts[,2] < 0, 2], col = "grey")
points(pred.loc[1], pred.loc[2], pch = 21, cex = 1.5, bg = "yellow")
text(pred.loc[1], pred.loc[2], pos = 2, cex = 1.1, labels = expression(bold('S')[0]))
points(1, 0, pch = 21, cex = 1.5, bg = "red")
dev.off()

------------------------------------------------------------------------------------------------------------------------

# Code for illustrating the screening effect, exponential model (Table 9.5)
X <- matrix(1,2,1)
u <- matrix(c(-1,1.1))
v <- matrix(c(0,0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist/4)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(-sqrt(2)/2,1.1))
v <- matrix(c(sqrt(2)/2,0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist/4)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(0,1.1))
v <- matrix(c(1,0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist/4)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(sqrt(2)/2,1.1))
v <- matrix(c(sqrt(2)/2,0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist/4)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(cos(pi/8),1.1))
v <- matrix(c(sin(pi/8),0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist/4)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(cos(pi/16),1.1))
v <- matrix(c(sin(pi/16),0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist/4)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(cos(pi/32),1.1))
v <- matrix(c(sin(pi/32),0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist/4)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(1,1.1))
v <- matrix(c(0,0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist/4)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist/4)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
--------------------------------------------------------------------------------------------------------------------

# Code for plotting Figure 9.4
fixedsite <- matrix(c(1.1,0),nrow=1,ncol=2,byrow=T)
predsite <- matrix(c(0,0),nrow=1,ncol=2,byrow=T)
transient <- matrix(c(-1,0,-sqrt(2)/2,sqrt(2)/2,0,1,sqrt(2)/2,sqrt(2)/2,cos(pi/8),sin(pi/8),
cos(pi/16),sin(pi/16),cos(pi/32),sin(pi/32),1,0),nrow=8,ncol=2,byrow=T)
plot(asp=1,transient[,1],transient[,2],xlim=c(-1.2,1.2),xlab="",ylab="",pch=16)
text(transient[,1]-0.1,transient[,2],labels=c(1:8))
points(fixedsite,pch=1)
points(predsite,pch=4)
#
---------------------------------------------------------------------------------------------------------------------

# Code for illustrating the screening effect, Gaussian model
X <- matrix(1,2,1)
u <- matrix(c(-1,1.1))
v <- matrix(c(0,0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist^2/16)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(-sqrt(2)/2,1.1))
v <- matrix(c(sqrt(2)/2,0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist^2/16)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(0,1.1))
v <- matrix(c(1,0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist^2/16)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(sqrt(2)/2,1.1))
v <- matrix(c(sqrt(2)/2,0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist^2/16)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(cos(pi/8),1.1))
v <- matrix(c(sin(pi/8),0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist^2/16)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(cos(pi/16),1.1))
v <- matrix(c(sin(pi/16),0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist^2/16)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(cos(pi/32),1.1))
v <- matrix(c(sin(pi/32),0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist^2/16)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#
u <- matrix(c(1,1.1))
v <- matrix(c(0,0))
R <- diag(2)
for(i in 1:2){
 for(j in 1:2){
  dist <- sqrt((u[i,1]-u[j,1])^2+(v[i,1]-v[j,1])^2)
R[i,j] <- exp(-3*dist^2/16)
}}
u01 <- 0
v01 <- 0
r0 <- matrix(0,2,1)
for(i in 1:2){
 dist <- sqrt((u[i,1]-u01)^2+(v[i,1]-v01)^2)
r0[i,1] <- exp(-3*dist^2/16)
}
okwts <- solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R)+t(r0)%*%solve(R)%*%(diag(2)-X%*%solve(t(X)%*%solve(R)%*%X)%*%t(X)%*%solve(R))
okpev <- 1-t(r0)%*%solve(R)%*%r0+(1-t(r0)%*%solve(R)%*%X)%*%solve(t(X)%*%solve(R)%*%X)%*%(1-t(X)%*%solve(R)%*%r0)
#

# Code for comparing quality of predictions based on ML and REML estimation (Table 9.6)
# Set up a small example of n=K^2 sites on a KxK square grid with unit spacing
K <- 8
locx <- matrix(rep(1:K,K),ncol=1,byrow=T)
locy <- matrix(rep(1:K,each=K),ncol=1,byrow=T)
n <- K^2
theta <- 0.7

# Define a covariance structure on the grid: an exponential covariance function with correlation theta at
# unit distance
Gtrue <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  Gtrue[i,j] <- theta^d
}}
Gtrueaug <- matrix(0,n+2,n+2)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  Gtrueaug[i,j] <- theta^d
}}
for(i in 1:n){
 d1 <- sqrt((locx[i,1]-4.5)^2+(locy[i,1]-4.5)^2)
 d2 <- sqrt((locx[i,1]-4.5)^2+(locy[i,1]-1)^2)
 Gtrueaug[i,65] <- theta^d1
 Gtrueaug[65,i] <- Gtrueaug[i,65]
 Gtrueaug[i,66] <- theta^d2
 Gtrueaug[66,i] <- Gtrueaug[i,66]
}
Gtrueaug[65,65] <- 1
Gtrueaug[66,66] <- 1
Gtrueaug[65,66] <- theta^3.5
Gtrueaug[66,65] <- theta^3.5

# Set up the model matrix: I, constant but unknown mean; II, row effects; III, row and column effects
int1 <- matrix(1,n,1)
int1aug <- matrix(1,n+2,1)
set.seed(806)
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

mle <- matrix(0,100,16)
remle <- matrix(0,100,14)
countml1 <- 0
countml2 <- 0
countml3 <- 0
countreml1 <- 0
countreml2 <- 0
countreml3 <- 0

for(k in 1:100){
y1aug <- int1aug+t(chol(Gtrueaug))%*%rnorm(n+2)
y1 <- y1aug[1:64,1]
# y1 <- t(chol(Gtrue))%*%rnorm(n)
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
 if(hold[i,2]==maxL){mle[k,3] <- (hold[i,4]/n)/log(1/hold[i,1])}
}
maxLR <- max(hold[,3])
for(i in 1:999){
 if(hold[i,3]==maxLR){remle[k,1] <- hold[i,1]}
 if(hold[i,3]==maxLR){remle[k,2] <- hold[i,4]/(n-p)}
 if(hold[i,3]==maxLR){remle[k,3] <- (hold[i,4]/(n-p))/log(1/hold[i,1])}
}
thetahatml <- mle[k,1]
sigma2hatml <- mle[k,2]
Ghatml <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  Ghatml[i,j] <- thetahatml^d
}}
g0hat1 <- matrix(0,n,1)
g0hat2 <- matrix(0,n,1)
for(i in 1:n){
 d1 <- sqrt((locx[i,1]-4.5)^2+(locy[i,1]-4.5)^2)
 d2 <- sqrt((locx[i,1]-4.5)^2+(locy[i,1]-1)^2)
 g0hat1[i,1] <- thetahatml^d1
 g0hat2[i,1] <- thetahatml^d2
}
Ghatinvml <- solve(Ghatml)
betahatml <- solve(t(X)%*%Ghatinvml%*%X)%*%t(X)%*%Ghatinvml%*%y1
covbetahatml <- sigma2hatml*solve(t(X)%*%Ghatinvml%*%X)
OK1 <- betahatml[1,1]+t(g0hat1)%*%Ghatinvml%*%(y1-X%*%betahatml[1,1])
OK2 <- betahatml[1,1]+t(g0hat2)%*%Ghatinvml%*%(y1-X%*%betahatml[1,1])
mle[k,4] <- betahatml[1,1]-qt(.90,n-p)*sqrt(covbetahatml[1,1])
mle[k,5] <- betahatml[1,1]+qt(.90,n-p)*sqrt(covbetahatml[1,1])
mle[k,6] <- mle[k,5]-mle[k,4]
mle[k,7] <- OK1
mle[k,8] <- OK2
covOK1 <- sigma2hatml*(1-t(g0hat1)%*%Ghatinvml%*%g0hat1+(1-t(g0hat1)%*%Ghatinvml%*%X)%*%solve(t(X)%*%Ghatinvml%*%X)%*%(1-t(X)%*%Ghatinvml%*%g0hat1))
covOK2 <- sigma2hatml*(1-t(g0hat2)%*%Ghatinvml%*%g0hat2+(1-t(g0hat2)%*%Ghatinvml%*%X)%*%solve(t(X)%*%Ghatinvml%*%X)%*%(1-t(X)%*%Ghatinvml%*%g0hat2))
mle[k,9] <- OK1-qt(.90,n-p)*sqrt(covOK1)
mle[k,10] <- OK1+qt(.90,n-p)*sqrt(covOK1)
mle[k,11] <- mle[k,10]-mle[k,9]
mle[k,12] <- OK2-qt(.90,n-p)*sqrt(covOK2)
mle[k,13] <- OK2+qt(.90,n-p)*sqrt(covOK2)
mle[k,14] <- mle[k,13]-mle[k,12]
mle[k,15] <- y1aug[65,1]
mle[k,16] <- y1aug[66,1]
if(mle[k,4]<1 && mle[k,5]>1){countml1 <- countml1+1}
if(mle[k,9]<y1aug[65,1] && mle[k,10]>y1aug[65,1]){countml2 <- countml2+1}
if(mle[k,12]<y1aug[66,1] && mle[k,13]>y1aug[66,1]){countml3 <- countml3+1}
thetahatreml <- remle[k,1]
sigma2hatreml <- remle[k,2]
Ghatreml <- matrix(0,n,n)
for(i in 1:n){
 for(j in 1:n){
  d <- sqrt((locx[i,1]-locx[j,1])^2+(locy[i,1]-locy[j,1])^2)
  Ghatreml[i,j] <- thetahatreml^d
}}
g0hatreml1 <- matrix(0,n,1)
g0hatreml2 <- matrix(0,n,1)
for(i in 1:n){
 d1 <- sqrt((locx[i,1]-4.5)^2+(locy[i,1]-4.5)^2)
 d2 <- sqrt((locx[i,1]-4.5)^2+(locy[i,1]-1)^2)
 g0hatreml1[i,1] <- thetahatreml^d1
 g0hatreml2[i,1] <- thetahatreml^d2
}
Ghatinvreml <- solve(Ghatreml)
betahatreml <- solve(t(X)%*%Ghatinvreml%*%X)%*%t(X)%*%Ghatinvreml%*%y1
covbetahatreml <- sigma2hatreml*solve(t(X)%*%Ghatinvreml%*%X)
remle[k,4] <- betahatreml[1,1]-qt(.90,n-p)*sqrt(covbetahatreml[1,1])
remle[k,5] <- betahatreml[1,1]+qt(.90,n-p)*sqrt(covbetahatreml[1,1])
remle[k,6] <- remle[k,5]-remle[k,4]
OKreml1 <- betahatreml[1,1]+t(g0hatreml1)%*%Ghatinvreml%*%(y1-X%*%betahatreml[1,1])
OKreml2 <- betahatreml[1,1]+t(g0hatreml2)%*%Ghatinvreml%*%(y1-X%*%betahatreml[1,1])
remle[k,4] <- betahatreml[1,1]-qt(.90,n-p)*sqrt(covbetahatreml[1,1])
remle[k,5] <- betahatreml[1,1]+qt(.90,n-p)*sqrt(covbetahatreml[1,1])
remle[k,6] <- remle[k,5]-remle[k,4]
remle[k,7] <- OKreml1
remle[k,8] <- OKreml2
covOKreml1 <- sigma2hatreml*(1-t(g0hatreml1)%*%Ghatinvreml%*%g0hatreml1+(1-t(g0hatreml1)%*%Ghatinvreml%*%X)%*%solve(t(X)%*%Ghatinvreml%*%X)%*%(1-t(X)%*%Ghatinvreml%*%g0hatreml1))
covOKreml2 <- sigma2hatreml*(1-t(g0hatreml2)%*%Ghatinvreml%*%g0hatreml2+(1-t(g0hatreml2)%*%Ghatinvreml%*%X)%*%solve(t(X)%*%Ghatinvreml%*%X)%*%(1-t(X)%*%Ghatinvreml%*%g0hatreml2))
remle[k,9] <- OKreml1-qt(.90,n-p)*sqrt(covOKreml1)
remle[k,10] <- OKreml1+qt(.90,n-p)*sqrt(covOKreml1)
remle[k,11] <- remle[k,10]-remle[k,9]
remle[k,12] <- OKreml2-qt(.90,n-p)*sqrt(covOKreml2)
remle[k,13] <- OKreml2+qt(.90,n-p)*sqrt(covOKreml2)
remle[k,14] <- remle[k,13]-remle[k,12]
if(remle[k,4]<1 && remle[k,5]>1){countreml1 <- countreml1+1}
if(remle[k,9]<y1aug[65,1] && remle[k,10]>y1aug[65,1]){countreml2 <- countreml2+1}
if(remle[k,12]<y1aug[66,1] && remle[k,13]>y1aug[66,1]){countreml3 <- countreml3+1}
}


