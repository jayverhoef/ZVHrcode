temp1 <- seq(0, 1/3, length.out = 6)[2:5]
x1 <- rep(temp1, 4)
y1 <- rep(temp1, rep(4, 4))


temp2 <- seq(0, 1/2, length.out = 9)[2:8]
x2 <- rep(temp2, 7)
y2 <- rep(temp2, rep(7, 7))


temp3 <- seq(0, 1, length.out = 15)[2:14]
x3 <- rep(temp3, 13)
y3 <- rep(temp3, rep(13, 13))

temp4 <- seq(0+0.05, 1-0.05, length.out = 4)
x4 <- rep(temp4, 4)
y4 <- rep(temp4, rep(4, 4))


temp5 <- seq(0+0.05, 1-0.05, length.out = 7)
x5 <- rep(temp5, 7)
y5 <- rep(temp5, rep(7, 7))


temp6 <- seq(0, 1, length.out = 15)[2:14]
x6 <- rep(temp3, 13)
y6 <- rep(temp3, rep(13, 13))


par(mfrow=c(2,3)) 
plot(x1, y1, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
     xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
plot(x2, y2, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
     xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
plot(x3, y3, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
     xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
plot(x4, y4, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
     xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
plot(x5, y5, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
     xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)
plot(x6, y6, yaxt='n', xaxt='n', xaxs = "i", yaxs = "i", ann=FALSE, 
     xlim = c(0,1), ylim = c(0,1), pch = 20, cex = 2)

