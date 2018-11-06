## Remove the strong lines and get the continumm

SNe <- read.table("../Data/dbf4.txt",skip = 20)
colnames(SNe) <- c("SN","Lambda","F")
D6_1 <- SNe[SNe[,1]=="D6-1",]
plot(D6_1[,2],D6_1[,3],"l")

## Linear interpolation
Obj <- D6_1
LambdaSeq <- Obj[,2]
Lambdamin <- min(LambdaSeq)
Lambdamax <- max(LambdaSeq)
FSeq <- Obj[,3]
FluxFun <- splinefun(LambdaSeq,FSeq,method ="fmm")
TestLambda <- seq(from = Lambdamin,to = Lambdamax,length.out = 2^11)
TestF <- FluxFun(TestLambda)
lines(TestLambda,TestF,"l",col = 2)

library(SpecMatch)
Rcpp::sourceCpp('~/Dropbox/project/Wang/SpecMatch/src/WTtools.cpp')
res <- StarletWT(TestF)
lines(TestLambda, StarletRC(res),col = 2)
res[12,] <-  0
res[11,] <- 0
res[10,]<-0
res[9,]<-0
res[8,]<-0
res[7,]<-0
res[6,]<-0
#
res[5,]<-0
res[4,]<- 0
res[3,] <- 0
res[2,]<-0
res[1,] <- 0
plot(D6_1[,2],D6_1[,3],"l")
lines(TestLambda,StarletRC(res),"l",lwd = 2,col=3)
lines(TestLambda,rep(0,2048),lwd = 2)
res$C[res$C > M]



###################################################
resraw <- StarletWT(TestF)
sigma <- MultiResSuppStarlet(resraw,2)
res <- HardThreshold(resraw,sigma)
resraw[6:12,] = 0
res[6:12,] = 0
plot(TestLambda,StarletRC(resraw),"l")
lines(TestLambda,StarletRC(res),col = 2,lwd = 2)



resraw <- StarletWT(TestF)
sigma <- MultiResSuppStarlet(resraw,2)
res <- HardThreshold(resraw,sigma)
plot(TestLambda,StarletRC(resraw),"l")
lines(TestLambda,StarletRC(res),col = 2,lwd = 2)
