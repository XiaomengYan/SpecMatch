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
TestLambda0 <- seq(from = Lambdamin,to = Lambdamax,length.out = 2^11)
TestF0 <- FluxFun(TestLambda)
lines(TestLambda,TestF0,"l",col = 2)

library(SpecMatch)
Rcpp::sourceCpp('~/Dropbox/project/Wang/SpecMatch/src/WTtools.cpp')


###################################################
resraw <- StarletWT(TestF0)

sigma <- MultiResSuppStarlet(TestF0,10)
res <- HardThreshold(resraw,3*sigma)
resraw[6:12,] = 0
res[6:12,] = 0
plot(TestLambda0,StarletRC(resraw),"l",lwd = 2)
lines(TestLambda0,StarletRC(res),col = 2,lwd = 2)
lines(TemInterpLambda,TemplateF_denoise,col=3)



resraw <- StarletWT(TestF)
sigma <- MultiResSuppStarlet(resraw,2)
res <- HardThreshold(resraw,sigma)
plot(TestLambda,StarletRC(resraw),"l")
lines(TestLambda,StarletRC(res),col = 2,lwd = 2)


#####################################################
library(SpecMatch)
CrossCorrelation(D6_1[10:3000,2:3],D6_1[,2:3])

