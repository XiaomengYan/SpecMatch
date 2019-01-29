SpectrumNormalization <- function(Specdf){
  maxfluxVec <- max(Specdf[,2])
  Specdf[,2] <- Specdf[,2]/maxfluxVec
  return(Specdf)
}

DeterStartEnd <- function(TestSignalWave, TemplateSignalWave){
  Startpoint <- max(min(TestSignalWave),min(TemplateSignalWave))
  Endpoint <- min(max(TestSignalWave),max(TemplateSignalWave))
  return(c(Startpoint,Endpoint))
}

######################
## Cross Correlation : Calculate the similarities between two spectrum.
## Input:
##        TestSpectra: dataframe(first column: wavelength, second column: flux)
##        TemplateSpectra: dataframe(first column: wavelength, second column: flux)


CrossCorrelation<- function(TestSpectra,TemplateSpectra){
  TemplateLambda <- TemplateSpectra[,1]
  TemplateF <- TemplateSpectra[,2]
  TemplateLen <- length(TemplateLambda)
  TemplateJ <- floor(log(TemplateLen,base = 2))
  commonLambda <- seq(from = 3000, to = 7000, length.out = 4000)
  TestLambda <- TestSpectra[,1]
  TestF <- TestSpectra[,2]
  TestLen <- length(TestLambda)
  TestJ <- floor(log(TestLen,base = 2))
  ## TemplateSpectra
  ###  Interpolation of TemplateSpectra
  TemplateFun <- splinefun(TemplateLambda,TemplateF,method = "fmm")
  TemInterpLambda <- seq(from = min(TemplateLambda),to = max(TemplateLambda),length.out = 2^TemplateJ)
  TemInterpF <- TemplateFun(TemInterpLambda)
  ### Denoise TemplateSpectra
  TemplateC0 <- StarletWT(TemInterpF)
  TemplateSigma <- MultiResSuppStarlet(TemInterpF,5)
  TemplateC <- HardThreshold(TemplateC0,3*TemplateSigma)
  TemplateC[floor((TemplateJ+1)/2):(TemplateJ+1),] = 0
  TemplateF_denoise <-StarletRC(TemplateC)## Remove noise and continumm

  TemplateDenoiseFun <- splinefun(TemInterpLambda,TemplateF_denoise,method = "fmm")
  TemInterpLambda_comm <- commonLambda[commonLambda>=min(TemInterpLambda)&commonLambda<=max(TemInterpLambda)]
  TemplateF_denoise_comm <- TemplateDenoiseFun(TemInterpLambda_comm)
  ## TestSpectra
  ### Interpolation of TestSpectra
  TestFun <- splinefun(TestLambda,TestF,method = "fmm")
  TestInterpLambda <- seq(from = min(TestLambda),to = max(TestLambda),length.out = 2^TestJ)
  TestInterpF <- TestFun(TestInterpLambda)
  ### Denoise TestSpectra
  TestC0 <- StarletWT(TestInterpF)
  TestSigma <- MultiResSuppStarlet(TestInterpF,5)
  TestC <- HardThreshold(TestC0,3*TestSigma)
  TestC[floor((TestJ+1)/2):TestJ+1,] <- 0
  TestF_denoise <- StarletRC(TestC)

  TestDenoiseFun <- splinefun(TestInterpLambda,TestF_denoise,method = "fmm")
  TestInterpLambda_comm <- commonLambda[commonLambda>=min(TestInterpLambda)&commonLambda<=max(TestInterpLambda)]
  TestF_denoise_comm <- TestDenoiseFun(TestInterpLambda_comm)

  ### Now Template (TemInterpLambda,Template_denoise)
  ###     Test (TestInterpLambda,TestF_denoise)
  ############################################################################
  ## Cross Correlation
  StartEnd <- DeterStartEnd(TestSignalWave = TestInterpLambda_comm,TemplateSignalWave = TemInterpLambda_comm)
  if(StartEnd[1]>=StartEnd[2])
    {print("There is no overlap between test signal and template signal")}
  else{
    N <- 10000
    Z <- seq(from = -0.01,to = 0.01,length.out = N)
    Crossfun <- rep(0,N)
    TestSeqWave <- log(TestInterpLambda_comm[which(TestInterpLambda_comm==StartEnd[1]):which(TestInterpLambda_comm==StartEnd[2])])
    TestSeqF <- TestF_denoise_comm[which(TestInterpLambda_comm==StartEnd[1]):which(TestInterpLambda_comm==StartEnd[2])]
    TestSeqFun <- splinefun(TestSeqWave,TestSeqF,method = "fmm")
    TemplateSeqWave <- log(TemInterpLambda_comm[which(TemInterpLambda_comm==StartEnd[1]):which(TemInterpLambda_comm==StartEnd[2])])
    TemplateSeqF <- TemplateF_denoise_comm[which(TemInterpLambda_comm==StartEnd[1]):which(TemInterpLambda_comm==StartEnd[2])]
    #TemplateSeqFun <- splinefun(TemplateSeqWave,TemplateSeqF,method = "fmm")
    for(i in 1:N){
      Crossfun[i] <- cor(TestSeqFun(TestSeqWave+log(1+Z[i])),TemplateSeqF)
    }
    tmp <- Z[which.max(Crossfun)]
  }##else
  return(c(max(Crossfun),tmp))
}














