######################
## Cross Correlation : Calculate the similarities between two spectrum.
## Input:
##        TestSpectra: dataframe(first column: wavelength, second column: flux)
##        TemplateSpectra: dataframe(first column: wavelength, second column: flux)


CrossCorrelation <- function(TestSpectra,TemplateSpectra){
  TemplateLambda <- TemplateSpectra[,1]
  TemplateF <- TemplateSpectra[,2]
  TemplateLen <- length(TemplateLambda)
  TemplateJ <- floor(log(TemplateLen,base = 2))
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
  ## TestSpectra
  ### Interpolation of TestSpectra
  TestFun <- splinefun(TestLambda,TestF,method = "fmm")
  TestInterpLambda <- seq(from = min(TestLambda),to = max(TestLambda),length.out = 2^TestJ)
  TestInterpF <- TemplateFun(TestInterpLambda)
  ### Denoise TestSpectra
  TestC0 <- StarletWT(TestInterpF)
  TestSigma <- MultiResSuppStarlet(TestInterpF,5)
  TestC <- HardThreshold(TestC0,TestSigma)
  TestC[floor((TestJ+1)/2):TestJ+1,] <- 0
  TestF_denoise <- StarletRC(TestC)

  ### Now Template (TemInterpLambda,Template_denoise)
  ###     Test (TestInterpLambda,TestF_denoise)
  ############################################################################
  ## Cross Correlation
  StartEnd <- DeterStartEnd(TestSignalWave = TestLambda,TemplateSignalWave = TemInterpLambda)
  if(StartEnd[1]>=StartEnd[2]){print("There is no overlap between test signal and template signal")}
  else{
    LambdaSeq <- TemInterpLambda[StartEnd[1]:StartEnd[2]]
    TemplateSeq <- TemplateF_denoise[StartEnd[1]:StartEnd[2]]
    TestSeq <- rep(0,StartEnd[2]-StartEnd[1]+1)
    TestFun_denoise <- splinefun(TestInterpLambda,TestF_denoise,method = "fmm")
    for(i in 1:(StartEnd[2]-StartEnd[1]+1)){
      TestSeq[i] = TestFun_denoise(LambdaSeq[i])
    }##for
  }##else
  return(CrossCorr(TestSignal = TestSeq,TemplateSignal = TemplateSeq))
}


