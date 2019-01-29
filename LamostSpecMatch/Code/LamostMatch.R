library(R.utils)
library(FITSio)
# Match three fast speed white dwarf with the spectrum in LAMOST

SNe <- read.table("../LamostSpecMatch/Data/dbf4.txt",skip = 20)
D6_1 <- SNe[SNe[,1]== "D6-1",c(2,3)]
Z_1 <- 1293000/299792458
D6_1 <- D6_1[D6_1[,1]>=3000&D6_1[,1]<=7000,]
D6_1[,1] <- D6_1[,1]/(1+Z_1)
D6_1[,2] <- 1/max(D6_1[,2])*D6_1[,2]

D6_2 <- SNe[SNe[,1]== "D6-2",c(2,3)]
D6_2 <- D6_2[D6_2[,1]>=3000&D6_2[,1]<=7000,]
Z_2 <- 894000/299792458
D6_2[,1] <- D6_2[,1]/(1+Z_2)
D6_2[,2] <- 1/max(D6_2[,2])*D6_2[,2]

D6_3 <- SNe[SNe[,1]== "D6-3",c(2,3)]
D6_3 <- D6_3[D6_3[,1]>=3000&D6_3[,1]<=7000,]
Z_3 <- 1247000/299792458
D6_3[,1] <- D6_3[,1]/(1+Z_3)
D6_3[,2]<- 1/max(D6_3[,2])*D6_3[,2]

source('../LamostSpecMatch/SpecMatch/R/CrossCorrelation.R', echo=TRUE)
Rcpp::sourceCpp('../LamostSpecMatch/SpecMatch/src/CrossCorr.cpp')
Rcpp::sourceCpp('../LamostSpecMatch/SpecMatch/src/WTtools.cpp')


## Unzip the file and replace the original .fits.gz with .fits file

eachpath <- "/Volumes/XiaomengYan/Lamost_all/fits/EG024121N042539M01/"
allfile <- list.files(eachpath)
library(R.utils)
for (i in 1:length(allfile)) {
  gunzip(paste0(eachpath,allfile[i]))
}

## Read file and normalization
##
fitfile <- list.files(eachpath)
Output <- "/Volumes/XiaomengYan/Lamost_all/fits/summary.txt"
t1 <- Sys.time()
for (i in 1:length(fitfile)) {
  # read .fits file one by one
  print(i)
  testSpectfile <- readFITS(paste0(eachpath,fitfile[i]))
  Z <- as.numeric(testSpect$hdr[which(testSpect$hdr=="Z")+1]) # get redshift
  Specdf <- testSpectfile$imDat # get spectrum dataframe
  testSpectrum <- SpectrumNormalization(Specdf[,c(3,1)]) # normalization
  testSpectrum[,1] <- testSpectrum[,1]/(1+Z)# Deredshift
  # CrossCorrelation with three hypervelocity white dwarfs
  r1 <- CrossCorrelation(TestSpectra = testSpectrum,D6_1)
  r2 <- CrossCorrelation(TestSpectra = testSpectrum,D6_2)
  r3 <- CrossCorrelation(TestSpectra = testSpectrum,D6_3)
  write.table(x = data.frame(Name = fitfile[i], r1 = r1, r2 = r2, r3 = r3 ),file = Output,append = TRUE,quote = FALSE,row.names = FALSE,col.names = FALSE)
}
time <- Sys.time() - t1


# Plot
i = 1
print(i)
testSpectfile <- readFITS(paste0(eachpath,fitfile[i]))
Z <- as.numeric(testSpect$hdr[which(testSpect$hdr=="Z")+1]) # get redshift
Specdf <- testSpectfile$imDat # get spectrum dataframe
testSpectrum <- SpectrumNormalization(Specdf[,c(3,1)]) # normalization
testSpectrum[,1] <- testSpectrum[,1]/(1+Z)# Deredshift
plot(testSpectrum[,1],testSpectrum[,2],"l",ylab = "normalized flux",xlab = "wavelength")
lines(D6_1[,1],D6_1[,2],col=2)
lines(D6_2[,1],D6_2[,2],col=3)
lines(D6_3[,1],D6_3[,2],col=4)

# remove the redshift/ Fast Fourier transform cross-correlation:

CrossCorrelation(D6_1,D6_2)
CrossCorrelation(D6_3,D6_2)
CrossCorrelation(D6_3,D6_1)
