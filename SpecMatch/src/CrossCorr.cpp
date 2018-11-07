#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double CrossCorr(arma::vec TestSignal, arma::vec TemplateSignal){
  //After determine the StartEnd point and intepolate towards the TestSignal. Now TestSignal and TemplateSignal
  // are aligned. Apply the general Cross Correlation
  int N = TestSignal.size(); // The size of TestSignal and TemplateSignal(cut-off) are the same.
  double Testbar = mean(TestSignal);
  double Templatebar = mean(TemplateSignal);
  double Testsigma = stddev(TestSignal);
  double Templatesigma = stddev(TemplateSignal);
  double tmp = 0;
  for(int n=0;n<N;n++){
    // Rcpp::Rcout<<TestSignal.at(n)-Testbar<<std::endl;
    // Rcpp::Rcout<<TemplateSignal.at(n)-Templatebar<<std::endl;
    tmp = tmp + (TestSignal.at(n)-Testbar)*(TemplateSignal.at(n)-Templatebar);
  }//for
  return tmp/((N-1)*Testsigma*Templatesigma);
}



//[[Rcpp::export]]
double test(arma::vec hi){
  return mean(hi);
}

//[[Rcpp::export]]
arma::vec DeterStartEnd(arma::vec TestSignalWave, arma::vec TemplateSignalWave){
  // To determine where to start the cross-correlation respect to the Template Signal
  double Templatemin = TemplateSignalWave.min();
  double Templatemax = TemplateSignalWave.max();
  int N_template = TemplateSignalWave.size();
  double Testmin = TestSignalWave.min();
  double Testmax = TestSignalWave.max();
  arma::vec StartEnd(2,fill::zeros);
  //Determine the start point
  if(Testmin > Templatemin && Testmin<Templatemax){
    for(int n=0;n<N_template;n++){
      if(Testmin >TemplateSignalWave.at(n)&&Testmin <= TemplateSignalWave.at(n+1)){StartEnd.at(0) = n+2;break;}//if
    }//for
  }else if(Testmin <= Templatemin){
    StartEnd.at(0) = 1;
  }else{
    StartEnd.at(0) = N_template+1;
  }
  //Determine the end point
  if(Testmax>Templatemin && Testmax<Templatemax){
    for(int n=0;n<N_template;n++){
      if(Testmax>TemplateSignalWave.at(n)&&Testmax <= TemplateSignalWave.at(n+1)){StartEnd.at(1) = n+2;break;}//if
    }//for
  }else if(Testmax <= Templatemin){
    StartEnd.at(1) = 0;
  }else{
    StartEnd.at(1) = N_template;
  }
  return StartEnd;
}
