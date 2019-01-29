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





