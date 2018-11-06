#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
using namespace Rcpp;
using namespace arma;


//[[Rcpp::export]]
arma::mat StarletWT(arma::vec X){
// Starlet Second Generation WT
  int N = X.size(); // Get the length of the input signal X
  int J = floor(log2(N)); // Get the total scale of the wavelet transform
  arma::mat C(J+1,N,fill::zeros),W(J,N,fill::zeros);
  arma::mat V(J+1,N,fill::zeros);
  for(int i = 0 ;i <N;i++){
    C.at(0,i) = X.at(i);
  }
  arma::vec h_1(5,fill::ones),h_1D(5,fill::ones);
  h_1[1] = 4;
  h_1[2] = 6;
  h_1[3] = 4;
  h_1D = 1.0/16 * h_1;
  double e = 0.0;
  for(int j = 0;j< J;j++){
    // Level
    for(int l = 0; l< N;l++){
      //position
      double sum = 0.0;
      for(int k = -2;k<3;k++){
        int lpk = l+pow(2,j)*k;
        while(lpk<0||lpk>(N-1)){
          if(lpk<0){
            lpk = -1*lpk;
          }
          if(lpk > (N-1)){
            lpk = 2*(N-1)-lpk;
          }
        }
        e = C.at(j,lpk)*h_1D[k+2];
        sum = sum + e;
      }
      C.at(j+1,l) = sum;
    }//position
// Second Transform towards c_j+1
    for(int l = 0; l< N;l++){
      double sum = 0.0;
      for(int k = -2;k<3;k++){
        int lpk = l+pow(2,j)*k;
        while(lpk<0||lpk>(N-1)){
          if(lpk<0){
            lpk = -1*lpk;
          }
          if(lpk>(N-1)){
            lpk = 2*(N-1)-lpk;
          }
        }
        e = C.at(j+1,lpk)*h_1D[k+2];
        sum = sum + e;
      }
      V.at(j+1,l) = sum;
    }//position
    W.row(j) = C.row(j)-V.row(j+1);
  }//level
  for(int j = 0; j< J+1;j++){
    if(j<J){
      C.row(j) = W.row(j);
    }else{
      C.row(j) = C.row(j);
    }
  }
  return(C);
}

//[[Rcpp::export]]
arma::mat StarletRC(arma::mat X){
  // Reconstruction Algorithm for the Starlet Second Generation WT
  int N = X.n_cols;
  int J = X.n_rows;
  arma::vec h_1(5,fill::ones),h_1D(5,fill::ones);
  h_1[1] = 4;
  h_1[2] = 6;
  h_1[3] = 4;
  h_1D = 1.0/16 * h_1;
  arma::mat V(1,N,fill::zeros);
  arma::mat C(J,N);//J_1 the Number of J + 1
  C.row(J-1) = X.row(J-1);
  double e = 0.0;
  for(int j = J-1;j > 0 ;j--){
    for(int l = 0;l< N;l++){
      double sum = 0.0;
      for(int k = -2;k<3;k++){
        double lpk = l+pow(2,j-1)*k;
        while(lpk<0||lpk>(N-1)){
          if(lpk<0){
            lpk = -1*lpk;
          }
          if(lpk>(N-1)){
            lpk = 2*(N-1)-lpk;
          }
        }//while
        e = C.at(j,lpk)*h_1D(k+2);
        sum = sum + e;
    }
      V.at(0,l) = sum;
  }
    C.row(j-1) = V + X.row(j-1);
}
  return(C.row(0));
}

//[[Rcpp::export]]
Rcpp::List PMT(arma::vec X){
  //pyramidal median transform
  int N = X.size();
  int J = floor(log2(N));
  arma::mat C(J+1,N,fill::zeros);
  arma::mat C_t(J,N,fill::zeros);
  arma::mat C_med(J,N,fill::zeros);
  arma::mat W(J,N,fill::zeros);
  arma::vec window(7,fill::zeros);
  arma::vec V(J+1,fill::zeros);
  V.at(0) = N;//Vector to record the number
               //of element in each level. The 0 element is the original signal lengt
  for(int i = 0 ;i <N;i++){
    C.at(0,i) = X.at(i);
  }
  for(int j=0;j<J;j++){
    for(int l=0;l<N;l++){
      for(int i = l-3;i<l+4;i++){
        int h = i;
        while(h<0||h>N-1){
          if(h<0){h = -1* h;}
          if(h>(N-1)){h = 2*(N-1)-h;}
        }
        window.at(i+3-l) = C.at(j,h);
      }//i
      C_med(j,l) = median(window);
    }//l
    // Extract the odd term of C_med(j+1,:) to get C(j+1,:)
    V.at(j+1) = floor(N/2.0)+N%2;
    for(int k = 0;k< V.at(j+1);k++){
      C(j+1,k) = C_med.at(j,2*k);
    }//k
    //Interpolation of C_{j+1} to the size of C_{j}
    int V_j = V.at(j+1);
    int t = 2*V_j - N%2;
    // If N is odd
    if(N%2==1){
      for(int m = 0;m<t;m++){
        for(int k = 0;k<V_j;k++){
          if(m==2*k){C_t.at(j+1,m) = C(j+1,k);}
          if(m==(2*k+1)&&m<N){C_t.at(j+1,m) = (C.at(j+1,k)+C.at(j+1,k+1))/2.0;}
        }//k
      }//m
    }//if
    // If N is even
    if(N%2==0){
      for(int m = 0;m<(t-1);m++){
        for(int k = 0;k<V_j;k++){
          if(m==2*k){C_t.at(j,m) = C(j+1,k);}
          if(m==(2*k+1)&&m<N){C_t.at(j,m) = (C.at(j+1,k)+C.at(j+1,k+1))/2.0;}
        }//k
      }//m
      //The last term is equal to the term before the last term
      C_t.at(j,t-1) = C_t.at(j,t-2);
    }//if
    N = V_j;
    W.row(j) = C.row(j) - C_t.row(j);
  }//j
  for(int j = 0;j<J;j++){
    C.row(j) = W.row(j);
  }
  return Rcpp::List::create(Rcpp::Named("C") = C,Rcpp::Named("NumberIndex") = V);
}

//[[Rcpp::export]]
arma::mat PMTRC(arma::mat C,arma::vec V){
  //reconstruction of prymdial median transform
  //We use linear intorpolation
  //Input: coefficient matrix X = {w1,...wJ,cJ}; V = Number of element in each scale
  int N = C.n_cols;
  int J = C.n_rows; //(log2(N)+1)
  arma::mat C_t(J,N,fill::zeros);
  for(int j = J-1;j>0;j--){
    int V_j_1 = V.at(j-1);
    N = V.at(j)*2 - V_j_1%2;
    //Interpolation of C_{j+1} to the size of C_{j}
    // If N is odd
    if(V_j_1%2==1){
      for(int m = 0;m<N;m++){
        for(int k = 0;k<V.at(j);k++){
          if(m==2*k){C_t.at(j,m) = C(j,k);}
          if(m==(2*k+1)&&m<N){C_t.at(j,m) = (C.at(j,k)+C.at(j,k+1))/2.0;}
        }//k
      }//m
    }//if
    // If N is even
    if(V_j_1%2==0){
      for(int m = 0;m<N;m++){
        for(int k = 0;k<V.at(j);k++){
          if(m==2*k){C_t.at(j,m) = C(j,k);}
          if(m==(2*k+1)&&m<N){C_t.at(j,m) = (C.at(j,k)+C.at(j,k+1))/2.0;}
        }//k
      }//m
      //The last term is equal to the term before the last term
      C_t.at(j,N-1) = C_t.at(j,N-2);
    }//if
    C.row(j-1) = C.row(j-1) + C_t.row(j);
  }//j
  return C.row(0);
}



//////////////////////////////////////////////////////////////////////////////
///// The Multiresolution support ///////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
// //[[Rcpp::export]]
// arma::mat MultiResSuppPMT(arma::vec X,int SimNo){
//   //Get the multiresolution support of wavelet coefficient C based on simulation
//   //Input: The origianl signal
//   //Output: Simulation size
//   Rcpp::List C(SimNo);
//   int N = X.size();
//   int J = floor(log2(N));
//   arma::mat C_n;
//   arma::mat M(J+1,N,fill::zeros);
//   arma::vec RandDraw(N,fill::zeros);
//   Rcpp::List PMT_res = PMT(X);
//   arma::mat C_signal = PMT_res.at(0);
//   arma::vec V_signal = PMT_res.at(1);
//   for(int n=0;n<SimNo;n++){
//     RandDraw = randn(N);
//     arma::mat C_n = PMT(RandDraw).at(0);
//     C.at(n) = C_n;
//     //Rcpp::Rcout<<C_n<<std::endl;
//   }
//   for(int j =0;j<J+1;j++){
//     //position
//     int V_j = V_signal.at(j);
//     //Rcpp::Rcout<< V_j <<std::endl;
//     for(int n =0;n< V_j;n++){
//       //scale
//       arma::vec tmp(SimNo,fill::zeros);
//       for(int i = 0;i<SimNo;i++){
//         arma::mat C_i = C.at(i);
//         tmp.at(i) = C_i.at(j,n);
//       }//i
//       //double absCjn = fabs(C_signal.at(j,n));
//       //if( absCjn > 3*stddev(tmp)){M.at(j,n) = 1;}
//       M.at(j,n) = stddev(tmp);
//       //Rcpp::Rcout<<tmp<<std::endl;
//     }//j
//   }//n
//   return M;
// }



//[[Rcpp::export]]
double MAD(arma::rowvec X){
  double Med = median(X);
  arma::vec MADres(X.size(),fill::zeros);
  for(int j = 0;j<X.size();j++){
    MADres.at(j) = fabs(X.at(j) - Med);
  }
  double res = median(MADres);
  return res;
}


//[[Rcpp::export]]
arma::vec MultiResSuppStarlet(arma::vec X,int NoSimu){
    //Get the multiresolution support of wavelet coefficient C based on simulation
    //Input: The origianl signal
    //Output: noise standard deviation of each level
  int N = X.size();
  int J = floor(log2(N));
  arma::mat sigmaEMat(NoSimu,J,fill::zeros);// The noise standard deviation generated by simulation and transform
  arma::mat ResSignal = StarletWT(X);
  arma::vec sigmaJ(J,fill::zeros);// The noise standard deviation of coefficients at each level j
  double sigmaS;//The noise standard deviation of original signal MED(|w_1|)/0.6745
  for(int n = 0;n<NoSimu;n++){
    arma::vec RandDraw = randn(N);
    arma::mat ResStarlet = StarletWT(RandDraw);
    for(int j=0; j< J; j++){
      sigmaEMat.at(n,j) = stddev(ResStarlet.row(j));
    }//j
  }//n
  arma::rowvec w1 = ResSignal.row(0);
  double sigmaE1 = mean(sigmaEMat.col(0));
  sigmaS = MAD(w1)/0.6745;
  for(int j = 0;j<J;j++){
    sigmaJ.at(j) = sigmaS/sigmaE1*mean(sigmaEMat.col(j));
  }
  return sigmaJ;
}





// //[[Rcpp::export]]
// arma::mat MultiResSuppStarlet(arma::vec X,int SimNo){
//   //Get the multiresolution support of wavelet coefficient C based on simulation
//   //Input: The origianl signal
//   //Output: Simulation size
//   Rcpp::List C(SimNo);
//   int N = X.size();
//   int J = floor(log2(N));
//   arma::mat C_n;
//   arma::mat M(J+1,N,fill::zeros);
//   arma::vec RandDraw(N,fill::zeros);
//   arma::mat C_signal = StarletWT(X);
//   for(int n=0;n<SimNo;n++){
//     RandDraw = randn(N);
//     arma::mat C_n= StarletWT(RandDraw);
//     C.at(n) = C_n;
//   }
//   for(int n =0;n< N;n++){
//     //position
//     for(int j =0;j<J+1;j++){
//       //scale
//       arma::vec tmp(SimNo,fill::zeros);
//       for(int i = 0;i<SimNo;i++){
//         arma::mat C_i = C.at(i);
//         tmp.at(i) = C_i.at(j,n);
//       }//i
//       //double absCjn = fabs(C_signal.at(j,n));
//       M.at(j,n) = stddev(tmp);
//      // if( absCjn > 3*stddev(tmp)){M.at(j,n) = 1;}
//     }//j
//   }//n
//   return M;
// }


///////////////////////////////////////////////////////////////////////
///////// Threshold////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
arma::mat HardThreshold(arma::mat C, arma::rowvec M){
  int N = C.n_cols;
  int J = C.n_rows;
  for(int n = 0;n<N;n++){
    for(int j = 0;j<J;j++){
      if(fabs(C.at(j,n)) < 3*M.at(j)){C.at(j,n) = 0;}
    }
  }
  return C;
}

//[[Rcpp::export]]
int sign(double x){
  int i=0;
  if(x>0){i = 1;}
  return i;
}

//[[Rcpp::export]]
arma::mat SoftThreshold(arma::mat C, arma::vec M){
  int N = C.n_cols;
  int J = C.n_rows;
  for(int n = 0;n<N;n++){
    for(int j = 0;j<J;j++){
      if(fabs(C.at(j,n)) < 3*M.at(j)){C.at(j,n) = 0;}
      else{C.at(j,n) = sign(C.at(j,n))*fabs(C.at(j,n))-M.at(j);}
    }
  }
  return C;
}
















