#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

//--------------------------------------------------------------
// Gamma function
//--------------------------------------------------------------
double Gamma(double xx)  {
  double cof[7],stp,x,tmp,ser;
  int j;
  cof[1]=76.18009173;
  cof[2]=-86.50532033;
  cof[3]=24.01409822;
  cof[4]=-1.231739516;
  cof[5]=0.120858003e-2;
  cof[6]=-0.536382e-5;
  stp=2.50662827465;

  x=xx-1.0;
  tmp=x+5.5;
  tmp=(x+0.5)*log(tmp)-tmp;
  ser=1.0;
  for (j=1; j<7; j++) {
    x=x+1.0;
    ser=ser+cof[j]/x;
  }
  return (exp(tmp+log(stp*ser)));
}

//--------------------------------------------------------------
//Beta function
//--------------------------------------------------------------
double Beta(double a, double b){
  return (Gamma(a)*Gamma(b)/Gamma(a+b));
}

//--------------------------------------------------------------
//F^(-1) for j = 2, ..., s
//--------------------------------------------------------------
double invF(int s, int j, double target, double by){
  // Description of Input
  // s      - p in paper; dimension of design space
  // j      - Desired dimension to evaluate
  // target - Desired output of F^(-1)
  // by     - Approximation step-size
  double constant = M_PI/Beta(0.5,(s-j+1.0)/2);
  double newtarg = target/constant;
  double runsum = 0.0;
  double index = 0.0;
  while(runsum < newtarg){
    runsum += pow(sin(M_PI*index),s-j)*by;
    index += by;
  }
  return(index);
}

//--------------------------------------------------------------
//Rosenblatt transformation from hypercube [0,1]^p to ball B^p
//--------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix CtoBp(NumericMatrix& D, double by, int num_proc){
  // Description of Input
  // D        - Low-discrepancy points on [0,1]^p
  // by       - Approximation step-size
  // num_proc - Number of processors to use in parallel computing
  int s = D.ncol();
  NumericVector b(s);
  NumericVector svect(s); //cumulative
  NumericVector cvect(s);
  // #ifdef _OPENMP
  //   omp_set_num_threads(num_proc);
  // #endif

  NumericMatrix ret(D.nrow(),D.ncol());

  // #pragma omp parallel for
  for (int i=0;i<D.nrow();i++){
    //compute new b vector
    // Rcout << "i: " << i << std::endl;
    for (int j=0;j<D.ncol();j++){
      // Rcout << "j: " << j << std::endl;
      if(j==0){
        // Rcout << "Here 0" << std::endl;
        b(0)=pow(D(i,0),1/double(s));
        svect(0)=1.0; //doesn't matter
        cvect(0)=1.0; //doesn't matter
      }
      else if(j==(D.ncol()-1)){
        // Rcout << "Here 1" << std::endl;
        b(j)=invF(s,j+1,D(i,j),by);
        svect(j)=svect(j-1)*sin(2*M_PI*b(j));
        cvect(j)=cos(2*M_PI*b(j));
      }
      else{
        // Rcout << "Here 2" << std::endl;
        b(j)=invF(s,j+1,D(i,j),by);
        svect(j)=svect(j-1)*sin(M_PI*b(j));
        cvect(j)=cos(M_PI*b(j));
      }
      }

    //...then compute the transform x
    for (int j=0;j<D.ncol();j++){
      // Rcout << "Here 3" << std::endl;
      if (j==(D.ncol()-1)){
        ret(i,j)=b(0)*svect(D.ncol()-1);
      }
      else{
        ret(i,j)=b(0)*svect(j)*cvect(j+1);
      }
    }
  }
  return (ret);
}
