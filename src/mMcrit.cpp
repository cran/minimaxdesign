# include <omp.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
#include <RcppArmadillo.h>
#include "header.h"

using namespace std;
using namespace Rcpp;

//--------------------------------------------------------------
//Computes the minimax criterion and corresponding farthest point for R plotting
//--------------------------------------------------------------
// [[Rcpp::export]]
List mMcritPt(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_evalpts){
  // Description of Inputs:
  // Rcpp_point     - Current design
  // Rcpp_evalpts   - Evaluation points for minimax criterion

  int dim_num = Rcpp_point.ncol(); //dimension of data points
  int point_num = Rcpp_point.nrow(); //number of data points
  int eval_num = Rcpp_evalpts.nrow(); // number of evaluation points for mM
  double ret = 0; //to track mM
  NumericVector far_pt(dim_num); //to record farthest point in design region
  NumericVector des_pt(dim_num); //to record farthest point in design region
//  NumericVector far_dpt(dim_num); // to record the design point closest to far_pt
//  NumericVector tmp_dpt(dim_num); //temporary container
  double dst;
  double tmp; //for temporary storage
  int idx;

  for (int i=0;i<eval_num;i++){ // for each evaluation point
    dst = dim_num; //i.e., arb. large
    NumericVector tmpvec = Rcpp_evalpts(i,_);
    for (int j=0;j<point_num;j++){ // for each design point
       //find the closest distance to a design point
//       tmp = sum(pow(Rcpp_evalpts(i,_)-Rcpp_point(j,_),2));
       tmp = sum(pow(tmpvec-Rcpp_point(j,_),2));
       if (dst > tmp){
         dst = tmp;
         idx = j;
       }
    }
    //find the maximum of such distances
    if (ret < dst){
      ret = dst;
      far_pt = Rcpp_evalpts(i,_);
      des_pt = Rcpp_point(idx,_);
    }
  }
  return List::create(Named("mM_dist")=sqrt(ret), Named("far_pt")=far_pt, Named("des_pt")=des_pt);
}

//--------------------------------------------------------------
//Returns the minimum distance for all design points
//--------------------------------------------------------------
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
NumericVector mMcrit_allpts(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_evalpts){
  // Description of Inputs:
  // Rcpp_point     - Current design
  // Rcpp_evalpts   - Evaluation points for minimax criterion

  int dim_num = Rcpp_point.ncol(); //dimension of data points
  int point_num = Rcpp_point.nrow(); //number of data points
  int eval_num = Rcpp_evalpts.nrow(); // number of evaluation points for mM
  double dst;
  NumericVector ret(point_num);
  for (int i=0; i<point_num; i++){
    ret(i) = 0.0;
  }

  #pragma omp parallel for
  for (int i=0;i<eval_num;i++){ // for each evaluation point
    double dst = DBL_MAX; //i.e., arb. large
    double tmp;
    int pt_idx;
    for (int j=0;j<point_num;j++){ // for each design point
       double tmp = sum(pow(Rcpp_evalpts(i,_)-Rcpp_point(j,_),2));
       if (tmp < dst){
         // Update distance and idx
         dst = tmp;
         pt_idx = j;
       }
    }
    //Update distance of closest point
    if (ret(pt_idx)<dst){
      ret(pt_idx)=dst;
    }
  }

  for (int i=0; i<point_num; i++){
    ret(i) = pow( ret(i), 0.5);
  }
  return (ret);
}

//--------------------------------------------------------------
//Computes minimax criterion for PSO (internal use)
//--------------------------------------------------------------
double mMcrit(arma::mat& point, NumericMatrix& Rcpp_evalpts){
  // Description of Inputs:
  // Rcpp_point     - Current design
  // Rcpp_evalpts   - Evaluation points for minimax criterion
  int dim_num = point.n_cols; //dimension of data points
  int point_num = point.n_rows; //number of data points
  int eval_num = Rcpp_evalpts.nrow(); // number of evaluation points for mM
  double ret = 0.0; //to track mM
  double dst;

  for (int i=0;i<eval_num;i++){ // for each evaluation point
    dst = DBL_MAX; //i.e., arb. large
    for (int j=0;j<point_num;j++){ // for each design point
       double tmp = 0.0;
       //find the closest distance to a design point
       for (int k=0; k<dim_num; k++){
          tmp += pow(Rcpp_evalpts(i,k) - point(j,k), 2.0);
       }
       if (tmp < dst){
         dst = tmp;
       }
    }
    //find the maximum of such distances
    if (dst > ret){
      ret = dst;
    }
  }
  return (sqrt(ret));
}

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
double mMcrit_proj(NumericMatrix& Rcpp_pts, NumericMatrix& Rcpp_evalpts, NumericMatrix& indices){
  //mMcrit_proj for use in subfunctions
  int dim_num = Rcpp_pts.ncol(); //dimension of data points
  int point_num = Rcpp_pts.nrow(); //number of data points
  int eval_num = Rcpp_evalpts.nrow(); // number of evaluation points for mM
  int ord = Rcpp_evalpts.ncol(); //Order
  int num_idx = indices.nrow(); //Number of indices
  arma::mat point(Rcpp_pts.begin(),point_num,dim_num,false);
  vector<double> dst(num_idx);

  #pragma omp parallel for
  for (int i=0; i<num_idx; i++){
    uvec idx(ord);
    for (int j=0; j<ord; j++){
      idx(j) = indices(i,j)-1;
    }
    arma::mat proj_pts = point.cols(idx);
    dst[i] = mMcrit(proj_pts,Rcpp_evalpts);
  }

  double ret = 0.0;
  for (int i=0; i<num_idx; i++){
    // Rcout << dst[i] << " " << endl;
    ret = max(ret, dst[i]);
  }

  return(ret);
}

//double mMcrit_proj(arma::mat& point, NumericMatrix& Rcpp_evalpts, arma::mat& theta, double ord){
//  //mMcrit_proj for use in subfunctions
//  //ord is q in formula
//  int dim_num = point.n_cols; //dimension of data points
//  int point_num = point.n_rows; //number of data points
//  int eval_num = Rcpp_evalpts.nrow(); // number of evaluation points for mM
//  int theta_num = theta.n_rows; //number of theta points
//
//  double tmp = 0; //to return
//  double tmp2 = 0;
//  double tmp3 = 0;
//  double tmp4 = 0;
//
//  for (int l=0;l<theta_num;l++){
//    tmp2 = 0;
//    for (int i=0;i<eval_num;i++){
//      tmp3 = 0;
//      for (int j=0;j<point_num;j++){
//        tmp4 = 0;
//        for (int m=0;m<dim_num;m++){
////          bool test = (theta(l,m)==0.0);
////          cout << test << endl;
//          if (theta(l,m)>0.0){
//            tmp4 += theta(l,m)*pow(point(j,m)-Rcpp_evalpts(i,m),2);
//          }
//        }
////        tmp3 += 1/tmp4; //tmp4 is theta norm of difference
//        tmp3 += 1/pow(tmp4,ord); //tmp4 is theta norm of difference
//      }
//      tmp2 += 1/(tmp3/point_num);
//    }
//    if (tmp < pow(tmp2/eval_num,1/(2*ord))){
//      tmp = pow(tmp2/eval_num,1/(2*ord));
//    }
//
//  }
//  return(tmp);
//}

//--------------------------------------------------------------
//Computes the finite-sample clustering approximation of minimax criterion
//--------------------------------------------------------------
// [[Rcpp::export]]
double kmeansobj(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_evalpts, int p){
  // Description of Inputs:
  // Rcpp_point     - Current design
  // Rcpp_evalpts   - Evaluation points for minimax criterion
  // p              - q in paper; approximation coefficient for minimax criterion
    int point_num = Rcpp_evalpts.nrow();
    int cluster_num = Rcpp_point.nrow();
    int dim_num = Rcpp_point.ncol();
    NumericVector cluster(point_num);

    double ret = 0;
    double point_energy = 0;
    double point_energy_min;
    double tmp_energy;

    for ( int j = 0; j < point_num; j++ ) // for each data point
    {
      //Assign clusters
      point_energy_min = DBL_MAX;

      for ( int k = 0; k < cluster_num; k++ ){ // for each cluster
        point_energy=0;
        for (int m=0;m<dim_num;m++){
          point_energy += pow(Rcpp_evalpts(j,m)-Rcpp_point(k,m),2);
        }
        if ( point_energy < point_energy_min ){
          point_energy_min = point_energy;
          cluster(j) = k;
        }
      }
      //Compute objective
      tmp_energy = 0;
      for (int m=0;m<dim_num;m++){
        tmp_energy += pow(Rcpp_evalpts(j,m)-Rcpp_point(cluster(j),m),2);
      }
      ret += pow(tmp_energy,p/2);

    }
    return (ret);
}
