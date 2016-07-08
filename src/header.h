#include <RcppArmadillo.h>
#ifdef _OPENMP
# include <omp.h>
#endif
# include <iostream>
# include <cmath>

using namespace std;
using namespace arma;
using namespace Rcpp;

//Transformation functions (subfunc.cpp)
NumericMatrix CtoA(NumericMatrix& D, double by, int num_proc);
NumericMatrix CtoB2(NumericMatrix& D, double by, int num_proc);

//Subfunctions functions (subfunc.cpp)
void printArr(int arr[],int len);
Rcpp::NumericVector iarrToRvec(int arr[], int len);
Rcpp::NumericVector darrToRvec(double arr[], int len);
Rcpp::NumericMatrix armamatToRmat(arma::mat& matr);
void printRcppMat(Rcpp::NumericMatrix matr);

//Main clustering subfunctions (subfunc.cpp)
uvec find(int arr[], int len, int seek);
uvec find(arma::rowvec& arr, int seek);
NumericVector cq_agd(NumericMatrix& pts, double q, double eps);
void kmeansreg (arma::mat& point, arma::mat& cluster_center, arma::rowvec& cluster, arma::rowvec& cluster_energy,
              double p, double pw, double inn_tol, int inn_itmax);
void kmeansreg (arma::mat& point, arma::mat& cluster_center, arma::rowvec& cluster, arma::rowvec& cluster_energy,
              double p, double pw, int it_max, double inn_tol, int inn_itmax);

//Minimax calculation functions (mMcrit.cpp)
List mMcritPt(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_evalpts);
double mMcrit(arma::mat& point, NumericMatrix& Rcpp_evalpts);
double kmeansobj(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_evalpts, int p);

// //Objects (subfunc.cpp)
// class objective
// {
//     // Description:
//     // This is a wrapper function for computing C^q centers in minimax clustering
// public:
//     objective (const arma::mat& x, const double y)
//     {	D = x;
// 		p = y;
// 	}
//
//     double operator() (const column_vector& c) const
//     {        // return function value with parameters
//         int num_pts = D.n_rows;
//       	int dim_num = D.n_cols;
//       	double ret = 0; //value to return
//       	double tmp = 0; // dummy variable for eucl distance of each point
//       	for (int i=0;i<num_pts;i++){ // for each point in D
//       		tmp = 0;
//       		for (int j=0;j<dim_num;j++){
//       			tmp += pow(D(i,j)-c(j),2);
//       		}
//       		ret += pow(tmp,double(p)/2);
//       	}
//       //  cout << pow(ret,1/p) << endl;
//       	return (pow(ret,1/p));
//     }
//
// private:
//   arma::mat D;
// 	double p;
// };
//
// class objective_grad
// {
//     // Description:
//     // This is a wrapper function for computing the gradients of C^q centers in minimax clustering
// public:
//     objective_grad(const arma::mat& x, const double y)
//     {	D = x;
// 		p = y;
// 	}
//
//     const column_vector operator() (const column_vector& c) const
//     {        // return function value with parameters
//        int num_pts = D.n_rows;
//       	int dim_num = D.n_cols;
//       	column_vector ret = 0*dlib::ones_matrix<double>(dim_num,1); //value to return
//       	column_vector tmpvect(dim_num);
//         arma::rowvec tmprowvec(dim_num);
//
//       	double tmp = 0; // dummy variables
//       	double tmp2 = 0;
//
//       	for (int i=0;i<num_pts;i++){ // for each point in D
//       		tmp = 0;
//       		for (int j=0;j<dim_num;j++){
//       			tmp += pow(D(i,j)-c(j),2);
//       		}
//       //    cout << "I died here 0" << endl;
//           tmprowvec = D.row(i);
//           rowToCol(tmprowvec,tmpvect);
//       		ret += pow(tmp,double(p-2)/2)*(tmpvect-c);
//       //    cout << "I died here 1" << endl;
//       		tmp2 += pow(tmp,double(p)/2);
//       //    cout << "I died here 2" << endl;
//       	}
//       //  cout << "ret: " << endl << ret << endl;
//       //  cout << "tmp2: " << endl << tmp2 << endl;
//       //  cout << (-1*ret*pow(tmp2,1/p-1)) << endl;
//       	return (-1*ret*pow(tmp2,1/p-1));
//     }
//
// private:
//   arma::mat D;
// 	double p;
// };
