#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstring>
#include <time.h>
#include <omp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

#include "header.h"

//--------------------------------------------------------------
//Enforcing box constraints of [0,1]^p
//--------------------------------------------------------------
void enfBox(arma::mat& des, double lb, double ub){
  // Description of Inputs:
  // des    - Current design
  // lb     - Lower bound to enforce
  // ub     - Upper bound to enforce
  int cluster_num = des.n_rows;
  int dim_num = des.n_cols;
  for (int i=0;i<cluster_num;i++){
    for (int j=0;j<dim_num;j++){
      	des(i,j) = min(max(des(i,j),lb),ub);
	  }
  }
}

//--------------------------------------------------------------
//Jitters each design point by no more than tol; used for minimax post-processing
//--------------------------------------------------------------
void jitter(arma::mat& des, arma::cube& des_part, double tol, double lb, double ub){
  // Description of Inputs:
  // des      - Current design
  // des_part - Container to return jittered particles
  // tol      - Maximum distance for jittering
  // lb       - Lower bound to enforce
  // ub       - Upper bound to enforce
  int part_num = des_part.n_slices;//number of particles
  int cluster_num = des.n_rows;
  int dim_num = des.n_cols;
  for (int i=0;i<cluster_num;i++){
    for (int j=0;j<dim_num;j++){
      for (int k=0;k<(part_num-1);k++){
  			des_part(i,j,k) = min(max(des(i,j) + 2*tol*randu() - tol,lb),ub); //don't let it jitter outside [0,1]^p
		  }
	  }
  }

  //Add the original design
  for (int i=0; i<cluster_num; i++){
    for (int j=0; j<dim_num; j++){
      des_part(i,j,part_num-1) = des(i,j);
    }
  }

}

//--------------------------------------------------------------
//Post-processing function using the minimax criterion
//--------------------------------------------------------------
// [[Rcpp::plugins(openmp)]]
double mMcritPSO(arma::mat& des, NumericMatrix& Rcpp_evalpts, int part_num,
    int it_max, int it_lim,
    double it_tol, double tol, double lb, double ub,
    double w, double c1, double c2){
  // Description of Inputs:
  // des          - Global-best design from minimax clustering
  // Rcpp_evalpts - Points for minimax post-processing
  // part_num     - Number of PSO particles for minimax clustering
  // it_max       - Maximum number of iterations for minimax clustering
  // tol          - Jitter tolerance in minimax post-processing
  // lb           - Lower bound of design
  // ub           - Upper bound of design
  part_num ++; //one particle for initial particle from kmeans
  int cluster_num = des.n_rows;
  int dim_num = des.n_cols;
  arma::cube des_vel = zeros(cluster_num,dim_num,part_num);//current velocity
  arma::cube des_part = zeros(cluster_num,dim_num,part_num);//particles
  arma::cube des_lbes = zeros(cluster_num,dim_num,part_num);//lbes positions
    //global best position will be saved in des
  arma::mat Dtmp(cluster_num,dim_num); //temporary container for computation

  arma::vec cur_obj = DBL_MAX*ones(part_num); //current objectives
  arma::vec lbes_obj = DBL_MAX*ones(part_num); //local best objectives
  double gbes_obj = DBL_MAX;
  double prev_gbes_obj = DBL_MAX;

  //Jitter points to make new particles
  jitter(des,des_part,tol,lb,ub);

  //Do PSO with mMcrit (true minimax) function
  clock_t t = clock(); //timing
  int it_num = 0;
  int it_same = 0;

  while ((it_num<it_max)&&(it_same<it_lim))
  {

    Rcout << "Post-processing iteration " << it_num << " ... " << endl;

    //Update minimax
    #pragma omp parallel for
    for (int i=0;i<part_num;++i){
      cur_obj(i) = mMcrit(des_part.slice(i), Rcpp_evalpts);
//      cout << "cur_obj(i): " << cur_obj(i) << endl;
    };

    //Update lbes and gbes
    bool ch_flg = false;
    for (int i=0;i<part_num;i++){
      if (cur_obj(i)<lbes_obj(i)){
        des_lbes.slice(i) = des_part.slice(i);
        lbes_obj(i) = cur_obj(i);
//        cout << "lbes_obj changed to:" << endl;
//        cout << lbes_obj.t << endl;
      }
      if (cur_obj(i)<gbes_obj){
        //Update gbes
        des = des_part.slice(i);
        prev_gbes_obj = gbes_obj;
        gbes_obj = cur_obj(i);
        Rcout <<  "gbes_obj changed to:" << endl;
        Rcout << gbes_obj << endl;

        //update flag if threshold exceeded
        if (prev_gbes_obj - gbes_obj > it_tol){
          ch_flg = true;
        }

      }

    }
    //Update or reset it_same
    if (ch_flg){
      it_same = 0;
    }
    else{
      it_same++;
    }

//    //Update velocities
//    double w = 0.72;//inertia constant (from Merwe and Engelbrecht, 2003)
//    double c1 = 1.49;//acceleration constants
//    double c2 = 1.49;

    for (int i=0; i<part_num; i++){
      Dtmp = des_part.slice(i);
      //update velocity
      des_vel.slice(i) = w*des_vel.slice(i) + c1*(randu(cluster_num,dim_num)%(des_lbes.slice(i)-Dtmp)) + c2*(randu(cluster_num,dim_num)%(des-Dtmp));
      //update particle
      des_part.slice(i) += des_vel.slice(i);
      enfBox(des_part.slice(i),lb,ub);
    }

//    cout << "Time :" << ((float)clock()-t)/CLOCKS_PER_SEC << endl;
//    cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

    it_num ++;
  }

  return(gbes_obj);
}

//--------------------------------------------------------------
//Test function: mMc (without PSO)
//--------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix kmeansreg(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_cluster_center,
            double p, double pw, int it_max, double inn_tol, int num_proc)
{
  // Description of Inputs:
  // Rcpp_point          - Clustering data for minimax clustering
  // Rcpp_cluster_center - Initial cluster centers
  // p                   - q in paper; approximation coefficient for minimax criterion
  // it_max              - Maximum number of iterations for minimax clustering
  // num_proc            - Number of processors to use in parallel computing

  int dim_num = Rcpp_point.ncol(); //dimension of cluster data points
  int point_num = Rcpp_point.nrow(); //number of cluster data points
	int cluster_num = Rcpp_cluster_center.nrow(); //number of desired clusters
  int inn_itmax = 1e4; //maximum number of agd iteration

  //cluster energy matrix
  arma::rowvec cluster_energy(cluster_num);
  cluster_energy.zeros(); //fill with zeros
  //keeps track of cluster assignment (each row for a particle)
  arma::rowvec cluster(point_num);
  cluster.ones(); //fill with -1
  cluster = -1*cluster;

  //clustering points
  arma::mat point(Rcpp_point.begin(),point_num,dim_num,false);
  //current cluster centers (or design points)
  arma::mat cluster_center(Rcpp_cluster_center.begin(),cluster_num,dim_num,false);

//  omp_set_num_threads(num_proc);

  //PSO iterations
  Rcout << "-------------------------------------------------" << endl;
  Rcout << "Minimax clustering ... " << endl;
  Rcout << "-------------------------------------------------" << endl;

  kmeansreg(point,cluster_center,cluster,cluster_energy,p, pw,it_max,inn_tol,inn_itmax);

  //wrap to Rcpp classes for output
  NumericMatrix ret_cluster_center(cluster_num,dim_num);

  ret_cluster_center = armamatToRmat(cluster_center);

  return (ret_cluster_center);

  }

//--------------------------------------------------------------
//Main function for mMc-PSO
//--------------------------------------------------------------
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::export]]
List kmeanspso(NumericMatrix& Rcpp_point, NumericMatrix& Rcpp_evalpts, NumericMatrix& Rcpp_cluster_center,
double p, double pw,
double w, double c1, double c2,
int mM_part_num, int it_max, int mM_it_max,
int it_lim, int mM_it_lim,
double it_tol, double mM_it_tol,
double inn_tol, int inn_itmax,
int num_proc, double tol, double lb, double ub)
{
  // Description of Inputs:
  // Rcpp_point          - Clustering data for minimax clustering
  // Rcpp_evalpts        - Evaluation points for minimax post-processing
  // Rcpp_cluster_center - Initial cluster centers
  // p                   - q in paper; approximation coefficient for minimax criterion
  // pw                  - l_pw norm
  // it_max              - Maximum number of iterations for minimax clustering
  // mM_it_max           - Maximum number of iterations for minimax post-processing
  // mM_part_num         - Number of PSO particles to use for minimax post-processing
  // num_proc            - Number of processors to use in parallel computing
  // tol                 - Jitter tolerance
  // lb                  - Lower bound for design
  // ub                  - Upper bound for design

  int it_num; //keeps track of iteration number
	int dim_num = Rcpp_point.ncol(); //dimension of cluster data points
	int point_num = Rcpp_point.nrow(); //number of cluster data points
	int cluster_num = Rcpp_cluster_center.nrow(); //number of desired clusters
  int part_num = Rcpp_cluster_center.ncol()/dim_num; //part_num - number of particles

  //cluster energy matrix
  arma::mat cluster_energy(part_num,cluster_num);
  arma::mat cluster_tot(part_num,1); //totals the energy
  //keeps track of cluster assignment (each row for a particle)
  arma::mat cluster = -1*ones(part_num,point_num);
  //clustering points
  arma::mat point(Rcpp_point.begin(),point_num,dim_num,false);
  //current cluster centers (or design points)
  arma::cube cluster_center(Rcpp_cluster_center.begin(),cluster_num,dim_num,part_num);
  //current velocity
  arma::cube cluster_vel = zeros(cluster_num,dim_num,part_num);
  //lbest positions
  arma::cube cluster_lbes = zeros(cluster_num,dim_num,part_num);
  arma::mat lbes_obj = DBL_MAX*ones(1,part_num);
  //gbest positions
  arma::mat cluster_gbes = zeros(cluster_num,dim_num);
  double gbes_obj = DBL_MAX;
  double prev_gbes_obj = DBL_MAX;
  arma::rowvec assign_gbes(point_num);
  arma::rowvec energy_gbes(cluster_num);

  arma::mat Dtmp(cluster_num,dim_num); //temporary matrix for computation

  omp_set_num_threads(num_proc);

//  double w = 0.72;//inertia constant (from Merwe and Engelbrecht, 2003)
//  double c1 = 1.49;//acceleration constants
//  double c2 = 1.49;

//
//  Data checks.
//

  clock_t t = clock();

  it_num = 0;
  int it_same = 0;

//PSO iterations
//  cout << "-------------------------------------------------" << endl;
//  cout << "Step 1: PSO clustering ... " << endl;
//  cout << "-------------------------------------------------" << endl;


  while ((it_num<it_max)&&(it_same<it_lim))
  {
    Rcout << "PSO clustering iteration " << it_num << " ... " << endl;

    //Update partition, centroids and energy
//    cout << "kmeanspso...1" << endl;
    #pragma omp parallel for
    for (int i=0;i<part_num;++i){
//      cout << "kmeanspso...1.1" << endl;
      arma::rowvec tmp1 = cluster.row(i);
      arma::rowvec tmp2 = cluster_energy.row(i);
//      cout << "kmeanspso...1.2" << endl;
      kmeansreg(point,cluster_center.slice(i),tmp1,tmp2,p,pw,inn_tol,inn_itmax);
//      cout << "kmeanspso...1.3" << endl;
      cluster.row(i) = tmp1;
      cluster_energy.row(i) = tmp2;
    };
//    cout << "kmeanspso...2" << endl;

    //Update lbes and gbes
    cluster_tot = sum(cluster_energy,1);
//    cout << cluster_tot << endl;
    if (p > 2){
      cluster_tot = pow(cluster_tot,1/p);
    }

    bool ch_flg = false;
    for (int i=0;i<part_num;i++){
      if (cluster_tot(i)<lbes_obj(i)){
        cluster_lbes.slice(i) = cluster_center.slice(i);
        lbes_obj(i) = cluster_tot(i);
//        cout << "lbes_obj changed to:" << endl;
//        cout << lbes_obj << endl;
      }
      if (cluster_tot(i)<gbes_obj){
        //Update gbes
        cluster_gbes = cluster_center.slice(i);
        prev_gbes_obj = gbes_obj;
        gbes_obj = cluster_tot(i);
        assign_gbes = cluster.row(i); //update optimal assignment
        energy_gbes = cluster_energy.row(i); //update optimal energy
        Rcout <<  "gbes_obj changed to:" << endl;
        Rcout << gbes_obj << endl;

        //Update flag if improvement exceeds threshold
        if (prev_gbes_obj - gbes_obj < it_tol){
          ch_flg = true;
        }

      }

    }

    //Update or reset it_same
    if (ch_flg){
      it_same = 0;
    }
    else{
      it_same++;
    }
//    cout << "and here..." << endl;

//    //Update centroids using velocities
//    w = 0.72;//inertia constant (from Merwe and Engelbrecht, 2003)
//    c1 = 1.49;//acceleration constants
//    c2 = 1.49;

    for (int i=0; i<part_num; i++){
      Dtmp = cluster_center.slice(i);
      //update velocities
      cluster_vel.slice(i) = w*cluster_vel.slice(i) + c1*(randu(cluster_num,dim_num)%(cluster_lbes.slice(i)-Dtmp)) + c2*(randu(cluster_num,dim_num)%(cluster_gbes-Dtmp));
      //update centroids
      cluster_center.slice(i) += cluster_vel.slice(i);
    }

    //Update iteration count
    it_num ++;

    // cout << "Time :" << ((float)clock()-t)/CLOCKS_PER_SEC << endl;
    // cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
  }

//  cout << "--------------------------------------------------" << endl;
//  cout << "Clustering on global solution..." << endl;
//  cout << "--------------------------------------------------" << endl;
  kmeansreg(point,cluster_gbes,assign_gbes,energy_gbes,p,pw,1000,inn_itmax);

  //Post-process by doing PSO on mMcrit
//  cout << "-------------------------------------------------" << endl;
//  cout << "Step 2: PSO post-processing ... " << endl;
//  cout << "-------------------------------------------------" << endl;

  //(point, cluster_center, cluster, cluster_energy, p, it_max, int inn_itmax)
  gbes_obj = mMcritPSO(cluster_gbes,Rcpp_evalpts,mM_part_num,mM_it_max,mM_it_lim,mM_it_tol,tol,lb,ub,w,c1,c2);
//  cout << "--------------------------------------------------" << endl;
//  cout << "Final mM estimate:" << gbes_obj << endl;
//  cout << "--------------------------------------------------" << endl;

  //wrap to Rcpp classes for output
  NumericMatrix ret_cluster_center(cluster_num,dim_num);
//  NumericMatrix ret_lbes_obj(1,cluster_num);

  ret_cluster_center = armamatToRmat(cluster_gbes);
//  ret_lbes_obj = armamatToRmat(lbes_obj);

  return List::create(Named("gbes_centers")=ret_cluster_center, Named("gbes_obj")=gbes_obj);

  }
