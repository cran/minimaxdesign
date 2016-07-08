#include <iostream>
#include <cmath>
#ifdef _OPENMP
# include <omp.h>
#endif
#include <RcppArmadillo.h>
#include "header.h"

using namespace std;
using namespace arma;
using namespace Rcpp;


//*****************************************************
//Transformation functions
//*****************************************************
//--------------------------------------------------------------
//Simplex: Transforms [0,1]^p to A^p
//--------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix CtoA(NumericMatrix& D, double by, int num_proc) {
  // Description of Input
  // D        - Low-discrepancy points on [0,1]^p
  // by       - Approximation step-size
  // num_proc - Number of processors to use in parallel computing
  NumericMatrix ret(D.nrow(),D.ncol());
  double tmp = 0; //temporary variable for computation
  for(int i=0;i<D.nrow();i++){
    tmp = 1; //reset tmp
    for(int j=D.ncol()-1;j>-1;j--){
      tmp = pow(D(i,j),(1/double(j+1)))*tmp;
      ret(i,j) = tmp;
    }
  }
  return(ret);
}

//--------------------------------------------------------------
//Ball: Transforms [0,1]^2 to B^2
//--------------------------------------------------------------
// [[Rcpp::export]]
NumericMatrix CtoB2(NumericMatrix& D, double by, int num_proc) {
  // Description of Input
  // D        - Low-discrepancy points on [0,1]^p
  // by       - Approximation step-size
  // num_proc - Number of processors to use in parallel computing
  NumericMatrix ret(D.nrow(),D.ncol());
  for(int i=0;i<D.nrow();i++){
    ret(i,0)=sqrt(D(i,0))*cos(2*M_PI*D(i,1));
    ret(i,1)=sqrt(D(i,0))*sin(2*M_PI*D(i,1));
  }
  return(ret);
}

//*****************************************************
//Conversion and print functions
//*****************************************************
void printArr(int arr[],int len){
  for (int i=0;i < len;i++){
    Rcout << arr[i] << endl;
  }
}

void printRcppMat(Rcpp::NumericMatrix& matr){
  for (int i=0;i<matr.nrow();i++){
    for (int j=0; j<matr.ncol();j++){
      Rcout << matr(i,j) << '\t';
    }
    Rcout << endl;
  }
}

//--------------------------------------------------------------
//Wrapping integer array in NumericVector
//--------------------------------------------------------------
Rcpp::NumericVector iarrToRvec(int arr[], int len){
  Rcpp::NumericVector ret(len);
  for (int i=0;i<len;i++){
    ret[i] = arr[i];
  }
  return ret;
}

//--------------------------------------------------------------
//Wraps double array in NumericVector
//--------------------------------------------------------------
Rcpp::NumericVector darrToRvec(double arr[], int len){
  //wraps double array in NumericVector
  Rcpp::NumericVector ret(len);
  for (int i=0;i<len;i++){
    ret[i] = arr[i];
  }
  return ret;
}

//--------------------------------------------------------------
//Wraps double array in NumericVector
//--------------------------------------------------------------
Rcpp::NumericMatrix armamatToRmat(arma::mat& matr){
  //wraps double array in NumericVector
  Rcpp::NumericMatrix ret(matr.n_rows,matr.n_cols);
  for (int i=0;i<matr.n_rows;i++){
    for (int j=0;j<matr.n_cols;j++){
      ret(i,j) = matr(i,j);
      //      cout << ret[i,j] << '\t';
    }
    //    cout << endl;
  }
  //  printRcppMat(ret);
  return ret;
}

//*****************************************************
//Sub-functions for kmeansreg
//*****************************************************
//--------------------------------------------------------------
//AGD for computing C_q centers
//--------------------------------------------------------------
arma::rowvec cq_agd(arma::mat& pts, double q, double eps, int it_max) {
  int n = pts.n_rows;
  int p = pts.n_cols;
  int it_num = 0;
  double beta = 0.0;
  double dst = 0.0;
  vector<double> maxvec(n,0.0); //Records the maximum distance within the group

  //Compute beta
  //  for (int i=0; i<n; i++){
  //    for (int j=(i+1); j<n; j++){
  //      dst = 0.0;
  //      for (int k=0; k<p; k++){
  //        dst += pow(pts(i,k) - pts(j,k), 2.0);
  //      }
  //      dst = pow(dst,0.5);
  //      if (dst > maxvec[i]){
  //        maxvec[i] = dst;
  //      }
  //      if (dst > maxvec[j]){
  //        maxvec[j] = dst;
  //      }
  //    }
  //  }
  //
  //  for (int i=0; i<n; i++){
  //    beta += pow(maxvec[i], q/2-1);
  //  }
  //  beta = beta*(q-1)/((double) n);

  beta = q-1;

  //  cout << "Here 0" << endl;

  //Initialize objects
  vector<double> curz(p,0.0); //Current intermediate point
  for (int i=0; i<n; i++){
    for (int j=0; j<p; j++){
      curz[j] += (1.0/(double)n)*pts(i,j);
    }
  }

  //  cout << "Here 1" << endl;

  vector<double> curu(p,0.0); //Current point
  double lambdat = 1.0; // lambda_t
  double lambdatp1 = 1.0; //lambda_{t+1}
  for (int j=0; j<p; j++){
    curu[j] = curz[j];
  }

  vector<double> prevu(p,0.0); //Previous point
  double gammat = (1.0 - lambdat)/lambdatp1; // gamma_t

  vector<double> wts(n,0.0); //weights for computing gradient
  vector<double> tmp(p,0.0);
  bool flag = true; //flag for error tolerance

  //  cout << "Here 2" << endl;

  //Do Nesterov updates until convergence
  while ( (flag) && (it_num<it_max) ){
    //Update lambda and gamma
    lambdat = lambdatp1;
    lambdatp1 = (1.0 + pow(1.0 + 4.0*pow(lambdat,2.0), 0.5 )) / 2.0;
    gammat = (1.0 - lambdat)/lambdatp1;

    //Update wts, u and z
    for (int i=0; i<n; i++){//Update wts
      wts[i] = 0.0;
      for (int j=0; j<p; j++){
        wts[i] += pow( curz[j]-pts(i,j), 2.0 );
      }
      wts[i] = (1.0/(double)n) * (1.0/beta) * pow(wts[i], q/2.0 - 1.0);
    }

    for (int j=0; j<p; j++){//Update prevu and prevz
      prevu[j] = curu[j];
      curu[j] = curz[j];
    }

    for (int i=0; i<n; i++){//Update prevu
      for (int j=0; j<p; j++){
        curu[j] -= wts[i] * (curz[j] - pts(i,j));
        //                  tmp[j] -= wts[i] * (curz[j] - pts(i,j));
      }
    }

    for (int j=0; j<p; j++){
      curz[j] = (1-gammat)*curu[j] + gammat*prevu[j];
    }

    //Update flags
    double dst = 0.0;
    for (int j=0; j<p; j++){
      dst += pow(curu[j]-prevu[j],2.0);
    }
    dst = pow(dst, 0.5);
    if (dst < eps){
      flag = false;
    }

    //Increment
    it_num++;
  }

  //  cout << "Here 3" << endl;

  arma::rowvec ret(p);
  for (int j=0; j<p; j++){
    ret(j) = curz[j];
  }

  return(ret);
}

//--------------------------------------------------------------
//Searches for all indices with value seek
//--------------------------------------------------------------
uvec find(int arr[], int len, int seek)
{
  uvec ret(len);
  ret.fill(len+1); // initialize at len+1's
  int count = 0; //count index
  for (int i = 0; i < len; ++i)
  {
    if (arr[i] == seek) {
      ret(count) = i;
      count += 1;
    }
  }
  if (count == 0){
    //      cout << "I found nothing" << endl;
    return(ret);
  }
  else {return (ret(span(0,count-1)));}
}

//--------------------------------------------------------------
//Searches for all indices with value seek
//--------------------------------------------------------------
uvec find(arma::rowvec& arr, int seek)
{
  int len = arr.n_elem;
  bool empty = true;
  uvec ret(len);
  ret.fill(0); // initialize at 0's
  int count = 0; //count index
  for (int i = 0; i < len; ++i){
    if (arr(i) == seek) {
      ret(count) = i;
      count += 1;
      empty = false;
    }
  }
  if (empty){
    uvec emptyret(1);
    emptyret << UINT_MAX;
    return(emptyret);
  }
  else{
    return(ret(span(0,count-1)));
  }
}

//--------------------------------------------------------------
//Searches for all indices with value seek
//--------------------------------------------------------------
void find(arma::mat& arr, arma::mat& point, int seek, list<arma::mat>& L)
{
  //arr - cluster indices
  //point - data points
  //seek - index to search
  //L - list to return
  bool empty;
  int count = 0;
  list<uvec> idxList;
  //  cout << "arr: " << arr.n_rows << ", " << arr.n_cols << endl;
  for (int i = 0; i < arr.n_rows; ++i){ //for each theta value
    empty = true;
    uvec idx(arr.n_cols);
    idx.fill(0);
    count = 0;
    for (int j = 0; j < arr.n_cols; ++j){ //for each cluster index
      if (arr(i,j) == seek){
        idx(count) = j;
        count += 1;
        empty = false;
      }
    }
    if (empty){//add empty index vector
      uvec emptyidx(1);
      emptyidx << UINT_MAX;
      idxList.push_back(emptyidx);
    }
    else{
      idxList.push_back(idx(span(0,count-1)));
    }
  }

  list<uvec>::iterator k;
  for (k=idxList.begin(); k!=idxList.end(); ++k){ //for each theta value
    if ((*k)(1)==UINT_MAX){
      arma::mat emptyMat(1,1);
      emptyMat << -1;
      L.push_back(emptyMat);
    }
    else{
      L.push_back(point.rows(*k));
    }
  }

}

//--------------------------------------------------------------
//One step of modified K-means (not parallel)
//--------------------------------------------------------------
// [[Rcpp::plugins(openmp)]]
void kmeansreg (arma::mat& point, arma::mat& cluster_center, arma::rowvec& cluster, arma::rowvec& cluster_energy,
                double p, double pw, double inn_tol, int inn_itmax)
{
  // Description of Inputs:
  // point          - Clustering data for minimax clustering
  // cluster_center - Initial cluster centers
  // cluster        - Keeps track of which points are in which cluster
  // cluster_energy - Keeps track of clustering objective in each cluster
  // p              - q in paper; approximation coefficient for minimax criterion
  // it_max         - Maximum iteration for AGD

  int dim_num = point.n_cols; //dimension of data points
  int point_num = point.n_rows; //number of data points
  int cluster_num = cluster_center.n_rows; //number of desired clusters
  //int cluster[point_num]; //keeps track of cluster assignment
  //  memset(cluster, 0, sizeof cluster);
  //  arma::mat cluster_energy(1,cluster_num);

  int j;
  int k;
  int k2;
  double point_energy = 0;
  double point_energy_min;
  double tmp_energy;
  int swap;

  arma::mat Dtmp; //temporary matrix object to compute centroid
  arma::rowvec rowtmp(dim_num);

  //
  //  If cluster is uninitialized, then assign clusters
  //

  //  cout << "kmeansreg...1" << endl;

  if (cluster(0) < 0){
    for ( j = 0; j < point_num; j++ ) // for each data point
    {
      point_energy_min = DBL_MAX;
      rowtmp = point.row(j);

      for ( k = 0; k < cluster_num; k++ ) // for each cluster
      {
        point_energy=0;
        for (int m=0;m<dim_num;m++){
          // //          point_energy += pow(abs(rowtmp(m)-cluster_center(k,m)),pw);
          point_energy += pow(rowtmp(m)-cluster_center(k,m),2);
        }
        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster[j] = k;
        }

      }
    }
  }


  //
  //  #1:
  //  Determine the C^q centers of the clusters by AGD
  //

  //  cout << "kmeansreg...2" << endl;

  for (k=0;k<cluster_num;k++){ // for each cluster
    uvec idx = find(cluster,k);
    //    cout << "kmeansreg...2: " << k << ",1" << endl;
    if (idx[1]>double(point_num)){
      //cluster is empty, assign random center
      cluster_center.row(k) = randu(1,dim_num);
    }
    else{
      //      cout << "kmeansreg...2: " << k << ",2" << endl;
      //      cout << "idx: " << idx.size() << ", " << idx[0] << endl;
      Dtmp = point.rows(idx);
      //      cout << "kmeansreg...2: " << k << ",3" << endl;
      cluster_center.row(k) = cq_agd(Dtmp,p,inn_tol,inn_itmax);
      //      cout << "kmeansreg...2: " << k << ",4" << endl;
    }
  }

  //
  //  #2:
  //  Assign each point to the cluster of its nearest center.
  //

  //  cout << "kmeansreg...3" << endl;

  swap = 0;

  for ( j = 0; j < point_num; j++ ) // for each data point
  {
    point_energy_min = DBL_MAX;
    k = cluster(j);
    rowtmp = point.row(j);

    for ( k2 = 0; k2 < cluster_num; k2++ ) // for each cluster
    {
      point_energy=0;
      for (int m=0;m<dim_num;m++){
        // //          point_energy += pow(abs(rowtmp(m)-cluster_center(k2,m)), pw);
        point_energy += pow(rowtmp(m)-cluster_center(k2,m), 2);
      }

      if ( point_energy < point_energy_min )
      {
        point_energy_min = point_energy;
        cluster(j) = k2;
      }
    }

    if ( k != cluster(j) )
    {
      swap = swap + 1;
    }
  }

  //
  //  #3:
  //  Determine the total energy of the current clustering with new centroids.
  //

  cluster_energy.fill(0); //Reset all values in cluster_energy

  for ( j = 0; j < point_num; j++ ){//for each point
    k = cluster(j);
    tmp_energy = 0;
    for (int m=0;m<dim_num;m++){
      // //        tmp_energy += pow(abs(point(j,m)-cluster_center(k,m)),pw);
      tmp_energy += pow(point(j,m)-cluster_center(k,m),2);
    }
    cluster_energy(k) += pow(tmp_energy,p/2.0);
  }
}

//--------------------------------------------------------------
//Multiple steps of modified K-means (parallel)
//--------------------------------------------------------------
// [[Rcpp::plugins(openmp)]]
void kmeansreg (arma::mat& point, arma::mat& cluster_center, arma::rowvec& cluster, arma::rowvec& cluster_energy,
                double p, double pw, int it_max, double inn_tol, int inn_itmax)
{
  // Description of Inputs:
  // point          - Clustering data for minimax clustering
  // cluster_center - Initial cluster centers
  // cluster        - Keeps track of which points are in which cluster
  // cluster_energy - Keeps track of clustering objective in each cluster
  // p              - q in paper; approximation coefficient for minimax criterion
  // it_max         - Maximum number of iterations for minimax clustering

  int it_num; //keeps track of iteration number
  int dim_num = point.n_cols; //dimension of data points
  int point_num = point.n_rows; //number of data points
  int cluster_num = cluster_center.n_rows; //number of desired clusters
  //int cluster[point_num]; //keeps track of cluster assignment
  //  memset(cluster, 0, sizeof cluster);
  //  arma::mat cluster_energy(1,cluster_num);

  int it_same = 0; //keeps track of same criterion energy
  bool ind_stop = TRUE;
  double min_sum = DBL_MAX;
  arma::rowvec min_energy;
  arma::rowvec min_cluster(cluster.n_rows);
  arma::mat min_center(cluster_num,dim_num);

  int j;
  int k;
  int k2;
  double point_energy = 0;
  double point_energy_min;
  double tmp_energy;
  int swap;

  arma::rowvec rowtmp(dim_num);
  arma::rowvec tmprowvec(dim_num);

  //
  //  If cluster is uninitialized, then assign clusters
  //

  //  cout << "test 0" << endl;
  if (cluster(0) < 0){
    for ( j = 0; j < point_num; j++ ) // for each data point
    {
      point_energy_min = DBL_MAX;
      rowtmp = point.row(j);

      for ( k = 0; k < cluster_num; k++ ) // for each cluster
      {
        point_energy=0;
        for (int m=0;m<dim_num;m++){
          // //          point_energy += pow(abs(rowtmp(m)-cluster_center(k,m)),pw);
          point_energy += pow(rowtmp(m)-cluster_center(k,m),2);
        }
        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster[j] = k;
        }

      }
    }
  }


  //  cout << "test 1" << endl;
  it_num = 0;

  while ((it_num<it_max)&&(ind_stop))
  {
    Rcout << "Post-process (cluster): iteration " << it_num << " ... " << endl;
    //  Increment iteration count
    it_num = it_num + 1;
    //    cout << "test 1.1" << endl;

    //
    //  #1:
    //  Determine C^q-centers of the clusters by L-BFGS
    //

#pragma omp parallel for
    for (k=0;k<cluster_num;k++){ // for each cluster
      uvec idx = find(cluster,k);
      if (idx[1]==UINT_MAX){
        //cluster is empty, assign random center
        cluster_center.row(k) = randu(1,dim_num);
      }
      else{
        arma::mat Dtmp = point.rows(idx);
        //update cluster centers
        cluster_center.row(k) = cq_agd(Dtmp,p,inn_tol,inn_itmax);
      }
    }
    //    cout << "test 1.2" << endl;


    //
    //  #2:
    //  Assign each point to the cluster of its nearest center.
    //

    swap = 0;

    for ( j = 0; j < point_num; j++ ) // for each data point
    {
      point_energy_min = DBL_MAX;
      k = cluster(j);
      rowtmp = point.row(j);

      for ( k2 = 0; k2 < cluster_num; k2++ ) // for each cluster
      {
        point_energy=0;
        for (int m=0;m<dim_num;m++){
          // //          point_energy += pow(abs(rowtmp(m)-cluster_center(k2,m)),pw);
          point_energy += pow(rowtmp(m)-cluster_center(k2,m),2);
        }

        if ( point_energy < point_energy_min )
        {
          point_energy_min = point_energy;
          cluster(j) = k2;
        }
      }

      if ( k != cluster(j) )
      {
        swap = swap + 1;
      }
    }

    //    cout << "test 1.3" << endl;
    // Rcout << "Swaps: " << swap << endl;
    // // Rcout << "Min. objective: " << sum(cluster_energy) << endl;
    // Rcout << "it_same: " << it_same << endl;

    //
    //  Terminate if no points were swapped.
    //
    if (1 < it_num){
      if (swap == 0){
        break;
      }
    }

    //
    //  #3:
    //  Determine the total energy of the current clustering with new centroids.
    //

    cluster_energy.fill(0); //Reset all values in cluster_energy

    for ( j = 0; j < point_num; j++ ){//for each point
      k = cluster(j);
      tmp_energy = 0;
      for (int m=0;m<dim_num;m++){
        // //        tmp_energy += pow(abs(point(j,m)-cluster_center(k,m)),pw);
        tmp_energy += pow(point(j,m)-cluster_center(k,m),2);
      }
      cluster_energy(k) += pow(tmp_energy,p/2.0);
    }

    //
    //  #4:
    //  Stop if objective does not change for a given number of times
    //

    //Update minimum energy
    if (sum(cluster_energy)<min_sum){
      min_sum = sum(cluster_energy);
      min_energy = cluster_energy;
      min_cluster = cluster;
      min_center = cluster_center;
      it_same = 0; //reset count
    }
    else{
      it_same++;
    }

    //Stop if min_energy does not change for iterations
    if (it_same >= 5){
      ind_stop = FALSE;
    }

  }

  //Finally, return the minimum values
  cluster_energy = min_energy;
  cluster = min_cluster;
  cluster_center = min_center;

}
