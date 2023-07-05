//
// Created by bgregor on 11/1/22.
//
#include "CondMutualInformation.h"

#include <cmath>
#include <array>
#include <numeric>
#include <functional>
#include <list>
#include <set>
#include <map>

#include <Rmath.h>
#include "R.h"
#include "Rinternals.h"
#include <R_ext/Random.h>

#include "nanoflann.hpp"

namespace CaDrA {

CondMutualInformation::CondMutualInformation(const int mK) : MutualInformationBase(mK) {}

CondMutualInformation::~CondMutualInformation() {}

double CondMutualInformation::compute(const ArrayXd &x, const ArrayXd& y, const ArrayXd& z) {
  // This implements the CMI algorithm described in: https://doi.org/10.1016/j.eswa.2012.05.014
  // 
  // Alkiviadis Tsimpiris, Ioannis Vlachos, Dimitris Kugiumtzis,
  // Nearest neighbor estimate of conditional mutual information in feature selection,
  // Expert Systems with Applications,
  // Volume 39, Issue 16, 2012, Pages 12697-12708
  //
  // Compute the conditional mutual information.
  // CMI = psi(m_k) - mean( psi(n_xz[i] + psy(n_yz[i])-psi(n_z[i]) ) for all elements i
  auto N = x.size() ;
  
  // Get the 3-dimensional distance, and then re-use that for the xz, yz, and calculations, so it
  // proceeds pretty much exactly as before.
  
  Array3col tmp_mat(N, 3) ;
  
  tmp_mat.col(0) = scale(x) ;
  tmp_mat.col(1) = scale(y) ;
  tmp_mat.col(2) = scale(z) ;
  
  // Calculating the distances also calculates the number of neighbors, so
  // calc_distances2 returns both.
  vector<double> dists = calc_distances3d(N, tmp_mat).first ;
  
  // Get the digamma_f values...These take single vector arguments
  // so map the 3D temp array without a copy.
  MapArrayConst x_scale(tmp_mat.col(0).data(), N) ;
  MapArrayConst y_scale(tmp_mat.col(1).data(), N) ;
  MapArrayConst z_scale(tmp_mat.col(2).data(), N) ;
  double xz_digamma_sum = sum_digamma_from_neighbors(x_scale, z_scale, dists) ;
  double yz_digamma_sum = sum_digamma_from_neighbors(y_scale, z_scale, dists) ;
  double z_digamma_sum = MutualInformationBase::sum_digamma_from_neighbors(z_scale, dists) ;
  
  // mutual info computation
  double mi = digamma_f(m_k)- (xz_digamma_sum + yz_digamma_sum - z_digamma_sum) / N;
  
  // Can't return less than 0.
  return std::max(0.0,mi) ;
}

double CondMutualInformation::sum_digamma_from_neighbors(MapArrayConst &vec1, MapArrayConst &vec2, const vector<double> &dists) {
  // Sum of digamma_f functions over neighbor counts for 2D.
  long N = dists.size() ;
  double sum = 0.0 ;
  
  // KD-Tree for this vector
  Array2col tmp_mat(N, 2) ;
  tmp_mat.col(0) = vec1 ;
  tmp_mat.col(1) = vec2 ;
  nanoflann::KDTreeEigenMatrixAdaptor<Array2col ,-1,metric_Chebyshev> vec_tree(2, tmp_mat, 10) ;
  
  std::vector<std::pair<Eigen::Index, double>> ret_matches;
  array<double,2> pt ;
  for (long i = 0 ; i < N ; ++i) {
    pt[0] = tmp_mat(i,0) ;
    pt[1] = tmp_mat(i,1) ;
    double tmp = vec_tree.index->radiusSearch(pt.data(), dists[i] , ret_matches , nanoflann::SearchParams(10));
    sum += digamma_f(tmp) ;
    ret_matches.clear() ;
  }
  return sum ;
}


pair<vector<double>,vector<long>>  CondMutualInformation::calc_distances3d(const long N, 
                                                                       const Array3col &tmp_mat) const {
  // Calculate the Chebyshev distances and numbers of neighbors for the 3D array tmp_mat.
  kd_tree_3d mat_index(tmp_mat.cols(),std::cref(tmp_mat),20) ;
  // We want N neighbors in addition to the point itself so
  // add 1 to the # of neighbors.
  int real_k = m_k + 1  ;
  
  // Chebyshev distance
  vector<double> dists(N) ;
  // Number of neighbors
  vector<long> neighbors(N) ;
  
  // a query point.
  array<double,3> query_pt ;
  for (long i = 0 ; i < N ; ++i) {
    // store indexes and distances
    vector<Eigen::Index> ret_indexes(real_k, 0.0);
    vector<double> out_dists_sqr(real_k,0.0);
    
    for (long j = 0 ; j < 3 ; ++j)
      query_pt[j] = tmp_mat(i, j);
    
    neighbors[i] = mat_index.index->knnSearch(&query_pt[0], real_k,
                                              &ret_indexes[0], &out_dists_sqr[0]) ;
    out_dists_sqr.resize(neighbors[i]) ;
    dists[i] = out_dists_sqr.back() ;
  }
  return make_pair(dists,neighbors) ;
}

} // CaDrA
