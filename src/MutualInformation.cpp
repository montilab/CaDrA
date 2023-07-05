//
// Created by bgregor on 11/1/22.
//
#include "MutualInformation.h"

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

MutualInformation::MutualInformation(const int mK) : MutualInformationBase(mK) {}

MutualInformation::~MutualInformation() {}

 

double MutualInformation::compute(const ArrayXd &x, const ArrayXd& y) {
  // compute mutual information based on k-nearest neighbors
  // This implements the algorithm described in: https://doi.org/10.1103/PhysRevE.69.066138
  // Alexander Kraskov, Harald Stogbauer, and Peter Grassberger
  // Phys. Rev. E 69, 066138 ?? Published 23 June 2004; Erratum Phys. Rev. E 83, 019903 (2011)
  
  // Creat a KD-tree using nanoflann and its Eigen matrix adapter
  // Store the size of the vectors for easy reference
  
  auto N = x.size() ;
  
  Array2col  tmp_mat(N, 2) ;
  tmp_mat.col(0) = scale(x) ;
  tmp_mat.col(1) = scale(y) ;
  // Map the double array pointer to an Eigen vector without a copy.
  MapArrayConst x_scale(tmp_mat.col(0).data(), N) ;
  MapArrayConst y_scale(tmp_mat.col(1).data(), N) ;
  
  vector<double> dists  = calc_distances2d(N, tmp_mat).first;
  
  double x_digamma_sum = sum_digamma_from_neighbors(x_scale, dists) ;
  double y_digamma_sum = sum_digamma_from_neighbors(y_scale, dists) ;
  
  // mutual info computation
  double mi = digamma_f(N) + digamma_f(m_k) - (x_digamma_sum + y_digamma_sum) / N;
  
  // Can't return less than 0.
  return std::max(0.0,mi) ;
}



pair<vector<double>,vector<long>> MutualInformation::calc_distances2d(const long N, 
                                                                      const Array2col &tmp_mat) const {
  // Calculate the Chebyshev distances and numbers of neighbors for the 2D array tmp_mat.
  kd_tree_2d mat_index(tmp_mat.cols(),std::cref(tmp_mat),20) ;
  // We want N neighbors in addition to the point itself so
  // add 1 to the # of neighbors.
  int real_k = m_k + 1 ;
  
  // Chebyshev distance
  vector<double> dists(N) ;
  // Number of neighbors
  vector<long> neighbors(N) ;
  
  // a query point.
  array<double,2> query_pt ;
  for (long i = 0 ; i < N ; ++i) {
    // store indexes and distances
    vector<Eigen::Index> ret_indexes(real_k, 0.0);
    vector<double> out_dists_sqr(real_k,0.0);
    
    query_pt[0] = tmp_mat(i, 0);
    query_pt[1] = tmp_mat(i, 1);
    
    neighbors[i] = mat_index.index->knnSearch(&query_pt[0], real_k,
                                              &ret_indexes[0], &out_dists_sqr[0]) ;
    out_dists_sqr.resize(neighbors[i]) ;
    auto max_dist= std::nextafter(*max_element(std::begin(out_dists_sqr), std::end(out_dists_sqr)),0.0)  ; // following sklearn
    dists[i] = max_dist ;
  }
  return make_pair(dists,neighbors) ;
}


} // CaDrA
