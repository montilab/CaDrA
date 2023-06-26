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

 

double MutualInformation::compute_c(const ArrayXd &x, const ArrayXd& y) {
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
    
  return std::max(0.0,std::min(mi,1.0)) ;
}

double MutualInformation::compute_d(const ArrayXd &x, const ArrayXi &y) {
  // Compute the mutual information for a continuous vector x and a
  // discrete vector y.
  // This implements the algorithm described in: https://doi.org/10.1371/journal.pone.0087357
  // Ross BC (2014) Mutual Information between Discrete and Continuous Data Sets. PLoS ONE 9(2): e87357. 
  
  auto N = x.size() ;
  ArrayXd x_scale = scale(x) ;
  // Make a kdtree for x_scale, it'll be needed later.
  kd_tree_1d xscale_index_tree(1, x_scale, 10);
  
  // Chebyshev distance
  ArrayXd dists(N) ;
  // Number of neighbors
  ArrayXd neighbors(N) ;
  
  // Get the unique labels in y by creating an STL set
  std::set<int> unique_labels{y.data(), y.data() + y.size()};

  // For each unique label store its count.
  vector<double> label_counts ;
  
  // Store the k's used for each label
  vector<vector<double>> k_all ;
  
  // Store the neighbor counts m for each label
  vector<double> m_all ;
  // master vector of all indices that are not of unique labels.
  vector<int> all_indices ;
  
  for (const auto label : unique_labels) {
    int count = (y == label).count();
    if (count > 1) {
      // Only process the non-unique labels
      label_counts.push_back(count) ;
      // Adjust real_k as necessary.
      auto real_k = min(m_k, count) + 1 ;
      vector<double> tmp = {static_cast<double>(count),static_cast<double>(real_k)};
      k_all.push_back(tmp) ;
      // Store the indices of y where they match
      // the current label.
      vector<int> label_indices(count);
      int j = 0 ;
      // Eigen doesn't seem to have a convenient way to do this
      // so just use a for loop.
      for (int i = 0; i < N; ++i) {
        if (y[i] == label) {
          label_indices[j] = i ; ++j ;
          all_indices.push_back(i) ;
        }
      }
      // Use Eigen 3.4's method of providing a vector of indices
      // to produce a sub-vector.
      ArrayXd masked_x = x_scale(label_indices) ;
      // Make a lookup tree for the points of this label.
      kd_tree_1d label_index_tree(1, masked_x, 10);
      // Get all of the distances for each point for this label.
      for (int i = 0; i < count; ++i) {
        vector<Eigen::Index> ret_indexes(real_k, 0.0);
        vector<double> out_dists(real_k,0.0);
        
        double query_pt[1] = {masked_x[i]};
        // out_dists stores the distances for this label. Get the max one.
        auto neighbors = label_index_tree.index->knnSearch(query_pt, 
                                                           real_k,
                                                           &ret_indexes[0], 
                                                           &out_dists[0]);
        out_dists.resize(neighbors) ;
        // The last one is out_dists is the furthest distance.
        auto max_dist = out_dists.back() ;
        
        std::vector<std::pair<Eigen::Index, double>> ret_matches;
        m_all.push_back(xscale_index_tree.index->radiusSearch(query_pt, 
                                                              max_dist, 
                                                              ret_matches , 
                                                              nanoflann::SearchParams(10))) ;
      }
    }
  }
  double N_mod = all_indices.size() ;
  double digamma_N = digamma_f(N) ;
  // Calculate the mean of digammas over the count of samples for each label.
  double digamma_labels = std::accumulate(label_counts.begin(),label_counts.end(),0.0,
                                          [&](double m, double n){ return m + 
                                            digamma_f(n) * n / N_mod ; } ) ;
  // Get the same for the k's used at each label.
  double digamma_k = std::accumulate(k_all.begin(),k_all.end(),0.0,
                                     [&](double m, vector<double> &n){ return m +
                                       digamma_f(std::max(n[1] - 1.0, 1.0)) * n[0] / N_mod; } ) ;
  double digamma_m = std::accumulate(m_all.begin(),m_all.end(), 0.0,  
                                     [&](double m, double n){return m + digamma_f(n); });
  digamma_m = digamma_m / N_mod ;
  
  // mutual info computation
  // Matlab:   4.602 - 4.3965 + 0.9228 - 1.0961
  double mi = digamma_N - digamma_labels + digamma_k - digamma_m ;
  return std::max(0.0,std::min(mi,1.0)) ;
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
