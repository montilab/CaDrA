//
// Created by bgregor on 11/1/22.
//
#include "MutualInformationBase.h"

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

#include <iostream>

namespace CaDrA {

MutualInformationBase::MutualInformationBase(const int mK) : m_k(mK) {
  // https://cran.r-project.org/doc/manuals/R-exts.html#Random-numbers
  GetRNGstate();
}

MutualInformationBase::~MutualInformationBase() {
  PutRNGstate();
}

int MutualInformationBase::get_k() const {
  return m_k;
}


void MutualInformationBase::set_k(int mK) {
  m_k = mK;
}
 

ArrayXd MutualInformationBase::scale(const ArrayXd &x, bool add_noise) const {
  // Scale the vector x by its standard deviation. Add a bit of random 
  // noise.
  auto size_v = x.size() ;
  double mean_x = x.mean() ;
  double sum = 0.0 ;
  for (auto i = 0 ; i < size_v ; i++) {
    double tmp = x[i] - mean_x ;
    sum += tmp * tmp ;
  }
  double std_dev = std::sqrt(sum / (size_v-1)) ;
  ArrayXd x_scale = x  / std_dev ;
  
  // add a wee bit of noise as suggested in
  // Kraskov et. al. 2004.
  if (add_noise) {
    double mean_xs = x_scale.mean() ;
    // Get a uniform value from the R RNG
    // https://cran.r-project.org/doc/manuals/R-exts.html#Random-numbers
    for (auto i = 0 ; i < x_scale.size() ; i++) {
      x_scale[i] += 1e-10 * mean_xs * unif_rand();
    }
  }
  return x_scale;
}

double MutualInformationBase::sum_digamma_from_neighbors(MapArrayConst &vec, const vector<double> &dists) {
  // This one is called from mutual_information_cc and cond_mutual_information for the neighbors
  // for a single vector.
  long N = dists.size() ;
  double sum = 0.0 ;
  
  // KD-Tree for this vector
  nanoflann::KDTreeEigenMatrixAdaptor<MapArrayConst,-1,metric_Chebyshev> vec_tree(1, vec, 10) ;
  
  std::vector<std::pair<Eigen::Index, double>> ret_matches;
  for (long i = 0 ; i < N ; ++i) {
    double pt = vec(i) ; // avoids type issues with the compiler and the radiusSearch.
    double tmp = vec_tree.index->radiusSearch(&pt, dists[i] , ret_matches , nanoflann::SearchParams(10));
    sum += digamma_f(tmp) ;
    ret_matches.clear() ;
  }
  return sum ;
}


// A wrapper for R's digamma function. This wrapper exists because
// the initial development was done using a standalone C++ program
// that used the Boost library's digamma function.
double MutualInformationBase::digamma_f(const double x) const {
    return digamma(x) ;
}

// placeholder
double MutualInformationBase::compute(){ return 0.0 ;}

} // CaDrA
