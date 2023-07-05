//
// Created by bgregor on 11/1/22.
//

#ifndef REVEALER_MUTUALINFORMATIONDISCRETE_H
#define REVEALER_MUTUALINFORMATIONDISCRETE_H

#include <eigen3/Eigen/Core>
#include <vector>

#include "MutualInformationBase.h"
#include "nanoflann.hpp"
#include "ChebyshevMetric.h"

using namespace Eigen;
using namespace std ;

namespace CaDrA {
 

class MutualInformationDiscrete : public MutualInformationBase {
  // Methods to compute mutual information for
  // continuous and discrete variables.
public:
  // Constructor
  // Initialize with a neighbor size.
  MutualInformationDiscrete(const int k) ;
  MutualInformationDiscrete() = delete ;
  
  // Destructor
  virtual ~MutualInformationDiscrete() ;
  
  //  Computes the Mutual Information of a continuous and a discrete variable
  virtual double compute(const ArrayXd &x, const ArrayXi& y) ;
  
protected:
  // Calculate distances and nearest neighbors in 2D.
  pair<vector<double>,vector<long>>  calc_distances2d(const long N, const Array2col &tmp_mat) const;
};

} // CaDrA

#endif //REVEALER_MUTUALINFORMATIONDISCRETE_H
