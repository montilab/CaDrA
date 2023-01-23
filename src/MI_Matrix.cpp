//
// Created by bgregor on 9/19/22.
//

#include "MI_Matrix.h"

#include <cmath>


#include <algorithm>
#include <array>


#include <boost/math/special_functions/digamma.hpp>
#include "nanoflann.hpp"
#include <iostream>
#include <cmath>
#include <numeric>
#include <functional>
#include <list>
#include <limits>
#include <stdexcept>

#include "MutualInformation.h"

int mutual_inf_cc_vec(const double *input_x, const double *input_y, const int n_elems,  const int k, double *mi) {
    // Case all continuous input_y or mixed discrete/continuous
    // input_x - the vector of continuous data. size N.
    // input_y - input vector. size N
    // n_elems - N
    // mi - the return value. Must be pre-allocated to size 1!!
    // k - the neighborhood size. 3 is the recommended number.

    CaDrA::MutualInformation mut_inf(k) ;

    // Convert to Eigen matrices
    CaDrA::MapArrayConst x_eig(input_x,n_elems) ;
    CaDrA::MapArrayConst y_eig(input_y,n_elems) ;
    *mi = mut_inf.mutual_information_cc(x_eig, y_eig);
    return 0 ;
}

int mutual_inf_cd_vec(const double *input_x, const int *input_y, const int n_elems, const int k, double *mi ) {
    // Case all discrete input_y where the input vector is all integers.
    // This can use the faster algorithm for continuous-discrete calculations.
    // input_x - the vector of continuous data. size N.
    // input_y - input vector. size N
    // n_elems - N
    // mi - the return value. Must be pre-allocated to size 1!!
    // k - the neighborhood size. 3 is the recommended number.
    CaDrA::MutualInformation mut_inf(k) ;

    // Convert to Eigen matrices
    CaDrA::MapArrayConst x_eig(input_x,n_elems) ;
    CaDrA::MapArrayIConst y_eig(input_y,n_elems) ;
    mi[0] = mut_inf.mutual_information_cd(x_eig, y_eig);
    return 0 ;
}

int cond_mutual_inf_vec(const double *input_x,  const double *input_y, const double *input_z, const int n_elems, const int k, double *mi) {
    // Conditional mutual information for a single vector of x,y,and z
    // input_x, input_y, input_z - input vectors, all of size n_elems.
    // k - number of nearest neighbors
    // mi - return value.
    CaDrA::MutualInformation mut_inf(k) ;
    // Convert to Eigen matrices
    CaDrA::MapArrayConst x_eig(input_x,n_elems) ;
    CaDrA::MapArrayConst y_eig(input_y,n_elems) ;
    CaDrA::MapArrayConst z_eig(input_z,n_elems) ;
    mi[0] = mut_inf.cond_mutual_information_ccc(x_eig, y_eig, z_eig);
    return 0 ;
}


