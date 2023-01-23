//
// Created by bgregor on 9/19/22.
//

#ifndef REVEALER_MI_MATRIX_H
#define REVEALER_MI_MATRIX_H

// This provides C-callable wrappers for the mutual information implementation
// which is in C++. The functions here are quite simple, they wrap the incoming
// arrays as Eigen Array types to call the

extern "C" {
        // and for when the input score is just 1 vector
        int mutual_inf_cc_vec(const double *input_x, const double *input_y, const int n_elems,  const int k, double *mi) ;

        // case where input_y is discrete integers
        int mutual_inf_cd_vec(const double *input_x, const int *input_y, const int n_elems, const int k, double *mi ) ;

        //  But - there's only 1 algorithm for the conditional mutual information, so there's
        // only a need for 1 function for this.
        // Conditional mutual information.
       int cond_mutual_inf_vec(const double *input_x,  const double *input_y, const double *input_z, const int n_elems, const int k, double *mi) ;
}
#endif //REVEALER_MI_MATRIX_H
