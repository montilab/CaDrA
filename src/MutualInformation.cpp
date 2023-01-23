//
// Created by bgregor on 11/1/22.
//

#include "MutualInformation.h"
#include <cmath>
#include <array>
#include "nanoflann.hpp"
#include <iostream>
#include <numeric>
#include <functional>
#include <list>
#include <set>
#include <map>

#ifndef BOOST_DIGAMMA
#include <Rmath.h>
#include "R.h"
#include "Rinternals.h"
#else
#include <boost/math/special_functions/digamma.hpp>
#endif

namespace CaDrA {

    double MutualInformation::mutual_information_cc(const ArrayXd &x, const ArrayXd& y) {
        // compute mutual information based on k-nearest neighbors
        // Creat a KD-tree using nanoflann and its Eigen matrix adapter
        // Store the size of the vectors for easy reference
        auto N = x.size() ;

        Array2col  tmp_mat(N, 2) ;
        tmp_mat.col(0) = scale(x) ;
        tmp_mat.col(1) = scale(y) ;
        // Pointers to those scaled values so that no extra copies are made.
        //const double *x_scale =  tmp_mat.col(0).data();
        //const double *y_scale =  tmp_mat.col(1).data() ;
        // Map the double array pointer to an Eigen vector without a copy.
        MapArrayConst x_scale(tmp_mat.col(0).data(), N) ;
        MapArrayConst y_scale(tmp_mat.col(1).data(), N) ;

        vector<double> dists  = calc_distances(N, tmp_mat).first;

        double x_digamma_sum = sum_digamma_from_neighbors(x_scale, dists) ;
        double y_digamma_sum = sum_digamma_from_neighbors(y_scale, dists) ;

        // mutual info computation
        double mi = MutualInformation::digamma_f(N)
                    + MutualInformation::digamma_f(m_k)
                    - (x_digamma_sum + y_digamma_sum) / N;

        return std::max(0.0,mi) ;
    }

    double MutualInformation::sum_digamma_from_neighbors(MapArrayConst &vec, const vector<double> &dists) {
        // This one is called from mutual_information_cc and cond_mutual_information for the neighbors
        // for a single vector.
        size_t N = dists.size() ;
        double sum = 0.0 ;

        // KD-Tree for this vector
        nanoflann::KDTreeEigenMatrixAdaptor<MapArrayConst,-1,metric_Chebyshev> vec_tree(1, vec, 10) ;

        std::vector<std::pair<long, double>> ret_matches;
        for (size_t i = 0 ; i < N ; ++i) {
            double pt = vec(i) ; // avoids type issues with the compiler and the radiusSearch.
            double tmp = vec_tree.index->radiusSearch(&pt, dists[i] , ret_matches , nanoflann::SearchParams(10));
            sum += MutualInformation::digamma_f(tmp) ;
            ret_matches.clear() ;
        }
        return sum ;
    }

    double MutualInformation::sum_digamma_from_neighbors(MapArrayConst &vec1, MapArrayConst &vec2, const vector<double> &dists) {
        // Sum of digamma_f functions over neighbor counts for 2D.
        size_t N = dists.size() ;
        double sum = 0.0 ;

        // KD-Tree for this vector
        Array2col tmp_mat(N, 2) ;
        tmp_mat.col(0) = vec1 ;
        tmp_mat.col(1) = vec2 ;
        nanoflann::KDTreeEigenMatrixAdaptor<Array2col ,-1,metric_Chebyshev> vec_tree(2, tmp_mat, 10) ;

        std::vector<std::pair<long, double>> ret_matches;
        array<double,2> pt ;
        for (size_t i = 0 ; i < N ; ++i) {
            pt[0] = tmp_mat(i,0) ;
            pt[1] = tmp_mat(i,1) ;
            double tmp = vec_tree.index->radiusSearch(pt.data(), dists[i] , ret_matches , nanoflann::SearchParams(10));
            sum += MutualInformation::digamma_f(tmp) ;
            ret_matches.clear() ;
        }
        return sum ;
    }



    ArrayXd MutualInformation::scale(const ArrayXd &x) const {// Center and scale the x data
        auto size_v = x.size() ;
        double mean_x = x.mean() ;
        auto x_p = (x - mean_x) ;
        double std_dev = std::sqrt((x_p).square().sum()/size_v);
        ArrayXd x_scale = x_p  / std_dev ;
        return x_scale;
    }

    double MutualInformation::cond_mutual_information_ccc(const ArrayXd &x, const ArrayXd& y, const ArrayXd& z) {
        // Compute the conditional mutual information.
        // CMI = psi(m_k) - mean( psi(n_xz[i] + psy(n_yz[i])-psi(n_z[i]) ) for all elements i
        auto N = x.size() ;

        // Get the 3-dimensional distance, and then re-use that for the xz, yz, and calculations, so it
        // proceeds pretty much exactly as before.

        Array3col tmp_mat(N, 3) ;

        tmp_mat.col(0) = scale(x) ;
        tmp_mat.col(1) = scale(y) ;
        tmp_mat.col(2) = scale(z) ;

        // Map the double array pointer to an Eigen vector without a copy.
        MapArrayConst x_scale(tmp_mat.col(0).data(), N) ;
        MapArrayConst y_scale(tmp_mat.col(1).data(), N) ;
        MapArrayConst z_scale(tmp_mat.col(2).data(), N) ;

        // Calculating the distances also calculates the number of neighbors, so
        // calc_distances returns both.
        vector<double> dists = calc_distances(N, tmp_mat).first ;

        // Get the digamma_f values...
        double xz_digamma_sum = sum_digamma_from_neighbors(x_scale, z_scale, dists) ;
        double yz_digamma_sum = sum_digamma_from_neighbors(y_scale, z_scale, dists) ;
        double z_digamma_sum = sum_digamma_from_neighbors(z_scale, dists) ;

        // mutual info computation
        double mi = MutualInformation::digamma_f(m_k)
                    - (xz_digamma_sum + yz_digamma_sum - z_digamma_sum) / N;

        return std::max(0.0,mi) ;
    }

    pair<vector<double>,vector<size_t>>  MutualInformation::calc_distances(const size_t N, const Array<double, -1, 3> &tmp_mat) const {
        // Calculate the Chebyshev distances and numbers of neighbors for the 2D array tmp_mat.
        kd_tree_3d mat_index(tmp_mat.cols(),std::cref(tmp_mat),10) ;
        // We want N neighbors in addition to the point itself so
        // add 1 to the # of neighbors.
        int real_k = m_k + 1  ;

        // Chebyshev distance
        vector<double> dists(N) ;
        // Number of neighbors
        vector<size_t> neighbors(N) ;

        // a query point.
        array<double,3> query_pt ;
        for (size_t i = 0 ; i < N ; ++i) {
            // store indexes and distances
            vector<long> ret_indexes(real_k, 0.0);
            vector<double> out_dists_sqr(real_k,0.0);

            for (size_t j = 0 ; j < 3 ; ++j)
                query_pt[j] = tmp_mat(i, j);

            neighbors[i] = mat_index.index->knnSearch(&query_pt[0], real_k,
                                                      &ret_indexes[0], &out_dists_sqr[0]) ;
            out_dists_sqr.resize(neighbors[i]) ;
            dists[i] = out_dists_sqr.back() ;
        }
        return make_pair(dists,neighbors) ;
    }

    double MutualInformation::cond_mutual_information_cdd(const ArrayXd &x, const ArrayXi& y, const ArrayXi& z) {
        // Implement conditional mutual information between continuous x and discrete y & z.
        // Quick implementation: convert the y & z arrays over to double precision, call the
        // cond_mutual_information_ccc function.

      return cond_mutual_information_ccc(x, y.cast<double>(), z.cast<double>()) ;
    }


    double MutualInformation::mutual_information_cd(const ArrayXd &x, const ArrayXi &y) {
        // Compute the mutual information for a continuous vector x and a
        // discrete vector y.
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
            if (count > 0) {
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
                    vector<long> ret_indexes(real_k, 0.0);
                    vector<double> out_dists(real_k, 0.0);
                    int idx = label_indices[i];
                    double query_pt[1] = {masked_x[i]};
                    // out_dists stores the distances for this label. Get the max one.
                    auto neighbors = label_index_tree.index->knnSearch(&query_pt[0], real_k,
                                                                       &ret_indexes[0], &out_dists[0]);
                    out_dists.resize(neighbors) ;
                    // The last one is out_dists is the furthest distance.
                    auto max_dist = out_dists.back() ;

                    std::vector<std::pair<long, double>> ret_matches;
                    m_all.push_back(xscale_index_tree.index->radiusSearch(query_pt, max_dist, ret_matches , nanoflann::SearchParams(10))) ;
                }
            }
        }
        double N_mod = all_indices.size() ;
        double digamma_N = digamma_f(N) ;
        // Calculate the mean of digammas over the count of samples for each label.
        double digamma_labels = std::accumulate(label_counts.begin(),label_counts.end(),0.0,
                                                [N_mod](double m, double n){ return m +
                                                        MutualInformation::digamma_f(n) * n / N_mod ; } ) ;
        // Get the same for the k's used at each label.
        double digamma_k = std::accumulate(k_all.begin(),k_all.end(),0.0,
                                           [N_mod](double m, vector<double> &n){ return m +
                                                   MutualInformation::digamma_f(std::max(n[1] - 1.0, 1.0)) * n[0] / N_mod; } ) ;
        double sum_m = std::accumulate(m_all.begin(),m_all.end(),0.0) ;
        double digamma_m = std::accumulate(m_all.begin(),m_all.end(),
                                                   0.0,  [](double m, double n){ return m +
                        MutualInformation::digamma_f(n); });
        digamma_m = digamma_m / N_mod ;

        // mutual info computation
        // Matlab:   4.602 - 4.3965 + 0.9228 - 1.0961
        double mi = digamma_N - digamma_labels + digamma_k - digamma_m ;
        return mi;
    }
    // Construct this object with the neighborhood size k.
    MutualInformation::MutualInformation(const int mK) : m_k(mK) {}

    int MutualInformation::get_k() const {
        return m_k;
    }

    void MutualInformation::set_k(int mK) {
        m_k = mK;
    }

    pair<vector<double>,vector<size_t>> MutualInformation::calc_distances(const size_t N, const Array<double, -1, 2> &tmp_mat) const {
        // Calculate the Chebyshev distances and numbers of neighbors for the 2D array tmp_mat.
        kd_tree_2d mat_index(tmp_mat.cols(),std::cref(tmp_mat),10) ;
        // We want N neighbors in addition to the point itself so
        // add 1 to the # of neighbors.
        int real_k = m_k + 1  ;

        // Chebyshev distance
        vector<double> dists(N) ;
        // Number of neighbors
        vector<size_t> neighbors(N) ;

        // a query point.
        array<double,2> query_pt ;
        for (size_t i = 0 ; i < N ; ++i) {
            // store indexes and distances
            vector<long> ret_indexes(real_k, 0.0);
            vector<double> out_dists_sqr(real_k,0.0);

            query_pt[0] = tmp_mat(i, 0);
            query_pt[1] = tmp_mat(i, 1);

            neighbors[i] = mat_index.index->knnSearch(&query_pt[0], real_k,
                                                      &ret_indexes[0], &out_dists_sqr[0]) ;
            out_dists_sqr.resize(neighbors[i]) ;
            auto max_dist= out_dists_sqr.back() ;//+ std::numeric_limits<double>::epsilon();
            dists[i] = max_dist ; //abs(x_scale(idx) - x_scale(i));
        }
        return make_pair(dists,neighbors) ;
    }

    // For development use the boost digamma_f funtion. For
    // R production this can be switched to call the digamma_f
    // function that comes with R.
    inline double MutualInformation::digamma_f(const double x) {
#ifndef BOOST_DIGAMMA
        // This is the R math library digamma function.
        return digamma(x) ;
#else
        // for development/testing outside of R,
        // use the Boost library's digamma function.
        return boost::math::digamma(x) ;
#endif
    }
    
} // CaDrA