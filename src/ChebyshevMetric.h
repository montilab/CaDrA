//
// Created by bgregor on 10/3/22.
//

#ifndef REVEALER_CHEBYSHEVMETRIC_H
#define REVEALER_CHEBYSHEVMETRIC_H

#include "nanoflann.hpp"

#include <vector>

namespace CaDrA {
/** Chebyshev distance metric
 *
 *  T Type of the elements (e.g. double, float, uint8_t)
 *  DataSource Source of the data, i.e. where the vectors are stored
 *  _DistanceType Type of distance variables (must be signed)
 *  AccessorType Type of the arguments with which the data can be
 * accessed (e.g. float, double, int64_t, T*)
 */
    template<
            class T, class DataSource, typename _DistanceType = T,
            typename AccessorType = uint32_t>
    struct Chebyshev_Adaptor {
        // This gets marked as unused here but it gets used when the
        // template is instantiated by nanoflann. Don't remove it.
        using ElementType = T;
        using DistanceType = _DistanceType;

        const DataSource &data_source;

        Chebyshev_Adaptor(const DataSource &_data_source)
                : data_source(_data_source) {
        }

        inline DistanceType evalMetric(
                const T *a, const AccessorType b_idx, size_t size) const {
            // Chebyshev distance: max difference between points in a and b
            DistanceType max_diff = std::numeric_limits<DistanceType>::min();
            for (size_t i = 0; i < size; ++i) {
                double diff = std::abs(a[i] - data_source.kdtree_get_pt(b_idx, i));
                if (diff > max_diff) {
                    max_diff = diff ;
                }
            }
           return max_diff ;
        }

        template<typename U, typename V>
        inline DistanceType accum_dist(const U a, const V b, const size_t) const {
            return std::abs(a - b);
        }
    };


    struct metric_Chebyshev : public nanoflann::Metric
    {
        template <class T, class DataSource, typename AccessorType = uint32_t>
        struct traits
        {
            using distance_t = Chebyshev_Adaptor<T, DataSource, T, AccessorType>;
        };
    };


}

#endif //REVEALER_CHEBYSHEVMETRIC_H
