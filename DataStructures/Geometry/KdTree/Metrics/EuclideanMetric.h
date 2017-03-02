#pragma once

#include <cstdint>

#include <nanoflann.hpp>

namespace kdtree {

// This class tells a kd-tree to use the Euclidean metric to compute distances.
template <typename PointSetAdapterT>
using EuclideanMetric = nanoflann::L2_Simple_Adaptor<int, PointSetAdapterT, int64_t>;

}
