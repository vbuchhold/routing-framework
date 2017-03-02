#pragma once

#include <nanoflann.hpp>

namespace kdtree {

// This class tells a kd-tree to use the Manhattan metric to compute distances.
template <typename PointSetAdapterT>
using ManhattanMetric = nanoflann::L1_Adaptor<int, PointSetAdapterT, int>;

}
