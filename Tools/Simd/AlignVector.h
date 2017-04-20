#pragma once

#include <vector>

#include <boost/align/aligned_allocator.hpp>

#include "Tools/CompilerSpecific.h"

// A std::vector whose elements are properly aligned for use with SIMD instructions.
template <typename T>
using AlignVector = std::vector<T, boost::alignment::aligned_allocator<T, MIN_ALIGNMENT>>;
