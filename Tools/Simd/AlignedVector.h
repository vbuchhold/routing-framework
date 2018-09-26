#pragma once

#include <vector>

#include <boost/align/aligned_allocator.hpp>

#include "Tools/CompilerSpecific.h"
#include "Tools/MachineSpecs.h"

// A std::vector whose elements are properly aligned for use with SIMD instructions.
template <typename T>
using AlignedVector = std::vector<
    T, boost::alignment::aligned_allocator<T, std::max(MIN_ALIGNMENT, CACHE_LINE_SIZE)>>;
