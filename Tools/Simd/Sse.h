#pragma once

#include <cassert>

#include <emmintrin.h>
#include <smmintrin.h>

// A wrapper around the intrinsic _mm_extract_epi32. Some compilers require its second argument to
// be a constant expression. This wrapper does exactly the same thing as _mm_extract_epi32, but
// allows the second argument to be a non-constant expression.
int mm_extract_epi32(__m128i a, const int index) {
  switch (index) {
    case 0:
      return _mm_extract_epi32(a, 0);
    case 1:
      return _mm_extract_epi32(a, 1);
    case 2:
      return _mm_extract_epi32(a, 2);
    case 3:
      return _mm_extract_epi32(a, 3);
    default:
      assert(false);
      return 0;
  }
}

// A wrapper around the intrinsic _mm_insert_epi32. Some compilers require its third argument to be
// a constant expression. This wrapper does exactly the same thing as _mm_insert_epi32, but allows
// the third argument to be a non-constant expression.
__m128i mm_insert_epi32(__m128i a, int i, const int index) {
  switch (index) {
    case 0:
      return _mm_insert_epi32(a, i, 0);
    case 1:
      return _mm_insert_epi32(a, i, 1);
    case 2:
      return _mm_insert_epi32(a, i, 2);
    case 3:
      return _mm_insert_epi32(a, i, 3);
    default:
      assert(false);
      return __m128i();
  }
}
