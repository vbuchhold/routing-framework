#pragma once

#include <cassert>

#include <immintrin.h>

// A wrapper around the intrinsic _mm256_extract_epi32. Some compilers require its second argument
// to be a constant expression. This wrapper does exactly the same thing as _mm256_extract_epi32,
// but allows the second argument to be a non-constant expression.
int mm256_extract_epi32(__m256i a, const int index) {
  switch (index) {
    case 0:
      return _mm256_extract_epi32(a, 0);
    case 1:
      return _mm256_extract_epi32(a, 1);
    case 2:
      return _mm256_extract_epi32(a, 2);
    case 3:
      return _mm256_extract_epi32(a, 3);
    case 4:
      return _mm256_extract_epi32(a, 4);
    case 5:
      return _mm256_extract_epi32(a, 5);
    case 6:
      return _mm256_extract_epi32(a, 6);
    case 7:
      return _mm256_extract_epi32(a, 7);
    default:
      assert(false);
      return 0;
  }
}

// A wrapper around the intrinsic _mm256_insert_epi32. Some compilers require its third argument to
// be a constant expression. This wrapper does exactly the same thing as _mm256_insert_epi32, but
// allows the third argument to be a non-constant expression.
__m256i mm256_insert_epi32(__m256i a, int i, const int index) {
  switch (index) {
    case 0:
      return _mm256_insert_epi32(a, i, 0);
    case 1:
      return _mm256_insert_epi32(a, i, 1);
    case 2:
      return _mm256_insert_epi32(a, i, 2);
    case 3:
      return _mm256_insert_epi32(a, i, 3);
    case 4:
      return _mm256_insert_epi32(a, i, 4);
    case 5:
      return _mm256_insert_epi32(a, i, 5);
    case 6:
      return _mm256_insert_epi32(a, i, 6);
    case 7:
      return _mm256_insert_epi32(a, i, 7);
    default:
      assert(false);
      return __m256i();
  }
}
