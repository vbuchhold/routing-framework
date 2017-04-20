#pragma once

// Hints to the compiler that the specified condition is expected to be true.
#ifdef __GNUC__ // GNU compiler.
# define LIKELY(cond) __builtin_expect(!!(cond), 1)
#endif
#ifndef LIKELY  // Default definition.
# define LIKELY(cond) (cond)
#endif

// Hints to the compiler that the specified condition is expected to be false.
#ifdef __GNUC__  // GNU compiler.
# define UNLIKELY(cond) __builtin_expect(!!(cond), 0)
#endif
#ifndef UNLIKELY // Default definition.
# define UNLIKELY(cond) (cond)
#endif

// The weakest alignment a type shall have for use with SIMD instructions.
#if defined __AVX2__
constexpr int MIN_ALIGNMENT = 32;
#elif defined __SSE4_2__
constexpr int MIN_ALIGNMENT = 16;
#else
constexpr int MIN_ALIGNMENT = 1;
#endif
