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
