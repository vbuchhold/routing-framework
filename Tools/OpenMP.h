#pragma once

// On platforms supporting OpenMP, this file includes the omp.h header file. On other platforms,
// this file provides stubs for the runtime library routines defined in the OpenMP API.

#ifdef _OPENMP

#include <omp.h>

#else

int omp_get_num_threads(void)
{
  return 1;
}

int omp_get_thread_num(void)
{
  return 0;
}

#endif
