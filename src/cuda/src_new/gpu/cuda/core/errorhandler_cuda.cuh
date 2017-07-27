/*
*   Errorhandling for CUDA
*/
#pragma once
#include <stdio.h>

#define CUDA_ERRCHK(ans) { cuda_assert((ans), __FILE__, __LINE__); }
inline void cuda_assert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define CUDA_ERRCHK_KERNEL() { CUDA_ERRCHK(cudaPeekAtLastError()); CUDA_ERRCHK(cudaDeviceSynchronize()); }
