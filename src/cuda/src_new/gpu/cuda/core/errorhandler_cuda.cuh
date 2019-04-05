/*
*   Errorhandling for CUDA
*/
#pragma once
#include "common/errorhandler.h"


#define CUDA_ERRCHK_ALWAYS(ans) { cuda_assert((ans), __FILE__, __LINE__); }
inline void cuda_assert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#define CUDA_ERRCHK_KERNEL_ALWAYS() { CUDA_ERRCHK_ALWAYS(cudaPeekAtLastError()); CUDA_ERRCHK_ALWAYS(cudaDeviceSynchronize()); }


//NDEBUG is defined implicitly at compile-time for RELEASE build
//(i.e. cmake is invoked without flags or -DCMAKE_BUILD_TYPE=RELEASE
//
//cmake -DCMAKE_BUILD_TYPE=DEBUG .. enables all error checking but disables all
//host-side concurrency and makes the GPU function calls to execute sequentially
#ifdef NDEBUG 
    #define CUDA_ERRCHK(ans) { ans; }
    #define CUDA_ERRCHK_KERNEL() {}
#else
    #define CUDA_ERRCHK(ans) { CUDA_ERRCHK_ALWAYS(ans); }
    #define CUDA_ERRCHK_KERNEL() { CUDA_ERRCHK(cudaPeekAtLastError()); CUDA_ERRCHK(cudaDeviceSynchronize()); }
#endif
