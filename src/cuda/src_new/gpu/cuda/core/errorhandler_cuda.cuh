/*
*   Errorhandling for CUDA
*/
#pragma once

//#define NO_GPU_ERRCHK

#ifdef NO_GPU_ERRCHK
    #define CUDA_ERRCHK(ans) { ans; }
    #define CUDA_ERRCHK_KERNEL() {}

    #define NDEBUG //Disable asserts
    #include "common/errorhandler.h"
#else
    #include "common/errorhandler.h"
    #define CUDA_ERRCHK(ans) { cuda_assert((ans), __FILE__, __LINE__); }
    inline void cuda_assert(cudaError_t code, const char *file, int line, bool abort=true)
    {
       if (code != cudaSuccess) {
          fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
          if (abort) exit(code);
       }
    }

    #define CUDA_ERRCHK_KERNEL() { CUDA_ERRCHK(cudaPeekAtLastError()); CUDA_ERRCHK(cudaDeviceSynchronize()); }
#endif
