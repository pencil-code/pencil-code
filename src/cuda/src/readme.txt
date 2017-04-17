##Initial drafts##

* All device modules should have postfix .cu or .cuh.

* All host modules should have postfix .cpp or .h (even though they contain C;
this because CUDA compiler generates essentially c++ code and we would need 
function prototype declarations for CUDA and C code to be compatible with each other)

* All CPU helper functions calling kernels should be named *_cuda, such as bouncond_cuda() so that
it's clear that the function doesn't do the job on the CPU.

* Constants in dconsts.cuh are implicitly static so they have a global scope, use them 
in kernels instead of passing nx,ny etc as function parameters. Device constants are defined as d_*NAME_IN_CAPS*.

* If some module containing only __device__ functions needs to be compiled separately, remember
to declare all constants needed from constants.cuh with extern classifier (extern __constant__ d_NX;)
In the worst case the program actually compiles and some very mysterious stuff starts happening.

* All modules that are separately compiled (that contain only __device__ functions) need to be
listed in the makefile

* unit_tests.cu can be used to check program validity (though it's very cumbersome to write new tests
and perhaps it's enough to check that everything looks ok from the slides)

* blockDim.x should always be a multiple of 32 and all threadIdx.x within a warp should access sequential 
indices in the computational domain

* #pragma once can be used as an include guard insead of header definitions and it looks nicer imo, but doesn't really matter and whatever floats your boat.




