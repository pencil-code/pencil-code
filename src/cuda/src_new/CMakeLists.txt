


#General
project (ASTAROTH_CODE)
cmake_minimum_required (VERSION 2.8)
cmake_policy (SET CMP0023 NEW)
set (CMAKE_CXX_STANDARD 11)

#Debug and verification
set (CMAKE_VERBOSE_MAKEFILE OFF)
set (CUDA_VERBOSE_BUILD OFF)

#Options
option (USE_CUDA "Determines whether we compile the CUDA code" ON)
find_package(CUDA REQUIRED)
set (CUDA_SEPARABLE_COMPILATION ON)
set (CUDA_PROPAGATE_HOST_FLAGS ON)

#Subdirectories
include_directories(.)#So that the modules can include their parent dir (f.ex. #include "common/stuff.h" instead of "../common/stuff")
add_subdirectory (common)

#CUDA flags
set (WARNING_FLAGS "-Wall,-Wextra,-Werror,-Wno-error=unused-parameter,-Wno-error=unknown-pragmas")
set (DEBUG_FLAGS "-g")
set (GENERIC_FLAGS "-std=c++11")
#set (ARCH_FLAGS -gencode arch=compute_20,code=sm_20)#Should not be hardcoded
#set (OTHER_FLAGS "-lineinfo,--ptxas-options=-v")#Any other CUDA flags, delimited by ;

set (CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};--compiler-options;${WARNING_FLAGS},${DEBUG_FLAGS};${GENERIC_FLAGS};${ARCH_FLAGS};${OTHER_FLAGS})


#Host code flags
set (CXX_WARNING_FLAGS "")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARNING_FLAGS} ${DEBUG_FLAGS}")


#CUDA sources (for separate compilation defined here since CMake refuses to play along otherwise)
set (CUDA_SRC_DIR gpu/cuda)

#Implementations
set (CUDA_CORE ${CUDA_SRC_DIR}/core/cuda_core.cu)

set (CUDA_GENERIC
     ${CUDA_SRC_DIR}/generic/slice_cuda_generic.cu 
     ${CUDA_SRC_DIR}/generic/boundcond_cuda_generic.cu  
     ${CUDA_SRC_DIR}/generic/rk3_cuda_generic.cu
     ${CUDA_SRC_DIR}/generic/collectiveops_cuda_generic.cu
     ${CUDA_SRC_DIR}/cuda_generic.cu)

set (CUDA_MODULES ${CUDA_CORE} ${CUDA_GENERIC})

CUDA_ADD_LIBRARY (cuda_modules ${CUDA_MODULES})

#Definitions
#add_definitions(-DSRC_DIR=${CMAKE_SOURCE_DIR})

#Link libraries
cuda_add_executable(ac_run gpu_astaroth_new.cc)
target_link_libraries (ac_run common cuda_modules)









































