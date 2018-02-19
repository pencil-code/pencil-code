#
# CMakeLists.txt for generating the makefile for Astaroth.
#   Usage: mkdir build && cd build && cmake <optional flags> ..
#   
#   For example: cmake -DUSE_DOUBLE_PRECISION=ON ..
# 
#   If you want to see the exact flags used during compilation, run
#   "make -j VERBOSE=1"

#-------------------General---------------------------------------------------#

project (ASTAROTH)
cmake_minimum_required (VERSION 3.5)
cmake_policy (SET CMP0023 NEW)


#-------------------Set user options with default values---------------------#

#Usage f.ex. cmake -DBUILD_DEBUG=ON ..
option(BUILD_DEBUG "Builds the program with debug symbols and extensive error checking" OFF)
option(BUILD_STANDALONE "Builds standalone Astaroth" OFF) #Changed default to off for PC
option(CUDA_BUILD_LEGACY "Builds GPU code for older GPUs (Fermi)" OFF)
option(USE_DOUBLE_PRECISION "Generates double precision code" OFF)


#-------------------Determine build type--------------------------------------#

#Available types (case-sensitive):
#RELEASE         (best performance)
#DEBUG           (w/ debug information, non-concurrent kernels)
if (BUILD_DEBUG)
    set(CMAKE_BUILD_TYPE DEBUG)
else ()
    set(CMAKE_BUILD_TYPE RELEASE)
endif()
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})


#----------------------Find packages------------------------------------------#

#C/C++
# the compiler used for C files 
message( STATUS "CMAKE_C_COMPILER: " ${CMAKE_C_COMPILER} )
message( STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER} )
message( STATUS "CMAKE_CXX_COMPILER_ID: " ${CMAKE_CXX_COMPILER_ID} )

#CUDA
find_package(CUDA)
if (NOT CUDA_FOUND)
    #find_package(CUDA REQUIRED) gives a confusing error message if it fails, 
    #therefore we print the reason here explicitly
    message(FATAL_ERROR "CUDA not found")
endif()

#OpenMP
find_package(OpenMP)
if (NOT OPENMP_FOUND)
    message(WARNING "OpenMP not found. All host-side concurrency disabled \
                    (lower performance).")
else ()
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")    
endif()


#----------------------Compilation settings-----------------------------------#

#Debug and verification
#set(CMAKE_VERBOSE_MAKEFILE OFF)
#set(CXX_VERBOSE_BUILD OFF)
#set(CUDA_VERBOSE_BUILD OFF)
#include(CTest)
#add_test(ac_test ac_run)
#find_program(MEMORYCHECK_COMMAND valgrind)
#set(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full" )

#CUDA settings
set(CUDA_SEPARABLE_COMPILATION ON)
set(CUDA_PROPAGATE_HOST_FLAGS ON)
#set(CUDA_BUILD_CUBIN ON) #Requires that we're compiling for only one architecture


#----------------------Setup defines------------------------------------------#

if (USE_DOUBLE_PRECISION)
	add_definitions(-DUSE_DOUBLE_PRECISION)
endif()
#add_definitions(-DUSE_DOUBLE_PRECISION)
#add_definitions(-DSRC_DIR=${CMAKE_SOURCE_DIR})


#----------------------Setup CXX compilation flags----------------------------#

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")

    set(CXX_WARNING_FLAGS "-Wall -Wextra -Werror -Wno-error=unused-parameter\
                           -Wno-error=unused-function -Wno-error=unknown-pragmas -Wno-error=strict-overflow")# -Wdouble-promotion -Wfloat-conversion")

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}\
                             -O3 -march=native -ffast-math -funroll-loops -D_FORCE_INLINES")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}\
                             -O0 -g")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARNING_FLAGS}")

elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    #Allow Intel-compiler specific flags

    set(CXX_WARNING_FLAGS "-Wall")

    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}\
                             -O3 -use_fast_math -D_FORCE_INLINES")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}\
                             -O0 -g")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARNING_FLAGS}")

endif()

message("CXX_FLAGS: " ${CMAKE_CXX_FLAGS})


#----------------------Setup CUDA compilation flags----------------------------#

#Generate code for the default architecture (Kepler)
set(CUDA_ARCH_FLAGS -gencode arch=compute_35,code=sm_35;
                    -gencode arch=compute_37,code=sm_37)

#Enable support for Pascal and Volta GPUs
if (CUDA_VERSION_MAJOR GREATER 8 OR CUDA_VERSION_MAJOR EQUAL 8)
    set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS};
                        -gencode arch=compute_60,code=sm_60)
    if (CUDA_VERSION_MAJOR GREATER 9 OR CUDA_VERSION_MAJOR EQUAL 9)
        set(CUDA_ARCH_FLAGS ${CUDA_ARCH_FLAGS};
                        -gencode arch=compute_70,code=sm_70)
    else()
        message(WARNING "CUDA 9.0 or greater not available, \
                         cannot generate code for cc 7.0 GPUs")
    endif()
else()
    message(WARNING "CUDA 8.0 or greater not available, \
                     cannot generate code for cc 6.0 GPUs")
endif()

#Generate code for older GPUs
if (CUDA_BUILD_LEGACY)
    set(CUDA_ARCH_FLAGS -gencode arch=compute_20,code=sm_20)
endif()

#Additional CUDA optimization flags
if (CMAKE_BUILD_TYPE MATCHES RELEASE)
    #Doesn't set any additional flags, see CUDA_NVCC_FLAGS_DEBUG how to add more
    set(CUDA_NVCC_FLAGS_RELEASE ${CUDA_NVCC_FLAGS_RELEASE}) 
endif()

#Additional CUDA debug flags
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    #The debug flags must be set inside this if clause, since either CMake 3.5
    #or nvcc 7.5 is bugged:
    #CMake converts these into empty strings when doing RELEASE build, but nvcc
    #7.5 fails to parse empty flags.
    set(CUDA_NVCC_FLAGS_DEBUG ${CUDA_NVCC_FLAGS_DEBUG};
                               --device-debug; 
                               --generate-line-info;
                               --ptxas-options=-v)
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};${CUDA_ARCH_FLAGS}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    #NOTE: nvcc does not support -std=c++11 with intel compilers :( 
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};${CUDA_ARCH_FLAGS}")
endif()

message("CUDA_NVCC_FLAGS: " ${CUDA_NVCC_FLAGS})


#----------------------Setup core subdirectories------------------------------#

#Include root directory (.) so that the following modules can include their
#parent dir (f.ex. #include "common/stuff.h" instead of "../common/stuff")
include_directories(.)
add_subdirectory (common)
add_subdirectory (gpu)

#CUDA sources
set(CUDA_SRC_DIR gpu/cuda)
set(CUDA_CORE ${CUDA_SRC_DIR}/core/cuda_core.cu)
set(CUDA_GENERIC
     ${CUDA_SRC_DIR}/generic/slice_cuda_generic.cu 
     ${CUDA_SRC_DIR}/generic/boundcond_cuda_generic.cu  
     ${CUDA_SRC_DIR}/generic/rk3_cuda_generic.cu
     ${CUDA_SRC_DIR}/generic/collectiveops_cuda_generic.cu
     ${CUDA_SRC_DIR}/cuda_generic.cu)
set(CUDA_MODULES ${CUDA_CORE} ${CUDA_GENERIC})


#----------------------Link---------------------------------------------------#

#Create shared library of the GPU module (astaroth_core)
#Use -fPIC if -fpic not supported. Some quick non-scientific tests:
#Without fpic: 4.94 user, 4.04 system, 0:09.88 elapsed
#With fpic: 4.96 user, 4.02 system, 0:09.90 elapsed
#With fPIC: 4.94 user, 4.05 system, 0:10.23 elapsed
CUDA_ADD_LIBRARY(astaroth_core SHARED gpu ${CUDA_MODULES} common OPTIONS --compiler-options "-fpic")

if (BUILD_STANDALONE)
    #Add additional subdirectories
    add_subdirectory (utils)
    add_subdirectory (cpu/model)
    #Create executable and link
    cuda_add_executable(ac_run main.cc)
    target_link_libraries(ac_run utils gpu common cpu_model astaroth_core)
endif()




























