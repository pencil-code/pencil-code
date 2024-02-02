## CMake settings
# V3.9 required for first-class CUDA support
# V3.17 required for the FindCUDAToolkit package
# V3.18 required for CMAKE_CUDA_ARCHITECTURES
cmake_minimum_required(VERSION 3.19)

## Options
option(DOUBLE_PRECISION   "Generates double precision code."                    ON)
option(BUILD_SAMPLES      "Builds projects in samples subdirectory."            ON)
option(BUILD_MHD_SAMPLES  "Builds MHD samples. Has no effect if BUILD_SAMPLES=OFF." ON)
option(MPI_ENABLED        "Enables additional functions for MPI communciation." ON)
option(USE_CUDA_AWARE_MPI "Uses GPUDirect RDMA for direct GPU-GPU communication instead of routing communication through host memory" ON)
option(MULTIGPU_ENABLED   "Enables multi-GPU on a single node. Uses peer-to-peer communication instead of MPI. Affects Legacy & Node layers only." ON)
option(VERBOSE            "Enables various status and warning messages"         OFF)
option(BUILD_UTILS        "Builds the utility library"                          ON)
option(SINGLEPASS_INTEGRATION "Perform integration in a single pass. Improves performance by ~20% but may introduce slightly larger floating-point arithmetic error than the conventional approach" OFF)
option(USE_DISTRIBUTED_IO "Use distributed IO (one file per process) instead of collective" ON)
option(USE_HIP     "Use HIP as the underlying GPGPU library instead of CUDA" OFF)
option(BUILD_SHARED_LIBS "Build Astaroth as a collection of shared libraries instead of statically built modules" OFF)
option(USE_OMP "Enables use of openmp" ON)
option(USE_PERFSTUBS "Build with perfstubs" OFF)
option(USE_VENDORED_PERFSTUBS "Use the vendored perfstubs in 3rd_party/perfstubs" ON)


if (USE_HIP)
    find_program(CMAKE_HIP_COMPILER hipcc REQUIRED)
    set(CMAKE_C_COMPILER ${CMAKE_HIP_COMPILER})
    set(CMAKE_CXX_COMPILER ${CMAKE_HIP_COMPILER})
endif()



## Project settings
project(astaroth C CXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

## Project-wide compilation flags
# -Werror TODO readd
#set(COMMON_FLAGS "-DOMPI_SKIP_MPICXX -Wall -Wextra -Wdouble-promotion -Wfloat-conversion -Wshadow -Wno-unused-result") #-DOMPI_SKIP_MPICXX is to force OpenMPI to use the C interface
set(COMMON_FLAGS "-DOMPI_SKIP_MPICXX -Wall -Wdouble-promotion -Wfloat-conversion -Wshadow -Wno-unused-result") #-DOMPI_SKIP_MPICXX is to force OpenMPI to use the C interface
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # For IDE highlighting (clangd)

## Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release") # Default
endif()
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

## Options
option(DOUBLE_PRECISION   "Generates double precision code."                        OFF)
option(BUILD_SAMPLES      "Builds projects in samples subdirectory."                ON)
option(BUILD_MHD_SAMPLES  "Builds MHD samples. Has no effect if BUILD_SAMPLES=OFF." ON)
option(MPI_ENABLED        "Enables additional functions for MPI communciation."     OFF)
option(USE_CUDA_AWARE_MPI "Uses GPUDirect RDMA for direct GPU-GPU communication instead of routing communication through host memory" ON)
option(MULTIGPU_ENABLED   "Enables multi-GPU on a single node. Uses peer-to-peer communication instead of MPI. Affects Legacy & Node layers only." ON)
option(VERBOSE            "Enables various status and warning messages"             OFF)
option(BUILD_UTILS        "Builds the utility library"                              ON)
option(BUILD_SHARED       "Creates shared library"                                  OFF)
option(PACKED_DATA_TRANSFERS    "Enables kernel for packed data transfers"          OFF)
option(SINGLEPASS_INTEGRATION "Perform integration in a single pass. Improves performance by ~20% but may introduce slightly larger floating-point arithmetic error than the conventional approach" ON)
option(ADJACENT_VERTEX_BUFFERS "Allocate vertex buffers such that they form a single massive" OFF)

## CUDA/HIP
if (USE_HIP)
    add_definitions(-DAC_USE_HIP=1)

    set(GPU_TARGETS "gfx908" CACHE STRING "GPU target")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --offload-arch=gfx908 --offload-arch=gfx90a")

    # Workaround for finding HIP on Triton
    list(APPEND CMAKE_PREFIX_PATH /opt/rocm/hip /opt/rocm)

    # Rocm
    #include_directories(/opt/rocm/include)
    #link_directories(/opt/rocm/lib)

    # Roctracer
    include_directories(/opt/rocm/roctracer/include)
    link_directories(/opt/rocm/roctracer/lib)

    # Rocrand
    include_directories(/opt/rocm/hiprand/include)
    include_directories(/opt/rocm/rocrand/include)
    link_directories(/opt/rocm/hiprand/lib)
    link_directories(/opt/rocm/rocrand/lib)

    find_package(hip)
else ()
    enable_language(CUDA)

    # Commented out CUDA host compiler flags: seem to be ignored anyways
    #string (REPLACE " " "," CUDA_COMMON_FLAGS "${COMMON_FLAGS}")
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -forward-unknown-to-host-compiler -Wall -Wextra -Wdouble-promotion -Wfloat-conversion -Wshadow --compiler-options=${CUDA_COMMON_FLAGS}")

    # To inspect ptx, run
    # cmake .. && make -j && nvcc -DAC_DOUBLE_PRECISION=1 --generate-code=arch=compute_60,code=[compute_60,sm_60] -c ../src/core/kernels/kernels.cu -I acc-runtime/api/ -I ../include/ -I ../acc-runtime/api/ -I ../acc-runtime -ptx -Xptxas=-Werror,--warn-on-double-precision-use
    # less kernels.ptx
    # Additional compilation info: --resource-usage

    # Not supported with older CUDA
    # set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --Werror all-warnings")
    if (NOT DOUBLE_PRECISION)
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xptxas=--warn-on-double-precision-use") # -Werror,
    endif()
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xptxas -dlcm=ca") # Cache global loads in L1
    #set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xptxas -dlcm=cg") # Cache global loads in L2 (default on >= NVIDIA Fermi)

    # Set device code architecture
    if (NOT CUDA_ARCHITECTURES)
        set(CMAKE_CUDA_ARCHITECTURES 60 61 70 80) # Default
    else ()
        set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCHITECTURES}) # User-specified
    endif()
endif()
if(USE_OMP)
find_package(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()
include(FetchContent)

## Global flags
if (DOUBLE_PRECISION)
    add_definitions(-DAC_DOUBLE_PRECISION=1)
    add_definitions(-DDOUBLE_PRECISION)
endif ()
if (VERBOSE)
    add_definitions(-DAC_VERBOSE=1)
endif ()
if (MPI_ENABLED)
    find_package(MPI REQUIRED COMPONENTS C)
    include_directories(${MPI_C_INCLUDE_DIRS})
    add_definitions(-DAC_MPI_ENABLED=1)
endif ()
if (USE_CUDA_AWARE_MPI)
    add_definitions(-DUSE_CUDA_AWARE_MPI=1)
endif()
if (MULTIGPU_ENABLED)
  add_definitions(-DAC_MULTIGPU_ENABLED=1)
endif()
if (SINGLEPASS_INTEGRATION)
    add_definitions(-DAC_SINGLEPASS_INTEGRATION)
endif()
if (USE_DISTRIBUTED_IO)
    add_definitions(-DUSE_DISTRIBUTED_IO=1)
endif()
if (USE_PERFSTUBS)
    add_definitions(-DUSE_PERFSTUBS=1)
    include(${PROJECT_SOURCE_DIR}/cmake/configure_perfstubs.cmake)
endif()

if (ADJACENT_VERTEX_BUFFERS)
    add_definitions(-DAC_ADJACENT_VERTEX_BUFFERS=1)
else ()
    add_definitions(-DAC_ADJACENT_VERTEX_BUFFERS=0)
endif ()

add_definitions(-DAC_DEFAULT_CONFIG="${CMAKE_SOURCE_DIR}/config/astaroth.conf")

## Includes
include_directories(./submodule/include)                                  # Library headers
include_directories(./submodule/build)
include_directories(./submodule/build/acc-runtime/api)
include_directories(./..)
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}) # CUDA headers
include_directories(./NUMA-scheduler/include)
include_directories(./submodule/build/src/core)
include_directories(./submodule/build/src/core/kernels)
include_directories(./submodule/build/src/utils)

link_directories(./submodule/build/src/core)
link_directories(./submodule/build/src/core/kernels)
link_directories(./submodule/build/src/utils)

add_subdirectory(NUMA-scheduler)
add_compile_options(-DMODPRE=${MODULE_PREFIX} -DMODIN=${MODULE_INFIX} -DMODSUF=${MODULE_SUFFIX})
add_compile_options(-DPACKED_DATA_TRANSFERS=1)
add_compile_options(-DAC_ADJACENT_VERTEX_BUFFERS=1)
add_definitions(-DAC_DOUBLE_PRECISION=1)
add_library(astaroth SHARED gpu_astaroth.cc)
target_link_libraries(astaroth astaroth_core kernels astaroth_utils mt)
