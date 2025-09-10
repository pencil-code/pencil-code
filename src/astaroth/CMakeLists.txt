## CMake settings
# V3.9 required for first-class CUDA support
# V3.17 required for the FindCUDAToolkit package
# V3.18 required for CMAKE_CUDA_ARCHITECTURES
cmake_minimum_required(VERSION 3.19)

## Options
option(DOUBLE_PRECISION       "Generates double precision code."                           ON)
option(USE_HIP                "Use HIP as the underlying GPGPU library instead of CUDA"    OFF)
option(USE_OMP                "Enables use of openmp"                                      ON)
option(USE_DARDEL             "Target machine is Dardel(PDC, Stockholm"                    OFF)
option(PACKED_DATA_TRANSFERS  "Enables kernel for packed data transfers"                   OFF)
option(SINGLEPASS_INTEGRATION "Enables single pass integration"                            OFF)
option(RUNTIME_COMPILATION    "Enables runtime compilation"                                OFF)
option(TRANSPILATION          "Tmp flag to make interface compatible with transpiled code" OFF)
option(TRAINING	              "Using torchfort"						   OFF)
option(CPU_BUILD              "CPU-only build"                                             OFF)

## Project settings
project(astaroth C CXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

## Project-wide compilation flags
# -Werror TODO readd
#set(COMMON_FLAGS "-DOMPI_SKIP_MPICXX -Wall -Wextra -Wdouble-promotion -Wfloat-conversion -Wshadow -Wno-unused-result") #-DOMPI_SKIP_MPICXX is to force OpenMPI to use the C interface
#set(COMMON_FLAGS "-DOMPI_SKIP_MPICXX -Wall -Wno-unused-variable -Wno-unused-parameter -Wno-unused-function -Wno-#pragma-messages -Wno-extra-tokens -Wdouble-promotion -Wfloat-conversion -Wshadow -Wno-unused-result") #-DOMPI_SKIP_MPICXX is to force OpenMPI to use the C interface
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${COMMON_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COMMON_FLAGS}")
set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # For IDE highlighting (clangd)

if (NOT CPU_BUILD AND DEFINED ENV{AC_CPU_BUILD})
	set(CPU_BUILD $ENV{AC_CPU_BUILD})
	message("Set CPU_BUILD=${CPU_BUILD} from ENV value")
endif()

if(CPU_BUILD)
	add_definitions(-DAC_CPU_BUILD=1)
else()
	add_definitions(-DAC_CPU_BUILD=0)
endif()

## Build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release") # Default
endif()
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})


find_package(CUDAToolkit QUIET)
find_package(HIP QUIET)
if(CPU_BUILD)
elseif (NOT USE_HIP AND CUDAToolkit_FOUND)
    message(STATUS "CUDA found")
    enable_language(CUDA)
    add_compile_options(-DDEVICE_ENABLED)
    add_compile_options(-DCUDA_ENABLED)
    include_directories(SYSTEM ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}) # Silence warnings originating from CUDA headers

    ##TP: do not set GPU architectures to speedup debug compilations
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    else()
    	set(CMAKE_CUDA_ARCHITECTURES "60;61;70;80")
    endif()
    message(STATUS "Enabled CUDA architectures (Astaroth): " ${CMAKE_CUDA_ARCHITECTURES})
elseif (USE_HIP OR HIP_FOUND)
    set(USE_HIP ON)
    message(STATUS "HIP found")
    enable_language(HIP)
    add_compile_options(-DDEVICE_ENABLED)
    add_compile_options(-DHIP_ENABLED)
    include_directories(SYSTEM ${HIP_INCLUDE_DIRS}) # Silence warnings originating from HIP headers

    ##TP: do not set GPU architectures to speedup debug compilations
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    else()
    	set(CMAKE_HIP_ARCHITECTURES "gfx90a;gfx908")
    endif()
    message(STATUS "Enabled HIP architectures (Astaroth): " ${CMAKE_HIP_ARCHITECTURES})
else()
    message(ERROR "Did not find CUDA or HIP")
endif()

if (USE_HIP)
    add_definitions(-DAC_USE_HIP=1)
else()
    add_definitions(-DAC_USE_HIP=0)
endif()

#if (USE_OMP)
#find_package(OpenMP)
#if (OPENMP_FOUND)
#set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()
#endif()
include(FetchContent)

## Global flags
if (DOUBLE_PRECISION)
  add_definitions(-DAC_DOUBLE_PRECISION=1)
endif()

if (SINGLEPASS_INTEGRATION)
  add_definitions(-DSINGLEPASS=1)
endif()

if (USE_DARDEL)
  set(MPI_INCLUDE /opt/cray/pe/mpich/8.1.28/ofi/crayclang/17.0/include)
  set(MPI_LIB /opt/cray/pe/mpich/8.1.28/ofi/crayclang/17.0/lib)
  include_directories(${MPI_INCLUDE})
  link_directories(${MPI_LIB})
elseif(CMAKE_CXX_COMPILER MATCHES "mpic\\+\\+$")
    message(STATUS "Detected mpic++ as the C++ compiler so skipping linking MPI")
else()
  set(MPI_CXX_SKIP_MPICXX ON) # Skip MPI C++ bindings
  find_package(MPI REQUIRED COMPONENTS C)
  include_directories(${MPI_C_INCLUDE_DIRS})
endif()

add_definitions(-DAC_MPI_ENABLED=1)
add_definitions(-DAC_DEFAULT_CONFIG="${CMAKE_SOURCE_DIR}/config/astaroth.conf")

## Includes
include_directories(./submodule/include)                                  # Library headers
include_directories(./submodule/build/acc-runtime/api)
include_directories(./submodule/acc-runtime/api)

link_directories(./submodule/build/src/core)

add_compile_options(-DMODPRE=${MODULE_PREFIX} -DMODIN=${MODULE_INFIX} -DMODSUF=${MODULE_SUFFIX})

if (DOUBLE_PRECISION)
  add_definitions(-DAC_DOUBLE_PRECISION=1)
  set(PREC "dbl")
else()
  add_definitions(-DAC_DOUBLE_PRECISION=0)
  set(PREC "sgl")
endif()

if (RUNTIME_COMPILATION)
  add_definitions(-DAC_RUNTIME_COMPILATION=1)
else()
  add_definitions(-DAC_RUNTIME_COMPILATION=0)
endif()

if (TRANSPILATION)
  add_definitions(-DTRANSPILATION=1)
else()
  add_definitions(-DTRANSPILATION=0)
endif()

if (TRAINING)
  add_definitions(-DTRAINING=1)
else()
  add_definitions(-DTRAINING=0)
endif()

if (PACKED_DATA_TRANSFERS)
  add_compile_options(-DPACKED_DATA_TRANSFERS=1)
  add_library(astaroth_${PREC} SHARED gpu_astaroth.cc loadStore.cc)
else()

  if (TRAINING)	
    add_library(astaroth_${PREC} SHARED gpu_astaroth.cc torch.cc)
    include_directories(/opt/torchfort/include/)
    link_directories(/opt/torchfort/lib/)
    target_link_libraries(astaroth_${PREC} /opt/torchfort/lib/libtorchfort.so)
  else()
    add_library(astaroth_${PREC} SHARED gpu_astaroth.cc)
  endif()

endif()

target_link_libraries(astaroth_${PREC} astaroth_core)
