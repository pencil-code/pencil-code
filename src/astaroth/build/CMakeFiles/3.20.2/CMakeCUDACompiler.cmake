set(CMAKE_CUDA_COMPILER "/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/bin/nvcc")
set(CMAKE_CUDA_HOST_COMPILER "")
set(CMAKE_CUDA_HOST_LINK_LAUNCHER "/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/bin/g++")
set(CMAKE_CUDA_COMPILER_ID "NVIDIA")
set(CMAKE_CUDA_COMPILER_VERSION "11.7.64")
set(CMAKE_CUDA_DEVICE_LINKER "/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/bin/nvlink")
set(CMAKE_CUDA_FATBINARY "/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/bin/fatbinary")
set(CMAKE_CUDA_STANDARD_COMPUTED_DEFAULT "17")
set(CMAKE_CUDA_COMPILE_FEATURES "cuda_std_03;cuda_std_11;cuda_std_14;cuda_std_17")
set(CMAKE_CUDA03_COMPILE_FEATURES "cuda_std_03")
set(CMAKE_CUDA11_COMPILE_FEATURES "cuda_std_11")
set(CMAKE_CUDA14_COMPILE_FEATURES "cuda_std_14")
set(CMAKE_CUDA17_COMPILE_FEATURES "cuda_std_17")
set(CMAKE_CUDA20_COMPILE_FEATURES "")
set(CMAKE_CUDA23_COMPILE_FEATURES "")

set(CMAKE_CUDA_PLATFORM_ID "Linux")
set(CMAKE_CUDA_SIMULATE_ID "GNU")
set(CMAKE_CUDA_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_CUDA_SIMULATE_VERSION "11.3")



set(CMAKE_CUDA_COMPILER_ENV_VAR "CUDACXX")
set(CMAKE_CUDA_HOST_COMPILER_ENV_VAR "CUDAHOSTCXX")

set(CMAKE_CUDA_COMPILER_LOADED 1)
set(CMAKE_CUDA_COMPILER_ID_RUN 1)
set(CMAKE_CUDA_SOURCE_FILE_EXTENSIONS cu)
set(CMAKE_CUDA_LINKER_PREFERENCE 15)
set(CMAKE_CUDA_LINKER_PREFERENCE_PROPAGATES 1)

set(CMAKE_CUDA_SIZEOF_DATA_PTR "8")
set(CMAKE_CUDA_COMPILER_ABI "ELF")
set(CMAKE_CUDA_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_CUDA_LIBRARY_ARCHITECTURE "")

if(CMAKE_CUDA_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CUDA_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CUDA_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CUDA_COMPILER_ABI}")
endif()

if(CMAKE_CUDA_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_CUDA_COMPILER_TOOLKIT_ROOT "/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4")
set(CMAKE_CUDA_COMPILER_TOOLKIT_LIBRARY_ROOT "/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4")
set(CMAKE_CUDA_COMPILER_LIBRARY_ROOT "/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4")

set(CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES "/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/targets/x86_64-linux/include")

set(CMAKE_CUDA_HOST_IMPLICIT_LINK_LIBRARIES "")
set(CMAKE_CUDA_HOST_IMPLICIT_LINK_DIRECTORIES "/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/targets/x86_64-linux/lib/stubs;/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/targets/x86_64-linux/lib")
set(CMAKE_CUDA_HOST_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")

set(CMAKE_CUDA_IMPLICIT_INCLUDE_DIRECTORIES "/appl/spack/v018/install-tree/gcc-11.3.0/openmpi-4.1.4-w2aekq/include;/appl/spack/v018/install-tree/gcc-11.3.0/intel-oneapi-mkl-2022.1.0-p2wpb4/mkl/2022.1.0/include;/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/include/c++/11.3.0;/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/include/c++/11.3.0/x86_64-pc-linux-gnu;/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/include/c++/11.3.0/backward;/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/lib/gcc/x86_64-pc-linux-gnu/11.3.0/include;/usr/local/include;/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/include;/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/lib/gcc/x86_64-pc-linux-gnu/11.3.0/include-fixed;/usr/include")
set(CMAKE_CUDA_IMPLICIT_LINK_LIBRARIES "stdc++;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES "/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/targets/x86_64-linux/lib/stubs;/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/targets/x86_64-linux/lib;/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/lib64;/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/lib/gcc/x86_64-pc-linux-gnu/11.3.0;/lib64;/usr/lib64;/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/lib64;/appl/spack/v018/install-tree/gcc-11.3.0/openmpi-4.1.4-w2aekq/lib;/appl/spack/v018/install-tree/gcc-11.3.0/intel-oneapi-mkl-2022.1.0-p2wpb4/mkl/2022.1.0/lib/intel64;/appl/spack/v018/install-tree/gcc-11.3.0/intel-oneapi-mkl-2022.1.0-p2wpb4/lib;/appl/spack/v018/install-tree/gcc-8.5.0/gcc-11.3.0-i44hho/lib;/appl/spack/v018/install-tree/gcc-11.3.0/cuda-11.7.0-zucvj4/lib64/stubs")
set(CMAKE_CUDA_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")

set(CMAKE_CUDA_RUNTIME_LIBRARY_DEFAULT "STATIC")

set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_MT "")
