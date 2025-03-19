
USE_HIP=$(<use_hip)
CMAKE_PREC=$(<cmake_prec)
mkdir -p build
echo $USE_HIP
echo $CMAKE_PREC
cd build
cmake -DMULTIGPU_ENABLED=ON -DDSL_MODULE_DIR=gputest -D$CMAKE_PREC -DUSE_HIP=$USE_HIP -DUSE_OMP=ON -DMPI_ENABLED=ON -DBUILD_SHARED_LIBS=ON -DSINGLEPASS_INTEGRATION=OFF -DUSE_CUDA_AWARE_MPI=OFF -DUSE_EXTERNAL_DECOMP=ON -DLINEAR_PROC_MAPPING=ON ..
make -j
cd ..
cp ./build/libastaroth.so libastaroth.so
cp libastaroth.so astaroth_sgl.so
cp libastaroth.so astaroth_dbl.so
mkdir -p scheduler-build
cd scheduler-build
cmake -DUSE_OPENMP=ON ../NUMA-scheduler
make -j
cd ..
cp -r ./NUMA-scheduler/include .


