echo $DSL_MODULE
echo $CMAKE_PREC
mkdir -p build
cd build
cmake -DMULTIGPU_ENABLED=ON -DDSL_MODULE_DIR=$DSL_MODULE -D$CMAKE_PREC  -DUSE_OMP=ON -DMPI_ENABLED=ON -DBUILD_SHARED_LIBS=ON -DSINGLEPASS_INTEGRATION=OFF -DUSE_CUDA_AWARE_MPI=OFF -DUSE_EXTERNAL_DECOMP=ON -DLINEAR_PROC_MAPPING=ON ..
make -j
cd ..
cp ./build/libastaroth.so libastaroth.so
cp libastaroth.so astaroth_sgl.so
cp libastaroth.so astaroth_dbl.so
cp ./build/NUMA-scheduler/src/mod/mt.mod ../.


