echo $DSL_MODULE
mkdir -p build
cd build
cmake -DMULTIGPU_ENABLED=ON -DDSL_MODULE_DIR=$DSL_MODULE -DUSE_OMP=ON -DMPI_ENABLED=ON -DBUILD_SHARED_LIBS=ON -DSINGLEPASS_INTEGRATION=OFF -DUSE_CUDA_AWARE_MPI=OFF -DUSE_EXTERNAL_DECOMP=ON -DLINEAR_PROC_MAPPING=ON ..
make -j
cd ..
cp ./build/libastaroth_sgl.so libastaroth_sgl.so
cp libastaroth_sgl.so astaroth_sgl.so

