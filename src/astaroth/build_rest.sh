CC -O1  -fPIC --debug -I.. -Isubmodule/include -Isubmodule/build -Isubmodule/build/acc-runtime/api -g -G -lm -DMODPRE=__ -DMODIN=_MOD_ -DMODSUF= -DAC_DOUBLE_PRECISION=0 -DPACKED_DATA_TRANSFERS=1 -DAC_MPI_ENABLED=1 -DAC_USE_HIP=1 -c loadStore.cc

CC -O1  -fPIC --debug -I.. -Isubmodule/include -Isubmodule/build -Isubmodule/build/acc-runtime/api -g -G -lm -DMODPRE=__ -DMODIN=_MOD_ -DMODSUF= -DAC_DOUBLE_PRECISION=0 -DPACKED_DATA_TRANSFERS=1 -DAC_MPI_ENABLED=1 -DAC_USE_HIP=1 -c gpu_astaroth.cc


