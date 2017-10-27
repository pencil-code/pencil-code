CUSOURCES = gpu/cuda/*.cu gpu/cuda/generic/*.cu gpu/cuda/core/*.cu
CCSOURCES = gpu/*.cc common/*.cc
MAIN_SRC = gpu_astaroth_new.cc
MAIN_OBJ = gpu_astaroth_new.o

include makefile.depend

OBJECTS = $(CUSOURCES:.cu=.o) $(CCSOURCES:.cc=.o) $(MAIN:.cc=.o)

# Settings for taito-gpu
CUFLAGS = -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 
CCFLAGS = -shared -Xcompiler -fPIC --debug
# Settings for Omer's machine
#CFLAGS = -gencode arch=compute_52,code=sm_52 --shared --relocatable-device-code=true --compiler-options '-fPIC'

#Importing environ variables into preprocessor
ENVIRON = -D MODPRE=${MODULE_PREFIX} -D MODIN=${MODULE_INFIX} -D MODSUF=${MODULE_SUFFIX}

%.o: %.cc
        nvcc $(CCFLAGS) $(ENVIRON) -dc $< -o $@

%.o: %.cu
        nvcc $(CUFLAGS) $(ENVIRON) -dc $< -o $@

astaroth_sgl.so: $(MAIN_OBJ) $(OBJECTS) 
#defines_PC.h defines_dims_PC.h
        nvcc $(CCFLAGS) $(CUFLAGS) -o astaroth_sgl.so $(OBJECTS)
        cp -f astaroth_sgl.so ../..

clean:
        @rm -f *.o
        @rm -f *.so
        @rm -f $(CHEADERS)

