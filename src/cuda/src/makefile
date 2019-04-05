
#Note that you only have to include the main file (compute.o)
#and modules that are *separately compiled*. Everything else is 
#automatically included in the compilation.

hostobjs = 
#model_collectiveops.o

deviceobjs = diagnostics.o integrators_v5.o collectiveops.o copyHalosConcur.o timestep.o gpu_astaroth.o copyouterhalostodevice.o copyinternalhalostohost.o copy_halos.o
#boundcond.o diff.o shear.o coriolis.o forcing.o continuity.o hydro.o

hostcpps = initutils.cpp io.cpp model_collectiveops.cpp
devicecus = slice.cu collectiveops.cu diagnostics.cu integrators_v5.cu boundcond.cu timestep.cu diff.cu shear.cu coriolis.cu forcing.cu continuity.cu hydro.cu gpu_astaroth.cu copyouterhalostodevice.cu copyinternalhalostohost.cu copy_halos.cu slice.cu

#Flags: 
#
#-g -G compiles debuggable code
#
#-O2 -use_fast_math uses all optimizations and replaces
#math functions with less precise, but faster versions
#
#--maxrregcount=64 sets maximum registers used by a 
#CUDA thread; Basically low register count increases
#occupancy in register limited programs, but causes
#register spilling to L1.
#
#-Xptxas -dlcm=ca caches global memory accesses to
#L1, but may cause overfetch because transactions
#from L1 are done in 128-byte cache lines, but
#transactions from L2/SGRAM are done as 32-byte 
#blocks.
# Settings for taito-gpu
CFLAGS = -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -shared -Xcompiler -fPIC --debug
# Settings for Omer's machine
#CFLAGS = -gencode arch=compute_52,code=sm_52 --shared --relocatable-device-code=true --compiler-options '-fPIC'

#Importing environ variables into preprocessor
ENVIRON = -D MODPRE=${MODULE_PREFIX} -D MODIN=${MODULE_INFIX} -D MODSUF=${MODULE_SUFFIX}

#Searches also src/verification for additional includes
VPATH = verification

#Include macro flags and their user defined values
#include makefile.local

#-rdc=true?
#all: hostobjs deviceobjs main
	#nvcc  -g -G -Xptxas -O0 $(CFLAGS) $(hostobjs) $(deviceobjs) $(main) -o ac_run_gdb
	#nvcc -g -lineinfo --ptxas-options=-v $(CFLAGS) $(hostobjs) $(deviceobjs) $(main) -o ac_run_gdb
	#nvcc -g $(CFLAGS) $(hostobjs) $(deviceobjs) $(main) -o ac_run_gdb	

%.o: %.cpp
	nvcc $(CFLAGS) $(ENVIRON) -dc $< -o $@

%.o: %.cu
	nvcc ${MFLAGS} $(CFLAGS) $(ENVIRON) -dc $< -o $@

#main: $(main)
#%.o: %.cu
	#nvcc -g -G -Xptxas -O0 ${INIT_CONFS} ${RUN_CONFS} $(CFLAGS) -dc $< -o $@
	#nvcc -g -lineinfo --ptxas-options=-v ${INIT_CONFS} ${RUN_CONFS} $(CFLAGS) -dc $< -o $@
	#nvcc -g ${INIT_CONFS} ${RUN_CONFS} $(CFLAGS) -dc $< -o $@

astaroth_sgl.so: $(hostobjs) $(deviceobjs) defines_PC.h defines_dims_PC.h
	nvcc $(CFLAGS) ${MFLAGS} -o astaroth_sgl.so $(hostobjs) $(deviceobjs)
	cp -f astaroth_sgl.so ..

include makefile.depend
hostobjs: $(hostobjs)
deviceobjs: $(deviceobjs)

clean:
	@rm -f *.o
	@rm -f *.so
	@rm -f $(CHEADERS)

