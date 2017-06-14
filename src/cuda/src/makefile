
#Note that you only have to include the main file (compute.o)
#and modules that are *separately compiled*. Everything else is 
#automatically included in the compilation.
hostobjs = initutils.o io.o model_collectiveops.o 
deviceobjs = slice.o collectiveops.o diagnostics.o integrators_v5s.o boundcond.o timestep.o diff.o shear.o coriolis.o forcing.o continuity.o hydro.o gpu_astaroth.o copyouterhalostodevice.o copyinternalhalostohost.o copy_halos.o copyHalosConcur.o

hostcpps = initutils.cpp io.cpp model_collectiveops.cpp
devicecus = collectiveops.cu diagnostics.cu integrators_v5s.cu boundcond.cu timestep.cu diff.cu shear.cu coriolis.cu forcing.cu continuity.cu hydro.cu gpu_astaroth.cu copyouterhalostodevice.cu copyinternalhalostohost.cu copy_halos.cu slice.cu
#deviceobjs = collectiveops.o diagnostics.o integrators_debug.o boundcond.o timestep.o diff.o shear.o coriolis.o forcing.o continuity.o hydro.o
#deviceobjs = collectiveops.o diagnostics.o integrators_ytraversing3.o boundcond.o timestep.o diff.o shear.o coriolis.o forcing.o continuity.o hydro.o
#deviceobjs = collectiveops.o diagnostics.o integrators_Th38x10.o boundcond.o timestep.o diff.o shear.o coriolis.o forcing.o continuity.o hydro.o
#main = compute.o

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
CFLAGS = -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -shared -Xcompiler -fPIC
# Settings for Omer's machine
#CFLAGS = -gencode arch=compute_52,code=sm_52 --shared --relocatable-device-code=true --compiler-options '-fPIC'

#Searches also src/verification for additional includes
VPATH = verification

#Include initial configuration
include init.conf

#Include run configuration
include run.conf

#Include macro flags and their user defined values
include makefile.local

#-rdc=true?
#all: hostobjs deviceobjs main
	#nvcc  -g -G -Xptxas -O0 $(CFLAGS) $(hostobjs) $(deviceobjs) $(main) -o ac_run_gdb
	#nvcc -g -lineinfo --ptxas-options=-v $(CFLAGS) $(hostobjs) $(deviceobjs) $(main) -o ac_run_gdb
	#nvcc -g $(CFLAGS) $(hostobjs) $(deviceobjs) $(main) -o ac_run_gdb	

hostobjs: $(hostobjs)
%.o: %.cpp
	nvcc -g ${INIT_CONFS} ${RUN_CONFS} $(CFLAGS) -dc $< -o $@

deviceobjs: $(deviceobjs)
%.o: %.cu
	nvcc -g ${INIT_CONFS} ${RUN_CONFS} ${MFLAGS} $(CFLAGS) -dc $< -o $@

#main: $(main)
#%.o: %.cu
	#nvcc -g -G -Xptxas -O0 ${INIT_CONFS} ${RUN_CONFS} $(CFLAGS) -dc $< -o $@
	#nvcc -g -lineinfo --ptxas-options=-v ${INIT_CONFS} ${RUN_CONFS} $(CFLAGS) -dc $< -o $@
	#nvcc -g ${INIT_CONFS} ${RUN_CONFS} $(CFLAGS) -dc $< -o $@
astaroth_sgl.so: $(hostobjs) $(deviceobjs) defines.h
	nvcc -g $(CFLAGS) ${INIT_CONFS} ${RUN_CONFS} ${MFLAGS} -o astaroth_sgl.so $(hostobjs) $(deviceobjs)
	cp -f astaroth_sgl.so ..

clean:
	@rm -f *.o
	@rm -f runme
	@rm data/*.ac
	@rm data/*.dat

