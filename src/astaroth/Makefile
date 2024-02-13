# Settings for taito-gpu
#-gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 
#-gencode arch=compute_60,code=sm_60 --expt-relaxed-constexpr
OPTFLAGS = -O1 
#-funroll-loops -march=native 
# -gencode arch=compute_60,code=sm_60
#  Original

# MV's adaptation for TIARA cluster. Compiles with: 
#  1) cuda/9.0   3) openmpi/2.0.1_ic16.0            5) hdf5/1.8.16_openmpi_2.0.1_ic16.0 7) gcc/5.3.0
#  2) intel/2016 4) fftw/2.1.5_openmpi_2.0.1_ic16.0 6) cuda/9.1
#CCFLAGS = -ccbin gcc -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. 
##CCFLAGS = -ccbin icpc -Xcompiler -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. 


#Importing environ variables into preprocessor
ENVIRON = -DMODPRE=${MODULE_PREFIX} -DMODIN=${MODULE_INFIX} -DMODSUF=${MODULE_SUFFIX}

PRECISION =
include PC_modulesources.h
ifeq ($(PRECISION),DOUBLE)
  ENVIRON += -DAC_DOUBLE_PRECISION=1 -DDOUBLE_PRECISION
  PREC=dbl
#  OPTFLAGS += -fsingle-PRECision-constant
  CMAKE_PREC = 'DOUBLE_PRECISION=ON'
else
  ENVIRON += -DAC_DOUBLE_PRECISION=0
  PREC=sgl
  CMAKE_PREC = 'DOUBLE_PRECISION=OFF'
endif

SOURCES = gpu_astaroth.cc

DATA_TRANSFERS = PACKED
ifeq ($(DATA_TRANSFERS),PACKED)
  ENVIRON += -DPACKED_DATA_TRANSFERS=1
  CMAKE_PACKED = 'PACKED_DATA_TRANSFERS=ON'
  SOURCES += loadStore.cc
else
  ENVIRON += -DPACKED_DATA_TRANSFERS=0
  CMAKE_PACKED = 'PACKED_DATA_TRANSFERS=OFF'
endif

OBJECTS = $(SOURCES:.cc=.o) 

CCFLAGS = $(OPTFLAGS) -std=c++11 -shared -Xcompiler -fPIC --debug -I.. -Isubmodule/include -Isubmodule/build -g -G -lm
#CCFLAGS = -std=c++11 -Xcompiler --debug -I. -I../.. -I.. -Isubmodule -Isubmodule/src -Isubmodule/include -Isubmodule/src/standalone -g -G -lm
CCFLAGS = $(OPTFLAGS) -fPIC --debug -I.. -Isubmodule/include -Isubmodule/build -Isubmodule/build/acc-runtime/api -g -lm

PCHEADER_DIR=.
DSL_MODULE_DIR=Pencil
#DSL_MODULE_DIR=samples/mhd

IS_CRAY = $(shell ftn --version | grep -i Cray)
ifeq ($(IS_CRAY),)
#GCC with CUDA
#export CUDA_NVCC_FLAGS='-shared -Xcompiler -fPIC'
  CPP = nvcc
  MODULE_NAME_ORDER=
  CMAKE_USE_HIP=OFF
else
#CRAY with HIP
  CPP = hipcc
  CCFLAGS += -DAC_USE_HIP=1
  CMAKE_USE_HIP=ON
  MODULE_NAME_ORDER=CRAY
  MODULE_INFIX=$
  MODULE_SUFFIX=_
#CC -xhip
endif
SINGLEPASS_INTEGRATION=OFF

include Makefile.depend
include Makefile.extern

astaroth_$(PREC).so: PC_modulesources.h $(MODULESOURCES) $(PCHEADER_DIR)/PC_moduleflags.h $(PC_HEADERS) $(OBJECTS) #astaroth_libs
	$(CPP) $(CCFLAGS) -shared -o astaroth_$(PREC).so $(OBJECTS) -L submodule/build/src/core -L submodule/build/src/core/kernels -L submodule/build/src/utils -lastaroth_core -lkernels -lastaroth_utils -Isubmodule/include

#astaroth_$(PREC).so: PC_modulesources.h $(MODULESOURCES) $(PCHEADER_DIR)/PC_moduleflags.h $(PC_HEADERS) $(OBJECTS) astaroth_libs
#	mkdir -p build && cd build && \
#	cmake -DMULTIGPU_ENABLED=OFF -DDSL_MODULE_DIR=../acc-runtime/$(DSL_MODULE_DIR) -DUSE_OMP=ON -DMPI_ENABLED=OFF -DBUILD_SHARED_LIBS=ON -DSINGLEPASS_INTEGRATION=OFF -DUSE_CUDA_AWARE_MPI=OFF -DUSE_EXTERNAL_DECOMP=ON -DLINEAR_PROC_MAPPING=ON -DBUILD_SAMPLES=OFF -DBUILD_STANDALONE=OFF ..
#	make -j 
#	cd ..
#	cp ./build/libastaroth_$(PREC).so astaroth_$(PREC).so

astaroth_libs: $(PC_HEADERS)
ifeq ($(SINGLEPASS_INTEGRATION),OFF)
	@cp submodule/acc-runtime/$(DSL_MODULE_DIR)/solve.ac submodule/acc-runtime/$(DSL_MODULE_DIR)/solve.ac.sv
	@sed -i -f solve.sed submodule/acc-runtime/$(DSL_MODULE_DIR)/solve.ac; rm -f sedtmp
endif
	@mkdir -p submodule/build && cd submodule/build && \
	cmake -DMULTIGPU_ENABLED=OFF -DDSL_MODULE_DIR=../acc-runtime/$(DSL_MODULE_DIR) -DUSE_OMP=ON -DMPI_ENABLED=OFF -DBUILD_SHARED=ON -DSINGLEPASS_INTEGRATION=$(SINGLEPASS_INTEGRATION) -DUSE_CUDA_AWARE_MPI=OFF\
              -D$(CMAKE_PACKED) -D$(CMAKE_PREC) -DBUILD_SAMPLES=OFF  -DADJACENT_VERTEX_BUFFERS=ON -DUSE_HIP=$(CMAKE_USE_HIP) -DLINEAR_PROC_MAPPING=ON -DUSE_EXTERNAL_DECOMP=ON -DVERBOSE=OFF -DCMAKE_BUILD_TYPE=DEBUG -DOPTIMIZE_MEM_ACCESSES=ON -DAUTOTUNE_TASKGRAPHS=OFF -DBUILD_STANDALONE=OFF .. \
	&& make -I../../.. VERBOSE=0 # Build the core Astaroth libraries
	@echo PRECISION=$(CMAKE_PREC)
	@echo ENVIRON: $(ENVIRON)

#	@echo PRECISION=$(CMAKE_PREC)
#-DCMAKE_CXX_COMPILER=/appl/spack/install-tree/gcc-4.8.5/gcc-8.3.0-qzmzn5/bin/c++ .. \
#export DSL_MODULE=pc_mhd_solver && \

%.o: %.cc
	$(CPP) $(CCFLAGS) $(ENVIRON) -c $< -o $@

$(PC_HEADERS): .sentinel

.sentinel: ../scripts/phys_modules2c PC_modulesources.h $(MODULESOURCES)
	@rm -f $(PCHEADER_DIR)/PC_modulepars.h submodule/acc-runtime/$(DSL_MODULE_DIR)/PC_modulepardecs.h 
#diagnostics/PC_modulediags.h
	@echo '// automatically generated; do not edit!' > $(PCHEADER_DIR)/PC_modulepars.h
	@echo '// automatically generated; do not edit!' > submodule/acc-runtime/$(DSL_MODULE_DIR)/PC_modulepardecs.h
	@export DSL_MODULE_DIR=$(DSL_MODULE_DIR); \
        for file in $(MODULESOURCES); do echo processing $$file ...; if [[ $$file =~ eos[_.] ]]; \
	then ../scripts/phys_modules2c par $$file equationofstate; \
        else ../scripts/phys_modules2c par $$file; fi; done
	@touch .sentinel
#
#@echo '// automatically generated; do not edit!' > diagnostics/PC_modulediags.h
#@rm -f diagnostics/PC_modulediags_init.h
#@sed -e's/_diags_/_diags_init_/' <  diagnostics/PC_modulediags.h > diagnostics/PC_modulediags_init.h
#        @for file in $(MODULESOURCES); do echo processing $$file ...; if [[ $$file =~ eos[_.] ]]; \
#        then ../scripts/phys_modules2c par $$file equationofstate; ../scripts/phys_modules2c diag $$file equationofstate;\
#        else ../scripts/phys_modules2c par $$file; ../scripts/phys_modules2c diag $$file; fi; done
#
clean:
	@rm -f $(OBJECTS)
	@rm -f *.so .sentinel
	@if [[ -d submodule/build ]]; then rm -f -rI submodule/build; fi # Removes the files generated with cmake for the Astaroth code
#
# Deep cleaning, e.g. necessary with compiler change
#
cleann: clean
	@rm -f ../cparam_c.h ../defines_cdata.h ../diagnostics_c.h ../forcing_c.h ../forcing_pars_c.h ../sub_c.h ../mpicomm_c.h ../boundcond_c.h
	@rm -f submodule/acc-runtime/$(DSL_MODULE_DIR)/PC_modulepardecs.h
	@rm -f $(PC_HEADER_DIR)/PC_modulepars.h

#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < submodule/acc-runtime/$(DSL_MODULE_DIR)/PC_modulepardecs.h | xargs rm -f
#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < $(PC_HEADER_DIR)/PC_modulepars.h | xargs rm -f
#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < diagnostics/PC_modulediags.h | xargs rm -f
#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < diagnostics/PC_modulediags_init.h | xargs rm -f
