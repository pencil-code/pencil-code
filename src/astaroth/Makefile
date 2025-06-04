SHELL=/bin/bash

#Importing environ variables into preprocessor
#ENVIRON = -DMODPRE=${MODULE_PREFIX} -DMODIN=${MODULE_INFIX} -DMODSUF=${MODULE_SUFFIX}

include PC_modulesources.h

ifeq ($(PRECISION),DOUBLE)
  ENVIRON = -DDOUBLE_PRECISION=1
  PREC=dbl
#  OPTFLAGS += -fsingle-PRECision-constant
  CMAKE_PREC = ON
else
  ENVIRON = -DDOUBLE_PRECISION=0
  PREC=sgl
  CMAKE_PREC = OFF
endif

include Makefile.vendor
include Makefile.target
include Makefile.module_info

ifeq ($(GPU_VENDOR),amd)
  USE_HIP = ON
else
  USE_HIP = OFF
endif

ifeq ($(MPICOMM),nompicomm)
  CMAKE_MPI = OFF
  $(error interfaced version requires MPI)
else
  CMAKE_MPI = ON
  $(info MPI enabled)
endif

SOURCES = gpu_astaroth.cc 

DATA_TRANSFERS = 
#PACKED
ifeq ($(DATA_TRANSFERS),PACKED)
  #ENVIRON += -DPACKED_DATA_TRANSFERS=1
  CMAKE_PACKED = ON
  SOURCES += loadStore.cc loadStore.h
else
  #ENVIRON += -DPACKED_DATA_TRANSFERS=0
  CMAKE_PACKED = OFF
endif

CMAKE_SINGLEPASS = OFF
OBJECTS = $(SOURCES:.cc=.o) 

#DSL_MODULE_DIR=samples/gputest
DSL_MODULE_DIR = ../../DSL/local
CWD = $(PWD)

ifeq ($(CMAKE_SINGLEPASS),ON)
  DSL_MODULE_FILE = solve_single.ac
else
  DSL_MODULE_FILE = solve_two.ac
endif

EXT_KERNELS =
ifneq ($(filter shock,$(MODULESOURCES)),)
  EXT_KERNELS += $(wildcard $(PENCIL_HOME)/src/astaroth/DSL/shock/*.ac)
endif

# C headers in submodule
AC_HEADERS=submodule/include/astaroth.h submodule/acc-runtime/api/math_utils.h

# C headers in pencil-code/src 
CHEADERS = ../cparam_c.h ../cdata_c.h ../boundcond_c.h ../defines_cdata.h ../forcing_c.h ../sub_c.h ../mpicomm_c.h

# Pencil Code related C headers in . (pencil-code/src/astaroth)
PCHEADER_DIR = .
PC_HEADERS = $(PCHEADER_DIR)/PC_moduleflags.h $(PCHEADER_DIR)/PC_module_parfuncs.h $(PCHEADER_DIR)/PC_modulepars.h 

RUNTIME_COMPILATION      ?= off
TRAINING                 ?= notraining
CMAKE_TRAINING_VALUE     ?= off
TRANSPILATION            ?= off
INLINING                 ?= off
READ_OVERRIDES           ?= off
ELIMINATE_CONDITIONALS   ?= $(TRANSPILATION)
OPTIMIZE_INPUT_PARAMS    ?= on

ifeq ($(RUNTIME_COMPILATION),on) 
	CMAKE_BUILD_TYPE ?= RelWithDebInfo
endif
ifeq ($(RUNTIME_COMPILATION),ON) 
	CMAKE_BUILD_TYPE ?= RelWithDebInfo
endif

ifeq ($(TRAINING),training_torchfort) 
	SUPPRESS_COMPILER_WARNINGS = on
	CMAKE_TRAINING_VALUE = on
        SOURCES += torch.cc
endif

CMAKE_BUILD_TYPE          ?= Debug
SUPPRESS_COMPILER_WARNINGS ?= off

astaroth_$(PREC).so: $(SOURCES) forcing.h $(AC_HEADERS) $(PC_HEADERS) $(CHEADERS) astaroth_libs
	mkdir -p build && \
	echo USE_HIP=$(USE_HIP) && \
	echo PRECISION=$(CMAKE_PREC) && \
	cd build && export LDFLAGS=""  && \
	cmake -DDOUBLE_PRECISION=$(CMAKE_PREC) -DRUNTIME_COMPILATION=$(RUNTIME_COMPILATION) -DTRANSPILATION=$(TRANSPILATION) -DUSE_HIP=$(USE_HIP) -DSINGLEPASS_INTEGRATION=$(CMAKE_SINGLEPASS) -DTRAINING=$(CMAKE_TRAINING_VALUE) .. && \
	make -j && \
	cd .. && \
	mv ./build/libastaroth_$(PREC).so libastaroth_$(PREC).so

#	mkdir -p scheduler-build && \
#	cd scheduler-build && \
#	cmake -DUSE_OPENMP=ON ../NUMA-scheduler && \
#	make -j && \
#	cd .. && \
#	cp -r ./NUMA-scheduler/include .

# -DRUNTIME_COMPILE=$(PRODUCTION)

CMAKE_OPTIONS = -DMULTIGPU_ENABLED=OFF -DMPI_ENABLED=$(CMAKE_MPI) -DBUILD_SHARED_LIBS=ON -DSINGLEPASS_INTEGRATION=$(CMAKE_SINGLEPASS) -DMAX_THREADS_PER_BLOCK=512 \
                -DPACKED_DATA_TRANSFERS=$(CMAKE_PACKED) -DDOUBLE_PRECISION=$(CMAKE_PREC) -DBUILD_SAMPLES=OFF -DUSE_HIP=$(USE_HIP) -DSUPPRESS_COMPILER_WARNINGS=$(SUPPRESS_COMPILER_WARNINGS) \
		-DVERBOSE=off -DOPTIMIZE_INPUT_PARAMS=$(OPTIMIZE_INPUT_PARAMS) -DOPTIMIZE_MEM_ACCESSES=ON -DOPTIMIZE_ARRAYS=$(RUNTIME_COMPILATION) -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) \
	        -DBUILD_STANDALONE=OFF -DOPTIMIZE_FIELDS=$(TRANSPILATION) -DALLOW_DEAD_VARIABLES=$(RUNTIME_COMPILATION) -DELIMINATE_CONDITIONALS=$(ELIMINATE_CONDITIONALS) -DBUILD_UTILS=OFF \
		-DCMAKE_POSITION_INDEPENDENT_CODE=ON -DDSL_MODULE_FILE=$(DSL_MODULE_FILE) -DINLINING=$(INLINING) -DBUILD_ACM=OFF
#
# External sources are provided with absolute pathnames to avoid any trouble with symbolic links.
#
astaroth_libs: $(PC_HEADERS)
	@ext_srcs="`ls -m $$PENCIL_HOME/src/astaroth/DSL/*/*.h \
                          $$PENCIL_HOME/src/astaroth/DSL/*.h \
		          $$PENCIL_HOME/src/astaroth/submodule/acc-runtime/stdlib/*.h \
			  $$PENCIL_HOME/src/astaroth/submodule/acc-runtime/stdlib/*/*.h \
                          $(CWD)/DSL/local/*.h \
			  |tr "\n" " "|sed -e's/, / /g'` $(EXT_KERNELS)" && \
	echo "const char* cmake_options = \" -DEXT_SOURCES=$$ext_srcs -DDSL_MODULE_DIR=../$(DSL_MODULE_DIR) $(CMAKE_OPTIONS)\";" > cmake_options.h && \
	mkdir -p submodule/build  && cd submodule/build && export LDFLAGS=""  &&  \
	cmake  $(CMAKE_OPTIONS) -DEXT_SOURCES=$$ext_srcs -DDSL_MODULE_DIR=../acc-runtime/$(DSL_MODULE_DIR) -DRUNTIME_COMPILATION=$(RUNTIME_COMPILATION) -DREAD_OVERRIDES=$(READ_OVERRIDES) ..
	cd submodule/build && make -j -I../../.. VERBOSE=0 # Build the core Astaroth library
	@echo PRECISION=$(CMAKE_PREC)

#        CMAKE_PREFIX_PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/24.1/comm_libs/mpi/lib /opt/nvidia/hpc_sdk/Linux_x86_64/24.1/comm_libs/mpi/include" \
update_astaroth_libs:
	cd submodule/build && export LDFLAGS="" && \
	make -I../../.. VERBOSE=0 # Build the core Astaroth library with existing Makefile.
	@echo PRECISION=$(CMAKE_PREC)

include Makefile.extern

$(PC_HEADERS): .sentinel

.sentinel: ../scripts/phys_modules2c PC_modulesources.h $(MODULESOURCES)
	@rm -f $(PCHEADER_DIR)/PC_modulepars.h $(PCHEADER_DIR)/PC_module*funcs.h submodule/acc-runtime/$(DSL_MODULE_DIR)/PC_modulepardecs.h 
	@echo '// automatically generated; do not edit!' > $(PCHEADER_DIR)/PC_modulepars.h
	@echo '// automatically generated; do not edit!' > submodule/acc-runtime/$(DSL_MODULE_DIR)/PC_modulepardecs.h
ifeq ($(TRANSPILATION),on)
	@export DSL_MODULE_DIR=$(DSL_MODULE_DIR); \
        for file in $(MODULESOURCES); do echo processing $$file ...; if [[ $$file =~ eos[_.] ]]; \
	then ../scripts/phys_modules2c_transpilation par $$file equationofstate; \
        else ../scripts/phys_modules2c_transpilation par $$file; fi; done
	@touch .sentinel
else
	@export DSL_MODULE_DIR=$(DSL_MODULE_DIR); \
        for file in $(MODULESOURCES); do echo processing $$file ...; if [[ $$file =~ eos[_.] ]]; \
	then ../scripts/phys_modules2c par $$file equationofstate; \
        else ../scripts/phys_modules2c par $$file; fi; done
	@touch .sentinel
endif
#
clean:
	@rm -f $(OBJECTS)
	@rm -f *.so .sentinel cmake_options.h
	@if [[ -d submodule/build ]]; then rm -f -rI submodule/build; fi # Removes the files generated with cmake for the Astaroth code
	@rm -rf build
	@rm -rf scheduler-build
	@rm -rf include
#
# Deep cleaning, e.g. necessary with compiler change
#
cleann: clean
	@rm -f ../cparam_c.h ../defines_cdata.h ../forcing_c.h ../forcing_pars_c.h ../sub_c.h ../mpicomm_c.h ../boundcond_c.h ../cdata_c.h
	@rm -f submodule/acc-runtime/$(DSL_MODULE_DIR)/PC_modulepardecs.h
	@rm -f $(PC_HEADER_DIR)/PC_modulepars.h $(PCHEADER_DIR)/PC_module_parfuncs.h
#
#--------------------------------------------------------------------------------------------------------------------------------------------
OPTFLAGS = -O1 
#-funroll-loops -march=native 
# -gencode arch=compute_60,code=sm_60
#  Original

# MV's adaptation for TIARA cluster. Compiles with: 
#  1) cuda/9.0   3) openmpi/2.0.1_ic16.0            5) hdf5/1.8.16_openmpi_2.0.1_ic16.0 7) gcc/5.3.0
#  2) intel/2016 4) fftw/2.1.5_openmpi_2.0.1_ic16.0 6) cuda/9.1
#CCFLAGS = -ccbin gcc -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. 
##CCFLAGS = -ccbin icpc -Xcompiler -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. 

#CCFLAGS = $(OPTFLAGS) -std=c++11 -shared -Xcompiler -fPIC --debug -I.. -Isubmodule/include -Isubmodule/build -g -G -lm
CCFLAGS = $(OPTFLAGS) -Xcompiler -shared -fPIC --debug -I.. -Isubmodule/include -Isubmodule/build -Isubmodule/build/acc-runtime/api -g -G -lm
#CCFLAGS = -std=c++11 -Xcompiler --debug -I. -I../.. -I.. -Isubmodule -Isubmodule/src -Isubmodule/include -Isubmodule/src/standalone -g -G -lm

#export CUDA_NVCC_FLAGS='-shared -Xcompiler -fPIC'
#export CUDA_NVCC_FLAGS='-shared -fPIC'

#astaroth_$(PREC).so: PC_modulesources.h $(MODULESOURCES) $(PCHEADER_DIR)/PC_moduleflags.h $(PC_HEADERS) astaroth_libs $(OBJECTS)
#	$(CPP) $(CCFLAGS) $(CUFLAGS) -shared -o astaroth_$(PREC).so $(OBJECTS) -L submodule/build/src/core -L submodule/build/src/core/kernels -L submodule/build/src/utils -lastaroth_core -lkernels -Isubmodule/include
#-lastaroth_utils

#-DCMAKE_CXX_COMPILER=/appl/spack/install-tree/gcc-4.8.5/gcc-8.3.0-qzmzn5/bin/c++ .. \

ifeq ($(USE_HIP),OFF)
#GCC with CUDA
  CPP = nvcc
else
#CRAY with HIP
  CPP = hipcc
endif

%.o: %.cc
	$(CPP) $(CCFLAGS) $(ENVIRON) -c $< -o $@

#@echo '// automatically generated; do not edit!' > diagnostics/PC_modulediags.h
#@rm -f diagnostics/PC_modulediags_init.h
#@sed -e's/_diags_/_diags_init_/' <  diagnostics/PC_modulediags.h > diagnostics/PC_modulediags_init.h
#        @for file in $(MODULESOURCES); do echo processing $$file ...; if [[ $$file =~ eos[_.] ]]; \
#        then ../scripts/phys_modules2c par $$file equationofstate; ../scripts/phys_modules2c diag $$file equationofstate;\
#        else ../scripts/phys_modules2c par $$file; ../scripts/phys_modules2c diag $$file; fi; done

#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < submodule/acc-runtime/$(DSL_MODULE_DIR)/PC_modulepardecs.h | xargs rm -f
#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < $(PC_HEADER_DIR)/PC_modulepars.h | xargs rm -f
#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < diagnostics/PC_modulediags.h | xargs rm -f
#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < diagnostics/PC_modulediags_init.h | xargs rm -f

