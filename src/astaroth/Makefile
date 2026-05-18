SHELL=/bin/bash

#Importing environ variables into preprocessor
#ENVIRON = -DMODPRE=${MODULE_PREFIX} -DMODIN=${MODULE_INFIX} -DMODSUF=${MODULE_SUFFIX}
#
# Get recent target.
#       
target := $(firstword $(MAKECMDGOALS))
#
# Avoid missing parts when cleaning
#
ifeq ($(filter clean cleann,$(target)),)
include PC_modulesources.h
include Makefile.module_info
endif

ifeq ($(PRECISION),DOUBLE)
  ENVIRON = -DDOUBLE_PRECISION=1
  PREC=dbl
  CMAKE_PREC = ON
else
ifeq ($(PRECISION),double)
  ENVIRON = -DDOUBLE_PRECISION=1
  PREC=dbl
  CMAKE_PREC = ON
else
  CMAKE_PREC = ON
  ENVIRON = -DDOUBLE_PRECISION=0
  PREC=sgl
  CMAKE_PREC = OFF
endif
endif

ifeq ($(MPICOMM),nompicomm)
  CMAKE_MPI = OFF
  $(error interfaced version requires MPI)
else
  CMAKE_MPI = ON
  $(info MPI enabled)
endif

SOURCES = gpu_astaroth.cc 


CMAKE_SINGLEPASS = OFF
OBJECTS = $(SOURCES:.cc=.o) 

#DSL_MODULE_DIR=samples/gputest
DSL_MODULE_DIR = DSL/local
CWD = $(PWD)

DSL_MODULE_FILE = solve.ac

EXT_KERNELS =
ifneq ($(filter shock,$(MODULESOURCES)),)
  EXT_KERNELS += $(wildcard $(PENCIL_HOME)/src/astaroth/DSL/shock/*.ac)
endif

# C headers in pencil-code/src 
CHEADERS = ../cparam_c.h ../boundcond_c.h ../forcing_c.h ../sub_c.h #../mpicomm_c.h

# Pencil Code related C headers in . (pencil-code/src/astaroth)
PCHEADER_DIR = .
PC_HEADERS = $(PCHEADER_DIR)/PC_moduleflags.h $(PCHEADER_DIR)/PC_module_parfuncs.h $(PCHEADER_DIR)/PC_modulepars.h 

RUNTIME_COMPILATION      ?= off
TRAINING                 ?= notraining
CMAKE_TRAINING_VALUE     ?= off
TRANSPILATION            ?= off
GENERATE_DSL_CODE        ?= off
MODIFY_SOURCE_CODE       ?= off
INLINING                 ?= off
READ_OVERRIDES           ?= off
USE_HEFFTE               ?= off
SUPPORT_OLD_ARCHITECTURES ?= ON
GPU_FFT                  ?= VENDOR
GPU_FFT_ENABLED              ?= ON

ifeq ($(GPU_FFT),HEFFTE)
	USE_HEFFTE = on
endif

FULL_DSL_OPTIMIZATIONS = off
ifeq ($(TRANSPILATION),on)
  FULL_DSL_OPTIMIZATIONS = on
endif

ifeq ($(READ_OVERRIDES),on)
  FULL_DSL_OPTIMIZATIONS = on
endif

ELIMINATE_CONDITIONALS   ?= $(TRANSPILATION)
CPU_BUILD                ?= off
OPTIMIZE_INPUT_PARAMS    ?= on
USE_ONLY_OPTIMIZED_KERNELS ?= on

ifeq ($(RUNTIME_COMPILATION),on) 
	CMAKE_BUILD_TYPE ?= RelWithDebInfo
endif
ifeq ($(RUNTIME_COMPILATION),ON) 
	CMAKE_BUILD_TYPE ?= RelWithDebInfo
endif

#ifeq ($(TRAINING),training_torchfort) 
ifneq ($(TORCHFORT_PATH),)
	SUPPRESS_COMPILER_WARNINGS = on
endif

ifeq ($(TRAINING),training_torchfort) 
	CMAKE_TRAINING_VALUE = on
  SOURCES += torch.cc
endif

CMAKE_BUILD_TYPE          ?= Debug
SUPPRESS_COMPILER_WARNINGS ?= off

#TP: mpicc does not work nicely with CMake.
#    Better to let CMake choose the compiler itself
ifeq ($(CC),mpicc)
  override CC = 
endif
ifeq ($(CXX),mpicc)
  override CXX = 
endif

astaroth_$(PREC).so: $(SOURCES) forcing.h $(PC_HEADERS) $(CHEADERS) astaroth_libs
	mkdir -p build && \
	echo PRECISION=$(CMAKE_PREC) && \
	cd build && export LDFLAGS=""  && \
	cmake -DCPU_BUILD=$(CPU_BUILD) -DDOUBLE_PRECISION=$(CMAKE_PREC) -DRUNTIME_COMPILATION=$(RUNTIME_COMPILATION) -DTRANSPILATION=$(TRANSPILATION) -DSINGLEPASS_INTEGRATION=$(CMAKE_SINGLEPASS) -DTRAINING=$(CMAKE_TRAINING_VALUE) -DTORCHFORT_PATH=$(TORCHFORT_PATH) .. && \
	make -j && \
	cd .. && \
	mv ./build/libastaroth_$(PREC).so libastaroth_$(PREC).so


CMAKE_OPTIONS = -DCPU_BUILD=$(CPU_BUILD) -DSUPPORT_OLD_ARCHITECTURES=$(SUPPORT_OLD_ARCHITECTURES) -DUSE_HEFFTE=$(USE_HEFFTE) -DFFT_ENABLED=$(GPU_FFT_ENABLED) -DMULTIGPU_ENABLED=OFF -DMPI_ENABLED=$(CMAKE_MPI) -DBUILD_SHARED_LIBS=ON -DSINGLEPASS_INTEGRATION=$(CMAKE_SINGLEPASS) -DMAX_THREADS_PER_BLOCK=512 \
                -DDOUBLE_PRECISION=$(CMAKE_PREC) -DBUILD_SAMPLES=OFF -DSUPPRESS_COMPILER_WARNINGS=$(SUPPRESS_COMPILER_WARNINGS) \
		-DVERBOSE=off -DUSE_ONLY_OPTIMIZED_KERNELS=$(USE_ONLY_OPTIMIZED_KERNELS) -DOPTIMIZE_INPUT_PARAMS=$(OPTIMIZE_INPUT_PARAMS) -DOPTIMIZE_MEM_ACCESSES=ON -DOPTIMIZE_ARRAYS=$(FULL_DSL_OPTIMIZATIONS) -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) \
	        -DBUILD_STANDALONE=OFF -DOPTIMIZE_FIELDS=$(TRANSPILATION) -DALLOW_DEAD_VARIABLES=$(FULL_DSL_OPTIMIZATIONS) -DBUILD_UTILS=OFF \
		-DCMAKE_POSITION_INDEPENDENT_CODE=ON -DELIMINATE_CONDITIONALS=$(ELIMINATE_CONDITIONALS) -DDSL_MODULE_FILE=$(DSL_MODULE_FILE) -DINLINING=$(INLINING) -DBUILD_ACM=OFF -DPROFILING_ENABLED=OFF

ifdef GPU_ARCH
CMAKE_OPTIONS += -DGPU_ARCH=$(GPU_ARCH)
endif
#
# External sources are provided with absolute pathnames to avoid any trouble with symbolic links.
#
astaroth_libs: $(PC_HEADERS) DSL/local/mhdsolver.ac
ifeq ($(TRANSPILATION),on)
	@if ! cmp -s $(PENCIL_HOME)/src/astaroth/DSL/$(DSL_MODULE_FILE) DSL/local/$(DSL_MODULE_FILE); then \
		cp $(PENCIL_HOME)/src/astaroth/DSL/$(DSL_MODULE_FILE) DSL/local/; \
	fi
endif
ifeq ($(READ_OVERRIDES),on)
	cp submodule/build/runtime_build/overrides.h submodule/build
endif
	#TP: if PC_nghost.h is found to be missing we turn to the default of Astaroth i.e. STENCIL_ORDER = 6
	@if [ ! -f DSL/local/PC_nghost.h ]; then \
		echo "hostdefine STENCIL_ORDER (6)" > DSL/local/PC_nghost.h; \
	fi

	@ext_srcs="`ls -m $$PENCIL_HOME/src/astaroth/DSL/*/*.h \
                          $$PENCIL_HOME/src/astaroth/DSL/*.h \
		          $$PENCIL_HOME/src/astaroth/submodule/acc-runtime/stdlib/*.h \
			  $$PENCIL_HOME/src/astaroth/submodule/acc-runtime/stdlib/*/*.h \
                          $(CWD)/DSL/local/*.h \
			  |tr "\n" " "|sed -e's/, /;/g'` $(EXT_KERNELS)" && \
	echo "const char* cmake_options = \" -DDSL_EXT_SOURCES=\\\"$$ext_srcs\\\" -DDSL_MODULE_DIR=../../../$(DSL_MODULE_DIR) $(CMAKE_OPTIONS)\";" > cmake_options.h && \
	mkdir -p submodule/build  && cd submodule/build && export LDFLAGS=""  &&  \
	cmake  $(CMAKE_OPTIONS) -DDSL_EXT_SOURCES="$$ext_srcs" -DDSL_MODULE_DIR=../acc-runtime/../../$(DSL_MODULE_DIR) -DRUNTIME_COMPILATION=$(RUNTIME_COMPILATION) -DREAD_OVERRIDES=$(READ_OVERRIDES) $(AC_HOME)
	cd submodule/build && sed -i -e's/SHELL *= *\/bin\/sh/SHELL = \/bin\/bash/' Makefile && make -j -I../../.. VERBOSE=0 # Build the core Astaroth library
	@echo PRECISION=$(CMAKE_PREC)

#        CMAKE_PREFIX_PATH="/opt/nvidia/hpc_sdk/Linux_x86_64/24.1/comm_libs/mpi/lib /opt/nvidia/hpc_sdk/Linux_x86_64/24.1/comm_libs/mpi/include" \
update_astaroth_libs:
#	cd submodule/build && export LDFLAGS="" && \
#	make -I../../.. VERBOSE=0 # Build the core Astaroth library with existing Makefile.
#	@echo PRECISION=$(CMAKE_PREC)

include Makefile.extern

$(PC_HEADERS): .sentinel

.sentinel: ../scripts/phys_modules2c PC_modulesources.h $(MODULESOURCES)
	@rm -f $(PCHEADER_DIR)/PC_modulepars.h $(PCHEADER_DIR)/PC_module*funcs.h $(DSL_MODULE_DIR)/PC_modulepardecs.h 
	@echo '// automatically generated; do not edit!' > $(PCHEADER_DIR)/PC_modulepars.h
	@echo '// automatically generated; do not edit!' > $(DSL_MODULE_DIR)/PC_modulepardecs.h
ifeq ($(TRANSPILATION),on)
	@export DSL_MODULE_DIR=$(DSL_MODULE_DIR); \
        for file in $(MODULESOURCES); do echo processing $$file ...; if [[ $$file =~ eos[_.] ]]; \
	then ../scripts/phys_modules2c_transpilation par $$file equationofstate; \
        else ../scripts/phys_modules2c_transpilation par $$file; fi; done
else
	@export DSL_MODULE_DIR=$(DSL_MODULE_DIR); \
        for file in $(MODULESOURCES); do echo processing $$file ...; if [[ $$file =~ eos[_.] ]]; \
	then ../scripts/phys_modules2c par $$file equationofstate; \
        else ../scripts/phys_modules2c par $$file; fi; done
endif
	@touch .sentinel
#
DSL/local/mhdsolver.ac: .sentinel
ifeq ($(TRANSPILATION),on)
        ifeq ($(GENERATE_DSL_CODE),on)
	  rm -f PC_modulesources_old.h
        endif
        ifeq ($(MODIFY_SOURCE_CODE),on)
	  python3 transpile.py --modify_source_code
        else
	  python3 transpile.py
        endif
	cp PC_modulesources.h PC_modulesources_old.h
else
	@touch DSL/local/mhdsolver.ac
endif

clean:
	@rm -f $(OBJECTS)
	@rm -f *.so .sentinel cmake_options.h
	@if [[ -d submodule/build ]]; then rm -f -r submodule/build; fi # Removes the files generated with cmake for the Astaroth code
	@rm -rf build
	@rm -rf scheduler-build
	@rm -rf include
	rm -f ac_build.log
#2026-May-18/Kishore: not sure if below belong here or in cleann
# stuff generated by fortran-parser
# @Touko, is it possible to get fortran-parser to write all its output in a dedicated subfolder? It seems difficult to manually keep track of all its files like below.
	rm -f res-*.txt res.txt save.txt macros.h ac_declarations.h var_declares.h
# transpilation
	rm -f PC_modulesources_old.h
#
# Deep cleaning, e.g. necessary with compiler change
#
cleann: clean
	@rm -f ../cparam_c.h ../defines_cdata.h ../forcing_c.h ../forcing_pars_c.h ../sub_c.h ../mpicomm_c.h ../boundcond_c.h ../cdata_c.h
	@rm -f $(DSL_MODULE_DIR)/PC_modulepardecs.h
	@rm -f $(PC_HEADER_DIR)/PC_modulepars.h $(PCHEADER_DIR)/PC_module_parfuncs.h
#
#--------------------------------------------------------------------------------------------------------------------------------------------
OPTFLAGS = -O1 
CCFLAGS = $(OPTFLAGS) -Xcompiler -shared -fPIC --debug -I.. -Isubmodule/include -Isubmodule/build -Isubmodule/build/acc-runtime/api -g -G -lm
%.o: %.cc
	$(CPP) $(CCFLAGS) $(ENVIRON) -c $< -o $@


