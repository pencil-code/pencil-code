SOURCES = gpu_astaroth.cc loadStore.cc 
OBJECTS = $(SOURCES:.cc=.o) 

# Settings for taito-gpu
CUFLAGS = -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 
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

# Settings for Omer's machine
#CFLAGS = -gencode arch=compute_52,code=sm_52 --shared --relocatable-device-code=true --compiler-options '-fPIC'

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
ENVIRON += -DPACKED_DATA_TRANSFERS=ON

CCFLAGS = $(OPTFLAGS) -std=c++11 -shared -Xcompiler -fPIC --debug -I.. -Isubmodule/include -Isubmodule/build -g -G -lm
#CCFLAGS = -std=c++11 -Xcompiler --debug -I. -I../.. -I.. -Isubmodule -Isubmodule/src -Isubmodule/include -Isubmodule/src/standalone -g -G -lm
CCFLAGS = $(OPTFLAGS) -shared -Xcompiler -fPIC --debug -I.. -Isubmodule/include -Isubmodule/build -Isubmodule/build/acc-runtime/api -g -G -lm

PCHEADER_DIR=.
#export CUDA_NVCC_FLAGS='-shared -Xcompiler -fPIC'
#export CUDA_NVCC_FLAGS='-shared -fPIC'

astaroth_$(PREC).so: PC_modulesources.h $(MODULESOURCES) $(PCHEADER_DIR)/PC_moduleflags.h $(PC_HEADERS) astaroth_libs $(OBJECTS)
	nvcc $(CCFLAGS) $(CUFLAGS) -o astaroth_$(PREC).so $(OBJECTS) -L submodule/build/src/core -L submodule/build/src/core/kernels -L submodule/build/src/utils -lastaroth_core -lkernels -lastaroth_utils -Isubmodule/include

astaroth_libs:
	mkdir -p submodule/build && cd submodule/build && \
	cmake -DMULTIGPU_ENABLED=OFF -DBUILD_DEBUG=OFF -DBUILD_SHARED=ON -DAC_VERBOSE=OFF \
              -DPACKED_DATA_TRANSFERS=ON -D$(CMAKE_PREC) -DBUILD_STANDALONE=OFF -DBUILD_SAMPLES=OFF .. \
	&& make -I../../.. VERBOSE=0 # Build the core Astaroth library

#-DCMAKE_CXX_COMPILER=/appl/spack/install-tree/gcc-4.8.5/gcc-8.3.0-qzmzn5/bin/c++ .. \
#export DSL_MODULE=pc_mhd_solver && \

include Makefile.depend
include Makefile.extern

%.o: %.cc
	nvcc $(CCFLAGS) $(ENVIRON) -dc $< -o $@

$(PC_HEADERS): .sentinel

.sentinel: ../scripts/phys_modules2c PC_modulesources.h $(MODULESOURCES)
	@rm -f $(PCHEADER_DIR)/PC_modulepars.h $(PCHEADER_DIR)/PC_modulepardecs.h 
#diagnostics/PC_modulediags.h
	@echo '// automatically generated; do not edit!' > $(PCHEADER_DIR)/PC_modulepars.h
	@echo '// automatically generated; do not edit!' > $(PCHEADER_DIR)/PC_modulepardecs.h
#@echo '// automatically generated; do not edit!' > diagnostics/PC_modulediags.h
	@for file in $(MODULESOURCES); do echo processing $$file ...; if [[ $$file =~ eos[_.] ]]; \
	then ../scripts/phys_modules2c par $$file equationofstate; \
        else ../scripts/phys_modules2c par $$file; fi; done
#@rm -f diagnostics/PC_modulediags_init.h
#@sed -e's/_diags_/_diags_init_/' <  diagnostics/PC_modulediags.h > diagnostics/PC_modulediags_init.h
	@touch .sentinel

#        @for file in $(MODULESOURCES); do echo processing $$file ...; if [[ $$file =~ eos[_.] ]]; \
#        then ../scripts/phys_modules2c par $$file equationofstate; ../scripts/phys_modules2c diag $$file equationofstate;\
#        else ../scripts/phys_modules2c par $$file; ../scripts/phys_modules2c diag $$file; fi; done

clean:
	@rm -f $(OBJECTS)
	@rm -f *.so .sentinel
	rm -rI submodule/build # Removes the files generated with cmake for the Astaroth code

# Deep cleaning, e.g. necessary with compiler change
#
cleann: clean
	@rm -f ../cparam_c.h ../defines_cdata.h ../diagnostics_c.h ../forcing_c.h ../forcing_pars_c.h ../sub_c.h ../mpicomm_c.h ../boundcond_c.h
	@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < $(PCHEADER_DIR)/PC_modulepardecs.h | xargs rm -f
	@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < $(PCHEADER_DIR)/PC_modulepars.h | xargs rm -f

	#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < diagnostics/PC_modulediags.h | xargs rm -f
	#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < diagnostics/PC_modulediags_init.h | xargs rm -f
