ACSOURCES = $(wildcard submodule/src/core/*.cu) 
#submodule/src/standalone/config_loader.cc
ACHEADERS = $(wildcard submodule/*.h) $(wildcard submodule/src/core/kernels/*.cuh)
ACOBJECTS = $(ACSOURCES:.cu=.o) $(ACSOURCES:.cc=.o)

CCSOURCES = loadStore.cc 
OBJECTS = $(CCSOURCES:.cc=.o) 
MAIN_SRC = gpu_astaroth.cc
MAIN_OBJ = gpu_astaroth.o

# Settings for taito-gpu
CUFLAGS = -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 
#-gencode arch=compute_60,code=sm_60 --expt-relaxed-constexpr
OPTFLAGS = -O1 
#-funroll-loops -march=native 
# -gencode arch=compute_60,code=sm_60
#  Original
CCFLAGS = -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. -Isubmodule -Isubmodule/src -Isubmodule/include -Isubmodule/src/standalone -Isubmodule/build -g -G -lm
#CCFLAGS = -std=c++11 -Xcompiler --debug -I. -I../.. -I.. -Isubmodule -Isubmodule/src -Isubmodule/include -Isubmodule/src/standalone -g -G -lm

# MV's adaptation for TIARA cluster. Compiles with: 
#  1) cuda/9.0   3) openmpi/2.0.1_ic16.0            5) hdf5/1.8.16_openmpi_2.0.1_ic16.0 7) gcc/5.3.0
#  2) intel/2016 4) fftw/2.1.5_openmpi_2.0.1_ic16.0 6) cuda/9.1
#CCFLAGS = -ccbin gcc -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. 
##CCFLAGS = -ccbin icpc -Xcompiler -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. 

# Settings for Omer's machine
#CFLAGS = -gencode arch=compute_52,code=sm_52 --shared --relocatable-device-code=true --compiler-options '-fPIC'

#Importing environ variables into preprocessor
ENVIRON = -DMODPRE=${MODULE_PREFIX} -DMODIN=${MODULE_INFIX} -DMODSUF=${MODULE_SUFFIX} -DPENCIL_ASTAROTH

PRECISION =
include PC_modulesources.h
ifeq ($(PRECISION),DOUBLE)
  ENVIRON += -DAC_DOUBLE_PRECISION=1
  PREC=dbl
#  OPTFLAGS += -fsingle-PRECision-constant
  CMAKE_PREC = 'DOUBLE_PRECISION=ON'
else
  ENVIRON += -DAC_DOUBLE_PRECISION=0
  PREC=sgl
  CMAKE_PREC = 'DOUBLE_PRECISION=OFF'
endif

PCHEADER_DIR=.
export CUDA_NVCC_FLAGS='-shared -Xcompiler -fPIC'

astaroth_$(PREC).so: PC_modulesources.h $(MODULESOURCES) $(PCHEADER_DIR)/PC_moduleflags.h $(PC_HEADERS) submodule/build/src/core/libastaroth_core.so  $(MAIN_OBJ) $(OBJECTS)
	nvcc $(CCFLAGS) $(CUFLAGS) -o astaroth_$(PREC).so $(OBJECTS) $(MAIN_OBJ) -L submodule/build/src/core -lastaroth_core

submodule/build/src/core/libastaroth_core.so: $(ACSOURCES) $(ACHEADERS)
	mkdir -p submodule/build && cd submodule/build && export DSL_MODULE=pc_mhd_solver && cmake -DAUTO_OPTIMIZATION=OFF -DMULTIGPU_ENABLED=OFF -DBUILD_DEBUG=OFF -DBUILD_SHARED=ON -DPACKED_DATA_TRANSFERS=ON -D$(CMAKE_PREC) -DBUILD_STANDALONE=OFF -DVERBOSE_PRINTING=OFF .. && make -I../../.. -j # VERBOSE=1 # Build the core Astaroth library

include Makefile.depend
include Makefile.extern

%.o: %.cc
	nvcc $(CCFLAGS) $(ENVIRON) -dc $< -o $@

%.o: %.cu
	nvcc $(CUFLAGS) $(CCFLAGS) $(OPTFLAGS) $(ENVIRON) -dc $< -o $@

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
	@rm -f $(OBJECTS) $(MAIN_OBJ)
	@rm -f *.so .sentinel
	rm -rI submodule/build # Removes the files generated with cmake for the Astaroth code
#	@rm -f $(CHEADERS)

# Deep cleaning, e.g. necessary with compiler change
#
cleann: clean
	@rm -f ../cparam_c.h ../defines_cdata.h ../diagnostics_c.h ../forcing_c.h ../forcing_pars_c.h ../sub_c.h ../mpicomm_c.h ../boundcond_c.h
	@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < $(PCHEADER_DIR)/PC_modulepardecs.h | xargs rm -f
	@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < $(PCHEADER_DIR)/PC_modulepars.h | xargs rm -f

	#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < diagnostics/PC_modulediags.h | xargs rm -f
	#@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < diagnostics/PC_modulediags_init.h | xargs rm -f
