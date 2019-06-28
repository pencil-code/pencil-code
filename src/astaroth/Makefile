CUSOURCES = #$(wildcard submodule/src/core/*.cu)  
CCSOURCES = #submodule/src/standalone/config_loader.cc
MAIN_SRC = gpu_astaroth.cc
MAIN_OBJ = gpu_astaroth.o

OBJECTS = $(CUSOURCES:.cu=.o) $(CCSOURCES:.cc=.o) 

# Settings for taito-gpu
CUFLAGS = -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 
#-gencode arch=compute_60,code=sm_60 --expt-relaxed-constexpr
OPTFLAGS = -O1 
#-funroll-loops -march=native 
# -gencode arch=compute_60,code=sm_60
#  Original
CCFLAGS = -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. -Isubmodule -Isubmodule/src -Isubmodule/include -Isubmodule/src/standalone -g -G -lm

# MV's adaptation for TIARA cluster. Compiles with: 
#  1) cuda/9.0   3) openmpi/2.0.1_ic16.0            5) hdf5/1.8.16_openmpi_2.0.1_ic16.0 7) gcc/5.3.0
#  2) intel/2016 4) fftw/2.1.5_openmpi_2.0.1_ic16.0 6) cuda/9.1
#CCFLAGS = -ccbin gcc -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. 
##CCFLAGS = -ccbin icpc -Xcompiler -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. 

# Settings for Omer's machine
#CFLAGS = -gencode arch=compute_52,code=sm_52 --shared --relocatable-device-code=true --compiler-options '-fPIC'

#Importing environ variables into preprocessor
ENVIRON = -D MODPRE=${MODULE_PREFIX} -D MODIN=${MODULE_INFIX} -D MODSUF=${MODULE_SUFFIX} -DPENCIL_ASTAROTH

PRECISION =
include PC_modulesources
ifeq ($(PRECISION),DOUBLE)
  ENVIRON += -DDOUBLE_PRECISION
  PREC=dbl
else ifeq ($(PRECISION),SINGLE)
  PREC=sgl
#  OPTFLAGS += -fsingle-PRECision-constant
else
  ENVIRON += -DDOUBLE_PRECISION
  PREC=dbl
endif

PCHEADER_DIR=.
astaroth_$(PREC).so: PC_modulesources $(MODULESOURCES) $(PCHEADER_DIR)/PC_moduleflags.h $(PC_HEADERS) $(MAIN_OBJ) $(OBJECTS)
	mkdir -p submodule/build && cd submodule/build && cmake -DMULTIGPU_ENABLED=ON -DDOUBLE_PRECISION=ON -DBUILD_STANDALONE=OFF .. && make -j VERBOSE=1 # Build the core Astaroth library
	nvcc $(CCFLAGS) $(CUFLAGS) -o astaroth_$(PREC).so $(OBJECTS) $(MAIN_OBJ)
	cp -f astaroth_$(PREC).so ..

include Makefile.depend

%.o: %.cc
	nvcc $(CCFLAGS) $(ENVIRON) -dc $< -o $@

%.o: %.cu
	nvcc $(CUFLAGS) $(CCFLAGS) $(OPTFLAGS) $(ENVIRON) -dc $< -o $@

$(PC_HEADERS): .sentinel

.sentinel: ../scripts/phys_modules2c PC_modulesources $(MODULESOURCES)
	@rm -f $(PCHEADER_DIR)/PC_modulepars.h $(PCHEADER_DIR)/PC_modulepardecs.h diagnostics/PC_modulediags.h
	@echo '// automatically generated; do not edit!' > $(PCHEADER_DIR)/PC_modulepars.h
	@echo '// automatically generated; do not edit!' > $(PCHEADER_DIR)/PC_modulepardecs.h
	@echo '// automatically generated; do not edit!' > diagnostics/PC_modulediags.h
	@for file in $(MODULESOURCES); do echo processing $$file ...; if [[ $$file =~ eos[_.] ]]; \
	then ../scripts/phys_modules2c par $$file equationofstate; ../scripts/phys_modules2c diag $$file equationofstate;\
        else ../scripts/phys_modules2c par $$file; ../scripts/phys_modules2c diag $$file; fi; done
	@rm -f diagnostics/PC_modulediags_init.h
	@sed -e's/_diags_/_diags_init_/' <  diagnostics/PC_modulediags.h > diagnostics/PC_modulediags_init.h
	@touch .sentinel

clean:
	@rm -f $(OBJECTS) $(MAIN_OBJ)
	@rm -f *.so .sentinel
	rm -rI submodule/build # Removes the files generated with cmake for the Astaroth code
#	@rm -f $(CHEADERS)

# Deep cleaning, e.g. necessary with compiler change
#
cleann: clean
	@rm -f ../cparam_c.h ../defines_cdata.h ../diagnostics_c.h ../forcing_c.h ../forcing_pars_c.h ../sub_c.h
	@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < diagnostics/PC_modulediags.h | xargs rm -f
	@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < diagnostics/PC_modulediags_init.h | xargs rm -f
	@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < $(PCHEADER_DIR)/PC_modulepardecs.h | xargs rm -f
	@sed -e'/\/\// d' -e's/#include *//' -e's/"//g' < $(PCHEADER_DIR)/PC_modulepars.h | xargs rm -f
