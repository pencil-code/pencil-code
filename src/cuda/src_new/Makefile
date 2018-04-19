CUSOURCES = $(wildcard gpu/cuda/*.cu) $(wildcard gpu/cuda/generic/*.cu) $(wildcard gpu/cuda/core/*.cu)
CCSOURCES = $(wildcard gpu/*.cc) $(wildcard common/*.cc) $(wildcard diagnostics/*.cc)
MAIN_SRC = gpu_astaroth.cc
MAIN_OBJ = gpu_astaroth.o

OBJECTS = $(CUSOURCES:.cu=.o) $(CCSOURCES:.cc=.o)

# Settings for taito-gpu
CUFLAGS = -gencode arch=compute_35,code=sm_35 -gencode arch=compute_37,code=sm_37 -gencode arch=compute_60,code=sm_60 --expt-relaxed-constexpr
OPTFLAGS = -O3 
#-funroll-loops -march=native 
# -gencode arch=compute_60,code=sm_60
CCFLAGS = -std=c++11 -shared -Xcompiler -fPIC --debug -I. -I../.. -I.. -g -G
# Settings for Omer's machine
#CFLAGS = -gencode arch=compute_52,code=sm_52 --shared --relocatable-device-code=true --compiler-options '-fPIC'

#Importing environ variables into preprocessor
ENVIRON = -D MODPRE=${MODULE_PREFIX} -D MODIN=${MODULE_INFIX} -D MODSUF=${MODULE_SUFFIX} -DGPU_ASTAROTH

DOUBLE_PRECISION=
include PC_modulesources
ifeq ($(DOUBLE_PRECISION),DOUBLE_PRECISION)
  ENVIRON += -DDOUBLE_PRECISION
  PREC=dbl
else
  PREC=sgl
#  OPTFLAGS += -fsingle-PRECision-constant
endif

astaroth_$(PREC).so: PC_modulesources $(MODULESOURCES) common/PC_moduleflags.h $(PC_HEADERS) $(MAIN_OBJ) $(OBJECTS)
	nvcc $(CCFLAGS) $(CUFLAGS) -o astaroth_$(PREC).so $(OBJECTS) $(MAIN_OBJ)
	cp -f astaroth_$(PREC).so ../..

include Makefile.depend

%.o: %.cc
	nvcc $(CCFLAGS) $(ENVIRON) -dc $< -o $@

%.o: %.cu
	nvcc $(CUFLAGS) $(CCFLAGS) $(OPTFLAGS) $(ENVIRON) -dc $< -o $@

$(PC_HEADERS): .sentinel

.sentinel: ../../scripts/phys_modules2c $(MODULESOURCES)
	@rm -f common/PC_modulepars.h common/PC_modulepardecs.h diagnostics/PC_modulediags.h
	@echo '// automatically generated; do not edit!' > common/PC_modulepars.h
	@echo '// automatically generated; do not edit!' > common/PC_modulepardecs.h
	@echo '// automatically generated; do not edit!' > diagnostics/PC_modulediags.h
	for file in $(MODULESOURCES); do echo file=$$file; if [[ $$file =~ eos[_.] ]]; \
	then echo $$file; ../../scripts/phys_modules2c par $$file equationofstate; ../../scripts/phys_modules2c diag $$file equationofstate;\
        else ../../scripts/phys_modules2c par $$file; ../../scripts/phys_modules2c diag $$file; fi; done
	rm -f diagnostics/PC_modulediags_init.h
	sed -e's/_diags_/_diags_init_/' <  diagnostics/PC_modulediags.h > diagnostics/PC_modulediags_init.h
	@touch .sentinel

clean:
	@rm -f $(OBJECTS) $(MAIN_OBJ)
	@rm -f *.so .sentinel
#	@rm -f $(CHEADERS)
