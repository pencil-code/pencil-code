DIR_COALA := ../../src/Coala

AR := ar
ARFLAGS := rcs

BUILD_DIR := build
MOD_DIR := $(BUILD_DIR)/mod
OBJ_DIR := $(BUILD_DIR)/obj

LIB := $(BUILD_DIR)/libcoala.a

SOURCES := \
	$(DIR_COALA)/precision.f90 \
	$(DIR_COALA)/progress_bar.f90 \
	$(DIR_COALA)/tension_module.f90 \
	$(DIR_COALA)/init_massgrid.f90 \
	$(DIR_COALA)/polynomials_legendre.f90 \
	$(DIR_COALA)/limiter.f90 \
	$(DIR_COALA)/reconstruction_g.f90 \
	$(DIR_COALA)/GQ_legendre_nodes_weights.f90 \
	$(DIR_COALA)/coag_term/coagflux_functions_GQ.f90 \
	$(DIR_COALA)/coag_term/coagintflux_functions_GQ.f90 \
	$(DIR_COALA)/collision_kernels.f90 \
	$(DIR_COALA)/functions_tabflux_tabintflux.f90 \
	$(DIR_COALA)/functions_flux_intflux.f90 \
	$(DIR_COALA)/solver_DG.f90 \
	$(DIR_COALA)/initial_condition_L2proj_GQ.f90 \
	$(DIR_COALA)/compute_solver_coag.f90 \
	$(DIR_COALA)/interface_coala_coag.f90

OBJECTS := $(patsubst $(DIR_COALA)/%.f90,$(OBJ_DIR)/%.o,$(SOURCES))

.PHONY: all clean

all: $(LIB)

$(LIB): $(OBJECTS)
	@mkdir -p $(dir $@)
	$(AR) $(ARFLAGS) $@ $^

# Fortran module dependencies

$(OBJ_DIR)/progress_bar.o: \
	$(OBJ_DIR)/precision.o

$(OBJ_DIR)/init_massgrid.o: \
	$(OBJ_DIR)/precision.o

$(OBJ_DIR)/polynomials_legendre.o: \
	$(OBJ_DIR)/precision.o

$(OBJ_DIR)/GQ_legendre_nodes_weights.o: \
	$(OBJ_DIR)/precision.o

$(OBJ_DIR)/collision_kernels.o: \
	$(OBJ_DIR)/precision.o

$(OBJ_DIR)/reconstruction_g.o: \
	$(OBJ_DIR)/precision.o \
	$(OBJ_DIR)/polynomials_legendre.o

$(OBJ_DIR)/coag_term/coagflux_functions_GQ.o: \
	$(OBJ_DIR)/precision.o \
	$(OBJ_DIR)/collision_kernels.o \
	$(OBJ_DIR)/polynomials_legendre.o

$(OBJ_DIR)/coag_term/coagintflux_functions_GQ.o: \
	$(OBJ_DIR)/precision.o \
	$(OBJ_DIR)/collision_kernels.o \
	$(OBJ_DIR)/polynomials_legendre.o

$(OBJ_DIR)/functions_tabflux_tabintflux.o: \
	$(OBJ_DIR)/precision.o \
	$(OBJ_DIR)/progress_bar.o \
	$(OBJ_DIR)/coag_term/coagflux_functions_GQ.o \
	$(OBJ_DIR)/coag_term/coagintflux_functions_GQ.o

$(OBJ_DIR)/functions_flux_intflux.o: \
	$(OBJ_DIR)/precision.o

$(OBJ_DIR)/limiter.o: \
	$(OBJ_DIR)/precision.o \
	$(OBJ_DIR)/polynomials_legendre.o \
	$(OBJ_DIR)/reconstruction_g.o \
	$(OBJ_DIR)/functions_flux_intflux.o

$(OBJ_DIR)/solver_DG.o: \
	$(OBJ_DIR)/precision.o \
	$(OBJ_DIR)/functions_flux_intflux.o \
	$(OBJ_DIR)/polynomials_legendre.o \
	$(OBJ_DIR)/limiter.o

$(OBJ_DIR)/compute_solver_coag.o: \
	$(OBJ_DIR)/precision.o \
	$(OBJ_DIR)/solver_DG.o \
	$(OBJ_DIR)/progress_bar.o

$(OBJ_DIR)/initial_condition_L2proj_GQ.o: \
	$(OBJ_DIR)/functions_flux_intflux.o \
	$(OBJ_DIR)/precision.o \
	$(OBJ_DIR)/polynomials_legendre.o \
	$(OBJ_DIR)/tension_module.o

$(OBJ_DIR)/interface_coala_coag.o: \
	$(OBJ_DIR)/precision.o \
	$(OBJ_DIR)/initial_condition_L2proj_GQ.o \
	$(OBJ_DIR)/limiter.o \
	$(OBJ_DIR)/compute_solver_coag.o

$(OBJ_DIR)/%.o: $(DIR_COALA)/%.f90
	@mkdir -p $(dir $@) $(MOD_DIR)
	$(FC) $(FFLAGS) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)
