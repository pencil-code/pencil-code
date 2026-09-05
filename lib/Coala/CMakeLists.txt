cmake_minimum_required(VERSION 3.20)

project(COALA LANGUAGES Fortran)

set(DIR_COALA ../../src/Coala)

add_library(coala STATIC
   ${DIR_COALA}/precision.f90
   ${DIR_COALA}/progress_bar.f90
   ${DIR_COALA}/tension_module.f90
   ${DIR_COALA}/init_massgrid.f90
   ${DIR_COALA}/polynomials_legendre.f90
   ${DIR_COALA}/limiter.f90
   ${DIR_COALA}/reconstruction_g.f90
   ${DIR_COALA}/GQ_legendre_nodes_weights.f90
   ${DIR_COALA}/coag_term/coagflux_functions_GQ.f90
   ${DIR_COALA}/coag_term/coagintflux_functions_GQ.f90
   ${DIR_COALA}/collision_kernels.f90
   ${DIR_COALA}/functions_tabflux_tabintflux.f90
   ${DIR_COALA}/functions_flux_intflux.f90
   ${DIR_COALA}/solver_DG.f90
   ${DIR_COALA}/initial_condition_L2proj_GQ.f90
   ${DIR_COALA}/compute_solver_coag.f90
   ${DIR_COALA}/interface_coala_coag.f90
)

set_target_properties(coala PROPERTIES
    Fortran_MODULE_DIRECTORY
        "${CMAKE_BINARY_DIR}/mod"
)
