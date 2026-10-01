!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
! This file declares all the integer tags used to allow variable
! numbers and types of records in varfiles and other datafiles.
!
! ***WARNING*** Fred 09/10/17 If persistent variables are generated in
! your init routine and depend on random seeds take care to force
! independence of processor as the default from start.f90 is processor
! dependent. See e.g. initialize_interstellar
!

! Persistent
integer, parameter :: id_block_PERSISTENT        = 2000

! Random Seeds
integer, parameter :: id_record_RANDOM_SEEDS     = 1     ! float(nseed)
integer, parameter :: id_record_RANDOM_SEEDS2    = 2     ! float(nseed)

! Iteration number
integer, parameter :: id_record_ITERATION_NUMBER = 100   ! int

! Interstellar
! deprecated:
integer, parameter :: id_record_ISM_T_NEXT_OLD   = 250
integer, parameter :: id_record_ISM_POS_NEXT_OLD = 251
integer, parameter :: id_record_ISM_BOLD_MASS    = 252
! currently active:
integer, parameter :: id_record_ISM_T_NEXT_SNI   = 253   ! float
integer, parameter :: id_record_ISM_T_NEXT_SNII  = 254   ! float
integer, parameter :: id_record_ISM_X_CLUSTER    = 255   ! float
integer, parameter :: id_record_ISM_Y_CLUSTER    = 256   ! float
integer, parameter :: id_record_ISM_Z_CLUSTER    = 260   ! float
integer, parameter :: id_record_ISM_T_CLUSTER    = 261   ! float
integer, parameter :: id_record_ISM_TOGGLE_SNI   = 257   ! bool
integer, parameter :: id_record_ISM_TOGGLE_SNII  = 258   ! bool
! deprecated:
integer, parameter :: id_record_ISM_SNRS         = 259
integer, parameter :: id_record_ISM_TOGGLE_OLD   = 1001
integer, parameter :: id_record_ISM_SNRS_OLD     = 1002

! Forcing
integer, parameter :: id_record_FORCING_LOCATION = 270   ! float(3,2)
integer, parameter :: id_record_FORCING_TSFORCE  = 271   ! float
integer, parameter :: id_record_FORCING_TORUS    = 272   ! type

! Hydro
integer, parameter :: id_record_HYDRO_TPHASE     = 280   ! float
integer, parameter :: id_record_HYDRO_PHASE1     = 281   ! float
integer, parameter :: id_record_HYDRO_PHASE2     = 282   ! float
integer, parameter :: id_record_HYDRO_TSFORCE    = 284   ! float
integer, parameter :: id_record_HYDRO_LOCATION   = 285   ! float(3)
integer, parameter :: id_record_HYDRO_AMPL       = 286   ! float
integer, parameter :: id_record_HYDRO_WAVENUMBER = 287   ! float
integer, parameter :: id_record_HYDRO_QVEC_GB    = 288   ! float(3)
integer, parameter :: id_record_HYDRO_AVEC_GB    = 289   ! float(3)

! Magnetic
integer, parameter :: id_record_MAGNETIC_PHASE   = 311   ! float
integer, parameter :: id_record_MAGNETIC_AMPL    = 312   ! float

! Shear
integer, parameter :: id_record_SHEAR_DELTA_Y    = 320   ! float

! Time stepping
integer, parameter :: id_record_TIME_STEP        = 330   ! float
integer, parameter :: id_record_EPS_RKF          = 331   ! float

! special/axionSU2back.f90
integer, parameter :: id_record_SPECIAL_LNKMIN0  = 340   ! float

! special/gravitational_waves_hTXk.f90
integer, parameter :: id_record_DT_GW            = 350   ! float

! special/backreact_infl.f90
integer, parameter :: id_record_LHEATING_ALWAYS  = 360   ! bool
integer, parameter :: id_record_LSOLVE_FOR_PHI   = 361   ! bool

! special/klein_gordon.f90
integer, parameter :: id_record_WALL_VEL         = 370   ! float
integer, parameter :: id_record_WALL_POS         = 371   ! float
