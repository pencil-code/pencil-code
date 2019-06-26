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
integer, parameter :: id_record_RANDOM_SEEDS     = 1

! Interstellar
! deprecated:
integer, parameter :: id_record_ISM_T_NEXT_OLD   = 250
integer, parameter :: id_record_ISM_POS_NEXT_OLD = 251
integer, parameter :: id_record_ISM_BOLD_MASS    = 252
! currently active:
integer, parameter :: id_record_ISM_T_NEXT_SNI   = 253
integer, parameter :: id_record_ISM_T_NEXT_SNII  = 254
integer, parameter :: id_record_ISM_X_CLUSTER    = 255
integer, parameter :: id_record_ISM_Y_CLUSTER    = 256
integer, parameter :: id_record_ISM_Z_CLUSTER    = 260
integer, parameter :: id_record_ISM_T_CLUSTER    = 261
integer, parameter :: id_record_ISM_TOGGLE_SNI   = 257
integer, parameter :: id_record_ISM_TOGGLE_SNII  = 258
! deprecated:
integer, parameter :: id_record_ISM_SNRS         = 259
integer, parameter :: id_record_ISM_TOGGLE_OLD   = 1001
integer, parameter :: id_record_ISM_SNRS_OLD     = 1002

! Forcing
integer, parameter :: id_record_FORCING_LOCATION = 270
integer, parameter :: id_record_FORCING_TSFORCE  = 271

! Hydro
integer, parameter :: id_record_HYDRO_TPHASE     = 280
integer, parameter :: id_record_HYDRO_PHASE1     = 281
integer, parameter :: id_record_HYDRO_PHASE2     = 282
integer, parameter :: id_record_HYDRO_TSFORCE    = 284
integer, parameter :: id_record_HYDRO_LOCATION   = 285
integer, parameter :: id_record_HYDRO_AMPL       = 286
integer, parameter :: id_record_HYDRO_WAVENUMBER = 287

! Magnetic
integer, parameter :: id_record_MAGNETIC_PHASE   = 311
integer, parameter :: id_record_MAGNETIC_AMPL    = 312

! Shear
integer, parameter :: id_record_SHEAR_DELTA_Y    = 320

! Special
integer, parameter :: id_record_SPECIAL_ILOAD    = 330

