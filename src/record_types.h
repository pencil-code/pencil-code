!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
! This file declares all the integer tags used to allow variable
! numbers and types of records in varfiles and other datafiles.

! Persistent
integer, parameter :: id_block_PERSISTENT        = 2000

! Random Seeds
integer, parameter :: id_record_RANDOM_SEEDS     = 1

! Interstellar
integer, parameter :: id_record_T_NEXT_SNI       = 250
integer, parameter :: id_record_POS_NEXT_SNII    = 251
integer, parameter :: id_record_BOLD_MASS        = 252

! Forcing
integer, parameter :: id_record_FORCING_LOCATION = 270
integer, parameter :: id_record_FORCING_TSFORCE  = 271

! Nohydro
integer, parameter :: id_record_NOHYDRO_TPHASE   = 280
integer, parameter :: id_record_NOHYDRO_PHASE1   = 281
integer, parameter :: id_record_NOHYDRO_PHASE2   = 282
integer, parameter :: id_record_NOHYDRO_TSFORCE  = 284
integer, parameter :: id_record_NOHYDRO_LOCATION = 285

! Magnetic
integer, parameter :: id_record_MAGNETIC_PHASE   = 311
integer, parameter :: id_record_MAGNETIC_AMPL    = 312

! Shear
integer, parameter :: id_record_DELTA_Y          = 320

! Interstellar
integer, parameter :: id_record_ISM_SN_TOGGLE    = 1001
integer, parameter :: id_record_ISM_SNRS         = 1002

