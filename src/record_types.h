!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
! This file declares all the integer tags used to allow variable
! numbers and types of records in varfiles and other datafiles.
!
integer, parameter :: id_block_PERSISTENT        = 2000
!
!
! Integers
!
integer, parameter :: id_record_RANDOM_SEEDS     = 1
!
! Reals
!
integer, parameter :: id_record_T_NEXT_SNI       = 250
integer, parameter :: id_record_BOLD_MASS        = 251
integer, parameter :: id_record_FORCING_LOCATION = 270
integer, parameter :: id_record_FORCING_TSFORCE  = 271
integer, parameter :: id_record_NOHYDRO_TPHASE   = 280
integer, parameter :: id_record_NOHYDRO_PHASE1   = 281
integer, parameter :: id_record_NOHYDRO_PHASE2   = 282
integer, parameter :: id_record_NOHYDRO_TSFORCE  = 284
integer, parameter :: id_record_NOHYDRO_LOCATION = 285
integer, parameter :: id_record_MAGNETIC_PHASE   = 311
integer, parameter :: id_record_MAGNETIC_AMPL    = 312
integer, parameter :: id_record_DELTA_Y          = 320
!
! Other
!
integer, parameter :: id_record_ISM_SN_TOGGLE    = 1001
integer, parameter :: id_record_ISM_SNRS         = 1002
!
