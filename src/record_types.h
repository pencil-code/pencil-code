!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
! This file declares all the integer tags used to allow variable
! numbers and types of records in varfiles and other datafiles.
!
integer, parameter :: id_block_PERSISTENT = 2000 
!
!
! Integers 
!
integer, parameter :: id_record_RANDOM_SEEDS = 1
!
! Reals
!
integer, parameter :: id_record_T_NEXT_SN    = 250 
integer, parameter :: id_record_FORCING_LOCATION = 270 
integer, parameter :: id_record_FORCING_TSFORCE  = 271 
integer, parameter :: id_record_MAGNETIC_BELTRAMI_PHASE = 311 
integer, parameter :: id_record_MAGNETIC_BELTRAMI_AMPL = 312 
!
! Other
!
integer, parameter :: id_record_ISM_SN_TOGGLE = 1001
integer, parameter :: id_record_ISM_SNRS = 1002
!
