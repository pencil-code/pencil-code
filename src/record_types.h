!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)
!
! This file declares all the integer tags used to allow variable
! numbers and types of records in varfiles and other datafiles.
!
integer, parameter :: id_block_PERSISTANT = 2000 
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
!
! Other
!
integer, parameter :: id_record_INTERSTELLAR_SN_TOGGLE = 1001
integer, parameter :: id_record_INTERSTELLAR_SNRS = 1002
!
