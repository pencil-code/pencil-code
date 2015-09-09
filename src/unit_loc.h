!  -*-f90-*-  (for emacs)    vim:set filetype=fortran:  (for vim)  
!
! Needed for switching between reading input data from external or internal files.
! If you use the preprocessor, remove all ! and set the appropriate compiler flag for preprocessor invocation;
! define f2003 if you want reading from internal files.
! Without preprocessor use, comment out 'integer, parameter :: unit=1' for reading from internal files 
! and 'character(LEN=*), allocatable :: unit' for reading from external files, respectively.
! In any case one of the two declarations must finally be present.
!
!#IFDEF f2003
     character(LEN=:), allocatable :: unit
!#ELSE
!      integer, parameter :: unit=1
!#ENDIF

