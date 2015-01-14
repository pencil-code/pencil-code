! Needed for switching between reading input data from external or internal files.
! If you use the preprocessor, remove all ! and set the appropriate compiler flag for preprocessor invocation;
! define f2003 if you want reading from internal files.
! Without preprocessor use, comment out 'integer :: unit' for reading from internal files 
! and 'character(LEN=*) :: unit' for reading from external files.
! In any case one of the two declarations must finally be present.
! 
!#IFDEF f2003
!  character(LEN=*) :: unit
!#ELSE
  integer :: unit
!#ENDIF
  intent(IN) :: unit
