module Global

!
!  A dummy module for the (lucky) case when we need no global variables.
!
  use Cparam

  implicit none

  contains

!***********************************************************************
    subroutine wglobal()
!
!  write global variables
!
    endsubroutine wglobal
!***********************************************************************
    subroutine rglobal()
!
!  read global variables
!
    endsubroutine rglobal
!***********************************************************************

endmodule Global
