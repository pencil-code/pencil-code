! $Id: noglobal.f90,v 1.2 2002-06-01 02:56:21 brandenb Exp $

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
