module Forcing

!! Dummy module for Naveir-Stokes equation without forcing

  use Cdata

  implicit none

!   real, dimension (mx,my,mz,3) :: fforce=0   !(forcing function)

  contains

!***********************************************************************
    subroutine addforce(df)
!
      use Cdata
!
!  add forcing in timestep()
!
      real, dimension (mx,my,mz,mvar) :: df
!
    endsubroutine addforce
!***********************************************************************
    subroutine forcing1
!
!  forcing function, using a set of precomputed wavevectors
!
    endsubroutine forcing1
!***********************************************************************
    subroutine forcing2
!
!  helical forcing function, using a set of precomputed wavevectors
!
    endsubroutine forcing2
!***********************************************************************

endmodule Forcing
