!  $Id: noinitcond_spec.f90,v 1.3 2003-04-04 16:47:16 anders Exp $
!
!  Substitute routines for vortex_solve.f90
!
module Initcond_spec

  implicit none

  contains

!***********************************************************************
    subroutine kepvor(f,xx,yy,zz,b_ell,q_ell,gamma,cs20,hh0)
!
      use Cdata
      use General

      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: q_ell, b_ell, hh0
      real :: gamma,cs20
    endsubroutine kepvor

end module Initcond_spec
