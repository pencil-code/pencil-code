!  $Id: noinitcond_spec.f90,v 1.2 2003-04-04 16:31:45 anders Exp $
!
!  Substitute routines for vortex_solve.f90
!
module Initcond_spec

  implicit none

  contains

!***********************************************************************
    subroutine kepvor(f,xx,yy,zz,a_ell,b_ell,xi0,gamma,cs20,hh0)
!
      use Cdata
      use General

      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: a_ell, b_ell, xi0, hh0
      real :: gamma,cs20
    endsubroutine kepvor

end module Initcond_spec
