!  $Id: noinitcond_spec.f90,v 1.1 2003-03-27 10:51:39 brandenb Exp $
!
!  Substitute routines for vortex_solve.f90
!
module Initcond_spec

  implicit none

  contains

!***********************************************************************
    subroutine special_ic(f,xx,yy,zz,a_ell,b_ell,xi0,gamma,cs20)
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: a_ell, b_ell, xi0
      real :: gamma,cs20
    endsubroutine special_ic

end module Initcond_spec
