!  $Id: noinitcond_spec.f90,v 1.6 2003-06-16 04:41:10 brandenb Exp $
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

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: q_ell, b_ell, hh0
      real :: gamma,cs20
      if(ip==0) print*,f,xx,yy,zz,b_ell,q_ell,gamma,cs20,hh0 ! compiler quiet
    endsubroutine kepvor
!***********************************************************************
    subroutine enthblob(f,xx,yy,zz,b_ell,q_ell,gamma,cs20,hh0)
!
      use Cdata
      use General

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
      real :: q_ell, b_ell, hh0
      real :: gamma,cs20
      if(ip==0) print*,f,xx,yy,zz,b_ell,q_ell,gamma,cs20,hh0 ! compiler quiet
    endsubroutine enthblob

end module Initcond_spec
