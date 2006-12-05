! $Id: deriv.f90,v 1.35 2006-12-05 19:32:16 ajohan Exp $

module Deriv

  use Messages

  implicit none

  private

  public :: der, der2, der3, der4, der5, der6, derij, der5i1j
  public :: der6_other
  public :: der_upwind1st

!debug  integer, parameter :: icount_der   = 1         !DERCOUNT
!debug  integer, parameter :: icount_der2  = 2         !DERCOUNT
!debug  integer, parameter :: icount_der4  = 3         !DERCOUNT
!debug  integer, parameter :: icount_der5  = 4         !DERCOUNT
!debug  integer, parameter :: icount_der6  = 5         !DERCOUNT
!debug  integer, parameter :: icount_derij = 6         !DERCOUNT
!debug  integer, parameter :: icount_der_upwind1st = 7 !DERCOUNT
!debug  integer, parameter :: icount_der_other = 8     !DERCOUNT

  interface der                 ! Overload the der function
    module procedure der_main   ! derivative of an 'mvar' variable
    module procedure der_other  ! derivative of another field
  endinterface

  contains

!***********************************************************************
    subroutine der_main(f,k,df,j)
!
!  calculate derivative df_k/dx_j
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!   1-oct-97/axel: coded
!  18-jul-98/axel: corrected mx -> my and mx -> mz in all y and z ders
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: j,k
!
!debug      if (loptimise_ders) der_call_count(k,icount_der,j,1) = & !DERCOUNT
!debug                            der_call_count(k,icount_der,j,1)+1 !DERCOUNT
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=(1./60)*dx_1(l1:l2)
          df=fac*(+ 45.0*(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k)) &
                  -  9.0*(f(l1+2:l2+2,m,n,k)-f(l1-2:l2-2,m,n,k)) &
                  +      (f(l1+3:l2+3,m,n,k)-f(l1-3:l2-3,m,n,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=(1./60)*dy_1(m)
          df=fac*(+ 45.0*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                  -  9.0*(f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) &
                  +      (f(l1:l2,m+3,n,k)-f(l1:l2,m-3,n,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=(1./60)*dz_1(n)
          df=fac*(+ 45.0*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k)) &
                  -  9.0*(f(l1:l2,m,n+2,k)-f(l1:l2,m,n-2,k)) &
                  +      (f(l1:l2,m,n+3,k)-f(l1:l2,m,n-3,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_main
!***********************************************************************
    subroutine der_other(f,df,j)
!
!  Along one pencil in NON f variable
!  calculate derivative of a scalar, get scalar
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!   26-nov-02/tony: coded - duplicate der_main but without k subscript
!                           then overload the der interface.
!   25-jun-04/tobi+wolf: adapted for non-equidistant grids

      use Cdata
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df,fac
      integer :: j
!
!debug      if (loptimise_ders) der_call_count(1,icount_der_other,j,1) = &
!debug                          der_call_count(1,icount_der_other,j,1) + 1
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=(1./60)*dx_1(l1:l2)
          df=fac*(+ 45.0*(f(l1+1:l2+1,m,n)-f(l1-1:l2-1,m,n)) &
                  -  9.0*(f(l1+2:l2+2,m,n)-f(l1-2:l2-2,m,n)) &
                  +      (f(l1+3:l2+3,m,n)-f(l1-3:l2-3,m,n)))
        else
          df=0.
          if (ip<=5) print*, 'der_other: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=(1./60)*dy_1(m)
          df=fac*(+ 45.0*(f(l1:l2,m+1,n)-f(l1:l2,m-1,n)) &
                  -  9.0*(f(l1:l2,m+2,n)-f(l1:l2,m-2,n)) &
                  +      (f(l1:l2,m+3,n)-f(l1:l2,m-3,n)))
        else
          df=0.
          if (ip<=5) print*, 'der_other: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=(1./60)*dz_1(n)
          df=fac*(+ 45.0*(f(l1:l2,m,n+1)-f(l1:l2,m,n-1)) &
                  -  9.0*(f(l1:l2,m,n+2)-f(l1:l2,m,n-2)) &
                  +      (f(l1:l2,m,n+3)-f(l1:l2,m,n-3)))
        else
          df=0.
          if (ip<=5) print*, 'der_other: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_other
!***********************************************************************
    subroutine der2(f,k,df2,j)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df2,fac,df
      integer :: j,k
!
!debug      if (loptimise_ders) der_call_count(k,icount_der2,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der2,j,1) + 1 !DERCOUNT
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=(1./180)*dx_1(l1:l2)**2
          df2=fac*(-490.0*f(l1  :l2  ,m,n,k) &
                   +270.0*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                   - 27.0*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                   +  2.0*(f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dx_tilde(l1:l2)*df
          endif
        else
          df2=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=(1./180)*dy_1(m)**2
          df2=fac*(-490.0*f(l1:l2,m  ,n,k) &
                   +270.0*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                   - 27.0*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                   +  2.0*(f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dy_tilde(m)*df
          endif
        else
          df2=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=(1./180)*dz_1(n)**2
          df2=fac*(-490.0*f(l1:l2,m,n  ,k) &
                   +270.0*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                   - 27.0*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                   +  2.0*(f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dz_tilde(n)*df
          endif
        else
          df2=0.
        endif
      endif

!
    endsubroutine der2
!***********************************************************************
    subroutine der3(f,k,df,j,ignoredx)
!
!  Calculate 3rd derivative of a scalar, get scalar
!
!  10-feb-06/anders: adapted from der5
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: j,k
      logical, optional :: ignoredx
      logical :: igndx
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der5,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der5,j,1) + 1 !DERCOUNT
!
      if (present(ignoredx)) then
        igndx = ignoredx
      else
        igndx = .false.
      endif

      if (.not. lequidist(j)) &
          call fatal_error('der3','NOT IMPLEMENTED for no equidistant grid')
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=(1.0/8)
          else
            fac=(1.0/8)*1./dx**3
          endif
          df=fac*(- 13.0*(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k)) &
                  +  8.0*(f(l1+2:l2+2,m,n,k)-f(l1-2:l2-2,m,n,k)) &
                  -  1.0*(f(l1+3:l2+3,m,n,k)-f(l1-3:l2-3,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=(1.0/8)
          else
            fac=(1.0/8)*1./dy**3
          endif
          df=fac*(- 13.0*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                  +  8.0*(f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) &
                  -  1.0*(f(l1:l2,m+3,n,k)-f(l1:l2,m-3,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=(1.0/8)
          else
            fac=(1.0/8)*1./dz**3
          endif
          df=fac*(- 13.0*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k)) &
                  +  8.0*(f(l1:l2,m,n+2,k)-f(l1:l2,m,n-2,k)) &
                  -  1.0*(f(l1:l2,m,n+3,k)-f(l1:l2,m,n-3,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der3
!***********************************************************************
    subroutine der4(f,k,df,j,ignoredx,upwind)
!
!  Calculate 4th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!
!   8-jul-02/wolf: coded
!   9-dec-03/nils: adapted from der6
!  10-feb-06/anders: corrected sign and factor
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df
      real :: fac
      integer :: j,k
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der4,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der4,j,1) + 1 !DERCOUNT
!
      if (.not. lequidist(j)) then
        call fatal_error('der4','NOT IMPLEMENTED for no equidistant grid')
      endif

      if (present(ignoredx)) then
        igndx = ignoredx
      else
        igndx = .false.
      endif
      if (present(upwind)) then
        upwnd = upwind
        call warning('der4','upwinding not implemented')
      else
        upwnd = .false.
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=(1.0/6)
          else
            fac=(1.0/6)*1/dx**4
          endif
          df=fac*(+ 56.0* f(l1:l2,m,n,k) &
                  - 39.0*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  + 12.0*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                  -      (f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=(1.0/6)
          else
            fac=(1.0/6)*1/dy**4
          endif
          df=fac*(+ 56.0* f(l1:l2,m  ,n,k) &
                  - 39.0*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  + 12.0*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                  -      (f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=(1.0/6)
          else
            fac=(1.0/6)*1/dz*4
          endif
          df=fac*(+ 56.0* f(l1:l2,m,n  ,k) &
                  - 39.0*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                  + 12.0*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                  -      (f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der4
!***********************************************************************
    subroutine der5(f,k,df,j,ignoredx)
!
!  Calculate 5th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!
!  29-oct-04/anders: adapted from der6
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: j,k
      logical, optional :: ignoredx
      logical :: igndx
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der5,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der5,j,1) + 1 !DERCOUNT
!
      if (present(ignoredx)) then
        igndx = ignoredx
      else
        igndx = .false.
      endif

      if (.not. lequidist(j)) &
          call fatal_error('der5','NOT IMPLEMENTED for no equidistant grid')
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=1.0
          else
            fac=1/dx**5
          endif
          df=fac*(+  2.5*(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k)) &
                  -  2.0*(f(l1+2:l2+2,m,n,k)-f(l1-2:l2-2,m,n,k)) &
                  +  0.5*(f(l1+3:l2+3,m,n,k)-f(l1-3:l2-3,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=1.0
          else
            fac=1/dy**5
          endif
          df=fac*(+  2.5*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                  -  2.0*(f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) &
                  +  0.5*(f(l1:l2,m+3,n,k)-f(l1:l2,m-3,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=1.0
          else
            fac=1/dz**5
          endif
          df=fac*(+  2.5*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k)) &
                  -  2.0*(f(l1:l2,m,n+2,k)-f(l1:l2,m,n-2,k)) &
                  +  0.5*(f(l1:l2,m,n+3,k)-f(l1:l2,m,n-3,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der5
!***********************************************************************
    subroutine der6(f,k,df,j,ignoredx,upwind)
!
!  Calculate 6th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!    The optional flag UPWIND is a variant thereof, which calculates
!  D^(6)*dx^5/60, which is the upwind correction of centered derivatives.
!
!   8-jul-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: j,k
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der6,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der6,j,1) + 1 !DERCOUNT
!
      if (present(ignoredx)) then
        igndx = ignoredx
      else
        igndx = .false.
      endif
      if (present(upwind)) then
        upwnd = upwind
      else
        upwnd = .false.
        if (.not. lequidist(j)) then
          call fatal_error('der6','NOT IMPLEMENTED for no equidistant grid')
        endif
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=1.0
          else if (upwnd) then
            fac=(1.0/60)*dx_1(l1:l2)
          else
            fac=1/dx**6
          endif
          df=fac*(- 20.0* f(l1:l2,m,n,k) &
                  + 15.0*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  -  6.0*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                  +      (f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=1.0
          else if (upwnd) then
            fac=(1.0/60)*dy_1(m)
          else
            fac=1/dy**6
          endif
          df=fac*(- 20.0* f(l1:l2,m  ,n,k) &
                  + 15.0*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  -  6.0*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                  +      (f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=1.
          else if (upwnd) then
            fac=(1.0/60)*dz_1(n)
          else
            fac=1/dz**6
          endif
          df=fac*(- 20.0* f(l1:l2,m,n  ,k) &
                  + 15.0*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                  -  6.0*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                  +      (f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der6
!***********************************************************************
    subroutine der6_other(f,df,j,ignoredx,upwind)
!
!  Calculate 6th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!    The optional flag UPWIND is a variant thereof, which calculates
!  D^(6)*dx^5/60, which is the upwind correction of centered derivatives.
!
!   8-jul-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df,fac
      integer :: j
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
!
      intent(in)  :: f,j,ignoredx
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der6,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der6,j,1) + 1 !DERCOUNT
!
      if (present(ignoredx)) then
        igndx = ignoredx
      else
        igndx = .false.
      endif
      if (present(upwind)) then
        upwnd = upwind
      else
        upwnd = .false.
        if (.not. lequidist(j)) then
          call fatal_error('der6_other','NOT IMPLEMENTED for no equidistant grid')
        endif
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=1.0
          else if (upwnd) then
            fac=(1.0/60)*dx_1(l1:l2)
          else
            fac=1/dx**6
          endif
          df=fac*(- 20.0* f(l1:l2,m,n) &
                  + 15.0*(f(l1+1:l2+1,m,n)+f(l1-1:l2-1,m,n)) &
                  -  6.0*(f(l1+2:l2+2,m,n)+f(l1-2:l2-2,m,n)) &
                  +      (f(l1+3:l2+3,m,n)+f(l1-3:l2-3,m,n)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=1.0
          else if (upwnd) then
            fac=(1.0/60)*dy_1(m)
          else
            fac=1/dy**6
          endif
          df=fac*(- 20.0* f(l1:l2,m  ,n) &
                  + 15.0*(f(l1:l2,m+1,n)+f(l1:l2,m-1,n)) &
                  -  6.0*(f(l1:l2,m+2,n)+f(l1:l2,m-2,n)) &
                  +      (f(l1:l2,m+3,n)+f(l1:l2,m-3,n)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=1.0
          else if (upwnd) then
            fac=(1.0/60)*dz_1(n)
          else
            fac=1/dz**6
          endif
          df=fac*(- 20.0* f(l1:l2,m,n  ) &
                  + 15.0*(f(l1:l2,m,n+1)+f(l1:l2,m,n-1)) &
                  -  6.0*(f(l1:l2,m,n+2)+f(l1:l2,m,n-2)) &
                  +      (f(l1:l2,m,n+3)+f(l1:l2,m,n-3)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der6_other
!***********************************************************************
    subroutine derij(f,k,df,i,j)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!   8-sep-01/axel: coded
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  14-nov-06/wolf: implemented bidiagonal scheme
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: i,j,k
!
!debug      if (loptimise_ders) der_call_count(k,icount_derij,i,j) = & !DERCOUNT
!debug                          der_call_count(k,icount_derij,i,j) + 1 !DERCOUNT
!

      if (lbidiagonal_derij) then
        !
        ! Use bidiagonal mixed-derivative operator, i.e.
        ! employ only the three neighbouring points on each of the four
        ! half-diagonals. This gives 6th-order mixed derivatives as the
        ! version below, but involves just 12 points instead of 36.
        !
        if ((i==1.and.j==2).or.(i==2.and.j==1)) then
          if (nxgrid/=1.and.nygrid/=1) then
            fac=(1./720.)*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
                        270.*( f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k)  &
                              +f(l1-1:l2-1,m-1,n,k)-f(l1+1:l2+1,m-1,n,k)) &
                       - 27.*( f(l1+2:l2+2,m+2,n,k)-f(l1-2:l2-2,m+2,n,k)  &
                              +f(l1-2:l2-2,m-2,n,k)-f(l1+2:l2+2,m-2,n,k)) &
                       +  2.*( f(l1+3:l2+3,m+3,n,k)-f(l1-3:l2-3,m+3,n,k)  &
                              +f(l1-3:l2-3,m-3,n,k)-f(l1+3:l2+3,m-3,n,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=(1./720.)*dy_1(m)*dz_1(n)
            df=fac*( &
                        270.*( f(l1:l2,m+1,n+1,k)-f(l1:l2,m+1,n-1,k)  &
                              +f(l1:l2,m-1,n-1,k)-f(l1:l2,m-1,n+1,k)) &
                       - 27.*( f(l1:l2,m+2,n+2,k)-f(l1:l2,m+2,n-2,k)  &
                              +f(l1:l2,m-2,n-2,k)-f(l1:l2,m-2,n+2,k)) &
                       +  2.*( f(l1:l2,m+3,n+3,k)-f(l1:l2,m+3,n-3,k)  &
                              +f(l1:l2,m-3,n-3,k)-f(l1:l2,m-3,n+3,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./720.)*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
                        270.*( f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k)  &
                              +f(l1-1:l2-1,m,n-1,k)-f(l1+1:l2+1,m,n-1,k)) &
                       - 27.*( f(l1+2:l2+2,m,n+2,k)-f(l1-2:l2-2,m,n+2,k)  &
                              +f(l1-2:l2-2,m,n-2,k)-f(l1+2:l2+2,m,n-2,k)) &
                       +  2.*( f(l1+3:l2+3,m,n+3,k)-f(l1-3:l2-3,m,n+3,k)  &
                              +f(l1-3:l2-3,m,n-3,k)-f(l1+3:l2+3,m,n-3,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or z-direction'
          endif
        endif

      else                      ! not using bidiagonal mixed derivatives
        !
        ! This is the old, straight-forward scheme
        !
        if ((i==1.and.j==2).or.(i==2.and.j==1)) then
          if (nxgrid/=1.and.nygrid/=1) then
            fac=(1./60.**2)*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
              45.*((45.*(f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k))  &
                    -9.*(f(l1+2:l2+2,m+1,n,k)-f(l1-2:l2-2,m+1,n,k))  &
                       +(f(l1+3:l2+3,m+1,n,k)-f(l1-3:l2-3,m+1,n,k))) &
                  -(45.*(f(l1+1:l2+1,m-1,n,k)-f(l1-1:l2-1,m-1,n,k))  &
                    -9.*(f(l1+2:l2+2,m-1,n,k)-f(l1-2:l2-2,m-1,n,k))  &
                       +(f(l1+3:l2+3,m-1,n,k)-f(l1-3:l2-3,m-1,n,k))))&
              -9.*((45.*(f(l1+1:l2+1,m+2,n,k)-f(l1-1:l2-1,m+2,n,k))  &
                    -9.*(f(l1+2:l2+2,m+2,n,k)-f(l1-2:l2-2,m+2,n,k))  &
                       +(f(l1+3:l2+3,m+2,n,k)-f(l1-3:l2-3,m+2,n,k))) &
                  -(45.*(f(l1+1:l2+1,m-2,n,k)-f(l1-1:l2-1,m-2,n,k))  &
                    -9.*(f(l1+2:l2+2,m-2,n,k)-f(l1-2:l2-2,m-2,n,k))  &
                       +(f(l1+3:l2+3,m-2,n,k)-f(l1-3:l2-3,m-2,n,k))))&
                 +((45.*(f(l1+1:l2+1,m+3,n,k)-f(l1-1:l2-1,m+3,n,k))  &
                    -9.*(f(l1+2:l2+2,m+3,n,k)-f(l1-2:l2-2,m+3,n,k))  &
                       +(f(l1+3:l2+3,m+3,n,k)-f(l1-3:l2-3,m+3,n,k))) &
                  -(45.*(f(l1+1:l2+1,m-3,n,k)-f(l1-1:l2-1,m-3,n,k))  &
                    -9.*(f(l1+2:l2+2,m-3,n,k)-f(l1-2:l2-2,m-3,n,k))  &
                       +(f(l1+3:l2+3,m-3,n,k)-f(l1-3:l2-3,m-3,n,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=(1./60.**2)*dy_1(m)*dz_1(n)
            df=fac*( &
              45.*((45.*(f(l1:l2,m+1,n+1,k)-f(l1:l2,m-1,n+1,k))  &
                    -9.*(f(l1:l2,m+2,n+1,k)-f(l1:l2,m-2,n+1,k))  &
                       +(f(l1:l2,m+3,n+1,k)-f(l1:l2,m-3,n+1,k))) &
                  -(45.*(f(l1:l2,m+1,n-1,k)-f(l1:l2,m-1,n-1,k))  &
                    -9.*(f(l1:l2,m+2,n-1,k)-f(l1:l2,m-2,n-1,k))  &
                       +(f(l1:l2,m+3,n-1,k)-f(l1:l2,m-3,n-1,k))))&
              -9.*((45.*(f(l1:l2,m+1,n+2,k)-f(l1:l2,m-1,n+2,k))  &
                    -9.*(f(l1:l2,m+2,n+2,k)-f(l1:l2,m-2,n+2,k))  &
                       +(f(l1:l2,m+3,n+2,k)-f(l1:l2,m-3,n+2,k))) &
                  -(45.*(f(l1:l2,m+1,n-2,k)-f(l1:l2,m-1,n-2,k))  &
                    -9.*(f(l1:l2,m+2,n-2,k)-f(l1:l2,m-2,n-2,k))  &
                       +(f(l1:l2,m+3,n-2,k)-f(l1:l2,m-3,n-2,k))))&
                 +((45.*(f(l1:l2,m+1,n+3,k)-f(l1:l2,m-1,n+3,k))  &
                    -9.*(f(l1:l2,m+2,n+3,k)-f(l1:l2,m-2,n+3,k))  &
                       +(f(l1:l2,m+3,n+3,k)-f(l1:l2,m-3,n+3,k))) &
                  -(45.*(f(l1:l2,m+1,n-3,k)-f(l1:l2,m-1,n-3,k))  &
                    -9.*(f(l1:l2,m+2,n-3,k)-f(l1:l2,m-2,n-3,k))  &
                       +(f(l1:l2,m+3,n-3,k)-f(l1:l2,m-3,n-3,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./60.**2)*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
              45.*((45.*(f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k))  &
                    -9.*(f(l1+2:l2+2,m,n+1,k)-f(l1-2:l2-2,m,n+1,k))  &
                       +(f(l1+3:l2+3,m,n+1,k)-f(l1-3:l2-3,m,n+1,k))) &
                  -(45.*(f(l1+1:l2+1,m,n-1,k)-f(l1-1:l2-1,m,n-1,k))  &
                    -9.*(f(l1+2:l2+2,m,n-1,k)-f(l1-2:l2-2,m,n-1,k))  &
                       +(f(l1+3:l2+3,m,n-1,k)-f(l1-3:l2-3,m,n-1,k))))&
              -9.*((45.*(f(l1+1:l2+1,m,n+2,k)-f(l1-1:l2-1,m,n+2,k))  &
                    -9.*(f(l1+2:l2+2,m,n+2,k)-f(l1-2:l2-2,m,n+2,k))  &
                       +(f(l1+3:l2+3,m,n+2,k)-f(l1-3:l2-3,m,n+2,k))) &
                  -(45.*(f(l1+1:l2+1,m,n-2,k)-f(l1-1:l2-1,m,n-2,k))  &
                    -9.*(f(l1+2:l2+2,m,n-2,k)-f(l1-2:l2-2,m,n-2,k))  &
                       +(f(l1+3:l2+3,m,n-2,k)-f(l1-3:l2-3,m,n-2,k))))&
                 +((45.*(f(l1+1:l2+1,m,n+3,k)-f(l1-1:l2-1,m,n+3,k))  &
                    -9.*(f(l1+2:l2+2,m,n+3,k)-f(l1-2:l2-2,m,n+3,k))  &
                       +(f(l1+3:l2+3,m,n+3,k)-f(l1-3:l2-3,m,n+3,k))) &
                  -(45.*(f(l1+1:l2+1,m,n-3,k)-f(l1-1:l2-1,m,n-3,k))  &
                    -9.*(f(l1+2:l2+2,m,n-3,k)-f(l1-2:l2-2,m,n-3,k))  &
                       +(f(l1+3:l2+3,m,n-3,k)-f(l1-3:l2-3,m,n-3,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or z-direction'
          endif
        endif

      endif                     ! bidiagonal derij
!
    endsubroutine derij
!***********************************************************************
    subroutine der5i1j(f,k,df,i,j)
!
!  Calculate 6th derivative with respect to two different directions.
!
!  05-dec-06/anders: adapted from derij
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: i,j,k
!
!debug      if (loptimise_ders) der_call_count(k,icount_derij,i,j) = & !DERCOUNT
!debug                          der_call_count(k,icount_derij,i,j) + 1 !DERCOUNT
!
      if ((i==1.and.j==1)) then
        if (nxgrid/=1) then
          call der6(f,k,df,j)
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in x-direction'
        endif
      elseif ((i==2.and.j==2)) then
        if (nygrid/=1) then
          call der6(f,k,df,j)
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in y-direction'
        endif
      elseif ((i==3.and.j==3)) then
        if (nzgrid/=1) then
          call der6(f,k,df,j)
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in z-direction'
        endif
      elseif ((i==1.and.j==2)) then
        if (nxgrid/=1.and.nygrid/=1) then
          fac=dx_1(l1:l2)**5*dy_1(m)
          df=fac*( &
            2.5*((45.*(f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k))  &
                  -9.*(f(l1+1:l2+1,m+2,n,k)-f(l1-1:l2-1,m+2,n,k))  &
                     +(f(l1+1:l2+1,m+3,n,k)-f(l1-1:l2-1,m+3,n,k))) &
                -(45.*(f(l1+1:l2+1,m-1,n,k)-f(l1-1:l2-1,m-1,n,k))  &
                  -9.*(f(l1+1:l2+1,m-2,n,k)-f(l1-1:l2-1,m-2,n,k))  &
                     +(f(l1+1:l2+1,m-3,n,k)-f(l1-1:l2-1,m-3,n,k))))&
           -2.0*((45.*(f(l1+2:l2+2,m+1,n,k)-f(l1-2:l2-2,m+1,n,k))  &
                  -9.*(f(l1+2:l2+2,m+2,n,k)-f(l1-2:l2-2,m+2,n,k))  &
                     +(f(l1+2:l2+2,m+3,n,k)-f(l1-2:l2-2,m+3,n,k))) &
                -(45.*(f(l1+2:l2+2,m-1,n,k)-f(l1-2:l2-2,m-1,n,k))  &
                  -9.*(f(l1+2:l2+2,m-2,n,k)-f(l1-2:l2-2,m-2,n,k))  &
                     +(f(l1+2:l2+2,m-3,n,k)-f(l1-2:l2-2,m-3,n,k))))&
           +0.5*((45.*(f(l1+3:l2+3,m+1,n,k)-f(l1-3:l2-3,m+1,n,k))  &
                  -9.*(f(l1+3:l2+3,m+2,n,k)-f(l1-3:l2-3,m+2,n,k))  &
                     +(f(l1+3:l2+3,m+3,n,k)-f(l1-3:l2-3,m+3,n,k))) &
                -(45.*(f(l1+3:l2+3,m-1,n,k)-f(l1-3:l2-3,m-1,n,k))  &
                  -9.*(f(l1+3:l2+3,m-2,n,k)-f(l1-3:l2-3,m-2,n,k))  &
                     +(f(l1+3:l2+3,m-3,n,k)-f(l1-3:l2-3,m-3,n,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in x- or y-direction'
        endif
      elseif ((i==2.and.j==1)) then
        if (nygrid/=1.and.nxgrid/=1) then
          fac=dy_1(m)**5*dx_1(l1:l2)
          df=fac*( &
            2.5*((45.*(f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k))  &
                  -9.*(f(l1+2:l2+2,m+1,n,k)-f(l1-2:l2-2,m+1,n,k))  &
                     +(f(l1+3:l2+3,m+1,n,k)-f(l1-3:l2-3,m+1,n,k))) &
                -(45.*(f(l1+1:l2+1,m-1,n,k)-f(l1-1:l2-1,m-1,n,k))  &
                  -9.*(f(l1+2:l2+2,m-1,n,k)-f(l1-2:l2-2,m-1,n,k))  &
                     +(f(l1+3:l2+3,m-1,n,k)-f(l1-3:l2-3,m-1,n,k))))&
           -2.0*((45.*(f(l1+1:l2+1,m+2,n,k)-f(l1-1:l2-1,m+2,n,k))  &
                  -9.*(f(l1+2:l2+2,m+2,n,k)-f(l1-2:l2-2,m+2,n,k))  &
                     +(f(l1+3:l2+3,m+2,n,k)-f(l1-3:l2-3,m+2,n,k))) &
                -(45.*(f(l1+1:l2+1,m-2,n,k)-f(l1-1:l2-1,m-2,n,k))  &
                  -9.*(f(l1+2:l2+2,m-2,n,k)-f(l1-2:l2-2,m-2,n,k))  &
                     +(f(l1+3:l2+3,m-2,n,k)-f(l1-3:l2-3,m-2,n,k))))&
           +0.5*((45.*(f(l1+1:l2+1,m+3,n,k)-f(l1-1:l2-1,m+3,n,k))  &
                  -9.*(f(l1+2:l2+2,m+3,n,k)-f(l1-2:l2-2,m+3,n,k))  &
                     +(f(l1+3:l2+3,m+3,n,k)-f(l1-3:l2-3,m+3,n,k))) &
                -(45.*(f(l1+1:l2+1,m-3,n,k)-f(l1-1:l2-1,m-3,n,k))  &
                  -9.*(f(l1+2:l2+2,m-3,n,k)-f(l1-2:l2-2,m-3,n,k))  &
                     +(f(l1+3:l2+3,m-3,n,k)-f(l1-3:l2-3,m-3,n,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in y- or x-direction'
        endif
      elseif ((i==1.and.j==3)) then
        if (nxgrid/=1.and.nzgrid/=1) then
          fac=dx_1(l1:l2)**5*dz_1(n)
          df=fac*( &
            2.5*((45.*(f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k))  &
                  -9.*(f(l1+1:l2+1,m,n+2,k)-f(l1-1:l2-1,m,n+2,k))  &
                     +(f(l1+1:l2+1,m,n+3,k)-f(l1-1:l2-1,m,n+3,k))) &
                -(45.*(f(l1+1:l2+1,m,n-1,k)-f(l1-1:l2-1,m,n-1,k))  &
                  -9.*(f(l1+1:l2+1,m,n-2,k)-f(l1-1:l2-1,m,n-2,k))  &
                     +(f(l1+1:l2+1,m,n-3,k)-f(l1-1:l2-1,m,n-3,k))))&
           -2.0*((45.*(f(l1+2:l2+2,m,n+1,k)-f(l1-2:l2-2,m,n+1,k))  &
                  -9.*(f(l1+2:l2+2,m,n+2,k)-f(l1-2:l2-2,m,n+2,k))  &
                     +(f(l1+2:l2+2,m,n+3,k)-f(l1-2:l2-2,m,n+3,k))) &
                -(45.*(f(l1+2:l2+2,m,n-1,k)-f(l1-2:l2-2,m,n-1,k))  &
                  -9.*(f(l1+2:l2+2,m,n-2,k)-f(l1-2:l2-2,m,n-2,k))  &
                     +(f(l1+2:l2+2,m,n-3,k)-f(l1-2:l2-2,m,n-3,k))))&
           +0.5*((45.*(f(l1+3:l2+3,m,n+1,k)-f(l1-3:l2-3,m,n+1,k))  &
                  -9.*(f(l1+3:l2+3,m,n+2,k)-f(l1-3:l2-3,m,n+2,k))  &
                     +(f(l1+3:l2+3,m,n+3,k)-f(l1-3:l2-3,m,n+3,k))) &
                -(45.*(f(l1+3:l2+3,m,n-1,k)-f(l1-3:l2-3,m,n-1,k))  &
                  -9.*(f(l1+3:l2+3,m,n-2,k)-f(l1-3:l2-3,m,n-2,k))  &
                     +(f(l1+3:l2+3,m,n-3,k)-f(l1-3:l2-3,m,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in x- or z-direction'
        endif
      elseif ((i==3.and.j==1)) then
        if (nzgrid/=1.and.nygrid/=1) then
          fac=dz_1(n)**5*dy_1(m)
          df=fac*( &
            2.5*((45.*(f(l1+1:l2+1,m,n+1,k)-f(l1+1-1:l2+1-1,m,n+1,k))  &
                  -9.*(f(l1+2:l2+2,m,n+1,k)-f(l1+2-1:l2+2-1,m,n+1,k))  &
                     +(f(l1+3:l2+3,m,n+1,k)-f(l1+3-1:l2+3-1,m,n+1,k))) &
                -(45.*(f(l1-1:l2-1,m,n+1,k)-f(l1-1-1:l2-1-1,m,n+1,k))  &
                  -9.*(f(l1-2:l2-2,m,n+1,k)-f(l1-2-1:l2-2-1,m,n+1,k))  &
                     +(f(l1-3:l2-3,m,n+1,k)-f(l1-3-1:l2-3-1,m,n+1,k))))&
           -2.0*((45.*(f(l1+1:l2+1,m,n+2,k)-f(l1+1-2:l2+1-2,m,n+2,k))  &
                  -9.*(f(l1+2:l2+2,m,n+2,k)-f(l1+2-2:l2+2-2,m,n+2,k))  &
                     +(f(l1+3:l2+3,m,n+2,k)-f(l1+3-2:l2+3-2,m,n+2,k))) &
                -(45.*(f(l1-1:l2-1,m,n+2,k)-f(l1-1-2:l2-1-2,m,n+2,k))  &
                  -9.*(f(l1-2:l2-2,m,n+2,k)-f(l1-2-2:l2-2-2,m,n+2,k))  &
                     +(f(l1-3:l2-3,m,n+2,k)-f(l1-3-2:l2-3-2,m,n+2,k))))&
           +0.5*((45.*(f(l1+1:l2+1,m,n+3,k)-f(l1+1-3:l2+1-3,m,n+3,k))  &
                  -9.*(f(l1+2:l2+2,m,n+3,k)-f(l1+2-3:l2+2-3,m,n+3,k))  &
                     +(f(l1+3:l2+3,m,n+3,k)-f(l1+3-3:l2+3-3,m,n+3,k))) &
                -(45.*(f(l1-1:l2-1,m,n+3,k)-f(l1-1-3:l2-1-3,m,n+3,k))  &
                  -9.*(f(l1-2:l2-2,m,n+3,k)-f(l1-2-3:l2-2-3,m,n+3,k))  &
                     +(f(l1-3:l2-3,m,n+3,k)-f(l1-3-3:l2-3-3,m,n+3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in z- or x-direction'
        endif
      elseif ((i==2.and.j==3)) then
        if (nygrid/=1.and.nzgrid/=1) then
          fac=dy_1(m)**5*dz_1(n)
          df=fac*( &
            2.5*((45.*(f(l1:l2,m+1,n+1,k)-f(l1:l2,m+1,n+1,k))  &
                  -9.*(f(l1:l2,m+1,n+2,k)-f(l1:l2,m+1,n+2,k))  &
                     +(f(l1:l2,m+1,n+3,k)-f(l1:l2,m+1,n+3,k))) &
                -(45.*(f(l1:l2,m+1,n-1,k)-f(l1:l2,m+1,n-1,k))  &
                  -9.*(f(l1:l2,m+1,n-2,k)-f(l1:l2,m+1,n-2,k))  &
                     +(f(l1:l2,m+1,n-3,k)-f(l1:l2,m+1,n-3,k))))&
           -2.0*((45.*(f(l1:l2,m+2,n+1,k)-f(l1:l2,m+2,n+1,k))  &
                  -9.*(f(l1:l2,m+2,n+2,k)-f(l1:l2,m+2,n+2,k))  &
                     +(f(l1:l2,m+2,n+3,k)-f(l1:l2,m+2,n+3,k))) &
                -(45.*(f(l1:l2,m+2,n-1,k)-f(l1:l2,m+2,n-1,k))  &
                  -9.*(f(l1:l2,m+2,n-2,k)-f(l1:l2,m+2,n-2,k))  &
                     +(f(l1:l2,m+2,n-3,k)-f(l1:l2,m+2,n-3,k))))&
           +0.5*((45.*(f(l1:l2,m+3,n+1,k)-f(l1:l2,m+3,n+1,k))  &
                  -9.*(f(l1:l2,m+3,n+2,k)-f(l1:l2,m+3,n+2,k))  &
                     +(f(l1:l2,m+3,n+3,k)-f(l1:l2,m+3,n+3,k))) &
                -(45.*(f(l1:l2,m+3,n-1,k)-f(l1:l2,m+3,n-1,k))  &
                  -9.*(f(l1:l2,m+3,n-2,k)-f(l1:l2,m+3,n-2,k))  &
                     +(f(l1:l2,m+3,n-3,k)-f(l1:l2,m+3,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in y- or z-direction'
        endif
      elseif ((i==3.and.j==2)) then
        if (nzgrid/=1.and.nygrid/=1) then
          fac=dz_1(n)**5*dy_1(m)
          df=fac*( &
            2.5*((45.*(f(l1:l2,m+1,n+1,k)-f(l1:l2,m+1,n+1,k))  &
                  -9.*(f(l1:l2,m+2,n+1,k)-f(l1:l2,m+2,n+1,k))  &
                     +(f(l1:l2,m+3,n+1,k)-f(l1:l2,m+3,n+1,k))) &
                -(45.*(f(l1:l2,m-1,n+1,k)-f(l1:l2,m-1,n+1,k))  &
                  -9.*(f(l1:l2,m-2,n+1,k)-f(l1:l2,m-2,n+1,k))  &
                     +(f(l1:l2,m-3,n+1,k)-f(l1:l2,m-3,n+1,k))))&
           -2.0*((45.*(f(l1:l2,m+1,n+2,k)-f(l1:l2,m+1,n+2,k))  &
                  -9.*(f(l1:l2,m+2,n+2,k)-f(l1:l2,m+2,n+2,k))  &
                     +(f(l1:l2,m+3,n+2,k)-f(l1:l2,m+3,n+2,k))) &
                -(45.*(f(l1:l2,m-1,n+2,k)-f(l1:l2,m-1,n+2,k))  &
                  -9.*(f(l1:l2,m-2,n+2,k)-f(l1:l2,m-2,n+2,k))  &
                     +(f(l1:l2,m-3,n+2,k)-f(l1:l2,m-3,n+2,k))))&
           +0.5*((45.*(f(l1:l2,m+1,n+3,k)-f(l1:l2,m+1,n+3,k))  &
                  -9.*(f(l1:l2,m+2,n+3,k)-f(l1:l2,m+2,n+3,k))  &
                     +(f(l1:l2,m+3,n+3,k)-f(l1:l2,m+3,n+3,k))) &
                -(45.*(f(l1:l2,m-1,n+3,k)-f(l1:l2,m-1,n+3,k))  &
                  -9.*(f(l1:l2,m-2,n+3,k)-f(l1:l2,m-2,n+3,k))  &
                     +(f(l1:l2,m-3,n+3,k)-f(l1:l2,m-3,n+3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in z- or y-direction'
        endif
      else
        print*, 'der5i1j: no such value for i,j=', i, j
        call fatal_error('der5i1j','')
      endif
!
    endsubroutine der5i1j
!***********************************************************************
    subroutine der_upwind1st(f,uu,k,df,j)
!
!  First order upwind derivative of variable
!
!  Useful for advecting non-logarithmic variables
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: df
      integer :: j,k,l
!
!debug      if (loptimise_ders) der_call_count(k,icount_der_upwind1st,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der_upwind1st,j,1) + 1 !DERCOUNT
!
      if (.not. lequidist(j)) then
        call fatal_error('der_upwind1st','NOT IMPLEMENTED for no equidistant grid')
      endif
!
      if (j == 1) then
        if (nxgrid /= 1) then
          do l=1,nx
            if (uu(3+l,1) > 0.) then
              df(l) = (f(3+l,m,n,k) - f(3+l-1,m,n,k))/dx
            else
              df(l) = (f(3+l+1,m,n,k) - f(3+l,m,n,k))/dx
            endif
          enddo
        else
          df=0.
          if (ip<=5) print*, 'der_upwind1st: Degenerate case in x-direction'
        endif
      elseif (j == 2) then
        if (nygrid /= 1) then
          do l=1,nx
            if (uu(l,2) > 0.) then
              df(l) = (f(3+l,m,n,k) - f(3+l,m-1,n,k))/dy
            else
              df(l) = (f(3+l,m+1,n,k) - f(3+l,m,n,k))/dy
            endif
          enddo
        else
          df=0.
          if (ip<=5) print*, 'der_upwind1st: Degenerate case in y-direction'
        endif
      elseif (j == 3) then
        if (nzgrid /= 1) then
          do l=1,nx
            if (uu(l,3) > 0.) then
              df(l) = (f(3+l,m,n,k) - f(3+l,m,n-1,k))/dz
            else
              df(l) = (f(3+l,m,n+1,k) - f(3+l,m,n,k))/dz
            endif
          enddo
        else
          df=0.
          if (ip<=5) print*, 'der_upwind1st: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_upwind1st
!***********************************************************************


endmodule Deriv
