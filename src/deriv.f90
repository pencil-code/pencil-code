! $Id: deriv.f90,v 1.6 2002-11-26 19:59:19 mee Exp $

module Deriv

  implicit none

  interface der                 ! Overload the der function
    module procedure der_main   ! derivative of an 'mvar' variable
    module procedure der_other  ! derivative of another field
  endinterface

  contains

!***********************************************************************
    subroutine der_main(f,k,df,j)
!
!  calculate derivative of a scalar, get scalar
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!   1-oct-97/axel: coded
!  18-jul-98/axel: corrected mx -> my and mx -> mz in all y and z ders
!   1-apr-01/axel+wolf: pencil formulation
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: df
      real :: fac
      integer :: j,k
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=1./(60.*dx)
          df=fac*(45.*(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k)) &
                  -9.*(f(l1+2:l2+2,m,n,k)-f(l1-2:l2-2,m,n,k)) &
                     +(f(l1+3:l2+3,m,n,k)-f(l1-3:l2-3,m,n,k)))
        else
          df=0.
          if (ip.le.10) print*, 'Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=1./(60.*dy)
          df=fac*(45.*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                  -9.*(f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) &
                     +(f(l1:l2,m+3,n,k)-f(l1:l2,m-3,n,k)))
        else
          df=0.
          if (ip.le.10) print*, 'Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=1./(60.*dz)
          df=fac*(45.*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k)) &
                  -9.*(f(l1:l2,m,n+2,k)-f(l1:l2,m,n-2,k)) &
                     +(f(l1:l2,m,n+3,k)-f(l1:l2,m,n-3,k)))
        else
          df=0.
          if (ip.le.10) print*, 'Degenerate case in z-direction'
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

      use Cdata
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df
      real :: fac
      integer :: j
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=1./(60.*dx)
          df=fac*(45.*(f(l1+1:l2+1,m,n)-f(l1-1:l2-1,m,n)) &
                  -9.*(f(l1+2:l2+2,m,n)-f(l1-2:l2-2,m,n)) &
                     +(f(l1+3:l2+3,m,n)-f(l1-3:l2-3,m,n)))
        else
          df=0.
          if (ip.le.10) print*, 'Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=1./(60.*dy)
          df=fac*(45.*(f(l1:l2,m+1,n)-f(l1:l2,m-1,n)) &
                  -9.*(f(l1:l2,m+2,n)-f(l1:l2,m-2,n)) &
                     +(f(l1:l2,m+3,n)-f(l1:l2,m-3,n)))
        else
          df=0.
          if (ip.le.10) print*, 'Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=1./(60.*dz)
          df=fac*(45.*(f(l1:l2,m,n+1)-f(l1:l2,m,n-1)) &
                  -9.*(f(l1:l2,m,n+2)-f(l1:l2,m,n-2)) &
                     +(f(l1:l2,m,n+3)-f(l1:l2,m,n-3)))
        else
          df=0.
          if (ip.le.10) print*, 'Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_other
!***********************************************************************
    subroutine der2(f,k,df,j)
!
!  calculate 2nd derivative of a scalar, get scalar
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: df
      real :: fac
      integer :: j,k
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=1./(180.*dx**2)
          df=fac*(-490.*f(l1:l2    ,m,n,k) &
                 +270.*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  -27.*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                   +2.*(f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=1./(180.*dy**2)
          df=fac*(-490.*f(l1:l2,m  ,n,k) &
                 +270.*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  -27.*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                   +2.*(f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=1./(180.*dz**2)
          df=fac*(-490.*f(l1:l2,m,n  ,k) &
                 +270.*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                  -27.*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                   +2.*(f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der2
!***********************************************************************
    subroutine der6(f,k,df,j,ignoredx)
!
!  Calculate 6th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!
!   8-jul-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: df
      real :: fac
      integer :: j,k
      logical, optional :: ignoredx
      logical :: igndx
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
      if (present(ignoredx)) then
        igndx = ignoredx
      else
        igndx = .false.
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then; fac=1.; else; fac=1./dx**6; endif
          df=fac*(-20.* f(l1:l2,m,n,k) &
                  +15.*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  - 6.*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                  +    (f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then; fac=1.; else; fac=1./dy**6; endif
          df=fac*(-20.* f(l1:l2,m  ,n,k) &
                  +15.*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  - 6.*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                  +    (f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then; fac=1.; else; fac=1./dz**6; endif
          df=fac*(-20.* f(l1:l2,m,n  ,k) &
                  +15.*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                  - 6.*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                  +    (f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der6
!***********************************************************************
    subroutine derij(f,k,df,i,j)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!   8-sep-01/axel: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (nx) :: df
      real :: fac
      integer :: i,j,k
!
      if ((i==1.and.j==2).or.(i==2.and.j==1)) then
        if (nxgrid/=1.and.nygrid/=1) then
          fac=1./(60.**2*dx*dy)
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
          if (ip.le.10) print*, 'Degenerate case in x-direction'
        endif
      elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
        if (nygrid/=1.and.nzgrid/=1) then
          fac=1./(60.**2*dy*dz)
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
          if (ip.le.10) print*, 'Degenerate case in y-direction'
        endif
      elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
        if (nzgrid/=1.and.nxgrid/=1) then
          fac=1./(60.**2*dz*dx)
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
          if (ip.le.10) print*, 'Degenerate case in z-direction'
        endif
      !else
      ! (don't waste any time if i=j)
      endif
!
    endsubroutine derij
!***********************************************************************


endmodule Deriv
