! $Id: deriv.f90,v 1.2 2002-06-01 02:56:21 brandenb Exp $

module Deriv

  implicit none

  contains

!***********************************************************************
    subroutine der(f,k,df,j)
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
        if (mx/=1) then
          fac=1./(60.*dx)
          df=fac*(45.*(f(5:mx-2,m,n,k)-f(3:mx-4,m,n,k)) &
                  -9.*(f(6:mx-1,m,n,k)-f(2:mx-5,m,n,k)) &
                     +(f(7:mx  ,m,n,k)-f(1:mx-6,m,n,k)))
        else
          df=0.
          if (ip.le.10) print*, 'Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (my/=1) then
          fac=1./(60.*dy)
          df=fac*(45.*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                  -9.*(f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) &
                     +(f(l1:l2,m+3,n,k)-f(l1:l2,m-3,n,k)))
        else
          df=0.
          if (ip.le.10) print*, 'Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (mz/=1) then
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
    endsubroutine der
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
        if (mx/=1) then
          fac=1./(180.*dx**2)
          df=fac*(-490.*f(4:mx-3,m,n,k) &
                 +270.*(f(5:mx-2,m,n,k)+f(3:mx-4,m,n,k)) &
                  -27.*(f(6:mx-1,m,n,k)+f(2:mx-5,m,n,k)) &
                   +2.*(f(7:mx  ,m,n,k)+f(1:mx-6,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (my/=1) then
          fac=1./(180.*dy**2)
          df=fac*(-490.*f(l1:l2,m  ,n,k) &
                 +270.*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  -27.*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                   +2.*(f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (mz/=1) then
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
        if (mx/=1.and.my/=1) then
          fac=1./(60.**2*dx*dy)
          df=fac*( &
            45.*((45.*(f(5:mx-2,m+1,n,k)-f(3:mx-4,m+1,n,k))  &
                  -9.*(f(6:mx-1,m+1,n,k)-f(2:mx-5,m+1,n,k))  &
                     +(f(7:mx  ,m+1,n,k)-f(1:mx-6,m+1,n,k))) &
                -(45.*(f(5:mx-2,m-1,n,k)-f(3:mx-4,m-1,n,k))  &
                  -9.*(f(6:mx-1,m-1,n,k)-f(2:mx-5,m-1,n,k))  &
                     +(f(7:mx  ,m-1,n,k)-f(1:mx-6,m-1,n,k))))&
            -9.*((45.*(f(5:mx-2,m+2,n,k)-f(3:mx-4,m+2,n,k))  &
                  -9.*(f(6:mx-1,m+2,n,k)-f(2:mx-5,m+2,n,k))  &
                     +(f(7:mx  ,m+2,n,k)-f(1:mx-6,m+2,n,k))) &
                -(45.*(f(5:mx-2,m-2,n,k)-f(3:mx-4,m-2,n,k))  &
                  -9.*(f(6:mx-1,m-2,n,k)-f(2:mx-5,m-2,n,k))  &
                     +(f(7:mx  ,m-2,n,k)-f(1:mx-6,m-2,n,k))))&
               +((45.*(f(5:mx-2,m+3,n,k)-f(3:mx-4,m+3,n,k))  &
                  -9.*(f(6:mx-1,m+3,n,k)-f(2:mx-5,m+3,n,k))  &
                     +(f(7:mx  ,m+3,n,k)-f(1:mx-6,m+3,n,k))) &
                -(45.*(f(5:mx-2,m-3,n,k)-f(3:mx-4,m-3,n,k))  &
                  -9.*(f(6:mx-1,m-3,n,k)-f(2:mx-5,m-3,n,k))  &
                     +(f(7:mx  ,m-3,n,k)-f(1:mx-6,m-3,n,k))))&
                 )
        else
          df=0.
          if (ip.le.10) print*, 'Degenerate case in x-direction'
        endif
      elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
        if (my/=1.and.mz/=1) then
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
        if (mz/=1.and.mx/=1) then
          fac=1./(60.**2*dz*dx)
          df=fac*( &
            45.*((45.*(f(5:mx-2,m,n+1,k)-f(3:mx-4,m,n+1,k))  &
                  -9.*(f(6:mx-1,m,n+1,k)-f(2:mx-5,m,n+1,k))  &
                     +(f(7:mx  ,m,n+1,k)-f(1:mx-6,m,n+1,k))) &
                -(45.*(f(5:mx-2,m,n-1,k)-f(3:mx-4,m,n-1,k))  &
                  -9.*(f(6:mx-1,m,n-1,k)-f(2:mx-5,m,n-1,k))  &
                     +(f(7:mx  ,m,n-1,k)-f(1:mx-6,m,n-1,k))))&
            -9.*((45.*(f(5:mx-2,m,n+2,k)-f(3:mx-4,m,n+2,k))  &
                  -9.*(f(6:mx-1,m,n+2,k)-f(2:mx-5,m,n+2,k))  &
                     +(f(7:mx  ,m,n+2,k)-f(1:mx-6,m,n+2,k))) &
                -(45.*(f(5:mx-2,m,n-2,k)-f(3:mx-4,m,n-2,k))  &
                  -9.*(f(6:mx-1,m,n-2,k)-f(2:mx-5,m,n-2,k))  &
                     +(f(7:mx  ,m,n-2,k)-f(1:mx-6,m,n-2,k))))&
               +((45.*(f(5:mx-2,m,n+3,k)-f(3:mx-4,m,n+3,k))  &
                  -9.*(f(6:mx-1,m,n+3,k)-f(2:mx-5,m,n+3,k))  &
                     +(f(7:mx  ,m,n+3,k)-f(1:mx-6,m,n+3,k))) &
                -(45.*(f(5:mx-2,m,n-3,k)-f(3:mx-4,m,n-3,k))  &
                  -9.*(f(6:mx-1,m,n-3,k)-f(2:mx-5,m,n-3,k))  &
                     +(f(7:mx  ,m,n-3,k)-f(1:mx-6,m,n-3,k))))&
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
