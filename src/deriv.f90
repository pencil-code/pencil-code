! $Id: deriv.f90,v 1.18 2004-07-03 02:13:13 theine Exp $

module Deriv
  
  use Mpicomm, only: stop_it

  implicit none

  interface der                 ! Overload the der function
    module procedure der_main   ! derivative of an 'mvar' variable
    module procedure der_other  ! derivative of another field
  endinterface

  
  contains

!***********************************************************************
    subroutine initialize_deriv()
!
!  should be moved to grid.f90
!
      use Cdata, only: dx,dy,dz
      use Cdata, only: dx_1,dy_1,dz_1

      dx_1=1./dx
      dy_1=1./dy
      dz_1=1./dz

    endsubroutine initialize_deriv
!***********************************************************************
    subroutine der_main(f,k,df,j)
!
!  calculate derivative df_k/dx_j
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!  28-may-04/bing+wolf: adapted for nonequidistant grid
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: df
      integer :: j,k

      call der_equidist(f,k,df,j)
      if (.not. lequidist(j)) then
         if (j==1) then
            df = df * dx * xiprim(l1:l2)
         elseif (j==2) then
            df = df * dy * psiprim(m)
         elseif (j==3) then
            df = df * dz * zetaprim(n)
         endif
      endif
!
    end subroutine der_main
!***********************************************************************
    subroutine der_equidist(f,k,df,j)
!
!  calculate derivative df_k/dx_j
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!   1-oct-97/axel: coded
!  18-jul-98/axel: corrected mx -> my and mx -> mz in all y and z ders
!   1-apr-01/axel+wolf: pencil formulation
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
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
          if (ip.le.10) print*, 'der_main: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=1./(60.*dy)
          df=fac*(45.*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                  -9.*(f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) &
                     +(f(l1:l2,m+3,n,k)-f(l1:l2,m-3,n,k)))
        else
          df=0.
          if (ip.le.10) print*, 'der_main: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=1./(60.*dz)
          df=fac*(45.*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k)) &
                  -9.*(f(l1:l2,m,n+2,k)-f(l1:l2,m,n-2,k)) &
                     +(f(l1:l2,m,n+3,k)-f(l1:l2,m,n-3,k)))
        else
          df=0.
          if (ip.le.10) print*, 'der_main: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_equidist
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
      if (.not. lequidist(j)) then
        call stop_it('der_other: NOT IMPLEMENTED for no equidistant grid')
      endif

      if (j==1) then
        if (nxgrid/=1) then
          fac=1./(60.*dx)
          df=fac*(45.*(f(l1+1:l2+1,m,n)-f(l1-1:l2-1,m,n)) &
                  -9.*(f(l1+2:l2+2,m,n)-f(l1-2:l2-2,m,n)) &
                     +(f(l1+3:l2+3,m,n)-f(l1-3:l2-3,m,n)))
        else
          df=0.
          if (ip.le.10) print*, 'der_other: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=1./(60.*dy)
          df=fac*(45.*(f(l1:l2,m+1,n)-f(l1:l2,m-1,n)) &
                  -9.*(f(l1:l2,m+2,n)-f(l1:l2,m-2,n)) &
                     +(f(l1:l2,m+3,n)-f(l1:l2,m-3,n)))
        else
          df=0.
          if (ip.le.10) print*, 'der_other: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=1./(60.*dz)
          df=fac*(45.*(f(l1:l2,m,n+1)-f(l1:l2,m,n-1)) &
                  -9.*(f(l1:l2,m,n+2)-f(l1:l2,m,n-2)) &
                     +(f(l1:l2,m,n+3)-f(l1:l2,m,n-3)))
        else
          df=0.
          if (ip.le.10) print*, 'der_other: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_other
!***********************************************************************
    subroutine der2(f,k,df2,j)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!
!  28-may-04/bing+wolf: adapted for nonequidistant grid
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: df2,df
      integer :: j,k

      intent(in)  :: f,k,j 
      intent(out) :: df2
      
      call der2_equidist(f,k,df2,j)
      if (.not. lequidist(j)) then
         call der_equidist(f,k,df,j)
         if (j==1) then
            df2 = df2*(dx*xiprim(l1:l2))**2+ df*xiprim2(l1:l2)*dx
         elseif (j==2) then
            df2 = df2*(dy*psiprim(m))**2   + df*psiprim2(m)*dy
         elseif (j==3) then
            df2 = df2*(dy*zetaprim(n))**2  + df*zetaprim2(n)*dz
         endif
      endif

    end subroutine der2
!***********************************************************************
    subroutine der2_equidist(f,k,df,j)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
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
    endsubroutine der2_equidist
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: df
      real :: fac
      integer :: j,k
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
      if (.not. lequidist(j)) then
        call stop_it('der6: NOT IMPLEMENTED for no equidistant grid')
      endif
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
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=1.
          else if (upwnd) then
            fac = 1./(60*dx)
          else
            fac=1./dx**6
          endif
          df=fac*(-20.* f(l1:l2,m,n,k) &
                  +15.*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  - 6.*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                  +    (f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=1.
          else if (upwnd) then
            fac = 1./(60*dy)
          else
            fac=1./dy**6
          endif
          df=fac*(-20.* f(l1:l2,m  ,n,k) &
                  +15.*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  - 6.*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                  +    (f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=1.
          else if (upwnd) then
            fac = 1./(60*dz)
          else
            fac=1./dz**6
          endif
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
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: df
      real :: fac
      integer :: j,k
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
      if (.not. lequidist(j)) then
        call stop_it('der4: NOT IMPLEMENTED for no equidistant grid')
      endif

      if (present(ignoredx)) then
        igndx = ignoredx
      else
        igndx = .false.
      endif
      if (present(upwind)) then
        upwnd = upwind
        print*, 'WARNING: upwinding not implemented for der4'
      else
        upwnd = .false.
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=1.
          else
            fac=1./dx**4
          endif
          df=fac*(-56.* f(l1:l2,m,n,k) &
                  +39.*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  -12.*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                  +    (f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=1.
          else
            fac=1./dy**4
          endif
          df=fac*(-56.* f(l1:l2,m  ,n,k) &
                  +39.*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  -12.*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                  +    (f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=1.
          else
            fac=1./dz*4
          endif
          df=fac*(-56.* f(l1:l2,m,n  ,k) &
                  +39.*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                  -12.*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                  +    (f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der4
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
          if (.not. lequidist(1)) then
            df =  df * dz * xiprim(l1:l2)
          endif
          if (.not. lequidist(2)) then
            df =  df * dz * psiprim(m)
          endif
        else
          df=0.
          if (ip.le.10) print*, 'derij: Degenerate case in x-direction'
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
          if (.not. lequidist(3)) then
            df =  df * dz * zetaprim(n)
          endif
          if (.not. lequidist(2)) then
            df =  df * dz * psiprim(m)
          endif
        else
          df=0.
          if (ip.le.10) print*, 'derij: Degenerate case in y-direction'
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
          if (.not. lequidist(3)) then
            df =  df * dz * zetaprim(n)
          endif 
          if (.not. lequidist(1)) then
            df =  df * dz * xiprim(l1:l2)
          endif
        else
          df=0.
          if (ip.le.10) print*, 'derij: Degenerate case in z-direction'
        endif
      !else
      ! (don't waste any time if i=j)
      endif
!
    endsubroutine derij
!***********************************************************************
    subroutine der_upwind1st(f,uu,k,df,j)
!
!  First order upwind derivative of variable
!
!  Useful for advecting non-logarithmic variables
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: df
      integer :: j,k,l
!
      if (.not. lequidist(j)) then
        call stop_it('der_upwind1st: NOT IMPLEMENTED for no equidistant grid')
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
          if (ip.le.10) print*, 'der_upwind1st: Degenerate case in x-direction'
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
          if (ip.le.10) print*, 'der_upwind1st: Degenerate case in y-direction'
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
          if (ip.le.10) print*, 'der_upwind1st: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_upwind1st
!***********************************************************************


endmodule Deriv
