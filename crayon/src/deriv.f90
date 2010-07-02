! $Id: deriv.f90 13582 2010-03-31 15:01:18Z sven.bingert $
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM integer, parameter :: nghost = 3
!
!***************************************************************
module Deriv
!
  use Messages, only: fatal_error, warning
  use Cdata
!
  implicit none
!
  private
!
  public :: initialize_deriv
  public :: der, der2, derij
!
  real :: der2_coef0, der2_coef1, der2_coef2, der2_coef3
!
  interface der                   ! Overload the der function
    module procedure der_main     ! derivative of an 'mvar' variable
    module procedure der_other    ! derivative of another field
  endinterface
!
  interface der2                  ! Overload the der function
    module procedure der2_main    ! derivative of an 'mvar' variable
    module procedure der2_other   ! derivative of another field
  endinterface
!
  interface derij                 ! Overload the der function
    module procedure derij_main   ! derivative of an 'mvar' variable
    module procedure derij_other  ! derivative of another field
  endinterface
!
  contains
!
!***********************************************************************
    subroutine initialize_deriv()
!
!  Initialize stencil coefficients
!
      select case (der2_type)
!
      case ('standard')
        der2_coef0=-490./180.; der2_coef1=270./180.
        der2_coef2=-27./180.; der2_coef3=2./180.
!
      case ('tuned1')
        der2_coef0=-0.75; der2_coef1=0.34375
        der2_coef2=0.125; der2_coef3=-0.09375
!
      case default
        write(unit=errormsg,fmt=*) &
             "der2_type doesn't exist"
        call fatal_error('initialize_deriv',errormsg)
!
      endselect
!
    endsubroutine initialize_deriv
!***********************************************************************
    subroutine der_main(f,k,df,j)
!
!  calculate derivative df_k/dx_j
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!   1-oct-97/axel: coded
!  18-jul-98/axel: corrected mx -> my and mx -> mz in all y and z ders
!   1-apr-01/axel+wolf: pencil formulation
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: j,k
!
      intent(in)  :: f,k,j
      intent(out) :: df
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
!
!  26-nov-02/tony: coded - duplicate der_main but without k subscript
!                          then overload the der interface.
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df,fac
      integer :: j
!
      intent(in)  :: f,j
      intent(out) :: df
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
    subroutine der2_main(f,k,df2,j)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df2,fac
      integer :: j,k
!
      intent(in)  :: f,k,j
      intent(out) :: df2
!
!debug      if (loptimise_ders) der_call_count(k,icount_der2,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der2,j,1) + 1 !DERCOUNT
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=dx_1(l1:l2)**2
          df2=fac*(der2_coef0*f (l1  :l2  ,m,n,k) &
                  +der2_coef1*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  +der2_coef2*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                  +der2_coef3*(f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
        else
          df2=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=dy_1(m)**2
          df2=fac*(der2_coef0*f(l1:l2,m  ,n,k) &
                  +der2_coef1*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  +der2_coef2*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                  +der2_coef3*(f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
        else
          df2=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=dz_1(n)**2
          df2=fac*(der2_coef0*f(l1:l2,m,n  ,k) &
                   +der2_coef1*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                   +der2_coef2*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                   +der2_coef3*(f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
        else
          df2=0.
        endif
      endif
!
    endsubroutine der2_main
!***********************************************************************
    subroutine der2_other(f,df2,j)
!
!  calculate 2nd derivative d^2f/dx_j^2 (of scalar f)
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df2,fac
      integer :: j
!
      intent(in)  :: f,j
      intent(out) :: df2
!
!debug      if (loptimise_ders) der_call_count(k,icount_der2,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der2,j,1) + 1 !DERCOUNT
!
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=(1./180)*dx_1(l1:l2)**2
          df2=fac*(-490.0*f(l1:l2,m,n) &
                   +270.0*(f(l1+1:l2+1,m,n)+f(l1-1:l2-1,m,n)) &
                   - 27.0*(f(l1+2:l2+2,m,n)+f(l1-2:l2-2,m,n)) &
                   +  2.0*(f(l1+3:l2+3,m,n)+f(l1-3:l2-3,m,n)))
        else
          df2=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=(1./180)*dy_1(m)**2
          df2=fac*(-490.0*f(l1:l2,m,n) &
                   +270.0*(f(l1:l2,m+1,n)+f(l1:l2,m-1,n)) &
                   - 27.0*(f(l1:l2,m+2,n)+f(l1:l2,m-2,n)) &
                   +  2.0*(f(l1:l2,m+3,n)+f(l1:l2,m-3,n)))
        else
          df2=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=(1./180)*dz_1(n)**2
          df2=fac*(-490.0*f(l1:l2,m,n) &
                   +270.0*(f(l1:l2,m,n+1)+f(l1:l2,m,n-1)) &
                   - 27.0*(f(l1:l2,m,n+2)+f(l1:l2,m,n-2)) &
                   +  2.0*(f(l1:l2,m,n+3)+f(l1:l2,m,n-3)))
        else
          df2=0.
        endif
      endif
!
    endsubroutine der2_other
!***********************************************************************
    subroutine derij_main(f,k,df,i,j)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!
!   8-sep-01/axel: coded
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  14-nov-06/wolf: implemented bidiagonal scheme
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
!
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
!
      endif                     ! bidiagonal derij
!
    endsubroutine derij_main
!***********************************************************************
    subroutine derij_other(f,df,i,j)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!   8-sep-01/axel: coded
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  14-nov-06/wolf: implemented bidiagonal scheme
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df,fac
      integer :: i,j
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
                        270.*( f(l1+1:l2+1,m+1,n)-f(l1-1:l2-1,m+1,n)  &
                              +f(l1-1:l2-1,m-1,n)-f(l1+1:l2+1,m-1,n)) &
                       - 27.*( f(l1+2:l2+2,m+2,n)-f(l1-2:l2-2,m+2,n)  &
                              +f(l1-2:l2-2,m-2,n)-f(l1+2:l2+2,m-2,n)) &
                       +  2.*( f(l1+3:l2+3,m+3,n)-f(l1-3:l2-3,m+3,n)  &
                              +f(l1-3:l2-3,m-3,n)-f(l1+3:l2+3,m-3,n)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=(1./720.)*dy_1(m)*dz_1(n)
            df=fac*( &
                        270.*( f(l1:l2,m+1,n+1)-f(l1:l2,m+1,n-1)  &
                              +f(l1:l2,m-1,n-1)-f(l1:l2,m-1,n+1)) &
                       - 27.*( f(l1:l2,m+2,n+2)-f(l1:l2,m+2,n-2)  &
                              +f(l1:l2,m-2,n-2)-f(l1:l2,m-2,n+2)) &
                       +  2.*( f(l1:l2,m+3,n+3)-f(l1:l2,m+3,n-3)  &
                              +f(l1:l2,m-3,n-3)-f(l1:l2,m-3,n+3)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./720.)*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
                        270.*( f(l1+1:l2+1,m,n+1)-f(l1-1:l2-1,m,n+1)  &
                              +f(l1-1:l2-1,m,n-1)-f(l1+1:l2+1,m,n-1)) &
                       - 27.*( f(l1+2:l2+2,m,n+2)-f(l1-2:l2-2,m,n+2)  &
                              +f(l1-2:l2-2,m,n-2)-f(l1+2:l2+2,m,n-2)) &
                       +  2.*( f(l1+3:l2+3,m,n+3)-f(l1-3:l2-3,m,n+3)  &
                              +f(l1-3:l2-3,m,n-3)-f(l1+3:l2+3,m,n-3)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or z-direction'
          endif
        endif
!
      else                      ! not using bidiagonal mixed derivatives
        !
        ! This is the old, straight-forward scheme
        !
        if ((i==1.and.j==2).or.(i==2.and.j==1)) then
          if (nxgrid/=1.and.nygrid/=1) then
            fac=(1./60.**2)*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
              45.*((45.*(f(l1+1:l2+1,m+1,n)-f(l1-1:l2-1,m+1,n))  &
                    -9.*(f(l1+2:l2+2,m+1,n)-f(l1-2:l2-2,m+1,n))  &
                       +(f(l1+3:l2+3,m+1,n)-f(l1-3:l2-3,m+1,n))) &
                  -(45.*(f(l1+1:l2+1,m-1,n)-f(l1-1:l2-1,m-1,n))  &
                    -9.*(f(l1+2:l2+2,m-1,n)-f(l1-2:l2-2,m-1,n))  &
                       +(f(l1+3:l2+3,m-1,n)-f(l1-3:l2-3,m-1,n))))&
              -9.*((45.*(f(l1+1:l2+1,m+2,n)-f(l1-1:l2-1,m+2,n))  &
                    -9.*(f(l1+2:l2+2,m+2,n)-f(l1-2:l2-2,m+2,n))  &
                       +(f(l1+3:l2+3,m+2,n)-f(l1-3:l2-3,m+2,n))) &
                  -(45.*(f(l1+1:l2+1,m-2,n)-f(l1-1:l2-1,m-2,n))  &
                    -9.*(f(l1+2:l2+2,m-2,n)-f(l1-2:l2-2,m-2,n))  &
                       +(f(l1+3:l2+3,m-2,n)-f(l1-3:l2-3,m-2,n))))&
                 +((45.*(f(l1+1:l2+1,m+3,n)-f(l1-1:l2-1,m+3,n))  &
                    -9.*(f(l1+2:l2+2,m+3,n)-f(l1-2:l2-2,m+3,n))  &
                       +(f(l1+3:l2+3,m+3,n)-f(l1-3:l2-3,m+3,n))) &
                  -(45.*(f(l1+1:l2+1,m-3,n)-f(l1-1:l2-1,m-3,n))  &
                    -9.*(f(l1+2:l2+2,m-3,n)-f(l1-2:l2-2,m-3,n))  &
                       +(f(l1+3:l2+3,m-3,n)-f(l1-3:l2-3,m-3,n))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=(1./60.**2)*dy_1(m)*dz_1(n)
            df=fac*( &
              45.*((45.*(f(l1:l2,m+1,n+1)-f(l1:l2,m-1,n+1))  &
                    -9.*(f(l1:l2,m+2,n+1)-f(l1:l2,m-2,n+1))  &
                       +(f(l1:l2,m+3,n+1)-f(l1:l2,m-3,n+1))) &
                  -(45.*(f(l1:l2,m+1,n-1)-f(l1:l2,m-1,n-1))  &
                    -9.*(f(l1:l2,m+2,n-1)-f(l1:l2,m-2,n-1))  &
                       +(f(l1:l2,m+3,n-1)-f(l1:l2,m-3,n-1))))&
              -9.*((45.*(f(l1:l2,m+1,n+2)-f(l1:l2,m-1,n+2))  &
                    -9.*(f(l1:l2,m+2,n+2)-f(l1:l2,m-2,n+2))  &
                       +(f(l1:l2,m+3,n+2)-f(l1:l2,m-3,n+2))) &
                  -(45.*(f(l1:l2,m+1,n-2)-f(l1:l2,m-1,n-2))  &
                    -9.*(f(l1:l2,m+2,n-2)-f(l1:l2,m-2,n-2))  &
                       +(f(l1:l2,m+3,n-2)-f(l1:l2,m-3,n-2))))&
                 +((45.*(f(l1:l2,m+1,n+3)-f(l1:l2,m-1,n+3))  &
                    -9.*(f(l1:l2,m+2,n+3)-f(l1:l2,m-2,n+3))  &
                       +(f(l1:l2,m+3,n+3)-f(l1:l2,m-3,n+3))) &
                  -(45.*(f(l1:l2,m+1,n-3)-f(l1:l2,m-1,n-3))  &
                    -9.*(f(l1:l2,m+2,n-3)-f(l1:l2,m-2,n-3))  &
                       +(f(l1:l2,m+3,n-3)-f(l1:l2,m-3,n-3))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./60.**2)*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
              45.*((45.*(f(l1+1:l2+1,m,n+1)-f(l1-1:l2-1,m,n+1))  &
                    -9.*(f(l1+2:l2+2,m,n+1)-f(l1-2:l2-2,m,n+1))  &
                       +(f(l1+3:l2+3,m,n+1)-f(l1-3:l2-3,m,n+1))) &
                  -(45.*(f(l1+1:l2+1,m,n-1)-f(l1-1:l2-1,m,n-1))  &
                    -9.*(f(l1+2:l2+2,m,n-1)-f(l1-2:l2-2,m,n-1))  &
                       +(f(l1+3:l2+3,m,n-1)-f(l1-3:l2-3,m,n-1))))&
              -9.*((45.*(f(l1+1:l2+1,m,n+2)-f(l1-1:l2-1,m,n+2))  &
                    -9.*(f(l1+2:l2+2,m,n+2)-f(l1-2:l2-2,m,n+2))  &
                       +(f(l1+3:l2+3,m,n+2)-f(l1-3:l2-3,m,n+2))) &
                  -(45.*(f(l1+1:l2+1,m,n-2)-f(l1-1:l2-1,m,n-2))  &
                    -9.*(f(l1+2:l2+2,m,n-2)-f(l1-2:l2-2,m,n-2))  &
                       +(f(l1+3:l2+3,m,n-2)-f(l1-3:l2-3,m,n-2))))&
                 +((45.*(f(l1+1:l2+1,m,n+3)-f(l1-1:l2-1,m,n+3))  &
                    -9.*(f(l1+2:l2+2,m,n+3)-f(l1-2:l2-2,m,n+3))  &
                       +(f(l1+3:l2+3,m,n+3)-f(l1-3:l2-3,m,n+3))) &
                  -(45.*(f(l1+1:l2+1,m,n-3)-f(l1-1:l2-1,m,n-3))  &
                    -9.*(f(l1+2:l2+2,m,n-3)-f(l1-2:l2-2,m,n-3))  &
                       +(f(l1+3:l2+3,m,n-3)-f(l1-3:l2-3,m,n-3))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or z-direction'
          endif
        endif
!
      endif                     ! bidiagonal derij
!
    endsubroutine derij_other
!***********************************************************************
endmodule Deriv
