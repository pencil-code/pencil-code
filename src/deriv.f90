! $Id$
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
  public :: initialize_deriv, finalize_deriv
  public :: der, der2, der3, der4, der5, der6, derij, der5i1j, der5_single
  public :: der6_other, der_pencil, der2_pencil, der6_pencil
  public :: der4i2j,der2i2j2k,der3i3j,der3i2j1k,der4i1j1k
  public :: deri_3d_inds
  public :: der_x,der2_x
  public :: der_z,der2_z
  public :: der_upwind1st
  public :: der_onesided_4_slice
  public :: der_onesided_4_slice_other
  public :: der2_minmod
  public :: heatflux_deriv_x
  public :: set_ghosts_for_onesided_ders
  public :: bval_from_neumann, bval_from_3rd, bval_from_4th
!
  real :: der2_coef0, der2_coef1, der2_coef2, der2_coef3
!
!debug  integer, parameter :: icount_der   = 1         !DERCOUNT
!debug  integer, parameter :: icount_der2  = 2         !DERCOUNT
!debug  integer, parameter :: icount_der4  = 3         !DERCOUNT
!debug  integer, parameter :: icount_der5  = 4         !DERCOUNT
!debug  integer, parameter :: icount_der6  = 5         !DERCOUNT
!debug  integer, parameter :: icount_derij = 6         !DERCOUNT
!debug  integer, parameter :: icount_der_upwind1st = 7 !DERCOUNT
!debug  integer, parameter :: icount_der_other = 8     !DERCOUNT
!
  interface der                 ! Overload the der function
    module procedure der_main   ! derivative of an 'mvar' variable
    module procedure der_other  ! derivative of another field
  endinterface
!
  interface der2                 ! Overload the der function
    module procedure der2_main   ! derivative of an 'mvar' variable
    module procedure der2_other  ! derivative of another field
  endinterface
!
  interface derij                 ! Overload the der function
    module procedure derij_main   ! derivative of an 'mvar' variable
    module procedure derij_other  ! derivative of another field
  endinterface
!
  interface  der_onesided_4_slice                ! Overload the der function
    module procedure  der_onesided_4_slice_main  ! derivative of an 'mvar' variable
    module procedure  der_onesided_4_slice_main_pt
    module procedure  der_onesided_4_slice_other ! derivative of another field
    module procedure  der_onesided_4_slice_other_pt
  endinterface
!
  interface bval_from_neumann
    module procedure bval_from_neumann_scl
    module procedure bval_from_neumann_arr
  endinterface
!
  interface bval_from_3rd
    module procedure bval_from_3rd_scl
    module procedure bval_from_3rd_arr
  endinterface
!
  interface bval_from_4th
    module procedure bval_from_4th_scl
    module procedure bval_from_4th_arr
  endinterface
!
  contains
!
!***********************************************************************
    subroutine initialize_deriv
!
!  Initialize stencil coefficients
!
!  23-sep-16/MR: added initialization and manipulation of l-offsets for complete
!                one-sided calculation of 2nd derivatives
!   5-jan-17/MR: removed offset manipulation as not needed; see set_ghosts_for_onesided_ders
!
      use General, only: indgen

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
    subroutine der_main(f, k, df, j, ignoredx)
!
!  calculate derivative df_k/dx_j
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!   1-oct-97/axel: coded
!  18-jul-98/axel: corrected mx -> my and mx -> mz in all y and z ders
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  21-feb-07/axel: added 1/r and 1/pomega factors for non-coord basis
!  20-sep-13/ccyang: added optional argument ignoredx
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: df
      integer, intent(in) :: j, k
      logical, intent(in), optional :: ignoredx
!
      real, parameter :: a = 1.0/60.0
      real, dimension(nx) :: fac
      real :: facs
      logical :: withdx
!
!debug      if (loptimise_ders) der_call_count(k,icount_der,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der,j,1)+1 !DERCOUNT
!
      if (present(ignoredx)) then
        withdx = .not. ignoredx
        if (ignoredx) then
          fac = a; facs = a
        endif
      else
        withdx = .true.
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (withdx) fac = a * dx_1(l1:l2)
          df=fac*(+ 45.0*(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k)) &
                  -  9.0*(f(l1+2:l2+2,m,n,k)-f(l1-2:l2-2,m,n,k)) &
                  +      (f(l1+3:l2+3,m,n,k)-f(l1-3:l2-3,m,n,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (withdx) facs = a * dy_1(m)
          df=facs*(+ 45.0*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                   -  9.0*(f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) &
                   +      (f(l1:l2,m+3,n,k)-f(l1:l2,m-3,n,k)))
          if (withdx) then
            if (lspherical_coords) then
              df = df * r1_mn
            elseif (lcylindrical_coords) then
              df = df * rcyl_mn1
            endif
          endif
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
!!          if (lpole(2) .and. lcoarse) then
          if (withdx) facs = a * dz_1(n)
          df=facs*(+ 45.0*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k)) &
                   -  9.0*(f(l1:l2,m,n+2,k)-f(l1:l2,m,n-2,k)) &
                   +      (f(l1:l2,m,n+3,k)-f(l1:l2,m,n-3,k)))
          if (withdx .and. lspherical_coords) df = df * r1_mn * sin1th(m)
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_main
!***********************************************************************
    subroutine der_x(f,df)
!
!  x derivative operating on an x-dependent 1-D array
!
!  23-jun-15/pete: adapted from der_z; note that f is not the f array!
!
      real, dimension (mx), intent(in)  :: f
      real, dimension (nx), intent(out) :: df
!
      real, dimension (nx) :: fac
!
      if (nxgrid/=1) then
        fac=(1./60)*dx_1(l1:l2)
        df=fac*(+ 45.0*(f(l1+1:l2+1)-f(l1-1:l2-1)) &
                -  9.0*(f(l1+2:l2+2)-f(l1-2:l2-2)) &
                +      (f(l1+3:l2+3)-f(l1-3:l2-3)))
      else
        df=0.
        if (ip<=5) print*, 'der_x: Degenerate case in x-direction'
      endif
!
    endsubroutine der_x
!***********************************************************************
    subroutine der2_x(f,df2)
!
!  Second x derivative operating on an x-dependent 1-D array
!
!  23-jun-15/pete: adapted from der2_z
!
      real, dimension (mx), intent(in)  :: f
      real, dimension (nx), intent(out) :: df2
!
      real, dimension (nx) :: fac
!
      if (nxgrid/=1) then
        fac=(1./180)*dx_1(l1:l2)**2
        df2=fac*(-490.0*f(l1:l2) &
                 +270.0*(f(l1+1:l2+1)+f(l1-1:l2-1)) &
                 - 27.0*(f(l1+2:l2+2)+f(l1-2:l2-2)) &
                 +  2.0*(f(l1+3:l2+3)+f(l1-3:l2-3)))
      else
        df2=0.
        if (ip<=5) print*, 'der2_x: Degenerate case in x-direction'
      endif
!
    endsubroutine der2_x
!***********************************************************************
    subroutine der_z(f,df)
!
!  z derivative operating on a z-dependent 1-D array
!
!   9-feb-07/axel: adapted from der_main; note that f is not the f array!
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(out) :: df
!
      real, dimension (nz) :: fac
!
      if (nzgrid/=1) then
!MR: coarse case/spherical missing!
        fac=(1./60)*dz_1(n1:n2)
        df=fac*(+ 45.0*(f(n1+1:n2+1)-f(n1-1:n2-1)) &
                -  9.0*(f(n1+2:n2+2)-f(n1-2:n2-2)) &
                +      (f(n1+3:n2+3)-f(n1-3:n2-3)))
      else
        df=0.
        if (ip<=5) print*, 'der_z: Degenerate case in z-direction'
      endif
!
    endsubroutine der_z
!***********************************************************************
    subroutine der2_z(f,df2)
!
!  z derivative operating on a z-dependent 1-D array
!
!   2-jan-10/axel: adapted from der_z and der_main
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(out) :: df2
!
      real, dimension (nz) :: fac
!
      if (nzgrid/=1) then
!MR: coarse case/spherical missing!
        fac=(1./180)*dz_1(n1:n2)**2
        df2=fac*(-490.0*f(n1:n2) &
                 +270.0*(f(n1+1:n2+1)+f(n1-1:n2-1)) &
                 - 27.0*(f(n1+2:n2+2)+f(n1-2:n2-2)) &
                 +  2.0*(f(n1+3:n2+3)+f(n1-3:n2-3)))
      else
        df2=0.
        if (ip<=5) print*, 'der2_z: Degenerate case in z-direction'
      endif
!
    endsubroutine der2_z
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
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  21-feb-07/axel: added 1/r and 1/pomega factors for non-coord basis
!  30-sep-16/MR: allow results dimensions > nx
!
      real, dimension (mx,my,mz) :: f
      real, dimension (:) :: df
      integer :: j
!
      intent(in)  :: f,j
      intent(out) :: df

      real, dimension (size(df)) :: fac
      real :: facs
      integer :: l1_,l2_,sdf
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
      else
        sdf=size(df)
        if (sdf>nx) then
          l1_=1; l2_=sdf
        else
          l1_=l1; l2_=l2 
        endif
        if (j==2) then
          if (nygrid/=1) then
            facs=(1./60.)*dy_1(m)
            df=facs*(+ 45.0*(f(l1_:l2_,m+1,n)-f(l1_:l2_,m-1,n)) &
                     -  9.0*(f(l1_:l2_,m+2,n)-f(l1_:l2_,m-2,n)) &
                     +      (f(l1_:l2_,m+3,n)-f(l1_:l2_,m-3,n)))
            if (lspherical_coords)   df=df*r1_mn
            if (lcylindrical_coords) df=df*rcyl_mn1
          else
            df=0.
            if (ip<=5) print*, 'der_other: Degenerate case in y-direction'
          endif
        elseif (j==3) then
          if (nzgrid/=1) then
!!            if (lpole(2) .and. lcoarse) then
            facs = (1./60.) * dz_1(n)
            df=facs*(+ 45.0*(f(l1_:l2_,m,n+1)-f(l1_:l2_,m,n-1)) &
                     -  9.0*(f(l1_:l2_,m,n+2)-f(l1_:l2_,m,n-2)) &
                     +      (f(l1_:l2_,m,n+3)-f(l1_:l2_,m,n-3)))
            if (lspherical_coords) df=df*r1_mn*sin1th(m)
          else
            df=0.
            if (ip<=5) print*, 'der_other: Degenerate case in z-direction'
          endif
        endif
      endif
!
    endsubroutine der_other
!***********************************************************************
    subroutine der_pencil(j,pencil,df)
!
!  Calculate first derivative of any x, y or z pencil.
!
!  01-nov-07/anders: adapted from der
!
      use General, only: itoa

      real, dimension (:) :: pencil,df
      integer :: j
!
      intent(in)  :: j, pencil
      intent(out) :: df
!
!  x-derivative
!
      if (j==1) then
        if (size(pencil)/=mx) &
          call fatal_error('der_pencil', &
                           'pencil must be of size mx for x derivative')
        df(l1:l2)=(1./60)*dx_1(l1:l2)*( &
            + 45.0*(pencil(l1+1:l2+1)-pencil(l1-1:l2-1)) &
            -  9.0*(pencil(l1+2:l2+2)-pencil(l1-2:l2-2)) &
            +      (pencil(l1+3:l2+3)-pencil(l1-3:l2-3)))
!
!  y-derivative
!
      else if (j==2) then
        if (size(pencil)/=my) &
          call fatal_error('der_pencil', &
                           'pencil must be of size my for y derivative')
        df(m1:m2)=(1./60)*dy_1(m1:m2)*( &
            + 45.0*(pencil(m1+1:m2+1)-pencil(m1-1:m2-1)) &
            -  9.0*(pencil(m1+2:m2+2)-pencil(m1-2:m2-2)) &
            +      (pencil(m1+3:m2+3)-pencil(m1-3:m2-3)))
        if (lspherical_coords) then
          df(m1:m2)=df(m1:m2)*r1_mn(lglob)
        elseif (lcylindrical_coords) then
          df(m1:m2)=df(m1:m2)*rcyl_mn1(lglob)
        endif
!
!  z-derivative
!
      else if (j==3) then
        if (size(pencil)/=mz) &
          call fatal_error('der_pencil', &
                           'pencil must be of size mz for z derivative')
!MR: coarse case missing!
        df(n1:n2)=(1./60)*dz_1(n1:n2)*( &
            + 45.0*(pencil(n1+1:n2+1)-pencil(n1-1:n2-1)) &
            -  9.0*(pencil(n1+2:n2+2)-pencil(n1-2:n2-2)) &
            +      (pencil(n1+3:n2+3)-pencil(n1-3:n2-3)))
        if (lspherical_coords) df(n1:n2)=df(n1:n2)*(r1_mn(lglob)*sin1th(m))
      else
        call fatal_error('der_pencil','no such direction j='//trim(itoa(j)))
      endif
!
    endsubroutine der_pencil
!***********************************************************************
    subroutine der2_main(f,k,df2,j,lwo_line_elem)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!  if lwo_line_elem=T no metric coefficents ar multiplied in the denominator;
!  default: lwo_line_elem=F
!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  23-sep-16/MR: introduced offset variables which can be manipulated for complete
!                one-sided calculation of 2nd derivatives
!  20-nov-16/MR: optional parameter lwo_line_elem added
!
      use General, only: loptest

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df2,fac,df
      real :: facs
      integer :: j,k
      logical, optional :: lwo_line_elem
!
      intent(in)  :: f,k,j,lwo_line_elem
      intent(out) :: df2
!
!debug      if (loptimise_ders) der_call_count(k,icount_der2,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der2,j,1) + 1 !DERCOUNT
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=dx_1(l1:l2)**2
          df2=fac*(der2_coef0* f(l1  :l2  ,m,n,k) &
                  +der2_coef1*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  +der2_coef2*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                  +der2_coef3*(f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dx_tilde(l1:l2)*df
          endif
        else
          df2=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          facs=dy_1(m)**2
          df2=facs*(der2_coef0* f(l1:l2,m  ,n,k) &
                   +der2_coef1*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                   +der2_coef2*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                   +der2_coef3*(f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
          if (.not.loptest(lwo_line_elem)) then
            if (lspherical_coords) then
              df2=df2*r2_mn
            elseif (lcylindrical_coords) then
              df2=df2*rcyl_mn2
            endif
          endif
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dy_tilde(m)*df
          endif
        else
          df2=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
!!          if (lpole(2) .and. lcoarse) then
          facs=dz_1(n)**2
          df2=facs*( der2_coef0* f(l1:l2,m,n  ,k) &
                    +der2_coef1*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                    +der2_coef2*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                    +der2_coef3*(f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)) )
          if (.not.loptest(lwo_line_elem).and.lspherical_coords) &
            df2=df2*r2_mn*sin2th(m)
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dz_tilde(n)*df !MR: coarse!?
          endif
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
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df2,fac,df
      real :: facs
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
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dx_tilde(l1:l2)*df
          endif
        else
          df2=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          facs=(1./180)*dy_1(m)**2
          df2=facs*(-490.0*f(l1:l2,m,n) &
                    +270.0*(f(l1:l2,m+1,n)+f(l1:l2,m-1,n)) &
                    - 27.0*(f(l1:l2,m+2,n)+f(l1:l2,m-2,n)) &
                    +  2.0*(f(l1:l2,m+3,n)+f(l1:l2,m-3,n)))
          if (lspherical_coords) then
            df2=df2*r2_mn
          elseif (lcylindrical_coords) then
            df2=df2*rcyl_mn2
          endif
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dy_tilde(m)*df
          endif
        else
          df2=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
!          if (lpole(2) .and. lcoarse) then
          facs=(1./180)*dz_1(n)**2
          df2=facs*(-490.0*f(l1:l2,m,n) &
                    +270.0*(f(l1:l2,m,n+1)+f(l1:l2,m,n-1)) &
                    - 27.0*(f(l1:l2,m,n+2)+f(l1:l2,m,n-2)) &
                    +  2.0*(f(l1:l2,m,n+3)+f(l1:l2,m,n-3)))
          if (lspherical_coords) df2=df2*r2_mn*sin2th(m)
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dz_tilde(n)*df  !MR: coarse?!
          endif
        else
          df2=0.
        endif
      endif
!
    endsubroutine der2_other
!***********************************************************************
    subroutine der2_pencil(j,pencil,df2)
!
!  Calculate 2nd derivative of any x, y or z pencil.
!
!  01-nov-07/anders: adapted from der2
!
      real, dimension (:) :: pencil,df2
      integer :: j
!
      intent(in)  :: j, pencil
      intent(out) :: df2
!
      if (j==1) then
!
!  x-derivative
!
        if (size(pencil)/=mx) &
          call fatal_error('der2_pencil','pencil must be of size mx for x derivative')
        df2=(1./180)*dx_1(l1:l2)**2*(-490.0*pencil(l1:l2) &
               +270.0*(pencil(l1+1:l2+1)+pencil(l1-1:l2-1)) &
               - 27.0*(pencil(l1+2:l2+2)+pencil(l1-2:l2-2)) &
               +  2.0*(pencil(l1+3:l2+3)+pencil(l1-3:l2-3)))
      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) &
          call fatal_error('der2_pencil','pencil must be of size my for y derivative')
!MR: spherical/cylindrical missing
        df2=(1./180)*dy_1(m1:m2)**2*(-490.0*pencil(m1:m2) &
               +270.0*(pencil(m1+1:m2+1)+pencil(m1-1:m2-1)) &
               - 27.0*(pencil(m1+2:m2+2)+pencil(m1-2:m2-2)) &
               +  2.0*(pencil(m1+3:m2+3)+pencil(m1-3:m2-3)))
      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) &
          call fatal_error('der2_pencil','pencil must be of size mz for z derivative')
!MR: spherical/coarse missing
        df2=(1./180)*dz_1(n1:n2)**2*(-490.0*pencil(n1:n2) &
               +270.0*(pencil(n1+1:n2+1)+pencil(n1-1:n2-1)) &
               - 27.0*(pencil(n1+2:n2+2)+pencil(n1-2:n2-2)) &
               +  2.0*(pencil(n1+3:n2+3)+pencil(n1-3:n2-3)))
      else
        if (lroot) print*, 'der2_pencil: no such direction j=', j
        call fatal_error('der2_pencil','')
      endif
!
    endsubroutine der2_pencil
!***********************************************************************
    subroutine der3(f,k,df,j,ignoredx)
!
!  Calculate 3rd derivative of a scalar, get scalar
!
!  10-feb-06/anders: adapted from der5
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
!
      if (.not. lequidist(j)) &
          call fatal_error('der3','NOT IMPLEMENTED for non-equidistant grid')
!
      if (lspherical_coords) &
           call fatal_error('der3','NOT IMPLEMENTED for spherical coordinates')
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
          if (lcylindrical_coords)   df=df*rcyl_mn1**3
!MR: spherical missing
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
!MR: spherical/coarse missing
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df
      real :: fac
      integer :: j,k
      logical, optional :: ignoredx,upwind
      logical :: igndx
!
      intent(in)  :: f,k,j,ignoredx,upwind
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der4,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der4,j,1) + 1 !DERCOUNT
!
      if (lspherical_coords) &
           call fatal_error('der4','NOT IMPLEMENTED for spherical coordinates')
!
      if (present(ignoredx)) then
        igndx = ignoredx
      else
        if (.not. lequidist(j)) &
            call fatal_error('der4','non-equidistant grid works only with IGNOREDX')
        igndx = .false.
      endif
!
      if (present(upwind)) &
        call warning('der4','upwinding not implemented')
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
          if (lcylindrical_coords)   df=df*rcyl_mn1**4
!MR: spherical missing
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=(1.0/6)
          else
            fac=(1.0/6)*1/dz**4
          endif
          df=fac*(+ 56.0* f(l1:l2,m,n  ,k) &
                  - 39.0*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                  + 12.0*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                  -      (f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
!MR: spherical/coarse missing
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
!
      if (.not. lequidist(j)) &
          call fatal_error('der5','NOT IMPLEMENTED for non-equidistant grid')
!
      if (lspherical_coords) &
           call fatal_error('der5','NOT IMPLEMENTED for spherical coordinates')
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
          if (lcylindrical_coords) df=df*rcyl_mn1**5
!MR: spherical missing
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
!MR: spherical/coarse missing
        else
          df=0.
        endif
      endif
!
    endsubroutine der5
!***********************************************************************
    subroutine der6(f,k,df,j,ignoredx,upwind)
!
!  Calculate 6th derivative of a scalar, get scalar.
!  Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!  The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!  The optional flag UPWIND is a variant thereof, which calculates
!  D^(6)*dx^5/60, which is the upwind correction of centered derivatives.
!
!   8-jul-02/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: j,k
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
      real :: facs
!
      intent(in)  :: f,k,j,ignoredx,upwind
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
!
      if (present(upwind)) then
        if (.not. lequidist(j).and..not.lignore_nonequi) &
          call fatal_error('der6','upwind cannot be used with '//&
                           'non-equidistant grid.')
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
            fac=(1.0/60)*dx_1(l1:l2)
          else
            fac=dx_1(l1:l2)**6
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
            facs=1.
          else if (upwnd) then
            facs=(1.0/60)*dy_1(m)
          else
            facs=dy_1(m)**6
          endif
          df=facs*(- 20.0* f(l1:l2,m  ,n,k) &
                   + 15.0*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                   -  6.0*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                   +      (f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
          if ((.not.igndx) .and. (.not.upwnd)) then
            if (lspherical_coords) df = df * r1_mn**6 
            if (lcylindrical_coords) df = df * rcyl_mn1**6 
          endif
        else
          df=0.
        endif        
      elseif (j==3) then
        if (nzgrid/=1) then
!!          if (lpole(2) .and. lcoarse) then
          if (igndx) then
            facs=1.
          else if (upwnd) then
            facs=(1.0/60)*dz_1(n)
          else
            facs=dz_1(n)**6
          endif
          df=facs*(- 20.0* f(l1:l2,m,n  ,k) &
                   + 15.0*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                   -  6.0*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                   +      (f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
          if ((.not.igndx) .and. (.not.upwnd) .and. lspherical_coords) &
            df = df * (r1_mn * sin1th(m))**6
        else
          df=0.
        endif
      endif
!
    endsubroutine der6
!***********************************************************************
    subroutine der2_minmod(f,j,delfk,delfkp1,delfkm1,k)
!
! Calculates gradient of a scalar along the direction j but
! get the derivatives at the point i-1,i,i+1
!
      intent(in) :: f,k,j
      intent(out) :: delfk,delfkp1,delfkm1
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: delfk,delfkp1,delfkm1,fac
      real, dimension (nx,-1:1) :: delf
      real, dimension (0:nx+1) :: delfx
      integer :: j,k
      integer :: i,ii,ix
      real :: tmp1,tmp2,tmp3
! x-component
      select case (k)
        case(1)
          do i=l1-1,l2+1
            ix=i-nghost
            tmp1=f(i,m,n,j)-f(i-1,m,n,j)
            tmp2=(f(i+1,m,n,j)-f(i-1,m,n,j))/4.
            tmp3=f(i+1,m,n,j)-f(i,m,n,j)
            delfx(ix) = minmod(tmp1,tmp2,tmp3)
          enddo
          do i=1,nx;do ii=-1,1
            delf(i,ii) = delfx(i+ii)
          enddo;enddo
          fac = dx_1(l1:l2)
! y-component
        case(2)
          do i=l1,l2
            ix=i-nghost
            do ii=-1,1
              tmp1=f(i,m+ii,n,j)-f(i,m+ii-1,n,j)
              tmp2=(f(i,m+ii+1,n,j)-f(i,m+ii-1,n,j))/4.
              tmp3=f(i,m+ii+1,n,j)-f(i,m+ii,n,j)
              delf(ix,ii) = minmod(tmp1,tmp2,tmp3)*dy_1(i)
            enddo
          enddo
          fac = dy_1(m)
          if (lspherical_coords) then
            fac = fac*r1_mn
          elseif (lcylindrical_coords) then
            fac = fac*rcyl_mn1
          endif
! z-component
        case(3)
          do i=l1,l2
            ix=i-nghost
            do ii=-1,1
              tmp1=f(i,m,n+ii,j)-f(i,m,n+ii-1,j)
              tmp2=(f(i,m,n+ii+1,j)-f(i,m,n+ii-1,j))/4.
              tmp3=f(i,m,n+ii+1,j)-f(i,m,n+ii,j)
              delf(ix,ii) = minmod(tmp1,tmp2,tmp3)
            enddo
          enddo
          fac = dz_1(n)
!MR: coarse missing
          if (lspherical_coords) fac = fac*r1_mn*sin1th(m)
        case default
          call fatal_error('deriv:der2_minmod','wrong component')
        endselect
        delfkm1 = delf(:,-1)*fac
        delfk   = delf(:, 0)*fac
        delfkp1 = delf(:, 1)*fac
!
    endsubroutine der2_minmod
!***********************************************************************
    real function minmod(a,b,c)
!
!  DOCUMENT ME!
!
      real :: a,b,c
!
      if (((a>0) .and. (b>0) .and. (c>0))) then
        minmod=max(a,b,c)
      elseif (((a<0) .and. (b<0) .and. (c<0))) then
        minmod=min(a,b,c)
      else
        minmod=0.0
      endif
!
    endfunction minmod
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
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df,fac
      integer :: j
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
!
      intent(in)  :: f,j,ignoredx,upwind
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
        if (.not. lequidist(j).and..not.lignore_nonequi) &
          call fatal_error('der6_other','upwind cannot be used with '//&
              'non-equidistant grid.')
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
            fac=(1.0/60)*dx_1(l1:l2)
          else
            fac=dx_1(l1:l2)**6
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
            fac=1.
          else if (upwnd) then
            fac=(1.0/60)*dy_1(m)
          else
            fac=dy_1(m)**6
          endif
          df=fac*(- 20.0* f(l1:l2,m  ,n) &
                  + 15.0*(f(l1:l2,m+1,n)+f(l1:l2,m-1,n)) &
                  -  6.0*(f(l1:l2,m+2,n)+f(l1:l2,m-2,n)) &
                  +      (f(l1:l2,m+3,n)+f(l1:l2,m-3,n)))
!MR: spherical/cylndrical missing
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
            fac=dz_1(n)**6
          endif
          df=fac*(- 20.0* f(l1:l2,m,n  ) &
                  + 15.0*(f(l1:l2,m,n+1)+f(l1:l2,m,n-1)) &
                  -  6.0*(f(l1:l2,m,n+2)+f(l1:l2,m,n-2)) &
                  +      (f(l1:l2,m,n+3)+f(l1:l2,m,n-3)))
!MR: spherical/coarse missing
        else
          df=0.
        endif
      endif
!
    endsubroutine der6_other
!***********************************************************************
    subroutine der6_pencil(j,pencil,df6,ignoredx,upwind)
!
!  Calculate 6th derivative of any x, y or z pencil.
!
!  20-jul-20/wlyra: adapted from der2_pencil
!
      real, dimension (:) :: pencil,df6
      real, dimension (nx) :: facx
      real, dimension (ny) :: facy
      real, dimension (nz) :: facz
      integer :: j
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
!
      intent(in)  :: j, pencil,ignoredx,upwind
      intent(out) :: df6

!
      if (present(ignoredx)) then
        igndx = ignoredx
      else
        igndx = .false.
      endif
!     
      if (present(upwind)) then
        if (.not. lequidist(j).and..not.lignore_nonequi) then
          call fatal_error('der6','upwind cannot be used with '//&
              'non-equidistant grid.')
        endif
        upwnd = upwind
      else
        upwnd = .false.
      endif     
!
!  x-derivative
!
      if (j==1) then
        if (size(pencil)/=mx) &
          call fatal_error('der6_pencil','pencil must be of size mx for x derivative')
        if (igndx) then
          facx=1.
        elseif (upwnd) then
          facx=(1.0/60)*dx_1(l1:l2)
        else
          facx=dx_1(l1:l2)**6
        endif
        df6=facx*(- 20.0* pencil(l1:l2) &
                  + 15.0*(pencil(l1+1:l2+1)+pencil(l1-1:l2-1)) &
                  -  6.0*(pencil(l1+2:l2+2)+pencil(l1-2:l2-2)) &
                  +      (pencil(l1+3:l2+3)+pencil(l1-3:l2-3)))
      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) &
          call fatal_error('der6_pencil','pencil must be of size my for y derivative')
        if (igndx) then
          facy=1.
        else if (upwnd) then
          facy=(1.0/60)*dy_1(m1:m2)
        else
          facy=dy_1(m1:m2)**6
        endif
        df6=facy*(- 20.0* pencil(m1:m2) &
                  + 15.0*(pencil(m1+1:m2+1)+pencil(m1-1:m2-1)) &
                  -  6.0*(pencil(m1+2:m2+2)+pencil(m1-2:m2-2)) &
                  +      (pencil(m1+3:m2+3)+pencil(m1-3:m2-3)))
!MR: no spherical/cylindrical
      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) &
          call fatal_error('der6_pencil','pencil must be of size mz for z derivative')
        if (igndx) then
          facz=1.
        else if (upwnd) then
          facz=(1.0/60)*dz_1(n1:n2)
        else
          facz=dz_1(n1:n2)**6
        endif
        df6=facz*(- 20.0* pencil(n1:n2) &
                  + 15.0*(pencil(n1+1:n2+1)+pencil(n1-1:n2-1)) &
                  -  6.0*(pencil(n1+2:n2+2)+pencil(n1-2:n2-2)) &
                  +      (pencil(n1+3:n2+3)+pencil(n1-3:n2-3)))
!MR: no spherical/coarse
      else
        if (lroot) print*, 'der6_pencil: no such direction j=', j
        call fatal_error('der6_pencil','')
      endif
!
    endsubroutine der6_pencil
!***********************************************************************
    real function der5_single(f,j,dc1)
!
!  computes 5th order derivative of function given by f at position j
!
!   3-oct-12/MR: coded
!
      real, dimension(:),  intent(in) :: f, dc1
      integer           ,  intent(in) :: j
!
      real :: fac
!
      if (size(f)/=1) then
        fac=dc1(j)**5
        der5_single=fac*(+  2.5*(f(j+1)-f(j-1)) &
                         -  2.0*(f(j+2)-f(j-2)) &
                         +  0.5*(f(j+3)-f(j-3)))
      else
        der5_single=0.
      endif
!
    endfunction der5_single
!***********************************************************************
    subroutine derij_main(f,k,df,i,j,lwo_line_elem)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!  if lwo_line_elem=T no metric coefficents ar multiplied in the denominator;
!  default: lwo_line_elem=F
!  !!! only for equidistant grids !!!
!
!   8-sep-01/axel: coded
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  14-nov-06/wolf: implemented bidiagonal scheme
!  20-nov-16/MR: optional parameter lwo_line_elem added
!
      use General, only: loptest
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: i,j,k
      logical, optional :: lwo_line_elem
!
      intent(in) :: f,k,i,j
      intent(out) :: df
!
! crash if this is called with i=j
!
!      if (i.eq.j) call fatal_error('derij_main','i=j, no derivative calculated')
!
!debug      if (loptimise_ders) der_call_count(k,icount_derij,i,j) = & !DERCOUNT
!debug                          der_call_count(k,icount_derij,i,j) + 1 !DERCOUNT
!
      if (lbidiagonal_derij) then   !!.and..not.lcoarse) then
        !
        ! Use bidiagonal mixed-derivative operator, i.e.
        ! employ only the three neighbouring points on each of the four
        ! half-diagonals. This gives 6th-order mixed derivatives as the
        ! version below, but involves just 12 points instead of 36.
        !
        if (i+j==3) then
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
        elseif (i+j==5) then
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
        if (i+j==3) then
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
        elseif (i+j==5) then
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
!  Spherical polars. The comments about "minus extra terms" refer to the
!  presence of extra terms that are being evaluated later in gij_etc.
!
      if (loptest(lwo_line_elem)) return

      if (lspherical_coords) then
        if (i+j==3) df=df*r1_mn                     !(minus extra terms)
        if (i==1.and.j==3 .or. i==3.and.j==1) df=df*r1_mn*sin1th(m) !(minus extra terms)
        if (i+j==5) df=df*r2_mn*sin1th(m)           !(minus extra terms)
      elseif (lcylindrical_coords) then
        if ( i==2 .or. j==2 ) df=df*rcyl_mn1
      endif
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
      intent(in)  :: f,i,j
      intent(out) :: df
!
! crash if this is called with i=j
!
!      if (i.eq.j) call fatal_error('derij_main','i=j, no derivative calculated')
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
        if (i+j==3) then
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
        elseif (i+j==5) then
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
        if (i+j==3) then
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
        elseif (i+j==5) then
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
!  Spherical polars. The comments about "minus extra terms" refer to the
!  presence of extra terms that are being evaluated later in gij_etc.
!
      if (lspherical_coords) then
        if (i+j==3) df=df*r1_mn !(minus extra terms)
        if ((i==1.and.j==3) .or. (i==3.and.j==1)) df=df*r1_mn*sin1th(m) !(minus extra terms)               
        if (i+j==5) df=df*r2_mn*sin1th(m) !(minus extra terms)
      elseif (lcylindrical_coords) then
        if ( i+j==3 .or. i+j==5 ) df=df*rcyl_mn1
      endif
!
    endsubroutine derij_other
!***********************************************************************
    subroutine der5i1j(f,k,df,i,j)
!
!  Calculate 6th derivative with respect to two different directions.
!
!  05-dec-06/anders: adapted from derij
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: i,j,k
!
      intent(in) :: f,k,i,j
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_derij,i,j) = & !DERCOUNT
!debug                          der_call_count(k,icount_derij,i,j) + 1 !DERCOUNT
!
!
      if (headtt) then
        if (lspherical_coords.or.lcylindrical_coords) &
             call warning('der5i1j','MISSING CURVATURE TERMS for non-cartesian coordinates')
      endif
!
      if (i==j) then
        if (i==1.and.nxgrid==1 .or. i==2.and.nygrid==1 .or. i==3.and.nzgrid==1) then    
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in '//coornames(i)//'-direction'
        else
          call der6(f,k,df,j)
        endif
      elseif ((i==1.and.j==2)) then
        if (nxgrid/=1.and.nygrid/=1) then
          fac=dx_1(l1:l2)**5*1/60.0*dy_1(m)
          df=fac*( &
            2.5*((45.*(f(l1+1:l2+1,m+1,n,k)-f(l1+1:l2+1,m-1,n,k))  &
                  -9.*(f(l1+1:l2+1,m+2,n,k)-f(l1+1:l2+1,m-2,n,k))  &
                     +(f(l1+1:l2+1,m+3,n,k)-f(l1+1:l2+1,m-3,n,k))) &
                -(45.*(f(l1-1:l2-1,m+1,n,k)-f(l1-1:l2-1,m-1,n,k))  &
                  -9.*(f(l1-1:l2-1,m+2,n,k)-f(l1-1:l2-1,m-2,n,k))  &
                     +(f(l1-1:l2-1,m+3,n,k)-f(l1-1:l2-1,m-3,n,k))))&
           -2.0*((45.*(f(l1+2:l2+2,m+1,n,k)-f(l1+2:l2+2,m-1,n,k))  &
                  -9.*(f(l1+2:l2+2,m+2,n,k)-f(l1+2:l2+2,m-2,n,k))  &
                     +(f(l1+2:l2+2,m+3,n,k)-f(l1+2:l2+2,m-3,n,k))) &
                -(45.*(f(l1-2:l2-2,m+1,n,k)-f(l1-2:l2-2,m-1,n,k))  &
                  -9.*(f(l1-2:l2-2,m+2,n,k)-f(l1-2:l2-2,m-2,n,k))  &
                     +(f(l1-2:l2-2,m+3,n,k)-f(l1-2:l2-2,m-3,n,k))))&
           +0.5*((45.*(f(l1+3:l2+3,m+1,n,k)-f(l1+3:l2+3,m-1,n,k))  &
                  -9.*(f(l1+3:l2+3,m+2,n,k)-f(l1+3:l2+3,m-2,n,k))  &
                     +(f(l1+3:l2+3,m+3,n,k)-f(l1+3:l2+3,m-3,n,k))) &
                -(45.*(f(l1-3:l2-3,m+1,n,k)-f(l1-3:l2-3,m-1,n,k))  &
                  -9.*(f(l1-3:l2-3,m+2,n,k)-f(l1-3:l2-3,m-2,n,k))  &
                     +(f(l1-3:l2-3,m+3,n,k)-f(l1-3:l2-3,m-3,n,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in x- or y-direction'
        endif
      elseif ((i==2.and.j==1)) then
        if (nygrid/=1.and.nxgrid/=1) then
          fac=dy_1(m)**5*1/60.0*dx_1(l1:l2)
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
          fac=dx_1(l1:l2)**5*1/60.0*dz_1(n)
          if (lspherical_coords) fac = fac*r1_mn*sin1th(m)
          df=fac*( &
            2.5*((45.*(f(l1+1:l2+1,m,n+1,k)-f(l1+1:l2+1,m,n-1,k))  &
                  -9.*(f(l1+1:l2+1,m,n+2,k)-f(l1+1:l2+1,m,n-2,k))  &
                     +(f(l1+1:l2+1,m,n+3,k)-f(l1+1:l2+1,m,n-3,k))) &
                -(45.*(f(l1-1:l2-1,m,n+1,k)-f(l1-1:l2-1,m,n-1,k))  &
                  -9.*(f(l1-1:l2-1,m,n+2,k)-f(l1-1:l2-1,m,n-2,k))  &
                     +(f(l1-1:l2-1,m,n+3,k)-f(l1-1:l2-1,m,n-3,k))))&
           -2.0*((45.*(f(l1+2:l2+2,m,n+1,k)-f(l1+2:l2+2,m,n-1,k))  &
                  -9.*(f(l1+2:l2+2,m,n+2,k)-f(l1+2:l2+2,m,n-2,k))  &
                     +(f(l1+2:l2+2,m,n+3,k)-f(l1+2:l2+2,m,n-3,k))) &
                -(45.*(f(l1-2:l2-2,m,n+1,k)-f(l1-2:l2-2,m,n-1,k))  &
                  -9.*(f(l1-2:l2-2,m,n+2,k)-f(l1-2:l2-2,m,n-2,k))  &
                     +(f(l1-2:l2-2,m,n+3,k)-f(l1-2:l2-2,m,n-3,k))))&
           +0.5*((45.*(f(l1+3:l2+3,m,n+1,k)-f(l1+3:l2+3,m,n-1,k))  &
                  -9.*(f(l1+3:l2+3,m,n+2,k)-f(l1+3:l2+3,m,n-2,k))  &
                     +(f(l1+3:l2+3,m,n+3,k)-f(l1+3:l2+3,m,n-3,k))) &
                -(45.*(f(l1-3:l2-3,m,n+1,k)-f(l1-3:l2-3,m,n-1,k))  &
                  -9.*(f(l1-3:l2-3,m,n+2,k)-f(l1-3:l2-3,m,n-2,k))  &
                     +(f(l1-3:l2-3,m,n+3,k)-f(l1-3:l2-3,m,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in x- or z-direction'
        endif
      elseif ((i==3.and.j==1)) then
        if (nzgrid/=1.and.nxgrid/=1) then
          fac=dz_1(n)**5*1/60.0*dx_1(l1:l2)
          df=fac*( &
            2.5*((45.*(f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k))  &
                  -9.*(f(l1+2:l2+2,m,n+1,k)-f(l1-2:l2-2,m,n+1,k))  &
                     +(f(l1+3:l2+3,m,n+1,k)-f(l1-3:l2-3,m,n+1,k))) &
                -(45.*(f(l1+1:l2+1,m,n-1,k)-f(l1-1:l2-1,m,n-1,k))  &
                  -9.*(f(l1+2:l2+2,m,n-1,k)-f(l1-2:l2-2,m,n-1,k))  &
                     +(f(l1+3:l2+3,m,n-1,k)-f(l1-3:l2-3,m,n-1,k))))&
           -2.0*((45.*(f(l1+1:l2+1,m,n+2,k)-f(l1-1:l2-1,m,n+2,k))  &
                  -9.*(f(l1+2:l2+2,m,n+2,k)-f(l1-2:l2-2,m,n+2,k))  &
                     +(f(l1+3:l2+3,m,n+2,k)-f(l1-3:l2-3,m,n+2,k))) &
                -(45.*(f(l1+1:l2+1,m,n-2,k)-f(l1-1:l2-1,m,n-2,k))  &
                  -9.*(f(l1+2:l2+2,m,n-2,k)-f(l1-2:l2-2,m,n-2,k))  &
                     +(f(l1+3:l2+3,m,n-2,k)-f(l1-3:l2-3,m,n-2,k))))&
           +0.5*((45.*(f(l1+1:l2+1,m,n+3,k)-f(l1-1:l2-1,m,n+3,k))  &
                  -9.*(f(l1+2:l2+2,m,n+3,k)-f(l1-2:l2-2,m,n+3,k))  &
                     +(f(l1+3:l2+3,m,n+3,k)-f(l1-3:l2-3,m,n+3,k))) &
                -(45.*(f(l1+1:l2+1,m,n-3,k)-f(l1-1:l2-1,m,n-3,k))  &
                  -9.*(f(l1+2:l2+2,m,n-3,k)-f(l1-2:l2-2,m,n-3,k))  &
                     +(f(l1+3:l2+3,m,n-3,k)-f(l1-3:l2-3,m,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in z- or x-direction'
        endif
      elseif ((i==2.and.j==3)) then
        if (nygrid/=1.and.nzgrid/=1) then
          fac=dy_1(m)**5*1/60.0*dz_1(n)
          df=fac*( &
            2.5*((45.*(f(l1:l2,m+1,n+1,k)-f(l1:l2,m+1,n-1,k))  &
                  -9.*(f(l1:l2,m+1,n+2,k)-f(l1:l2,m+1,n-2,k))  &
                     +(f(l1:l2,m+1,n+3,k)-f(l1:l2,m+1,n-3,k))) &
                -(45.*(f(l1:l2,m-1,n+1,k)-f(l1:l2,m-1,n-1,k))  &
                  -9.*(f(l1:l2,m-1,n+2,k)-f(l1:l2,m-1,n-2,k))  &
                     +(f(l1:l2,m-1,n+3,k)-f(l1:l2,m-1,n-3,k))))&
           -2.0*((45.*(f(l1:l2,m+2,n+1,k)-f(l1:l2,m+2,n-1,k))  &
                  -9.*(f(l1:l2,m+2,n+2,k)-f(l1:l2,m+2,n-2,k))  &
                     +(f(l1:l2,m+2,n+3,k)-f(l1:l2,m+2,n-3,k))) &
                -(45.*(f(l1:l2,m-2,n+1,k)-f(l1:l2,m-2,n-1,k))  &
                  -9.*(f(l1:l2,m-2,n+2,k)-f(l1:l2,m-2,n-2,k))  &
                     +(f(l1:l2,m-2,n+3,k)-f(l1:l2,m-2,n-3,k))))&
           +0.5*((45.*(f(l1:l2,m+3,n+1,k)-f(l1:l2,m+3,n-1,k))  &
                  -9.*(f(l1:l2,m+3,n+2,k)-f(l1:l2,m+3,n-2,k))  &
                     +(f(l1:l2,m+3,n+3,k)-f(l1:l2,m+3,n-3,k))) &
                -(45.*(f(l1:l2,m-3,n+1,k)-f(l1:l2,m-3,n-1,k))  &
                  -9.*(f(l1:l2,m-3,n+2,k)-f(l1:l2,m-3,n-2,k))  &
                     +(f(l1:l2,m-3,n+3,k)-f(l1:l2,m-3,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in y- or z-direction'
        endif
      elseif ((i==3.and.j==2)) then
        if (nzgrid/=1.and.nygrid/=1) then
          fac=dz_1(n)**5*1/60.0*dy_1(m)
          df=fac*( &
            2.5*((45.*(f(l1:l2,m+1,n+1,k)-f(l1:l2,m-1,n+1,k))  &
                  -9.*(f(l1:l2,m+2,n+1,k)-f(l1:l2,m-2,n+1,k))  &
                     +(f(l1:l2,m+3,n+1,k)-f(l1:l2,m-3,n+1,k))) &
                -(45.*(f(l1:l2,m+1,n-1,k)-f(l1:l2,m-1,n-1,k))  &
                  -9.*(f(l1:l2,m+2,n-1,k)-f(l1:l2,m-2,n-1,k))  &
                     +(f(l1:l2,m+3,n-1,k)-f(l1:l2,m-3,n-1,k))))&
           -2.0*((45.*(f(l1:l2,m+1,n+2,k)-f(l1:l2,m-1,n+2,k))  &
                  -9.*(f(l1:l2,m+2,n+2,k)-f(l1:l2,m-2,n+2,k))  &
                     +(f(l1:l2,m+3,n+2,k)-f(l1:l2,m-3,n+2,k))) &
                -(45.*(f(l1:l2,m+1,n-2,k)-f(l1:l2,m-1,n-2,k))  &
                  -9.*(f(l1:l2,m+2,n-2,k)-f(l1:l2,m-2,n-2,k))  &
                     +(f(l1:l2,m+3,n-2,k)-f(l1:l2,m-3,n-2,k))))&
           +0.5*((45.*(f(l1:l2,m+1,n+3,k)-f(l1:l2,m-1,n+3,k))  &
                  -9.*(f(l1:l2,m+2,n+3,k)-f(l1:l2,m-2,n+3,k))  &
                     +(f(l1:l2,m+3,n+3,k)-f(l1:l2,m-3,n+3,k))) &
                -(45.*(f(l1:l2,m+1,n-3,k)-f(l1:l2,m-1,n-3,k))  &
                  -9.*(f(l1:l2,m+2,n-3,k)-f(l1:l2,m-2,n-3,k))  &
                     +(f(l1:l2,m+3,n-3,k)-f(l1:l2,m-3,n-3,k))))&
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
    subroutine der4i2j(f,k,df,i,j)
!
!  Calculate 6th derivative with respect to two different directions.
!
!  02-apr-17/wlyra: adapted from der5i1j
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: i,j,k
!
      intent(in) :: f,k,i,j
      intent(out) :: df

!debug      if (loptimise_ders) der_call_count(k,icount_derij,i,j) = & !DERCOUNT
!debug                          der_call_count(k,icount_derij,i,j) + 1 !DERCOUNT
!
!
      if (headtt) then
        if (lspherical_coords.or.lcylindrical_coords) &
             call warning('der4i2j','MISSING CURVATURE TERMS for non-cartesian coordinates')
      endif
!
      if (i==j) then
        if (i==1.and.nxgrid==1 .or. i==2.and.nygrid==1 .or. i==3.and.nzgrid==1) then
          df=0.
          if (ip<=5) print*, 'der4i2j: Degenerate case in '//coornames(i)//'-direction'
        else
          call der6(f,k,df,j)
        endif
      elseif ((i==1.and.j==2)) then
        if (nxgrid/=1.and.nygrid/=1) then
          fac=(1./6.0*dx_1(l1:l2)**4) * (1./180.0*dy_1(m)**2)
          df=fac*( &
               56.0*( -490.* f(l1:l2,m  ,n,k)                            & 
                      +270.*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k))          &
                       -27.*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k))          &
                        +2.*(f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))         &
              -39.0*((-490.* f(l1+1:l2+1,m  ,n,k)                        &
                      +270.*(f(l1+1:l2+1,m+1,n,k)+f(l1+1:l2+1,m-1,n,k))  &
                       -27.*(f(l1+1:l2+1,m+2,n,k)+f(l1+1:l2+1,m-2,n,k))  &
                        +2.*(f(l1+1:l2+1,m+3,n,k)+f(l1+1:l2+1,m-3,n,k))) &
                    +(-490.* f(l1-1:l2-1,m  ,n,k)                        &
                      +270.*(f(l1-1:l2-1,m+1,n,k)+f(l1-1:l2-1,m-1,n,k))  &
                       -27.*(f(l1-1:l2-1,m+2,n,k)+f(l1-1:l2-1,m-2,n,k))  &
                        +2.*(f(l1-1:l2-1,m+3,n,k)+f(l1-1:l2-1,m-3,n,k))))&
              +12.0*((-490.* f(l1+2:l2+2,m  ,n,k)                        &
                      +270.*(f(l1+2:l2+2,m+1,n,k)+f(l1+2:l2+2,m-1,n,k))  &
                       -27.*(f(l1+2:l2+2,m+2,n,k)+f(l1+2:l2+2,m-2,n,k))  &
                        +2.*(f(l1+2:l2+2,m+3,n,k)+f(l1+2:l2+2,m-3,n,k))) &
                    +(-490.* f(l1-2:l2-2,m  ,n,k)                        &
                      +270.*(f(l1-2:l2-2,m+1,n,k)+f(l1-2:l2-2,m-1,n,k))  &
                       -27.*(f(l1-2:l2-2,m+2,n,k)+f(l1-2:l2-2,m-2,n,k))  &
                        +2.*(f(l1-2:l2-2,m+3,n,k)+f(l1-2:l2-2,m-3,n,k))))&
               -1.0*((-490.* f(l1+3:l2+3,m  ,n,k)                        &
                      +270.*(f(l1+3:l2+3,m+1,n,k)+f(l1+3:l2+3,m-1,n,k))  &
                       -27.*(f(l1+3:l2+3,m+2,n,k)+f(l1+3:l2+3,m-2,n,k))  &
                        +2.*(f(l1+3:l2+3,m+3,n,k)+f(l1+3:l2+3,m-3,n,k))) &
                    +(-490.* f(l1-3:l2-3,m  ,n,k)                        &
                    +  270.*(f(l1-3:l2-3,m+1,n,k)+f(l1-3:l2-3,m-1,n,k))  &
                       -27.*(f(l1-3:l2-3,m+2,n,k)+f(l1-3:l2-3,m-2,n,k))  &
                        +2.*(f(l1-3:l2-3,m+3,n,k)+f(l1-3:l2-3,m-3,n,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der4i2j: Degenerate case in x- or y-direction'
        endif
      elseif ((i==2.and.j==1)) then
        if (nygrid/=1.and.nxgrid/=1) then
          fac=(1./6.0*dy_1(m)**4) * (1./180.0*dx_1(l1:l2)**2)
          df=fac*( &
               56.0*( -490.* f(l1  :l2  ,m  ,n,k)                        & 
                      +270.*(f(l1+1:l2+1,m  ,n,k)+f(l1-1:l2-1,m  ,n,k))  &
                       -27.*(f(l1+2:l2+2,m  ,n,k)+f(l1-2:l2-2,m  ,n,k))  &
                        +2.*(f(l1+3:l2+3,m  ,n,k)+f(l1-3:l2-3,m  ,n,k))) &
              -39.0*((-490.* f(l1  :l2  ,m+1,n,k)                        &
                      +270.*(f(l1+1:l2+1,m+1,n,k)+f(l1-1:l2-1,m+1,n,k))  &
                       -27.*(f(l1+2:l2+2,m+1,n,k)+f(l1-2:l2-2,m+1,n,k))  &
                        +2.*(f(l1+3:l2+3,m+1,n,k)+f(l1-3:l2-3,m+1,n,k))) &
                    +(-490.* f(l1  :l2  ,m-1,n,k)                        &
                      +270.*(f(l1+1:l2+1,m-1,n,k)+f(l1-1:l2-1,m-1,n,k))  &
                       -27.*(f(l1+2:l2+2,m-1,n,k)+f(l1-2:l2-2,m-1,n,k))  &
                        +2.*(f(l1+3:l2+3,m-1,n,k)+f(l1-3:l2-3,m-1,n,k))))&
              +12.0*((-490.* f(l1  :l2  ,m+2,n,k)                        &
                      +270.*(f(l1+1:l2+1,m+2,n,k)+f(l1-1:l2-1,m+2,n,k))  &
                       -27.*(f(l1+2:l2+2,m+2,n,k)+f(l1-2:l2-2,m+2,n,k))  &
                        +2.*(f(l1+3:l2+3,m+2,n,k)+f(l1-3:l2-3,m+2,n,k))) &
                    +(-490.* f(l1  :l2  ,m-2,n,k)                        &
                      +270.*(f(l1+1:l2+1,m-2,n,k)+f(l1-1:l2-1,m-2,n,k))  &
                       -27.*(f(l1+2:l2+2,m-2,n,k)+f(l1-2:l2-2,m-2,n,k))  &
                        +2.*(f(l1+3:l2+3,m-2,n,k)+f(l1-3:l2-3,m-2,n,k))))&
               -1.0*((-490.* f(l1  :l2  ,m+3,n,k)                        &
                      +270.*(f(l1+1:l2+1,m+3,n,k)+f(l1-1:l2-1,m+3,n,k))  &
                       -27.*(f(l1+2:l2+2,m+3,n,k)+f(l1-2:l2-2,m+3,n,k))  &
                        +2.*(f(l1+3:l2+3,m+3,n,k)+f(l1-3:l2-3,m+3,n,k))) &
                    +(-490.* f(l1  :l2  ,m-3,n,k)                        &
                      +270.*(f(l1+1:l2+1,m-3,n,k)+f(l1-1:l2-1,m-3,n,k))  &
                       -27.*(f(l1+2:l2+2,m-3,n,k)+f(l1-2:l2-2,m-3,n,k))  &
                        +2.*(f(l1+3:l2+3,m-3,n,k)+f(l1-3:l2-3,m-3,n,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der4i2j: Degenerate case in y- or x-direction'
        endif
      elseif ((i==1.and.j==3)) then
        if (nxgrid/=1.and.nzgrid/=1) then
          fac=(1./6.0*dx_1(l1:l2)**4) * (1./180.0*dz_1(n)**2)
          if (lspherical_coords) fac = fac*(r1_mn*sin1th(m))**2
          df=fac*( &
               56.0*( -490.* f(l1  :l2  ,m,n  ,k)                        & 
                      +270.*(f(l1  :l2  ,m,n+1,k)+f(l1  :l2  ,m,n-1,k))  &
                       -27.*(f(l1  :l2  ,m,n+2,k)+f(l1  :l2  ,m,n-2,k))  &
                        +2.*(f(l1  :l2  ,m,n+3,k)+f(l1  :l2  ,m,n-3,k))) &
              -39.0*((-490.* f(l1+1:l2+1,m,n  ,k)                        &
                      +270.*(f(l1+1:l2+1,m,n+1,k)+f(l1+1:l2+1,m,n-1,k))  &
                       -27.*(f(l1+1:l2+1,m,n+2,k)+f(l1+1:l2+1,m,n-2,k))  &
                        +2.*(f(l1+1:l2+1,m,n+3,k)+f(l1+1:l2+1,m,n-3,k))) &
                    +(-490.* f(l1-1:l2-1,m,n  ,k)                        &
                      +270.*(f(l1-1:l2-1,m,n+1,k)+f(l1-1:l2-1,m,n-1,k))  &
                       -27.*(f(l1-1:l2-1,m,n+2,k)+f(l1-1:l2-1,m,n-2,k))  &
                        +2.*(f(l1-1:l2-1,m,n+3,k)+f(l1-1:l2-1,m,n-3,k))))&
              +12.0*((-490.* f(l1+2:l2+2,m,n  ,k)                        &
                      +270.*(f(l1+2:l2+2,m,n+1,k)+f(l1+2:l2+2,m,n-1,k))  &
                       -27.*(f(l1+2:l2+2,m,n+2,k)+f(l1+2:l2+2,m,n-2,k))  &
                        +2.*(f(l1+2:l2+2,m,n+3,k)+f(l1+2:l2+2,m,n-3,k))) &
                    +(-490.* f(l1-2:l2-2,m,n  ,k)                        &
                      +270.*(f(l1-2:l2-2,m,n+1,k)+f(l1-2:l2-2,m,n-1,k))  &
                       -27.*(f(l1-2:l2-2,m,n+2,k)+f(l1-2:l2-2,m,n-2,k))  &
                        +2.*(f(l1-2:l2-2,m,n+3,k)+f(l1-2:l2-2,m,n-3,k))))&
               -1.0*((-490.* f(l1+3:l2+3,m,n  ,k)                        &
                      +270.*(f(l1+3:l2+3,m,n+1,k)+f(l1+3:l2+3,m,n-1,k))  &
                       -27.*(f(l1+3:l2+3,m,n+2,k)+f(l1+3:l2+3,m,n-2,k))  &
                        +2.*(f(l1+3:l2+3,m,n+3,k)+f(l1+3:l2+3,m,n-3,k))) &
                    +(-490.* f(l1-3:l2-3,m,n  ,k)                        &
                      +270.*(f(l1-3:l2-3,m,n+1,k)+f(l1-3:l2-3,m,n-1,k))  &
                       -27.*(f(l1-3:l2-3,m,n+2,k)+f(l1-3:l2-3,m,n-2,k))  &
                        +2.*(f(l1-3:l2-3,m,n+3,k)+f(l1-3:l2-3,m,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der4i2j: Degenerate case in x- or z-direction'
        endif
      elseif ((i==3.and.j==1)) then
        if (nzgrid/=1.and.nxgrid/=1) then
          fac=(1./6.0*dz_1(n)**4) * (1./180.0*dx_1(l1:l2)**2)
          if (lspherical_coords) fac = fac*(r1_mn*sin1th(m))**4
          df=fac*( &
               56.0*( -490.* f(l1  :l2  ,m,n  ,k)                        & 
                      +270.*(f(l1+1:l2+1,m,n  ,k)+f(l1-1:l2-1,m,n  ,k))  &
                       -27.*(f(l1+2:l2+2,m,n  ,k)+f(l1-2:l2-2,m,n  ,k))  &
                        +2.*(f(l1+3:l2+3,m,n  ,k)+f(l1-3:l2-3,m,n  ,k))) &
              -39.0*((-490.* f(l1  :l2  ,m,n+1,k)                        &
                      +270.*(f(l1+1:l2+1,m,n+1,k)+f(l1-1:l2-1,m,n+1,k))  &
                       -27.*(f(l1+2:l2+2,m,n+1,k)+f(l1-2:l2-2,m,n+1,k))  &
                        +2.*(f(l1+3:l2+3,m,n+1,k)+f(l1-3:l2-3,m,n+1,k))) &
                    +(-490.* f(l1  :l2  ,m,n-1,k)                        &
                      +270.*(f(l1+1:l2+1,m,n-1,k)+f(l1-1:l2-1,m,n-1,k))  &
                       -27.*(f(l1+2:l2+2,m,n-1,k)+f(l1-2:l2-2,m,n-1,k))  &
                        +2.*(f(l1+3:l2+3,m,n-1,k)+f(l1-3:l2-3,m,n-1,k))))&
              +12.0*((-490.* f(l1  :l2  ,m,n+2,k)                        &
                      +270.*(f(l1+1:l2+1,m,n+2,k)+f(l1-1:l2-1,m,n+2,k))  &
                       -27.*(f(l1+2:l2+2,m,n+2,k)+f(l1-2:l2-2,m,n+2,k))  &
                        +2.*(f(l1+3:l2+3,m,n+2,k)+f(l1-3:l2-3,m,n+2,k))) &
                    +(-490.* f(l1  :l2  ,m,n-2,k)                        &
                      +270.*(f(l1+1:l2+1,m,n-2,k)+f(l1-1:l2-1,m,n-2,k))  &
                       -27.*(f(l1+2:l2+2,m,n-2,k)+f(l1-2:l2-2,m,n-2,k))  &
                        +2.*(f(l1+3:l2+3,m,n-2,k)+f(l1-3:l2-3,m,n-2,k))))&
               -1.0*((-490.* f(l1  :l2  ,m,n+3,k)                        &
                      +270.*(f(l1+1:l2+1,m,n+3,k)+f(l1-1:l2-1,m,n+3,k))  &
                       -27.*(f(l1+2:l2+2,m,n+3,k)+f(l1-2:l2-2,m,n+3,k))  &
                        +2.*(f(l1+3:l2+3,m,n+3,k)+f(l1-3:l2-3,m,n+3,k))) &
                    +(-490.* f(l1  :l2  ,m,n-3,k)                        &
                      +270.*(f(l1+1:l2+1,m,n-3,k)+f(l1-1:l2-1,m,n-3,k))  &
                       -27.*(f(l1+2:l2+2,m,n-3,k)+f(l1-2:l2-2,m,n-3,k))  &
                        +2.*(f(l1+3:l2+3,m,n-3,k)+f(l1-3:l2-3,m,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der4i2j: Degenerate case in z- or x-direction'
        endif
      elseif ((i==2.and.j==3)) then
        if (nygrid/=1.and.nzgrid/=1) then
          fac=(1./6.0*dy_1(m)**4) * (1./180.0*dz_1(n)**2)
          df=fac*( &
               56.0*( -490.* f(l1:l2,m  ,n  ,k)                      & 
                      +270.*(f(l1:l2,m  ,n+1,k)+f(l1:l2,m  ,n-1,k))  &
                       -27.*(f(l1:l2,m  ,n+2,k)+f(l1:l2,m  ,n-2,k))  &
                        +2.*(f(l1:l2,m  ,n+3,k)+f(l1:l2,m  ,n-3,k))) &
              -39.0*((-490.* f(l1:l2,m+1,n  ,k)                      &
                      +270.*(f(l1:l2,m+1,n+1,k)+f(l1:l2,m+1,n-1,k))  &
                       -27.*(f(l1:l2,m+1,n+2,k)+f(l1:l2,m+1,n-2,k))  &
                        +2.*(f(l1:l2,m+1,n+3,k)+f(l1:l2,m+1,n-3,k))) &
                    +(-490.* f(l1:l2,m-1,n  ,k)                      &
                      +270.*(f(l1:l2,m-1,n+1,k)+f(l1:l2,m-1,n-1,k))  &
                       -27.*(f(l1:l2,m-1,n+2,k)+f(l1:l2,m-1,n-2,k))  &
                        +2.*(f(l1:l2,m-1,n+3,k)+f(l1:l2,m-1,n-3,k))))&
              +12.0*((-490.* f(l1:l2,m+2,n  ,k)                      &
                      +270.*(f(l1:l2,m+2,n+1,k)+f(l1:l2,m+2,n-1,k))  &
                       -27.*(f(l1:l2,m+2,n+2,k)+f(l1:l2,m+2,n-2,k))  &
                        +2.*(f(l1:l2,m+2,n+3,k)+f(l1:l2,m+2,n-3,k))) &
                    +(-490.* f(l1:l2,m-2,n  ,k)                      &
                      +270.*(f(l1:l2,m-2,n+1,k)+f(l1:l2,m-2,n-1,k))  &
                       -27.*(f(l1:l2,m-2,n+2,k)+f(l1:l2,m-2,n-2,k))  &
                        +2.*(f(l1:l2,m-2,n+3,k)+f(l1:l2,m-2,n-3,k))))&
               -1.0*((-490.* f(l1:l2,m+3,n  ,k)                      &
                      +270.*(f(l1:l2,m+3,n+1,k)+f(l1:l2,m+3,n-1,k))  &
                       -27.*(f(l1:l2,m+3,n+2,k)+f(l1:l2,m+3,n-2,k))  &
                        +2.*(f(l1:l2,m+3,n+3,k)+f(l1:l2,m+3,n-3,k))) &
                    +(-490.* f(l1:l2,m-3,n  ,k)                      &
                      +270.*(f(l1:l2,m-3,n+1,k)+f(l1:l2,m-3,n-1,k))  &
                       -27.*(f(l1:l2,m-3,n+2,k)+f(l1:l2,m-3,n-2,k))  &
                        +2.*(f(l1:l2,m-3,n+3,k)+f(l1:l2,m-3,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der4i2j: Degenerate case in y- or z-direction'
        endif
      elseif ((i==3.and.j==2)) then
        if (nzgrid/=1.and.nygrid/=1) then
          fac=(1./6.0*dz_1(n)**4) * (1./180.0*dy_1(m)**2)
          df=fac*( &
               56.0*( -490.* f(l1:l2,m  ,n  ,k)                      & 
                      +270.*(f(l1:l2,m+1,n  ,k)+f(l1:l2,m-1,n  ,k))  &
                       -27.*(f(l1:l2,m+2,n  ,k)+f(l1:l2,m-2,n  ,k))  &
                        +2.*(f(l1:l2,m+3,n  ,k)+f(l1:l2,m-3,n  ,k))) &
              -39.0*((-490.* f(l1:l2,m  ,n+1,k)                      &
                      +270.*(f(l1:l2,m+1,n+1,k)+f(l1:l2,m-1,n+1,k))  &
                       -27.*(f(l1:l2,m+2,n+1,k)+f(l1:l2,m-2,n+1,k))  &
                        +2.*(f(l1:l2,m+3,n+1,k)+f(l1:l2,m-3,n+1,k))) &
                    +(-490.* f(l1:l2,m  ,n-1,k)                      &
                      +270.*(f(l1:l2,m+1,n-1,k)+f(l1:l2,m-1,n-1,k))  &
                       -27.*(f(l1:l2,m+2,n-1,k)+f(l1:l2,m-2,n-1,k))  &
                        +2.*(f(l1:l2,m+3,n-1,k)+f(l1:l2,m-3,n-1,k))))&
              +12.0*((-490.* f(l1:l2,m  ,n+2,k)                      &
                      +270.*(f(l1:l2,m+1,n+2,k)+f(l1:l2,m-1,n+2,k))  &
                       -27.*(f(l1:l2,m+2,n+2,k)+f(l1:l2,m-2,n+2,k))  &
                        +2.*(f(l1:l2,m+3,n+2,k)+f(l1:l2,m-3,n+2,k))) &
                    +(-490.* f(l1:l2,m  ,n-2,k)                      &
                      +270.*(f(l1:l2,m+1,n-2,k)+f(l1:l2,m-1,n-2,k))  &
                       -27.*(f(l1:l2,m+2,n-2,k)+f(l1:l2,m-2,n-2,k))  &
                        +2.*(f(l1:l2,m+3,n-2,k)+f(l1:l2,m-3,n-2,k))))&
               -1.0*((-490.* f(l1:l2,m  ,n+3,k)                      &
                      +270.*(f(l1:l2,m+1,n+3,k)+f(l1:l2,m-1,n+3,k))  &
                       -27.*(f(l1:l2,m+2,n+3,k)+f(l1:l2,m-2,n+3,k))  &
                        +2.*(f(l1:l2,m+3,n+3,k)+f(l1:l2,m-3,n+3,k))) &
                    +(-490.* f(l1:l2,m  ,n-3,k)                      &
                      +270.*(f(l1:l2,m+1,n-3,k)+f(l1:l2,m-1,n-3,k))  &
                       -27.*(f(l1:l2,m+2,n-3,k)+f(l1:l2,m-2,n-3,k))  &
                        +2.*(f(l1:l2,m+3,n-3,k)+f(l1:l2,m-3,n-3,k))))&
                        )
        else
          df=0.
          if (ip<=5) print*, 'der4i2j: Degenerate case in z- or y-direction'
        endif
      else
        print*, 'der4i2j: no such value for i,j=', i, j
        call fatal_error('der4i2j','')
      endif
!
    endsubroutine der4i2j
!***********************************************************************
    subroutine der2i2j2k(f,k,df)
!      
!  Mixed 6th derivative of der2x(der2y(der2z(f))). Worked out symbolically
!  in python. Result as spit from the python routine.       
!
!  02-apr-17/wlyra: coded
!
      real, dimension (mx,my,mz,mfarray),intent(in) :: f
      real, dimension (nx) :: fac
      integer,intent(in) :: k
      real, dimension(nx), intent(out) :: df
!
      if (headtt) then
        if (lspherical_coords.or.lcylindrical_coords) &
             call warning('der2i2j2k','MISSING CURVATURE TERMS for non-cartesian coordinates')
      endif
!
!  MR: cases i=j/=k etc. missing
!
      if (nxgrid/=1.and.nzgrid/=1.and.nygrid/=1) then
        fac=1./180.0**3*(dx_1(l1:l2)*dy_1(m)*dz_1(n))**2
        df = fac*(& 
             ( -117649000.0 * f( l1:l2 , m , n , k ))+&
             ( 64827000.0 * f( l1:l2 , m , n - 1 , k ))+&
             ( -6482700.0 * f( l1:l2 , m , n - 2 , k ))+&
             ( 480200.0 * f( l1:l2 , m , n - 3 , k ))+&
             ( 64827000.0 * f( l1:l2 , m , n + 1 , k ))+&
             ( -6482700.0 * f( l1:l2 , m , n + 2 , k ))+&
             ( 480200.0 * f( l1:l2 , m , n + 3 , k ))+&
             ( 64827000.0 * f( l1:l2 , m - 1 , n , k ))+&
             ( -35721000.0 * f( l1:l2 , m - 1 , n - 1 , k ))+&
             ( 3572100.0 * f( l1:l2 , m - 1 , n - 2 , k ))+&
             ( -264600.0 * f( l1:l2 , m - 1 , n - 3 , k ))+&
             ( -35721000.0 * f( l1:l2 , m - 1 , n + 1 , k ))+&
             ( 3572100.0 * f( l1:l2 , m - 1 , n + 2 , k ))+&
             ( -264600.0 * f( l1:l2 , m - 1 , n + 3 , k ))+&
             ( -6482700.0 * f( l1:l2 , m - 2 , n , k ))+&
             ( 3572100.0 * f( l1:l2 , m - 2 , n - 1 , k ))+&
             ( -357210.0 * f( l1:l2 , m - 2 , n - 2 , k ))+&
             ( 26460.0 * f( l1:l2 , m - 2 , n - 3 , k ))+&
             ( 3572100.0 * f( l1:l2 , m - 2 , n + 1 , k ))+&
             ( -357210.0 * f( l1:l2 , m - 2 , n + 2 , k ))+&
             ( 26460.0 * f( l1:l2 , m - 2 , n + 3 , k ))+&
             ( 480200.0 * f( l1:l2 , m - 3 , n , k ))+&
             ( -264600.0 * f( l1:l2 , m - 3 , n - 1 , k ))+&
             ( 26460.0 * f( l1:l2 , m - 3 , n - 2 , k ))+&
             ( -1960.0 * f( l1:l2 , m - 3 , n - 3 , k ))+&
             ( -264600.0 * f( l1:l2 , m - 3 , n + 1 , k ))+&
             ( 26460.0 * f( l1:l2 , m - 3 , n + 2 , k ))+&
             ( -1960.0 * f( l1:l2 , m - 3 , n + 3 , k ))+&
             ( 64827000.0 * f( l1:l2 , m + 1 , n , k ))+&
             ( -35721000.0 * f( l1:l2 , m + 1 , n - 1 , k ))+&
             ( 3572100.0 * f( l1:l2 , m + 1 , n - 2 , k ))+&
             ( -264600.0 * f( l1:l2 , m + 1 , n - 3 , k ))+&
             ( -35721000.0 * f( l1:l2 , m + 1 , n + 1 , k ))+&
             ( 3572100.0 * f( l1:l2 , m + 1 , n + 2 , k ))+&
             ( -264600.0 * f( l1:l2 , m + 1 , n + 3 , k ))+&
             ( -6482700.0 * f( l1:l2 , m + 2 , n , k ))+&
             ( 3572100.0 * f( l1:l2 , m + 2 , n - 1 , k ))+&
             ( -357210.0 * f( l1:l2 , m + 2 , n - 2 , k ))+&
             ( 26460.0 * f( l1:l2 , m + 2 , n - 3 , k ))+&
             ( 3572100.0 * f( l1:l2 , m + 2 , n + 1 , k ))+&
             ( -357210.0 * f( l1:l2 , m + 2 , n + 2 , k ))+&
             ( 26460.0 * f( l1:l2 , m + 2 , n + 3 , k ))+&
             ( 480200.0 * f( l1:l2 , m + 3 , n , k ))+&
             ( -264600.0 * f( l1:l2 , m + 3 , n - 1 , k ))+&
             ( 26460.0 * f( l1:l2 , m + 3 , n - 2 , k ))+&
             ( -1960.0 * f( l1:l2 , m + 3 , n - 3 , k ))+&
             ( -264600.0 * f( l1:l2 , m + 3 , n + 1 , k ))+&
             ( 26460.0 * f( l1:l2 , m + 3 , n + 2 , k ))+&
             ( -1960.0 * f( l1:l2 , m + 3 , n + 3 , k ))+&
             ( 64827000.0 * f( l1-1:l2-1 , m , n , k ))+&
             ( -35721000.0 * f( l1-1:l2-1 , m , n - 1 , k ))+&
             ( 3572100.0 * f( l1-1:l2-1 , m , n - 2 , k ))+&
             ( -264600.0 * f( l1-1:l2-1 , m , n - 3 , k ))+&
             ( -35721000.0 * f( l1-1:l2-1 , m , n + 1 , k ))+&
             ( 3572100.0 * f( l1-1:l2-1 , m , n + 2 , k ))+&
             ( -264600.0 * f( l1-1:l2-1 , m , n + 3 , k ))+&
             ( -35721000.0 * f( l1-1:l2-1 , m - 1 , n , k ))+&
             ( 19683000.0 * f( l1-1:l2-1 , m - 1 , n - 1 , k ))+&
             ( -1968300.0 * f( l1-1:l2-1 , m - 1 , n - 2 , k ))+&
             ( 145800.0 * f( l1-1:l2-1 , m - 1 , n - 3 , k ))+&
             ( 19683000.0 * f( l1-1:l2-1 , m - 1 , n + 1 , k ))+&
             ( -1968300.0 * f( l1-1:l2-1 , m - 1 , n + 2 , k ))+&
             ( 145800.0 * f( l1-1:l2-1 , m - 1 , n + 3 , k ))+&
             ( 3572100.0 * f( l1-1:l2-1 , m - 2 , n , k ))+&
             ( -1968300.0 * f( l1-1:l2-1 , m - 2 , n - 1 , k ))+&
             ( 196830.0 * f( l1-1:l2-1 , m - 2 , n - 2 , k ))+&
             ( -14580.0 * f( l1-1:l2-1 , m - 2 , n - 3 , k ))+&
             ( -1968300.0 * f( l1-1:l2-1 , m - 2 , n + 1 , k ))+&
             ( 196830.0 * f( l1-1:l2-1 , m - 2 , n + 2 , k ))+&
             ( -14580.0 * f( l1-1:l2-1 , m - 2 , n + 3 , k ))+&
             ( -264600.0 * f( l1-1:l2-1 , m - 3 , n , k ))+&
             ( 145800.0 * f( l1-1:l2-1 , m - 3 , n - 1 , k ))+&
             ( -14580.0 * f( l1-1:l2-1 , m - 3 , n - 2 , k ))+&
             ( 1080.0 * f( l1-1:l2-1 , m - 3 , n - 3 , k ))+&
             ( 145800.0 * f( l1-1:l2-1 , m - 3 , n + 1 , k ))+&
             ( -14580.0 * f( l1-1:l2-1 , m - 3 , n + 2 , k ))+&
             ( 1080.0 * f( l1-1:l2-1 , m - 3 , n + 3 , k ))+&
             ( -35721000.0 * f( l1-1:l2-1 , m + 1 , n , k ))+&
             ( 19683000.0 * f( l1-1:l2-1 , m + 1 , n - 1 , k ))+&
             ( -1968300.0 * f( l1-1:l2-1 , m + 1 , n - 2 , k ))+&
             ( 145800.0 * f( l1-1:l2-1 , m + 1 , n - 3 , k ))+&
             ( 19683000.0 * f( l1-1:l2-1 , m + 1 , n + 1 , k ))+&
             ( -1968300.0 * f( l1-1:l2-1 , m + 1 , n + 2 , k ))+&
             ( 145800.0 * f( l1-1:l2-1 , m + 1 , n + 3 , k ))+&
             ( 3572100.0 * f( l1-1:l2-1 , m + 2 , n , k ))+&
             ( -1968300.0 * f( l1-1:l2-1 , m + 2 , n - 1 , k ))+&
             ( 196830.0 * f( l1-1:l2-1 , m + 2 , n - 2 , k ))+&
             ( -14580.0 * f( l1-1:l2-1 , m + 2 , n - 3 , k ))+&
             ( -1968300.0 * f( l1-1:l2-1 , m + 2 , n + 1 , k ))+&
             ( 196830.0 * f( l1-1:l2-1 , m + 2 , n + 2 , k ))+&
             ( -14580.0 * f( l1-1:l2-1 , m + 2 , n + 3 , k ))+&
             ( -264600.0 * f( l1-1:l2-1 , m + 3 , n , k ))+&
             ( 145800.0 * f( l1-1:l2-1 , m + 3 , n - 1 , k ))+&
             ( -14580.0 * f( l1-1:l2-1 , m + 3 , n - 2 , k ))+&
             ( 1080.0 * f( l1-1:l2-1 , m + 3 , n - 3 , k ))+&
             ( 145800.0 * f( l1-1:l2-1 , m + 3 , n + 1 , k ))+&
             ( -14580.0 * f( l1-1:l2-1 , m + 3 , n + 2 , k ))+&
             ( 1080.0 * f( l1-1:l2-1 , m + 3 , n + 3 , k ))+&
             ( -6482700.0 * f( l1-2:l2-2 , m , n , k ))+&
             ( 3572100.0 * f( l1-2:l2-2 , m , n - 1 , k ))+&
             ( -357210.0 * f( l1-2:l2-2 , m , n - 2 , k ))+&
             ( 26460.0 * f( l1-2:l2-2 , m , n - 3 , k ))+&
             ( 3572100.0 * f( l1-2:l2-2 , m , n + 1 , k ))+&
             ( -357210.0 * f( l1-2:l2-2 , m , n + 2 , k ))+&
             ( 26460.0 * f( l1-2:l2-2 , m , n + 3 , k ))+&
             ( 3572100.0 * f( l1-2:l2-2 , m - 1 , n , k ))+&
             ( -1968300.0 * f( l1-2:l2-2 , m - 1 , n - 1 , k ))+&
             ( 196830.0 * f( l1-2:l2-2 , m - 1 , n - 2 , k ))+&
             ( -14580.0 * f( l1-2:l2-2 , m - 1 , n - 3 , k ))+&
             ( -1968300.0 * f( l1-2:l2-2 , m - 1 , n + 1 , k ))+&
             ( 196830.0 * f( l1-2:l2-2 , m - 1 , n + 2 , k ))+&
             ( -14580.0 * f( l1-2:l2-2 , m - 1 , n + 3 , k ))+&
             ( -357210.0 * f( l1-2:l2-2 , m - 2 , n , k ))+&
             ( 196830.0 * f( l1-2:l2-2 , m - 2 , n - 1 , k ))+&
             ( -19683.0 * f( l1-2:l2-2 , m - 2 , n - 2 , k ))+&
             ( 1458.0 * f( l1-2:l2-2 , m - 2 , n - 3 , k ))+&
             ( 196830.0 * f( l1-2:l2-2 , m - 2 , n + 1 , k ))+&
             ( -19683.0 * f( l1-2:l2-2 , m - 2 , n + 2 , k ))+&
             ( 1458.0 * f( l1-2:l2-2 , m - 2 , n + 3 , k ))+&
             ( 26460.0 * f( l1-2:l2-2 , m - 3 , n , k ))+&
             ( -14580.0 * f( l1-2:l2-2 , m - 3 , n - 1 , k ))+&
             ( 1458.0 * f( l1-2:l2-2 , m - 3 , n - 2 , k ))+&
             ( -108.0 * f( l1-2:l2-2 , m - 3 , n - 3 , k ))+&
             ( -14580.0 * f( l1-2:l2-2 , m - 3 , n + 1 , k ))+&
             ( 1458.0 * f( l1-2:l2-2 , m - 3 , n + 2 , k ))+&
             ( -108.0 * f( l1-2:l2-2 , m - 3 , n + 3 , k ))+&
             ( 3572100.0 * f( l1-2:l2-2 , m + 1 , n , k ))+&
             ( -1968300.0 * f( l1-2:l2-2 , m + 1 , n - 1 , k ))+&
             ( 196830.0 * f( l1-2:l2-2 , m + 1 , n - 2 , k ))+&
             ( -14580.0 * f( l1-2:l2-2 , m + 1 , n - 3 , k ))+&
             ( -1968300.0 * f( l1-2:l2-2 , m + 1 , n + 1 , k ))+&
             ( 196830.0 * f( l1-2:l2-2 , m + 1 , n + 2 , k ))+&
             ( -14580.0 * f( l1-2:l2-2 , m + 1 , n + 3 , k ))+&
             ( -357210.0 * f( l1-2:l2-2 , m + 2 , n , k ))+&
             ( 196830.0 * f( l1-2:l2-2 , m + 2 , n - 1 , k ))+&
             ( -19683.0 * f( l1-2:l2-2 , m + 2 , n - 2 , k ))+&
             ( 1458.0 * f( l1-2:l2-2 , m + 2 , n - 3 , k ))+&
             ( 196830.0 * f( l1-2:l2-2 , m + 2 , n + 1 , k ))+&
             ( -19683.0 * f( l1-2:l2-2 , m + 2 , n + 2 , k ))+&
             ( 1458.0 * f( l1-2:l2-2 , m + 2 , n + 3 , k ))+&
             ( 26460.0 * f( l1-2:l2-2 , m + 3 , n , k ))+&
             ( -14580.0 * f( l1-2:l2-2 , m + 3 , n - 1 , k ))+&
             ( 1458.0 * f( l1-2:l2-2 , m + 3 , n - 2 , k ))+&
             ( -108.0 * f( l1-2:l2-2 , m + 3 , n - 3 , k ))+&
             ( -14580.0 * f( l1-2:l2-2 , m + 3 , n + 1 , k ))+&
             ( 1458.0 * f( l1-2:l2-2 , m + 3 , n + 2 , k ))+&
             ( -108.0 * f( l1-2:l2-2 , m + 3 , n + 3 , k ))+&
             ( 480200.0 * f( l1-3:l2-3 , m , n , k ))+&
             ( -264600.0 * f( l1-3:l2-3 , m , n - 1 , k ))+&
             ( 26460.0 * f( l1-3:l2-3 , m , n - 2 , k ))+&
             ( -1960.0 * f( l1-3:l2-3 , m , n - 3 , k ))+&
             ( -264600.0 * f( l1-3:l2-3 , m , n + 1 , k ))+&
             ( 26460.0 * f( l1-3:l2-3 , m , n + 2 , k ))+&
             ( -1960.0 * f( l1-3:l2-3 , m , n + 3 , k ))+&
             ( -264600.0 * f( l1-3:l2-3 , m - 1 , n , k ))+&
             ( 145800.0 * f( l1-3:l2-3 , m - 1 , n - 1 , k ))+&
             ( -14580.0 * f( l1-3:l2-3 , m - 1 , n - 2 , k ))+&
             ( 1080.0 * f( l1-3:l2-3 , m - 1 , n - 3 , k ))+&
             ( 145800.0 * f( l1-3:l2-3 , m - 1 , n + 1 , k ))+&
             ( -14580.0 * f( l1-3:l2-3 , m - 1 , n + 2 , k ))+&
             ( 1080.0 * f( l1-3:l2-3 , m - 1 , n + 3 , k ))+&
             ( 26460.0 * f( l1-3:l2-3 , m - 2 , n , k ))+&
             ( -14580.0 * f( l1-3:l2-3 , m - 2 , n - 1 , k ))+&
             ( 1458.0 * f( l1-3:l2-3 , m - 2 , n - 2 , k ))+&
             ( -108.0 * f( l1-3:l2-3 , m - 2 , n - 3 , k ))+&
             ( -14580.0 * f( l1-3:l2-3 , m - 2 , n + 1 , k ))+&
             ( 1458.0 * f( l1-3:l2-3 , m - 2 , n + 2 , k ))+&
             ( -108.0 * f( l1-3:l2-3 , m - 2 , n + 3 , k ))+&
             ( -1960.0 * f( l1-3:l2-3 , m - 3 , n , k ))+&
             ( 1080.0 * f( l1-3:l2-3 , m - 3 , n - 1 , k ))+&
             ( -108.0 * f( l1-3:l2-3 , m - 3 , n - 2 , k ))+&
             ( 8.0 * f( l1-3:l2-3 , m - 3 , n - 3 , k ))+&
             ( 1080.0 * f( l1-3:l2-3 , m - 3 , n + 1 , k ))+&
             ( -108.0 * f( l1-3:l2-3 , m - 3 , n + 2 , k ))+&
             ( 8.0 * f( l1-3:l2-3 , m - 3 , n + 3 , k ))+&
             ( -264600.0 * f( l1-3:l2-3 , m + 1 , n , k ))+&
             ( 145800.0 * f( l1-3:l2-3 , m + 1 , n - 1 , k ))+&
             ( -14580.0 * f( l1-3:l2-3 , m + 1 , n - 2 , k ))+&
             ( 1080.0 * f( l1-3:l2-3 , m + 1 , n - 3 , k ))+&
             ( 145800.0 * f( l1-3:l2-3 , m + 1 , n + 1 , k ))+&
             ( -14580.0 * f( l1-3:l2-3 , m + 1 , n + 2 , k ))+&
             ( 1080.0 * f( l1-3:l2-3 , m + 1 , n + 3 , k ))+&
             ( 26460.0 * f( l1-3:l2-3 , m + 2 , n , k ))+&
             ( -14580.0 * f( l1-3:l2-3 , m + 2 , n - 1 , k ))+&
             ( 1458.0 * f( l1-3:l2-3 , m + 2 , n - 2 , k ))+&
             ( -108.0 * f( l1-3:l2-3 , m + 2 , n - 3 , k ))+&
             ( -14580.0 * f( l1-3:l2-3 , m + 2 , n + 1 , k ))+&
             ( 1458.0 * f( l1-3:l2-3 , m + 2 , n + 2 , k ))+&
             ( -108.0 * f( l1-3:l2-3 , m + 2 , n + 3 , k ))+&
             ( -1960.0 * f( l1-3:l2-3 , m + 3 , n , k ))+&
             ( 1080.0 * f( l1-3:l2-3 , m + 3 , n - 1 , k ))+&
             ( -108.0 * f( l1-3:l2-3 , m + 3 , n - 2 , k ))+&
             ( 8.0 * f( l1-3:l2-3 , m + 3 , n - 3 , k ))+&
             ( 1080.0 * f( l1-3:l2-3 , m + 3 , n + 1 , k ))+&
             ( -108.0 * f( l1-3:l2-3 , m + 3 , n + 2 , k ))+&
             ( 8.0 * f( l1-3:l2-3 , m + 3 , n + 3 , k ))+&
             ( 64827000.0 * f( l1+1:l2+1 , m , n , k ))+&
             ( -35721000.0 * f( l1+1:l2+1 , m , n - 1 , k ))+&
             ( 3572100.0 * f( l1+1:l2+1 , m , n - 2 , k ))+&
             ( -264600.0 * f( l1+1:l2+1 , m , n - 3 , k ))+&
             ( -35721000.0 * f( l1+1:l2+1 , m , n + 1 , k ))+&
             ( 3572100.0 * f( l1+1:l2+1 , m , n + 2 , k ))+&
             ( -264600.0 * f( l1+1:l2+1 , m , n + 3 , k ))+&
             ( -35721000.0 * f( l1+1:l2+1 , m - 1 , n , k ))+&
             ( 19683000.0 * f( l1+1:l2+1 , m - 1 , n - 1 , k ))+&
             ( -1968300.0 * f( l1+1:l2+1 , m - 1 , n - 2 , k ))+&
             ( 145800.0 * f( l1+1:l2+1 , m - 1 , n - 3 , k ))+&
             ( 19683000.0 * f( l1+1:l2+1 , m - 1 , n + 1 , k ))+&
             ( -1968300.0 * f( l1+1:l2+1 , m - 1 , n + 2 , k ))+&
             ( 145800.0 * f( l1+1:l2+1 , m - 1 , n + 3 , k ))+&
             ( 3572100.0 * f( l1+1:l2+1 , m - 2 , n , k ))+&
             ( -1968300.0 * f( l1+1:l2+1 , m - 2 , n - 1 , k ))+&
             ( 196830.0 * f( l1+1:l2+1 , m - 2 , n - 2 , k ))+&
             ( -14580.0 * f( l1+1:l2+1 , m - 2 , n - 3 , k ))+&
             ( -1968300.0 * f( l1+1:l2+1 , m - 2 , n + 1 , k ))+&
             ( 196830.0 * f( l1+1:l2+1 , m - 2 , n + 2 , k ))+&
             ( -14580.0 * f( l1+1:l2+1 , m - 2 , n + 3 , k ))+&
             ( -264600.0 * f( l1+1:l2+1 , m - 3 , n , k ))+&
             ( 145800.0 * f( l1+1:l2+1 , m - 3 , n - 1 , k ))+&
             ( -14580.0 * f( l1+1:l2+1 , m - 3 , n - 2 , k ))+&
             ( 1080.0 * f( l1+1:l2+1 , m - 3 , n - 3 , k ))+&
             ( 145800.0 * f( l1+1:l2+1 , m - 3 , n + 1 , k ))+&
             ( -14580.0 * f( l1+1:l2+1 , m - 3 , n + 2 , k ))+&
             ( 1080.0 * f( l1+1:l2+1 , m - 3 , n + 3 , k ))+&
             ( -35721000.0 * f( l1+1:l2+1 , m + 1 , n , k ))+&
             ( 19683000.0 * f( l1+1:l2+1 , m + 1 , n - 1 , k ))+&
             ( -1968300.0 * f( l1+1:l2+1 , m + 1 , n - 2 , k ))+&
             ( 145800.0 * f( l1+1:l2+1 , m + 1 , n - 3 , k ))+&
             ( 19683000.0 * f( l1+1:l2+1 , m + 1 , n + 1 , k ))+&
             ( -1968300.0 * f( l1+1:l2+1 , m + 1 , n + 2 , k ))+&
             ( 145800.0 * f( l1+1:l2+1 , m + 1 , n + 3 , k ))+&
             ( 3572100.0 * f( l1+1:l2+1 , m + 2 , n , k ))+&
             ( -1968300.0 * f( l1+1:l2+1 , m + 2 , n - 1 , k ))+&
             ( 196830.0 * f( l1+1:l2+1 , m + 2 , n - 2 , k ))+&
             ( -14580.0 * f( l1+1:l2+1 , m + 2 , n - 3 , k ))+&
             ( -1968300.0 * f( l1+1:l2+1 , m + 2 , n + 1 , k ))+&
             ( 196830.0 * f( l1+1:l2+1 , m + 2 , n + 2 , k ))+&
             ( -14580.0 * f( l1+1:l2+1 , m + 2 , n + 3 , k ))+&
             ( -264600.0 * f( l1+1:l2+1 , m + 3 , n , k ))+&
             ( 145800.0 * f( l1+1:l2+1 , m + 3 , n - 1 , k ))+&
             ( -14580.0 * f( l1+1:l2+1 , m + 3 , n - 2 , k ))+&
             ( 1080.0 * f( l1+1:l2+1 , m + 3 , n - 3 , k ))+&
             ( 145800.0 * f( l1+1:l2+1 , m + 3 , n + 1 , k ))+&
             ( -14580.0 * f( l1+1:l2+1 , m + 3 , n + 2 , k ))+&
             ( 1080.0 * f( l1+1:l2+1 , m + 3 , n + 3 , k ))+&
             ( -6482700.0 * f( l1+2:l2+2 , m , n , k ))+&
             ( 3572100.0 * f( l1+2:l2+2 , m , n - 1 , k ))+&
             ( -357210.0 * f( l1+2:l2+2 , m , n - 2 , k ))+&
             ( 26460.0 * f( l1+2:l2+2 , m , n - 3 , k ))+&
             ( 3572100.0 * f( l1+2:l2+2 , m , n + 1 , k ))+&
             ( -357210.0 * f( l1+2:l2+2 , m , n + 2 , k ))+&
             ( 26460.0 * f( l1+2:l2+2 , m , n + 3 , k ))+&
             ( 3572100.0 * f( l1+2:l2+2 , m - 1 , n , k ))+&
             ( -1968300.0 * f( l1+2:l2+2 , m - 1 , n - 1 , k ))+&
             ( 196830.0 * f( l1+2:l2+2 , m - 1 , n - 2 , k ))+&
             ( -14580.0 * f( l1+2:l2+2 , m - 1 , n - 3 , k ))+&
             ( -1968300.0 * f( l1+2:l2+2 , m - 1 , n + 1 , k ))+&
             ( 196830.0 * f( l1+2:l2+2 , m - 1 , n + 2 , k ))+&
             ( -14580.0 * f( l1+2:l2+2 , m - 1 , n + 3 , k ))+&
             ( -357210.0 * f( l1+2:l2+2 , m - 2 , n , k ))+&
             ( 196830.0 * f( l1+2:l2+2 , m - 2 , n - 1 , k ))+&
             ( -19683.0 * f( l1+2:l2+2 , m - 2 , n - 2 , k ))+&
             ( 1458.0 * f( l1+2:l2+2 , m - 2 , n - 3 , k ))+&
             ( 196830.0 * f( l1+2:l2+2 , m - 2 , n + 1 , k ))+&
             ( -19683.0 * f( l1+2:l2+2 , m - 2 , n + 2 , k ))+&
             ( 1458.0 * f( l1+2:l2+2 , m - 2 , n + 3 , k ))+&
             ( 26460.0 * f( l1+2:l2+2 , m - 3 , n , k ))+&
             ( -14580.0 * f( l1+2:l2+2 , m - 3 , n - 1 , k ))+&
             ( 1458.0 * f( l1+2:l2+2 , m - 3 , n - 2 , k ))+&
             ( -108.0 * f( l1+2:l2+2 , m - 3 , n - 3 , k ))+&
             ( -14580.0 * f( l1+2:l2+2 , m - 3 , n + 1 , k ))+&
             ( 1458.0 * f( l1+2:l2+2 , m - 3 , n + 2 , k ))+&
             ( -108.0 * f( l1+2:l2+2 , m - 3 , n + 3 , k ))+&
             ( 3572100.0 * f( l1+2:l2+2 , m + 1 , n , k ))+&
             ( -1968300.0 * f( l1+2:l2+2 , m + 1 , n - 1 , k ))+&
             ( 196830.0 * f( l1+2:l2+2 , m + 1 , n - 2 , k ))+&
             ( -14580.0 * f( l1+2:l2+2 , m + 1 , n - 3 , k ))+&
             ( -1968300.0 * f( l1+2:l2+2 , m + 1 , n + 1 , k ))+&
             ( 196830.0 * f( l1+2:l2+2 , m + 1 , n + 2 , k ))+&
             ( -14580.0 * f( l1+2:l2+2 , m + 1 , n + 3 , k ))+&
             ( -357210.0 * f( l1+2:l2+2 , m + 2 , n , k ))+&
             ( 196830.0 * f( l1+2:l2+2 , m + 2 , n - 1 , k ))+&
             ( -19683.0 * f( l1+2:l2+2 , m + 2 , n - 2 , k ))+&
             ( 1458.0 * f( l1+2:l2+2 , m + 2 , n - 3 , k ))+&
             ( 196830.0 * f( l1+2:l2+2 , m + 2 , n + 1 , k ))+&
             ( -19683.0 * f( l1+2:l2+2 , m + 2 , n + 2 , k ))+&
             ( 1458.0 * f( l1+2:l2+2 , m + 2 , n + 3 , k ))+&
             ( 26460.0 * f( l1+2:l2+2 , m + 3 , n , k ))+&
             ( -14580.0 * f( l1+2:l2+2 , m + 3 , n - 1 , k ))+&
             ( 1458.0 * f( l1+2:l2+2 , m + 3 , n - 2 , k ))+&
             ( -108.0 * f( l1+2:l2+2 , m + 3 , n - 3 , k ))+&
             ( -14580.0 * f( l1+2:l2+2 , m + 3 , n + 1 , k ))+&
             ( 1458.0 * f( l1+2:l2+2 , m + 3 , n + 2 , k ))+&
             ( -108.0 * f( l1+2:l2+2 , m + 3 , n + 3 , k ))+&
             ( 480200.0 * f( l1+3:l2+3 , m , n , k ))+&
             ( -264600.0 * f( l1+3:l2+3 , m , n - 1 , k ))+&
             ( 26460.0 * f( l1+3:l2+3 , m , n - 2 , k ))+&
             ( -1960.0 * f( l1+3:l2+3 , m , n - 3 , k ))+&
             ( -264600.0 * f( l1+3:l2+3 , m , n + 1 , k ))+&
             ( 26460.0 * f( l1+3:l2+3 , m , n + 2 , k ))+&
             ( -1960.0 * f( l1+3:l2+3 , m , n + 3 , k ))+&
             ( -264600.0 * f( l1+3:l2+3 , m - 1 , n , k ))+&
             ( 145800.0 * f( l1+3:l2+3 , m - 1 , n - 1 , k ))+&
             ( -14580.0 * f( l1+3:l2+3 , m - 1 , n - 2 , k ))+&
             ( 1080.0 * f( l1+3:l2+3 , m - 1 , n - 3 , k ))+&
             ( 145800.0 * f( l1+3:l2+3 , m - 1 , n + 1 , k ))+&
             ( -14580.0 * f( l1+3:l2+3 , m - 1 , n + 2 , k ))+&
             ( 1080.0 * f( l1+3:l2+3 , m - 1 , n + 3 , k ))+&
             ( 26460.0 * f( l1+3:l2+3 , m - 2 , n , k ))+&
             ( -14580.0 * f( l1+3:l2+3 , m - 2 , n - 1 , k ))+&
             ( 1458.0 * f( l1+3:l2+3 , m - 2 , n - 2 , k ))+&
             ( -108.0 * f( l1+3:l2+3 , m - 2 , n - 3 , k ))+&
             ( -14580.0 * f( l1+3:l2+3 , m - 2 , n + 1 , k ))+&
             ( 1458.0 * f( l1+3:l2+3 , m - 2 , n + 2 , k ))+&
             ( -108.0 * f( l1+3:l2+3 , m - 2 , n + 3 , k ))+&
             ( -1960.0 * f( l1+3:l2+3 , m - 3 , n , k ))+&
             ( 1080.0 * f( l1+3:l2+3 , m - 3 , n - 1 , k ))+&
             ( -108.0 * f( l1+3:l2+3 , m - 3 , n - 2 , k ))+&
             ( 8.0 * f( l1+3:l2+3 , m - 3 , n - 3 , k ))+&
             ( 1080.0 * f( l1+3:l2+3 , m - 3 , n + 1 , k ))+&
             ( -108.0 * f( l1+3:l2+3 , m - 3 , n + 2 , k ))+&
             ( 8.0 * f( l1+3:l2+3 , m - 3 , n + 3 , k ))+&
             ( -264600.0 * f( l1+3:l2+3 , m + 1 , n , k ))+&
             ( 145800.0 * f( l1+3:l2+3 , m + 1 , n - 1 , k ))+&
             ( -14580.0 * f( l1+3:l2+3 , m + 1 , n - 2 , k ))+&
             ( 1080.0 * f( l1+3:l2+3 , m + 1 , n - 3 , k ))+&
             ( 145800.0 * f( l1+3:l2+3 , m + 1 , n + 1 , k ))+&
             ( -14580.0 * f( l1+3:l2+3 , m + 1 , n + 2 , k ))+&
             ( 1080.0 * f( l1+3:l2+3 , m + 1 , n + 3 , k ))+&
             ( 26460.0 * f( l1+3:l2+3 , m + 2 , n , k ))+&
             ( -14580.0 * f( l1+3:l2+3 , m + 2 , n - 1 , k ))+&
             ( 1458.0 * f( l1+3:l2+3 , m + 2 , n - 2 , k ))+&
             ( -108.0 * f( l1+3:l2+3 , m + 2 , n - 3 , k ))+&
             ( -14580.0 * f( l1+3:l2+3 , m + 2 , n + 1 , k ))+&
             ( 1458.0 * f( l1+3:l2+3 , m + 2 , n + 2 , k ))+&
             ( -108.0 * f( l1+3:l2+3 , m + 2 , n + 3 , k ))+&
             ( -1960.0 * f( l1+3:l2+3 , m + 3 , n , k ))+&
             ( 1080.0 * f( l1+3:l2+3 , m + 3 , n - 1 , k ))+&
             ( -108.0 * f( l1+3:l2+3 , m + 3 , n - 2 , k ))+&
             ( 8.0 * f( l1+3:l2+3 , m + 3 , n - 3 , k ))+&
             ( 1080.0 * f( l1+3:l2+3 , m + 3 , n + 1 , k ))+&
             ( -108.0 * f( l1+3:l2+3 , m + 3 , n + 2 , k ))+&
             ( 8.0 * f( l1+3:l2+3 , m + 3 , n + 3 , k ))&
             )
      else
        df=0.
      endif
      
    endsubroutine der2i2j2k
!***********************************************************************
    subroutine der3i3j(f,k,df,i,j)
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx), intent(out) :: df
      real, dimension (nx) :: fac
      integer, intent(in) :: k,i,j
!
      if (headtt) then
        if (lspherical_coords.or.lcylindrical_coords) &
             call warning('der3i3j','MISSING CURVATURE TERMS for non-cartesian coordinates')
      endif
!
      if (i==j) then
        if (i==1.and.nxgrid==1 .or. i==2.and.nygrid==1 .or. i==3.and.nzgrid==1) then
          df=0.
          if (ip<=5) print*, 'der3i3j: Degenerate case in '//coornames(i)//'-direction'
        else
          call der6(f,k,df,j)
        endif
      elseif (i+j==3) then
        xy: if (nxgrid/=1.and.nygrid/=1) then
          fac= 1./8.0**2 * dx_1(l1:l2)**3 * dy_1(m)**3
          df = fac*(&
               ( 169.0 * f( l1-1:l2-1 , m - 1 , n , k ))+&
               (-104.0 * f( l1-2:l2-2 , m - 1 , n , k ))+&
               ( 13.0 * f( l1-3:l2-3 , m - 1 , n , k ))+&
               (-169.0 * f( l1+1:l2+1 , m - 1 , n , k ))+&
               ( 104.0 * f( l1+2:l2+2 , m - 1 , n , k ))+&
               (-13.0 * f( l1+3:l2+3 , m - 1 , n , k ))+&
               (-104.0 * f( l1-1:l2-1 , m - 2 , n , k ))+&
               ( 64.0 * f( l1-2:l2-2 , m - 2 , n , k ))+&
               (-8.0 * f( l1-3:l2-3 , m - 2 , n , k ))+&
               ( 104.0 * f( l1+1:l2+1 , m - 2 , n , k ))+&
               (-64.0 * f( l1+2:l2+2 , m - 2 , n , k ))+&
               ( 8.0 * f( l1+3:l2+3 , m - 2 , n , k ))+&
               ( 13.0 * f( l1-1:l2-1 , m - 3 , n , k ))+&
               (-8.0 * f( l1-2:l2-2 , m - 3 , n , k ))+&
               ( 1.0 * f( l1-3:l2-3 , m - 3 , n , k ))+&
               (-13.0 * f( l1+1:l2+1 , m - 3 , n , k ))+&
               ( 8.0 * f( l1+2:l2+2 , m - 3 , n , k ))+&
               (-1.0 * f( l1+3:l2+3 , m - 3 , n , k ))+&
               (-169.0 * f( l1-1:l2-1 , m + 1 , n , k ))+&
               ( 104.0 * f( l1-2:l2-2 , m + 1 , n , k ))+&
               (-13.0 * f( l1-3:l2-3 , m + 1 , n , k ))+&
               ( 169.0 * f( l1+1:l2+1 , m + 1 , n , k ))+&
               (-104.0 * f( l1+2:l2+2 , m + 1 , n , k ))+&
               ( 13.0 * f( l1+3:l2+3 , m + 1 , n , k ))+&
               ( 104.0 * f( l1-1:l2-1 , m + 2 , n , k ))+&
               (-64.0 * f( l1-2:l2-2 , m + 2 , n , k ))+&
               ( 8.0 * f( l1-3:l2-3 , m + 2 , n , k ))+&
               (-104.0 * f( l1+1:l2+1 , m + 2 , n , k ))+&
               ( 64.0 * f( l1+2:l2+2 , m + 2 , n , k ))+&
               (-8.0 * f( l1+3:l2+3 , m + 2 , n , k ))+&
               (-13.0 * f( l1-1:l2-1 , m + 3 , n , k ))+&
               ( 8.0 * f( l1-2:l2-2 , m + 3 , n , k ))+&
               (-1.0 * f( l1-3:l2-3 , m + 3 , n , k ))+&
               ( 13.0 * f( l1+1:l2+1 , m + 3 , n , k ))+&
               (-8.0 * f( l1+2:l2+2 , m + 3 , n , k ))+&
               ( 1.0 * f( l1+3:l2+3 , m + 3 , n , k ))&
               )
        else
          df=0.0
        endif xy
!
      elseif (i+j==4) then
        xz: if (nxgrid/=1.and.nzgrid/=1) then
          fac= 1./8.0**2 * dx_1(l1:l2)**3 * dz_1(n)**3
          df = fac*(& 
               ( 169.0 * f( l1-1:l2-1 , m , n - 1 , k ))+&
               (-104.0 * f( l1-2:l2-2 , m , n - 1 , k ))+&
               ( 13.0 * f( l1-3:l2-3 , m , n - 1 , k ))+&
               (-169.0 * f( l1+1:l2+1 , m , n - 1 , k ))+&
               ( 104.0 * f( l1+2:l2+2 , m , n - 1 , k ))+&
               (-13.0 * f( l1+3:l2+3 , m , n - 1 , k ))+&
               (-104.0 * f( l1-1:l2-1 , m , n - 2 , k ))+&
               ( 64.0 * f( l1-2:l2-2 , m , n - 2 , k ))+&
               (-8.0 * f( l1-3:l2-3 , m , n - 2 , k ))+&
               ( 104.0 * f( l1+1:l2+1 , m , n - 2 , k ))+&
               (-64.0 * f( l1+2:l2+2 , m , n - 2 , k ))+&
               ( 8.0 * f( l1+3:l2+3 , m , n - 2 , k ))+&
               ( 13.0 * f( l1-1:l2-1 , m , n - 3 , k ))+&
               (-8.0 * f( l1-2:l2-2 , m , n - 3 , k ))+&
               ( 1.0 * f( l1-3:l2-3 , m , n - 3 , k ))+&
               (-13.0 * f( l1+1:l2+1 , m , n - 3 , k ))+&
               ( 8.0 * f( l1+2:l2+2 , m , n - 3 , k ))+&
               (-1.0 * f( l1+3:l2+3 , m , n - 3 , k ))+&
               (-169.0 * f( l1-1:l2-1 , m , n + 1 , k ))+&
               ( 104.0 * f( l1-2:l2-2 , m , n + 1 , k ))+&
               (-13.0 * f( l1-3:l2-3 , m , n + 1 , k ))+&
               ( 169.0 * f( l1+1:l2+1 , m , n + 1 , k ))+&
               (-104.0 * f( l1+2:l2+2 , m , n + 1 , k ))+&
               ( 13.0 * f( l1+3:l2+3 , m , n + 1 , k ))+&
               ( 104.0 * f( l1-1:l2-1 , m , n + 2 , k ))+&
               (-64.0 * f( l1-2:l2-2 , m , n + 2 , k ))+&
               ( 8.0 * f( l1-3:l2-3 , m , n + 2 , k ))+&
               (-104.0 * f( l1+1:l2+1 , m , n + 2 , k ))+&
               ( 64.0 * f( l1+2:l2+2 , m , n + 2 , k ))+&
               (-8.0 * f( l1+3:l2+3 , m , n + 2 , k ))+&
               (-13.0 * f( l1-1:l2-1 , m , n + 3 , k ))+&
               ( 8.0 * f( l1-2:l2-2 , m , n + 3 , k ))+&
               (-1.0 * f( l1-3:l2-3 , m , n + 3 , k ))+&
               ( 13.0 * f( l1+1:l2+1 , m , n + 3 , k ))+&
               (-8.0 * f( l1+2:l2+2 , m , n + 3 , k ))+&
               ( 1.0 * f( l1+3:l2+3 , m , n + 3 , k ))&
               )
        else
          df=0.0
        endif xz
      elseif (i+j==5) then
        yz: if (nygrid/=1.and.nzgrid/=1) then
          fac= 1./8.0**2 * dy_1(m)**3 * dz_1(n)**3
          df = fac*(& 
               ( 169.0 * f( l1:l2 , m - 1 , n - 1 , k ))+&
               (-104.0 * f( l1:l2 , m - 2 , n - 1 , k ))+&
               ( 13.0 * f( l1:l2 , m - 3 , n - 1 , k ))+&
               (-169.0 * f( l1:l2 , m + 1 , n - 1 , k ))+&
               ( 104.0 * f( l1:l2 , m + 2 , n - 1 , k ))+&
               (-13.0 * f( l1:l2 , m + 3 , n - 1 , k ))+&
               (-104.0 * f( l1:l2 , m - 1 , n - 2 , k ))+&
               ( 64.0 * f( l1:l2 , m - 2 , n - 2 , k ))+&
               (-8.0 * f( l1:l2 , m - 3 , n - 2 , k ))+&
               ( 104.0 * f( l1:l2 , m + 1 , n - 2 , k ))+&
               (-64.0 * f( l1:l2 , m + 2 , n - 2 , k ))+&
               ( 8.0 * f( l1:l2 , m + 3 , n - 2 , k ))+&
               ( 13.0 * f( l1:l2 , m - 1 , n - 3 , k ))+&
               (-8.0 * f( l1:l2 , m - 2 , n - 3 , k ))+&
               ( 1.0 * f( l1:l2 , m - 3 , n - 3 , k ))+&
               (-13.0 * f( l1:l2 , m + 1 , n - 3 , k ))+&
               ( 8.0 * f( l1:l2 , m + 2 , n - 3 , k ))+&
               (-1.0 * f( l1:l2 , m + 3 , n - 3 , k ))+&
               (-169.0 * f( l1:l2 , m - 1 , n + 1 , k ))+&
               ( 104.0 * f( l1:l2 , m - 2 , n + 1 , k ))+&
               (-13.0 * f( l1:l2 , m - 3 , n + 1 , k ))+&
               ( 169.0 * f( l1:l2 , m + 1 , n + 1 , k ))+&
               (-104.0 * f( l1:l2 , m + 2 , n + 1 , k ))+&
               ( 13.0 * f( l1:l2 , m + 3 , n + 1 , k ))+&
               ( 104.0 * f( l1:l2 , m - 1 , n + 2 , k ))+&
               (-64.0 * f( l1:l2 , m - 2 , n + 2 , k ))+&
               ( 8.0 * f( l1:l2 , m - 3 , n + 2 , k ))+&
               (-104.0 * f( l1:l2 , m + 1 , n + 2 , k ))+&
               ( 64.0 * f( l1:l2 , m + 2 , n + 2 , k ))+&
               (-8.0 * f( l1:l2 , m + 3 , n + 2 , k ))+&
               (-13.0 * f( l1:l2 , m - 1 , n + 3 , k ))+&
               ( 8.0 * f( l1:l2 , m - 2 , n + 3 , k ))+&
               (-1.0 * f( l1:l2 , m - 3 , n + 3 , k ))+&
               ( 13.0 * f( l1:l2 , m + 1 , n + 3 , k ))+&
               (-8.0 * f( l1:l2 , m + 2 , n + 3 , k ))+&
               ( 1.0 * f( l1:l2 , m + 3 , n + 3 , k ))&
               )
        else
          df=0.0
        endif yz
!
      endif
!               
    endsubroutine der3i3j
!***********************************************************************          
    subroutine der3i2j1k(f,ik,df,i,j,k)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx), intent(out) :: df
      real, dimension (nx) :: fac
      integer, intent(in) :: ik,i,j,k
!
      if (headtt) then
        if (lspherical_coords.or.lcylindrical_coords) &
             call warning('der3i2j1k','MISSING CURVATURE TERMS for non-cartesian coordinates')
      endif
!
      if (i==j.and.j==k) then
        if (i==1.and.nxgrid==1 .or. i==2.and.nygrid==1 .or. i==3.and.nzgrid==1) then
          df=0.
          if (ip<=5) print*, 'der3i2j1k: Degenerate case in '//coornames(i)//'-direction'
        else
          call der6(f,ik,df,j)
        endif
!  MR: cases i=j/=k etc. missing
      elseif (i==1.and.j==2.and.k==3) then
        xyz: if (nxgrid/=1.and.nygrid/=1.and.nzgrid/=1) then
          fac= (1./8.0*dx_1(l1:l2)**3) * (1./180.0*dy_1(m)**2) * (1/60.0*dz_1(n))
          df = fac*(&
               ( 286650.0 * f( l1-1:l2-1 , m , n - 1 , ik ))+&
               ( -176400.0 * f( l1-2:l2-2 , m , n - 1 , ik ))+&
               ( 22050.0 * f( l1-3:l2-3 , m , n - 1 , ik ))+&
               ( -286650.0 * f( l1+1:l2+1 , m , n - 1 , ik ))+&
               ( 176400.0 * f( l1+2:l2+2 , m , n - 1 , ik ))+&
               ( -22050.0 * f( l1+3:l2+3 , m , n - 1 , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m - 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1-2:l2-2 , m - 1 , n - 1 , ik ))+&
               ( -12150.0 * f( l1-3:l2-3 , m - 1 , n - 1 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m - 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1+2:l2+2 , m - 1 , n - 1 , ik ))+&
               ( 12150.0 * f( l1+3:l2+3 , m - 1 , n - 1 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 2 , n - 1 , ik ))+&
               ( -9720.0 * f( l1-2:l2-2 , m - 2 , n - 1 , ik ))+&
               ( 1215.0 * f( l1-3:l2-3 , m - 2 , n - 1 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 2 , n - 1 , ik ))+&
               ( 9720.0 * f( l1+2:l2+2 , m - 2 , n - 1 , ik ))+&
               ( -1215.0 * f( l1+3:l2+3 , m - 2 , n - 1 , ik ))+&
               ( -1170.0 * f( l1-1:l2-1 , m - 3 , n - 1 , ik ))+&
               ( 720.0 * f( l1-2:l2-2 , m - 3 , n - 1 , ik ))+&
               ( -90.0 * f( l1-3:l2-3 , m - 3 , n - 1 , ik ))+&
               ( 1170.0 * f( l1+1:l2+1 , m - 3 , n - 1 , ik ))+&
               ( -720.0 * f( l1+2:l2+2 , m - 3 , n - 1 , ik ))+&
               ( 90.0 * f( l1+3:l2+3 , m - 3 , n - 1 , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m + 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1-2:l2-2 , m + 1 , n - 1 , ik ))+&
               ( -12150.0 * f( l1-3:l2-3 , m + 1 , n - 1 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m + 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1+2:l2+2 , m + 1 , n - 1 , ik ))+&
               ( 12150.0 * f( l1+3:l2+3 , m + 1 , n - 1 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m + 2 , n - 1 , ik ))+&
               ( -9720.0 * f( l1-2:l2-2 , m + 2 , n - 1 , ik ))+&
               ( 1215.0 * f( l1-3:l2-3 , m + 2 , n - 1 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m + 2 , n - 1 , ik ))+&
               ( 9720.0 * f( l1+2:l2+2 , m + 2 , n - 1 , ik ))+&
               ( -1215.0 * f( l1+3:l2+3 , m + 2 , n - 1 , ik ))+&
               ( -1170.0 * f( l1-1:l2-1 , m + 3 , n - 1 , ik ))+&
               ( 720.0 * f( l1-2:l2-2 , m + 3 , n - 1 , ik ))+&
               ( -90.0 * f( l1-3:l2-3 , m + 3 , n - 1 , ik ))+&
               ( 1170.0 * f( l1+1:l2+1 , m + 3 , n - 1 , ik ))+&
               ( -720.0 * f( l1+2:l2+2 , m + 3 , n - 1 , ik ))+&
               ( 90.0 * f( l1+3:l2+3 , m + 3 , n - 1 , ik ))+&
               ( -57330.0 * f( l1-1:l2-1 , m , n - 2 , ik ))+&
               ( 35280.0 * f( l1-2:l2-2 , m , n - 2 , ik ))+&
               ( -4410.0 * f( l1-3:l2-3 , m , n - 2 , ik ))+&
               ( 57330.0 * f( l1+1:l2+1 , m , n - 2 , ik ))+&
               ( -35280.0 * f( l1+2:l2+2 , m , n - 2 , ik ))+&
               ( 4410.0 * f( l1+3:l2+3 , m , n - 2 , ik ))+&
               ( 31590.0 * f( l1-1:l2-1 , m - 1 , n - 2 , ik ))+&
               ( -19440.0 * f( l1-2:l2-2 , m - 1 , n - 2 , ik ))+&
               ( 2430.0 * f( l1-3:l2-3 , m - 1 , n - 2 , ik ))+&
               ( -31590.0 * f( l1+1:l2+1 , m - 1 , n - 2 , ik ))+&
               ( 19440.0 * f( l1+2:l2+2 , m - 1 , n - 2 , ik ))+&
               ( -2430.0 * f( l1+3:l2+3 , m - 1 , n - 2 , ik ))+&
               ( -3159.0 * f( l1-1:l2-1 , m - 2 , n - 2 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m - 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1-3:l2-3 , m - 2 , n - 2 , ik ))+&
               ( 3159.0 * f( l1+1:l2+1 , m - 2 , n - 2 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m - 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1+3:l2+3 , m - 2 , n - 2 , ik ))+&
               ( 234.0 * f( l1-1:l2-1 , m - 3 , n - 2 , ik ))+&
               ( -144.0 * f( l1-2:l2-2 , m - 3 , n - 2 , ik ))+&
               ( 18.0 * f( l1-3:l2-3 , m - 3 , n - 2 , ik ))+&
               ( -234.0 * f( l1+1:l2+1 , m - 3 , n - 2 , ik ))+&
               ( 144.0 * f( l1+2:l2+2 , m - 3 , n - 2 , ik ))+&
               ( -18.0 * f( l1+3:l2+3 , m - 3 , n - 2 , ik ))+&
               ( 31590.0 * f( l1-1:l2-1 , m + 1 , n - 2 , ik ))+&
               ( -19440.0 * f( l1-2:l2-2 , m + 1 , n - 2 , ik ))+&
               ( 2430.0 * f( l1-3:l2-3 , m + 1 , n - 2 , ik ))+&
               ( -31590.0 * f( l1+1:l2+1 , m + 1 , n - 2 , ik ))+&
               ( 19440.0 * f( l1+2:l2+2 , m + 1 , n - 2 , ik ))+&
               ( -2430.0 * f( l1+3:l2+3 , m + 1 , n - 2 , ik ))+&
               ( -3159.0 * f( l1-1:l2-1 , m + 2 , n - 2 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m + 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1-3:l2-3 , m + 2 , n - 2 , ik ))+&
               ( 3159.0 * f( l1+1:l2+1 , m + 2 , n - 2 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m + 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1+3:l2+3 , m + 2 , n - 2 , ik ))+&
               ( 234.0 * f( l1-1:l2-1 , m + 3 , n - 2 , ik ))+&
               ( -144.0 * f( l1-2:l2-2 , m + 3 , n - 2 , ik ))+&
               ( 18.0 * f( l1-3:l2-3 , m + 3 , n - 2 , ik ))+&
               ( -234.0 * f( l1+1:l2+1 , m + 3 , n - 2 , ik ))+&
               ( 144.0 * f( l1+2:l2+2 , m + 3 , n - 2 , ik ))+&
               ( -18.0 * f( l1+3:l2+3 , m + 3 , n - 2 , ik ))+&
               ( 6370.0 * f( l1-1:l2-1 , m , n - 3 , ik ))+&
               ( -3920.0 * f( l1-2:l2-2 , m , n - 3 , ik ))+&
               ( 490.0 * f( l1-3:l2-3 , m , n - 3 , ik ))+&
               ( -6370.0 * f( l1+1:l2+1 , m , n - 3 , ik ))+&
               ( 3920.0 * f( l1+2:l2+2 , m , n - 3 , ik ))+&
               ( -490.0 * f( l1+3:l2+3 , m , n - 3 , ik ))+&
               ( -3510.0 * f( l1-1:l2-1 , m - 1 , n - 3 , ik ))+&
               ( 2160.0 * f( l1-2:l2-2 , m - 1 , n - 3 , ik ))+&
               ( -270.0 * f( l1-3:l2-3 , m - 1 , n - 3 , ik ))+&
               ( 3510.0 * f( l1+1:l2+1 , m - 1 , n - 3 , ik ))+&
               ( -2160.0 * f( l1+2:l2+2 , m - 1 , n - 3 , ik ))+&
               ( 270.0 * f( l1+3:l2+3 , m - 1 , n - 3 , ik ))+&
               ( 351.0 * f( l1-1:l2-1 , m - 2 , n - 3 , ik ))+&
               ( -216.0 * f( l1-2:l2-2 , m - 2 , n - 3 , ik ))+&
               ( 27.0 * f( l1-3:l2-3 , m - 2 , n - 3 , ik ))+&
               ( -351.0 * f( l1+1:l2+1 , m - 2 , n - 3 , ik ))+&
               ( 216.0 * f( l1+2:l2+2 , m - 2 , n - 3 , ik ))+&
               ( -27.0 * f( l1+3:l2+3 , m - 2 , n - 3 , ik ))+&
               ( -26.0 * f( l1-1:l2-1 , m - 3 , n - 3 , ik ))+&
               ( 16.0 * f( l1-2:l2-2 , m - 3 , n - 3 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m - 3 , n - 3 , ik ))+&
               ( 26.0 * f( l1+1:l2+1 , m - 3 , n - 3 , ik ))+&
               ( -16.0 * f( l1+2:l2+2 , m - 3 , n - 3 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m - 3 , n - 3 , ik ))+&
               ( -3510.0 * f( l1-1:l2-1 , m + 1 , n - 3 , ik ))+&
               ( 2160.0 * f( l1-2:l2-2 , m + 1 , n - 3 , ik ))+&
               ( -270.0 * f( l1-3:l2-3 , m + 1 , n - 3 , ik ))+&
               ( 3510.0 * f( l1+1:l2+1 , m + 1 , n - 3 , ik ))+&
               ( -2160.0 * f( l1+2:l2+2 , m + 1 , n - 3 , ik ))+&
               ( 270.0 * f( l1+3:l2+3 , m + 1 , n - 3 , ik ))+&
               ( 351.0 * f( l1-1:l2-1 , m + 2 , n - 3 , ik ))+&
               ( -216.0 * f( l1-2:l2-2 , m + 2 , n - 3 , ik ))+&
               ( 27.0 * f( l1-3:l2-3 , m + 2 , n - 3 , ik ))+&
               ( -351.0 * f( l1+1:l2+1 , m + 2 , n - 3 , ik ))+&
               ( 216.0 * f( l1+2:l2+2 , m + 2 , n - 3 , ik ))+&
               ( -27.0 * f( l1+3:l2+3 , m + 2 , n - 3 , ik ))+&
               ( -26.0 * f( l1-1:l2-1 , m + 3 , n - 3 , ik ))+&
               ( 16.0 * f( l1-2:l2-2 , m + 3 , n - 3 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m + 3 , n - 3 , ik ))+&
               ( 26.0 * f( l1+1:l2+1 , m + 3 , n - 3 , ik ))+&
               ( -16.0 * f( l1+2:l2+2 , m + 3 , n - 3 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m + 3 , n - 3 , ik ))+&
               ( -286650.0 * f( l1-1:l2-1 , m , n + 1 , ik ))+&
               ( 176400.0 * f( l1-2:l2-2 , m , n + 1 , ik ))+&
               ( -22050.0 * f( l1-3:l2-3 , m , n + 1 , ik ))+&
               ( 286650.0 * f( l1+1:l2+1 , m , n + 1 , ik ))+&
               ( -176400.0 * f( l1+2:l2+2 , m , n + 1 , ik ))+&
               ( 22050.0 * f( l1+3:l2+3 , m , n + 1 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m - 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1-2:l2-2 , m - 1 , n + 1 , ik ))+&
               ( 12150.0 * f( l1-3:l2-3 , m - 1 , n + 1 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m - 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1+2:l2+2 , m - 1 , n + 1 , ik ))+&
               ( -12150.0 * f( l1+3:l2+3 , m - 1 , n + 1 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m - 2 , n + 1 , ik ))+&
               ( 9720.0 * f( l1-2:l2-2 , m - 2 , n + 1 , ik ))+&
               ( -1215.0 * f( l1-3:l2-3 , m - 2 , n + 1 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m - 2 , n + 1 , ik ))+&
               ( -9720.0 * f( l1+2:l2+2 , m - 2 , n + 1 , ik ))+&
               ( 1215.0 * f( l1+3:l2+3 , m - 2 , n + 1 , ik ))+&
               ( 1170.0 * f( l1-1:l2-1 , m - 3 , n + 1 , ik ))+&
               ( -720.0 * f( l1-2:l2-2 , m - 3 , n + 1 , ik ))+&
               ( 90.0 * f( l1-3:l2-3 , m - 3 , n + 1 , ik ))+&
               ( -1170.0 * f( l1+1:l2+1 , m - 3 , n + 1 , ik ))+&
               ( 720.0 * f( l1+2:l2+2 , m - 3 , n + 1 , ik ))+&
               ( -90.0 * f( l1+3:l2+3 , m - 3 , n + 1 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m + 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1-2:l2-2 , m + 1 , n + 1 , ik ))+&
               ( 12150.0 * f( l1-3:l2-3 , m + 1 , n + 1 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m + 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1+2:l2+2 , m + 1 , n + 1 , ik ))+&
               ( -12150.0 * f( l1+3:l2+3 , m + 1 , n + 1 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 2 , n + 1 , ik ))+&
               ( 9720.0 * f( l1-2:l2-2 , m + 2 , n + 1 , ik ))+&
               ( -1215.0 * f( l1-3:l2-3 , m + 2 , n + 1 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 2 , n + 1 , ik ))+&
               ( -9720.0 * f( l1+2:l2+2 , m + 2 , n + 1 , ik ))+&
               ( 1215.0 * f( l1+3:l2+3 , m + 2 , n + 1 , ik ))+&
               ( 1170.0 * f( l1-1:l2-1 , m + 3 , n + 1 , ik ))+&
               ( -720.0 * f( l1-2:l2-2 , m + 3 , n + 1 , ik ))+&
               ( 90.0 * f( l1-3:l2-3 , m + 3 , n + 1 , ik ))+&
               ( -1170.0 * f( l1+1:l2+1 , m + 3 , n + 1 , ik ))+&
               ( 720.0 * f( l1+2:l2+2 , m + 3 , n + 1 , ik ))+&
               ( -90.0 * f( l1+3:l2+3 , m + 3 , n + 1 , ik ))+&
               ( 57330.0 * f( l1-1:l2-1 , m , n + 2 , ik ))+&
               ( -35280.0 * f( l1-2:l2-2 , m , n + 2 , ik ))+&
               ( 4410.0 * f( l1-3:l2-3 , m , n + 2 , ik ))+&
               ( -57330.0 * f( l1+1:l2+1 , m , n + 2 , ik ))+&
               ( 35280.0 * f( l1+2:l2+2 , m , n + 2 , ik ))+&
               ( -4410.0 * f( l1+3:l2+3 , m , n + 2 , ik ))+&
               ( -31590.0 * f( l1-1:l2-1 , m - 1 , n + 2 , ik ))+&
               ( 19440.0 * f( l1-2:l2-2 , m - 1 , n + 2 , ik ))+&
               ( -2430.0 * f( l1-3:l2-3 , m - 1 , n + 2 , ik ))+&
               ( 31590.0 * f( l1+1:l2+1 , m - 1 , n + 2 , ik ))+&
               ( -19440.0 * f( l1+2:l2+2 , m - 1 , n + 2 , ik ))+&
               ( 2430.0 * f( l1+3:l2+3 , m - 1 , n + 2 , ik ))+&
               ( 3159.0 * f( l1-1:l2-1 , m - 2 , n + 2 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m - 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1-3:l2-3 , m - 2 , n + 2 , ik ))+&
               ( -3159.0 * f( l1+1:l2+1 , m - 2 , n + 2 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m - 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1+3:l2+3 , m - 2 , n + 2 , ik ))+&
               ( -234.0 * f( l1-1:l2-1 , m - 3 , n + 2 , ik ))+&
               ( 144.0 * f( l1-2:l2-2 , m - 3 , n + 2 , ik ))+&
               ( -18.0 * f( l1-3:l2-3 , m - 3 , n + 2 , ik ))+&
               ( 234.0 * f( l1+1:l2+1 , m - 3 , n + 2 , ik ))+&
               ( -144.0 * f( l1+2:l2+2 , m - 3 , n + 2 , ik ))+&
               ( 18.0 * f( l1+3:l2+3 , m - 3 , n + 2 , ik ))+&
               ( -31590.0 * f( l1-1:l2-1 , m + 1 , n + 2 , ik ))+&
               ( 19440.0 * f( l1-2:l2-2 , m + 1 , n + 2 , ik ))+&
               ( -2430.0 * f( l1-3:l2-3 , m + 1 , n + 2 , ik ))+&
               ( 31590.0 * f( l1+1:l2+1 , m + 1 , n + 2 , ik ))+&
               ( -19440.0 * f( l1+2:l2+2 , m + 1 , n + 2 , ik ))+&
               ( 2430.0 * f( l1+3:l2+3 , m + 1 , n + 2 , ik ))+&
               ( 3159.0 * f( l1-1:l2-1 , m + 2 , n + 2 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m + 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1-3:l2-3 , m + 2 , n + 2 , ik ))+&
               ( -3159.0 * f( l1+1:l2+1 , m + 2 , n + 2 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m + 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1+3:l2+3 , m + 2 , n + 2 , ik ))+&
               ( -234.0 * f( l1-1:l2-1 , m + 3 , n + 2 , ik ))+&
               ( 144.0 * f( l1-2:l2-2 , m + 3 , n + 2 , ik ))+&
               ( -18.0 * f( l1-3:l2-3 , m + 3 , n + 2 , ik ))+&
               ( 234.0 * f( l1+1:l2+1 , m + 3 , n + 2 , ik ))+&
               ( -144.0 * f( l1+2:l2+2 , m + 3 , n + 2 , ik ))+&
               ( 18.0 * f( l1+3:l2+3 , m + 3 , n + 2 , ik ))+&
               ( -6370.0 * f( l1-1:l2-1 , m , n + 3 , ik ))+&
               ( 3920.0 * f( l1-2:l2-2 , m , n + 3 , ik ))+&
               ( -490.0 * f( l1-3:l2-3 , m , n + 3 , ik ))+&
               ( 6370.0 * f( l1+1:l2+1 , m , n + 3 , ik ))+&
               ( -3920.0 * f( l1+2:l2+2 , m , n + 3 , ik ))+&
               ( 490.0 * f( l1+3:l2+3 , m , n + 3 , ik ))+&
               ( 3510.0 * f( l1-1:l2-1 , m - 1 , n + 3 , ik ))+&
               ( -2160.0 * f( l1-2:l2-2 , m - 1 , n + 3 , ik ))+&
               ( 270.0 * f( l1-3:l2-3 , m - 1 , n + 3 , ik ))+&
               ( -3510.0 * f( l1+1:l2+1 , m - 1 , n + 3 , ik ))+&
               ( 2160.0 * f( l1+2:l2+2 , m - 1 , n + 3 , ik ))+&
               ( -270.0 * f( l1+3:l2+3 , m - 1 , n + 3 , ik ))+&
               ( -351.0 * f( l1-1:l2-1 , m - 2 , n + 3 , ik ))+&
               ( 216.0 * f( l1-2:l2-2 , m - 2 , n + 3 , ik ))+&
               ( -27.0 * f( l1-3:l2-3 , m - 2 , n + 3 , ik ))+&
               ( 351.0 * f( l1+1:l2+1 , m - 2 , n + 3 , ik ))+&
               ( -216.0 * f( l1+2:l2+2 , m - 2 , n + 3 , ik ))+&
               ( 27.0 * f( l1+3:l2+3 , m - 2 , n + 3 , ik ))+&
               ( 26.0 * f( l1-1:l2-1 , m - 3 , n + 3 , ik ))+&
               ( -16.0 * f( l1-2:l2-2 , m - 3 , n + 3 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m - 3 , n + 3 , ik ))+&
               ( -26.0 * f( l1+1:l2+1 , m - 3 , n + 3 , ik ))+&
               ( 16.0 * f( l1+2:l2+2 , m - 3 , n + 3 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m - 3 , n + 3 , ik ))+&
               ( 3510.0 * f( l1-1:l2-1 , m + 1 , n + 3 , ik ))+&
               ( -2160.0 * f( l1-2:l2-2 , m + 1 , n + 3 , ik ))+&
               ( 270.0 * f( l1-3:l2-3 , m + 1 , n + 3 , ik ))+&
               ( -3510.0 * f( l1+1:l2+1 , m + 1 , n + 3 , ik ))+&
               ( 2160.0 * f( l1+2:l2+2 , m + 1 , n + 3 , ik ))+&
               ( -270.0 * f( l1+3:l2+3 , m + 1 , n + 3 , ik ))+&
               ( -351.0 * f( l1-1:l2-1 , m + 2 , n + 3 , ik ))+&
               ( 216.0 * f( l1-2:l2-2 , m + 2 , n + 3 , ik ))+&
               ( -27.0 * f( l1-3:l2-3 , m + 2 , n + 3 , ik ))+&
               ( 351.0 * f( l1+1:l2+1 , m + 2 , n + 3 , ik ))+&
               ( -216.0 * f( l1+2:l2+2 , m + 2 , n + 3 , ik ))+&
               ( 27.0 * f( l1+3:l2+3 , m + 2 , n + 3 , ik ))+&
               ( 26.0 * f( l1-1:l2-1 , m + 3 , n + 3 , ik ))+&
               ( -16.0 * f( l1-2:l2-2 , m + 3 , n + 3 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m + 3 , n + 3 , ik ))+&
               ( -26.0 * f( l1+1:l2+1 , m + 3 , n + 3 , ik ))+&
               ( 16.0 * f( l1+2:l2+2 , m + 3 , n + 3 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m + 3 , n + 3 , ik ))&
               )
        else
          df=0.0
        endif xyz
!        
      elseif (i==1.and.j==3.and.k==2) then
        xzy: if (nxgrid/=1.and.nygrid/=1.and.nzgrid/=1) then
          fac= (1./8.0*dx_1(l1:l2)**3) * (1./180.0*dz_1(n)**2) * (1/60.0*dy_1(m))
          df = fac*(&
               ( 286650.0 * f( l1-1:l2-1 , m - 1 , n , ik ))+&
               ( -176400.0 * f( l1-2:l2-2 , m - 1 , n , ik ))+&
               ( 22050.0 * f( l1-3:l2-3 , m - 1 , n , ik ))+&
               ( -286650.0 * f( l1+1:l2+1 , m - 1 , n , ik ))+&
               ( 176400.0 * f( l1+2:l2+2 , m - 1 , n , ik ))+&
               ( -22050.0 * f( l1+3:l2+3 , m - 1 , n , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m - 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1-2:l2-2 , m - 1 , n - 1 , ik ))+&
               ( -12150.0 * f( l1-3:l2-3 , m - 1 , n - 1 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m - 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1+2:l2+2 , m - 1 , n - 1 , ik ))+&
               ( 12150.0 * f( l1+3:l2+3 , m - 1 , n - 1 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 1 , n - 2 , ik ))+&
               ( -9720.0 * f( l1-2:l2-2 , m - 1 , n - 2 , ik ))+&
               ( 1215.0 * f( l1-3:l2-3 , m - 1 , n - 2 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 1 , n - 2 , ik ))+&
               ( 9720.0 * f( l1+2:l2+2 , m - 1 , n - 2 , ik ))+&
               ( -1215.0 * f( l1+3:l2+3 , m - 1 , n - 2 , ik ))+&
               ( -1170.0 * f( l1-1:l2-1 , m - 1 , n - 3 , ik ))+&
               ( 720.0 * f( l1-2:l2-2 , m - 1 , n - 3 , ik ))+&
               ( -90.0 * f( l1-3:l2-3 , m - 1 , n - 3 , ik ))+&
               ( 1170.0 * f( l1+1:l2+1 , m - 1 , n - 3 , ik ))+&
               ( -720.0 * f( l1+2:l2+2 , m - 1 , n - 3 , ik ))+&
               ( 90.0 * f( l1+3:l2+3 , m - 1 , n - 3 , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m - 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1-2:l2-2 , m - 1 , n + 1 , ik ))+&
               ( -12150.0 * f( l1-3:l2-3 , m - 1 , n + 1 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m - 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1+2:l2+2 , m - 1 , n + 1 , ik ))+&
               ( 12150.0 * f( l1+3:l2+3 , m - 1 , n + 1 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 1 , n + 2 , ik ))+&
               ( -9720.0 * f( l1-2:l2-2 , m - 1 , n + 2 , ik ))+&
               ( 1215.0 * f( l1-3:l2-3 , m - 1 , n + 2 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 1 , n + 2 , ik ))+&
               ( 9720.0 * f( l1+2:l2+2 , m - 1 , n + 2 , ik ))+&
               ( -1215.0 * f( l1+3:l2+3 , m - 1 , n + 2 , ik ))+&
               ( -1170.0 * f( l1-1:l2-1 , m - 1 , n + 3 , ik ))+&
               ( 720.0 * f( l1-2:l2-2 , m - 1 , n + 3 , ik ))+&
               ( -90.0 * f( l1-3:l2-3 , m - 1 , n + 3 , ik ))+&
               ( 1170.0 * f( l1+1:l2+1 , m - 1 , n + 3 , ik ))+&
               ( -720.0 * f( l1+2:l2+2 , m - 1 , n + 3 , ik ))+&
               ( 90.0 * f( l1+3:l2+3 , m - 1 , n + 3 , ik ))+&
               ( -57330.0 * f( l1-1:l2-1 , m - 2 , n , ik ))+&
               ( 35280.0 * f( l1-2:l2-2 , m - 2 , n , ik ))+&
               ( -4410.0 * f( l1-3:l2-3 , m - 2 , n , ik ))+&
               ( 57330.0 * f( l1+1:l2+1 , m - 2 , n , ik ))+&
               ( -35280.0 * f( l1+2:l2+2 , m - 2 , n , ik ))+&
               ( 4410.0 * f( l1+3:l2+3 , m - 2 , n , ik ))+&
               ( 31590.0 * f( l1-1:l2-1 , m - 2 , n - 1 , ik ))+&
               ( -19440.0 * f( l1-2:l2-2 , m - 2 , n - 1 , ik ))+&
               ( 2430.0 * f( l1-3:l2-3 , m - 2 , n - 1 , ik ))+&
               ( -31590.0 * f( l1+1:l2+1 , m - 2 , n - 1 , ik ))+&
               ( 19440.0 * f( l1+2:l2+2 , m - 2 , n - 1 , ik ))+&
               ( -2430.0 * f( l1+3:l2+3 , m - 2 , n - 1 , ik ))+&
               ( -3159.0 * f( l1-1:l2-1 , m - 2 , n - 2 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m - 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1-3:l2-3 , m - 2 , n - 2 , ik ))+&
               ( 3159.0 * f( l1+1:l2+1 , m - 2 , n - 2 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m - 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1+3:l2+3 , m - 2 , n - 2 , ik ))+&
               ( 234.0 * f( l1-1:l2-1 , m - 2 , n - 3 , ik ))+&
               ( -144.0 * f( l1-2:l2-2 , m - 2 , n - 3 , ik ))+&
               ( 18.0 * f( l1-3:l2-3 , m - 2 , n - 3 , ik ))+&
               ( -234.0 * f( l1+1:l2+1 , m - 2 , n - 3 , ik ))+&
               ( 144.0 * f( l1+2:l2+2 , m - 2 , n - 3 , ik ))+&
               ( -18.0 * f( l1+3:l2+3 , m - 2 , n - 3 , ik ))+&
               ( 31590.0 * f( l1-1:l2-1 , m - 2 , n + 1 , ik ))+&
               ( -19440.0 * f( l1-2:l2-2 , m - 2 , n + 1 , ik ))+&
               ( 2430.0 * f( l1-3:l2-3 , m - 2 , n + 1 , ik ))+&
               ( -31590.0 * f( l1+1:l2+1 , m - 2 , n + 1 , ik ))+&
               ( 19440.0 * f( l1+2:l2+2 , m - 2 , n + 1 , ik ))+&
               ( -2430.0 * f( l1+3:l2+3 , m - 2 , n + 1 , ik ))+&
               ( -3159.0 * f( l1-1:l2-1 , m - 2 , n + 2 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m - 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1-3:l2-3 , m - 2 , n + 2 , ik ))+&
               ( 3159.0 * f( l1+1:l2+1 , m - 2 , n + 2 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m - 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1+3:l2+3 , m - 2 , n + 2 , ik ))+&
               ( 234.0 * f( l1-1:l2-1 , m - 2 , n + 3 , ik ))+&
               ( -144.0 * f( l1-2:l2-2 , m - 2 , n + 3 , ik ))+&
               ( 18.0 * f( l1-3:l2-3 , m - 2 , n + 3 , ik ))+&
               ( -234.0 * f( l1+1:l2+1 , m - 2 , n + 3 , ik ))+&
               ( 144.0 * f( l1+2:l2+2 , m - 2 , n + 3 , ik ))+&
               ( -18.0 * f( l1+3:l2+3 , m - 2 , n + 3 , ik ))+&
               ( 6370.0 * f( l1-1:l2-1 , m - 3 , n , ik ))+&
               ( -3920.0 * f( l1-2:l2-2 , m - 3 , n , ik ))+&
               ( 490.0 * f( l1-3:l2-3 , m - 3 , n , ik ))+&
               ( -6370.0 * f( l1+1:l2+1 , m - 3 , n , ik ))+&
               ( 3920.0 * f( l1+2:l2+2 , m - 3 , n , ik ))+&
               ( -490.0 * f( l1+3:l2+3 , m - 3 , n , ik ))+&
               ( -3510.0 * f( l1-1:l2-1 , m - 3 , n - 1 , ik ))+&
               ( 2160.0 * f( l1-2:l2-2 , m - 3 , n - 1 , ik ))+&
               ( -270.0 * f( l1-3:l2-3 , m - 3 , n - 1 , ik ))+&
               ( 3510.0 * f( l1+1:l2+1 , m - 3 , n - 1 , ik ))+&
               ( -2160.0 * f( l1+2:l2+2 , m - 3 , n - 1 , ik ))+&
               ( 270.0 * f( l1+3:l2+3 , m - 3 , n - 1 , ik ))+&
               ( 351.0 * f( l1-1:l2-1 , m - 3 , n - 2 , ik ))+&
               ( -216.0 * f( l1-2:l2-2 , m - 3 , n - 2 , ik ))+&
               ( 27.0 * f( l1-3:l2-3 , m - 3 , n - 2 , ik ))+&
               ( -351.0 * f( l1+1:l2+1 , m - 3 , n - 2 , ik ))+&
               ( 216.0 * f( l1+2:l2+2 , m - 3 , n - 2 , ik ))+&
               ( -27.0 * f( l1+3:l2+3 , m - 3 , n - 2 , ik ))+&
               ( -26.0 * f( l1-1:l2-1 , m - 3 , n - 3 , ik ))+&
               ( 16.0 * f( l1-2:l2-2 , m - 3 , n - 3 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m - 3 , n - 3 , ik ))+&
               ( 26.0 * f( l1+1:l2+1 , m - 3 , n - 3 , ik ))+&
               ( -16.0 * f( l1+2:l2+2 , m - 3 , n - 3 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m - 3 , n - 3 , ik ))+&
               ( -3510.0 * f( l1-1:l2-1 , m - 3 , n + 1 , ik ))+&
               ( 2160.0 * f( l1-2:l2-2 , m - 3 , n + 1 , ik ))+&
               ( -270.0 * f( l1-3:l2-3 , m - 3 , n + 1 , ik ))+&
               ( 3510.0 * f( l1+1:l2+1 , m - 3 , n + 1 , ik ))+&
               ( -2160.0 * f( l1+2:l2+2 , m - 3 , n + 1 , ik ))+&
               ( 270.0 * f( l1+3:l2+3 , m - 3 , n + 1 , ik ))+&
               ( 351.0 * f( l1-1:l2-1 , m - 3 , n + 2 , ik ))+&
               ( -216.0 * f( l1-2:l2-2 , m - 3 , n + 2 , ik ))+&
               ( 27.0 * f( l1-3:l2-3 , m - 3 , n + 2 , ik ))+&
               ( -351.0 * f( l1+1:l2+1 , m - 3 , n + 2 , ik ))+&
               ( 216.0 * f( l1+2:l2+2 , m - 3 , n + 2 , ik ))+&
               ( -27.0 * f( l1+3:l2+3 , m - 3 , n + 2 , ik ))+&
               ( -26.0 * f( l1-1:l2-1 , m - 3 , n + 3 , ik ))+&
               ( 16.0 * f( l1-2:l2-2 , m - 3 , n + 3 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m - 3 , n + 3 , ik ))+&
               ( 26.0 * f( l1+1:l2+1 , m - 3 , n + 3 , ik ))+&
               ( -16.0 * f( l1+2:l2+2 , m - 3 , n + 3 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m - 3 , n + 3 , ik ))+&
               ( -286650.0 * f( l1-1:l2-1 , m + 1 , n , ik ))+&
               ( 176400.0 * f( l1-2:l2-2 , m + 1 , n , ik ))+&
               ( -22050.0 * f( l1-3:l2-3 , m + 1 , n , ik ))+&
               ( 286650.0 * f( l1+1:l2+1 , m + 1 , n , ik ))+&
               ( -176400.0 * f( l1+2:l2+2 , m + 1 , n , ik ))+&
               ( 22050.0 * f( l1+3:l2+3 , m + 1 , n , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m + 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1-2:l2-2 , m + 1 , n - 1 , ik ))+&
               ( 12150.0 * f( l1-3:l2-3 , m + 1 , n - 1 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m + 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1+2:l2+2 , m + 1 , n - 1 , ik ))+&
               ( -12150.0 * f( l1+3:l2+3 , m + 1 , n - 1 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 1 , n - 2 , ik ))+&
               ( 9720.0 * f( l1-2:l2-2 , m + 1 , n - 2 , ik ))+&
               ( -1215.0 * f( l1-3:l2-3 , m + 1 , n - 2 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 1 , n - 2 , ik ))+&
               ( -9720.0 * f( l1+2:l2+2 , m + 1 , n - 2 , ik ))+&
               ( 1215.0 * f( l1+3:l2+3 , m + 1 , n - 2 , ik ))+&
               ( 1170.0 * f( l1-1:l2-1 , m + 1 , n - 3 , ik ))+&
               ( -720.0 * f( l1-2:l2-2 , m + 1 , n - 3 , ik ))+&
               ( 90.0 * f( l1-3:l2-3 , m + 1 , n - 3 , ik ))+&
               ( -1170.0 * f( l1+1:l2+1 , m + 1 , n - 3 , ik ))+&
               ( 720.0 * f( l1+2:l2+2 , m + 1 , n - 3 , ik ))+&
               ( -90.0 * f( l1+3:l2+3 , m + 1 , n - 3 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m + 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1-2:l2-2 , m + 1 , n + 1 , ik ))+&
               ( 12150.0 * f( l1-3:l2-3 , m + 1 , n + 1 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m + 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1+2:l2+2 , m + 1 , n + 1 , ik ))+&
               ( -12150.0 * f( l1+3:l2+3 , m + 1 , n + 1 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 1 , n + 2 , ik ))+&
               ( 9720.0 * f( l1-2:l2-2 , m + 1 , n + 2 , ik ))+&
               ( -1215.0 * f( l1-3:l2-3 , m + 1 , n + 2 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 1 , n + 2 , ik ))+&
               ( -9720.0 * f( l1+2:l2+2 , m + 1 , n + 2 , ik ))+&
               ( 1215.0 * f( l1+3:l2+3 , m + 1 , n + 2 , ik ))+&
               ( 1170.0 * f( l1-1:l2-1 , m + 1 , n + 3 , ik ))+&
               ( -720.0 * f( l1-2:l2-2 , m + 1 , n + 3 , ik ))+&
               ( 90.0 * f( l1-3:l2-3 , m + 1 , n + 3 , ik ))+&
               ( -1170.0 * f( l1+1:l2+1 , m + 1 , n + 3 , ik ))+&
               ( 720.0 * f( l1+2:l2+2 , m + 1 , n + 3 , ik ))+&
               ( -90.0 * f( l1+3:l2+3 , m + 1 , n + 3 , ik ))+&
               ( 57330.0 * f( l1-1:l2-1 , m + 2 , n , ik ))+&
               ( -35280.0 * f( l1-2:l2-2 , m + 2 , n , ik ))+&
               ( 4410.0 * f( l1-3:l2-3 , m + 2 , n , ik ))+&
               ( -57330.0 * f( l1+1:l2+1 , m + 2 , n , ik ))+&
               ( 35280.0 * f( l1+2:l2+2 , m + 2 , n , ik ))+&
               ( -4410.0 * f( l1+3:l2+3 , m + 2 , n , ik ))+&
               ( -31590.0 * f( l1-1:l2-1 , m + 2 , n - 1 , ik ))+&
               ( 19440.0 * f( l1-2:l2-2 , m + 2 , n - 1 , ik ))+&
               ( -2430.0 * f( l1-3:l2-3 , m + 2 , n - 1 , ik ))+&
               ( 31590.0 * f( l1+1:l2+1 , m + 2 , n - 1 , ik ))+&
               ( -19440.0 * f( l1+2:l2+2 , m + 2 , n - 1 , ik ))+&
               ( 2430.0 * f( l1+3:l2+3 , m + 2 , n - 1 , ik ))+&
               ( 3159.0 * f( l1-1:l2-1 , m + 2 , n - 2 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m + 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1-3:l2-3 , m + 2 , n - 2 , ik ))+&
               ( -3159.0 * f( l1+1:l2+1 , m + 2 , n - 2 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m + 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1+3:l2+3 , m + 2 , n - 2 , ik ))+&
               ( -234.0 * f( l1-1:l2-1 , m + 2 , n - 3 , ik ))+&
               ( 144.0 * f( l1-2:l2-2 , m + 2 , n - 3 , ik ))+&
               ( -18.0 * f( l1-3:l2-3 , m + 2 , n - 3 , ik ))+&
               ( 234.0 * f( l1+1:l2+1 , m + 2 , n - 3 , ik ))+&
               ( -144.0 * f( l1+2:l2+2 , m + 2 , n - 3 , ik ))+&
               ( 18.0 * f( l1+3:l2+3 , m + 2 , n - 3 , ik ))+&
               ( -31590.0 * f( l1-1:l2-1 , m + 2 , n + 1 , ik ))+&
               ( 19440.0 * f( l1-2:l2-2 , m + 2 , n + 1 , ik ))+&
               ( -2430.0 * f( l1-3:l2-3 , m + 2 , n + 1 , ik ))+&
               ( 31590.0 * f( l1+1:l2+1 , m + 2 , n + 1 , ik ))+&
               ( -19440.0 * f( l1+2:l2+2 , m + 2 , n + 1 , ik ))+&
               ( 2430.0 * f( l1+3:l2+3 , m + 2 , n + 1 , ik ))+&
               ( 3159.0 * f( l1-1:l2-1 , m + 2 , n + 2 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m + 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1-3:l2-3 , m + 2 , n + 2 , ik ))+&
               ( -3159.0 * f( l1+1:l2+1 , m + 2 , n + 2 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m + 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1+3:l2+3 , m + 2 , n + 2 , ik ))+&
               ( -234.0 * f( l1-1:l2-1 , m + 2 , n + 3 , ik ))+&
               ( 144.0 * f( l1-2:l2-2 , m + 2 , n + 3 , ik ))+&
               ( -18.0 * f( l1-3:l2-3 , m + 2 , n + 3 , ik ))+&
               ( 234.0 * f( l1+1:l2+1 , m + 2 , n + 3 , ik ))+&
               ( -144.0 * f( l1+2:l2+2 , m + 2 , n + 3 , ik ))+&
               ( 18.0 * f( l1+3:l2+3 , m + 2 , n + 3 , ik ))+&
               ( -6370.0 * f( l1-1:l2-1 , m + 3 , n , ik ))+&
               ( 3920.0 * f( l1-2:l2-2 , m + 3 , n , ik ))+&
               ( -490.0 * f( l1-3:l2-3 , m + 3 , n , ik ))+&
               ( 6370.0 * f( l1+1:l2+1 , m + 3 , n , ik ))+&
               ( -3920.0 * f( l1+2:l2+2 , m + 3 , n , ik ))+&
               ( 490.0 * f( l1+3:l2+3 , m + 3 , n , ik ))+&
               ( 3510.0 * f( l1-1:l2-1 , m + 3 , n - 1 , ik ))+&
               ( -2160.0 * f( l1-2:l2-2 , m + 3 , n - 1 , ik ))+&
               ( 270.0 * f( l1-3:l2-3 , m + 3 , n - 1 , ik ))+&
               ( -3510.0 * f( l1+1:l2+1 , m + 3 , n - 1 , ik ))+&
               ( 2160.0 * f( l1+2:l2+2 , m + 3 , n - 1 , ik ))+&
               ( -270.0 * f( l1+3:l2+3 , m + 3 , n - 1 , ik ))+&
               ( -351.0 * f( l1-1:l2-1 , m + 3 , n - 2 , ik ))+&
               ( 216.0 * f( l1-2:l2-2 , m + 3 , n - 2 , ik ))+&
               ( -27.0 * f( l1-3:l2-3 , m + 3 , n - 2 , ik ))+&
               ( 351.0 * f( l1+1:l2+1 , m + 3 , n - 2 , ik ))+&
               ( -216.0 * f( l1+2:l2+2 , m + 3 , n - 2 , ik ))+&
               ( 27.0 * f( l1+3:l2+3 , m + 3 , n - 2 , ik ))+&
               ( 26.0 * f( l1-1:l2-1 , m + 3 , n - 3 , ik ))+&
               ( -16.0 * f( l1-2:l2-2 , m + 3 , n - 3 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m + 3 , n - 3 , ik ))+&
               ( -26.0 * f( l1+1:l2+1 , m + 3 , n - 3 , ik ))+&
               ( 16.0 * f( l1+2:l2+2 , m + 3 , n - 3 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m + 3 , n - 3 , ik ))+&
               ( 3510.0 * f( l1-1:l2-1 , m + 3 , n + 1 , ik ))+&
               ( -2160.0 * f( l1-2:l2-2 , m + 3 , n + 1 , ik ))+&
               ( 270.0 * f( l1-3:l2-3 , m + 3 , n + 1 , ik ))+&
               ( -3510.0 * f( l1+1:l2+1 , m + 3 , n + 1 , ik ))+&
               ( 2160.0 * f( l1+2:l2+2 , m + 3 , n + 1 , ik ))+&
               ( -270.0 * f( l1+3:l2+3 , m + 3 , n + 1 , ik ))+&
               ( -351.0 * f( l1-1:l2-1 , m + 3 , n + 2 , ik ))+&
               ( 216.0 * f( l1-2:l2-2 , m + 3 , n + 2 , ik ))+&
               ( -27.0 * f( l1-3:l2-3 , m + 3 , n + 2 , ik ))+&
               ( 351.0 * f( l1+1:l2+1 , m + 3 , n + 2 , ik ))+&
               ( -216.0 * f( l1+2:l2+2 , m + 3 , n + 2 , ik ))+&
               ( 27.0 * f( l1+3:l2+3 , m + 3 , n + 2 , ik ))+&
               ( 26.0 * f( l1-1:l2-1 , m + 3 , n + 3 , ik ))+&
               ( -16.0 * f( l1-2:l2-2 , m + 3 , n + 3 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m + 3 , n + 3 , ik ))+&
               ( -26.0 * f( l1+1:l2+1 , m + 3 , n + 3 , ik ))+&
               ( 16.0 * f( l1+2:l2+2 , m + 3 , n + 3 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m + 3 , n + 3 , ik ))&
               )
        else
          df=0.0
        endif xzy
!        
      elseif (i==2.and.j==1.and.k==3) then
        yxz: if (nxgrid/=1.and.nygrid/=1.and.nzgrid/=1) then
          fac= (1./8.0*dy_1(m)**3) * (1./180.0*dx_1(l1:l2)**2) * (1/60.0*dz_1(n))
          df = fac*(&
               ( 286650.0 * f( l1:l2 , m - 1 , n - 1 , ik ))+&
               ( -176400.0 * f( l1:l2 , m - 2 , n - 1 , ik ))+&
               ( 22050.0 * f( l1:l2 , m - 3 , n - 1 , ik ))+&
               ( -286650.0 * f( l1:l2 , m + 1 , n - 1 , ik ))+&
               ( 176400.0 * f( l1:l2 , m + 2 , n - 1 , ik ))+&
               ( -22050.0 * f( l1:l2 , m + 3 , n - 1 , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m - 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1-1:l2-1 , m - 2 , n - 1 , ik ))+&
               ( -12150.0 * f( l1-1:l2-1 , m - 3 , n - 1 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m + 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1-1:l2-1 , m + 2 , n - 1 , ik ))+&
               ( 12150.0 * f( l1-1:l2-1 , m + 3 , n - 1 , ik ))+&
               ( 15795.0 * f( l1-2:l2-2 , m - 1 , n - 1 , ik ))+&
               ( -9720.0 * f( l1-2:l2-2 , m - 2 , n - 1 , ik ))+&
               ( 1215.0 * f( l1-2:l2-2 , m - 3 , n - 1 , ik ))+&
               ( -15795.0 * f( l1-2:l2-2 , m + 1 , n - 1 , ik ))+&
               ( 9720.0 * f( l1-2:l2-2 , m + 2 , n - 1 , ik ))+&
               ( -1215.0 * f( l1-2:l2-2 , m + 3 , n - 1 , ik ))+&
               ( -1170.0 * f( l1-3:l2-3 , m - 1 , n - 1 , ik ))+&
               ( 720.0 * f( l1-3:l2-3 , m - 2 , n - 1 , ik ))+&
               ( -90.0 * f( l1-3:l2-3 , m - 3 , n - 1 , ik ))+&
               ( 1170.0 * f( l1-3:l2-3 , m + 1 , n - 1 , ik ))+&
               ( -720.0 * f( l1-3:l2-3 , m + 2 , n - 1 , ik ))+&
               ( 90.0 * f( l1-3:l2-3 , m + 3 , n - 1 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m - 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1+1:l2+1 , m - 2 , n - 1 , ik ))+&
               ( -12150.0 * f( l1+1:l2+1 , m - 3 , n - 1 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m + 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1+1:l2+1 , m + 2 , n - 1 , ik ))+&
               ( 12150.0 * f( l1+1:l2+1 , m + 3 , n - 1 , ik ))+&
               ( 15795.0 * f( l1+2:l2+2 , m - 1 , n - 1 , ik ))+&
               ( -9720.0 * f( l1+2:l2+2 , m - 2 , n - 1 , ik ))+&
               ( 1215.0 * f( l1+2:l2+2 , m - 3 , n - 1 , ik ))+&
               ( -15795.0 * f( l1+2:l2+2 , m + 1 , n - 1 , ik ))+&
               ( 9720.0 * f( l1+2:l2+2 , m + 2 , n - 1 , ik ))+&
               ( -1215.0 * f( l1+2:l2+2 , m + 3 , n - 1 , ik ))+&
               ( -1170.0 * f( l1+3:l2+3 , m - 1 , n - 1 , ik ))+&
               ( 720.0 * f( l1+3:l2+3 , m - 2 , n - 1 , ik ))+&
               ( -90.0 * f( l1+3:l2+3 , m - 3 , n - 1 , ik ))+&
               ( 1170.0 * f( l1+3:l2+3 , m + 1 , n - 1 , ik ))+&
               ( -720.0 * f( l1+3:l2+3 , m + 2 , n - 1 , ik ))+&
               ( 90.0 * f( l1+3:l2+3 , m + 3 , n - 1 , ik ))+&
               ( -57330.0 * f( l1:l2 , m - 1 , n - 2 , ik ))+&
               ( 35280.0 * f( l1:l2 , m - 2 , n - 2 , ik ))+&
               ( -4410.0 * f( l1:l2 , m - 3 , n - 2 , ik ))+&
               ( 57330.0 * f( l1:l2 , m + 1 , n - 2 , ik ))+&
               ( -35280.0 * f( l1:l2 , m + 2 , n - 2 , ik ))+&
               ( 4410.0 * f( l1:l2 , m + 3 , n - 2 , ik ))+&
               ( 31590.0 * f( l1-1:l2-1 , m - 1 , n - 2 , ik ))+&
               ( -19440.0 * f( l1-1:l2-1 , m - 2 , n - 2 , ik ))+&
               ( 2430.0 * f( l1-1:l2-1 , m - 3 , n - 2 , ik ))+&
               ( -31590.0 * f( l1-1:l2-1 , m + 1 , n - 2 , ik ))+&
               ( 19440.0 * f( l1-1:l2-1 , m + 2 , n - 2 , ik ))+&
               ( -2430.0 * f( l1-1:l2-1 , m + 3 , n - 2 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m - 1 , n - 2 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m - 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1-2:l2-2 , m - 3 , n - 2 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m + 1 , n - 2 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m + 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1-2:l2-2 , m + 3 , n - 2 , ik ))+&
               ( 234.0 * f( l1-3:l2-3 , m - 1 , n - 2 , ik ))+&
               ( -144.0 * f( l1-3:l2-3 , m - 2 , n - 2 , ik ))+&
               ( 18.0 * f( l1-3:l2-3 , m - 3 , n - 2 , ik ))+&
               ( -234.0 * f( l1-3:l2-3 , m + 1 , n - 2 , ik ))+&
               ( 144.0 * f( l1-3:l2-3 , m + 2 , n - 2 , ik ))+&
               ( -18.0 * f( l1-3:l2-3 , m + 3 , n - 2 , ik ))+&
               ( 31590.0 * f( l1+1:l2+1 , m - 1 , n - 2 , ik ))+&
               ( -19440.0 * f( l1+1:l2+1 , m - 2 , n - 2 , ik ))+&
               ( 2430.0 * f( l1+1:l2+1 , m - 3 , n - 2 , ik ))+&
               ( -31590.0 * f( l1+1:l2+1 , m + 1 , n - 2 , ik ))+&
               ( 19440.0 * f( l1+1:l2+1 , m + 2 , n - 2 , ik ))+&
               ( -2430.0 * f( l1+1:l2+1 , m + 3 , n - 2 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m - 1 , n - 2 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m - 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1+2:l2+2 , m - 3 , n - 2 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m + 1 , n - 2 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m + 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1+2:l2+2 , m + 3 , n - 2 , ik ))+&
               ( 234.0 * f( l1+3:l2+3 , m - 1 , n - 2 , ik ))+&
               ( -144.0 * f( l1+3:l2+3 , m - 2 , n - 2 , ik ))+&
               ( 18.0 * f( l1+3:l2+3 , m - 3 , n - 2 , ik ))+&
               ( -234.0 * f( l1+3:l2+3 , m + 1 , n - 2 , ik ))+&
               ( 144.0 * f( l1+3:l2+3 , m + 2 , n - 2 , ik ))+&
               ( -18.0 * f( l1+3:l2+3 , m + 3 , n - 2 , ik ))+&
               ( 6370.0 * f( l1:l2 , m - 1 , n - 3 , ik ))+&
               ( -3920.0 * f( l1:l2 , m - 2 , n - 3 , ik ))+&
               ( 490.0 * f( l1:l2 , m - 3 , n - 3 , ik ))+&
               ( -6370.0 * f( l1:l2 , m + 1 , n - 3 , ik ))+&
               ( 3920.0 * f( l1:l2 , m + 2 , n - 3 , ik ))+&
               ( -490.0 * f( l1:l2 , m + 3 , n - 3 , ik ))+&
               ( -3510.0 * f( l1-1:l2-1 , m - 1 , n - 3 , ik ))+&
               ( 2160.0 * f( l1-1:l2-1 , m - 2 , n - 3 , ik ))+&
               ( -270.0 * f( l1-1:l2-1 , m - 3 , n - 3 , ik ))+&
               ( 3510.0 * f( l1-1:l2-1 , m + 1 , n - 3 , ik ))+&
               ( -2160.0 * f( l1-1:l2-1 , m + 2 , n - 3 , ik ))+&
               ( 270.0 * f( l1-1:l2-1 , m + 3 , n - 3 , ik ))+&
               ( 351.0 * f( l1-2:l2-2 , m - 1 , n - 3 , ik ))+&
               ( -216.0 * f( l1-2:l2-2 , m - 2 , n - 3 , ik ))+&
               ( 27.0 * f( l1-2:l2-2 , m - 3 , n - 3 , ik ))+&
               ( -351.0 * f( l1-2:l2-2 , m + 1 , n - 3 , ik ))+&
               ( 216.0 * f( l1-2:l2-2 , m + 2 , n - 3 , ik ))+&
               ( -27.0 * f( l1-2:l2-2 , m + 3 , n - 3 , ik ))+&
               ( -26.0 * f( l1-3:l2-3 , m - 1 , n - 3 , ik ))+&
               ( 16.0 * f( l1-3:l2-3 , m - 2 , n - 3 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m - 3 , n - 3 , ik ))+&
               ( 26.0 * f( l1-3:l2-3 , m + 1 , n - 3 , ik ))+&
               ( -16.0 * f( l1-3:l2-3 , m + 2 , n - 3 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m + 3 , n - 3 , ik ))+&
               ( -3510.0 * f( l1+1:l2+1 , m - 1 , n - 3 , ik ))+&
               ( 2160.0 * f( l1+1:l2+1 , m - 2 , n - 3 , ik ))+&
               ( -270.0 * f( l1+1:l2+1 , m - 3 , n - 3 , ik ))+&
               ( 3510.0 * f( l1+1:l2+1 , m + 1 , n - 3 , ik ))+&
               ( -2160.0 * f( l1+1:l2+1 , m + 2 , n - 3 , ik ))+&
               ( 270.0 * f( l1+1:l2+1 , m + 3 , n - 3 , ik ))+&
               ( 351.0 * f( l1+2:l2+2 , m - 1 , n - 3 , ik ))+&
               ( -216.0 * f( l1+2:l2+2 , m - 2 , n - 3 , ik ))+&
               ( 27.0 * f( l1+2:l2+2 , m - 3 , n - 3 , ik ))+&
               ( -351.0 * f( l1+2:l2+2 , m + 1 , n - 3 , ik ))+&
               ( 216.0 * f( l1+2:l2+2 , m + 2 , n - 3 , ik ))+&
               ( -27.0 * f( l1+2:l2+2 , m + 3 , n - 3 , ik ))+&
               ( -26.0 * f( l1+3:l2+3 , m - 1 , n - 3 , ik ))+&
               ( 16.0 * f( l1+3:l2+3 , m - 2 , n - 3 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m - 3 , n - 3 , ik ))+&
               ( 26.0 * f( l1+3:l2+3 , m + 1 , n - 3 , ik ))+&
               ( -16.0 * f( l1+3:l2+3 , m + 2 , n - 3 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m + 3 , n - 3 , ik ))+&
               ( -286650.0 * f( l1:l2 , m - 1 , n + 1 , ik ))+&
               ( 176400.0 * f( l1:l2 , m - 2 , n + 1 , ik ))+&
               ( -22050.0 * f( l1:l2 , m - 3 , n + 1 , ik ))+&
               ( 286650.0 * f( l1:l2 , m + 1 , n + 1 , ik ))+&
               ( -176400.0 * f( l1:l2 , m + 2 , n + 1 , ik ))+&
               ( 22050.0 * f( l1:l2 , m + 3 , n + 1 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m - 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1-1:l2-1 , m - 2 , n + 1 , ik ))+&
               ( 12150.0 * f( l1-1:l2-1 , m - 3 , n + 1 , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m + 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1-1:l2-1 , m + 2 , n + 1 , ik ))+&
               ( -12150.0 * f( l1-1:l2-1 , m + 3 , n + 1 , ik ))+&
               ( -15795.0 * f( l1-2:l2-2 , m - 1 , n + 1 , ik ))+&
               ( 9720.0 * f( l1-2:l2-2 , m - 2 , n + 1 , ik ))+&
               ( -1215.0 * f( l1-2:l2-2 , m - 3 , n + 1 , ik ))+&
               ( 15795.0 * f( l1-2:l2-2 , m + 1 , n + 1 , ik ))+&
               ( -9720.0 * f( l1-2:l2-2 , m + 2 , n + 1 , ik ))+&
               ( 1215.0 * f( l1-2:l2-2 , m + 3 , n + 1 , ik ))+&
               ( 1170.0 * f( l1-3:l2-3 , m - 1 , n + 1 , ik ))+&
               ( -720.0 * f( l1-3:l2-3 , m - 2 , n + 1 , ik ))+&
               ( 90.0 * f( l1-3:l2-3 , m - 3 , n + 1 , ik ))+&
               ( -1170.0 * f( l1-3:l2-3 , m + 1 , n + 1 , ik ))+&
               ( 720.0 * f( l1-3:l2-3 , m + 2 , n + 1 , ik ))+&
               ( -90.0 * f( l1-3:l2-3 , m + 3 , n + 1 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m - 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1+1:l2+1 , m - 2 , n + 1 , ik ))+&
               ( 12150.0 * f( l1+1:l2+1 , m - 3 , n + 1 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m + 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1+1:l2+1 , m + 2 , n + 1 , ik ))+&
               ( -12150.0 * f( l1+1:l2+1 , m + 3 , n + 1 , ik ))+&
               ( -15795.0 * f( l1+2:l2+2 , m - 1 , n + 1 , ik ))+&
               ( 9720.0 * f( l1+2:l2+2 , m - 2 , n + 1 , ik ))+&
               ( -1215.0 * f( l1+2:l2+2 , m - 3 , n + 1 , ik ))+&
               ( 15795.0 * f( l1+2:l2+2 , m + 1 , n + 1 , ik ))+&
               ( -9720.0 * f( l1+2:l2+2 , m + 2 , n + 1 , ik ))+&
               ( 1215.0 * f( l1+2:l2+2 , m + 3 , n + 1 , ik ))+&
               ( 1170.0 * f( l1+3:l2+3 , m - 1 , n + 1 , ik ))+&
               ( -720.0 * f( l1+3:l2+3 , m - 2 , n + 1 , ik ))+&
               ( 90.0 * f( l1+3:l2+3 , m - 3 , n + 1 , ik ))+&
               ( -1170.0 * f( l1+3:l2+3 , m + 1 , n + 1 , ik ))+&
               ( 720.0 * f( l1+3:l2+3 , m + 2 , n + 1 , ik ))+&
               ( -90.0 * f( l1+3:l2+3 , m + 3 , n + 1 , ik ))+&
               ( 57330.0 * f( l1:l2 , m - 1 , n + 2 , ik ))+&
               ( -35280.0 * f( l1:l2 , m - 2 , n + 2 , ik ))+&
               ( 4410.0 * f( l1:l2 , m - 3 , n + 2 , ik ))+&
               ( -57330.0 * f( l1:l2 , m + 1 , n + 2 , ik ))+&
               ( 35280.0 * f( l1:l2 , m + 2 , n + 2 , ik ))+&
               ( -4410.0 * f( l1:l2 , m + 3 , n + 2 , ik ))+&
               ( -31590.0 * f( l1-1:l2-1 , m - 1 , n + 2 , ik ))+&
               ( 19440.0 * f( l1-1:l2-1 , m - 2 , n + 2 , ik ))+&
               ( -2430.0 * f( l1-1:l2-1 , m - 3 , n + 2 , ik ))+&
               ( 31590.0 * f( l1-1:l2-1 , m + 1 , n + 2 , ik ))+&
               ( -19440.0 * f( l1-1:l2-1 , m + 2 , n + 2 , ik ))+&
               ( 2430.0 * f( l1-1:l2-1 , m + 3 , n + 2 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m - 1 , n + 2 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m - 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1-2:l2-2 , m - 3 , n + 2 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m + 1 , n + 2 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m + 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1-2:l2-2 , m + 3 , n + 2 , ik ))+&
               ( -234.0 * f( l1-3:l2-3 , m - 1 , n + 2 , ik ))+&
               ( 144.0 * f( l1-3:l2-3 , m - 2 , n + 2 , ik ))+&
               ( -18.0 * f( l1-3:l2-3 , m - 3 , n + 2 , ik ))+&
               ( 234.0 * f( l1-3:l2-3 , m + 1 , n + 2 , ik ))+&
               ( -144.0 * f( l1-3:l2-3 , m + 2 , n + 2 , ik ))+&
               ( 18.0 * f( l1-3:l2-3 , m + 3 , n + 2 , ik ))+&
               ( -31590.0 * f( l1+1:l2+1 , m - 1 , n + 2 , ik ))+&
               ( 19440.0 * f( l1+1:l2+1 , m - 2 , n + 2 , ik ))+&
               ( -2430.0 * f( l1+1:l2+1 , m - 3 , n + 2 , ik ))+&
               ( 31590.0 * f( l1+1:l2+1 , m + 1 , n + 2 , ik ))+&
               ( -19440.0 * f( l1+1:l2+1 , m + 2 , n + 2 , ik ))+&
               ( 2430.0 * f( l1+1:l2+1 , m + 3 , n + 2 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m - 1 , n + 2 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m - 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1+2:l2+2 , m - 3 , n + 2 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m + 1 , n + 2 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m + 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1+2:l2+2 , m + 3 , n + 2 , ik ))+&
               ( -234.0 * f( l1+3:l2+3 , m - 1 , n + 2 , ik ))+&
               ( 144.0 * f( l1+3:l2+3 , m - 2 , n + 2 , ik ))+&
               ( -18.0 * f( l1+3:l2+3 , m - 3 , n + 2 , ik ))+&
               ( 234.0 * f( l1+3:l2+3 , m + 1 , n + 2 , ik ))+&
               ( -144.0 * f( l1+3:l2+3 , m + 2 , n + 2 , ik ))+&
               ( 18.0 * f( l1+3:l2+3 , m + 3 , n + 2 , ik ))+&
               ( -6370.0 * f( l1:l2 , m - 1 , n + 3 , ik ))+&
               ( 3920.0 * f( l1:l2 , m - 2 , n + 3 , ik ))+&
               ( -490.0 * f( l1:l2 , m - 3 , n + 3 , ik ))+&
               ( 6370.0 * f( l1:l2 , m + 1 , n + 3 , ik ))+&
               ( -3920.0 * f( l1:l2 , m + 2 , n + 3 , ik ))+&
               ( 490.0 * f( l1:l2 , m + 3 , n + 3 , ik ))+&
               ( 3510.0 * f( l1-1:l2-1 , m - 1 , n + 3 , ik ))+&
               ( -2160.0 * f( l1-1:l2-1 , m - 2 , n + 3 , ik ))+&
               ( 270.0 * f( l1-1:l2-1 , m - 3 , n + 3 , ik ))+&
               ( -3510.0 * f( l1-1:l2-1 , m + 1 , n + 3 , ik ))+&
               ( 2160.0 * f( l1-1:l2-1 , m + 2 , n + 3 , ik ))+&
               ( -270.0 * f( l1-1:l2-1 , m + 3 , n + 3 , ik ))+&
               ( -351.0 * f( l1-2:l2-2 , m - 1 , n + 3 , ik ))+&
               ( 216.0 * f( l1-2:l2-2 , m - 2 , n + 3 , ik ))+&
               ( -27.0 * f( l1-2:l2-2 , m - 3 , n + 3 , ik ))+&
               ( 351.0 * f( l1-2:l2-2 , m + 1 , n + 3 , ik ))+&
               ( -216.0 * f( l1-2:l2-2 , m + 2 , n + 3 , ik ))+&
               ( 27.0 * f( l1-2:l2-2 , m + 3 , n + 3 , ik ))+&
               ( 26.0 * f( l1-3:l2-3 , m - 1 , n + 3 , ik ))+&
               ( -16.0 * f( l1-3:l2-3 , m - 2 , n + 3 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m - 3 , n + 3 , ik ))+&
               ( -26.0 * f( l1-3:l2-3 , m + 1 , n + 3 , ik ))+&
               ( 16.0 * f( l1-3:l2-3 , m + 2 , n + 3 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m + 3 , n + 3 , ik ))+&
               ( 3510.0 * f( l1+1:l2+1 , m - 1 , n + 3 , ik ))+&
               ( -2160.0 * f( l1+1:l2+1 , m - 2 , n + 3 , ik ))+&
               ( 270.0 * f( l1+1:l2+1 , m - 3 , n + 3 , ik ))+&
               ( -3510.0 * f( l1+1:l2+1 , m + 1 , n + 3 , ik ))+&
               ( 2160.0 * f( l1+1:l2+1 , m + 2 , n + 3 , ik ))+&
               ( -270.0 * f( l1+1:l2+1 , m + 3 , n + 3 , ik ))+&
               ( -351.0 * f( l1+2:l2+2 , m - 1 , n + 3 , ik ))+&
               ( 216.0 * f( l1+2:l2+2 , m - 2 , n + 3 , ik ))+&
               ( -27.0 * f( l1+2:l2+2 , m - 3 , n + 3 , ik ))+&
               ( 351.0 * f( l1+2:l2+2 , m + 1 , n + 3 , ik ))+&
               ( -216.0 * f( l1+2:l2+2 , m + 2 , n + 3 , ik ))+&
               ( 27.0 * f( l1+2:l2+2 , m + 3 , n + 3 , ik ))+&
               ( 26.0 * f( l1+3:l2+3 , m - 1 , n + 3 , ik ))+&
               ( -16.0 * f( l1+3:l2+3 , m - 2 , n + 3 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m - 3 , n + 3 , ik ))+&
               ( -26.0 * f( l1+3:l2+3 , m + 1 , n + 3 , ik ))+&
               ( 16.0 * f( l1+3:l2+3 , m + 2 , n + 3 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m + 3 , n + 3 , ik ))&
               )
        else
          df=0.0
        endif yxz
!        
      elseif (i==2.and.j==3.and.k==1) then
        yzx: if (nxgrid/=1.and.nygrid/=1.and.nzgrid/=1) then
          fac= (1./8.0*dy_1(m)**3) * (1./180.0*dz_1(n)**2) * (1/60.0*dx_1(l1:l2))
          df = fac*(&
               ( 286650.0 * f( l1-1:l2-1 , m - 1 , n , ik ))+&
               ( -176400.0 * f( l1-1:l2-1 , m - 2 , n , ik ))+&
               ( 22050.0 * f( l1-1:l2-1 , m - 3 , n , ik ))+&
               ( -286650.0 * f( l1-1:l2-1 , m + 1 , n , ik ))+&
               ( 176400.0 * f( l1-1:l2-1 , m + 2 , n , ik ))+&
               ( -22050.0 * f( l1-1:l2-1 , m + 3 , n , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m - 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1-1:l2-1 , m - 2 , n - 1 , ik ))+&
               ( -12150.0 * f( l1-1:l2-1 , m - 3 , n - 1 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m + 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1-1:l2-1 , m + 2 , n - 1 , ik ))+&
               ( 12150.0 * f( l1-1:l2-1 , m + 3 , n - 1 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 1 , n - 2 , ik ))+&
               ( -9720.0 * f( l1-1:l2-1 , m - 2 , n - 2 , ik ))+&
               ( 1215.0 * f( l1-1:l2-1 , m - 3 , n - 2 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 1 , n - 2 , ik ))+&
               ( 9720.0 * f( l1-1:l2-1 , m + 2 , n - 2 , ik ))+&
               ( -1215.0 * f( l1-1:l2-1 , m + 3 , n - 2 , ik ))+&
               ( -1170.0 * f( l1-1:l2-1 , m - 1 , n - 3 , ik ))+&
               ( 720.0 * f( l1-1:l2-1 , m - 2 , n - 3 , ik ))+&
               ( -90.0 * f( l1-1:l2-1 , m - 3 , n - 3 , ik ))+&
               ( 1170.0 * f( l1-1:l2-1 , m + 1 , n - 3 , ik ))+&
               ( -720.0 * f( l1-1:l2-1 , m + 2 , n - 3 , ik ))+&
               ( 90.0 * f( l1-1:l2-1 , m + 3 , n - 3 , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m - 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1-1:l2-1 , m - 2 , n + 1 , ik ))+&
               ( -12150.0 * f( l1-1:l2-1 , m - 3 , n + 1 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m + 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1-1:l2-1 , m + 2 , n + 1 , ik ))+&
               ( 12150.0 * f( l1-1:l2-1 , m + 3 , n + 1 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 1 , n + 2 , ik ))+&
               ( -9720.0 * f( l1-1:l2-1 , m - 2 , n + 2 , ik ))+&
               ( 1215.0 * f( l1-1:l2-1 , m - 3 , n + 2 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 1 , n + 2 , ik ))+&
               ( 9720.0 * f( l1-1:l2-1 , m + 2 , n + 2 , ik ))+&
               ( -1215.0 * f( l1-1:l2-1 , m + 3 , n + 2 , ik ))+&
               ( -1170.0 * f( l1-1:l2-1 , m - 1 , n + 3 , ik ))+&
               ( 720.0 * f( l1-1:l2-1 , m - 2 , n + 3 , ik ))+&
               ( -90.0 * f( l1-1:l2-1 , m - 3 , n + 3 , ik ))+&
               ( 1170.0 * f( l1-1:l2-1 , m + 1 , n + 3 , ik ))+&
               ( -720.0 * f( l1-1:l2-1 , m + 2 , n + 3 , ik ))+&
               ( 90.0 * f( l1-1:l2-1 , m + 3 , n + 3 , ik ))+&
               ( -57330.0 * f( l1-2:l2-2 , m - 1 , n , ik ))+&
               ( 35280.0 * f( l1-2:l2-2 , m - 2 , n , ik ))+&
               ( -4410.0 * f( l1-2:l2-2 , m - 3 , n , ik ))+&
               ( 57330.0 * f( l1-2:l2-2 , m + 1 , n , ik ))+&
               ( -35280.0 * f( l1-2:l2-2 , m + 2 , n , ik ))+&
               ( 4410.0 * f( l1-2:l2-2 , m + 3 , n , ik ))+&
               ( 31590.0 * f( l1-2:l2-2 , m - 1 , n - 1 , ik ))+&
               ( -19440.0 * f( l1-2:l2-2 , m - 2 , n - 1 , ik ))+&
               ( 2430.0 * f( l1-2:l2-2 , m - 3 , n - 1 , ik ))+&
               ( -31590.0 * f( l1-2:l2-2 , m + 1 , n - 1 , ik ))+&
               ( 19440.0 * f( l1-2:l2-2 , m + 2 , n - 1 , ik ))+&
               ( -2430.0 * f( l1-2:l2-2 , m + 3 , n - 1 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m - 1 , n - 2 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m - 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1-2:l2-2 , m - 3 , n - 2 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m + 1 , n - 2 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m + 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1-2:l2-2 , m + 3 , n - 2 , ik ))+&
               ( 234.0 * f( l1-2:l2-2 , m - 1 , n - 3 , ik ))+&
               ( -144.0 * f( l1-2:l2-2 , m - 2 , n - 3 , ik ))+&
               ( 18.0 * f( l1-2:l2-2 , m - 3 , n - 3 , ik ))+&
               ( -234.0 * f( l1-2:l2-2 , m + 1 , n - 3 , ik ))+&
               ( 144.0 * f( l1-2:l2-2 , m + 2 , n - 3 , ik ))+&
               ( -18.0 * f( l1-2:l2-2 , m + 3 , n - 3 , ik ))+&
               ( 31590.0 * f( l1-2:l2-2 , m - 1 , n + 1 , ik ))+&
               ( -19440.0 * f( l1-2:l2-2 , m - 2 , n + 1 , ik ))+&
               ( 2430.0 * f( l1-2:l2-2 , m - 3 , n + 1 , ik ))+&
               ( -31590.0 * f( l1-2:l2-2 , m + 1 , n + 1 , ik ))+&
               ( 19440.0 * f( l1-2:l2-2 , m + 2 , n + 1 , ik ))+&
               ( -2430.0 * f( l1-2:l2-2 , m + 3 , n + 1 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m - 1 , n + 2 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m - 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1-2:l2-2 , m - 3 , n + 2 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m + 1 , n + 2 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m + 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1-2:l2-2 , m + 3 , n + 2 , ik ))+&
               ( 234.0 * f( l1-2:l2-2 , m - 1 , n + 3 , ik ))+&
               ( -144.0 * f( l1-2:l2-2 , m - 2 , n + 3 , ik ))+&
               ( 18.0 * f( l1-2:l2-2 , m - 3 , n + 3 , ik ))+&
               ( -234.0 * f( l1-2:l2-2 , m + 1 , n + 3 , ik ))+&
               ( 144.0 * f( l1-2:l2-2 , m + 2 , n + 3 , ik ))+&
               ( -18.0 * f( l1-2:l2-2 , m + 3 , n + 3 , ik ))+&
               ( 6370.0 * f( l1-3:l2-3 , m - 1 , n , ik ))+&
               ( -3920.0 * f( l1-3:l2-3 , m - 2 , n , ik ))+&
               ( 490.0 * f( l1-3:l2-3 , m - 3 , n , ik ))+&
               ( -6370.0 * f( l1-3:l2-3 , m + 1 , n , ik ))+&
               ( 3920.0 * f( l1-3:l2-3 , m + 2 , n , ik ))+&
               ( -490.0 * f( l1-3:l2-3 , m + 3 , n , ik ))+&
               ( -3510.0 * f( l1-3:l2-3 , m - 1 , n - 1 , ik ))+&
               ( 2160.0 * f( l1-3:l2-3 , m - 2 , n - 1 , ik ))+&
               ( -270.0 * f( l1-3:l2-3 , m - 3 , n - 1 , ik ))+&
               ( 3510.0 * f( l1-3:l2-3 , m + 1 , n - 1 , ik ))+&
               ( -2160.0 * f( l1-3:l2-3 , m + 2 , n - 1 , ik ))+&
               ( 270.0 * f( l1-3:l2-3 , m + 3 , n - 1 , ik ))+&
               ( 351.0 * f( l1-3:l2-3 , m - 1 , n - 2 , ik ))+&
               ( -216.0 * f( l1-3:l2-3 , m - 2 , n - 2 , ik ))+&
               ( 27.0 * f( l1-3:l2-3 , m - 3 , n - 2 , ik ))+&
               ( -351.0 * f( l1-3:l2-3 , m + 1 , n - 2 , ik ))+&
               ( 216.0 * f( l1-3:l2-3 , m + 2 , n - 2 , ik ))+&
               ( -27.0 * f( l1-3:l2-3 , m + 3 , n - 2 , ik ))+&
               ( -26.0 * f( l1-3:l2-3 , m - 1 , n - 3 , ik ))+&
               ( 16.0 * f( l1-3:l2-3 , m - 2 , n - 3 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m - 3 , n - 3 , ik ))+&
               ( 26.0 * f( l1-3:l2-3 , m + 1 , n - 3 , ik ))+&
               ( -16.0 * f( l1-3:l2-3 , m + 2 , n - 3 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m + 3 , n - 3 , ik ))+&
               ( -3510.0 * f( l1-3:l2-3 , m - 1 , n + 1 , ik ))+&
               ( 2160.0 * f( l1-3:l2-3 , m - 2 , n + 1 , ik ))+&
               ( -270.0 * f( l1-3:l2-3 , m - 3 , n + 1 , ik ))+&
               ( 3510.0 * f( l1-3:l2-3 , m + 1 , n + 1 , ik ))+&
               ( -2160.0 * f( l1-3:l2-3 , m + 2 , n + 1 , ik ))+&
               ( 270.0 * f( l1-3:l2-3 , m + 3 , n + 1 , ik ))+&
               ( 351.0 * f( l1-3:l2-3 , m - 1 , n + 2 , ik ))+&
               ( -216.0 * f( l1-3:l2-3 , m - 2 , n + 2 , ik ))+&
               ( 27.0 * f( l1-3:l2-3 , m - 3 , n + 2 , ik ))+&
               ( -351.0 * f( l1-3:l2-3 , m + 1 , n + 2 , ik ))+&
               ( 216.0 * f( l1-3:l2-3 , m + 2 , n + 2 , ik ))+&
               ( -27.0 * f( l1-3:l2-3 , m + 3 , n + 2 , ik ))+&
               ( -26.0 * f( l1-3:l2-3 , m - 1 , n + 3 , ik ))+&
               ( 16.0 * f( l1-3:l2-3 , m - 2 , n + 3 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m - 3 , n + 3 , ik ))+&
               ( 26.0 * f( l1-3:l2-3 , m + 1 , n + 3 , ik ))+&
               ( -16.0 * f( l1-3:l2-3 , m + 2 , n + 3 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m + 3 , n + 3 , ik ))+&
               ( -286650.0 * f( l1+1:l2+1 , m - 1 , n , ik ))+&
               ( 176400.0 * f( l1+1:l2+1 , m - 2 , n , ik ))+&
               ( -22050.0 * f( l1+1:l2+1 , m - 3 , n , ik ))+&
               ( 286650.0 * f( l1+1:l2+1 , m + 1 , n , ik ))+&
               ( -176400.0 * f( l1+1:l2+1 , m + 2 , n , ik ))+&
               ( 22050.0 * f( l1+1:l2+1 , m + 3 , n , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m - 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1+1:l2+1 , m - 2 , n - 1 , ik ))+&
               ( 12150.0 * f( l1+1:l2+1 , m - 3 , n - 1 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m + 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1+1:l2+1 , m + 2 , n - 1 , ik ))+&
               ( -12150.0 * f( l1+1:l2+1 , m + 3 , n - 1 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 1 , n - 2 , ik ))+&
               ( 9720.0 * f( l1+1:l2+1 , m - 2 , n - 2 , ik ))+&
               ( -1215.0 * f( l1+1:l2+1 , m - 3 , n - 2 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 1 , n - 2 , ik ))+&
               ( -9720.0 * f( l1+1:l2+1 , m + 2 , n - 2 , ik ))+&
               ( 1215.0 * f( l1+1:l2+1 , m + 3 , n - 2 , ik ))+&
               ( 1170.0 * f( l1+1:l2+1 , m - 1 , n - 3 , ik ))+&
               ( -720.0 * f( l1+1:l2+1 , m - 2 , n - 3 , ik ))+&
               ( 90.0 * f( l1+1:l2+1 , m - 3 , n - 3 , ik ))+&
               ( -1170.0 * f( l1+1:l2+1 , m + 1 , n - 3 , ik ))+&
               ( 720.0 * f( l1+1:l2+1 , m + 2 , n - 3 , ik ))+&
               ( -90.0 * f( l1+1:l2+1 , m + 3 , n - 3 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m - 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1+1:l2+1 , m - 2 , n + 1 , ik ))+&
               ( 12150.0 * f( l1+1:l2+1 , m - 3 , n + 1 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m + 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1+1:l2+1 , m + 2 , n + 1 , ik ))+&
               ( -12150.0 * f( l1+1:l2+1 , m + 3 , n + 1 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 1 , n + 2 , ik ))+&
               ( 9720.0 * f( l1+1:l2+1 , m - 2 , n + 2 , ik ))+&
               ( -1215.0 * f( l1+1:l2+1 , m - 3 , n + 2 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 1 , n + 2 , ik ))+&
               ( -9720.0 * f( l1+1:l2+1 , m + 2 , n + 2 , ik ))+&
               ( 1215.0 * f( l1+1:l2+1 , m + 3 , n + 2 , ik ))+&
               ( 1170.0 * f( l1+1:l2+1 , m - 1 , n + 3 , ik ))+&
               ( -720.0 * f( l1+1:l2+1 , m - 2 , n + 3 , ik ))+&
               ( 90.0 * f( l1+1:l2+1 , m - 3 , n + 3 , ik ))+&
               ( -1170.0 * f( l1+1:l2+1 , m + 1 , n + 3 , ik ))+&
               ( 720.0 * f( l1+1:l2+1 , m + 2 , n + 3 , ik ))+&
               ( -90.0 * f( l1+1:l2+1 , m + 3 , n + 3 , ik ))+&
               ( 57330.0 * f( l1+2:l2+2 , m - 1 , n , ik ))+&
               ( -35280.0 * f( l1+2:l2+2 , m - 2 , n , ik ))+&
               ( 4410.0 * f( l1+2:l2+2 , m - 3 , n , ik ))+&
               ( -57330.0 * f( l1+2:l2+2 , m + 1 , n , ik ))+&
               ( 35280.0 * f( l1+2:l2+2 , m + 2 , n , ik ))+&
               ( -4410.0 * f( l1+2:l2+2 , m + 3 , n , ik ))+&
               ( -31590.0 * f( l1+2:l2+2 , m - 1 , n - 1 , ik ))+&
               ( 19440.0 * f( l1+2:l2+2 , m - 2 , n - 1 , ik ))+&
               ( -2430.0 * f( l1+2:l2+2 , m - 3 , n - 1 , ik ))+&
               ( 31590.0 * f( l1+2:l2+2 , m + 1 , n - 1 , ik ))+&
               ( -19440.0 * f( l1+2:l2+2 , m + 2 , n - 1 , ik ))+&
               ( 2430.0 * f( l1+2:l2+2 , m + 3 , n - 1 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m - 1 , n - 2 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m - 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1+2:l2+2 , m - 3 , n - 2 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m + 1 , n - 2 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m + 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1+2:l2+2 , m + 3 , n - 2 , ik ))+&
               ( -234.0 * f( l1+2:l2+2 , m - 1 , n - 3 , ik ))+&
               ( 144.0 * f( l1+2:l2+2 , m - 2 , n - 3 , ik ))+&
               ( -18.0 * f( l1+2:l2+2 , m - 3 , n - 3 , ik ))+&
               ( 234.0 * f( l1+2:l2+2 , m + 1 , n - 3 , ik ))+&
               ( -144.0 * f( l1+2:l2+2 , m + 2 , n - 3 , ik ))+&
               ( 18.0 * f( l1+2:l2+2 , m + 3 , n - 3 , ik ))+&
               ( -31590.0 * f( l1+2:l2+2 , m - 1 , n + 1 , ik ))+&
               ( 19440.0 * f( l1+2:l2+2 , m - 2 , n + 1 , ik ))+&
               ( -2430.0 * f( l1+2:l2+2 , m - 3 , n + 1 , ik ))+&
               ( 31590.0 * f( l1+2:l2+2 , m + 1 , n + 1 , ik ))+&
               ( -19440.0 * f( l1+2:l2+2 , m + 2 , n + 1 , ik ))+&
               ( 2430.0 * f( l1+2:l2+2 , m + 3 , n + 1 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m - 1 , n + 2 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m - 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1+2:l2+2 , m - 3 , n + 2 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m + 1 , n + 2 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m + 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1+2:l2+2 , m + 3 , n + 2 , ik ))+&
               ( -234.0 * f( l1+2:l2+2 , m - 1 , n + 3 , ik ))+&
               ( 144.0 * f( l1+2:l2+2 , m - 2 , n + 3 , ik ))+&
               ( -18.0 * f( l1+2:l2+2 , m - 3 , n + 3 , ik ))+&
               ( 234.0 * f( l1+2:l2+2 , m + 1 , n + 3 , ik ))+&
               ( -144.0 * f( l1+2:l2+2 , m + 2 , n + 3 , ik ))+&
               ( 18.0 * f( l1+2:l2+2 , m + 3 , n + 3 , ik ))+&
               ( -6370.0 * f( l1+3:l2+3 , m - 1 , n , ik ))+&
               ( 3920.0 * f( l1+3:l2+3 , m - 2 , n , ik ))+&
               ( -490.0 * f( l1+3:l2+3 , m - 3 , n , ik ))+&
               ( 6370.0 * f( l1+3:l2+3 , m + 1 , n , ik ))+&
               ( -3920.0 * f( l1+3:l2+3 , m + 2 , n , ik ))+&
               ( 490.0 * f( l1+3:l2+3 , m + 3 , n , ik ))+&
               ( 3510.0 * f( l1+3:l2+3 , m - 1 , n - 1 , ik ))+&
               ( -2160.0 * f( l1+3:l2+3 , m - 2 , n - 1 , ik ))+&
               ( 270.0 * f( l1+3:l2+3 , m - 3 , n - 1 , ik ))+&
               ( -3510.0 * f( l1+3:l2+3 , m + 1 , n - 1 , ik ))+&
               ( 2160.0 * f( l1+3:l2+3 , m + 2 , n - 1 , ik ))+&
               ( -270.0 * f( l1+3:l2+3 , m + 3 , n - 1 , ik ))+&
               ( -351.0 * f( l1+3:l2+3 , m - 1 , n - 2 , ik ))+&
               ( 216.0 * f( l1+3:l2+3 , m - 2 , n - 2 , ik ))+&
               ( -27.0 * f( l1+3:l2+3 , m - 3 , n - 2 , ik ))+&
               ( 351.0 * f( l1+3:l2+3 , m + 1 , n - 2 , ik ))+&
               ( -216.0 * f( l1+3:l2+3 , m + 2 , n - 2 , ik ))+&
               ( 27.0 * f( l1+3:l2+3 , m + 3 , n - 2 , ik ))+&
               ( 26.0 * f( l1+3:l2+3 , m - 1 , n - 3 , ik ))+&
               ( -16.0 * f( l1+3:l2+3 , m - 2 , n - 3 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m - 3 , n - 3 , ik ))+&
               ( -26.0 * f( l1+3:l2+3 , m + 1 , n - 3 , ik ))+&
               ( 16.0 * f( l1+3:l2+3 , m + 2 , n - 3 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m + 3 , n - 3 , ik ))+&
               ( 3510.0 * f( l1+3:l2+3 , m - 1 , n + 1 , ik ))+&
               ( -2160.0 * f( l1+3:l2+3 , m - 2 , n + 1 , ik ))+&
               ( 270.0 * f( l1+3:l2+3 , m - 3 , n + 1 , ik ))+&
               ( -3510.0 * f( l1+3:l2+3 , m + 1 , n + 1 , ik ))+&
               ( 2160.0 * f( l1+3:l2+3 , m + 2 , n + 1 , ik ))+&
               ( -270.0 * f( l1+3:l2+3 , m + 3 , n + 1 , ik ))+&
               ( -351.0 * f( l1+3:l2+3 , m - 1 , n + 2 , ik ))+&
               ( 216.0 * f( l1+3:l2+3 , m - 2 , n + 2 , ik ))+&
               ( -27.0 * f( l1+3:l2+3 , m - 3 , n + 2 , ik ))+&
               ( 351.0 * f( l1+3:l2+3 , m + 1 , n + 2 , ik ))+&
               ( -216.0 * f( l1+3:l2+3 , m + 2 , n + 2 , ik ))+&
               ( 27.0 * f( l1+3:l2+3 , m + 3 , n + 2 , ik ))+&
               ( 26.0 * f( l1+3:l2+3 , m - 1 , n + 3 , ik ))+&
               ( -16.0 * f( l1+3:l2+3 , m - 2 , n + 3 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m - 3 , n + 3 , ik ))+&
               ( -26.0 * f( l1+3:l2+3 , m + 1 , n + 3 , ik ))+&
               ( 16.0 * f( l1+3:l2+3 , m + 2 , n + 3 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m + 3 , n + 3 , ik ))&
               )
        else
          df=0.0
        endif yzx
!        
      elseif (i==3.and.j==1.and.k==2) then
        zxy: if (nxgrid/=1.and.nygrid/=1.and.nzgrid/=1) then
          fac= (1./8.0*dz_1(n)**3) * (1./180.0*dx_1(l1:l2)**2) * (1/60.0*dy_1(m))
          df = fac*(&
               ( 286650.0 * f( l1:l2 , m - 1 , n - 1 , ik ))+&
               ( -176400.0 * f( l1:l2 , m - 1 , n - 2 , ik ))+&
               ( 22050.0 * f( l1:l2 , m - 1 , n - 3 , ik ))+&
               ( -286650.0 * f( l1:l2 , m - 1 , n + 1 , ik ))+&
               ( 176400.0 * f( l1:l2 , m - 1 , n + 2 , ik ))+&
               ( -22050.0 * f( l1:l2 , m - 1 , n + 3 , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m - 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1-1:l2-1 , m - 1 , n - 2 , ik ))+&
               ( -12150.0 * f( l1-1:l2-1 , m - 1 , n - 3 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m - 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1-1:l2-1 , m - 1 , n + 2 , ik ))+&
               ( 12150.0 * f( l1-1:l2-1 , m - 1 , n + 3 , ik ))+&
               ( 15795.0 * f( l1-2:l2-2 , m - 1 , n - 1 , ik ))+&
               ( -9720.0 * f( l1-2:l2-2 , m - 1 , n - 2 , ik ))+&
               ( 1215.0 * f( l1-2:l2-2 , m - 1 , n - 3 , ik ))+&
               ( -15795.0 * f( l1-2:l2-2 , m - 1 , n + 1 , ik ))+&
               ( 9720.0 * f( l1-2:l2-2 , m - 1 , n + 2 , ik ))+&
               ( -1215.0 * f( l1-2:l2-2 , m - 1 , n + 3 , ik ))+&
               ( -1170.0 * f( l1-3:l2-3 , m - 1 , n - 1 , ik ))+&
               ( 720.0 * f( l1-3:l2-3 , m - 1 , n - 2 , ik ))+&
               ( -90.0 * f( l1-3:l2-3 , m - 1 , n - 3 , ik ))+&
               ( 1170.0 * f( l1-3:l2-3 , m - 1 , n + 1 , ik ))+&
               ( -720.0 * f( l1-3:l2-3 , m - 1 , n + 2 , ik ))+&
               ( 90.0 * f( l1-3:l2-3 , m - 1 , n + 3 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m - 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1+1:l2+1 , m - 1 , n - 2 , ik ))+&
               ( -12150.0 * f( l1+1:l2+1 , m - 1 , n - 3 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m - 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1+1:l2+1 , m - 1 , n + 2 , ik ))+&
               ( 12150.0 * f( l1+1:l2+1 , m - 1 , n + 3 , ik ))+&
               ( 15795.0 * f( l1+2:l2+2 , m - 1 , n - 1 , ik ))+&
               ( -9720.0 * f( l1+2:l2+2 , m - 1 , n - 2 , ik ))+&
               ( 1215.0 * f( l1+2:l2+2 , m - 1 , n - 3 , ik ))+&
               ( -15795.0 * f( l1+2:l2+2 , m - 1 , n + 1 , ik ))+&
               ( 9720.0 * f( l1+2:l2+2 , m - 1 , n + 2 , ik ))+&
               ( -1215.0 * f( l1+2:l2+2 , m - 1 , n + 3 , ik ))+&
               ( -1170.0 * f( l1+3:l2+3 , m - 1 , n - 1 , ik ))+&
               ( 720.0 * f( l1+3:l2+3 , m - 1 , n - 2 , ik ))+&
               ( -90.0 * f( l1+3:l2+3 , m - 1 , n - 3 , ik ))+&
               ( 1170.0 * f( l1+3:l2+3 , m - 1 , n + 1 , ik ))+&
               ( -720.0 * f( l1+3:l2+3 , m - 1 , n + 2 , ik ))+&
               ( 90.0 * f( l1+3:l2+3 , m - 1 , n + 3 , ik ))+&
               ( -57330.0 * f( l1:l2 , m - 2 , n - 1 , ik ))+&
               ( 35280.0 * f( l1:l2 , m - 2 , n - 2 , ik ))+&
               ( -4410.0 * f( l1:l2 , m - 2 , n - 3 , ik ))+&
               ( 57330.0 * f( l1:l2 , m - 2 , n + 1 , ik ))+&
               ( -35280.0 * f( l1:l2 , m - 2 , n + 2 , ik ))+&
               ( 4410.0 * f( l1:l2 , m - 2 , n + 3 , ik ))+&
               ( 31590.0 * f( l1-1:l2-1 , m - 2 , n - 1 , ik ))+&
               ( -19440.0 * f( l1-1:l2-1 , m - 2 , n - 2 , ik ))+&
               ( 2430.0 * f( l1-1:l2-1 , m - 2 , n - 3 , ik ))+&
               ( -31590.0 * f( l1-1:l2-1 , m - 2 , n + 1 , ik ))+&
               ( 19440.0 * f( l1-1:l2-1 , m - 2 , n + 2 , ik ))+&
               ( -2430.0 * f( l1-1:l2-1 , m - 2 , n + 3 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m - 2 , n - 1 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m - 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1-2:l2-2 , m - 2 , n - 3 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m - 2 , n + 1 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m - 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1-2:l2-2 , m - 2 , n + 3 , ik ))+&
               ( 234.0 * f( l1-3:l2-3 , m - 2 , n - 1 , ik ))+&
               ( -144.0 * f( l1-3:l2-3 , m - 2 , n - 2 , ik ))+&
               ( 18.0 * f( l1-3:l2-3 , m - 2 , n - 3 , ik ))+&
               ( -234.0 * f( l1-3:l2-3 , m - 2 , n + 1 , ik ))+&
               ( 144.0 * f( l1-3:l2-3 , m - 2 , n + 2 , ik ))+&
               ( -18.0 * f( l1-3:l2-3 , m - 2 , n + 3 , ik ))+&
               ( 31590.0 * f( l1+1:l2+1 , m - 2 , n - 1 , ik ))+&
               ( -19440.0 * f( l1+1:l2+1 , m - 2 , n - 2 , ik ))+&
               ( 2430.0 * f( l1+1:l2+1 , m - 2 , n - 3 , ik ))+&
               ( -31590.0 * f( l1+1:l2+1 , m - 2 , n + 1 , ik ))+&
               ( 19440.0 * f( l1+1:l2+1 , m - 2 , n + 2 , ik ))+&
               ( -2430.0 * f( l1+1:l2+1 , m - 2 , n + 3 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m - 2 , n - 1 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m - 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1+2:l2+2 , m - 2 , n - 3 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m - 2 , n + 1 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m - 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1+2:l2+2 , m - 2 , n + 3 , ik ))+&
               ( 234.0 * f( l1+3:l2+3 , m - 2 , n - 1 , ik ))+&
               ( -144.0 * f( l1+3:l2+3 , m - 2 , n - 2 , ik ))+&
               ( 18.0 * f( l1+3:l2+3 , m - 2 , n - 3 , ik ))+&
               ( -234.0 * f( l1+3:l2+3 , m - 2 , n + 1 , ik ))+&
               ( 144.0 * f( l1+3:l2+3 , m - 2 , n + 2 , ik ))+&
               ( -18.0 * f( l1+3:l2+3 , m - 2 , n + 3 , ik ))+&
               ( 6370.0 * f( l1:l2 , m - 3 , n - 1 , ik ))+&
               ( -3920.0 * f( l1:l2 , m - 3 , n - 2 , ik ))+&
               ( 490.0 * f( l1:l2 , m - 3 , n - 3 , ik ))+&
               ( -6370.0 * f( l1:l2 , m - 3 , n + 1 , ik ))+&
               ( 3920.0 * f( l1:l2 , m - 3 , n + 2 , ik ))+&
               ( -490.0 * f( l1:l2 , m - 3 , n + 3 , ik ))+&
               ( -3510.0 * f( l1-1:l2-1 , m - 3 , n - 1 , ik ))+&
               ( 2160.0 * f( l1-1:l2-1 , m - 3 , n - 2 , ik ))+&
               ( -270.0 * f( l1-1:l2-1 , m - 3 , n - 3 , ik ))+&
               ( 3510.0 * f( l1-1:l2-1 , m - 3 , n + 1 , ik ))+&
               ( -2160.0 * f( l1-1:l2-1 , m - 3 , n + 2 , ik ))+&
               ( 270.0 * f( l1-1:l2-1 , m - 3 , n + 3 , ik ))+&
               ( 351.0 * f( l1-2:l2-2 , m - 3 , n - 1 , ik ))+&
               ( -216.0 * f( l1-2:l2-2 , m - 3 , n - 2 , ik ))+&
               ( 27.0 * f( l1-2:l2-2 , m - 3 , n - 3 , ik ))+&
               ( -351.0 * f( l1-2:l2-2 , m - 3 , n + 1 , ik ))+&
               ( 216.0 * f( l1-2:l2-2 , m - 3 , n + 2 , ik ))+&
               ( -27.0 * f( l1-2:l2-2 , m - 3 , n + 3 , ik ))+&
               ( -26.0 * f( l1-3:l2-3 , m - 3 , n - 1 , ik ))+&
               ( 16.0 * f( l1-3:l2-3 , m - 3 , n - 2 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m - 3 , n - 3 , ik ))+&
               ( 26.0 * f( l1-3:l2-3 , m - 3 , n + 1 , ik ))+&
               ( -16.0 * f( l1-3:l2-3 , m - 3 , n + 2 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m - 3 , n + 3 , ik ))+&
               ( -3510.0 * f( l1+1:l2+1 , m - 3 , n - 1 , ik ))+&
               ( 2160.0 * f( l1+1:l2+1 , m - 3 , n - 2 , ik ))+&
               ( -270.0 * f( l1+1:l2+1 , m - 3 , n - 3 , ik ))+&
               ( 3510.0 * f( l1+1:l2+1 , m - 3 , n + 1 , ik ))+&
               ( -2160.0 * f( l1+1:l2+1 , m - 3 , n + 2 , ik ))+&
               ( 270.0 * f( l1+1:l2+1 , m - 3 , n + 3 , ik ))+&
               ( 351.0 * f( l1+2:l2+2 , m - 3 , n - 1 , ik ))+&
               ( -216.0 * f( l1+2:l2+2 , m - 3 , n - 2 , ik ))+&
               ( 27.0 * f( l1+2:l2+2 , m - 3 , n - 3 , ik ))+&
               ( -351.0 * f( l1+2:l2+2 , m - 3 , n + 1 , ik ))+&
               ( 216.0 * f( l1+2:l2+2 , m - 3 , n + 2 , ik ))+&
               ( -27.0 * f( l1+2:l2+2 , m - 3 , n + 3 , ik ))+&
               ( -26.0 * f( l1+3:l2+3 , m - 3 , n - 1 , ik ))+&
               ( 16.0 * f( l1+3:l2+3 , m - 3 , n - 2 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m - 3 , n - 3 , ik ))+&
               ( 26.0 * f( l1+3:l2+3 , m - 3 , n + 1 , ik ))+&
               ( -16.0 * f( l1+3:l2+3 , m - 3 , n + 2 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m - 3 , n + 3 , ik ))+&
               ( -286650.0 * f( l1:l2 , m + 1 , n - 1 , ik ))+&
               ( 176400.0 * f( l1:l2 , m + 1 , n - 2 , ik ))+&
               ( -22050.0 * f( l1:l2 , m + 1 , n - 3 , ik ))+&
               ( 286650.0 * f( l1:l2 , m + 1 , n + 1 , ik ))+&
               ( -176400.0 * f( l1:l2 , m + 1 , n + 2 , ik ))+&
               ( 22050.0 * f( l1:l2 , m + 1 , n + 3 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m + 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1-1:l2-1 , m + 1 , n - 2 , ik ))+&
               ( 12150.0 * f( l1-1:l2-1 , m + 1 , n - 3 , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m + 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1-1:l2-1 , m + 1 , n + 2 , ik ))+&
               ( -12150.0 * f( l1-1:l2-1 , m + 1 , n + 3 , ik ))+&
               ( -15795.0 * f( l1-2:l2-2 , m + 1 , n - 1 , ik ))+&
               ( 9720.0 * f( l1-2:l2-2 , m + 1 , n - 2 , ik ))+&
               ( -1215.0 * f( l1-2:l2-2 , m + 1 , n - 3 , ik ))+&
               ( 15795.0 * f( l1-2:l2-2 , m + 1 , n + 1 , ik ))+&
               ( -9720.0 * f( l1-2:l2-2 , m + 1 , n + 2 , ik ))+&
               ( 1215.0 * f( l1-2:l2-2 , m + 1 , n + 3 , ik ))+&
               ( 1170.0 * f( l1-3:l2-3 , m + 1 , n - 1 , ik ))+&
               ( -720.0 * f( l1-3:l2-3 , m + 1 , n - 2 , ik ))+&
               ( 90.0 * f( l1-3:l2-3 , m + 1 , n - 3 , ik ))+&
               ( -1170.0 * f( l1-3:l2-3 , m + 1 , n + 1 , ik ))+&
               ( 720.0 * f( l1-3:l2-3 , m + 1 , n + 2 , ik ))+&
               ( -90.0 * f( l1-3:l2-3 , m + 1 , n + 3 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m + 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1+1:l2+1 , m + 1 , n - 2 , ik ))+&
               ( 12150.0 * f( l1+1:l2+1 , m + 1 , n - 3 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m + 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1+1:l2+1 , m + 1 , n + 2 , ik ))+&
               ( -12150.0 * f( l1+1:l2+1 , m + 1 , n + 3 , ik ))+&
               ( -15795.0 * f( l1+2:l2+2 , m + 1 , n - 1 , ik ))+&
               ( 9720.0 * f( l1+2:l2+2 , m + 1 , n - 2 , ik ))+&
               ( -1215.0 * f( l1+2:l2+2 , m + 1 , n - 3 , ik ))+&
               ( 15795.0 * f( l1+2:l2+2 , m + 1 , n + 1 , ik ))+&
               ( -9720.0 * f( l1+2:l2+2 , m + 1 , n + 2 , ik ))+&
               ( 1215.0 * f( l1+2:l2+2 , m + 1 , n + 3 , ik ))+&
               ( 1170.0 * f( l1+3:l2+3 , m + 1 , n - 1 , ik ))+&
               ( -720.0 * f( l1+3:l2+3 , m + 1 , n - 2 , ik ))+&
               ( 90.0 * f( l1+3:l2+3 , m + 1 , n - 3 , ik ))+&
               ( -1170.0 * f( l1+3:l2+3 , m + 1 , n + 1 , ik ))+&
               ( 720.0 * f( l1+3:l2+3 , m + 1 , n + 2 , ik ))+&
               ( -90.0 * f( l1+3:l2+3 , m + 1 , n + 3 , ik ))+&
               ( 57330.0 * f( l1:l2 , m + 2 , n - 1 , ik ))+&
               ( -35280.0 * f( l1:l2 , m + 2 , n - 2 , ik ))+&
               ( 4410.0 * f( l1:l2 , m + 2 , n - 3 , ik ))+&
               ( -57330.0 * f( l1:l2 , m + 2 , n + 1 , ik ))+&
               ( 35280.0 * f( l1:l2 , m + 2 , n + 2 , ik ))+&
               ( -4410.0 * f( l1:l2 , m + 2 , n + 3 , ik ))+&
               ( -31590.0 * f( l1-1:l2-1 , m + 2 , n - 1 , ik ))+&
               ( 19440.0 * f( l1-1:l2-1 , m + 2 , n - 2 , ik ))+&
               ( -2430.0 * f( l1-1:l2-1 , m + 2 , n - 3 , ik ))+&
               ( 31590.0 * f( l1-1:l2-1 , m + 2 , n + 1 , ik ))+&
               ( -19440.0 * f( l1-1:l2-1 , m + 2 , n + 2 , ik ))+&
               ( 2430.0 * f( l1-1:l2-1 , m + 2 , n + 3 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m + 2 , n - 1 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m + 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1-2:l2-2 , m + 2 , n - 3 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m + 2 , n + 1 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m + 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1-2:l2-2 , m + 2 , n + 3 , ik ))+&
               ( -234.0 * f( l1-3:l2-3 , m + 2 , n - 1 , ik ))+&
               ( 144.0 * f( l1-3:l2-3 , m + 2 , n - 2 , ik ))+&
               ( -18.0 * f( l1-3:l2-3 , m + 2 , n - 3 , ik ))+&
               ( 234.0 * f( l1-3:l2-3 , m + 2 , n + 1 , ik ))+&
               ( -144.0 * f( l1-3:l2-3 , m + 2 , n + 2 , ik ))+&
               ( 18.0 * f( l1-3:l2-3 , m + 2 , n + 3 , ik ))+&
               ( -31590.0 * f( l1+1:l2+1 , m + 2 , n - 1 , ik ))+&
               ( 19440.0 * f( l1+1:l2+1 , m + 2 , n - 2 , ik ))+&
               ( -2430.0 * f( l1+1:l2+1 , m + 2 , n - 3 , ik ))+&
               ( 31590.0 * f( l1+1:l2+1 , m + 2 , n + 1 , ik ))+&
               ( -19440.0 * f( l1+1:l2+1 , m + 2 , n + 2 , ik ))+&
               ( 2430.0 * f( l1+1:l2+1 , m + 2 , n + 3 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m + 2 , n - 1 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m + 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1+2:l2+2 , m + 2 , n - 3 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m + 2 , n + 1 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m + 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1+2:l2+2 , m + 2 , n + 3 , ik ))+&
               ( -234.0 * f( l1+3:l2+3 , m + 2 , n - 1 , ik ))+&
               ( 144.0 * f( l1+3:l2+3 , m + 2 , n - 2 , ik ))+&
               ( -18.0 * f( l1+3:l2+3 , m + 2 , n - 3 , ik ))+&
               ( 234.0 * f( l1+3:l2+3 , m + 2 , n + 1 , ik ))+&
               ( -144.0 * f( l1+3:l2+3 , m + 2 , n + 2 , ik ))+&
               ( 18.0 * f( l1+3:l2+3 , m + 2 , n + 3 , ik ))+&
               ( -6370.0 * f( l1:l2 , m + 3 , n - 1 , ik ))+&
               ( 3920.0 * f( l1:l2 , m + 3 , n - 2 , ik ))+&
               ( -490.0 * f( l1:l2 , m + 3 , n - 3 , ik ))+&
               ( 6370.0 * f( l1:l2 , m + 3 , n + 1 , ik ))+&
               ( -3920.0 * f( l1:l2 , m + 3 , n + 2 , ik ))+&
               ( 490.0 * f( l1:l2 , m + 3 , n + 3 , ik ))+&
               ( 3510.0 * f( l1-1:l2-1 , m + 3 , n - 1 , ik ))+&
               ( -2160.0 * f( l1-1:l2-1 , m + 3 , n - 2 , ik ))+&
               ( 270.0 * f( l1-1:l2-1 , m + 3 , n - 3 , ik ))+&
               ( -3510.0 * f( l1-1:l2-1 , m + 3 , n + 1 , ik ))+&
               ( 2160.0 * f( l1-1:l2-1 , m + 3 , n + 2 , ik ))+&
               ( -270.0 * f( l1-1:l2-1 , m + 3 , n + 3 , ik ))+&
               ( -351.0 * f( l1-2:l2-2 , m + 3 , n - 1 , ik ))+&
               ( 216.0 * f( l1-2:l2-2 , m + 3 , n - 2 , ik ))+&
               ( -27.0 * f( l1-2:l2-2 , m + 3 , n - 3 , ik ))+&
               ( 351.0 * f( l1-2:l2-2 , m + 3 , n + 1 , ik ))+&
               ( -216.0 * f( l1-2:l2-2 , m + 3 , n + 2 , ik ))+&
               ( 27.0 * f( l1-2:l2-2 , m + 3 , n + 3 , ik ))+&
               ( 26.0 * f( l1-3:l2-3 , m + 3 , n - 1 , ik ))+&
               ( -16.0 * f( l1-3:l2-3 , m + 3 , n - 2 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m + 3 , n - 3 , ik ))+&
               ( -26.0 * f( l1-3:l2-3 , m + 3 , n + 1 , ik ))+&
               ( 16.0 * f( l1-3:l2-3 , m + 3 , n + 2 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m + 3 , n + 3 , ik ))+&
               ( 3510.0 * f( l1+1:l2+1 , m + 3 , n - 1 , ik ))+&
               ( -2160.0 * f( l1+1:l2+1 , m + 3 , n - 2 , ik ))+&
               ( 270.0 * f( l1+1:l2+1 , m + 3 , n - 3 , ik ))+&
               ( -3510.0 * f( l1+1:l2+1 , m + 3 , n + 1 , ik ))+&
               ( 2160.0 * f( l1+1:l2+1 , m + 3 , n + 2 , ik ))+&
               ( -270.0 * f( l1+1:l2+1 , m + 3 , n + 3 , ik ))+&
               ( -351.0 * f( l1+2:l2+2 , m + 3 , n - 1 , ik ))+&
               ( 216.0 * f( l1+2:l2+2 , m + 3 , n - 2 , ik ))+&
               ( -27.0 * f( l1+2:l2+2 , m + 3 , n - 3 , ik ))+&
               ( 351.0 * f( l1+2:l2+2 , m + 3 , n + 1 , ik ))+&
               ( -216.0 * f( l1+2:l2+2 , m + 3 , n + 2 , ik ))+&
               ( 27.0 * f( l1+2:l2+2 , m + 3 , n + 3 , ik ))+&
               ( 26.0 * f( l1+3:l2+3 , m + 3 , n - 1 , ik ))+&
               ( -16.0 * f( l1+3:l2+3 , m + 3 , n - 2 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m + 3 , n - 3 , ik ))+&
               ( -26.0 * f( l1+3:l2+3 , m + 3 , n + 1 , ik ))+&
               ( 16.0 * f( l1+3:l2+3 , m + 3 , n + 2 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m + 3 , n + 3 , ik ))&
               )
        else
          df=0.0
        endif zxy
!        
      elseif (i==3.and.j==2.and.k==1) then
        zyx: if (nxgrid/=1.and.nygrid/=1.and.nzgrid/=1) then
          fac= (1./8.0*dz_1(n)**3) * (1./180.0*dy_1(m)**2) * (1/60.0*dx_1(l1:l2))
          df = fac*(&
               ( 286650.0 * f( l1-1:l2-1 , m , n - 1 , ik ))+&
               ( -176400.0 * f( l1-1:l2-1 , m , n - 2 , ik ))+&
               ( 22050.0 * f( l1-1:l2-1 , m , n - 3 , ik ))+&
               ( -286650.0 * f( l1-1:l2-1 , m , n + 1 , ik ))+&
               ( 176400.0 * f( l1-1:l2-1 , m , n + 2 , ik ))+&
               ( -22050.0 * f( l1-1:l2-1 , m , n + 3 , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m - 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1-1:l2-1 , m - 1 , n - 2 , ik ))+&
               ( -12150.0 * f( l1-1:l2-1 , m - 1 , n - 3 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m - 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1-1:l2-1 , m - 1 , n + 2 , ik ))+&
               ( 12150.0 * f( l1-1:l2-1 , m - 1 , n + 3 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 2 , n - 1 , ik ))+&
               ( -9720.0 * f( l1-1:l2-1 , m - 2 , n - 2 , ik ))+&
               ( 1215.0 * f( l1-1:l2-1 , m - 2 , n - 3 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m - 2 , n + 1 , ik ))+&
               ( 9720.0 * f( l1-1:l2-1 , m - 2 , n + 2 , ik ))+&
               ( -1215.0 * f( l1-1:l2-1 , m - 2 , n + 3 , ik ))+&
               ( -1170.0 * f( l1-1:l2-1 , m - 3 , n - 1 , ik ))+&
               ( 720.0 * f( l1-1:l2-1 , m - 3 , n - 2 , ik ))+&
               ( -90.0 * f( l1-1:l2-1 , m - 3 , n - 3 , ik ))+&
               ( 1170.0 * f( l1-1:l2-1 , m - 3 , n + 1 , ik ))+&
               ( -720.0 * f( l1-1:l2-1 , m - 3 , n + 2 , ik ))+&
               ( 90.0 * f( l1-1:l2-1 , m - 3 , n + 3 , ik ))+&
               ( -157950.0 * f( l1-1:l2-1 , m + 1 , n - 1 , ik ))+&
               ( 97200.0 * f( l1-1:l2-1 , m + 1 , n - 2 , ik ))+&
               ( -12150.0 * f( l1-1:l2-1 , m + 1 , n - 3 , ik ))+&
               ( 157950.0 * f( l1-1:l2-1 , m + 1 , n + 1 , ik ))+&
               ( -97200.0 * f( l1-1:l2-1 , m + 1 , n + 2 , ik ))+&
               ( 12150.0 * f( l1-1:l2-1 , m + 1 , n + 3 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m + 2 , n - 1 , ik ))+&
               ( -9720.0 * f( l1-1:l2-1 , m + 2 , n - 2 , ik ))+&
               ( 1215.0 * f( l1-1:l2-1 , m + 2 , n - 3 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 2 , n + 1 , ik ))+&
               ( 9720.0 * f( l1-1:l2-1 , m + 2 , n + 2 , ik ))+&
               ( -1215.0 * f( l1-1:l2-1 , m + 2 , n + 3 , ik ))+&
               ( -1170.0 * f( l1-1:l2-1 , m + 3 , n - 1 , ik ))+&
               ( 720.0 * f( l1-1:l2-1 , m + 3 , n - 2 , ik ))+&
               ( -90.0 * f( l1-1:l2-1 , m + 3 , n - 3 , ik ))+&
               ( 1170.0 * f( l1-1:l2-1 , m + 3 , n + 1 , ik ))+&
               ( -720.0 * f( l1-1:l2-1 , m + 3 , n + 2 , ik ))+&
               ( 90.0 * f( l1-1:l2-1 , m + 3 , n + 3 , ik ))+&
               ( -57330.0 * f( l1-2:l2-2 , m , n - 1 , ik ))+&
               ( 35280.0 * f( l1-2:l2-2 , m , n - 2 , ik ))+&
               ( -4410.0 * f( l1-2:l2-2 , m , n - 3 , ik ))+&
               ( 57330.0 * f( l1-2:l2-2 , m , n + 1 , ik ))+&
               ( -35280.0 * f( l1-2:l2-2 , m , n + 2 , ik ))+&
               ( 4410.0 * f( l1-2:l2-2 , m , n + 3 , ik ))+&
               ( 31590.0 * f( l1-2:l2-2 , m - 1 , n - 1 , ik ))+&
               ( -19440.0 * f( l1-2:l2-2 , m - 1 , n - 2 , ik ))+&
               ( 2430.0 * f( l1-2:l2-2 , m - 1 , n - 3 , ik ))+&
               ( -31590.0 * f( l1-2:l2-2 , m - 1 , n + 1 , ik ))+&
               ( 19440.0 * f( l1-2:l2-2 , m - 1 , n + 2 , ik ))+&
               ( -2430.0 * f( l1-2:l2-2 , m - 1 , n + 3 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m - 2 , n - 1 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m - 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1-2:l2-2 , m - 2 , n - 3 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m - 2 , n + 1 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m - 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1-2:l2-2 , m - 2 , n + 3 , ik ))+&
               ( 234.0 * f( l1-2:l2-2 , m - 3 , n - 1 , ik ))+&
               ( -144.0 * f( l1-2:l2-2 , m - 3 , n - 2 , ik ))+&
               ( 18.0 * f( l1-2:l2-2 , m - 3 , n - 3 , ik ))+&
               ( -234.0 * f( l1-2:l2-2 , m - 3 , n + 1 , ik ))+&
               ( 144.0 * f( l1-2:l2-2 , m - 3 , n + 2 , ik ))+&
               ( -18.0 * f( l1-2:l2-2 , m - 3 , n + 3 , ik ))+&
               ( 31590.0 * f( l1-2:l2-2 , m + 1 , n - 1 , ik ))+&
               ( -19440.0 * f( l1-2:l2-2 , m + 1 , n - 2 , ik ))+&
               ( 2430.0 * f( l1-2:l2-2 , m + 1 , n - 3 , ik ))+&
               ( -31590.0 * f( l1-2:l2-2 , m + 1 , n + 1 , ik ))+&
               ( 19440.0 * f( l1-2:l2-2 , m + 1 , n + 2 , ik ))+&
               ( -2430.0 * f( l1-2:l2-2 , m + 1 , n + 3 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m + 2 , n - 1 , ik ))+&
               ( 1944.0 * f( l1-2:l2-2 , m + 2 , n - 2 , ik ))+&
               ( -243.0 * f( l1-2:l2-2 , m + 2 , n - 3 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m + 2 , n + 1 , ik ))+&
               ( -1944.0 * f( l1-2:l2-2 , m + 2 , n + 2 , ik ))+&
               ( 243.0 * f( l1-2:l2-2 , m + 2 , n + 3 , ik ))+&
               ( 234.0 * f( l1-2:l2-2 , m + 3 , n - 1 , ik ))+&
               ( -144.0 * f( l1-2:l2-2 , m + 3 , n - 2 , ik ))+&
               ( 18.0 * f( l1-2:l2-2 , m + 3 , n - 3 , ik ))+&
               ( -234.0 * f( l1-2:l2-2 , m + 3 , n + 1 , ik ))+&
               ( 144.0 * f( l1-2:l2-2 , m + 3 , n + 2 , ik ))+&
               ( -18.0 * f( l1-2:l2-2 , m + 3 , n + 3 , ik ))+&
               ( 6370.0 * f( l1-3:l2-3 , m , n - 1 , ik ))+&
               ( -3920.0 * f( l1-3:l2-3 , m , n - 2 , ik ))+&
               ( 490.0 * f( l1-3:l2-3 , m , n - 3 , ik ))+&
               ( -6370.0 * f( l1-3:l2-3 , m , n + 1 , ik ))+&
               ( 3920.0 * f( l1-3:l2-3 , m , n + 2 , ik ))+&
               ( -490.0 * f( l1-3:l2-3 , m , n + 3 , ik ))+&
               ( -3510.0 * f( l1-3:l2-3 , m - 1 , n - 1 , ik ))+&
               ( 2160.0 * f( l1-3:l2-3 , m - 1 , n - 2 , ik ))+&
               ( -270.0 * f( l1-3:l2-3 , m - 1 , n - 3 , ik ))+&
               ( 3510.0 * f( l1-3:l2-3 , m - 1 , n + 1 , ik ))+&
               ( -2160.0 * f( l1-3:l2-3 , m - 1 , n + 2 , ik ))+&
               ( 270.0 * f( l1-3:l2-3 , m - 1 , n + 3 , ik ))+&
               ( 351.0 * f( l1-3:l2-3 , m - 2 , n - 1 , ik ))+&
               ( -216.0 * f( l1-3:l2-3 , m - 2 , n - 2 , ik ))+&
               ( 27.0 * f( l1-3:l2-3 , m - 2 , n - 3 , ik ))+&
               ( -351.0 * f( l1-3:l2-3 , m - 2 , n + 1 , ik ))+&
               ( 216.0 * f( l1-3:l2-3 , m - 2 , n + 2 , ik ))+&
               ( -27.0 * f( l1-3:l2-3 , m - 2 , n + 3 , ik ))+&
               ( -26.0 * f( l1-3:l2-3 , m - 3 , n - 1 , ik ))+&
               ( 16.0 * f( l1-3:l2-3 , m - 3 , n - 2 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m - 3 , n - 3 , ik ))+&
               ( 26.0 * f( l1-3:l2-3 , m - 3 , n + 1 , ik ))+&
               ( -16.0 * f( l1-3:l2-3 , m - 3 , n + 2 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m - 3 , n + 3 , ik ))+&
               ( -3510.0 * f( l1-3:l2-3 , m + 1 , n - 1 , ik ))+&
               ( 2160.0 * f( l1-3:l2-3 , m + 1 , n - 2 , ik ))+&
               ( -270.0 * f( l1-3:l2-3 , m + 1 , n - 3 , ik ))+&
               ( 3510.0 * f( l1-3:l2-3 , m + 1 , n + 1 , ik ))+&
               ( -2160.0 * f( l1-3:l2-3 , m + 1 , n + 2 , ik ))+&
               ( 270.0 * f( l1-3:l2-3 , m + 1 , n + 3 , ik ))+&
               ( 351.0 * f( l1-3:l2-3 , m + 2 , n - 1 , ik ))+&
               ( -216.0 * f( l1-3:l2-3 , m + 2 , n - 2 , ik ))+&
               ( 27.0 * f( l1-3:l2-3 , m + 2 , n - 3 , ik ))+&
               ( -351.0 * f( l1-3:l2-3 , m + 2 , n + 1 , ik ))+&
               ( 216.0 * f( l1-3:l2-3 , m + 2 , n + 2 , ik ))+&
               ( -27.0 * f( l1-3:l2-3 , m + 2 , n + 3 , ik ))+&
               ( -26.0 * f( l1-3:l2-3 , m + 3 , n - 1 , ik ))+&
               ( 16.0 * f( l1-3:l2-3 , m + 3 , n - 2 , ik ))+&
               ( -2.0 * f( l1-3:l2-3 , m + 3 , n - 3 , ik ))+&
               ( 26.0 * f( l1-3:l2-3 , m + 3 , n + 1 , ik ))+&
               ( -16.0 * f( l1-3:l2-3 , m + 3 , n + 2 , ik ))+&
               ( 2.0 * f( l1-3:l2-3 , m + 3 , n + 3 , ik ))+&
               ( -286650.0 * f( l1+1:l2+1 , m , n - 1 , ik ))+&
               ( 176400.0 * f( l1+1:l2+1 , m , n - 2 , ik ))+&
               ( -22050.0 * f( l1+1:l2+1 , m , n - 3 , ik ))+&
               ( 286650.0 * f( l1+1:l2+1 , m , n + 1 , ik ))+&
               ( -176400.0 * f( l1+1:l2+1 , m , n + 2 , ik ))+&
               ( 22050.0 * f( l1+1:l2+1 , m , n + 3 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m - 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1+1:l2+1 , m - 1 , n - 2 , ik ))+&
               ( 12150.0 * f( l1+1:l2+1 , m - 1 , n - 3 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m - 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1+1:l2+1 , m - 1 , n + 2 , ik ))+&
               ( -12150.0 * f( l1+1:l2+1 , m - 1 , n + 3 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 2 , n - 1 , ik ))+&
               ( 9720.0 * f( l1+1:l2+1 , m - 2 , n - 2 , ik ))+&
               ( -1215.0 * f( l1+1:l2+1 , m - 2 , n - 3 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m - 2 , n + 1 , ik ))+&
               ( -9720.0 * f( l1+1:l2+1 , m - 2 , n + 2 , ik ))+&
               ( 1215.0 * f( l1+1:l2+1 , m - 2 , n + 3 , ik ))+&
               ( 1170.0 * f( l1+1:l2+1 , m - 3 , n - 1 , ik ))+&
               ( -720.0 * f( l1+1:l2+1 , m - 3 , n - 2 , ik ))+&
               ( 90.0 * f( l1+1:l2+1 , m - 3 , n - 3 , ik ))+&
               ( -1170.0 * f( l1+1:l2+1 , m - 3 , n + 1 , ik ))+&
               ( 720.0 * f( l1+1:l2+1 , m - 3 , n + 2 , ik ))+&
               ( -90.0 * f( l1+1:l2+1 , m - 3 , n + 3 , ik ))+&
               ( 157950.0 * f( l1+1:l2+1 , m + 1 , n - 1 , ik ))+&
               ( -97200.0 * f( l1+1:l2+1 , m + 1 , n - 2 , ik ))+&
               ( 12150.0 * f( l1+1:l2+1 , m + 1 , n - 3 , ik ))+&
               ( -157950.0 * f( l1+1:l2+1 , m + 1 , n + 1 , ik ))+&
               ( 97200.0 * f( l1+1:l2+1 , m + 1 , n + 2 , ik ))+&
               ( -12150.0 * f( l1+1:l2+1 , m + 1 , n + 3 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m + 2 , n - 1 , ik ))+&
               ( 9720.0 * f( l1+1:l2+1 , m + 2 , n - 2 , ik ))+&
               ( -1215.0 * f( l1+1:l2+1 , m + 2 , n - 3 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 2 , n + 1 , ik ))+&
               ( -9720.0 * f( l1+1:l2+1 , m + 2 , n + 2 , ik ))+&
               ( 1215.0 * f( l1+1:l2+1 , m + 2 , n + 3 , ik ))+&
               ( 1170.0 * f( l1+1:l2+1 , m + 3 , n - 1 , ik ))+&
               ( -720.0 * f( l1+1:l2+1 , m + 3 , n - 2 , ik ))+&
               ( 90.0 * f( l1+1:l2+1 , m + 3 , n - 3 , ik ))+&
               ( -1170.0 * f( l1+1:l2+1 , m + 3 , n + 1 , ik ))+&
               ( 720.0 * f( l1+1:l2+1 , m + 3 , n + 2 , ik ))+&
               ( -90.0 * f( l1+1:l2+1 , m + 3 , n + 3 , ik ))+&
               ( 57330.0 * f( l1+2:l2+2 , m , n - 1 , ik ))+&
               ( -35280.0 * f( l1+2:l2+2 , m , n - 2 , ik ))+&
               ( 4410.0 * f( l1+2:l2+2 , m , n - 3 , ik ))+&
               ( -57330.0 * f( l1+2:l2+2 , m , n + 1 , ik ))+&
               ( 35280.0 * f( l1+2:l2+2 , m , n + 2 , ik ))+&
               ( -4410.0 * f( l1+2:l2+2 , m , n + 3 , ik ))+&
               ( -31590.0 * f( l1+2:l2+2 , m - 1 , n - 1 , ik ))+&
               ( 19440.0 * f( l1+2:l2+2 , m - 1 , n - 2 , ik ))+&
               ( -2430.0 * f( l1+2:l2+2 , m - 1 , n - 3 , ik ))+&
               ( 31590.0 * f( l1+2:l2+2 , m - 1 , n + 1 , ik ))+&
               ( -19440.0 * f( l1+2:l2+2 , m - 1 , n + 2 , ik ))+&
               ( 2430.0 * f( l1+2:l2+2 , m - 1 , n + 3 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m - 2 , n - 1 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m - 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1+2:l2+2 , m - 2 , n - 3 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m - 2 , n + 1 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m - 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1+2:l2+2 , m - 2 , n + 3 , ik ))+&
               ( -234.0 * f( l1+2:l2+2 , m - 3 , n - 1 , ik ))+&
               ( 144.0 * f( l1+2:l2+2 , m - 3 , n - 2 , ik ))+&
               ( -18.0 * f( l1+2:l2+2 , m - 3 , n - 3 , ik ))+&
               ( 234.0 * f( l1+2:l2+2 , m - 3 , n + 1 , ik ))+&
               ( -144.0 * f( l1+2:l2+2 , m - 3 , n + 2 , ik ))+&
               ( 18.0 * f( l1+2:l2+2 , m - 3 , n + 3 , ik ))+&
               ( -31590.0 * f( l1+2:l2+2 , m + 1 , n - 1 , ik ))+&
               ( 19440.0 * f( l1+2:l2+2 , m + 1 , n - 2 , ik ))+&
               ( -2430.0 * f( l1+2:l2+2 , m + 1 , n - 3 , ik ))+&
               ( 31590.0 * f( l1+2:l2+2 , m + 1 , n + 1 , ik ))+&
               ( -19440.0 * f( l1+2:l2+2 , m + 1 , n + 2 , ik ))+&
               ( 2430.0 * f( l1+2:l2+2 , m + 1 , n + 3 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m + 2 , n - 1 , ik ))+&
               ( -1944.0 * f( l1+2:l2+2 , m + 2 , n - 2 , ik ))+&
               ( 243.0 * f( l1+2:l2+2 , m + 2 , n - 3 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m + 2 , n + 1 , ik ))+&
               ( 1944.0 * f( l1+2:l2+2 , m + 2 , n + 2 , ik ))+&
               ( -243.0 * f( l1+2:l2+2 , m + 2 , n + 3 , ik ))+&
               ( -234.0 * f( l1+2:l2+2 , m + 3 , n - 1 , ik ))+&
               ( 144.0 * f( l1+2:l2+2 , m + 3 , n - 2 , ik ))+&
               ( -18.0 * f( l1+2:l2+2 , m + 3 , n - 3 , ik ))+&
               ( 234.0 * f( l1+2:l2+2 , m + 3 , n + 1 , ik ))+&
               ( -144.0 * f( l1+2:l2+2 , m + 3 , n + 2 , ik ))+&
               ( 18.0 * f( l1+2:l2+2 , m + 3 , n + 3 , ik ))+&
               ( -6370.0 * f( l1+3:l2+3 , m , n - 1 , ik ))+&
               ( 3920.0 * f( l1+3:l2+3 , m , n - 2 , ik ))+&
               ( -490.0 * f( l1+3:l2+3 , m , n - 3 , ik ))+&
               ( 6370.0 * f( l1+3:l2+3 , m , n + 1 , ik ))+&
               ( -3920.0 * f( l1+3:l2+3 , m , n + 2 , ik ))+&
               ( 490.0 * f( l1+3:l2+3 , m , n + 3 , ik ))+&
               ( 3510.0 * f( l1+3:l2+3 , m - 1 , n - 1 , ik ))+&
               ( -2160.0 * f( l1+3:l2+3 , m - 1 , n - 2 , ik ))+&
               ( 270.0 * f( l1+3:l2+3 , m - 1 , n - 3 , ik ))+&
               ( -3510.0 * f( l1+3:l2+3 , m - 1 , n + 1 , ik ))+&
               ( 2160.0 * f( l1+3:l2+3 , m - 1 , n + 2 , ik ))+&
               ( -270.0 * f( l1+3:l2+3 , m - 1 , n + 3 , ik ))+&
               ( -351.0 * f( l1+3:l2+3 , m - 2 , n - 1 , ik ))+&
               ( 216.0 * f( l1+3:l2+3 , m - 2 , n - 2 , ik ))+&
               ( -27.0 * f( l1+3:l2+3 , m - 2 , n - 3 , ik ))+&
               ( 351.0 * f( l1+3:l2+3 , m - 2 , n + 1 , ik ))+&
               ( -216.0 * f( l1+3:l2+3 , m - 2 , n + 2 , ik ))+&
               ( 27.0 * f( l1+3:l2+3 , m - 2 , n + 3 , ik ))+&
               ( 26.0 * f( l1+3:l2+3 , m - 3 , n - 1 , ik ))+&
               ( -16.0 * f( l1+3:l2+3 , m - 3 , n - 2 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m - 3 , n - 3 , ik ))+&
               ( -26.0 * f( l1+3:l2+3 , m - 3 , n + 1 , ik ))+&
               ( 16.0 * f( l1+3:l2+3 , m - 3 , n + 2 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m - 3 , n + 3 , ik ))+&
               ( 3510.0 * f( l1+3:l2+3 , m + 1 , n - 1 , ik ))+&
               ( -2160.0 * f( l1+3:l2+3 , m + 1 , n - 2 , ik ))+&
               ( 270.0 * f( l1+3:l2+3 , m + 1 , n - 3 , ik ))+&
               ( -3510.0 * f( l1+3:l2+3 , m + 1 , n + 1 , ik ))+&
               ( 2160.0 * f( l1+3:l2+3 , m + 1 , n + 2 , ik ))+&
               ( -270.0 * f( l1+3:l2+3 , m + 1 , n + 3 , ik ))+&
               ( -351.0 * f( l1+3:l2+3 , m + 2 , n - 1 , ik ))+&
               ( 216.0 * f( l1+3:l2+3 , m + 2 , n - 2 , ik ))+&
               ( -27.0 * f( l1+3:l2+3 , m + 2 , n - 3 , ik ))+&
               ( 351.0 * f( l1+3:l2+3 , m + 2 , n + 1 , ik ))+&
               ( -216.0 * f( l1+3:l2+3 , m + 2 , n + 2 , ik ))+&
               ( 27.0 * f( l1+3:l2+3 , m + 2 , n + 3 , ik ))+&
               ( 26.0 * f( l1+3:l2+3 , m + 3 , n - 1 , ik ))+&
               ( -16.0 * f( l1+3:l2+3 , m + 3 , n - 2 , ik ))+&
               ( 2.0 * f( l1+3:l2+3 , m + 3 , n - 3 , ik ))+&
               ( -26.0 * f( l1+3:l2+3 , m + 3 , n + 1 , ik ))+&
               ( 16.0 * f( l1+3:l2+3 , m + 3 , n + 2 , ik ))+&
               ( -2.0 * f( l1+3:l2+3 , m + 3 , n + 3 , ik ))&
               )
        else
          df=0.0
        endif zyx
!        
      endif
!
    endsubroutine der3i2j1k
!***********************************************************************
    subroutine der4i1j1k(f,ik,df,i,j,k)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx), intent(out) :: df
      real, dimension (nx) :: fac
      integer, intent(in) :: ik,i,j,k
!
      if (headtt) then
        if (lspherical_coords.or.lcylindrical_coords) &
             call warning('der4i1j1k','MISSING CURVATURE TERMS for non-cartesian coordinates')
      endif
!
      if (i==j.and.j==k) then
        if (i==1.and.nxgrid==1 .or. i==2.and.nygrid==1 .or. i==3.and.nzgrid==1) then
          df=0.
          if (ip<=5) print*, 'der4i1j1k: Degenerate case in '//coornames(i)//'-direction'                             
        else
          call der6(f,ik,df,j)
        endif
!  MR: cases i=j/=k etc. missing
      elseif ((i==1.and.j==2.and.k==3).or.(i==1.and.j==3.and.k==2)) then
        xyz: if (nxgrid/=1.and.nygrid/=1.and.nzgrid/=1) then
          fac= (1./6.0*dx_1(l1:l2)**4) * (1./60.0*dy_1(m)) * (1/60.0*dz_1(n))
          df = fac*(&
               ( 113400.0 * f( l1:l2 , m - 1 , n - 1 , ik ))+&
               ( -78975.0 * f( l1-1:l2-1 , m - 1 , n - 1 , ik ))+&
               ( 24300.0 * f( l1-2:l2-2 , m - 1 , n - 1 , ik ))+&
               ( -2025.0 * f( l1-3:l2-3 , m - 1 , n - 1 , ik ))+&
               ( -78975.0 * f( l1+1:l2+1 , m - 1 , n - 1 , ik ))+&
               ( 24300.0 * f( l1+2:l2+2 , m - 1 , n - 1 , ik ))+&
               ( -2025.0 * f( l1+3:l2+3 , m - 1 , n - 1 , ik ))+&
               ( -22680.0 * f( l1:l2 , m - 2 , n - 1 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 2 , n - 1 , ik ))+&
               ( -4860.0 * f( l1-2:l2-2 , m - 2 , n - 1 , ik ))+&
               ( 405.0 * f( l1-3:l2-3 , m - 2 , n - 1 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m - 2 , n - 1 , ik ))+&
               ( -4860.0 * f( l1+2:l2+2 , m - 2 , n - 1 , ik ))+&
               ( 405.0 * f( l1+3:l2+3 , m - 2 , n - 1 , ik ))+&
               ( 2520.0 * f( l1:l2 , m - 3 , n - 1 , ik ))+&
               ( -1755.0 * f( l1-1:l2-1 , m - 3 , n - 1 , ik ))+&
               ( 540.0 * f( l1-2:l2-2 , m - 3 , n - 1 , ik ))+&
               ( -45.0 * f( l1-3:l2-3 , m - 3 , n - 1 , ik ))+&
               ( -1755.0 * f( l1+1:l2+1 , m - 3 , n - 1 , ik ))+&
               ( 540.0 * f( l1+2:l2+2 , m - 3 , n - 1 , ik ))+&
               ( -45.0 * f( l1+3:l2+3 , m - 3 , n - 1 , ik ))+&
               ( -113400.0 * f( l1:l2 , m + 1 , n - 1 , ik ))+&
               ( 78975.0 * f( l1-1:l2-1 , m + 1 , n - 1 , ik ))+&
               ( -24300.0 * f( l1-2:l2-2 , m + 1 , n - 1 , ik ))+&
               ( 2025.0 * f( l1-3:l2-3 , m + 1 , n - 1 , ik ))+&
               ( 78975.0 * f( l1+1:l2+1 , m + 1 , n - 1 , ik ))+&
               ( -24300.0 * f( l1+2:l2+2 , m + 1 , n - 1 , ik ))+&
               ( 2025.0 * f( l1+3:l2+3 , m + 1 , n - 1 , ik ))+&
               ( 22680.0 * f( l1:l2 , m + 2 , n - 1 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 2 , n - 1 , ik ))+&
               ( 4860.0 * f( l1-2:l2-2 , m + 2 , n - 1 , ik ))+&
               ( -405.0 * f( l1-3:l2-3 , m + 2 , n - 1 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m + 2 , n - 1 , ik ))+&
               ( 4860.0 * f( l1+2:l2+2 , m + 2 , n - 1 , ik ))+&
               ( -405.0 * f( l1+3:l2+3 , m + 2 , n - 1 , ik ))+&
               ( -2520.0 * f( l1:l2 , m + 3 , n - 1 , ik ))+&
               ( 1755.0 * f( l1-1:l2-1 , m + 3 , n - 1 , ik ))+&
               ( -540.0 * f( l1-2:l2-2 , m + 3 , n - 1 , ik ))+&
               ( 45.0 * f( l1-3:l2-3 , m + 3 , n - 1 , ik ))+&
               ( 1755.0 * f( l1+1:l2+1 , m + 3 , n - 1 , ik ))+&
               ( -540.0 * f( l1+2:l2+2 , m + 3 , n - 1 , ik ))+&
               ( 45.0 * f( l1+3:l2+3 , m + 3 , n - 1 , ik ))+&
               ( -22680.0 * f( l1:l2 , m - 1 , n - 2 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 1 , n - 2 , ik ))+&
               ( -4860.0 * f( l1-2:l2-2 , m - 1 , n - 2 , ik ))+&
               ( 405.0 * f( l1-3:l2-3 , m - 1 , n - 2 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m - 1 , n - 2 , ik ))+&
               ( -4860.0 * f( l1+2:l2+2 , m - 1 , n - 2 , ik ))+&
               ( 405.0 * f( l1+3:l2+3 , m - 1 , n - 2 , ik ))+&
               ( 4536.0 * f( l1:l2 , m - 2 , n - 2 , ik ))+&
               ( -3159.0 * f( l1-1:l2-1 , m - 2 , n - 2 , ik ))+&
               ( 972.0 * f( l1-2:l2-2 , m - 2 , n - 2 , ik ))+&
               ( -81.0 * f( l1-3:l2-3 , m - 2 , n - 2 , ik ))+&
               ( -3159.0 * f( l1+1:l2+1 , m - 2 , n - 2 , ik ))+&
               ( 972.0 * f( l1+2:l2+2 , m - 2 , n - 2 , ik ))+&
               ( -81.0 * f( l1+3:l2+3 , m - 2 , n - 2 , ik ))+&
               ( -504.0 * f( l1:l2 , m - 3 , n - 2 , ik ))+&
               ( 351.0 * f( l1-1:l2-1 , m - 3 , n - 2 , ik ))+&
               ( -108.0 * f( l1-2:l2-2 , m - 3 , n - 2 , ik ))+&
               ( 9.0 * f( l1-3:l2-3 , m - 3 , n - 2 , ik ))+&
               ( 351.0 * f( l1+1:l2+1 , m - 3 , n - 2 , ik ))+&
               ( -108.0 * f( l1+2:l2+2 , m - 3 , n - 2 , ik ))+&
               ( 9.0 * f( l1+3:l2+3 , m - 3 , n - 2 , ik ))+&
               ( 22680.0 * f( l1:l2 , m + 1 , n - 2 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 1 , n - 2 , ik ))+&
               ( 4860.0 * f( l1-2:l2-2 , m + 1 , n - 2 , ik ))+&
               ( -405.0 * f( l1-3:l2-3 , m + 1 , n - 2 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m + 1 , n - 2 , ik ))+&
               ( 4860.0 * f( l1+2:l2+2 , m + 1 , n - 2 , ik ))+&
               ( -405.0 * f( l1+3:l2+3 , m + 1 , n - 2 , ik ))+&
               ( -4536.0 * f( l1:l2 , m + 2 , n - 2 , ik ))+&
               ( 3159.0 * f( l1-1:l2-1 , m + 2 , n - 2 , ik ))+&
               ( -972.0 * f( l1-2:l2-2 , m + 2 , n - 2 , ik ))+&
               ( 81.0 * f( l1-3:l2-3 , m + 2 , n - 2 , ik ))+&
               ( 3159.0 * f( l1+1:l2+1 , m + 2 , n - 2 , ik ))+&
               ( -972.0 * f( l1+2:l2+2 , m + 2 , n - 2 , ik ))+&
               ( 81.0 * f( l1+3:l2+3 , m + 2 , n - 2 , ik ))+&
               ( 504.0 * f( l1:l2 , m + 3 , n - 2 , ik ))+&
               ( -351.0 * f( l1-1:l2-1 , m + 3 , n - 2 , ik ))+&
               ( 108.0 * f( l1-2:l2-2 , m + 3 , n - 2 , ik ))+&
               ( -9.0 * f( l1-3:l2-3 , m + 3 , n - 2 , ik ))+&
               ( -351.0 * f( l1+1:l2+1 , m + 3 , n - 2 , ik ))+&
               ( 108.0 * f( l1+2:l2+2 , m + 3 , n - 2 , ik ))+&
               ( -9.0 * f( l1+3:l2+3 , m + 3 , n - 2 , ik ))+&
               ( 2520.0 * f( l1:l2 , m - 1 , n - 3 , ik ))+&
               ( -1755.0 * f( l1-1:l2-1 , m - 1 , n - 3 , ik ))+&
               ( 540.0 * f( l1-2:l2-2 , m - 1 , n - 3 , ik ))+&
               ( -45.0 * f( l1-3:l2-3 , m - 1 , n - 3 , ik ))+&
               ( -1755.0 * f( l1+1:l2+1 , m - 1 , n - 3 , ik ))+&
               ( 540.0 * f( l1+2:l2+2 , m - 1 , n - 3 , ik ))+&
               ( -45.0 * f( l1+3:l2+3 , m - 1 , n - 3 , ik ))+&
               ( -504.0 * f( l1:l2 , m - 2 , n - 3 , ik ))+&
               ( 351.0 * f( l1-1:l2-1 , m - 2 , n - 3 , ik ))+&
               ( -108.0 * f( l1-2:l2-2 , m - 2 , n - 3 , ik ))+&
               ( 9.0 * f( l1-3:l2-3 , m - 2 , n - 3 , ik ))+&
               ( 351.0 * f( l1+1:l2+1 , m - 2 , n - 3 , ik ))+&
               ( -108.0 * f( l1+2:l2+2 , m - 2 , n - 3 , ik ))+&
               ( 9.0 * f( l1+3:l2+3 , m - 2 , n - 3 , ik ))+&
               ( 56.0 * f( l1:l2 , m - 3 , n - 3 , ik ))+&
               ( -39.0 * f( l1-1:l2-1 , m - 3 , n - 3 , ik ))+&
               ( 12.0 * f( l1-2:l2-2 , m - 3 , n - 3 , ik ))+&
               ( -1.0 * f( l1-3:l2-3 , m - 3 , n - 3 , ik ))+&
               ( -39.0 * f( l1+1:l2+1 , m - 3 , n - 3 , ik ))+&
               ( 12.0 * f( l1+2:l2+2 , m - 3 , n - 3 , ik ))+&
               ( -1.0 * f( l1+3:l2+3 , m - 3 , n - 3 , ik ))+&
               ( -2520.0 * f( l1:l2 , m + 1 , n - 3 , ik ))+&
               ( 1755.0 * f( l1-1:l2-1 , m + 1 , n - 3 , ik ))+&
               ( -540.0 * f( l1-2:l2-2 , m + 1 , n - 3 , ik ))+&
               ( 45.0 * f( l1-3:l2-3 , m + 1 , n - 3 , ik ))+&
               ( 1755.0 * f( l1+1:l2+1 , m + 1 , n - 3 , ik ))+&
               ( -540.0 * f( l1+2:l2+2 , m + 1 , n - 3 , ik ))+&
               ( 45.0 * f( l1+3:l2+3 , m + 1 , n - 3 , ik ))+&
               ( 504.0 * f( l1:l2 , m + 2 , n - 3 , ik ))+&
               ( -351.0 * f( l1-1:l2-1 , m + 2 , n - 3 , ik ))+&
               ( 108.0 * f( l1-2:l2-2 , m + 2 , n - 3 , ik ))+&
               ( -9.0 * f( l1-3:l2-3 , m + 2 , n - 3 , ik ))+&
               ( -351.0 * f( l1+1:l2+1 , m + 2 , n - 3 , ik ))+&
               ( 108.0 * f( l1+2:l2+2 , m + 2 , n - 3 , ik ))+&
               ( -9.0 * f( l1+3:l2+3 , m + 2 , n - 3 , ik ))+&
               ( -56.0 * f( l1:l2 , m + 3 , n - 3 , ik ))+&
               ( 39.0 * f( l1-1:l2-1 , m + 3 , n - 3 , ik ))+&
               ( -12.0 * f( l1-2:l2-2 , m + 3 , n - 3 , ik ))+&
               ( 1.0 * f( l1-3:l2-3 , m + 3 , n - 3 , ik ))+&
               ( 39.0 * f( l1+1:l2+1 , m + 3 , n - 3 , ik ))+&
               ( -12.0 * f( l1+2:l2+2 , m + 3 , n - 3 , ik ))+&
               ( 1.0 * f( l1+3:l2+3 , m + 3 , n - 3 , ik ))+&
               ( -113400.0 * f( l1:l2 , m - 1 , n + 1 , ik ))+&
               ( 78975.0 * f( l1-1:l2-1 , m - 1 , n + 1 , ik ))+&
               ( -24300.0 * f( l1-2:l2-2 , m - 1 , n + 1 , ik ))+&
               ( 2025.0 * f( l1-3:l2-3 , m - 1 , n + 1 , ik ))+&
               ( 78975.0 * f( l1+1:l2+1 , m - 1 , n + 1 , ik ))+&
               ( -24300.0 * f( l1+2:l2+2 , m - 1 , n + 1 , ik ))+&
               ( 2025.0 * f( l1+3:l2+3 , m - 1 , n + 1 , ik ))+&
               ( 22680.0 * f( l1:l2 , m - 2 , n + 1 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m - 2 , n + 1 , ik ))+&
               ( 4860.0 * f( l1-2:l2-2 , m - 2 , n + 1 , ik ))+&
               ( -405.0 * f( l1-3:l2-3 , m - 2 , n + 1 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 2 , n + 1 , ik ))+&
               ( 4860.0 * f( l1+2:l2+2 , m - 2 , n + 1 , ik ))+&
               ( -405.0 * f( l1+3:l2+3 , m - 2 , n + 1 , ik ))+&
               ( -2520.0 * f( l1:l2 , m - 3 , n + 1 , ik ))+&
               ( 1755.0 * f( l1-1:l2-1 , m - 3 , n + 1 , ik ))+&
               ( -540.0 * f( l1-2:l2-2 , m - 3 , n + 1 , ik ))+&
               ( 45.0 * f( l1-3:l2-3 , m - 3 , n + 1 , ik ))+&
               ( 1755.0 * f( l1+1:l2+1 , m - 3 , n + 1 , ik ))+&
               ( -540.0 * f( l1+2:l2+2 , m - 3 , n + 1 , ik ))+&
               ( 45.0 * f( l1+3:l2+3 , m - 3 , n + 1 , ik ))+&
               ( 113400.0 * f( l1:l2 , m + 1 , n + 1 , ik ))+&
               ( -78975.0 * f( l1-1:l2-1 , m + 1 , n + 1 , ik ))+&
               ( 24300.0 * f( l1-2:l2-2 , m + 1 , n + 1 , ik ))+&
               ( -2025.0 * f( l1-3:l2-3 , m + 1 , n + 1 , ik ))+&
               ( -78975.0 * f( l1+1:l2+1 , m + 1 , n + 1 , ik ))+&
               ( 24300.0 * f( l1+2:l2+2 , m + 1 , n + 1 , ik ))+&
               ( -2025.0 * f( l1+3:l2+3 , m + 1 , n + 1 , ik ))+&
               ( -22680.0 * f( l1:l2 , m + 2 , n + 1 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m + 2 , n + 1 , ik ))+&
               ( -4860.0 * f( l1-2:l2-2 , m + 2 , n + 1 , ik ))+&
               ( 405.0 * f( l1-3:l2-3 , m + 2 , n + 1 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 2 , n + 1 , ik ))+&
               ( -4860.0 * f( l1+2:l2+2 , m + 2 , n + 1 , ik ))+&
               ( 405.0 * f( l1+3:l2+3 , m + 2 , n + 1 , ik ))+&
               ( 2520.0 * f( l1:l2 , m + 3 , n + 1 , ik ))+&
               ( -1755.0 * f( l1-1:l2-1 , m + 3 , n + 1 , ik ))+&
               ( 540.0 * f( l1-2:l2-2 , m + 3 , n + 1 , ik ))+&
               ( -45.0 * f( l1-3:l2-3 , m + 3 , n + 1 , ik ))+&
               ( -1755.0 * f( l1+1:l2+1 , m + 3 , n + 1 , ik ))+&
               ( 540.0 * f( l1+2:l2+2 , m + 3 , n + 1 , ik ))+&
               ( -45.0 * f( l1+3:l2+3 , m + 3 , n + 1 , ik ))+&
               ( 22680.0 * f( l1:l2 , m - 1 , n + 2 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m - 1 , n + 2 , ik ))+&
               ( 4860.0 * f( l1-2:l2-2 , m - 1 , n + 2 , ik ))+&
               ( -405.0 * f( l1-3:l2-3 , m - 1 , n + 2 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 1 , n + 2 , ik ))+&
               ( 4860.0 * f( l1+2:l2+2 , m - 1 , n + 2 , ik ))+&
               ( -405.0 * f( l1+3:l2+3 , m - 1 , n + 2 , ik ))+&
               ( -4536.0 * f( l1:l2 , m - 2 , n + 2 , ik ))+&
               ( 3159.0 * f( l1-1:l2-1 , m - 2 , n + 2 , ik ))+&
               ( -972.0 * f( l1-2:l2-2 , m - 2 , n + 2 , ik ))+&
               ( 81.0 * f( l1-3:l2-3 , m - 2 , n + 2 , ik ))+&
               ( 3159.0 * f( l1+1:l2+1 , m - 2 , n + 2 , ik ))+&
               ( -972.0 * f( l1+2:l2+2 , m - 2 , n + 2 , ik ))+&
               ( 81.0 * f( l1+3:l2+3 , m - 2 , n + 2 , ik ))+&
               ( 504.0 * f( l1:l2 , m - 3 , n + 2 , ik ))+&
               ( -351.0 * f( l1-1:l2-1 , m - 3 , n + 2 , ik ))+&
               ( 108.0 * f( l1-2:l2-2 , m - 3 , n + 2 , ik ))+&
               ( -9.0 * f( l1-3:l2-3 , m - 3 , n + 2 , ik ))+&
               ( -351.0 * f( l1+1:l2+1 , m - 3 , n + 2 , ik ))+&
               ( 108.0 * f( l1+2:l2+2 , m - 3 , n + 2 , ik ))+&
               ( -9.0 * f( l1+3:l2+3 , m - 3 , n + 2 , ik ))+&
               ( -22680.0 * f( l1:l2 , m + 1 , n + 2 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m + 1 , n + 2 , ik ))+&
               ( -4860.0 * f( l1-2:l2-2 , m + 1 , n + 2 , ik ))+&
               ( 405.0 * f( l1-3:l2-3 , m + 1 , n + 2 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 1 , n + 2 , ik ))+&
               ( -4860.0 * f( l1+2:l2+2 , m + 1 , n + 2 , ik ))+&
               ( 405.0 * f( l1+3:l2+3 , m + 1 , n + 2 , ik ))+&
               ( 4536.0 * f( l1:l2 , m + 2 , n + 2 , ik ))+&
               ( -3159.0 * f( l1-1:l2-1 , m + 2 , n + 2 , ik ))+&
               ( 972.0 * f( l1-2:l2-2 , m + 2 , n + 2 , ik ))+&
               ( -81.0 * f( l1-3:l2-3 , m + 2 , n + 2 , ik ))+&
               ( -3159.0 * f( l1+1:l2+1 , m + 2 , n + 2 , ik ))+&
               ( 972.0 * f( l1+2:l2+2 , m + 2 , n + 2 , ik ))+&
               ( -81.0 * f( l1+3:l2+3 , m + 2 , n + 2 , ik ))+&
               ( -504.0 * f( l1:l2 , m + 3 , n + 2 , ik ))+&
               ( 351.0 * f( l1-1:l2-1 , m + 3 , n + 2 , ik ))+&
               ( -108.0 * f( l1-2:l2-2 , m + 3 , n + 2 , ik ))+&
               ( 9.0 * f( l1-3:l2-3 , m + 3 , n + 2 , ik ))+&
               ( 351.0 * f( l1+1:l2+1 , m + 3 , n + 2 , ik ))+&
               ( -108.0 * f( l1+2:l2+2 , m + 3 , n + 2 , ik ))+&
               ( 9.0 * f( l1+3:l2+3 , m + 3 , n + 2 , ik ))+&
               ( -2520.0 * f( l1:l2 , m - 1 , n + 3 , ik ))+&
               ( 1755.0 * f( l1-1:l2-1 , m - 1 , n + 3 , ik ))+&
               ( -540.0 * f( l1-2:l2-2 , m - 1 , n + 3 , ik ))+&
               ( 45.0 * f( l1-3:l2-3 , m - 1 , n + 3 , ik ))+&
               ( 1755.0 * f( l1+1:l2+1 , m - 1 , n + 3 , ik ))+&
               ( -540.0 * f( l1+2:l2+2 , m - 1 , n + 3 , ik ))+&
               ( 45.0 * f( l1+3:l2+3 , m - 1 , n + 3 , ik ))+&
               ( 504.0 * f( l1:l2 , m - 2 , n + 3 , ik ))+&
               ( -351.0 * f( l1-1:l2-1 , m - 2 , n + 3 , ik ))+&
               ( 108.0 * f( l1-2:l2-2 , m - 2 , n + 3 , ik ))+&
               ( -9.0 * f( l1-3:l2-3 , m - 2 , n + 3 , ik ))+&
               ( -351.0 * f( l1+1:l2+1 , m - 2 , n + 3 , ik ))+&
               ( 108.0 * f( l1+2:l2+2 , m - 2 , n + 3 , ik ))+&
               ( -9.0 * f( l1+3:l2+3 , m - 2 , n + 3 , ik ))+&
               ( -56.0 * f( l1:l2 , m - 3 , n + 3 , ik ))+&
               ( 39.0 * f( l1-1:l2-1 , m - 3 , n + 3 , ik ))+&
               ( -12.0 * f( l1-2:l2-2 , m - 3 , n + 3 , ik ))+&
               ( 1.0 * f( l1-3:l2-3 , m - 3 , n + 3 , ik ))+&
               ( 39.0 * f( l1+1:l2+1 , m - 3 , n + 3 , ik ))+&
               ( -12.0 * f( l1+2:l2+2 , m - 3 , n + 3 , ik ))+&
               ( 1.0 * f( l1+3:l2+3 , m - 3 , n + 3 , ik ))+&
               ( 2520.0 * f( l1:l2 , m + 1 , n + 3 , ik ))+&
               ( -1755.0 * f( l1-1:l2-1 , m + 1 , n + 3 , ik ))+&
               ( 540.0 * f( l1-2:l2-2 , m + 1 , n + 3 , ik ))+&
               ( -45.0 * f( l1-3:l2-3 , m + 1 , n + 3 , ik ))+&
               ( -1755.0 * f( l1+1:l2+1 , m + 1 , n + 3 , ik ))+&
               ( 540.0 * f( l1+2:l2+2 , m + 1 , n + 3 , ik ))+&
               ( -45.0 * f( l1+3:l2+3 , m + 1 , n + 3 , ik ))+&
               ( -504.0 * f( l1:l2 , m + 2 , n + 3 , ik ))+&
               ( 351.0 * f( l1-1:l2-1 , m + 2 , n + 3 , ik ))+&
               ( -108.0 * f( l1-2:l2-2 , m + 2 , n + 3 , ik ))+&
               ( 9.0 * f( l1-3:l2-3 , m + 2 , n + 3 , ik ))+&
               ( 351.0 * f( l1+1:l2+1 , m + 2 , n + 3 , ik ))+&
               ( -108.0 * f( l1+2:l2+2 , m + 2 , n + 3 , ik ))+&
               ( 9.0 * f( l1+3:l2+3 , m + 2 , n + 3 , ik ))+&
               ( 56.0 * f( l1:l2 , m + 3 , n + 3 , ik ))+&
               ( -39.0 * f( l1-1:l2-1 , m + 3 , n + 3 , ik ))+&
               ( 12.0 * f( l1-2:l2-2 , m + 3 , n + 3 , ik ))+&
               ( -1.0 * f( l1-3:l2-3 , m + 3 , n + 3 , ik ))+&
               ( -39.0 * f( l1+1:l2+1 , m + 3 , n + 3 , ik ))+&
               ( 12.0 * f( l1+2:l2+2 , m + 3 , n + 3 , ik ))+&
               ( -1.0 * f( l1+3:l2+3 , m + 3 , n + 3 , ik ))&
               )
        else
          df=0.0
        endif xyz
!        
      elseif ((i==2.and.j==1.and.k==3).or.(i==2.and.j==3.and.k==1)) then
        yxz: if (nxgrid/=1.and.nygrid/=1.and.nzgrid/=1) then
          fac= (1./6.0*dy_1(m)**4) * (1./60.0*dx_1(l1:l2)) * (1/60.0*dz_1(n))
          df = fac*(&
               ( 113400.0 * f( l1-1:l2-1 , m , n - 1 , ik ))+&
               ( -78975.0 * f( l1-1:l2-1 , m - 1 , n - 1 , ik ))+&
               ( 24300.0 * f( l1-1:l2-1 , m - 2 , n - 1 , ik ))+&
               ( -2025.0 * f( l1-1:l2-1 , m - 3 , n - 1 , ik ))+&
               ( -78975.0 * f( l1-1:l2-1 , m + 1 , n - 1 , ik ))+&
               ( 24300.0 * f( l1-1:l2-1 , m + 2 , n - 1 , ik ))+&
               ( -2025.0 * f( l1-1:l2-1 , m + 3 , n - 1 , ik ))+&
               ( -22680.0 * f( l1-2:l2-2 , m , n - 1 , ik ))+&
               ( 15795.0 * f( l1-2:l2-2 , m - 1 , n - 1 , ik ))+&
               ( -4860.0 * f( l1-2:l2-2 , m - 2 , n - 1 , ik ))+&
               ( 405.0 * f( l1-2:l2-2 , m - 3 , n - 1 , ik ))+&
               ( 15795.0 * f( l1-2:l2-2 , m + 1 , n - 1 , ik ))+&
               ( -4860.0 * f( l1-2:l2-2 , m + 2 , n - 1 , ik ))+&
               ( 405.0 * f( l1-2:l2-2 , m + 3 , n - 1 , ik ))+&
               ( 2520.0 * f( l1-3:l2-3 , m , n - 1 , ik ))+&
               ( -1755.0 * f( l1-3:l2-3 , m - 1 , n - 1 , ik ))+&
               ( 540.0 * f( l1-3:l2-3 , m - 2 , n - 1 , ik ))+&
               ( -45.0 * f( l1-3:l2-3 , m - 3 , n - 1 , ik ))+&
               ( -1755.0 * f( l1-3:l2-3 , m + 1 , n - 1 , ik ))+&
               ( 540.0 * f( l1-3:l2-3 , m + 2 , n - 1 , ik ))+&
               ( -45.0 * f( l1-3:l2-3 , m + 3 , n - 1 , ik ))+&
               ( -113400.0 * f( l1+1:l2+1 , m , n - 1 , ik ))+&
               ( 78975.0 * f( l1+1:l2+1 , m - 1 , n - 1 , ik ))+&
               ( -24300.0 * f( l1+1:l2+1 , m - 2 , n - 1 , ik ))+&
               ( 2025.0 * f( l1+1:l2+1 , m - 3 , n - 1 , ik ))+&
               ( 78975.0 * f( l1+1:l2+1 , m + 1 , n - 1 , ik ))+&
               ( -24300.0 * f( l1+1:l2+1 , m + 2 , n - 1 , ik ))+&
               ( 2025.0 * f( l1+1:l2+1 , m + 3 , n - 1 , ik ))+&
               ( 22680.0 * f( l1+2:l2+2 , m , n - 1 , ik ))+&
               ( -15795.0 * f( l1+2:l2+2 , m - 1 , n - 1 , ik ))+&
               ( 4860.0 * f( l1+2:l2+2 , m - 2 , n - 1 , ik ))+&
               ( -405.0 * f( l1+2:l2+2 , m - 3 , n - 1 , ik ))+&
               ( -15795.0 * f( l1+2:l2+2 , m + 1 , n - 1 , ik ))+&
               ( 4860.0 * f( l1+2:l2+2 , m + 2 , n - 1 , ik ))+&
               ( -405.0 * f( l1+2:l2+2 , m + 3 , n - 1 , ik ))+&
               ( -2520.0 * f( l1+3:l2+3 , m , n - 1 , ik ))+&
               ( 1755.0 * f( l1+3:l2+3 , m - 1 , n - 1 , ik ))+&
               ( -540.0 * f( l1+3:l2+3 , m - 2 , n - 1 , ik ))+&
               ( 45.0 * f( l1+3:l2+3 , m - 3 , n - 1 , ik ))+&
               ( 1755.0 * f( l1+3:l2+3 , m + 1 , n - 1 , ik ))+&
               ( -540.0 * f( l1+3:l2+3 , m + 2 , n - 1 , ik ))+&
               ( 45.0 * f( l1+3:l2+3 , m + 3 , n - 1 , ik ))+&
               ( -22680.0 * f( l1-1:l2-1 , m , n - 2 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 1 , n - 2 , ik ))+&
               ( -4860.0 * f( l1-1:l2-1 , m - 2 , n - 2 , ik ))+&
               ( 405.0 * f( l1-1:l2-1 , m - 3 , n - 2 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m + 1 , n - 2 , ik ))+&
               ( -4860.0 * f( l1-1:l2-1 , m + 2 , n - 2 , ik ))+&
               ( 405.0 * f( l1-1:l2-1 , m + 3 , n - 2 , ik ))+&
               ( 4536.0 * f( l1-2:l2-2 , m , n - 2 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m - 1 , n - 2 , ik ))+&
               ( 972.0 * f( l1-2:l2-2 , m - 2 , n - 2 , ik ))+&
               ( -81.0 * f( l1-2:l2-2 , m - 3 , n - 2 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m + 1 , n - 2 , ik ))+&
               ( 972.0 * f( l1-2:l2-2 , m + 2 , n - 2 , ik ))+&
               ( -81.0 * f( l1-2:l2-2 , m + 3 , n - 2 , ik ))+&
               ( -504.0 * f( l1-3:l2-3 , m , n - 2 , ik ))+&
               ( 351.0 * f( l1-3:l2-3 , m - 1 , n - 2 , ik ))+&
               ( -108.0 * f( l1-3:l2-3 , m - 2 , n - 2 , ik ))+&
               ( 9.0 * f( l1-3:l2-3 , m - 3 , n - 2 , ik ))+&
               ( 351.0 * f( l1-3:l2-3 , m + 1 , n - 2 , ik ))+&
               ( -108.0 * f( l1-3:l2-3 , m + 2 , n - 2 , ik ))+&
               ( 9.0 * f( l1-3:l2-3 , m + 3 , n - 2 , ik ))+&
               ( 22680.0 * f( l1+1:l2+1 , m , n - 2 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 1 , n - 2 , ik ))+&
               ( 4860.0 * f( l1+1:l2+1 , m - 2 , n - 2 , ik ))+&
               ( -405.0 * f( l1+1:l2+1 , m - 3 , n - 2 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m + 1 , n - 2 , ik ))+&
               ( 4860.0 * f( l1+1:l2+1 , m + 2 , n - 2 , ik ))+&
               ( -405.0 * f( l1+1:l2+1 , m + 3 , n - 2 , ik ))+&
               ( -4536.0 * f( l1+2:l2+2 , m , n - 2 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m - 1 , n - 2 , ik ))+&
               ( -972.0 * f( l1+2:l2+2 , m - 2 , n - 2 , ik ))+&
               ( 81.0 * f( l1+2:l2+2 , m - 3 , n - 2 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m + 1 , n - 2 , ik ))+&
               ( -972.0 * f( l1+2:l2+2 , m + 2 , n - 2 , ik ))+&
               ( 81.0 * f( l1+2:l2+2 , m + 3 , n - 2 , ik ))+&
               ( 504.0 * f( l1+3:l2+3 , m , n - 2 , ik ))+&
               ( -351.0 * f( l1+3:l2+3 , m - 1 , n - 2 , ik ))+&
               ( 108.0 * f( l1+3:l2+3 , m - 2 , n - 2 , ik ))+&
               ( -9.0 * f( l1+3:l2+3 , m - 3 , n - 2 , ik ))+&
               ( -351.0 * f( l1+3:l2+3 , m + 1 , n - 2 , ik ))+&
               ( 108.0 * f( l1+3:l2+3 , m + 2 , n - 2 , ik ))+&
               ( -9.0 * f( l1+3:l2+3 , m + 3 , n - 2 , ik ))+&
               ( 2520.0 * f( l1-1:l2-1 , m , n - 3 , ik ))+&
               ( -1755.0 * f( l1-1:l2-1 , m - 1 , n - 3 , ik ))+&
               ( 540.0 * f( l1-1:l2-1 , m - 2 , n - 3 , ik ))+&
               ( -45.0 * f( l1-1:l2-1 , m - 3 , n - 3 , ik ))+&
               ( -1755.0 * f( l1-1:l2-1 , m + 1 , n - 3 , ik ))+&
               ( 540.0 * f( l1-1:l2-1 , m + 2 , n - 3 , ik ))+&
               ( -45.0 * f( l1-1:l2-1 , m + 3 , n - 3 , ik ))+&
               ( -504.0 * f( l1-2:l2-2 , m , n - 3 , ik ))+&
               ( 351.0 * f( l1-2:l2-2 , m - 1 , n - 3 , ik ))+&
               ( -108.0 * f( l1-2:l2-2 , m - 2 , n - 3 , ik ))+&
               ( 9.0 * f( l1-2:l2-2 , m - 3 , n - 3 , ik ))+&
               ( 351.0 * f( l1-2:l2-2 , m + 1 , n - 3 , ik ))+&
               ( -108.0 * f( l1-2:l2-2 , m + 2 , n - 3 , ik ))+&
               ( 9.0 * f( l1-2:l2-2 , m + 3 , n - 3 , ik ))+&
               ( 56.0 * f( l1-3:l2-3 , m , n - 3 , ik ))+&
               ( -39.0 * f( l1-3:l2-3 , m - 1 , n - 3 , ik ))+&
               ( 12.0 * f( l1-3:l2-3 , m - 2 , n - 3 , ik ))+&
               ( -1.0 * f( l1-3:l2-3 , m - 3 , n - 3 , ik ))+&
               ( -39.0 * f( l1-3:l2-3 , m + 1 , n - 3 , ik ))+&
               ( 12.0 * f( l1-3:l2-3 , m + 2 , n - 3 , ik ))+&
               ( -1.0 * f( l1-3:l2-3 , m + 3 , n - 3 , ik ))+&
               ( -2520.0 * f( l1+1:l2+1 , m , n - 3 , ik ))+&
               ( 1755.0 * f( l1+1:l2+1 , m - 1 , n - 3 , ik ))+&
               ( -540.0 * f( l1+1:l2+1 , m - 2 , n - 3 , ik ))+&
               ( 45.0 * f( l1+1:l2+1 , m - 3 , n - 3 , ik ))+&
               ( 1755.0 * f( l1+1:l2+1 , m + 1 , n - 3 , ik ))+&
               ( -540.0 * f( l1+1:l2+1 , m + 2 , n - 3 , ik ))+&
               ( 45.0 * f( l1+1:l2+1 , m + 3 , n - 3 , ik ))+&
               ( 504.0 * f( l1+2:l2+2 , m , n - 3 , ik ))+&
               ( -351.0 * f( l1+2:l2+2 , m - 1 , n - 3 , ik ))+&
               ( 108.0 * f( l1+2:l2+2 , m - 2 , n - 3 , ik ))+&
               ( -9.0 * f( l1+2:l2+2 , m - 3 , n - 3 , ik ))+&
               ( -351.0 * f( l1+2:l2+2 , m + 1 , n - 3 , ik ))+&
               ( 108.0 * f( l1+2:l2+2 , m + 2 , n - 3 , ik ))+&
               ( -9.0 * f( l1+2:l2+2 , m + 3 , n - 3 , ik ))+&
               ( -56.0 * f( l1+3:l2+3 , m , n - 3 , ik ))+&
               ( 39.0 * f( l1+3:l2+3 , m - 1 , n - 3 , ik ))+&
               ( -12.0 * f( l1+3:l2+3 , m - 2 , n - 3 , ik ))+&
               ( 1.0 * f( l1+3:l2+3 , m - 3 , n - 3 , ik ))+&
               ( 39.0 * f( l1+3:l2+3 , m + 1 , n - 3 , ik ))+&
               ( -12.0 * f( l1+3:l2+3 , m + 2 , n - 3 , ik ))+&
               ( 1.0 * f( l1+3:l2+3 , m + 3 , n - 3 , ik ))+&
               ( -113400.0 * f( l1-1:l2-1 , m , n + 1 , ik ))+&
               ( 78975.0 * f( l1-1:l2-1 , m - 1 , n + 1 , ik ))+&
               ( -24300.0 * f( l1-1:l2-1 , m - 2 , n + 1 , ik ))+&
               ( 2025.0 * f( l1-1:l2-1 , m - 3 , n + 1 , ik ))+&
               ( 78975.0 * f( l1-1:l2-1 , m + 1 , n + 1 , ik ))+&
               ( -24300.0 * f( l1-1:l2-1 , m + 2 , n + 1 , ik ))+&
               ( 2025.0 * f( l1-1:l2-1 , m + 3 , n + 1 , ik ))+&
               ( 22680.0 * f( l1-2:l2-2 , m , n + 1 , ik ))+&
               ( -15795.0 * f( l1-2:l2-2 , m - 1 , n + 1 , ik ))+&
               ( 4860.0 * f( l1-2:l2-2 , m - 2 , n + 1 , ik ))+&
               ( -405.0 * f( l1-2:l2-2 , m - 3 , n + 1 , ik ))+&
               ( -15795.0 * f( l1-2:l2-2 , m + 1 , n + 1 , ik ))+&
               ( 4860.0 * f( l1-2:l2-2 , m + 2 , n + 1 , ik ))+&
               ( -405.0 * f( l1-2:l2-2 , m + 3 , n + 1 , ik ))+&
               ( -2520.0 * f( l1-3:l2-3 , m , n + 1 , ik ))+&
               ( 1755.0 * f( l1-3:l2-3 , m - 1 , n + 1 , ik ))+&
               ( -540.0 * f( l1-3:l2-3 , m - 2 , n + 1 , ik ))+&
               ( 45.0 * f( l1-3:l2-3 , m - 3 , n + 1 , ik ))+&
               ( 1755.0 * f( l1-3:l2-3 , m + 1 , n + 1 , ik ))+&
               ( -540.0 * f( l1-3:l2-3 , m + 2 , n + 1 , ik ))+&
               ( 45.0 * f( l1-3:l2-3 , m + 3 , n + 1 , ik ))+&
               ( 113400.0 * f( l1+1:l2+1 , m , n + 1 , ik ))+&
               ( -78975.0 * f( l1+1:l2+1 , m - 1 , n + 1 , ik ))+&
               ( 24300.0 * f( l1+1:l2+1 , m - 2 , n + 1 , ik ))+&
               ( -2025.0 * f( l1+1:l2+1 , m - 3 , n + 1 , ik ))+&
               ( -78975.0 * f( l1+1:l2+1 , m + 1 , n + 1 , ik ))+&
               ( 24300.0 * f( l1+1:l2+1 , m + 2 , n + 1 , ik ))+&
               ( -2025.0 * f( l1+1:l2+1 , m + 3 , n + 1 , ik ))+&
               ( -22680.0 * f( l1+2:l2+2 , m , n + 1 , ik ))+&
               ( 15795.0 * f( l1+2:l2+2 , m - 1 , n + 1 , ik ))+&
               ( -4860.0 * f( l1+2:l2+2 , m - 2 , n + 1 , ik ))+&
               ( 405.0 * f( l1+2:l2+2 , m - 3 , n + 1 , ik ))+&
               ( 15795.0 * f( l1+2:l2+2 , m + 1 , n + 1 , ik ))+&
               ( -4860.0 * f( l1+2:l2+2 , m + 2 , n + 1 , ik ))+&
               ( 405.0 * f( l1+2:l2+2 , m + 3 , n + 1 , ik ))+&
               ( 2520.0 * f( l1+3:l2+3 , m , n + 1 , ik ))+&
               ( -1755.0 * f( l1+3:l2+3 , m - 1 , n + 1 , ik ))+&
               ( 540.0 * f( l1+3:l2+3 , m - 2 , n + 1 , ik ))+&
               ( -45.0 * f( l1+3:l2+3 , m - 3 , n + 1 , ik ))+&
               ( -1755.0 * f( l1+3:l2+3 , m + 1 , n + 1 , ik ))+&
               ( 540.0 * f( l1+3:l2+3 , m + 2 , n + 1 , ik ))+&
               ( -45.0 * f( l1+3:l2+3 , m + 3 , n + 1 , ik ))+&
               ( 22680.0 * f( l1-1:l2-1 , m , n + 2 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m - 1 , n + 2 , ik ))+&
               ( 4860.0 * f( l1-1:l2-1 , m - 2 , n + 2 , ik ))+&
               ( -405.0 * f( l1-1:l2-1 , m - 3 , n + 2 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 1 , n + 2 , ik ))+&
               ( 4860.0 * f( l1-1:l2-1 , m + 2 , n + 2 , ik ))+&
               ( -405.0 * f( l1-1:l2-1 , m + 3 , n + 2 , ik ))+&
               ( -4536.0 * f( l1-2:l2-2 , m , n + 2 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m - 1 , n + 2 , ik ))+&
               ( -972.0 * f( l1-2:l2-2 , m - 2 , n + 2 , ik ))+&
               ( 81.0 * f( l1-2:l2-2 , m - 3 , n + 2 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m + 1 , n + 2 , ik ))+&
               ( -972.0 * f( l1-2:l2-2 , m + 2 , n + 2 , ik ))+&
               ( 81.0 * f( l1-2:l2-2 , m + 3 , n + 2 , ik ))+&
               ( 504.0 * f( l1-3:l2-3 , m , n + 2 , ik ))+&
               ( -351.0 * f( l1-3:l2-3 , m - 1 , n + 2 , ik ))+&
               ( 108.0 * f( l1-3:l2-3 , m - 2 , n + 2 , ik ))+&
               ( -9.0 * f( l1-3:l2-3 , m - 3 , n + 2 , ik ))+&
               ( -351.0 * f( l1-3:l2-3 , m + 1 , n + 2 , ik ))+&
               ( 108.0 * f( l1-3:l2-3 , m + 2 , n + 2 , ik ))+&
               ( -9.0 * f( l1-3:l2-3 , m + 3 , n + 2 , ik ))+&
               ( -22680.0 * f( l1+1:l2+1 , m , n + 2 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m - 1 , n + 2 , ik ))+&
               ( -4860.0 * f( l1+1:l2+1 , m - 2 , n + 2 , ik ))+&
               ( 405.0 * f( l1+1:l2+1 , m - 3 , n + 2 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 1 , n + 2 , ik ))+&
               ( -4860.0 * f( l1+1:l2+1 , m + 2 , n + 2 , ik ))+&
               ( 405.0 * f( l1+1:l2+1 , m + 3 , n + 2 , ik ))+&
               ( 4536.0 * f( l1+2:l2+2 , m , n + 2 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m - 1 , n + 2 , ik ))+&
               ( 972.0 * f( l1+2:l2+2 , m - 2 , n + 2 , ik ))+&
               ( -81.0 * f( l1+2:l2+2 , m - 3 , n + 2 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m + 1 , n + 2 , ik ))+&
               ( 972.0 * f( l1+2:l2+2 , m + 2 , n + 2 , ik ))+&
               ( -81.0 * f( l1+2:l2+2 , m + 3 , n + 2 , ik ))+&
               ( -504.0 * f( l1+3:l2+3 , m , n + 2 , ik ))+&
               ( 351.0 * f( l1+3:l2+3 , m - 1 , n + 2 , ik ))+&
               ( -108.0 * f( l1+3:l2+3 , m - 2 , n + 2 , ik ))+&
               ( 9.0 * f( l1+3:l2+3 , m - 3 , n + 2 , ik ))+&
               ( 351.0 * f( l1+3:l2+3 , m + 1 , n + 2 , ik ))+&
               ( -108.0 * f( l1+3:l2+3 , m + 2 , n + 2 , ik ))+&
               ( 9.0 * f( l1+3:l2+3 , m + 3 , n + 2 , ik ))+&
               ( -2520.0 * f( l1-1:l2-1 , m , n + 3 , ik ))+&
               ( 1755.0 * f( l1-1:l2-1 , m - 1 , n + 3 , ik ))+&
               ( -540.0 * f( l1-1:l2-1 , m - 2 , n + 3 , ik ))+&
               ( 45.0 * f( l1-1:l2-1 , m - 3 , n + 3 , ik ))+&
               ( 1755.0 * f( l1-1:l2-1 , m + 1 , n + 3 , ik ))+&
               ( -540.0 * f( l1-1:l2-1 , m + 2 , n + 3 , ik ))+&
               ( 45.0 * f( l1-1:l2-1 , m + 3 , n + 3 , ik ))+&
               ( 504.0 * f( l1-2:l2-2 , m , n + 3 , ik ))+&
               ( -351.0 * f( l1-2:l2-2 , m - 1 , n + 3 , ik ))+&
               ( 108.0 * f( l1-2:l2-2 , m - 2 , n + 3 , ik ))+&
               ( -9.0 * f( l1-2:l2-2 , m - 3 , n + 3 , ik ))+&
               ( -351.0 * f( l1-2:l2-2 , m + 1 , n + 3 , ik ))+&
               ( 108.0 * f( l1-2:l2-2 , m + 2 , n + 3 , ik ))+&
               ( -9.0 * f( l1-2:l2-2 , m + 3 , n + 3 , ik ))+&
               ( -56.0 * f( l1-3:l2-3 , m , n + 3 , ik ))+&
               ( 39.0 * f( l1-3:l2-3 , m - 1 , n + 3 , ik ))+&
               ( -12.0 * f( l1-3:l2-3 , m - 2 , n + 3 , ik ))+&
               ( 1.0 * f( l1-3:l2-3 , m - 3 , n + 3 , ik ))+&
               ( 39.0 * f( l1-3:l2-3 , m + 1 , n + 3 , ik ))+&
               ( -12.0 * f( l1-3:l2-3 , m + 2 , n + 3 , ik ))+&
               ( 1.0 * f( l1-3:l2-3 , m + 3 , n + 3 , ik ))+&
               ( 2520.0 * f( l1+1:l2+1 , m , n + 3 , ik ))+&
               ( -1755.0 * f( l1+1:l2+1 , m - 1 , n + 3 , ik ))+&
               ( 540.0 * f( l1+1:l2+1 , m - 2 , n + 3 , ik ))+&
               ( -45.0 * f( l1+1:l2+1 , m - 3 , n + 3 , ik ))+&
               ( -1755.0 * f( l1+1:l2+1 , m + 1 , n + 3 , ik ))+&
               ( 540.0 * f( l1+1:l2+1 , m + 2 , n + 3 , ik ))+&
               ( -45.0 * f( l1+1:l2+1 , m + 3 , n + 3 , ik ))+&
               ( -504.0 * f( l1+2:l2+2 , m , n + 3 , ik ))+&
               ( 351.0 * f( l1+2:l2+2 , m - 1 , n + 3 , ik ))+&
               ( -108.0 * f( l1+2:l2+2 , m - 2 , n + 3 , ik ))+&
               ( 9.0 * f( l1+2:l2+2 , m - 3 , n + 3 , ik ))+&
               ( 351.0 * f( l1+2:l2+2 , m + 1 , n + 3 , ik ))+&
               ( -108.0 * f( l1+2:l2+2 , m + 2 , n + 3 , ik ))+&
               ( 9.0 * f( l1+2:l2+2 , m + 3 , n + 3 , ik ))+&
               ( 56.0 * f( l1+3:l2+3 , m , n + 3 , ik ))+&
               ( -39.0 * f( l1+3:l2+3 , m - 1 , n + 3 , ik ))+&
               ( 12.0 * f( l1+3:l2+3 , m - 2 , n + 3 , ik ))+&
               ( -1.0 * f( l1+3:l2+3 , m - 3 , n + 3 , ik ))+&
               ( -39.0 * f( l1+3:l2+3 , m + 1 , n + 3 , ik ))+&
               ( 12.0 * f( l1+3:l2+3 , m + 2 , n + 3 , ik ))+&
               ( -1.0 * f( l1+3:l2+3 , m + 3 , n + 3 , ik ))&
               )
        else
          df=0.0
        endif yxz
!        
      elseif ((i==3.and.j==1.and.k==2).or.(i==3.and.j==2.and.k==1)) then
        zxy: if (nxgrid/=1.and.nygrid/=1.and.nzgrid/=1) then
          fac= (1./6.0*dz_1(n)**4) * (1./60.0*dx_1(l1:l2)) * (1/60.0*dy_1(m))
          df = fac*(&               
               ( 113400.0 * f( l1-1:l2-1 , m - 1 , n , ik ))+&
               ( -78975.0 * f( l1-1:l2-1 , m - 1 , n - 1 , ik ))+&
               ( 24300.0 * f( l1-1:l2-1 , m - 1 , n - 2 , ik ))+&
               ( -2025.0 * f( l1-1:l2-1 , m - 1 , n - 3 , ik ))+&
               ( -78975.0 * f( l1-1:l2-1 , m - 1 , n + 1 , ik ))+&
               ( 24300.0 * f( l1-1:l2-1 , m - 1 , n + 2 , ik ))+&
               ( -2025.0 * f( l1-1:l2-1 , m - 1 , n + 3 , ik ))+&
               ( -22680.0 * f( l1-2:l2-2 , m - 1 , n , ik ))+&
               ( 15795.0 * f( l1-2:l2-2 , m - 1 , n - 1 , ik ))+&
               ( -4860.0 * f( l1-2:l2-2 , m - 1 , n - 2 , ik ))+&
               ( 405.0 * f( l1-2:l2-2 , m - 1 , n - 3 , ik ))+&
               ( 15795.0 * f( l1-2:l2-2 , m - 1 , n + 1 , ik ))+&
               ( -4860.0 * f( l1-2:l2-2 , m - 1 , n + 2 , ik ))+&
               ( 405.0 * f( l1-2:l2-2 , m - 1 , n + 3 , ik ))+&
               ( 2520.0 * f( l1-3:l2-3 , m - 1 , n , ik ))+&
               ( -1755.0 * f( l1-3:l2-3 , m - 1 , n - 1 , ik ))+&
               ( 540.0 * f( l1-3:l2-3 , m - 1 , n - 2 , ik ))+&
               ( -45.0 * f( l1-3:l2-3 , m - 1 , n - 3 , ik ))+&
               ( -1755.0 * f( l1-3:l2-3 , m - 1 , n + 1 , ik ))+&
               ( 540.0 * f( l1-3:l2-3 , m - 1 , n + 2 , ik ))+&
               ( -45.0 * f( l1-3:l2-3 , m - 1 , n + 3 , ik ))+&
               ( -113400.0 * f( l1+1:l2+1 , m - 1 , n , ik ))+&
               ( 78975.0 * f( l1+1:l2+1 , m - 1 , n - 1 , ik ))+&
               ( -24300.0 * f( l1+1:l2+1 , m - 1 , n - 2 , ik ))+&
               ( 2025.0 * f( l1+1:l2+1 , m - 1 , n - 3 , ik ))+&
               ( 78975.0 * f( l1+1:l2+1 , m - 1 , n + 1 , ik ))+&
               ( -24300.0 * f( l1+1:l2+1 , m - 1 , n + 2 , ik ))+&
               ( 2025.0 * f( l1+1:l2+1 , m - 1 , n + 3 , ik ))+&
               ( 22680.0 * f( l1+2:l2+2 , m - 1 , n , ik ))+&
               ( -15795.0 * f( l1+2:l2+2 , m - 1 , n - 1 , ik ))+&
               ( 4860.0 * f( l1+2:l2+2 , m - 1 , n - 2 , ik ))+&
               ( -405.0 * f( l1+2:l2+2 , m - 1 , n - 3 , ik ))+&
               ( -15795.0 * f( l1+2:l2+2 , m - 1 , n + 1 , ik ))+&
               ( 4860.0 * f( l1+2:l2+2 , m - 1 , n + 2 , ik ))+&
               ( -405.0 * f( l1+2:l2+2 , m - 1 , n + 3 , ik ))+&
               ( -2520.0 * f( l1+3:l2+3 , m - 1 , n , ik ))+&
               ( 1755.0 * f( l1+3:l2+3 , m - 1 , n - 1 , ik ))+&
               ( -540.0 * f( l1+3:l2+3 , m - 1 , n - 2 , ik ))+&
               ( 45.0 * f( l1+3:l2+3 , m - 1 , n - 3 , ik ))+&
               ( 1755.0 * f( l1+3:l2+3 , m - 1 , n + 1 , ik ))+&
               ( -540.0 * f( l1+3:l2+3 , m - 1 , n + 2 , ik ))+&
               ( 45.0 * f( l1+3:l2+3 , m - 1 , n + 3 , ik ))+&
               ( -22680.0 * f( l1-1:l2-1 , m - 2 , n , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 2 , n - 1 , ik ))+&
               ( -4860.0 * f( l1-1:l2-1 , m - 2 , n - 2 , ik ))+&
               ( 405.0 * f( l1-1:l2-1 , m - 2 , n - 3 , ik ))+&
               ( 15795.0 * f( l1-1:l2-1 , m - 2 , n + 1 , ik ))+&
               ( -4860.0 * f( l1-1:l2-1 , m - 2 , n + 2 , ik ))+&
               ( 405.0 * f( l1-1:l2-1 , m - 2 , n + 3 , ik ))+&
               ( 4536.0 * f( l1-2:l2-2 , m - 2 , n , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m - 2 , n - 1 , ik ))+&
               ( 972.0 * f( l1-2:l2-2 , m - 2 , n - 2 , ik ))+&
               ( -81.0 * f( l1-2:l2-2 , m - 2 , n - 3 , ik ))+&
               ( -3159.0 * f( l1-2:l2-2 , m - 2 , n + 1 , ik ))+&
               ( 972.0 * f( l1-2:l2-2 , m - 2 , n + 2 , ik ))+&
               ( -81.0 * f( l1-2:l2-2 , m - 2 , n + 3 , ik ))+&
               ( -504.0 * f( l1-3:l2-3 , m - 2 , n , ik ))+&
               ( 351.0 * f( l1-3:l2-3 , m - 2 , n - 1 , ik ))+&
               ( -108.0 * f( l1-3:l2-3 , m - 2 , n - 2 , ik ))+&
               ( 9.0 * f( l1-3:l2-3 , m - 2 , n - 3 , ik ))+&
               ( 351.0 * f( l1-3:l2-3 , m - 2 , n + 1 , ik ))+&
               ( -108.0 * f( l1-3:l2-3 , m - 2 , n + 2 , ik ))+&
               ( 9.0 * f( l1-3:l2-3 , m - 2 , n + 3 , ik ))+&
               ( 22680.0 * f( l1+1:l2+1 , m - 2 , n , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 2 , n - 1 , ik ))+&
               ( 4860.0 * f( l1+1:l2+1 , m - 2 , n - 2 , ik ))+&
               ( -405.0 * f( l1+1:l2+1 , m - 2 , n - 3 , ik ))+&
               ( -15795.0 * f( l1+1:l2+1 , m - 2 , n + 1 , ik ))+&
               ( 4860.0 * f( l1+1:l2+1 , m - 2 , n + 2 , ik ))+&
               ( -405.0 * f( l1+1:l2+1 , m - 2 , n + 3 , ik ))+&
               ( -4536.0 * f( l1+2:l2+2 , m - 2 , n , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m - 2 , n - 1 , ik ))+&
               ( -972.0 * f( l1+2:l2+2 , m - 2 , n - 2 , ik ))+&
               ( 81.0 * f( l1+2:l2+2 , m - 2 , n - 3 , ik ))+&
               ( 3159.0 * f( l1+2:l2+2 , m - 2 , n + 1 , ik ))+&
               ( -972.0 * f( l1+2:l2+2 , m - 2 , n + 2 , ik ))+&
               ( 81.0 * f( l1+2:l2+2 , m - 2 , n + 3 , ik ))+&
               ( 504.0 * f( l1+3:l2+3 , m - 2 , n , ik ))+&
               ( -351.0 * f( l1+3:l2+3 , m - 2 , n - 1 , ik ))+&
               ( 108.0 * f( l1+3:l2+3 , m - 2 , n - 2 , ik ))+&
               ( -9.0 * f( l1+3:l2+3 , m - 2 , n - 3 , ik ))+&
               ( -351.0 * f( l1+3:l2+3 , m - 2 , n + 1 , ik ))+&
               ( 108.0 * f( l1+3:l2+3 , m - 2 , n + 2 , ik ))+&
               ( -9.0 * f( l1+3:l2+3 , m - 2 , n + 3 , ik ))+&
               ( 2520.0 * f( l1-1:l2-1 , m - 3 , n , ik ))+&
               ( -1755.0 * f( l1-1:l2-1 , m - 3 , n - 1 , ik ))+&
               ( 540.0 * f( l1-1:l2-1 , m - 3 , n - 2 , ik ))+&
               ( -45.0 * f( l1-1:l2-1 , m - 3 , n - 3 , ik ))+&
               ( -1755.0 * f( l1-1:l2-1 , m - 3 , n + 1 , ik ))+&
               ( 540.0 * f( l1-1:l2-1 , m - 3 , n + 2 , ik ))+&
               ( -45.0 * f( l1-1:l2-1 , m - 3 , n + 3 , ik ))+&
               ( -504.0 * f( l1-2:l2-2 , m - 3 , n , ik ))+&
               ( 351.0 * f( l1-2:l2-2 , m - 3 , n - 1 , ik ))+&
               ( -108.0 * f( l1-2:l2-2 , m - 3 , n - 2 , ik ))+&
               ( 9.0 * f( l1-2:l2-2 , m - 3 , n - 3 , ik ))+&
               ( 351.0 * f( l1-2:l2-2 , m - 3 , n + 1 , ik ))+&
               ( -108.0 * f( l1-2:l2-2 , m - 3 , n + 2 , ik ))+&
               ( 9.0 * f( l1-2:l2-2 , m - 3 , n + 3 , ik ))+&
               ( 56.0 * f( l1-3:l2-3 , m - 3 , n , ik ))+&
               ( -39.0 * f( l1-3:l2-3 , m - 3 , n - 1 , ik ))+&
               ( 12.0 * f( l1-3:l2-3 , m - 3 , n - 2 , ik ))+&
               ( -1.0 * f( l1-3:l2-3 , m - 3 , n - 3 , ik ))+&
               ( -39.0 * f( l1-3:l2-3 , m - 3 , n + 1 , ik ))+&
               ( 12.0 * f( l1-3:l2-3 , m - 3 , n + 2 , ik ))+&
               ( -1.0 * f( l1-3:l2-3 , m - 3 , n + 3 , ik ))+&
               ( -2520.0 * f( l1+1:l2+1 , m - 3 , n , ik ))+&
               ( 1755.0 * f( l1+1:l2+1 , m - 3 , n - 1 , ik ))+&
               ( -540.0 * f( l1+1:l2+1 , m - 3 , n - 2 , ik ))+&
               ( 45.0 * f( l1+1:l2+1 , m - 3 , n - 3 , ik ))+&
               ( 1755.0 * f( l1+1:l2+1 , m - 3 , n + 1 , ik ))+&
               ( -540.0 * f( l1+1:l2+1 , m - 3 , n + 2 , ik ))+&
               ( 45.0 * f( l1+1:l2+1 , m - 3 , n + 3 , ik ))+&
               ( 504.0 * f( l1+2:l2+2 , m - 3 , n , ik ))+&
               ( -351.0 * f( l1+2:l2+2 , m - 3 , n - 1 , ik ))+&
               ( 108.0 * f( l1+2:l2+2 , m - 3 , n - 2 , ik ))+&
               ( -9.0 * f( l1+2:l2+2 , m - 3 , n - 3 , ik ))+&
               ( -351.0 * f( l1+2:l2+2 , m - 3 , n + 1 , ik ))+&
               ( 108.0 * f( l1+2:l2+2 , m - 3 , n + 2 , ik ))+&
               ( -9.0 * f( l1+2:l2+2 , m - 3 , n + 3 , ik ))+&
               ( -56.0 * f( l1+3:l2+3 , m - 3 , n , ik ))+&
               ( 39.0 * f( l1+3:l2+3 , m - 3 , n - 1 , ik ))+&
               ( -12.0 * f( l1+3:l2+3 , m - 3 , n - 2 , ik ))+&
               ( 1.0 * f( l1+3:l2+3 , m - 3 , n - 3 , ik ))+&
               ( 39.0 * f( l1+3:l2+3 , m - 3 , n + 1 , ik ))+&
               ( -12.0 * f( l1+3:l2+3 , m - 3 , n + 2 , ik ))+&
               ( 1.0 * f( l1+3:l2+3 , m - 3 , n + 3 , ik ))+&
               ( -113400.0 * f( l1-1:l2-1 , m + 1 , n , ik ))+&
               ( 78975.0 * f( l1-1:l2-1 , m + 1 , n - 1 , ik ))+&
               ( -24300.0 * f( l1-1:l2-1 , m + 1 , n - 2 , ik ))+&
               ( 2025.0 * f( l1-1:l2-1 , m + 1 , n - 3 , ik ))+&
               ( 78975.0 * f( l1-1:l2-1 , m + 1 , n + 1 , ik ))+&
               ( -24300.0 * f( l1-1:l2-1 , m + 1 , n + 2 , ik ))+&
               ( 2025.0 * f( l1-1:l2-1 , m + 1 , n + 3 , ik ))+&
               ( 22680.0 * f( l1-2:l2-2 , m + 1 , n , ik ))+&
               ( -15795.0 * f( l1-2:l2-2 , m + 1 , n - 1 , ik ))+&
               ( 4860.0 * f( l1-2:l2-2 , m + 1 , n - 2 , ik ))+&
               ( -405.0 * f( l1-2:l2-2 , m + 1 , n - 3 , ik ))+&
               ( -15795.0 * f( l1-2:l2-2 , m + 1 , n + 1 , ik ))+&
               ( 4860.0 * f( l1-2:l2-2 , m + 1 , n + 2 , ik ))+&
               ( -405.0 * f( l1-2:l2-2 , m + 1 , n + 3 , ik ))+&
               ( -2520.0 * f( l1-3:l2-3 , m + 1 , n , ik ))+&
               ( 1755.0 * f( l1-3:l2-3 , m + 1 , n - 1 , ik ))+&
               ( -540.0 * f( l1-3:l2-3 , m + 1 , n - 2 , ik ))+&
               ( 45.0 * f( l1-3:l2-3 , m + 1 , n - 3 , ik ))+&
               ( 1755.0 * f( l1-3:l2-3 , m + 1 , n + 1 , ik ))+&
               ( -540.0 * f( l1-3:l2-3 , m + 1 , n + 2 , ik ))+&
               ( 45.0 * f( l1-3:l2-3 , m + 1 , n + 3 , ik ))+&
               ( 113400.0 * f( l1+1:l2+1 , m + 1 , n , ik ))+&
               ( -78975.0 * f( l1+1:l2+1 , m + 1 , n - 1 , ik ))+&
               ( 24300.0 * f( l1+1:l2+1 , m + 1 , n - 2 , ik ))+&
               ( -2025.0 * f( l1+1:l2+1 , m + 1 , n - 3 , ik ))+&
               ( -78975.0 * f( l1+1:l2+1 , m + 1 , n + 1 , ik ))+&
               ( 24300.0 * f( l1+1:l2+1 , m + 1 , n + 2 , ik ))+&
               ( -2025.0 * f( l1+1:l2+1 , m + 1 , n + 3 , ik ))+&
               ( -22680.0 * f( l1+2:l2+2 , m + 1 , n , ik ))+&
               ( 15795.0 * f( l1+2:l2+2 , m + 1 , n - 1 , ik ))+&
               ( -4860.0 * f( l1+2:l2+2 , m + 1 , n - 2 , ik ))+&
               ( 405.0 * f( l1+2:l2+2 , m + 1 , n - 3 , ik ))+&
               ( 15795.0 * f( l1+2:l2+2 , m + 1 , n + 1 , ik ))+&
               ( -4860.0 * f( l1+2:l2+2 , m + 1 , n + 2 , ik ))+&
               ( 405.0 * f( l1+2:l2+2 , m + 1 , n + 3 , ik ))+&
               ( 2520.0 * f( l1+3:l2+3 , m + 1 , n , ik ))+&
               ( -1755.0 * f( l1+3:l2+3 , m + 1 , n - 1 , ik ))+&
               ( 540.0 * f( l1+3:l2+3 , m + 1 , n - 2 , ik ))+&
               ( -45.0 * f( l1+3:l2+3 , m + 1 , n - 3 , ik ))+&
               ( -1755.0 * f( l1+3:l2+3 , m + 1 , n + 1 , ik ))+&
               ( 540.0 * f( l1+3:l2+3 , m + 1 , n + 2 , ik ))+&
               ( -45.0 * f( l1+3:l2+3 , m + 1 , n + 3 , ik ))+&
               ( 22680.0 * f( l1-1:l2-1 , m + 2 , n , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 2 , n - 1 , ik ))+&
               ( 4860.0 * f( l1-1:l2-1 , m + 2 , n - 2 , ik ))+&
               ( -405.0 * f( l1-1:l2-1 , m + 2 , n - 3 , ik ))+&
               ( -15795.0 * f( l1-1:l2-1 , m + 2 , n + 1 , ik ))+&
               ( 4860.0 * f( l1-1:l2-1 , m + 2 , n + 2 , ik ))+&
               ( -405.0 * f( l1-1:l2-1 , m + 2 , n + 3 , ik ))+&
               ( -4536.0 * f( l1-2:l2-2 , m + 2 , n , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m + 2 , n - 1 , ik ))+&
               ( -972.0 * f( l1-2:l2-2 , m + 2 , n - 2 , ik ))+&
               ( 81.0 * f( l1-2:l2-2 , m + 2 , n - 3 , ik ))+&
               ( 3159.0 * f( l1-2:l2-2 , m + 2 , n + 1 , ik ))+&
               ( -972.0 * f( l1-2:l2-2 , m + 2 , n + 2 , ik ))+&
               ( 81.0 * f( l1-2:l2-2 , m + 2 , n + 3 , ik ))+&
               ( 504.0 * f( l1-3:l2-3 , m + 2 , n , ik ))+&
               ( -351.0 * f( l1-3:l2-3 , m + 2 , n - 1 , ik ))+&
               ( 108.0 * f( l1-3:l2-3 , m + 2 , n - 2 , ik ))+&
               ( -9.0 * f( l1-3:l2-3 , m + 2 , n - 3 , ik ))+&
               ( -351.0 * f( l1-3:l2-3 , m + 2 , n + 1 , ik ))+&
               ( 108.0 * f( l1-3:l2-3 , m + 2 , n + 2 , ik ))+&
               ( -9.0 * f( l1-3:l2-3 , m + 2 , n + 3 , ik ))+&
               ( -22680.0 * f( l1+1:l2+1 , m + 2 , n , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 2 , n - 1 , ik ))+&
               ( -4860.0 * f( l1+1:l2+1 , m + 2 , n - 2 , ik ))+&
               ( 405.0 * f( l1+1:l2+1 , m + 2 , n - 3 , ik ))+&
               ( 15795.0 * f( l1+1:l2+1 , m + 2 , n + 1 , ik ))+&
               ( -4860.0 * f( l1+1:l2+1 , m + 2 , n + 2 , ik ))+&
               ( 405.0 * f( l1+1:l2+1 , m + 2 , n + 3 , ik ))+&
               ( 4536.0 * f( l1+2:l2+2 , m + 2 , n , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m + 2 , n - 1 , ik ))+&
               ( 972.0 * f( l1+2:l2+2 , m + 2 , n - 2 , ik ))+&
               ( -81.0 * f( l1+2:l2+2 , m + 2 , n - 3 , ik ))+&
               ( -3159.0 * f( l1+2:l2+2 , m + 2 , n + 1 , ik ))+&
               ( 972.0 * f( l1+2:l2+2 , m + 2 , n + 2 , ik ))+&
               ( -81.0 * f( l1+2:l2+2 , m + 2 , n + 3 , ik ))+&
               ( -504.0 * f( l1+3:l2+3 , m + 2 , n , ik ))+&
               ( 351.0 * f( l1+3:l2+3 , m + 2 , n - 1 , ik ))+&
               ( -108.0 * f( l1+3:l2+3 , m + 2 , n - 2 , ik ))+&
               ( 9.0 * f( l1+3:l2+3 , m + 2 , n - 3 , ik ))+&
               ( 351.0 * f( l1+3:l2+3 , m + 2 , n + 1 , ik ))+&
               ( -108.0 * f( l1+3:l2+3 , m + 2 , n + 2 , ik ))+&
               ( 9.0 * f( l1+3:l2+3 , m + 2 , n + 3 , ik ))+&
               ( -2520.0 * f( l1-1:l2-1 , m + 3 , n , ik ))+&
               ( 1755.0 * f( l1-1:l2-1 , m + 3 , n - 1 , ik ))+&
               ( -540.0 * f( l1-1:l2-1 , m + 3 , n - 2 , ik ))+&
               ( 45.0 * f( l1-1:l2-1 , m + 3 , n - 3 , ik ))+&
               ( 1755.0 * f( l1-1:l2-1 , m + 3 , n + 1 , ik ))+&
               ( -540.0 * f( l1-1:l2-1 , m + 3 , n + 2 , ik ))+&
               ( 45.0 * f( l1-1:l2-1 , m + 3 , n + 3 , ik ))+&
               ( 504.0 * f( l1-2:l2-2 , m + 3 , n , ik ))+&
               ( -351.0 * f( l1-2:l2-2 , m + 3 , n - 1 , ik ))+&
               ( 108.0 * f( l1-2:l2-2 , m + 3 , n - 2 , ik ))+&
               ( -9.0 * f( l1-2:l2-2 , m + 3 , n - 3 , ik ))+&
               ( -351.0 * f( l1-2:l2-2 , m + 3 , n + 1 , ik ))+&
               ( 108.0 * f( l1-2:l2-2 , m + 3 , n + 2 , ik ))+&
               ( -9.0 * f( l1-2:l2-2 , m + 3 , n + 3 , ik ))+&
               ( -56.0 * f( l1-3:l2-3 , m + 3 , n , ik ))+&
               ( 39.0 * f( l1-3:l2-3 , m + 3 , n - 1 , ik ))+&
               ( -12.0 * f( l1-3:l2-3 , m + 3 , n - 2 , ik ))+&
               ( 1.0 * f( l1-3:l2-3 , m + 3 , n - 3 , ik ))+&
               ( 39.0 * f( l1-3:l2-3 , m + 3 , n + 1 , ik ))+&
               ( -12.0 * f( l1-3:l2-3 , m + 3 , n + 2 , ik ))+&
               ( 1.0 * f( l1-3:l2-3 , m + 3 , n + 3 , ik ))+&
               ( 2520.0 * f( l1+1:l2+1 , m + 3 , n , ik ))+&
               ( -1755.0 * f( l1+1:l2+1 , m + 3 , n - 1 , ik ))+&
               ( 540.0 * f( l1+1:l2+1 , m + 3 , n - 2 , ik ))+&
               ( -45.0 * f( l1+1:l2+1 , m + 3 , n - 3 , ik ))+&
               ( -1755.0 * f( l1+1:l2+1 , m + 3 , n + 1 , ik ))+&
               ( 540.0 * f( l1+1:l2+1 , m + 3 , n + 2 , ik ))+&
               ( -45.0 * f( l1+1:l2+1 , m + 3 , n + 3 , ik ))+&
               ( -504.0 * f( l1+2:l2+2 , m + 3 , n , ik ))+&
               ( 351.0 * f( l1+2:l2+2 , m + 3 , n - 1 , ik ))+&
               ( -108.0 * f( l1+2:l2+2 , m + 3 , n - 2 , ik ))+&
               ( 9.0 * f( l1+2:l2+2 , m + 3 , n - 3 , ik ))+&
               ( 351.0 * f( l1+2:l2+2 , m + 3 , n + 1 , ik ))+&
               ( -108.0 * f( l1+2:l2+2 , m + 3 , n + 2 , ik ))+&
               ( 9.0 * f( l1+2:l2+2 , m + 3 , n + 3 , ik ))+&
               ( 56.0 * f( l1+3:l2+3 , m + 3 , n , ik ))+&
               ( -39.0 * f( l1+3:l2+3 , m + 3 , n - 1 , ik ))+&
               ( 12.0 * f( l1+3:l2+3 , m + 3 , n - 2 , ik ))+&
               ( -1.0 * f( l1+3:l2+3 , m + 3 , n - 3 , ik ))+&
               ( -39.0 * f( l1+3:l2+3 , m + 3 , n + 1 , ik ))+&
               ( 12.0 * f( l1+3:l2+3 , m + 3 , n + 2 , ik ))+&
               ( -1.0 * f( l1+3:l2+3 , m + 3 , n + 3 , ik ))&
               )
        else
          df=0.0
        endif zxy
!        
      endif
!      
    endsubroutine der4i1j1k
!***********************************************************************
    subroutine der_upwind1st(f,uu,k,df,j)
!
!  First order upwind derivative of variable
!
!  Useful for advecting non-logarithmic variables
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: df
      integer :: j,k,l
!
      intent(in)  :: f,uu,k,j
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der_upwind1st,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der_upwind1st,j,1) + 1 !DERCOUNT
!
      if (lspherical_coords.or.lcylindrical_coords) &
          call fatal_error('der_upwind1st','NOT IMPLEMENTED for non-cartesian grid')
!
      if (j == 1) then
        if (nxgrid /= 1) then
          do l=1,nx
            if (uu(l,1) > 0.) then
              df(l) = (f(nghost+l  ,m,n,k) - f(nghost+l-1,m,n,k))*dx_1(nghost+l)
            else
              df(l) = (f(nghost+l+1,m,n,k) - f(nghost+l  ,m,n,k))*dx_1(nghost+l)
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
              df(l) = (f(nghost+l,m  ,n,k) - f(nghost+l,m-1,n,k))*dy_1(m)
            else
              df(l) = (f(nghost+l,m+1,n,k) - f(nghost+l,m  ,n,k))*dy_1(m)
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
              df(l) = (f(nghost+l,m,n,k  ) - f(nghost+l,m,n-1,k))*dz_1(n)
            else
              df(l) = (f(nghost+l,m,n+1,k) - f(nghost+l,m,n,k  ))*dz_1(n)
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
    subroutine der_onesided_4_slice_main(f,sgn,k,df,pos,j)
!
!   Calculate x/y/z-derivative on a yz/xz/xy-slice at gridpoint pos.
!   Uses a one-sided 4th order stencil.
!   sgn = +1 for forward difference, sgn = -1 for backwards difference.
!
!   Because of its original intended use in relation to solving
!   characteristic equations on boundaries (NSCBC), this sub should
!   return only PARTIAL derivatives, NOT COVARIANT. Applying the right
!   scaling factors and connection terms should instead be done when
!   solving the characteristic equations.
!
!   7-jul-08/arne: coded.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:,:) :: df
      real :: fac
      integer :: pos,k,sgn,j
!
      intent(in)  :: f,k,pos,sgn,j
      intent(out) :: df
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=1./12.*dx_1(pos)
          df = fac*(-sgn*25*f(pos,m1:m2,n1:n2,k)&
                  +sgn*48*f(pos+sgn*1,m1:m2,n1:n2,k)&
                  -sgn*36*f(pos+sgn*2,m1:m2,n1:n2,k)&
                  +sgn*16*f(pos+sgn*3,m1:m2,n1:n2,k)&
                  -sgn*3 *f(pos+sgn*4,m1:m2,n1:n2,k))
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in x-directder_onesided_4_sliceion'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=1./12.*dy_1(pos)
          df = fac*(-sgn*25*f(l1:l2,pos,n1:n2,k)&
                  +sgn*48*f(l1:l2,pos+sgn*1,n1:n2,k)&
                  -sgn*36*f(l1:l2,pos+sgn*2,n1:n2,k)&
                  +sgn*16*f(l1:l2,pos+sgn*3,n1:n2,k)&
                  -sgn*3 *f(l1:l2,pos+sgn*4,n1:n2,k))
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=1./12.*dz_1(pos)
          df = fac*(-sgn*25*f(l1:l2,m1:m2,pos,k)&
                  +sgn*48*f(l1:l2,m1:m2,pos+sgn*1,k)&
                  -sgn*36*f(l1:l2,m1:m2,pos+sgn*2,k)&
                  +sgn*16*f(l1:l2,m1:m2,pos+sgn*3,k)&
                  -sgn*3 *f(l1:l2,m1:m2,pos+sgn*4,k))
!MR: coarse missing
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in z-direction'
        endif
      endif
    endsubroutine der_onesided_4_slice_main
!***********************************************************************
    subroutine der_onesided_4_slice_main_pt(f,sgn,k,df,lll,mmm,nnn,j)
!
!  made using der_onesided_4_slice_main. One sided derivative is calculated
!  at one point (lll,mmm,nnn)
!
!  15-oct-09/Natalia: coded.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real  :: df
      real :: fac
      integer :: pos,lll,mmm,nnn,k,sgn,j
!
      intent(in)  :: f,k,lll,mmm,nnn,sgn,j
      intent(out) :: df
!
      if (j==1) then
       pos=lll
        if (nxgrid/=1) then
          fac=1./12.*dx_1(pos)
          df = fac*(-sgn*25*f(pos,mmm,nnn,k)&
                  +sgn*48*f(pos+sgn*1,mmm,nnn,k)&
                  -sgn*36*f(pos+sgn*2,mmm,nnn,k)&
                  +sgn*16*f(pos+sgn*3,mmm,nnn,k)&
                  -sgn*3 *f(pos+sgn*4,mmm,nnn,k))
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in x-directder_onesided_4_sliceion'
        endif
      elseif (j==2) then
       pos=mmm
        if (nygrid/=1) then
          fac=1./12.*dy_1(pos)
          df = fac*(-sgn*25*f(lll,pos,nnn,k)&
                  +sgn*48*f(lll,pos+sgn*1,nnn,k)&
                  -sgn*36*f(lll,pos+sgn*2,nnn,k)&
                  +sgn*16*f(lll,pos+sgn*3,nnn,k)&
                  -sgn*3 *f(lll,pos+sgn*4,nnn,k))
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in y-direction'
        endif
      elseif (j==3) then
       pos=nnn
        if (nzgrid/=1) then
          fac=1./12.*dz_1(pos)
          df = fac*(-sgn*25*f(lll,mmm,pos,k)&
                  +sgn*48*f(lll,mmm,pos+sgn*1,k)&
                  -sgn*36*f(lll,mmm,pos+sgn*2,k)&
                  +sgn*16*f(lll,mmm,pos+sgn*3,k)&
                  -sgn*3 *f(lll,mmm,pos+sgn*4,k))
!MR: coarse missing
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in z-direction'
        endif
      endif
    endsubroutine der_onesided_4_slice_main_pt
!***********************************************************************
    subroutine der_onesided_4_slice_other(f,sgn,df,pos,j)
!
!   Calculate x/y/z-derivative on a yz/xz/xy-slice at gridpoint pos.
!   Uses a one-sided 4th order stencil.
!   sgn = +1 for forward difference, sgn = -1 for backwards difference.
!
!   Because of its original intended use in relation to solving
!   characteristic equations on boundaries (NSCBC), this sub should
!   return only PARTIAL derivatives, NOT COVARIANT. Applying the right
!   scaling factors and connection terms should instead be done when
!   solving the characteristic equations.
!
!   7-jul-08/arne: coded.
!
      real, dimension (mx,my,mz) :: f
      real, dimension (:,:) :: df
      real :: fac
      integer :: pos,sgn,j
!
      intent(in)  :: f,pos,sgn,j
      intent(out) :: df
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=1./12.*dx_1(pos)
          df = fac*(-sgn*25*f(pos,m1:m2,n1:n2)&
                  +sgn*48*f(pos+sgn*1,m1:m2,n1:n2)&
                  -sgn*36*f(pos+sgn*2,m1:m2,n1:n2)&
                  +sgn*16*f(pos+sgn*3,m1:m2,n1:n2)&
                  -sgn*3 *f(pos+sgn*4,m1:m2,n1:n2))
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in x-directder_onesided_4_sliceion'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=1./12.*dy_1(pos)
          df = fac*(-sgn*25*(f(l1:l2,pos,n1:n2)-f(l1:l2,pos+sgn*1,n1:n2))&
                    +sgn*23*(f(l1:l2,pos+sgn*1,n1:n2)-f(l1:l2,pos+sgn*2,n1:n2))&
                    -sgn*13*(f(l1:l2,pos+sgn*2,n1:n2)-f(l1:l2,pos+sgn*3,n1:n2))&
                    +sgn*3*(f(l1:l2,pos+sgn*3,n1:n2)-f(l1:l2,pos+sgn*4,n1:n2)))
!
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=1./12.*dz_1(pos)
          df = fac*(-sgn*25*f(l1:l2,m1:m2,pos)&
                    +sgn*48*f(l1:l2,m1:m2,pos+sgn*1)&
                    -sgn*36*f(l1:l2,m1:m2,pos+sgn*2)&
                    +sgn*16*f(l1:l2,m1:m2,pos+sgn*3)&
                    -sgn*3 *f(l1:l2,m1:m2,pos+sgn*4))
!MR: coarse missing
!
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine der_onesided_4_slice_other
!***********************************************************************
    subroutine der_onesided_4_slice_other_pt(f,sgn,df,lll,mmm,nnn,j)
!
!  One sided derivative is calculated
!  at one point (lll,mmm,nnn).
!
!  15-oct-09/Natalia: coded.
!  15-oct-09/axel: changed file name to shorter version
!
      real, dimension (mx,my,mz) :: f
      real :: df
      real :: fac
      integer :: pos,lll,mmm,nnn,sgn,j
!
      intent(in)  :: f,lll,mmm,nnn,sgn,j
      intent(out) :: df
!
      if (j==1) then
        pos=lll
        if (nxgrid/=1) then
          fac=1./12.*dx_1(pos)
          df = fac*(-sgn*25*f(pos,mmm,nnn)&
                    +sgn*48*f(pos+sgn*1,mmm,nnn)&
                    -sgn*36*f(pos+sgn*2,mmm,nnn)&
                    +sgn*16*f(pos+sgn*3,mmm,nnn)&
                    -sgn*3 *f(pos+sgn*4,mmm,nnn))
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice_other_pt: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        pos=mmm
        if (nygrid/=1) then
          fac=1./12.*dy_1(pos)
          df = fac*(-sgn*25*(f(lll,pos,nnn)-f(lll,pos+sgn*1,nnn))&
                    +sgn*23*(f(lll,pos+sgn*1,nnn)-f(lll,pos+sgn*2,nnn))&
                    -sgn*13*(f(lll,pos+sgn*2,nnn)-f(lll,pos+sgn*3,nnn))&
                    +sgn* 3*(f(lll,pos+sgn*3,nnn)-f(lll,pos+sgn*4,nnn)))
!
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice_other_pt: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        pos=nnn
        if (nzgrid/=1) then
          fac=1./12.*dz_1(pos)
          df = fac*(-sgn*25*f(lll,mmm,pos)&
                    +sgn*48*f(lll,mmm,pos+sgn*1)&
                    -sgn*36*f(lll,mmm,pos+sgn*2)&
                    +sgn*16*f(lll,mmm,pos+sgn*3)&
                    -sgn*3 *f(lll,mmm,pos+sgn*4))
!
!MR: coarse missing
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice_other_pt: Degenerate case in z-direction'
        endif
      endif
    endsubroutine der_onesided_4_slice_other_pt
!***********************************************************************
    subroutine finalize_deriv
!
!  Dummy
!
    endsubroutine finalize_deriv
!***********************************************************************
    subroutine deri_3d_inds(f,df,inds,j,lignored,lnometric)
!
!  dummy routine for compatibility
!
!  26-mar-12/MR: coded
!
      real, dimension (mx,my,mz)          :: f
      real, dimension (nx)                :: df
      integer                             :: j
      logical,                   optional :: lignored, lnometric
      integer, dimension(nx)              :: inds
!
      intent(in)  :: f,j,inds,lignored,lnometric
      intent(out) :: df
!
      call fatal_error('deri_3d_inds','Upwinding not implemented for nonuniform grids')
!
! dummy computation to avoid compiler warnings of unused variables
!
      if (present(lignored).and.present(lnometric)) &
          df  = inds + f(l1:l2,1,1) + j
!
    endsubroutine deri_3d_inds
!************************************************************************
    logical function heatflux_deriv_x(f, inh, fac, topbot)
!
!   dummy routine
!
!  17-apr-12/MR: coded
!
      use General, only: keep_compiler_quiet
!
      real, dimension(mx,my,mz,mfarray), intent(IN):: f
      real, dimension(my,mz)           , intent(IN):: inh
      real                             , intent(IN):: fac
      integer                          , intent(IN):: topbot
!
      heatflux_deriv_x = .true.
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(inh)
      call keep_compiler_quiet(fac)
      call keep_compiler_quiet(topbot)
!
    endfunction heatflux_deriv_x
!***********************************************************************
    subroutine bval_from_neumann_scl(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val is a constant.
!
!  30-sep-16/MR: coded
!
      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      real :: val

      integer :: k

!MR: coarse missing
      if (topbot=='bot') then
        if (idir==1) then
          k=l1
          f(l1,:,:,j) = (-val*60.*dx + 360.*f(k+1,:,:,j) &
                                     - 450.*f(k+2,:,:,j) &
                                     + 400.*f(k+3,:,:,j) &
                                     - 225.*f(k+4,:,:,j) &
                                     +  72.*f(k+5,:,:,j) &
                                     -  10.*f(k+6,:,:,j) )/147.
        elseif (idir==2) then
          k=m1
          f(:,m1,:,j) = (-val*60.*dy + 360.*f(:,k+1,:,j) &
                                     - 450.*f(:,k+2,:,j) &
                                     + 400.*f(:,k+3,:,j) &
                                     - 225.*f(:,k+4,:,j) &
                                     +  72.*f(:,k+5,:,j) &
                                     -  10.*f(:,k+6,:,j) )/147.
        elseif (idir==3) then
          k=n1
          f(:,:,n1,j) = (-val*60.*dz + 360.*f(:,:,k+1,j) &
                                     - 450.*f(:,:,k+2,j) &
                                     + 400.*f(:,:,k+3,j) &
                                     - 225.*f(:,:,k+4,j) &
                                     +  72.*f(:,:,k+5,j) &
                                     -  10.*f(:,:,k+6,j) )/147.
        endif
      else
        if (idir==1) then
          k=l2
          f(l2,:,:,j) = (val*60.*dx + 360.*f(k-1,:,:,j) &
                                    - 450.*f(k-2,:,:,j) &
                                    + 400.*f(k-3,:,:,j) &
                                    - 225.*f(k-4,:,:,j) &
                                    +  72.*f(k-5,:,:,j) &
                                    -  10.*f(k-6,:,:,j) )/147.
        elseif (idir==2) then
          k=m2
          f(:,m2,:,j) = (val*60.*dy + 360.*f(:,k-1,:,j) &
                                    - 450.*f(:,k-2,:,j) &
                                    + 400.*f(:,k-3,:,j) &
                                    - 225.*f(:,k-4,:,j) &
                                    +  72.*f(:,k-5,:,j) &
                                    -  10.*f(:,k-6,:,j) )/147.
        elseif (idir==3) then
          k=n2
          f(:,:,n2,j) = (val*60.*dz + 360.*f(:,:,k-1,j) &
                                    - 450.*f(:,:,k-2,j) &
                                    + 400.*f(:,:,k-3,j) &
                                    - 225.*f(:,:,k-4,j) &
                                    +  72.*f(:,:,k-5,j) &
                                    -  10.*f(:,:,k-6,j) )/147.
        endif
      endif

    endsubroutine bval_from_neumann_scl
!***********************************************************************
    subroutine bval_from_3rd_scl(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the 3rd kind BC d f/d x_i = val*f employing
!  one-sided difference formulae. val is a constant.
!
!  30-sep-16/MR: coded
!
      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      real :: val

      integer :: k

!MR: coarse missing
      if (topbot=='bot') then
        if (idir==1) then
          k=l1
          f(l1,:,:,j) = (  360.*f(k+1,:,:,j) &
                         - 450.*f(k+2,:,:,j) &
                         + 400.*f(k+3,:,:,j) &
                         - 225.*f(k+4,:,:,j) &
                         +  72.*f(k+5,:,:,j) &
                         -  10.*f(k+6,:,:,j) )/(147.+val*60.*dx)
        elseif (idir==2) then
          k=m1
          f(:,m1,:,j) = (  360.*f(:,k+1,:,j) &
                         - 450.*f(:,k+2,:,j) &
                         + 400.*f(:,k+3,:,j) &
                         - 225.*f(:,k+4,:,j) &
                         +  72.*f(:,k+5,:,j) &
                         -  10.*f(:,k+6,:,j) )/(147.+val*60.*dy)
        elseif (idir==3) then
          k=n1
          f(:,:,n1,j) = (  360.*f(:,:,k+1,j) &
                         - 450.*f(:,:,k+2,j) &
                         + 400.*f(:,:,k+3,j) &
                         - 225.*f(:,:,k+4,j) &
                         +  72.*f(:,:,k+5,j) &
                         -  10.*f(:,:,k+6,j) )/(147.+val*60.*dz)
        endif
      else
        if (idir==1) then
          k=l2
          f(l2,:,:,j) = (  360.*f(k-1,:,:,j) &
                         - 450.*f(k-2,:,:,j) &
                         + 400.*f(k-3,:,:,j) &
                         - 225.*f(k-4,:,:,j) &
                         +  72.*f(k-5,:,:,j) &
                         -  10.*f(k-6,:,:,j) )/(147.-val*60.*dx)
        elseif (idir==2) then
          k=m2
          f(:,m2,:,j) = (  360.*f(:,k-1,:,j) &
                         - 450.*f(:,k-2,:,j) &
                         + 400.*f(:,k-3,:,j) &
                         - 225.*f(:,k-4,:,j) &
                         +  72.*f(:,k-5,:,j) &
                         -  10.*f(:,k-6,:,j) )/(147.-val*60.*dy)
        elseif (idir==3) then
          k=n2
          f(:,:,n2,j) = (+ 360.*f(:,:,k-1,j) &
                         - 450.*f(:,:,k-2,j) &
                         + 400.*f(:,:,k-3,j) &
                         - 225.*f(:,:,k-4,j) &
                         +  72.*f(:,:,k-5,j) &
                         -  10.*f(:,:,k-6,j) )/(147.-val*60.*dz)
        endif
      endif

    endsubroutine bval_from_3rd_scl
!***********************************************************************
    subroutine bval_from_4th_scl(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the 4th kind BC d^2 f/d x_i^2 = val*f employing
!  one-sided difference formulae. val is a constant.
!
!  30-dec-16/MR: coded
!
      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      real :: val

      integer :: k

!MR: coarse missing
      if (topbot=='bot') then
        if (idir==1) then
          k=l1
          f(l1,:,:,j) = (- 3132.*f(k+1,:,:,j) &
                         + 5265.*f(k+2,:,:,j) &
                         - 5080.*f(k+3,:,:,j) &
                         + 2970.*f(k+4,:,:,j) &
                         -  972.*f(k+5,:,:,j) & 
                         +  137.*f(k+6,:,:,j)  )/(-812.+val*180.*dx*dx)
        elseif (idir==2) then
          k=m1
          f(:,m1,:,j) = (- 3132.*f(:,k+1,:,j) &
                         + 5265.*f(:,k+2,:,j) &
                         - 5080.*f(:,k+3,:,j) &
                         + 2970.*f(:,k+4,:,j) &
                         -  972.*f(:,k+5,:,j) & 
                         +  137.*f(:,k+6,:,j)  )/(-812.+val*180.*dy*dy)
        elseif (idir==3) then
          k=n1
          f(:,:,n1,j) = (- 3132.*f(:,:,k+1,j) &
                         + 5265.*f(:,:,k+2,j) &
                         - 5080.*f(:,:,k+3,j) &
                         + 2970.*f(:,:,k+4,j) &
                         -  972.*f(:,:,k+5,j) & 
                         +  137.*f(:,:,k+6,j)  )/(-812.+val*180.*dz*dz)
        endif
      else
        if (idir==1) then
          k=l2
          f(l2,:,:,j) = (- 3132.*f(k-1,:,:,j) &
                         + 5265.*f(k-2,:,:,j) &
                         - 5080.*f(k-3,:,:,j) &
                         + 2970.*f(k-4,:,:,j) &
                         -  972.*f(k-5,:,:,j) & 
                         +  137.*f(k-6,:,:,j)  )/(-812.+val*180.*dx*dx)
        elseif (idir==2) then
          k=m2
          f(:,m2,:,j) = (- 3132.*f(:,k-1,:,j) &
                         + 5265.*f(:,k-2,:,j) &
                         - 5080.*f(:,k-3,:,j) &
                         + 2970.*f(:,k-4,:,j) &
                         -  972.*f(:,k-5,:,j) & 
                         +  137.*f(:,k-6,:,j)  )/(-812.+val*180.*dy*dy)
        elseif (idir==3) then
          k=n2
          f(:,:,n2,j) = (- 3132.*f(:,:,k-1,j) &
                         + 5265.*f(:,:,k-2,j) &
                         - 5080.*f(:,:,k-3,j) &
                         + 2970.*f(:,:,k-4,j) &
                         -  972.*f(:,:,k-5,j) & 
                         +  137.*f(:,:,k-6,j)  )/(-812.+val*180.*dz*dz)
        endif
      endif

    endsubroutine bval_from_4th_scl
!***********************************************************************
    subroutine set_ghosts_for_onesided_ders(f,topbot,j,idir,l2nd)

!  20-sep-16/MR: added optional parameter bval for boundary value.
!                added ghost value setting for having the second derivative
!                correct in one-sided formulation for first inner point.
!   5-jan-17/MR: simplified, as ghost zone values for 1st and 2nd derivative
!                also in 2nd ghost zone point equal, as Fred spotted.
!
      use General, only: loptest

      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      logical, optional :: l2nd

      integer :: k,off

      if (loptest(l2nd)) then
        off=2
      else
        off=3
      endif

      if (topbot=='bot') then
        if (idir==1) then

          do k=l1-1,l1-off,-1
            f(k,:,:,j)=   7.*(f(k+1,:,:,j)-f(k+6,:,:,j)) &
                        -21.*(f(k+2,:,:,j)-f(k+5,:,:,j)) &
                        +35.*(f(k+3,:,:,j)-f(k+4,:,:,j)) &
                            + f(k+7,:,:,j)
          enddo
        elseif (idir==2) then

          do k=m1-1,m1-off,-1
            f(:,k,:,j)=   7.*(f(:,k+1,:,j)-f(:,k+6,:,j)) &
                        -21.*(f(:,k+2,:,j)-f(:,k+5,:,j)) &
                        +35.*(f(:,k+3,:,j)-f(:,k+4,:,j)) &
                            + f(:,k+7,:,j)
          enddo
        elseif (idir==3) then

          do k=n1-1,n1-off,-1
            f(:,:,k,j)=   7.*(f(:,:,k+1,j)-f(:,:,k+6,j)) &
                        -21.*(f(:,:,k+2,j)-f(:,:,k+5,j)) &
                        +35.*(f(:,:,k+3,j)-f(:,:,k+4,j)) &
                            + f(:,:,k+7,j)
          enddo
        endif
      else
        if (idir==1) then
          do k=l2+1,l2+off
            f(k,:,:,j)=   7.*(f(k-1,:,:,j)-f(k-6,:,:,j)) &
                        -21.*(f(k-2,:,:,j)-f(k-5,:,:,j)) &
                        +35.*(f(k-3,:,:,j)-f(k-4,:,:,j)) &
                            + f(k-7,:,:,j)
          enddo
        elseif (idir==2) then
          do k=m2+1,m2+off
            f(:,k,:,j)=   7.*(f(:,k-1,:,j)-f(:,k-6,:,j)) &
                        -21.*(f(:,k-2,:,j)-f(:,k-5,:,j)) &
                        +35.*(f(:,k-3,:,j)-f(:,k-4,:,j)) &
                            + f(:,k-7,:,j)
          enddo
        elseif (idir==3) then
          do k=n2+1,n2+off
            f(:,:,k,j)=   7.*(f(:,:,k-1,j)-f(:,:,k-6,j)) &
                        -21.*(f(:,:,k-2,j)-f(:,:,k-5,j)) &
                        +35.*(f(:,:,k-3,j)-f(:,:,k-4,j)) &
                            + f(:,:,k-7,j)
          enddo
        endif
      endif

    endsubroutine set_ghosts_for_onesided_ders
!***********************************************************************
    subroutine bval_from_neumann_arr(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val depends on x,y.
!
!  30-sep-16/MR: coded
!
      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val

      integer :: k

!MR: coarse missing
      if (topbot=='bot') then
        if (idir==1) then
          k=l1
          f(l1,:,:,j) = (-val*60.*dx + 360.*f(k+1,:,:,j) &
                                     - 450.*f(k+2,:,:,j) &
                                     + 400.*f(k+3,:,:,j) &
                                     - 225.*f(k+4,:,:,j) &
                                     +  72.*f(k+5,:,:,j) &
                                     -  10.*f(k+6,:,:,j) )/147.
        elseif (idir==2) then
          k=m1
          f(:,m1,:,j) = (-val*60.*dy + 360.*f(:,k+1,:,j) &
                                     - 450.*f(:,k+2,:,j) &
                                     + 400.*f(:,k+3,:,j) &
                                     - 225.*f(:,k+4,:,j) &
                                     +  72.*f(:,k+5,:,j) &
                                     -  10.*f(:,k+6,:,j) )/147.
        elseif (idir==3) then
          k=n1
          f(:,:,n1,j) = (-val*60.*dz + 360.*f(:,:,k+1,j) &
                                     - 450.*f(:,:,k+2,j) &
                                     + 400.*f(:,:,k+3,j) &
                                     - 225.*f(:,:,k+4,j) &
                                     +  72.*f(:,:,k+5,j) &
                                     -  10.*f(:,:,k+6,j) )/147.
        endif
      else
        if (idir==1) then
          k=l2
          f(l2,:,:,j) = (val*60.*dx + 360.*f(k-1,:,:,j) &
                                    - 450.*f(k-2,:,:,j) &
                                    + 400.*f(k-3,:,:,j) &
                                    - 225.*f(k-4,:,:,j) &
                                    +  72.*f(k-5,:,:,j) &
                                    -  10.*f(k-6,:,:,j) )/147.
        elseif (idir==2) then
          k=m2
          f(:,m2,:,j) = (val*60.*dy + 360.*f(:,k-1,:,j) &
                                    - 450.*f(:,k-2,:,j) &
                                    + 400.*f(:,k-3,:,j) &
                                    - 225.*f(:,k-4,:,j) &
                                    +  72.*f(:,k-5,:,j) &
                                    -  10.*f(:,k-6,:,j) )/147.
        elseif (idir==3) then
          k=n2
          f(:,:,n2,j) = (val*60.*dz + 360.*f(:,:,k-1,j) &
                                    - 450.*f(:,:,k-2,j) &
                                    + 400.*f(:,:,k-3,j) &
                                    - 225.*f(:,:,k-4,j) &
                                    +  72.*f(:,:,k-5,j) &
                                    -  10.*f(:,:,k-6,j) )/147.
        endif
      endif

    endsubroutine bval_from_neumann_arr
!***********************************************************************
    subroutine bval_from_3rd_arr(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val depends on x,y.
!
!  30-sep-16/MR: coded
!  09-feb-17/Ivan: completed dummy routine
!
      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val
!
!MR: coarse missing
      integer :: k

      if (topbot=='bot') then
        if (idir==1) then
          k=l1
          f(l1,:,:,j) = (  360.*f(k+1,:,:,j) &
                         - 450.*f(k+2,:,:,j) &
                         + 400.*f(k+3,:,:,j) &
                         - 225.*f(k+4,:,:,j) &
                         +  72.*f(k+5,:,:,j) &
                         -  10.*f(k+6,:,:,j) )/(147.+val*60.*dx)
        elseif (idir==2) then
          k=m1
          f(:,m1,:,j) = (  360.*f(:,k+1,:,j) &
                         - 450.*f(:,k+2,:,j) &
                         + 400.*f(:,k+3,:,j) &
                         - 225.*f(:,k+4,:,j) &
                         +  72.*f(:,k+5,:,j) &
                         -  10.*f(:,k+6,:,j) )/(147.+val*60.*dy)
        elseif (idir==3) then
          k=n1
          f(:,:,n1,j) = (  360.*f(:,:,k+1,j) &
                         - 450.*f(:,:,k+2,j) &
                         + 400.*f(:,:,k+3,j) &
                         - 225.*f(:,:,k+4,j) &
                         +  72.*f(:,:,k+5,j) &
                         -  10.*f(:,:,k+6,j) )/(147.+val*60.*dz)
        endif
      else
        if (idir==1) then
          k=l2
          f(l2,:,:,j) = (  360.*f(k-1,:,:,j) &
                         - 450.*f(k-2,:,:,j) &
                         + 400.*f(k-3,:,:,j) &
                         - 225.*f(k-4,:,:,j) &
                         +  72.*f(k-5,:,:,j) &
                         -  10.*f(k-6,:,:,j) )/(147.-val*60.*dx)
        elseif (idir==2) then
          k=m2
          f(:,m2,:,j) = (  360.*f(:,k-1,:,j) &
                         - 450.*f(:,k-2,:,j) &
                         + 400.*f(:,k-3,:,j) &
                         - 225.*f(:,k-4,:,j) &
                         +  72.*f(:,k-5,:,j) &
                         -  10.*f(:,k-6,:,j) )/(147.-val*60.*dy)
        elseif (idir==3) then
          k=n2
          f(:,:,n2,j) = (+ 360.*f(:,:,k-1,j) &
                         - 450.*f(:,:,k-2,j) &
                         + 400.*f(:,:,k-3,j) &
                         - 225.*f(:,:,k-4,j) &
                         +  72.*f(:,:,k-5,j) &
                         -  10.*f(:,:,k-6,j) )/(147.-val*60.*dz)
        endif
      endif

    endsubroutine bval_from_3rd_arr
!***********************************************************************
    subroutine bval_from_4th_arr(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the 4th kind BC d^2 f/d x_i^2 = val*f employing
!  one-sided difference formulae. val is a constant.
!
!  09-feb-17/Ivan: coded
!
!MR: coarse missing
      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val

      integer :: k

      if (topbot=='bot') then
        if (idir==1) then
          k=l1
          f(l1,:,:,j) = (- 3132.*f(k+1,:,:,j) &
                         + 5265.*f(k+2,:,:,j) &
                         - 5080.*f(k+3,:,:,j) &
                         + 2970.*f(k+4,:,:,j) &
                         -  972.*f(k+5,:,:,j) & 
                         +  137.*f(k+6,:,:,j)  )/(-812.+val*180.*dx*dx)
        elseif (idir==2) then
          k=m1
          f(:,m1,:,j) = (- 3132.*f(:,k+1,:,j) &
                         + 5265.*f(:,k+2,:,j) &
                         - 5080.*f(:,k+3,:,j) &
                         + 2970.*f(:,k+4,:,j) &
                         -  972.*f(:,k+5,:,j) & 
                         +  137.*f(:,k+6,:,j)  )/(-812.+val*180.*dy*dy)
        elseif (idir==3) then
          k=n1
          f(:,:,n1,j) = (- 3132.*f(:,:,k+1,j) &
                         + 5265.*f(:,:,k+2,j) &
                         - 5080.*f(:,:,k+3,j) &
                         + 2970.*f(:,:,k+4,j) &
                         -  972.*f(:,:,k+5,j) & 
                         +  137.*f(:,:,k+6,j)  )/(-812.+val*180.*dz*dz)
        endif
      else
        if (idir==1) then
          k=l2
          f(l2,:,:,j) = (- 3132.*f(k-1,:,:,j) &
                         + 5265.*f(k-2,:,:,j) &
                         - 5080.*f(k-3,:,:,j) &
                         + 2970.*f(k-4,:,:,j) &
                         -  972.*f(k-5,:,:,j) & 
                         +  137.*f(k-6,:,:,j)  )/(-812.+val*180.*dx*dx)
        elseif (idir==2) then
          k=m2
          f(:,m2,:,j) = (- 3132.*f(:,k-1,:,j) &
                         + 5265.*f(:,k-2,:,j) &
                         - 5080.*f(:,k-3,:,j) &
                         + 2970.*f(:,k-4,:,j) &
                         -  972.*f(:,k-5,:,j) & 
                         +  137.*f(:,k-6,:,j)  )/(-812.+val*180.*dy*dy)
        elseif (idir==3) then
          k=n2
          f(:,:,n2,j) = (- 3132.*f(:,:,k-1,j) &
                         + 5265.*f(:,:,k-2,j) &
                         - 5080.*f(:,:,k-3,j) &
                         + 2970.*f(:,:,k-4,j) &
                         -  972.*f(:,:,k-5,j) & 
                         +  137.*f(:,:,k-6,j)  )/(-812.+val*180.*dz*dz)
        endif
      endif

    endsubroutine bval_from_4th_arr
!***********************************************************************
 endmodule Deriv
