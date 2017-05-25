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
  public :: der6_other, der_pencil, der2_pencil
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
      real, parameter :: a = 1.0 / 60.0
      real, dimension(nx) :: fac
      logical :: withdx
!
!debug      if (loptimise_ders) der_call_count(k,icount_der,j,1) = & !DERCOUNT
!debug                            der_call_count(k,icount_der,j,1)+1 !DERCOUNT
!
      if (present(ignoredx)) then
        withdx = .not. ignoredx
        if (ignoredx) fac = a
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
          if (withdx) fac = a * dy_1(m)
          df=fac*(+ 45.0*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                  -  9.0*(f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) &
                  +      (f(l1:l2,m+3,n,k)-f(l1:l2,m-3,n,k)))
          if (withdx .and. lspherical_coords  ) df = df * r1_mn
          if (withdx .and. lcylindrical_coords) df = df * rcyl_mn1
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (withdx) fac = a * dz_1(n)
          df=fac*(+ 45.0*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k)) &
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
      integer :: l1_, l2_, sdf
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
            fac=(1./60)*dy_1(m)
            df=fac*(+ 45.0*(f(l1_:l2_,m+1,n)-f(l1_:l2_,m-1,n)) &
                    -  9.0*(f(l1_:l2_,m+2,n)-f(l1_:l2_,m-2,n)) &
                    +      (f(l1_:l2_,m+3,n)-f(l1_:l2_,m-3,n)))
            if (lspherical_coords)     df=df*r1_mn
            if (lcylindrical_coords)   df=df*rcyl_mn1
          else
            df=0.
            if (ip<=5) print*, 'der_other: Degenerate case in y-direction'
          endif
        elseif (j==3) then
          if (nzgrid/=1) then
            fac=(1./60)*dz_1(n)
            df=fac*(+ 45.0*(f(l1_:l2_,m,n+1)-f(l1_:l2_,m,n-1)) &
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
      real, dimension (:) :: pencil,df
      integer :: j
!
      intent(in)  :: j, pencil
      intent(out) :: df
!
!  x-derivative
!
      if (j==1) then
        if (size(pencil)/=mx) then
          if (lroot) print*, 'der_pencil: pencil must be of size mx for x derivative'
          call fatal_error('der_pencil','')
        endif
        df(l1:l2)=(1./60)*dx_1(l1:l2)*( &
            + 45.0*(pencil(l1+1:l2+1)-pencil(l1-1:l2-1)) &
            -  9.0*(pencil(l1+2:l2+2)-pencil(l1-2:l2-2)) &
            +      (pencil(l1+3:l2+3)-pencil(l1-3:l2-3)))
      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) then
          if (lroot) print*, 'der_pencil: pencil must be of size my for y derivative'
          call fatal_error('der_pencil','')
        endif
        df(m1:m2)=(1./60)*dy_1(m1:m2)*( &
            + 45.0*(pencil(m1+1:m2+1)-pencil(m1-1:m2-1)) &
            -  9.0*(pencil(m1+2:m2+2)-pencil(m1-2:m2-2)) &
            +      (pencil(m1+3:m2+3)-pencil(m1-3:m2-3)))
      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) then
          if (lroot) print*, 'der_pencil: pencil must be of size mz for z derivative'
          call fatal_error('der_pencil','')
        endif
        df(n1:n2)=(1./60)*dz_1(n1:n2)*( &
            + 45.0*(pencil(n1+1:n2+1)-pencil(n1-1:n2-1)) &
            -  9.0*(pencil(n1+2:n2+2)-pencil(n1-2:n2-2)) &
            +      (pencil(n1+3:n2+3)-pencil(n1-3:n2-3)))
      else
        if (lroot) print*, 'der_pencil: no such direction j=', j
        call fatal_error('der_pencil','')
      endif
!
      if (lcylindrical_coords.or.lspherical_coords) &
           call fatal_error("der_pencil","Not implemented for non-cartesian")
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
          fac=dy_1(m)**2
          df2=fac*(der2_coef0* f(l1:l2,m   ,n,k) &
                  +der2_coef1*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  +der2_coef2*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                  +der2_coef3*(f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
          if (.not.loptest(lwo_line_elem)) then
            if (lspherical_coords)   df2=df2*r2_mn
            if (lcylindrical_coords) df2=df2*rcyl_mn2
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
          fac=dz_1(n)**2
          df2=fac*( der2_coef0* f(l1:l2,m,n    ,k) &
                   +der2_coef1*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                   +der2_coef2*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                   +der2_coef3*(f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
          if (.not.loptest(lwo_line_elem)) then
            if (lspherical_coords) df2=df2*r2_mn*sin2th(m)
          endif
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dz_tilde(n)*df
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
          fac=(1./180)*dy_1(m)**2
          df2=fac*(-490.0*f(l1:l2,m,n) &
                   +270.0*(f(l1:l2,m+1,n)+f(l1:l2,m-1,n)) &
                   - 27.0*(f(l1:l2,m+2,n)+f(l1:l2,m-2,n)) &
                   +  2.0*(f(l1:l2,m+3,n)+f(l1:l2,m-3,n)))
          if (lspherical_coords)     df2=df2*r2_mn
          if (lcylindrical_coords)   df2=df2*rcyl_mn2
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dy_tilde(m)*df
          endif
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
          if (lspherical_coords) df2=df2*r2_mn*sin2th(m)
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dz_tilde(n)*df
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
!  x-derivative
!
      if (j==1) then
        if (size(pencil)/=mx) then
          if (lroot) print*, 'der2_pencil: pencil must be of size mx for x derivative'
          call fatal_error('der2_pencil','')
        endif
        df2=(1./180)*dx_1(l1:l2)**2*(-490.0*pencil(l1:l2) &
               +270.0*(pencil(l1+1:l2+1)+pencil(l1-1:l2-1)) &
               - 27.0*(pencil(l1+2:l2+2)+pencil(l1-2:l2-2)) &
               +  2.0*(pencil(l1+3:l2+3)+pencil(l1-3:l2-3)))
      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) then
          if (lroot) &
              print*, 'der2_pencil: pencil must be of size my for y derivative'
          call fatal_error('der2_pencil','')
        endif
        df2=(1./180)*dy_1(m1:m2)**2*(-490.0*pencil(m1:m2) &
               +270.0*(pencil(m1+1:m2+1)+pencil(m1-1:m2-1)) &
               - 27.0*(pencil(m1+2:m2+2)+pencil(m1-2:m2-2)) &
               +  2.0*(pencil(m1+3:m2+3)+pencil(m1-3:m2-3)))
      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) then
          if (lroot) &
              print*, 'der2_pencil: pencil must be of size mz for z derivative'
          call fatal_error('der2_pencil','')
        endif
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
          if (lcylindrical_coords)   df=df*rcyl_mn1**5
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: j,k
      logical, optional :: ignoredx,upwind
      logical :: igndx,upwnd
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
        if (.not. lequidist(j)) then
          call fatal_error('der6','for non-equidistant grid only '//&
              'if dx is ignored.')
          igndx = .true.
        endif
        igndx = .false.
      endif
!
      if (present(upwind)) then
        if (.not. lequidist(j)) then
          call fatal_error('der6','upwind cannot be used with '//&
              'non-equidistant grid.')
        endif
        upwnd = upwind
      else
        upwnd = .false.
!        if ((.not.lcartesian_coords).and.(.not.igndx)) then
!DM: non cartesian grids should not necessarily use upwinding. Wlad do you disagree ?
!         if (.not.igndx) then
!          call fatal_error('der6','in non-cartesian coordinates '//&
!              'just works if upwinding is used')
!        endif
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
            fac=1.
          else if (upwnd) then
            fac=(1.0/60)*dy_1(m)
          else
            fac=dy_1(m)**6
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
            fac=dz_1(n)**6
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
    subroutine der2_minmod(f,j,delfk,delfkp1,delfkm1,k)
!
! Calculates gradient of a scalar along the direction j but
! get for the derivatives at the point i-1,i,i+1
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
          if (lspherical_coords) fac = fac*r1_mn
          if (lcylindrical_coords) fac = fac*rcyl_mn1
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
          if (lspherical_coords) fac = fac*r1_mn*sin1th(m)
        case default
          call fatal_error('deriv:der2_minmod','wrong component')
        endselect
        delfkm1(:) = delf(:,-1)*fac
        delfk(:) = delf(:,0)*fac
        delfkp1(:) = delf(:,1)*fac
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
        upwnd = upwind
      else
        upwnd = .false.
        if (.not. lequidist(j)) &
             call fatal_error('der6_other','NOT IMPLEMENTED for '//&
             'non equidistant grid')
        if (.not.lcartesian_coords) &
             call fatal_error('der6_other','in non-cartesian coordinates '//&
             'just works if upwinding is used')
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=1.
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
            fac=1.
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
            fac=1.
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
!  Spherical polars. The comments about "minus extra terms" refer to the
!  presence of extra terms that are being evaluated later in gij_etc.
!
      if (loptest(lwo_line_elem)) return

      if (lspherical_coords) then
        if ((i==1.and.j==2)) df=df*r1_mn
        if ((i==2.and.j==1)) df=df*r1_mn           !(minus extra terms)
        if ((i==1.and.j==3)) df=df*r1_mn*sin1th(m)
        if ((i==3.and.j==1)) df=df*r1_mn*sin1th(m) !(minus extra terms)
        if ((i==2.and.j==3)) df=df*r2_mn*sin1th(m)
        if ((i==3.and.j==2)) df=df*r2_mn*sin1th(m) !(minus extra terms)
      endif
!
      if (lcylindrical_coords) then
        if ( i+j==3 .or. i+j==5 ) df=df*rcyl_mn1
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
!  Spherical polars. The comments about "minus extra terms" refer to the
!  presence of extra terms that are being evaluated later in gij_etc.
!
      if (lspherical_coords) then
        if ((i==1.and.j==2)) df=df*r1_mn
        if ((i==2.and.j==1)) df=df*r1_mn !(minus extra terms)
        if ((i==1.and.j==3)) df=df*r1_mn*sin1th(m)
        if ((i==3.and.j==1)) df=df*r1_mn*sin1th(m) !(minus extra terms)
        if ((i==2.and.j==3)) df=df*r2_mn*sin1th(m)
        if ((i==3.and.j==2)) df=df*r2_mn*sin1th(m) !(minus extra terms)
      endif
!
      if (lcylindrical_coords) then
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
!debug      if (loptimise_ders) der_call_count(k,icount_derij,i,j) = & !DERCOUNT
!debug                          der_call_count(k,icount_derij,i,j) + 1 !DERCOUNT
!
!
      !if (lspherical_coords.or.lcylindrical_coords) &
      !    call fatal_error('der5i1j','NOT IMPLEMENTED for non-cartesian coordinates')
      if (headtt) then
        if (lspherical_coords.or.lcylindrical_coords) &
             call warning('der5i1j','MISSING CURVATURE TERMS for non-cartesian coordinates')
      endif
!
      df=0.0
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
        if (nzgrid/=1.and.nygrid/=1) then
          fac=dz_1(n)**5*1/60.0*dy_1(m)
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
      df=0.0
      if ((i==1.and.j==1)) then
        if (nxgrid/=1) then
          call der6(f,k,df,j)
        else
          df=0.
          if (ip<=5) print*, 'der4i2j: Degenerate case in x-direction'
        endif
      elseif ((i==2.and.j==2)) then
        if (nygrid/=1) then
          call der6(f,k,df,j)
        else
          df=0.
          if (ip<=5) print*, 'der4i2j: Degenerate case in y-direction'
        endif
      elseif ((i==3.and.j==3)) then
        if (nzgrid/=1) then
          call der6(f,k,df,j)
        else
          df=0.
          if (ip<=5) print*, 'der4i2j: Degenerate case in z-direction'
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
        if (nzgrid/=1.and.nygrid/=1) then
          fac=(1./6.0*dz_1(n)**4) * (1./180.0*dx_1(l1:l2)**2)
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
    subroutine der3i3j(f,ik,df,i,j)
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: ik,i,j
      df=0.0
    endsubroutine der3i3j
!***********************************************************************          
    subroutine der3i2j1k(f,ik,df,i,j,k)
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: ik,i,j,k
      df=0.0
    endsubroutine der3i2j1k
!***********************************************************************
    subroutine der4i1j1k(f,ik,df,i,j,k)
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: ik,i,j,k
      df=0.0
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
  AUTOMATIC else* AUTOMATIC C  df(l) = (f(nghost+l,m+1,n,k) - ************  *****$Id$
!
!** AUTOMATIC CPndif* AUTOMATICenddo* AUTOMATPARAM.INC GENERdf=0.* AUTOMATICif (ip<=5) print*, 'der_upwind1st: Degenerate case in y-direction'
! variablhe number oPARAed bj == 3) then* AUTOMATed bnzgrid /= 1***************  do l=1,nx* AUTOMATIC Ced buu(l,3) > 0.dule Deriv
!
  usATION ***********************  ***
! Declare (f,n-1eneratz_1(n** AUTOMATIC CPARAM.INC GENERATION *******************,n+_derblic :: initializ
  puiv, finalize_deriv
  puhe number of f array
! variables and auxiliary variables added by this module
!
! CPARAM integer, parameter :: nghzst = 3
!
!***********************he nu!* AUTendsubroutine ! CPARAM inte
!* public :: der_onesided_4_slice
  public :: der_onesided_4_slice_other* AUT public :: der_onesided_4_slice_main(f,sgn,k,df,pos,j)
!
!   Calculmetex/y/z-derivative on a yz/xz/xy-c ::  at ****point pos. :: sUses a one-
  pu 4th order stencil_from_sgn = +1 forcoefward difference,f0, der-_coef1backer2_scoef2, der2.ic :: sBecause of its original intenm_4ter  in rela
!
! to solving :: scharacteristic equger, ss
  boundaries (NSCBC), thisblic should :: sreturn only PARTIAL derded_ders, NOT COVARIANT. Applying**** right_coef0caler :factors and conn 3
!
! termser :: i instead be done w****coef0eter ::: icunt_der2  = 2         !D :: icou7-jul-08/arne: coded ::   implreal, dimensug  (mx,my,mz,mfarray) :: number o!DERCOUNT
!debug :,:ter :d: icount_derer ::ac       RCOUgerer :  pukx_derj         RCOUNt(in) er ::,k
  pu_other 8     !DERCOoutDERCOUNT = 8     ****==odule Deriv
!
  *****x****/ction
    module  ger,=1./12.*dx_1(pos** AUTOMATICdf =ger,*(-sgn*25*f' va,m1:m2,n1:n2,k)& Cdata
!
  impleld
+dure48er_oth endi1er  ! derivative of another field
 edure36erface
!
  2er  ! derivative of another field
  endi1the der func3nterface der2                 ! Overload  erface
!
  4er  ! derivativ** AUTOMATles and auxiliary variables added by this module
!
! CPmod
  public :: ger, parameter :: nghxst = 3
  ! Overload the der!
!*******************************==2********************yure der_main   ! derivative of an 'm$
!
 variable
    module procedure der_l1:l2
  puerivative of another field
  endinterfded_4_sli
!
  ine der2                 ! Overload the nction
    mod   module procedure der2_main   ! derivanction
    mod3r' variable
    module procedure der2_onction
    mod4ive of another field
  endinterface
!
  interface derij                 ! Overload the der function
    module ost = 3
!
!*******************************==**************************der_main   ! derivative of an 'm findinterface
!
  interface  der_onesided_4_r  ! dnt_de          ! Overload the der functionbval_from
!
  inedure  der_onesided_4_slice_main  ! de_scl
    modul 'mvar' variable
    module procedure  d_scl
    modul3lice_main_pt
    module procedure  der__scl
    modul4 another field
  endinterface
!
  interface derij                 ! Overload the der function
    module ublic :: der_x,der2_x
  public :: der_zer2_z
  public :: der_mod
  public :: heatf
  public :: der_onesided_4_slice
  public :: der_onesided_4_slice_other
  public :: der2_minmod
  public :: heatf_ptlux_deriv_x
lll,mmm,nnnublic :: made us    !  23-sep-16/MR: added in. One rom_4tesided_dersis ct_ghostsicounatval_ ann, b(: removed olic :: 15-oct-09/Nataliaerij = 6         !DERCOUNT
!debug  integer, parameter :: icount_der
!
  UNT
!debug  integer, parameter :: icount_d: removed ofer_other = 8     !DERCOUNT
!
  inte: removed ofe der                 ! Overload the der function
    modulepos=lll  module procedure der_main   ! derivative of an 'mvar' variable
    module procedure der_other/180.
!
ve of another field
  endinterface
!
  int"
        call fatal_error('irload the der functio"
        call fatal_error('initiaerivative of an rmsg)
!
      endselect
!
    endsur2_other  ! der"
       other field
  endinterface
!
  interface derij                 ! Overload the der function
    module procedure derij_main   ! derivative of an 'mvar' variable
    module procedure .093mmmocedure derij_other  ! derivative of another field
  endinterface
!
  interface  der_onesidll_slice      call fatal_error('initialize_+wolf:   module)
!
      endselect
!
    endsubrouwolf: adapte2ze_deriv
!***************************wolf: adapte3*******************************
    swolf: adapte4_main(f, k, df, j, ignoredx)
!
!  calculate derivative df_k/dx_j
!  accurate to 6th order, explicit, periodic!
  interface bval_from_neumann
    module procedure bval_f.093nnre bval_from_neumann_scl
    module procedure bval_from_neumann_arr
  endinterface
!
  inter removerom_3rd
    module procedure bval_fro  logical module procedure bval_from_3rd_arr
  endise_ders) der_c 'mvar' variable
    module proceduise_ders) der_cscl
    module procedure bval_from_ise_ders) der_cace
!
  contains
!
!***********************************************************************
    subroutine initialize_deriv
!
!  Initialize stencil coefficients
!
!  23-sep-16/MR: added in_pt
  public :: der_onesided_4_slice
  public :: der_onesided_4_slice_other
  public :: der2_minmod
  public :: hotherlux_der_x
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
!debug  integerter :: icount_der_upwind1st = 7 !DERCOUNT
!debug  integer, parameter :: icount_d_other = 8     !DERCOUNT
!
  inrface der                 ! Overload the der function
    module procedure der_main   ! derivative of an 'mvar' variable
    module procedure der_other  ! derivave of another field
  endinterface
!
  interface der                ! Overload the der function
    modu procedure der2_main   ! derivative of an 'mvar' varble
    module procedure der2_other  ! derivative ofnother field
  endinterface
!
  interface derij                 ! Overload the der function
    module procedure derij_main   ! derivative of an 'mvar' variable
    module procedure derij_other  ! derivative of another field
  endinterface
!
  interface  der_one(sided_4_slice    )-function
    module pro)ve of another field
 
  endi23************
x(f,df)
!
!  routine der_x(f,dand. lsph derivative operating edure1x-dependent 1-D arrand. lsproutine der_x(f,dded_4_s x derivative operating on anx-dependent 1-D arrintent(iroutine der_x(f,d_4_slic))lic er ! derivative of another field
    module procedure  der_onesided_4_slice_other_pt
  endinterface
!
  interface bval_from_neumann
    module procedure bval_from_neumann_scl
    module procedure bval_from_neumann_arr
  endinterface
!
  interface bval_fromx derivative operating on anal_from_3rd_scl
    modulerom der_z; note that f is no
  endinterface
!
  interfx derivative operating on an procedure bval_from_4th_stion'
      endif
!
    endsufrom_4th_arr
  endinterfanx) :: fac
!
      if (nxgrid/=1) then
        fac=(1./60)*dx_1(l1:l2)
        df=fac*(+ 45.0*(f(l1+1:l2+ublic :: der_x,der2_x
  public :: der_z,der2_z
  public :: der_            -  9.0*(f(c = a * dx_1(l1:l2)
          df=fac*(+ 45.0*(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k)) &
                  -  9.0*(f(!   5-jan-/MR: removed offset mans_for_onesided_ders
!
      use General, only: indgen

      sel :: icose (der2_type)
!
      case)+f(l1-2:l2-axel:ERCOnged file nameparamhorter verdebuds  ) df = df * r1_mn
          if (withdx .and. lc0./180.
        der2_coef2=-27./180.; der2_coef3=2./180.
!_other = 8     !DERCOUNT
!
  in_coef0=-0.75; der2_coef1=0.34375
        der2_coef2=0.125; der2_coef3=-00.09375
!
      case default
        write(unit=errormsg,fmt=*) &
             "der2_type doesn't exist"
     if (ip<=5) print*, 'der_x: Degenderiv',errormsg)
tion'
      endif
!
    endsubroutine initialize_d**********************************************************
    subroutine der2_x(f,ubroutine der_maif, k, df, j, ignoredx)
!
!  calculate derivative df_k/dx_j
!  accurate to 6th orde-490.0*f( function
    module procedurivative of an 'mvar' variable
    module procedure dl: coded
!  18-jul-98/axel: corrected mx -> my and mx -> mz in all y and z ders
!   1-apr-01/ax****wolf: penc
!
  olf: adapted for x derivative operating on an x-depe  endsubroutine d!
    endsubroued 1/from der_z; note that f is not the f*****************!
    endsubrou 20-s x derivative operating on an df2)
!
!  z derivz-depe!
    endsubroux
!
 (nx) :: fac
!
      if (nxgrid/=1) then
        fac=(1./60)*dx_1(l1:l2)
        df=fa       -  9.0*(f(n1+2:n2+2)-f(!
  interface bval_from_neumann
    module procedure bval_frtional :: ignoredx
!
      real, parameter :: a = 1.0 / 60.0
      real, dimension(nx) :: fac
      logical  if (ip<=5) print*, 'der_x: Degeneise_ders) der_cation'
      endif
!
    endsubroutRCOUNT
!debug   **********************************count(k,icount_d********
    subroutine der2_x(f,desent(ignoredx))    real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(out) :: df2
!
      real, dimension  initialize_deriv
!
!  Initialize stencil coefficients
!
!  23-sep-16/MR: adde490.0*f(
  public :: der_onesided_4_slice
  public :: der_onesided_4_slice_other
  public :: der2f  !Dize_DERCOic :: Dummyz,der2_z
  public :: hifts by explici public :: der_onesided_4_slice
  public :: der_onesided_4_slice_other
  public :: der2_mii_3d_inds(f &
 ad t,j,lignored,lnometriclic :: dstru  faster!
or compatibilitucti!  26-mar-12/MRerij = ds  ) df = df * r1_mn
          if ( non-coorr :: icount_der_upwind1st = nxr non-coord-coord baUNT
!debuter :: i, dimension (mx,my,mz) :: f
:: r       logicERCOn (:) :: df
      op-2)) Degenace.
!  2 5-jun-04/ parameter :: COUNT
!debuow results dimensi:: ad tr = 8     !DERCOUNT
!
  injr inteace.
!  25-jun-04/                 ! Overload the dcall fpe)
_error(
! Cverload t','URAM i    not implement+3)+or nonuniform_neums'lic ::lf: adn-equtger, paraavoiOUNTmpiler warnings  = unused variablesd the der fpresERCOace.
!  ).and.id/=1) th-jun-04/t) e of anotherdf  =dimen + re  der_1,1) + er = 8  z
  public :: derverload t2/tony: coded - duplicate der_main but without k subscript
!             
  pubdimensi funebug  heatflux expli_x(f, inh,ger,, topbotlic :: slf: adapted fl2+2)+f7-ap1-feb-07/axel: added 1/er  Gparaml, = 3 : keep_j,1) + 1_quiet         !DERCOUNT
!debu  integer, paramet)))
DERCOIN)r :: icount_der_upwind1st (rs for non-coord        if (sdf>inher2_coef1=270 
        endif
         if       if (sdf>nr, parameter :: i       endif
         if (j==2) then
   else
 ) :: fac
n)-f(l1-3:l2-3,m = .trueariablesrs) d x-direction'
     (f** AUTOM(l1_:l2_,m-1,n)) &
      inh             -  9.0*(f(l1_:l2_,m+fac             -  9.0*(f(l1_:l2_,m+else
     n)-f(l13:l2+3,m,n)-f(l1-3:l2-3,m2/tony: coded - duplicate der_main but without k subscript
!                          thebval_from_neumann_scl(f,else
     dir,vallic :: et_ghostss:: icUNT
!dey value 
   :: icN       BC d f/d x_i =nera emploter ridsal_from_4toef2, der2    mulae.endifi, bvconstant :: ico30-sep-16b-07/axel: added 1/r and 1/pomega   integer,*ter :: icountunt_der2 (LENproc(m)
      parameter :: icou     i 'der2_x: Degenval
              -  9k1_:l2_,m,f (l1_:l2=='bot'*******************   inction
    module   k=l1  ! derivati(l1,:,:,j*****-val*60 'mv + 31_mnf(k+l_coords) derivative operating r_z; note that f i 45(m)
   2      else
            df=0.
            if (ip+ 40(m)
   3      else
            df=0.
            if (ip<=225m)
   4 'der_other: Degenerate case in z-direction'
     72m)
   5      else
            df=0.
            if (ip<= 1(m)
   6_coords))/147ariables a********l1_:l2dule procedure d    m  if (lspheric:,ml_cods) df=df*r1_mn*yin1th(m)
 :,       else
            df=0.
            if (ip<=5) prinmens, 'd_other: Degenerate case in z-direction'
          emensif
      endif
      endif
!
    endsubroutine der_othemens!******************************************************mens***********
    subroutine der_pencil(j,pencil,df)
!
! mensalcute first derivative of any x, y ****************    nnov-07/anders: a:,n1ed from der
!
   z  real, dimemensio(:) :: pencil,df
      integer :: j
!
      intent(inin)  :j, pencil
      intent(out) :: df
!
!  x-derivative
!
!
    if (j==1) then
        if (size(pencil)/=mx) then
          if (lroot) print*, 'der_pencil: pencil must be of s size  for x derivative'
          call fatal_error('der_pepencilte first derivative***************l2_,m,n+3)-f(l1_:l2_,m,n-3)))
          2 if (lspherica, 'der_ot= (df*r1_mn*sin1th(m)
  -       else
            df=0.
            if (i<=5) print-, 'der_other: Degenerate case in z-direction'
         en-if
        endif
      endif
!
    endsubroutin der_other-!********************************************************-*************
    subroutine der_pencil(j,penci,df)
!
!  -alculate first derivative of any x, y or z pencil.
!
!  01-=(1./60)*dy_1(:,m :: j, 
              real, dimencil(m1:m2+1)-pencil(m1-1:m2-1)) &
            -  9.0*(p1./6 :: j, pencil
      intent(out) :: df
!
!  x-erivative
!
encil1+3:m2+3)-pencil(m1-3:m2-3)))
      else if (j==3)+2)-      if (lroot) print*, 'der_pencil: pencil ust be of sin
        if (lroot) print*, 'der_pencil: pencil must bnt*,il','')
        endif
        df(l1:l2)=(1./60)*dx_1(l1:lr_pencil','')
  :,n-2:l2
           1:l2+1)-pencil(l--1:l2-1)) &
            -  9.0*(pencil(l1+2l2+2)-pencil(l--2:l2-2)) &
            +      (pencil(l1+3l2+3)-pencil(l--3:l2-3)))
      else if (j==2) then
!
!  yderivative
!
 -      if (size(pencil)/=my) then
          f (lroot) prin-*, 'der_pencil: pencil must be of size my fr y derivative-
          call fatal_error('der_he nucil coefficients
!
 else
            df=l_coords)     df=df*r1_mn
            if (lcylindrical_coords)   df=df*rcyl_mn1
          else
    3rd df=0.
            if (ip<=5) print*, 'der_other: Degenerate case in y3rd kindion'
          endi*ff
        elseif (j==3) then
          if (nzgrid/=1) then
            fac=(1./60)*dz_1(n)
            df=fac*(+ 45.0*(f(l1_:l2_,m,n+1)-f(l1_:l2_,m,n-1)) &
                    -  9.0*(f(l1_:l2_,m,n+2)-f(l1_:l2_,m,n-2)) &
                    +      (f(l1_:l2_,m,n+3)-f(l1_:l2_,m,n-3)))
            if (lspherical_coords) df 1th(m)
          else
            df=0.
      <=5) print*, 'der_other: Degenerate case in z-d        endif
        endif
      endif
!
    e der_other
!***********************************************************
    subroutine der_pen,df)
!
!  Calculate fi(rst +df*r1_mn*sother field
   any x, y or z pencil.
!
!  01-nov-07/anders: adapted fromnsion (nxmension (:) :: pencil,df
      integer  intent(in)  :: j, pencil
      intent(out) ::erivative
!
      if (j==1) then
        if (simx) then
          if (lroot) print*, 'der_pencust be of size mx for x derivative'
          cror('der_pencil','')
             derycall_count(k,icount_der2l2)=(1./60)*dx_1(l1:l2)*( &
            + 45.0*(   if (nxgrl(l1-1:l2-1)) &
            -  9.0*(l2+2)-pencil(l1-2:l2-2)) &
            +      (l2+3)-pencil(l1-3:l2-3)))
      else if (j==2) derivative
!
        if (size(pencil)/=my) thenf (lroot) print*, 'der_pencil: pencil must be or y derivative'
                  derzother field_error('der_pencil','')
        endif
        df(m1:m2)=(1./60)*dy_1(m1:m2)*( &
  nsion (nx)cil(m1+1:m2+1)-pencil(m1-1:m2-1)) &
 -  9.0*(pencil(m1+2:m2+2)-pencil(m1-2:m2-2)) &
 +      (pencil(m1+3:m2+3)-pencil(m1-3:m2-3)))
if (j==3) then
!
!  z-derivative
!
        if ()/=mz) then
          if (lroot) print*, 'der_pil must be of size mz      =df*r1_mn*scall_count(k,icount_der2,j,1) + 1 !DERCOUNT
!r_pencil','')
        endif   if (nxgri60)*dz_1(n1:n2)*( &
            + 4(n1+1:n2+1)-pencil(n1-1:n2-1)) &
            - l(n1+2:n2+2)-pencil(n1-2:n2-2)) &
            +il(n1+3:n2+3)-pencil(n1-3:n2-3)))
      else
  root) print*, 'der_pencil: no such direction j=  call fatal_error('demn2
          ecoef3*(f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))al_coords.or.lspherical_coo       call fatal_error("der_pencil","Not implenon-cartesian")
!
    endsubroutine der_pencil
**************************************************
    subroutine der2_main(f,k,df2,j,lwo_line calculate 2nd derivative d^2f_k/dx_j^2
!  accu order, explicit, perimn2
          e                 +der2_c explicit construction -> x6.5 fas
!
!   l_coords)     df=df*r1_mn
            if (lcylindrical_coords)   df=df*rcyl_mn1
          else
    4th df=0.
            if (ip<=5) print*, 'der_other: Degenerate case in yh
!
jun-04/to^2        ^2 adapted for non-equidistant grids
!  23-sep-16/MR: introduced offset variables dec(1./60)*dz_1(n)
            df=fac*(+ 45.0*(f(l1_:l2_,m,n+1)-f(l1_:l2_,m,n-1)) &
                    -  9.0*(f(l1_:l2_,m,n+2)-f(l1_:l2_,m,n-2)) &
                    +      (f(l1_:l2_,m,n+3)-f(l1_:l2_,m,n-3)))
            if (lspherical_coords) df= 313*******       else
            df=0.
      + 526_other
, 'der_other: Degenerate case in z-d- 508    endif
        endif
      endif
!
    e+ 297(m)
   !***********************************-  9****************** _coords) df2=df2*r2_mn*sin2th137!
!  Calculate     -812      12
!
dxer_call_count(k,icount_der2,j,1) + 1 !DERCOUNT
!
      if (j==1) then
     , dimensiomension (:) :: pencil,df
      integerj
!
      in)  :: j, pencil
      intent(out) ::: df2
!
!d
!
      if (j==1) then
        if (si_count(k,i          if (lroot) print*, 'der_penc           size mx for  der_call_count(k,icount_der2,j,1) +pencil','')T
!
!
      if (j==1y2_coef3*(f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)))
          if (.not.lequidi       df2=fl(l1-1:l2-1)) &
            -  9.0*(         +27l(l1-2:l2-2)) &
            +      (n)) &
      l(l1-3:l2-3)))
      else if (j==2) )+f(l1-2:l2-!
        if (size(pencil)/=my) then(f(l1+3:l2+3rint*, 'de-3:l2-3,m,n)))
          if (.not.lequiive'
    T
!
!
      if (j==1z
                  +der2_coef1*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  +der2_coe, dimensioncil(m1+1:m2+1)-pencil(m1-1:m2-1)) &
 j
!
      incil(m1+2:m2+2)-pencil(m1-2:m2-2)) &
: df2
!
!deencil(m1+3:m2+3)-pencil(m1-3:m2-3)))
_count(k,icthen
!
!  z-derivative
!
        if (           n
           der_call_count(k,icount_der2,j,1) +  of size mT
!
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=r_pencil','')
        endif       df2=fa60)*dz_1(n1:n2)*( &
            + 4         +270pencil(n1-1:n2-1)) &
            - n)) &
       -pencil(n1-2:n2-2)) &
            +)+f(l1-2:l2-2)-pencil(n1-3:n2-3)))
      else
  (f(l1+3:l2+3,, 'der_pe-3:l2-3,m,n)))
          if (.not.lequid_error('then
            call der(f,df,j)
            df2=df2+dx_tilde(l1:l2)al_coords.or.lspherical_coo else
         tal_error("der_pencil","Not implef (j==2) then
 )
!
    endsubroutine der_pencil
    fac=(1./180*********************************490.0*f(l1:l2,mtine der2_main(f,k,df2,j,lwo_linef(l1:l2,m+1,n)+ deriva,m-1,n)) &
                   - 27.0*(f(l1it, pe2,n)+f(l1:l2,m-2,n)) &
                   +  2. explicit construction -> x6.5 fas*******2/tony: coded - duplicate der_main but without k subscript
!                          theset_*****s_foent(out) ::derthe             ifl2nd)
rids
 fac=(1./60)*adm_4tnteger ::parame     elsdebug r: Degenerate_from_r x derivativt) pri*****nerate sett     or ha      !DEseconesided_ders
x for x derivativecor= 3
g  ial_from_4t    if bug  oef1,irsdx_1nery: ind_from_5-jan-17b-07/s!debifi  inas       zonlyeratescil(l1stERCOU2)
        endif
        df2=(1.alsog  i    ncil(l1-1:lann, b    l1)+peFred spott= 6          Degenerate case inloptest1_:l2_,m      df=fac*(+ 45.0*(f(l1_:l2_,m,n+1)-f(l1_:l2_,m,n-1)) &
                    -  9.0*(f(l1_:l2_dimensionnteger :: j
!2nd1_:l2_,m,n-2)) &
   ,ofplicit             e( then***************off==(1./60)PARAM.INC GENatal3  real, dimension            +      (f(l1_:l2_,m,n+3)-f(l1_:l2_,m,n-3)eriv
!
  use M    -1,l1-off,-  if (lsphererick_coords=7)
          else
            df=0.
  -21   intent(in)  :: f,j
      intent(ou+3der_ndif
        endif
      endif
!
 -       !*******************************+l(m1+2:*************
    subroutine der_-1:m2-1alculate se
            df=0.
     +    7_coordsumber of f array
! variables  any x, y or z pen &
             m +27m.0*(pencil(m1+1:m2+1)+p:,encim1-1:m2mension (:) :: pencil,df
      intil(m1+in)  :: j, pencil
      intent(out      
!
      if (j==1) then
        ifil(m1-          if (lroot) print*, 'der_n
!
!  size mx for x derivative'
        ncil)pencil','')          if (lroot) &
      mens    print*, 'der2_pencil: pencil must be of s******** &
             n +27n     call fatal_error('deder2encil',''l(l1-1:l2-1)) &
            -  9180)*dz_l(l1-2:l2-2)) &
            +   ) &
    l(l1-3:l2-3)))
      else if (j=2+1)+pen!
        if (size(pencil)/=my)  - 27.0*rint*, 'der_pencil: pencil must b-2)) &
ive'
          +  2.0*(pencil(n1+3:n2+3)3)+penprint*, 'der2_pencil: pencil_error('der_pencil','')
        endif
        df(m1:      2+1,l2+ sizl(m1+1:m2+1)+pencil(m1-1:m2-cil(m1+1:m2+1)-pencil(m1-1:m2-1))il(m1+2ncil(m1+2:m2+2)-pencil(m1-2:m2-2)       encil(m1+3:m2+3)-pencil(m1-3:m2-3il(m1-3then
!
!  z-derivative
!
        n
!
!  n
          if (lroot) print*, 'dencil)/ of size m          if (lroot) &
       -      print*, 'der2_pencil: pencil must be of size mz fr z derivative'
eal,mdimension (nx) :: dfder2_pencil','')60)*dz_1(n1:n2)*( &
           180)*dz_1pencil(n1-1:n2-1)) &
          ) &
     -pencil(n1-2:n2-2)) &
         2+1)+penc)-pencil(n1-3:n2-3)))
      els - 27.0*(, 'der_pencil: no such direction-2)) &
 _error('      +  2.0*(pencil(n1+3:n2+3)+     l(n1-3:n2-3)))
      else
        if (lroot) prin*, 'der2_pencil:eal,nnoredx)) then
            call fatal_etal_error("der_pencil","Not iendif
!
   )
!
    endsubroutine der_pen***************************************************tine der2_main(f,k,df2,j,lwo_utine der3( derivative d^2f_k/dx_j^2
!  aate 3rd deit, peve of a scalar, get scalar
!
!  10-feb-06/anders: adapted from der5
!
      r explicit construction - if (j==1) then
        if (l_coords)     df=df*r1_mn
            if (lcylindrical_coords)   df=df*rcyl_mn1
          else
            ar(l1+            if (ip<=5) print*, 'der_other: Degenerate case in y-direction'
          endif
        elseif (j==3) then
          if (nzgriddepz
  s
  x,y         fac=(1./60)*dz_1(n)
            df=fac*(+ 45.0*(f(l1_:l2_,m,n+1)-f(l1_:l2_,m,n-1)) &
                    -  9.0*(f(l1_:l2_,m,n  df=fac*(+ 7 !DERCOf(l1_:l2_,m,n-2)) &
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
      real, dimension (:) :: pencil,df
      integer :: j
!
      intent(in)  :: j, pencil
      intent(out) :: df
!
!  x-derivative
!
      if (j==1) then
        if (size(pencil)/=mx) then
          if (lroot) print*, 'der_pencil: pencil must be of size mx for x derivative'
          call fatal_error('der_pencil','')
        endif
        df(l1:l2)=(1./60)*dx_1(l1:l2)*( &
            + 45.0*(pencil(l1+1:l2+1)-pencil(l1-1:l2-1)) &
            -  9.0*(pencil(l1+2:l2+2)-pencil(l1-2:l2-2)) &
            +      (pencil(l1+3:l2+3)-pencil(l1-3:l2-3)))
      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) then
          if (lroot) print*, 'der_pencil: pencil must be of size my for y derivative'
          call fatal_error('der_pencil','')
        endif
        df(m1:m2)=(1./60)*dy_1(m1:m2)*( &
            + 45.0*(pencil(m1+1:m2+1)-pencil(m1-1:m2-1)) &
            -  9.0*(pencil(m1+2:m2+2)-pencil(m1-2:m2-2)) &
            +      (pencil(m1+3:m2+3)-pencil(m1-3:m2-3)))
      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) then
          if (lroot) print*, 'der_pencil: pencil must be of size mz for z derivative'
          call fatal_error('der_pencil','')
        endif
        df(n1:n2)=(1./60)*dz_1(n1:n2)*( &
            + 45.0*(pencil(n1+1:n2+1)-pencil(n1-1:n2-1)) &
            -  9.0*(pencil(n1+2:n2+2)-pencil(n1-2:n2-2)) &
            +      (pencil(n1+3:n2+3)-pencil(n1-3:n2-3)))
      else
        if (lroot) print*, 'der_pencil: no such direction j=', j
        call fatal_error('der_pencil','')
      endif
!
      if (lcylindrical_coords.or.lspherical_coords) &
           call fatal_error("der_pencil","Not implemented for non-cartesian")
!
    endsubroutine der_pencil
!***********************************************************************
    subroutine der2_main(f,k,df2,j,lwo_line_elem)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!  ar real, dimension (nx) :: fac
!
      if (nxgrid/=1) then
        fac=(1./180)*dx_1(l1:l2)*     endif
    else
            fac=(,3:l21.0/8)*1./dy**3
          endif
          df=fac*(- 13.0*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                  +  8.0*(f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) &
       09-feb0*(pIvanerijdebu     df=0.
         e if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) then
          if (lroot) &
       l_mn1**3
        else
      rror('dxterr :: j
3:l2r = 8     !DEm added
!
      use General, only: loptest

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df2,fac,df
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
          fac=dy_1(m)**2
          df2=fac*(der2_coef0* f(l1:l2,m   ,n,k) &
                  +der2_coef1*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  +der2_coef2*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                  +der2_coef3*(f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)))
          if (.not.loptest(lwo_line_elem)) then
            if (lspherical_coords)   df2=df2*r2_mn
            if (lcylindrical_coords) df2=df2*rcyl_mn2
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
          fac=dz_1(n)**2
          df2=fac*( der2_coef0* f(l1:l2,m,n    ,k) &
                   +der2_coef1*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                   +der2_coef2*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                   +der2_coef3*(f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)))
          if (.not.loptest(lwo_line_elem)) then
            if (lspherical_coords) df2=df2*r2_mn*sin2th(m)
          endif
          if (.not.lequidist(j)) then
            call der(f,k,df,j)
            df2=df2+dz_tilde(n)*df
          endif
           - 39.0*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  + 12.0 else
          df2=0.
        endif
      endif
!
    endsubroutine der2_main
!********************** else
            fac=(1.0/8)*1./dy**3
          endif
          df=f_other(f,df2,j)
!
!  calculate 2nd derivative d^2f/dx_j^2 (of scalar f)
!  accurate to 6th order, expln
        if (nzgr                  -  1.0*(f(l1:l2,m+3,n,k)-f(l1:l2,m-3,n,k)))
          if (lcylindrical_coords)   df=df*rcyl_mn1**3
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac dimension (nx) :: df2,fac,df
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
          fac=(1./180)*dy_1(m)**2
          df2=fac*(-490.0*f(l1:l2,m,n) &
                   +270.0*(f(l1:l2,m+1,n)+f(l1:l2,m-1,n)) &
                   - 27.0*(f(l1:l2,m+2,n)+f(l1:l2,m-2,n)) &
                   +  2.0*(f(l1:l2,m+3,n)+f(l1:l2,m-3,n)))
          if (lspherical_coords)     df2=df2*r2_mn
          if (lcylindrical_coords)   df2=df2*rcyl_mn2
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dy_tilde(m)*df
          endif
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
          if (lspherical_coords) df2=df2*r2_mn*sin2th(m)
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dz_tilde(n)*df
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
      inte        - 39.0*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  + 12.0*endmodule Dxplic