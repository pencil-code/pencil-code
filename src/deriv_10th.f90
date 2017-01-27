! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM integer, parameter :: nghost = 5
!
!***************************************************************
module Deriv
!
  use Cdata
  use Messages, only: fatal_error, warning
!
  implicit none
!
  private
!
  public :: initialize_deriv, finalize_deriv
  public :: der, der2, der3, der4, der5, der6, der10, derij, der5i1j
  public :: der6_other, der_pencil, der2_pencil
  public :: deri_3d_inds
  public :: der_upwind1st, der_z, der2_z
  public :: der_onesided_4_slice
  public :: der_onesided_4_slice_other
  public :: der2_minmod
  public :: heatflux_deriv_x
  public :: set_ghosts_for_onesided_ders
  public :: bval_from_neumann, bval_from_3rd, bval_from_4th_scl
!
  real :: der2_coef0, der2_coef1, der2_coef2, der2_coef3, der2_coef4, der2_coef5
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
    module procedure  der_onesided_4_slice_other ! derivative of another field
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
        der2_coef0=-73766./25200.; der2_coef1=42000./25200.
        der2_coef2=-6000./25200.; der2_coef3=1000./25200.
        der2_coef4=-125./25200.; der2_coef5=8./25200.;
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
    subroutine der_main(f,k,df,j,ignoredx)
!
!  calculate derivative df_k/dx_j
!  accurate to 10th order, explicit, periodic
!
!   1-oct-97/axel: coded
!  18-jul-98/axel: corrected mx -> my and mx -> mz in all y and z ders
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  21-feb-07/axel: added 1/r and 1/pomega factors for non-coord basis
!  25-aug-09/axel: adapted from deriv
!  12-dec-10/axel: adapted also y and z derivatives
!
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      logical, intent(in), optional :: ignoredx
      integer :: j,k
!
      intent(in)  :: f,k,j
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der,j,1) = & !DERCOUNT
!debug                            der_call_count(k,icount_der,j,1)+1 !DERCOUNT
!
      if (present(ignoredx)) call fatal_error('der_main', 'optional argument ignoredx is not implemented. ')
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=(1./2520.)*dx_1(l1:l2)
          df=fac*(+2100.*(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k)) &
                   -600.*(f(l1+2:l2+2,m,n,k)-f(l1-2:l2-2,m,n,k)) &
                   +150.*(f(l1+3:l2+3,m,n,k)-f(l1-3:l2-3,m,n,k)) &
                   - 25.*(f(l1+4:l2+4,m,n,k)-f(l1-4:l2-4,m,n,k)) &
                   +  2.*(f(l1+5:l2+5,m,n,k)-f(l1-5:l2-5,m,n,k)))
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=(1./2520.)*dy_1(m)
          df=fac*(+2100.*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                  - 600.*(f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) &
                  + 150.*(f(l1:l2,m+3,n,k)-f(l1:l2,m-3,n,k)) &
                  -  25.*(f(l1:l2,m+4,n,k)-f(l1:l2,m-4,n,k)) &
                  +   2.*(f(l1:l2,m+5,n,k)-f(l1:l2,m-5,n,k)))
          if (lspherical_coords)     df=df*r1_mn
          if (lcylindrical_coords)   df=df*rcyl_mn1
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=(1./2520.)*dz_1(n)
          df=fac*(+2100.*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k)) &
                  - 600.*(f(l1:l2,m,n+2,k)-f(l1:l2,m,n-2,k)) &
                  + 150.*(f(l1:l2,m,n+3,k)-f(l1:l2,m,n-3,k)) &
                  -  25.*(f(l1:l2,m,n+4,k)-f(l1:l2,m,n-4,k)) &
                  +   2.*(f(l1:l2,m,n+5,k)-f(l1:l2,m,n-5,k)))
          if (lspherical_coords) df=df*r1_mn*sin1th(m)
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
!  accurate to 10th order, explicit, periodic
!
!  26-nov-02/tony: coded, duplicate der_main but without k subscript, overload
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  21-feb-07/axel: added 1/r and 1/pomega factors for non-coord basis
!  25-aug-09/axel: adapted from deriv
!  12-dec-10/axel: adapted also y and z derivatives
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
          fac=(1./2520.)*dx_1(l1:l2)
          df=fac*(+2100.*(f(l1+1:l2+1,m,n)-f(l1-1:l2-1,m,n)) &
                  - 600.*(f(l1+2:l2+2,m,n)-f(l1-2:l2-2,m,n)) &
                  + 150.*(f(l1+3:l2+3,m,n)-f(l1-3:l2-3,m,n)) &
                  -  25.*(f(l1+4:l2+4,m,n)-f(l1-4:l2-4,m,n)) &
                  +   2.*(f(l1+5:l2+5,m,n)-f(l1-5:l2-5,m,n)))
        else
          df=0.
          if (ip<=5) print*, 'der_other: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          fac=(1./2520.)*dy_1(m)
          df=fac*(+2100.*(f(l1:l2,m+1,n)-f(l1:l2,m-1,n)) &
                  - 600.*(f(l1:l2,m+2,n)-f(l1:l2,m-2,n)) &
                  + 150.*(f(l1:l2,m+3,n)-f(l1:l2,m-3,n)) &
                  -  25.*(f(l1:l2,m+4,n)-f(l1:l2,m-4,n)) &
                  +   2.*(f(l1:l2,m+5,n)-f(l1:l2,m-5,n)))
          if (lspherical_coords)     df=df*r1_mn
          if (lcylindrical_coords)   df=df*rcyl_mn1
        else
          df=0.
          if (ip<=5) print*, 'der_other: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=(1./2520.)*dz_1(n)
          df=fac*(+2100.*(f(l1:l2,m,n+1)-f(l1:l2,m,n-1)) &
                  - 600.*(f(l1:l2,m,n+2)-f(l1:l2,m,n-2)) &
                  + 150.*(f(l1:l2,m,n+3)-f(l1:l2,m,n-3)) &
                  -  25.*(f(l1:l2,m,n+4)-f(l1:l2,m,n-4)) &
                  +   2.*(f(l1:l2,m,n+5)-f(l1:l2,m,n-5)))
          if (lspherical_coords) df=df*r1_mn*sin1th(m)
        else
          df=0.
          if (ip<=5) print*, 'der_other: Degenerate case in z-direction'
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
!  25-aug-09/axel: adapted from deriv
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
        df(l1:l2)=(1./2520)*dx_1(l1:l2)*( &
            + 2100.0*(pencil(l1+1:l2+1)-pencil(l1-1:l2-1)) &
            - 168.0*(pencil(l1+2:l2+2)-pencil(l1-2:l2-2)) &
            +      (pencil(l1+3:l2+3)-pencil(l1-3:l2-3)))
      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) then
          if (lroot) print*, 'der_pencil: pencil must be of size my for y derivative'
          call fatal_error('der_pencil','')
        endif
        df(m1:m2)=(1./2520)*dy_1(m1:m2)*( &
            + 2100.0*(pencil(m1+1:m2+1)-pencil(m1-1:m2-1)) &
            - 168.0*(pencil(m1+2:m2+2)-pencil(m1-2:m2-2)) &
            +      (pencil(m1+3:m2+3)-pencil(m1-3:m2-3)))
      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) then
          if (lroot) print*, 'der_pencil: pencil must be of size mz for z derivative'
          call fatal_error('der_pencil','')
        endif
        df(n1:n2)=(1./2520)*dz_1(n1:n2)*( &
            + 2100.0*(pencil(n1+1:n2+1)-pencil(n1-1:n2-1)) &
            - 168.0*(pencil(n1+2:n2+2)-pencil(n1-2:n2-2)) &
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
!  accurate to 10th order, explicit, periodic
!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  25-aug-09/axel: adapted from deriv
!  12-dec-10/axel: adapted also y and z derivatives
!  20-nov-16/MR: optional parameter lwo_line_elem added
!
      use General, only: loptest

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df2,fac,df
      integer :: j,k
      logical, optional :: lwo_line_elem
!
      intent(in)  :: f,k,j
      intent(out) :: df2
!
!debug      if (loptimise_ders) der_call_count(k,icount_der2,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der2,j,1) + 1 !DERCOUNT
!
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=dx_1(l1:l2)**2
          df2=fac*(der2_coef0*f (l1  :l2  ,m,n,k) &
                  +der2_coef1*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  +der2_coef2*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
                  +der2_coef3*(f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)) &
                  +der2_coef4*(f(l1+4:l2+4,m,n,k)+f(l1-4:l2-4,m,n,k)) &
                  +der2_coef5*(f(l1+5:l2+5,m,n,k)+f(l1-5:l2-5,m,n,k)))
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
          df2=fac*(der2_coef0*f(l1:l2,m  ,n,k) &
                  +der2_coef1*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  +der2_coef2*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
                  +der2_coef3*(f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)) &
                  +der2_coef4*(f(l1:l2,m+4,n,k)+f(l1:l2,m-4,n,k)) &
                  +der2_coef5*(f(l1:l2,m+5,n,k)+f(l1:l2,m-5,n,k)))
          if (.not.loptest(lwo_line_elem)) then
            if (lspherical_coords)     df2=df2*r2_mn
            if (lcylindrical_coords)   df2=df2*rcyl_mn2
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
          df2=fac*(der2_coef0*f(l1:l2,m,n  ,k) &
                   +der2_coef1*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                   +der2_coef2*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
                   +der2_coef3*(f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)) &
                   +der2_coef4*(f(l1:l2,m,n+4,k)+f(l1:l2,m,n-4,k)) &
                   +der2_coef5*(f(l1:l2,m,n+5,k)+f(l1:l2,m,n-5,k)))
          if (lspherical_coords.and..not.loptest(lwo_line_elem)) df2=df2*r2_mn*sin2th(m)
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
!  accurate to 10th order, explicit, periodic
!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  25-aug-09/axel: adapted from deriv
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
!
    endsubroutine der2_other
!***********************************************************************
    subroutine der2_pencil(j,pencil,df2)
!
!  Calculate 2nd derivative of any x, y or z pencil.
!
!  01-nov-07/anders: adapted from der2
!  25-aug-09/axel: adapted from deriv
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
          if (lroot) print*, 'der2_pencil: pencil must be of size my for y derivative'
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
          if (lroot) print*, 'der2_pencil: pencil must be of size mz for z derivative'
          call fatal_error('der2_pencil','')
        endif
        df2(n1:n2)=(1./180)*dz_1(n1:n2)**2*(-490.0*pencil(n1:n2) &
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
!  25-aug-09/axel: copied from deriv, but not adapted yet
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
!  25-aug-09/axel: copied from deriv, but not adapted yet
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
!
      if (lspherical_coords) &
           call fatal_error('der4','NOT IMPLEMENTED for spherical coordinates')
!
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
!  25-aug-09/axel: copied from deriv, but not adapted yet
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
          call fatal_error('der5','NOT IMPLEMENTED for no equidistant grid')
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
!  25-aug-09/axel: copied from deriv, but not adapted yet
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
        if (.not. lequidist(j)) then
          call fatal_error('der6','NOT IMPLEMENTED for non-equidistant grid')
        endif
        igndx = .false.
      endif
!
      if (present(upwind)) then
        upwnd = upwind
      else
        upwnd = .false.
        if ((.not.lcartesian_coords).and.(.not.igndx)) then
          call fatal_error('der6','in non-cartesian coordinates '//&
               'just works if upwinding is used')
        endif
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=1.0/60
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
            fac=1.0/60
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
            fac=1.0/60
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
!  25-aug-09/axel: copied from deriv, but not adapted yet
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
        if (.not. lequidist(j)) &
             call fatal_error('der6_other','NOT IMPLEMENTED for '//&
             'non equidistant grid')
        if (.not.lcartesian_coords) &
             call fatal_error('der6_other','in non-cartesian coordinates '//&
             'just works if upwiding is used')
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
    subroutine der10(f,k,df,j,ignoredx,upwind)
!
!  Calculate 10th derivative of a scalar, get scalar
!    Used for hyperdiffusion that affects small wave numbers as little as
!  possible (useful for density).
!    The optional flag IGNOREDX is useful for numerical purposes, where
!  you want to affect the Nyquist scale in each direction, independent of
!  the ratios dx:dy:dz.
!    The optional flag UPWIND is a variant thereof, which calculates
!  D^(10)*dx^9/840, which is the upwind correction of centered derivatives.
!
! 30-oct-13/pete: adapted from der6
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
        if (.not. lequidist(j)) then
          call fatal_error('der10','NOT IMPLEMENTED for non-equidistant grid')
        endif
        igndx = .false.
      endif
!
      if (present(upwind)) then
        upwnd = upwind
      else
        upwnd = .false.
        if ((.not.lcartesian_coords).and.(.not.igndx)) then
          call fatal_error('der10','in non-cartesian coordinates '//&
               'just works if upwinding is used')
        endif
      endif
!
      if (j==1) then
        if (nxgrid/=1) then
          if (igndx) then
            fac=1.0/840
          else if (upwnd) then
            fac=(1.0/840)*dx_1(l1:l2)
          else
            fac=dx_1(l1:l2)**10
          endif
          df=fac*(- 252.0* f(l1:l2,m,n,k) &
              +210.0*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
              -120.0*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) &
              + 45.0*(f(l1+3:l2+3,m,n,k)+f(l1-3:l2-3,m,n,k)) &
              - 10.0*(f(l1+4:l2+4,m,n,k)+f(l1-4:l2-4,m,n,k)) &
              +      (f(l1+5:l2+5,m,n,k)+f(l1-5:l2-5,m,n,k)))
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=1.0/840
          else if (upwnd) then
            fac=(1.0/840)*dy_1(m)
          else
            fac=dy_1(m)**10
          endif
          df=fac*(-252.0* f(l1:l2,m  ,n,k) &
              +210.0*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
              -120.0*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) &
              + 45.0*(f(l1:l2,m+3,n,k)+f(l1:l2,m-3,n,k)) &
              - 10.0*(f(l1:l2,m+4,n,k)+f(l1:l2,m-4,n,k)) &
              +      (f(l1:l2,m+5,n,k)+f(l1:l2,m-5,n,k)))
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=1.0/840
          else if (upwnd) then
            fac=(1.0/840)*dz_1(n)
          else
            fac=dz_1(n)**10
          endif
          df=fac*(-252.0* f(l1:l2,m,n  ,k) &
              +210.0*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
              -120.0*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) &
              + 45.0*(f(l1:l2,m,n+3,k)+f(l1:l2,m,n-3,k)) &
              - 10.0*(f(l1:l2,m,n+4,k)+f(l1:l2,m,n-4,k)) &
              +      (f(l1:l2,m,n+5,k)+f(l1:l2,m,n-5,k)))
        else
          df=0.
        endif
      endif
!
    endsubroutine der10
!***********************************************************************
    subroutine derij_main(f,k,df,i,j,lwo_line_elem)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!
!   8-sep-01/axel: coded
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  14-nov-06/wolf: implemented bidiagonal scheme
!  25-aug-09/axel: adapted from deriv
!  20-nov-16/MR: optional parameter lwo_line_elem added
!
      use General, only: loptest

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: i,j,k
      logical, optional :: lwo_line_elem
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
            fac=(1./100800.)*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
                       42000.*( f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k)  &
                               +f(l1-1:l2-1,m-1,n,k)-f(l1+1:l2+1,m-1,n,k)) &
                      - 6000.*( f(l1+2:l2+2,m+2,n,k)-f(l1-2:l2-2,m+2,n,k)  &
                               +f(l1-2:l2-2,m-2,n,k)-f(l1+2:l2+2,m-2,n,k)) &
                      + 1000.*( f(l1+3:l2+3,m+3,n,k)-f(l1-3:l2-3,m+3,n,k)  &
                               +f(l1-3:l2-3,m-3,n,k)-f(l1+3:l2+3,m-3,n,k)) &
                      -  125.*( f(l1+4:l2+4,m+4,n,k)-f(l1-4:l2-4,m+4,n,k)  &
                               +f(l1-4:l2-4,m-4,n,k)-f(l1+4:l2+4,m-4,n,k)) &
                      +    8.*( f(l1+5:l2+5,m+5,n,k)-f(l1-5:l2-5,m+5,n,k)  &
                               +f(l1-5:l2-5,m-5,n,k)-f(l1+5:l2+5,m-5,n,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=(1./100800.)*dy_1(m)*dz_1(n)
            df=fac*( &
                       42000.*( f(l1:l2,m+1,n+1,k)-f(l1:l2,m+1,n-1,k)  &
                               +f(l1:l2,m-1,n-1,k)-f(l1:l2,m-1,n+1,k)) &
                      - 6000.*( f(l1:l2,m+2,n+2,k)-f(l1:l2,m+2,n-2,k)  &
                               +f(l1:l2,m-2,n-2,k)-f(l1:l2,m-2,n+2,k)) &
                      + 1000.*( f(l1:l2,m+3,n+3,k)-f(l1:l2,m+3,n-3,k)  &
                               +f(l1:l2,m-3,n-3,k)-f(l1:l2,m-3,n+3,k)) &
                      -  125.*( f(l1:l2,m+4,n+4,k)-f(l1:l2,m+4,n-4,k)  &
                               +f(l1:l2,m-4,n-4,k)-f(l1:l2,m-4,n+4,k)) &
                      +    8.*( f(l1:l2,m+5,n+5,k)-f(l1:l2,m+5,n-5,k)  &
                               +f(l1:l2,m-5,n-5,k)-f(l1:l2,m-5,n+5,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./100800.)*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
                       42000.*( f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k)  &
                               +f(l1-1:l2-1,m,n-1,k)-f(l1+1:l2+1,m,n-1,k)) &
                      - 6000.*( f(l1+2:l2+2,m,n+2,k)-f(l1-2:l2-2,m,n+2,k)  &
                               +f(l1-2:l2-2,m,n-2,k)-f(l1+2:l2+2,m,n-2,k)) &
                      + 1000.*( f(l1+3:l2+3,m,n+3,k)-f(l1-3:l2-3,m,n+3,k)  &
                               +f(l1-3:l2-3,m,n-3,k)-f(l1+3:l2+3,m,n-3,k)) &
                      -  125.*( f(l1+4:l2+4,m,n+4,k)-f(l1-4:l2-4,m,n+4,k)  &
                               +f(l1-4:l2-4,m,n-4,k)-f(l1+4:l2+4,m,n-4,k)) &
                      +    8.*( f(l1+5:l2+5,m,n+5,k)-f(l1-5:l2-5,m,n+5,k)  &
                               +f(l1-5:l2-5,m,n-5,k)-f(l1+5:l2+5,m,n-5,k)) &
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
            fac=(1./2520.**2)*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
              2100.*((2100.*(f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k))  &
                    -9.*(f(l1+2:l2+2,m+1,n,k)-f(l1-2:l2-2,m+1,n,k))  &
                       +(f(l1+3:l2+3,m+1,n,k)-f(l1-3:l2-3,m+1,n,k))) &
                  -(2100.*(f(l1+1:l2+1,m-1,n,k)-f(l1-1:l2-1,m-1,n,k))  &
                    -9.*(f(l1+2:l2+2,m-1,n,k)-f(l1-2:l2-2,m-1,n,k))  &
                       +(f(l1+3:l2+3,m-1,n,k)-f(l1-3:l2-3,m-1,n,k))))&
              -9.*((2100.*(f(l1+1:l2+1,m+2,n,k)-f(l1-1:l2-1,m+2,n,k))  &
                    -9.*(f(l1+2:l2+2,m+2,n,k)-f(l1-2:l2-2,m+2,n,k))  &
                       +(f(l1+3:l2+3,m+2,n,k)-f(l1-3:l2-3,m+2,n,k))) &
                  -(2100.*(f(l1+1:l2+1,m-2,n,k)-f(l1-1:l2-1,m-2,n,k))  &
                    -9.*(f(l1+2:l2+2,m-2,n,k)-f(l1-2:l2-2,m-2,n,k))  &
                       +(f(l1+3:l2+3,m-2,n,k)-f(l1-3:l2-3,m-2,n,k))))&
                 +((2100.*(f(l1+1:l2+1,m+3,n,k)-f(l1-1:l2-1,m+3,n,k))  &
                    -9.*(f(l1+2:l2+2,m+3,n,k)-f(l1-2:l2-2,m+3,n,k))  &
                       +(f(l1+3:l2+3,m+3,n,k)-f(l1-3:l2-3,m+3,n,k))) &
                  -(2100.*(f(l1+1:l2+1,m-3,n,k)-f(l1-1:l2-1,m-3,n,k))  &
                    -9.*(f(l1+2:l2+2,m-3,n,k)-f(l1-2:l2-2,m-3,n,k))  &
                       +(f(l1+3:l2+3,m-3,n,k)-f(l1-3:l2-3,m-3,n,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=(1./2520.**2)*dy_1(m)*dz_1(n)
            df=fac*( &
              2100.*((2100.*(f(l1:l2,m+1,n+1,k)-f(l1:l2,m-1,n+1,k))  &
                    -9.*(f(l1:l2,m+2,n+1,k)-f(l1:l2,m-2,n+1,k))  &
                       +(f(l1:l2,m+3,n+1,k)-f(l1:l2,m-3,n+1,k))) &
                  -(2100.*(f(l1:l2,m+1,n-1,k)-f(l1:l2,m-1,n-1,k))  &
                    -9.*(f(l1:l2,m+2,n-1,k)-f(l1:l2,m-2,n-1,k))  &
                       +(f(l1:l2,m+3,n-1,k)-f(l1:l2,m-3,n-1,k))))&
              -9.*((2100.*(f(l1:l2,m+1,n+2,k)-f(l1:l2,m-1,n+2,k))  &
                    -9.*(f(l1:l2,m+2,n+2,k)-f(l1:l2,m-2,n+2,k))  &
                       +(f(l1:l2,m+3,n+2,k)-f(l1:l2,m-3,n+2,k))) &
                  -(2100.*(f(l1:l2,m+1,n-2,k)-f(l1:l2,m-1,n-2,k))  &
                    -9.*(f(l1:l2,m+2,n-2,k)-f(l1:l2,m-2,n-2,k))  &
                       +(f(l1:l2,m+3,n-2,k)-f(l1:l2,m-3,n-2,k))))&
                 +((2100.*(f(l1:l2,m+1,n+3,k)-f(l1:l2,m-1,n+3,k))  &
                    -9.*(f(l1:l2,m+2,n+3,k)-f(l1:l2,m-2,n+3,k))  &
                       +(f(l1:l2,m+3,n+3,k)-f(l1:l2,m-3,n+3,k))) &
                  -(2100.*(f(l1:l2,m+1,n-3,k)-f(l1:l2,m-1,n-3,k))  &
                    -9.*(f(l1:l2,m+2,n-3,k)-f(l1:l2,m-2,n-3,k))  &
                       +(f(l1:l2,m+3,n-3,k)-f(l1:l2,m-3,n-3,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./2520.**2)*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
              2100.*((2100.*(f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k))  &
                    -9.*(f(l1+2:l2+2,m,n+1,k)-f(l1-2:l2-2,m,n+1,k))  &
                       +(f(l1+3:l2+3,m,n+1,k)-f(l1-3:l2-3,m,n+1,k))) &
                  -(2100.*(f(l1+1:l2+1,m,n-1,k)-f(l1-1:l2-1,m,n-1,k))  &
                    -9.*(f(l1+2:l2+2,m,n-1,k)-f(l1-2:l2-2,m,n-1,k))  &
                       +(f(l1+3:l2+3,m,n-1,k)-f(l1-3:l2-3,m,n-1,k))))&
              -9.*((2100.*(f(l1+1:l2+1,m,n+2,k)-f(l1-1:l2-1,m,n+2,k))  &
                    -9.*(f(l1+2:l2+2,m,n+2,k)-f(l1-2:l2-2,m,n+2,k))  &
                       +(f(l1+3:l2+3,m,n+2,k)-f(l1-3:l2-3,m,n+2,k))) &
                  -(2100.*(f(l1+1:l2+1,m,n-2,k)-f(l1-1:l2-1,m,n-2,k))  &
                    -9.*(f(l1+2:l2+2,m,n-2,k)-f(l1-2:l2-2,m,n-2,k))  &
                       +(f(l1+3:l2+3,m,n-2,k)-f(l1-3:l2-3,m,n-2,k))))&
                 +((2100.*(f(l1+1:l2+1,m,n+3,k)-f(l1-1:l2-1,m,n+3,k))  &
                    -9.*(f(l1+2:l2+2,m,n+3,k)-f(l1-2:l2-2,m,n+3,k))  &
                       +(f(l1+3:l2+3,m,n+3,k)-f(l1-3:l2-3,m,n+3,k))) &
                  -(2100.*(f(l1+1:l2+1,m,n-3,k)-f(l1-1:l2-1,m,n-3,k))  &
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
        if ((i==2.and.j==1)) df=df*r1_mn !(minus extra terms)
        if ((i==1.and.j==3)) df=df*r1_mn*sin1th(m)
        if ((i==3.and.j==1)) df=df*r1_mn*sin1th(m) !(minus extra terms)
        if ((i==2.and.j==3)) df=df*r2_mn*sin1th(m)
        if ((i==3.and.j==2)) df=df*r2_mn*sin1th(m) !(minus extra terms)
      endif
!
      if (lcylindrical_coords) then
        if ((i+j==3.or.i+j==5)) df=df*rcyl_mn1
      endif
!
    endsubroutine derij_main
!***********************************************************************
    subroutine derij_other(f,df,i,j)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 10th order, explicit, periodic
!
!   8-sep-01/axel: coded
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  14-nov-06/wolf: implemented bidiagonal scheme
!  25-aug-09/axel: adapted from deriv
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
            fac=(1./2520.**2)*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
              2100.*((2100.*(f(l1+1:l2+1,m+1,n)-f(l1-1:l2-1,m+1,n))  &
                    -9.*(f(l1+2:l2+2,m+1,n)-f(l1-2:l2-2,m+1,n))  &
                       +(f(l1+3:l2+3,m+1,n)-f(l1-3:l2-3,m+1,n))) &
                  -(2100.*(f(l1+1:l2+1,m-1,n)-f(l1-1:l2-1,m-1,n))  &
                    -9.*(f(l1+2:l2+2,m-1,n)-f(l1-2:l2-2,m-1,n))  &
                       +(f(l1+3:l2+3,m-1,n)-f(l1-3:l2-3,m-1,n))))&
              -9.*((2100.*(f(l1+1:l2+1,m+2,n)-f(l1-1:l2-1,m+2,n))  &
                    -9.*(f(l1+2:l2+2,m+2,n)-f(l1-2:l2-2,m+2,n))  &
                       +(f(l1+3:l2+3,m+2,n)-f(l1-3:l2-3,m+2,n))) &
                  -(2100.*(f(l1+1:l2+1,m-2,n)-f(l1-1:l2-1,m-2,n))  &
                    -9.*(f(l1+2:l2+2,m-2,n)-f(l1-2:l2-2,m-2,n))  &
                       +(f(l1+3:l2+3,m-2,n)-f(l1-3:l2-3,m-2,n))))&
                 +((2100.*(f(l1+1:l2+1,m+3,n)-f(l1-1:l2-1,m+3,n))  &
                    -9.*(f(l1+2:l2+2,m+3,n)-f(l1-2:l2-2,m+3,n))  &
                       +(f(l1+3:l2+3,m+3,n)-f(l1-3:l2-3,m+3,n))) &
                  -(2100.*(f(l1+1:l2+1,m-3,n)-f(l1-1:l2-1,m-3,n))  &
                    -9.*(f(l1+2:l2+2,m-3,n)-f(l1-2:l2-2,m-3,n))  &
                       +(f(l1+3:l2+3,m-3,n)-f(l1-3:l2-3,m-3,n))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==3).or.(i==3.and.j==2)) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=(1./2520.**2)*dy_1(m)*dz_1(n)
            df=fac*( &
              2100.*((2100.*(f(l1:l2,m+1,n+1)-f(l1:l2,m-1,n+1))  &
                    -9.*(f(l1:l2,m+2,n+1)-f(l1:l2,m-2,n+1))  &
                       +(f(l1:l2,m+3,n+1)-f(l1:l2,m-3,n+1))) &
                  -(2100.*(f(l1:l2,m+1,n-1)-f(l1:l2,m-1,n-1))  &
                    -9.*(f(l1:l2,m+2,n-1)-f(l1:l2,m-2,n-1))  &
                       +(f(l1:l2,m+3,n-1)-f(l1:l2,m-3,n-1))))&
              -9.*((2100.*(f(l1:l2,m+1,n+2)-f(l1:l2,m-1,n+2))  &
                    -9.*(f(l1:l2,m+2,n+2)-f(l1:l2,m-2,n+2))  &
                       +(f(l1:l2,m+3,n+2)-f(l1:l2,m-3,n+2))) &
                  -(2100.*(f(l1:l2,m+1,n-2)-f(l1:l2,m-1,n-2))  &
                    -9.*(f(l1:l2,m+2,n-2)-f(l1:l2,m-2,n-2))  &
                       +(f(l1:l2,m+3,n-2)-f(l1:l2,m-3,n-2))))&
                 +((2100.*(f(l1:l2,m+1,n+3)-f(l1:l2,m-1,n+3))  &
                    -9.*(f(l1:l2,m+2,n+3)-f(l1:l2,m-2,n+3))  &
                       +(f(l1:l2,m+3,n+3)-f(l1:l2,m-3,n+3))) &
                  -(2100.*(f(l1:l2,m+1,n-3)-f(l1:l2,m-1,n-3))  &
                    -9.*(f(l1:l2,m+2,n-3)-f(l1:l2,m-2,n-3))  &
                       +(f(l1:l2,m+3,n-3)-f(l1:l2,m-3,n-3))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./2520.**2)*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
              2100.*((2100.*(f(l1+1:l2+1,m,n+1)-f(l1-1:l2-1,m,n+1))  &
                    -9.*(f(l1+2:l2+2,m,n+1)-f(l1-2:l2-2,m,n+1))  &
                       +(f(l1+3:l2+3,m,n+1)-f(l1-3:l2-3,m,n+1))) &
                  -(2100.*(f(l1+1:l2+1,m,n-1)-f(l1-1:l2-1,m,n-1))  &
                    -9.*(f(l1+2:l2+2,m,n-1)-f(l1-2:l2-2,m,n-1))  &
                       +(f(l1+3:l2+3,m,n-1)-f(l1-3:l2-3,m,n-1))))&
              -9.*((2100.*(f(l1+1:l2+1,m,n+2)-f(l1-1:l2-1,m,n+2))  &
                    -9.*(f(l1+2:l2+2,m,n+2)-f(l1-2:l2-2,m,n+2))  &
                       +(f(l1+3:l2+3,m,n+2)-f(l1-3:l2-3,m,n+2))) &
                  -(2100.*(f(l1+1:l2+1,m,n-2)-f(l1-1:l2-1,m,n-2))  &
                    -9.*(f(l1+2:l2+2,m,n-2)-f(l1-2:l2-2,m,n-2))  &
                       +(f(l1+3:l2+3,m,n-2)-f(l1-3:l2-3,m,n-2))))&
                 +((2100.*(f(l1+1:l2+1,m,n+3)-f(l1-1:l2-1,m,n+3))  &
                    -9.*(f(l1+2:l2+2,m,n+3)-f(l1-2:l2-2,m,n+3))  &
                       +(f(l1+3:l2+3,m,n+3)-f(l1-3:l2-3,m,n+3))) &
                  -(2100.*(f(l1+1:l2+1,m,n-3)-f(l1-1:l2-1,m,n-3))  &
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
        if ((i==1.and.j==2)) df=df*rcyl_mn1
        if ((i==2.and.j==1)) df=df*rcyl_mn1
        if ((i==1.and.j==3)) df=df
        if ((i==3.and.j==1)) df=df
        if ((i==2.and.j==3)) df=df*rcyl_mn1
        if ((i==3.and.j==2)) df=df*rcyl_mn1
      endif
!
    endsubroutine derij_other
!***********************************************************************
    subroutine der5i1j(f,k,df,i,j)
!
!  Calculate 6th derivative with respect to two different directions.
!
!  05-dec-06/anders: adapted from derij
!  25-aug-09/axel: copied from deriv, but not adapted yet
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df,fac
      integer :: i,j,k
!
!debug      if (loptimise_ders) der_call_count(k,icount_derij,i,j) = & !DERCOUNT
!debug                          der_call_count(k,icount_derij,i,j) + 1 !DERCOUNT
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
          fac=dx_1(l1:l2)**5*1/2520.0*dy_1(m)
          df=fac*( &
            2.5*((2100.*(f(l1+1:l2+1,m+1,n,k)-f(l1+1:l2+1,m-1,n,k))  &
                  -9.*(f(l1+1:l2+1,m+2,n,k)-f(l1+1:l2+1,m-2,n,k))  &
                     +(f(l1+1:l2+1,m+3,n,k)-f(l1+1:l2+1,m-3,n,k))) &
                -(2100.*(f(l1-1:l2-1,m+1,n,k)-f(l1-1:l2-1,m-1,n,k))  &
                  -9.*(f(l1-1:l2-1,m+2,n,k)-f(l1-1:l2-1,m-2,n,k))  &
                     +(f(l1-1:l2-1,m+3,n,k)-f(l1-1:l2-1,m-3,n,k))))&
           -2.0*((2100.*(f(l1+2:l2+2,m+1,n,k)-f(l1+2:l2+2,m-1,n,k))  &
                  -9.*(f(l1+2:l2+2,m+2,n,k)-f(l1+2:l2+2,m-2,n,k))  &
                     +(f(l1+2:l2+2,m+3,n,k)-f(l1+2:l2+2,m-3,n,k))) &
                -(2100.*(f(l1-2:l2-2,m+1,n,k)-f(l1-2:l2-2,m-1,n,k))  &
                  -9.*(f(l1-2:l2-2,m+2,n,k)-f(l1-2:l2-2,m-2,n,k))  &
                     +(f(l1-2:l2-2,m+3,n,k)-f(l1-2:l2-2,m-3,n,k))))&
           +0.5*((2100.*(f(l1+3:l2+3,m+1,n,k)-f(l1+3:l2+3,m-1,n,k))  &
                  -9.*(f(l1+3:l2+3,m+2,n,k)-f(l1+3:l2+3,m-2,n,k))  &
                     +(f(l1+3:l2+3,m+3,n,k)-f(l1+3:l2+3,m-3,n,k))) &
                -(2100.*(f(l1-3:l2-3,m+1,n,k)-f(l1-3:l2-3,m-1,n,k))  &
                  -9.*(f(l1-3:l2-3,m+2,n,k)-f(l1-3:l2-3,m-2,n,k))  &
                     +(f(l1-3:l2-3,m+3,n,k)-f(l1-3:l2-3,m-3,n,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in x- or y-direction'
        endif
      elseif ((i==2.and.j==1)) then
        if (nygrid/=1.and.nxgrid/=1) then
          fac=dy_1(m)**5*1/2520.0*dx_1(l1:l2)
          df=fac*( &
            2.5*((2100.*(f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k))  &
                  -9.*(f(l1+2:l2+2,m+1,n,k)-f(l1-2:l2-2,m+1,n,k))  &
                     +(f(l1+3:l2+3,m+1,n,k)-f(l1-3:l2-3,m+1,n,k))) &
                -(2100.*(f(l1+1:l2+1,m-1,n,k)-f(l1-1:l2-1,m-1,n,k))  &
                  -9.*(f(l1+2:l2+2,m-1,n,k)-f(l1-2:l2-2,m-1,n,k))  &
                     +(f(l1+3:l2+3,m-1,n,k)-f(l1-3:l2-3,m-1,n,k))))&
           -2.0*((2100.*(f(l1+1:l2+1,m+2,n,k)-f(l1-1:l2-1,m+2,n,k))  &
                  -9.*(f(l1+2:l2+2,m+2,n,k)-f(l1-2:l2-2,m+2,n,k))  &
                     +(f(l1+3:l2+3,m+2,n,k)-f(l1-3:l2-3,m+2,n,k))) &
                -(2100.*(f(l1+1:l2+1,m-2,n,k)-f(l1-1:l2-1,m-2,n,k))  &
                  -9.*(f(l1+2:l2+2,m-2,n,k)-f(l1-2:l2-2,m-2,n,k))  &
                     +(f(l1+3:l2+3,m-2,n,k)-f(l1-3:l2-3,m-2,n,k))))&
           +0.5*((2100.*(f(l1+1:l2+1,m+3,n,k)-f(l1-1:l2-1,m+3,n,k))  &
                  -9.*(f(l1+2:l2+2,m+3,n,k)-f(l1-2:l2-2,m+3,n,k))  &
                     +(f(l1+3:l2+3,m+3,n,k)-f(l1-3:l2-3,m+3,n,k))) &
                -(2100.*(f(l1+1:l2+1,m-3,n,k)-f(l1-1:l2-1,m-3,n,k))  &
                  -9.*(f(l1+2:l2+2,m-3,n,k)-f(l1-2:l2-2,m-3,n,k))  &
                     +(f(l1+3:l2+3,m-3,n,k)-f(l1-3:l2-3,m-3,n,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in y- or x-direction'
        endif
      elseif ((i==1.and.j==3)) then
        if (nxgrid/=1.and.nzgrid/=1) then
          fac=dx_1(l1:l2)**5*1/2520.0*dz_1(n)
          df=fac*( &
            2.5*((2100.*(f(l1+1:l2+1,m,n+1,k)-f(l1+1:l2+1,m,n-1,k))  &
                  -9.*(f(l1+1:l2+1,m,n+2,k)-f(l1+1:l2+1,m,n-2,k))  &
                     +(f(l1+1:l2+1,m,n+3,k)-f(l1+1:l2+1,m,n-3,k))) &
                -(2100.*(f(l1-1:l2-1,m,n+1,k)-f(l1-1:l2-1,m,n-1,k))  &
                  -9.*(f(l1-1:l2-1,m,n+2,k)-f(l1-1:l2-1,m,n-2,k))  &
                     +(f(l1-1:l2-1,m,n+3,k)-f(l1-1:l2-1,m,n-3,k))))&
           -2.0*((2100.*(f(l1+2:l2+2,m,n+1,k)-f(l1+2:l2+2,m,n-1,k))  &
                  -9.*(f(l1+2:l2+2,m,n+2,k)-f(l1+2:l2+2,m,n-2,k))  &
                     +(f(l1+2:l2+2,m,n+3,k)-f(l1+2:l2+2,m,n-3,k))) &
                -(2100.*(f(l1-2:l2-2,m,n+1,k)-f(l1-2:l2-2,m,n-1,k))  &
                  -9.*(f(l1-2:l2-2,m,n+2,k)-f(l1-2:l2-2,m,n-2,k))  &
                     +(f(l1-2:l2-2,m,n+3,k)-f(l1-2:l2-2,m,n-3,k))))&
           +0.5*((2100.*(f(l1+3:l2+3,m,n+1,k)-f(l1+3:l2+3,m,n-1,k))  &
                  -9.*(f(l1+3:l2+3,m,n+2,k)-f(l1+3:l2+3,m,n-2,k))  &
                     +(f(l1+3:l2+3,m,n+3,k)-f(l1+3:l2+3,m,n-3,k))) &
                -(2100.*(f(l1-3:l2-3,m,n+1,k)-f(l1-3:l2-3,m,n-1,k))  &
                  -9.*(f(l1-3:l2-3,m,n+2,k)-f(l1-3:l2-3,m,n-2,k))  &
                     +(f(l1-3:l2-3,m,n+3,k)-f(l1-3:l2-3,m,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in x- or z-direction'
        endif
      elseif ((i==3.and.j==1)) then
        if (nzgrid/=1.and.nygrid/=1) then
          fac=dz_1(n)**5*1/2520.0*dy_1(m)
          df=fac*( &
            2.5*((2100.*(f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k))  &
                  -9.*(f(l1+2:l2+2,m,n+1,k)-f(l1-2:l2-2,m,n+1,k))  &
                     +(f(l1+3:l2+3,m,n+1,k)-f(l1-3:l2-3,m,n+1,k))) &
                -(2100.*(f(l1+1:l2+1,m,n-1,k)-f(l1-1:l2-1,m,n-1,k))  &
                  -9.*(f(l1+2:l2+2,m,n-1,k)-f(l1-2:l2-2,m,n-1,k))  &
                     +(f(l1+3:l2+3,m,n-1,k)-f(l1-3:l2-3,m,n-1,k))))&
           -2.0*((2100.*(f(l1+1:l2+1,m,n+2,k)-f(l1-1:l2-1,m,n+2,k))  &
                  -9.*(f(l1+2:l2+2,m,n+2,k)-f(l1-2:l2-2,m,n+2,k))  &
                     +(f(l1+3:l2+3,m,n+2,k)-f(l1-3:l2-3,m,n+2,k))) &
                -(2100.*(f(l1+1:l2+1,m,n-2,k)-f(l1-1:l2-1,m,n-2,k))  &
                  -9.*(f(l1+2:l2+2,m,n-2,k)-f(l1-2:l2-2,m,n-2,k))  &
                     +(f(l1+3:l2+3,m,n-2,k)-f(l1-3:l2-3,m,n-2,k))))&
           +0.5*((2100.*(f(l1+1:l2+1,m,n+3,k)-f(l1-1:l2-1,m,n+3,k))  &
                  -9.*(f(l1+2:l2+2,m,n+3,k)-f(l1-2:l2-2,m,n+3,k))  &
                     +(f(l1+3:l2+3,m,n+3,k)-f(l1-3:l2-3,m,n+3,k))) &
                -(2100.*(f(l1+1:l2+1,m,n-3,k)-f(l1-1:l2-1,m,n-3,k))  &
                  -9.*(f(l1+2:l2+2,m,n-3,k)-f(l1-2:l2-2,m,n-3,k))  &
                     +(f(l1+3:l2+3,m,n-3,k)-f(l1-3:l2-3,m,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in z- or x-direction'
        endif
      elseif ((i==2.and.j==3)) then
        if (nygrid/=1.and.nzgrid/=1) then
          fac=dy_1(m)**5*1/2520.0*dz_1(n)
          df=fac*( &
            2.5*((2100.*(f(l1:l2,m+1,n+1,k)-f(l1:l2,m+1,n-1,k))  &
                  -9.*(f(l1:l2,m+1,n+2,k)-f(l1:l2,m+1,n-2,k))  &
                     +(f(l1:l2,m+1,n+3,k)-f(l1:l2,m+1,n-3,k))) &
                -(2100.*(f(l1:l2,m-1,n+1,k)-f(l1:l2,m-1,n-1,k))  &
                  -9.*(f(l1:l2,m-1,n+2,k)-f(l1:l2,m-1,n-2,k))  &
                     +(f(l1:l2,m-1,n+3,k)-f(l1:l2,m-1,n-3,k))))&
           -2.0*((2100.*(f(l1:l2,m+2,n+1,k)-f(l1:l2,m+2,n-1,k))  &
                  -9.*(f(l1:l2,m+2,n+2,k)-f(l1:l2,m+2,n-2,k))  &
                     +(f(l1:l2,m+2,n+3,k)-f(l1:l2,m+2,n-3,k))) &
                -(2100.*(f(l1:l2,m-2,n+1,k)-f(l1:l2,m-2,n-1,k))  &
                  -9.*(f(l1:l2,m-2,n+2,k)-f(l1:l2,m-2,n-2,k))  &
                     +(f(l1:l2,m-2,n+3,k)-f(l1:l2,m-2,n-3,k))))&
           +0.5*((2100.*(f(l1:l2,m+3,n+1,k)-f(l1:l2,m+3,n-1,k))  &
                  -9.*(f(l1:l2,m+3,n+2,k)-f(l1:l2,m+3,n-2,k))  &
                     +(f(l1:l2,m+3,n+3,k)-f(l1:l2,m+3,n-3,k))) &
                -(2100.*(f(l1:l2,m-3,n+1,k)-f(l1:l2,m-3,n-1,k))  &
                  -9.*(f(l1:l2,m-3,n+2,k)-f(l1:l2,m-3,n-2,k))  &
                     +(f(l1:l2,m-3,n+3,k)-f(l1:l2,m-3,n-3,k))))&
                 )
        else
          df=0.
          if (ip<=5) print*, 'der5i1j: Degenerate case in y- or z-direction'
        endif
      elseif ((i==3.and.j==2)) then
        if (nzgrid/=1.and.nygrid/=1) then
          fac=dz_1(n)**5*1/2520.0*dy_1(m)
          df=fac*( &
            2.5*((2100.*(f(l1:l2,m+1,n+1,k)-f(l1:l2,m-1,n+1,k))  &
                  -9.*(f(l1:l2,m+2,n+1,k)-f(l1:l2,m-2,n+1,k))  &
                     +(f(l1:l2,m+3,n+1,k)-f(l1:l2,m-3,n+1,k))) &
                -(2100.*(f(l1:l2,m+1,n-1,k)-f(l1:l2,m-1,n-1,k))  &
                  -9.*(f(l1:l2,m+2,n-1,k)-f(l1:l2,m-2,n-1,k))  &
                     +(f(l1:l2,m+3,n-1,k)-f(l1:l2,m-3,n-1,k))))&
           -2.0*((2100.*(f(l1:l2,m+1,n+2,k)-f(l1:l2,m-1,n+2,k))  &
                  -9.*(f(l1:l2,m+2,n+2,k)-f(l1:l2,m-2,n+2,k))  &
                     +(f(l1:l2,m+3,n+2,k)-f(l1:l2,m-3,n+2,k))) &
                -(2100.*(f(l1:l2,m+1,n-2,k)-f(l1:l2,m-1,n-2,k))  &
                  -9.*(f(l1:l2,m+2,n-2,k)-f(l1:l2,m-2,n-2,k))  &
                     +(f(l1:l2,m+3,n-2,k)-f(l1:l2,m-3,n-2,k))))&
           +0.5*((2100.*(f(l1:l2,m+1,n+3,k)-f(l1:l2,m-1,n+3,k))  &
                  -9.*(f(l1:l2,m+2,n+3,k)-f(l1:l2,m-2,n+3,k))  &
                     +(f(l1:l2,m+3,n+3,k)-f(l1:l2,m-3,n+3,k))) &
                -(2100.*(f(l1:l2,m+1,n-3,k)-f(l1:l2,m-1,n-3,k))  &
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
      if (lspherical_coords.or.lcylindrical_coords) &
           call fatal_error('der5i1j','NOT IMPLEMENTED for non-cartesian coordinates')
!
    endsubroutine der5i1j
!***********************************************************************
    subroutine der_upwind1st(f,uu,k,df,j)
!
!  First order upwind derivative of variable
!  Useful for advecting non-logarithmic variables
!
!  25-aug-09/axel: copied from deriv, but not adapted yet
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
      if (.not. lequidist(j)) &
        call fatal_error('der_upwind1st','NOT IMPLEMENTED for no equidistant grid')
!
      if (lspherical_coords.or.lcylindrical_coords) &
           call fatal_error('der_upwind1st','NOT IMPLEMENTED for non-cartesian grid')
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
!  25-aug-09/axel: copied from deriv, but not adapted yet
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
                  -sgn*36*f(l1:l2,m1:m2,pos+sgn*1,k)&
                  +sgn*16*f(l1:l2,m1:m2,pos+sgn*1,k)&
                  -sgn*3 *f(l1:l2,m1:m2,pos+sgn*1,k))
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in z-direction'
        endif
      endif
    endsubroutine
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
!  25-aug-09/axel: copied from deriv, but not adapted yet
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
          df = fac*(-sgn*25*f(l1:l2,pos,n1:n2)&
                  +sgn*48*f(l1:l2,pos+sgn*1,n1:n2)&
                  -sgn*36*f(l1:l2,pos+sgn*2,n1:n2)&
                  +sgn*16*f(l1:l2,pos+sgn*3,n1:n2)&
                  -sgn*3 *f(l1:l2,pos+sgn*4,n1:n2))
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          fac=1./12.*dz_1(pos)
          df = fac*(-sgn*25*f(l1:l2,m1:m2,pos)&
                  +sgn*48*f(l1:l2,m1:m2,pos+sgn*1)&
                  -sgn*36*f(l1:l2,m1:m2,pos+sgn*1)&
                  +sgn*16*f(l1:l2,m1:m2,pos+sgn*1)&
                  -sgn*3 *f(l1:l2,m1:m2,pos+sgn*1))
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in z-direction'
        endif
      endif
    endsubroutine
!***********************************************************************
    subroutine der_z(f,df)
!
! Dummy routine.
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(out) :: df
!
      call fatal_error("deriv_10t","der_z not implemented yet")
!
! To avoid compiler warnings:
      df=f(n1:n2)
!
    endsubroutine der_z
!***********************************************************************
    subroutine der2_z(f,df2)
!
! Dummy routine.
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(out) :: df2
!
      call fatal_error("deriv_10th","der2_z not implemented yet")
!
! To avoid compiler warnings:
      df2=f(n1:n2)
!
    endsubroutine der2_z
!***********************************************************************
    subroutine der2_minmod(f,j,delfk,delfkp1,delfkm1,k)
!
!  Dummy routine
!
      intent(in) :: f,k,j
      intent(out) :: delfk,delfkp1,delfkm1
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: delfk,delfkp1,delfkm1
      integer :: j,k
!
      call fatal_error('der2_minmod','Not implemented for deriv_10th')
!
!  Fill with dummy values to keep compiler quiet
      delfk(:) = j; delfkp1(:) = k; delfkm1(:) = f(l1,m1,n1,1)
!
    endsubroutine der2_minmod
!***********************************************************************
    subroutine finalize_deriv()
!
! Dummy
!
    endsubroutine finalize_deriv
!***********************************************************************
    subroutine deri_3d_inds(f,df,inds,j,lignored,lnometric)
!
!  dummy routine for compatibility
!
!  26-mar-12/MR: coded
!
!      use General, only: keep_compiler_quiet

      real, dimension (mx,my,mz)          :: f
      real, dimension (nx)                :: df
      integer                             :: j
      logical,                   optional :: lignored, lnometric
      integer, dimension(nx)              :: inds
!
      intent(in)  :: f,j,inds,lignored,lnometric
      intent(out) :: df
!
!      call keep_compiler_quiet(df)
      call fatal_error('deri_3d_inds','Upwinding not implemented for nonuniform grids')
!
! dummy computation to avoid compiler warnings of unused variables
      if (present(lignored).and.present(lnometric)) &
          df  = inds + f(l1:l2,1,1) + j
!
    endsubroutine deri_3d_inds
!************************************************************************
    logical function heatflux_deriv_x( f, inh, fac, topbot )
!
!   dummy routine
!
!  17-apr-12/MR: coded
!
     real, dimension(mx,my,mz,mfarray), intent(IN):: f
     real, dimension(my,mz)           , intent(IN):: inh
     real                             , intent(IN):: fac
     integer                          , intent(IN):: topbot
!
     heatflux_deriv_x = .false.

    endfunction heatflux_deriv_x
!***********************************************************************
    subroutine set_ghosts_for_onesided_ders(f,topbot,j,idir,l2nd_)
!
!  Calculates the ghost point value. The coefficients are derived from two FD formulae:
!  1) derivative is evaluated at point 4 for the given grid -1 0 1 2 3 |4| 5 6 7 8 9
!  2) derivative is evaluated at point 4 for the other grid    0 1 2 3 |4| 5 6 7 8 9 10
!  the second expression is substituted into the first equation and then solved for f(i-1)
!  resulting in onesided formula for the ghost point.
!
!  24-jan-17/Ivan: coded.
!
      use General, only: loptest

      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      logical, optional :: l2nd_

      integer :: k,off

      !call fatal_error('set_ghosts_for_onesided_ders','Not implemented for 10th order.')
      if (loptest(l2nd)) then
        off=4
      else
        off=5
      endif

      if (topbot=='bot') then
        if (idir==1) then

          do k=l1-1,l1-off,-1
            f(k,:,:,j)=11*f(k+1,:,:,j) &
                      -55*f(k+2,:,:,j) &
                     +165*f(k+3,:,:,j) &
                     -330*f(k+4,:,:,j) &
                     +462*f(k+5,:,:,j) &
                     -462*f(k+6,:,:,j) &
                     +330*f(k+7,:,:,j) &
                     -165*f(k+8,:,:,j) &
                      +55*f(k+9,:,:,j) &
                      -11*f(k+10,:,:,j) &
                         +f(k+11,:,:,j)
          enddo
        elseif (idir==2) then

          do k=m1-1,m1-off,-1
            f(:,k,:,j)=11*f(:,k+1,:,j) &
                      -55*f(:,k+2,:,j) &
                     +165*f(:,k+3,:,j) &
                     -330*f(:,k+4,:,j) &
                     +462*f(:,k+5,:,j) &
                     -462*f(:,k+6,:,j) &
                     +330*f(:,k+7,:,j) &
                     -165*f(:,k+8,:,j) &
                      +55*f(:,k+9,:,j) &
                      -11*f(:,k+10,:,j) &
                         +f(:,k+11,:,j)
          enddo
        elseif (idir==3) then

          do k=n1-1,n1-off,-1
            f(:,:,k,j)=11*f(:,:,k+1,j) &
                      -55*f(:,:,k+2,j) &
                     +165*f(:,:,k+3,j) &
                     -330*f(:,:,k+4,j) &
                     +462*f(:,:,k+5,j) &
                     -462*f(:,:,k+6,j) &
                     +330*f(:,:,k+7,j) &
                     -165*f(:,:,k+8,j) &
                      +55*f(:,:,k+9,j) &
                      -11*f(:,:,k+10,j) &
                         +f(:,:,k+11,j)
          enddo
        endif
      else
        if (idir==1) then
          do k=l2+1,l2+off
            f(k,:,:,j)=11*f(k-1,:,:,j) &
                      -55*f(k-2,:,:,j) &
                     +165*f(k-3,:,:,j) &
                     -330*f(k-4,:,:,j) &
                     +462*f(k-5,:,:,j) &
                     -462*f(k-6,:,:,j) &
                     +330*f(k-7,:,:,j) &
                     -165*f(k-8,:,:,j) &
                      +55*f(k-9,:,:,j) &
                      -11*f(k-10,:,:,j) &
                         +f(k-11,:,:,j)
          enddo
        elseif (idir==2) then
          do k=m2+1,m2+off
            f(:,k,:,j)=11*f(:,k-1,:,j) &
                      -55*f(:,k-2,:,j) &
                     +165*f(:,k-3,:,j) &
                     -330*f(:,k-4,:,j) &
                     +462*f(:,k-5,:,j) &
                     -462*f(:,k-6,:,j) &
                     +330*f(:,k-7,:,j) &
                     -165*f(:,k-8,:,j) &
                      +55*f(:,k-9,:,j) &
                      -11*f(:,k-10,:,j) &
                         +f(:,k-11,:,j)
          enddo
        elseif (idir==3) then
          do k=n2+1,n2+off
            f(:,:,k,j)=11*f(:,:,k-1,j) &
                      -55*f(:,:,k-2,j) &
                     +165*f(:,:,k-3,j) &
                     -330*f(:,:,k-4,j) &
                     +462*f(:,:,k-5,j) &
                     -462*f(:,:,k-6,j) &
                     +330*f(:,:,k-7,j) &
                     -165*f(:,:,k-8,j) &
                      +55*f(:,:,k-9,j) &
                      -11*f(:,:,k-10,j) &
                         +f(:,:,k-11,j)
          enddo
        endif
      endif

    endsubroutine set_ghosts_for_onesided_ders
!***********************************************************************
    subroutine bval_from_neumann_scl(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val is a constant.
!
!  27-jan-17/Ivan: coded
!
      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      real :: val

      integer :: k

      if (topbot=='bot') then
        if (idir==1) then
          k=l1
          f(l1,:,:,j) = (-val*2520.*dx+  25200.*f(k+1,:,:,j) &
                                      -  56700.*f(k+2,:,:,j) &
                                      + 100800.*f(k+3,:,:,j) &
                                      - 132300.*f(k+4,:,:,j) &
                                      + 127008.*f(k+5,:,:,j) &
                                      -  88200.*f(k+6,:,:,j) &
                                      +  43200.*f(k+7,:,:,j) &
                                      -  14175.*f(k+8,:,:,j) &
                                      +   2800.*f(k+9,:,:,j) &
                                      -    252.*f(k+10,:,:,j))/7381.
        elseif (idir==2) then
          k=m1
          f(:,m1,:,j) = (-val*2520.*dy+  25200.*f(:,k+1,:,j) &
                                      -  56700.*f(:,k+2,:,j) &
                                      + 100800.*f(:,k+3,:,j) &
                                      - 132300.*f(:,k+4,:,j) &
                                      + 127008.*f(:,k+5,:,j) &
                                      -  88200.*f(:,k+6,:,j) &
                                      +  43200.*f(:,k+7,:,j) &
                                      -  14175.*f(:,k+8,:,j) &
                                      +   2800.*f(:,k+9,:,j) &
                                      -    252.*f(:,k+10,:,j))/7381.
        elseif (idir==3) then
          k=n1
          f(:,:,n1,j) = (-val*2520.*dz+  25200.*f(:,:,k+1,j) &
                                      -  56700.*f(:,:,k+2,j) &
                                      + 100800.*f(:,:,k+3,j) &
                                      - 132300.*f(:,:,k+4,j) &
                                      + 127008.*f(:,:,k+5,j) &
                                      -  88200.*f(:,:,k+6,j) &
                                      +  43200.*f(:,:,k+7,j) &
                                      -  14175.*f(:,:,k+8,j) &
                                      +   2800.*f(:,:,k+9,j) &
                                      -    252.*f(:,:,k+10,j))/7381.
        endif
      else
        if (idir==1) then
          k=l2
          f(l2,:,:,j) = (val*2520.*dx+  25200.*f(k-1,:,:,j) &
                                     -  56700.*f(k-2,:,:,j) &
                                     + 100800.*f(k-3,:,:,j) &
                                     - 132300.*f(k-4,:,:,j) &
                                     + 127008.*f(k-5,:,:,j) &
                                     -  88200.*f(k-6,:,:,j) &
                                     +  43200.*f(k-7,:,:,j) &
                                     -  14175.*f(k-8,:,:,j) &
                                     +   2800.*f(k-9,:,:,j) &
                                     -    252.*f(k-10,:,:,j))/7381.
        elseif (idir==2) then
          k=m2
          f(:,m2,:,j) = (val*2520.*dy+  25200.*f(:,k-1,:,j) &
                                     -  56700.*f(:,k-2,:,j) &
                                     + 100800.*f(:,k-3,:,j) &
                                     - 132300.*f(:,k-4,:,j) &
                                     + 127008.*f(:,k-5,:,j) &
                                     -  88200.*f(:,k-6,:,j) &
                                     +  43200.*f(:,k-7,:,j) &
                                     -  14175.*f(:,k-8,:,j) &
                                     +   2800.*f(:,k-9,:,j) &
                                     -    252.*f(:,k-10,:,j))/7381.
        elseif (idir==3) then
          k=n2
          f(:,:,n2,j) = (val*2520.*dz+  25200.*f(:,:,k-1,j) &
                                     -  56700.*f(:,:,k-2,j) &
                                     + 100800.*f(:,:,k-3,j) &
                                     - 132300.*f(:,:,k-4,j) &
                                     + 127008.*f(:,:,k-5,j) &
                                     -  88200.*f(:,:,k-6,j) &
                                     +  43200.*f(:,:,k-7,j) &
                                     -  14175.*f(:,:,k-8,j) &
                                     +   2800.*f(:,:,k-9,j) &
                                     -    252.*f(:,:,k-10,j))/7381.
        endif
      endif

    endsubroutine bval_from_neumann_scl
!***********************************************************************
    subroutine bval_from_neumann_arr(f,topbot,j,idir,val)
!
!  Dummy.
!
      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val

      call fatal_error('bval_from_neumann_arr','Not implemented for 10th order.')

    endsubroutine bval_from_neumann_arr
!***********************************************************************
    subroutine bval_from_3rd_scl(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the 3rd kind BC d f/d x_i = val*f employing
!  one-sided difference formulae. val is a constant.
!
!  27-jan-17/Ivan: coded
!
      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      real :: val

      integer :: k

      if (topbot=='bot') then
        if (idir==1) then
          k=l1
          f(l1,:,:,j) = (   25200.*f(k+1,:,:,j) &
                         -  56700.*f(k+2,:,:,j) &
                         + 100800.*f(k+3,:,:,j) &
                         - 132300.*f(k+4,:,:,j) &
                         + 127008.*f(k+5,:,:,j) &
                         -  88200.*f(k+6,:,:,j) &
                         +  43200.*f(k+7,:,:,j) &
                         -  14175.*f(k+8,:,:,j) &
                         +   2800.*f(k+9,:,:,j) &
                         -    252.*f(k+10,:,:,j))/(7381.+val*2520.*dx)
        elseif (idir==2) then
          k=m1
          f(:,m1,:,j) = (   25200.*f(:,k+1,:,j) &
                         -  56700.*f(:,k+2,:,j) &
                         + 100800.*f(:,k+3,:,j) &
                         - 132300.*f(:,k+4,:,j) &
                         + 127008.*f(:,k+5,:,j) &
                         -  88200.*f(:,k+6,:,j) &
                         +  43200.*f(:,k+7,:,j) &
                         -  14175.*f(:,k+8,:,j) &
                         +   2800.*f(:,k+9,:,j) &
                         -    252.*f(:,k+10,:,j))/(7381.+val*2520.*dy)
        elseif (idir==3) then
          k=n1
          f(:,:,n1,j) = (   25200.*f(:,:,k+1,j) &
                         -  56700.*f(:,:,k+2,j) &
                         + 100800.*f(:,:,k+3,j) &
                         - 132300.*f(:,:,k+4,j) &
                         + 127008.*f(:,:,k+5,j) &
                         -  88200.*f(:,:,k+6,j) &
                         +  43200.*f(:,:,k+7,j) &
                         -  14175.*f(:,:,k+8,j) &
                         +   2800.*f(:,:,k+9,j) &
                         -    252.*f(:,:,k+10,j))/(7381.+val*2520.*dz)
        endif
      else
        if (idir==1) then
          k=l2
          f(l2,:,:,j) = (   25200.*f(k-1,:,:,j) &
                         -  56700.*f(k-2,:,:,j) &
                         + 100800.*f(k-3,:,:,j) &
                         - 132300.*f(k-4,:,:,j) &
                         + 127008.*f(k-5,:,:,j) &
                         -  88200.*f(k-6,:,:,j) &
                         +  43200.*f(k-7,:,:,j) &
                         -  14175.*f(k-8,:,:,j) &
                         +   2800.*f(k-9,:,:,j) &
                         -    252.*f(k-10,:,:,j))/(7381.-val*2520.*dx)
        elseif (idir==2) then
          k=m2
          f(:,m2,:,j) = (   25200.*f(:,k-1,:,j) &
                         -  56700.*f(:,k-2,:,j) &
                         + 100800.*f(:,k-3,:,j) &
                         - 132300.*f(:,k-4,:,j) &
                         + 127008.*f(:,k-5,:,j) &
                         -  88200.*f(:,k-6,:,j) &
                         +  43200.*f(:,k-7,:,j) &
                         -  14175.*f(:,k-8,:,j) &
                         +   2800.*f(:,k-9,:,j) &
                         -    252.*f(:,k-10,:,j))/(7381.-val*2520.*dy)
        elseif (idir==3) then
          k=n2
          f(:,:,n2,j) = (   25200.*f(:,:,k-1,j) &
                         -  56700.*f(:,:,k-2,j) &
                         + 100800.*f(:,:,k-3,j) &
                         - 132300.*f(:,:,k-4,j) &
                         + 127008.*f(:,:,k-5,j) &
                         -  88200.*f(:,:,k-6,j) &
                         +  43200.*f(:,:,k-7,j) &
                         -  14175.*f(:,:,k-8,j) &
                         +   2800.*f(:,:,k-9,j) &
                         -    252.*f(:,:,k-10,j))/(7381.-val*2520.*dz)
        endif
      endif

    endsubroutine bval_from_3rd_scl
!***********************************************************************
    subroutine bval_from_3rd_arr(f,topbot,j,idir,val,func)
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
      external :: func
!
      call fatal_error('bval_from_3rd_arr','not implemented')

    endsubroutine bval_from_3rd_arr
!***********************************************************************
    subroutine bval_from_4th_scl(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the 4th kind BC d^2 f/d x_i^2 = val*f employing
!  one-sided difference formulae. val is a constant.
!
!  27-jan-17/Ivan: coded
!
      real, dimension(mx,my,mz,*) :: f
      character(LEN=3) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val

      integer :: k

      if (topbot=='bot') then
        if (idir==1) then
          k=l1
          f(l1,:,:,j) = (-  972200.*f(k+1,:,:,j) &
                         + 2754450.*f(k+2,:,:,j) &
                         - 5232800.*f(k+3,:,:,j) &
                         + 7088550.*f(k+4,:,:,j) &
                         - 6932016.*f(k+5,:,:,j) &
                         + 4872700.*f(k+6,:,:,j) &
                         - 2407200.*f(k+7,:,:,j) &
                         +  794925.*f(k+8,:,:,j) &
                         -  157800.*f(k+9,:,:,j) &
                         +   14258.*f(k+10,:,:,j))/(-177133.+val*25200.*dx*dx)
        elseif (idir==2) then
          k=m1
          f(:,m1,:,j) = (-  972200.*f(:,k+1,:,j) &
                         + 2754450.*f(:,k+2,:,j) &
                         - 5232800.*f(:,k+3,:,j) &
                         + 7088550.*f(:,k+4,:,j) &
                         - 6932016.*f(:,k+5,:,j) &
                         + 4872700.*f(:,k+6,:,j) &
                         - 2407200.*f(:,k+7,:,j) &
                         +  794925.*f(:,k+8,:,j) &
                         -  157800.*f(:,k+9,:,j) &
                         +   14258.*f(:,k+10,:,j))/(-177133.+val*25200.*dy*dy)
        elseif (idir==3) then
          k=n1
          f(:,:,n1,j) = (-  972200.*f(:,:,k+1,j) &
                         + 2754450.*f(:,:,k+2,j) &
                         - 5232800.*f(:,:,k+3,j) &
                         + 7088550.*f(:,:,k+4,j) &
                         - 6932016.*f(:,:,k+5,j) &
                         + 4872700.*f(:,:,k+6,j) &
                         - 2407200.*f(:,:,k+7,j) &
                         +  794925.*f(:,:,k+8,j) &
                         -  157800.*f(:,:,k+9,j) &
                         +   14258.*f(:,:,k+10,j))/(-177133.+val*25200.*dz*dz)
        endif
      else
        if (idir==1) then
          k=l2
          f(l2,:,:,j) = (-  972200.*f(k-1,:,:,j) &
                         + 2754450.*f(k-2,:,:,j) &
                         - 5232800.*f(k-3,:,:,j) &
                         + 7088550.*f(k-4,:,:,j) &
                         - 6932016.*f(k-5,:,:,j) &
                         + 4872700.*f(k-6,:,:,j) &
                         - 2407200.*f(k-7,:,:,j) &
                         +  794925.*f(k-8,:,:,j) &
                         -  157800.*f(k-9,:,:,j) &
                         +   14258.*f(k-10,:,:,j))/(-177133.+val*25200.*dx*dx)
        elseif (idir==2) then
          k=m2
          f(:,m2,:,j) = (-  972200.*f(:,k-1,:,j) &
                         + 2754450.*f(:,k-2,:,j) &
                         - 5232800.*f(:,k-3,:,j) &
                         + 7088550.*f(:,k-4,:,j) &
                         - 6932016.*f(:,k-5,:,j) &
                         + 4872700.*f(:,k-6,:,j) &
                         - 2407200.*f(:,k-7,:,j) &
                         +  794925.*f(:,k-8,:,j) &
                         -  157800.*f(:,k-9,:,j) &
                         +   14258.*f(:,k-10,:,j))/(-177133.+val*25200.*dy*dy)
        elseif (idir==3) then
          k=n2
          f(:,:,n2,j) = (-  972200.*f(:,:,k-1,j) &
                         + 2754450.*f(:,:,k-2,j) &
                         - 5232800.*f(:,:,k-3,j) &
                         + 7088550.*f(:,:,k-4,j) &
                         - 6932016.*f(:,:,k-5,j) &
                         + 4872700.*f(:,:,k-6,j) &
                         - 2407200.*f(:,:,k-7,j) &
                         +  794925.*f(:,:,k-8,j) &
                         -  157800.*f(:,:,k-9,j) &
                         +   14258.*f(:,:,k-10,j))/(-177133.+val*25200.*dz*dz)
        endif
      endif

    endsubroutine bval_from_4th_scl
!***********************************************************************
endmodule Deriv
