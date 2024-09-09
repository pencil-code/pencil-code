! $Id$
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM integer, parameter :: nghost = 2
!
!***************************************************************
module Deriv
!
  use Messages, only: fatal_error, warning, not_implemented
  use Cdata
  use General, only: keep_compiler_quiet
!
  implicit none
!
  private
!
  include 'deriv.h'
!
  real :: der2_coef0, der2_coef1, der2_coef2
!
  contains
!
!***********************************************************************
    subroutine initialize_deriv
!
!  Initialize stencil coefficients (dummy routine)
!
    endsubroutine initialize_deriv
!***********************************************************************
    subroutine calc_coeffs_1( grid, coeffs )
!
!  dummy
!
      real, dimension(-0:1), intent(in) :: grid
      real, dimension(-1:1), intent(out) :: coeffs
!
      if (lroot) print*,'calc_coeffs_1 is not evaluated'
!
  endsubroutine calc_coeffs_1
!***********************************************************************
    subroutine der_main(f, k, df, j, ignoredx)
!
!  calculate derivative df_k/dx_j
!  accurate to 4th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: df
      integer, intent(in) :: j, k
      logical, intent(in), optional :: ignoredx
!
      real, parameter :: a = 1.0/12.0
      real, dimension(nx) :: fac
      real :: facs
      logical :: withdx
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
          df=fac*(+  8.0*(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k)) &
                  -      (f(l1+2:l2+2,m,n,k)-f(l1-2:l2-2,m,n,k)) )
        else
          df=0.
          if (ip<=5) print*, 'der_main: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (withdx) facs = a * dy_1(m)
          df=facs*(+  8.0*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                   -      (f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) )
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
          if (lcoarse_mn) then
            if (withdx) facs = a * dz_1(n) * nphis1(m)
            df=facs*(+ 8.0*(f(l1:l2,m,ninds(+1,m,n),k)-f(l1:l2,m,ninds(-1,m,n),k)) &
                     -     (f(l1:l2,m,ninds(+2,m,n),k)-f(l1:l2,m,ninds(-2,m,n),k)) )
          else
            if (withdx) facs = a * dz_1(n)
            df=facs*(+  8.0*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k)) &
                     -      (f(l1:l2,m,n+2,k)-f(l1:l2,m,n-2,k)) )
          endif
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension (mx), intent(in)  :: f
      real, dimension (nx), intent(out) :: df
!
      real, dimension (nx) :: fac
!
      if (nxgrid/=1) then
        fac=(1./12)*dx_1(l1:l2)
        df=fac*(+  8.0*(f(l1+1:l2+1)-f(l1-1:l2-1)) &
                -      (f(l1+2:l2+2)-f(l1-2:l2-2)) )
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension (mx), intent(in)  :: f
      real, dimension (nx), intent(out) :: df2
!
      real, dimension (nx) :: fac
!
      if (nxgrid/=1) then
        fac=(1./12)*dx_1(l1:l2)**2
        df2=fac*(- 30.0*f(l1:l2) &
                 + 16.0*(f(l1+1:l2+1)+f(l1-1:l2-1)) &
                 -      (f(l1+2:l2+2)+f(l1-2:l2-2)) )
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(out) :: df
!
      real, dimension (nz) :: fac
!
      if (nzgrid/=1) then
!MR: coarse case/spherical missing!
        fac=(1./12)*dz_1(n1:n2)
        df=fac*(+  8.0*(f(n1+1:n2+1)-f(n1-1:n2-1)) &
                -      (f(n1+2:n2+2)-f(n1-2:n2-2)) )
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(out) :: df2
!
      real, dimension (nz) :: fac
!
      if (nzgrid/=1) then
!MR: coarse case/spherical missing!
        fac=(1./12)*dz_1(n1:n2)**2
        df2=fac*(- 30.0*f(n1:n2) &
                 + 16.0*(f(n1+1:n2+1)+f(n1-1:n2-1)) &
                 -      (f(n1+2:n2+2)+f(n1-2:n2-2)) )
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension (mx,my,mz), intent(in) :: f
      real, dimension (:), intent(out) :: df
      integer, intent(in) :: j
!
      real, dimension (size(df)) :: fac
      real :: facs
      integer :: l1_,l2_,sdf
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=(1./12)*dx_1(l1:l2)
          df=fac*(+  8.0*(f(l1+1:l2+1,m,n)-f(l1-1:l2-1,m,n)) &
                  -     (f(l1+2:l2+2,m,n)-f(l1-2:l2-2,m,n)) )
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
            facs=(1./12)*dy_1(m)
            df=facs*(+  8.0*(f(l1_:l2_,m+1,n)-f(l1_:l2_,m-1,n)) &
                     -      (f(l1_:l2_,m+2,n)-f(l1_:l2_,m-2,n)) )
            if (lspherical_coords)   df=df*r1_mn
            if (lcylindrical_coords) df=df*rcyl_mn1
          else
            df=0.
            if (ip<=5) print*, 'der_other: Degenerate case in y-direction'
          endif
        elseif (j==3) then
          if (nzgrid/=1) then
            if (lcoarse_mn) then
              facs = (1./12) * dz_1(n) * nphis1(m)
              df=facs*(+  8.0*(f(l1_:l2_,m,ninds(+1,m,n))-f(l1_:l2_,m,ninds(-1,m,n))) &
                       -      (f(l1_:l2_,m,ninds(+2,m,n))-f(l1_:l2_,m,ninds(-2,m,n))) )
            else
              facs = (1./12) * dz_1(n)
              df=facs*(+  8.0*(f(l1_:l2_,m,n+1)-f(l1_:l2_,m,n-1)) &
                       -      (f(l1_:l2_,m,n+2)-f(l1_:l2_,m,n-2)) )
            endif
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      use General, only: itoa

      integer, intent(in) :: j
      real, dimension (:), intent(in) :: pencil
      real, dimension (:), intent(out) :: df
!
!  x-derivative
!
      if (j==1) then
        if (size(pencil)/=mx) &
          call fatal_error('der_pencil', &
                           'pencil must be of size mx for x derivative')
        df(l1:l2)=(1./12)*dx_1(l1:l2)*( &
            +  8.0*(pencil(l1+1:l2+1)-pencil(l1-1:l2-1)) &
            -      (pencil(l1+2:l2+2)-pencil(l1-2:l2-2)) )
!
!  y-derivative
!
      else if (j==2) then
        if (size(pencil)/=my) &
          call fatal_error('der_pencil', &
                           'pencil must be of size my for y derivative')
        df(m1:m2)=(1./12)*dy_1(m1:m2)*( &
            +  8.0*(pencil(m1+1:m2+1)-pencil(m1-1:m2-1)) &
            -      (pencil(m1+2:m2+2)-pencil(m1-2:m2-2)) )
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
        df(n1:n2)=(1./12)*dz_1(n1:n2)*( &
            +  8.0*(pencil(n1+1:n2+1)-pencil(n1-1:n2-1)) &
            -      (pencil(n1+2:n2+2)-pencil(n1-2:n2-2)) )
        if (lspherical_coords) df(n1:n2)=df(n1:n2)*(r1_mn(lglob)*sin1th(m))
      else
        call fatal_error('der_pencil','no such direction j='//trim(itoa(j)))
      endif
!
    endsubroutine der_pencil
!***********************************************************************
  subroutine distr_der(arr,idir,der,order)
!
!  Calculates 1st or 2nd derivative of a 1D array (of vectors, so 2nd dim for components),
!  which is distributed across the procs of its dimension.
!  At the moment only for z-direction (idir=IZBEAM=3).
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
    real, dimension(:,:), intent(in) :: arr
    real, dimension(:,:), intent(out) :: der
    integer, intent(in) :: idir
    integer, intent(in), optional :: order
!
    call not_implemented('distr_der','for 4th order')
    call keep_compiler_quiet(idir)
    call keep_compiler_quiet(arr)
    call keep_compiler_quiet(der)
!
  endsubroutine distr_der
!***********************************************************************
    subroutine der2_main(f,k,df2,j,lwo_line_elem)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!  if lwo_line_elem=T no metric coefficents ar multiplied in the denominator;
!  default: lwo_line_elem=F
!
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      use General, only: loptest

      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df2
      integer, intent(in) :: j,k
      logical, intent(in), optional :: lwo_line_elem
!
      real, dimension (nx) :: fac, df
      real, parameter :: der2_coef0=-30./12, der2_coef1=16./12, der2_coef2=-1./12
      real :: facs
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=dx_1(l1:l2)**2
          df2=fac*(der2_coef0* f(l1  :l2  ,m,n,k) &
                  +der2_coef1*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  +der2_coef2*(f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) )
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
                   +der2_coef2*(f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) )
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
          if (lcoarse_mn) then
            facs=dz_1(n)**2*nphis2(m)
            df2=facs*( der2_coef0* f(l1:l2,m,n  ,k) &
                      +der2_coef1*(f(l1:l2,m,ninds(+1,m,n),k)+f(l1:l2,m,ninds(-1,m,n),k)) &
                      +der2_coef2*(f(l1:l2,m,ninds(+2,m,n),k)+f(l1:l2,m,ninds(-2,m,n),k)) )
          else
            facs=dz_1(n)**2
            df2=facs*( der2_coef0* f(l1:l2,m,n  ,k) &
                      +der2_coef1*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                      +der2_coef2*(f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) )
          endif
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension (mx,my,mz), intent(in) :: f
      real, dimension (nx), intent(out) :: df2
      integer, intent(in) :: j
!
      real :: facs
      real, dimension (nx) :: fac, df
!
      if (j==1) then
        if (nxgrid/=1) then
          fac=(1./12)*dx_1(l1:l2)**2
          df2=fac*(- 30.0*f(l1:l2,m,n) &
                   + 16.0*(f(l1+1:l2+1,m,n)+f(l1-1:l2-1,m,n)) &
                   -      (f(l1+2:l2+2,m,n)+f(l1-2:l2-2,m,n)) )
          if (.not.lequidist(j)) then
            call der(f,df,j)
            df2=df2+dx_tilde(l1:l2)*df
          endif
        else
          df2=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          facs=(1./12)*dy_1(m)**2
          df2=facs*(- 30.0*f(l1:l2,m,n) &
                    + 16.0*(f(l1:l2,m+1,n)+f(l1:l2,m-1,n)) &
                    -      (f(l1:l2,m+2,n)+f(l1:l2,m-2,n)) )
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
!          if (lcoarse_mn) then
          facs=(1./12)*dz_1(n)**2
          df2=facs*(- 30.0*f(l1:l2,m,n) &
                    + 16.0*(f(l1:l2,m,n+1)+f(l1:l2,m,n-1)) &
                    -      (f(l1:l2,m,n+2)+f(l1:l2,m,n-2)) )
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      integer, intent(in) :: j
      real, dimension (:), intent(in) :: pencil
      real, dimension (:), intent(out) :: df2
!
      if (j==1) then
!
!  x-derivative
!
        if (size(pencil)/=mx) &
          call fatal_error('der2_pencil','pencil must be of size mx for x derivative')
        df2=(1./12)*dx_1(l1:l2)**2*(- 30.0*pencil(l1:l2) &
               + 16.0*(pencil(l1+1:l2+1)+pencil(l1-1:l2-1)) &
               -      (pencil(l1+2:l2+2)+pencil(l1-2:l2-2)) )
      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) &
          call fatal_error('der2_pencil','pencil must be of size my for y derivative')
!MR: spherical/cylindrical missing
        df2=(1./12)*dy_1(m1:m2)**2*(- 30.0*pencil(m1:m2) &
               + 16.0*(pencil(m1+1:m2+1)+pencil(m1-1:m2-1)) &
               -      (pencil(m1+2:m2+2)+pencil(m1-2:m2-2)) )
      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) &
          call fatal_error('der2_pencil','pencil must be of size mz for z derivative')
!MR: spherical/coarse missing
        df2=(1./12)*dz_1(n1:n2)**2*(- 30.0*pencil(n1:n2) &
               + 16.0*(pencil(n1+1:n2+1)+pencil(n1-1:n2-1)) &
               -      (pencil(n1+2:n2+2)+pencil(n1-2:n2-2)) )
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: j, k
      logical, intent(in), optional :: ignoredx
!
      real, dimension (nx) :: fac
      logical :: igndx
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
            fac=(1.0/2)
          else
            fac=(1.0/2)*1./dx**3
          endif
          df=fac*(-  2.0*(f(l1+1:l2+1,m,n,k)-f(l1-1:l2-1,m,n,k)) &
                  +      (f(l1+2:l2+2,m,n,k)-f(l1-2:l2-2,m,n,k)) )
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=(1.0/2)
          else
            fac=(1.0/2)*1./dy**3
          endif
          df=fac*(-  2.0*(f(l1:l2,m+1,n,k)-f(l1:l2,m-1,n,k)) &
                  +      (f(l1:l2,m+2,n,k)-f(l1:l2,m-2,n,k)) )
          if (lcylindrical_coords)   df=df*rcyl_mn1**3
!MR: spherical missing
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=(1.0/2)
          else
            fac=(1.0/2)*1./dz**3
          endif
          df=fac*(-  2.0*(f(l1:l2,m,n+1,k)-f(l1:l2,m,n-1,k)) &
                  +      (f(l1:l2,m,n+2,k)-f(l1:l2,m,n-2,k)) )
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: j, k
      logical, intent(in), optional :: ignoredx, upwind
!
      logical :: igndx
      real :: fac
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
            fac=1.0
          else
            fac=1.0/dx**4
          endif
          df=fac*(+  6.0* f(l1:l2,m,n,k) &
                  -  4.0*(f(l1+1:l2+1,m,n,k)+f(l1-1:l2-1,m,n,k)) &
                  +      (f(l1+2:l2+2,m,n,k)+f(l1-2:l2-2,m,n,k)) )
        else
          df=0.
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if (igndx) then
            fac=1.0
          else
            fac=1.0/dy**4
          endif
          df=fac*(+  6.0* f(l1:l2,m  ,n,k) &
                  -  4.0*(f(l1:l2,m+1,n,k)+f(l1:l2,m-1,n,k)) &
                  +      (f(l1:l2,m+2,n,k)+f(l1:l2,m-2,n,k)) )
          if (lcylindrical_coords)   df=df*rcyl_mn1**4
!MR: spherical missing
        else
          df=0.
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if (igndx) then
            fac=1.0
          else
            fac=1.0/dz**4
          endif
          df=fac*(+  6.0* f(l1:l2,m,n  ,k) &
                  -  4.0*(f(l1:l2,m,n+1,k)+f(l1:l2,m,n-1,k)) &
                  +      (f(l1:l2,m,n+2,k)+f(l1:l2,m,n-2,k)) )
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
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: j,k
      logical, intent(in), optional :: ignoredx
!
      logical :: igndx
      real, dimension (nx) :: fac
!
      call fatal_error('deriv_4th','der5 not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine der5
!***********************************************************************
    subroutine der6_main(f,k,df,j,ignoredx,upwind)
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
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      use General, only: loptest

      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: j, k
      logical, intent(in), optional :: ignoredx, upwind
!
      call fatal_error('deriv_4th','der6 not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine der6_main
!***********************************************************************
    subroutine der2_minmod(f,j,delfk,delfkp1,delfkm1,k)
!
!  Dummy routine
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: delfk,delfkp1,delfkm1
      integer, intent(in) :: j, k
!
      call fatal_error('der2_minmod','Not implemented for deriv_2nd')
      call keep_compiler_quiet(delfk)
      call keep_compiler_quiet(delfkp1)
      call keep_compiler_quiet(delfkm1)
!
    endsubroutine der2_minmod
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
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz), intent(in) :: f
      real, dimension (nx),intent(out) :: df
      integer, intent(in) :: j
      logical, intent(in), optional :: ignoredx, upwind
!
      call fatal_error('deriv_4th','der6_other not implemented yet')
      call keep_compiler_quiet(df)
!
    endsubroutine der6_other
!***********************************************************************
    subroutine der6_pencil(j,pencil,df6,ignoredx,upwind)
!
!  Calculate 6th derivative of any x, y or z pencil.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (:), intent(in) :: pencil
      real, dimension (:), intent(out) :: df6
      integer, intent(in) :: j
      logical, intent(in), optional :: ignoredx, upwind
!
      call fatal_error('deriv_4th','der6_pencil not implemented yet')
      call keep_compiler_quiet(df6)
!
    endsubroutine der6_pencil
!***********************************************************************
    real function der5_single(f,j,dc1)
!
!  computes 5th order derivative of function given by f at position j
!
!  09-Sep-2024/PABourdin: adapted from 2nd-order
!
      real, dimension(:), intent(in) :: f, dc1
      integer, intent(in) :: j
!
      call fatal_error('deriv_4th','der6_pencil not implemented yet')
      call keep_compiler_quiet(der5_single)
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
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      use General, only: loptest
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: i, j, k
      logical, intent(in), optional :: lwo_line_elem
!
      real :: facs
      real, dimension (nx) :: fac
!
! crash if this is called with i=j
!
!      if (i.eq.j) call fatal_error('derij_main','i=j, no derivative calculated')
!
      if (lbidiagonal_derij .and. .not.lcoarse_mn) then
        !
        ! Use bidiagonal mixed-derivative operator, i.e.
        ! employ only the three neighbouring points on each of the four
        ! half-diagonals. This gives 6th-order mixed derivatives as the
        ! version below, but involves just 12 points instead of 36.
        !
        if (i+j==3) then
          if (nxgrid/=1.and.nygrid/=1) then
            fac=(1./144.)*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
                         64.*( f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k)  &
                              +f(l1-1:l2-1,m-1,n,k)-f(l1+1:l2+1,m-1,n,k)) &
                       -     ( f(l1+2:l2+2,m+2,n,k)-f(l1-2:l2-2,m+2,n,k)  &
                              +f(l1-2:l2-2,m-2,n,k)-f(l1+2:l2+2,m-2,n,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in x- or y-direction'
          endif
        elseif (i+j==5) then
          if (nygrid/=1.and.nzgrid/=1) then
            facs=(1./144.)*dy_1(m)*dz_1(n)
            df=facs*( &
                         64.*( f(l1:l2,m+1,n+1,k)-f(l1:l2,m+1,n-1,k)  &
                              +f(l1:l2,m-1,n-1,k)-f(l1:l2,m-1,n+1,k)) &
                       -     ( f(l1:l2,m+2,n+2,k)-f(l1:l2,m+2,n-2,k)  &
                              +f(l1:l2,m-2,n-2,k)-f(l1:l2,m-2,n+2,k)) &
                    )
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./144.)*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
                         64.*( f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k)  &
                              +f(l1-1:l2-1,m,n-1,k)-f(l1+1:l2+1,m,n-1,k)) &
                       -     ( f(l1+2:l2+2,m,n+2,k)-f(l1-2:l2-2,m,n+2,k)  &
                              +f(l1-2:l2-2,m,n-2,k)-f(l1+2:l2+2,m,n-2,k)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in x- or z-direction'
          endif
        endif
!
      else                      ! not using bidiagonal mixed derivatives
        !
        ! This is the old, straight-forward scheme
        !
        if (i+j==3) then
          if (nxgrid/=1.and.nygrid/=1) then
            fac=(1./144.)*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
               8.*(( 8.*(f(l1+1:l2+1,m+1,n,k)-f(l1-1:l2-1,m+1,n,k))  &
                    -   (f(l1+2:l2+2,m+1,n,k)-f(l1-2:l2-2,m+1,n,k))) &
                  -( 8.*(f(l1+1:l2+1,m-1,n,k)-f(l1-1:l2-1,m-1,n,k))  &
                    -   (f(l1+2:l2+2,m-1,n,k)-f(l1-2:l2-2,m-1,n,k))))&
              -   (( 8.*(f(l1+1:l2+1,m+2,n,k)-f(l1-1:l2-1,m+2,n,k))  &
                    -   (f(l1+2:l2+2,m+2,n,k)-f(l1-2:l2-2,m+2,n,k))) &
                  -( 8.*(f(l1+1:l2+1,m-2,n,k)-f(l1-1:l2-1,m-2,n,k))  &
                    -   (f(l1+2:l2+2,m-2,n,k)-f(l1-2:l2-2,m-2,n,k))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in x- or y-direction'
          endif
        elseif (i+j==5) then
          if (nygrid/=1.and.nzgrid/=1) then
            facs=(1./144.)*dy_1(m)*dz_1(n)
            if (lcoarse_mn) then
              facs=facs*nphis1(m)
              df=facs*( &
                 8.*(( 8.*(f(l1:l2,m+1,ninds(+1,m,n),k)-f(l1:l2,m-1,ninds(+1,m,n),k))  &
                      -   (f(l1:l2,m+2,ninds(+1,m,n),k)-f(l1:l2,m-2,ninds(+1,m,n),k))) &
                    -( 8.*(f(l1:l2,m+1,ninds(-1,m,n),k)-f(l1:l2,m-1,ninds(-1,m,n),k))  &
                      -   (f(l1:l2,m+2,ninds(-1,m,n),k)-f(l1:l2,m-2,ninds(-1,m,n),k))))&
                -   (( 8.*(f(l1:l2,m+1,ninds(+2,m,n),k)-f(l1:l2,m-1,ninds(+2,m,n),k))  &
                      -   (f(l1:l2,m+2,ninds(+2,m,n),k)-f(l1:l2,m-2,ninds(+2,m,n),k))) &
                    -( 8.*(f(l1:l2,m+1,ninds(-2,m,n),k)-f(l1:l2,m-1,ninds(-2,m,n),k))  &
                      -   (f(l1:l2,m+2,ninds(-2,m,n),k)-f(l1:l2,m-2,ninds(-2,m,n),k))))&
                    )
            else
              df=facs*( &
                 8.*(( 8.*(f(l1:l2,m+1,n+1,k)-f(l1:l2,m-1,n+1,k))  &
                      -   (f(l1:l2,m+2,n+1,k)-f(l1:l2,m-2,n+1,k))) &
                    -( 8.*(f(l1:l2,m+1,n-1,k)-f(l1:l2,m-1,n-1,k))  &
                      -   (f(l1:l2,m+2,n-1,k)-f(l1:l2,m-2,n-1,k))))&
                -   (( 8.*(f(l1:l2,m+1,n+2,k)-f(l1:l2,m-1,n+2,k))  &
                      -   (f(l1:l2,m+2,n+2,k)-f(l1:l2,m-2,n+2,k))) &
                    -( 8.*(f(l1:l2,m+1,n-2,k)-f(l1:l2,m-1,n-2,k))  &
                      -   (f(l1:l2,m+2,n-2,k)-f(l1:l2,m-2,n-2,k))))&
                    )
            endif
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./60.**2)*dz_1(n)*dx_1(l1:l2)
            if (lcoarse_mn) then
              fac=fac*nphis1(m)
              df=fac*( &
                 8.*(( 8.*(f(l1+1:l2+1,m,ninds(+1,m,n),k)-f(l1-1:l2-1,m,ninds(+1,m,n),k))  &
                      -   (f(l1+2:l2+2,m,ninds(+1,m,n),k)-f(l1-2:l2-2,m,ninds(+1,m,n),k))) &
                    -( 8.*(f(l1+1:l2+1,m,ninds(-1,m,n),k)-f(l1-1:l2-1,m,ninds(-1,m,n),k))  &
                      -   (f(l1+2:l2+2,m,ninds(-1,m,n),k)-f(l1-2:l2-2,m,ninds(-1,m,n),k))))&
                -   (( 8.*(f(l1+1:l2+1,m,ninds(+2,m,n),k)-f(l1-1:l2-1,m,ninds(+2,m,n),k))  &
                      -   (f(l1+2:l2+2,m,ninds(+2,m,n),k)-f(l1-2:l2-2,m,ninds(+2,m,n),k))) &
                    -( 8.*(f(l1+1:l2+1,m,ninds(-2,m,n),k)-f(l1-1:l2-1,m,ninds(-2,m,n),k))  &
                      -   (f(l1+2:l2+2,m,ninds(-2,m,n),k)-f(l1-2:l2-2,m,ninds(-2,m,n),k))))&
                     )
            else
              df=fac*( &
                 8.*(( 8.*(f(l1+1:l2+1,m,n+1,k)-f(l1-1:l2-1,m,n+1,k))  &
                      -   (f(l1+2:l2+2,m,n+1,k)-f(l1-2:l2-2,m,n+1,k))) &
                    -( 8.*(f(l1+1:l2+1,m,n-1,k)-f(l1-1:l2-1,m,n-1,k))  &
                      -   (f(l1+2:l2+2,m,n-1,k)-f(l1-2:l2-2,m,n-1,k))))&
                -   (( 8.*(f(l1+1:l2+1,m,n+2,k)-f(l1-1:l2-1,m,n+2,k))  &
                      -   (f(l1+2:l2+2,m,n+2,k)-f(l1-2:l2-2,m,n+2,k))) &
                    -( 8.*(f(l1+1:l2+1,m,n-2,k)-f(l1-1:l2-1,m,n-2,k))  &
                      -   (f(l1+2:l2+2,m,n-2,k)-f(l1-2:l2-2,m,n-2,k))))&
                     )
            endif
          else
            df=0.
            if (ip<=5) print*, 'derij_main: Degenerate case in x- or z-direction'
          endif
        endif
!
      endif
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
        if ( i==2 .or. j==2 ) df=df*rcyl_mn1        ! not correct for i=j
      endif
!
    endsubroutine derij_main
!***********************************************************************
    subroutine derij_other(f,df,i,j)
!
!  calculate 2nd derivative with respect to two different directions
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension (mx,my,mz), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: i, j
!
      real, dimension (nx) :: fac
!
! crash if this is called with i=j
!
!      if (i.eq.j) call fatal_error('derij_main','i=j, no derivative calculated')
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
            fac=(1./144.)*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
                         64.*( f(l1+1:l2+1,m+1,n)-f(l1-1:l2-1,m+1,n)  &
                              +f(l1-1:l2-1,m-1,n)-f(l1+1:l2+1,m-1,n)) &
                       -     ( f(l1+2:l2+2,m+2,n)-f(l1-2:l2-2,m+2,n)  &
                              +f(l1-2:l2-2,m-2,n)-f(l1+2:l2+2,m-2,n)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif (i+j==5) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=(1./144.)*dy_1(m)*dz_1(n)
            df=fac*( &
                         64.*( f(l1:l2,m+1,n+1)-f(l1:l2,m+1,n-1)  &
                              +f(l1:l2,m-1,n-1)-f(l1:l2,m-1,n+1)) &
                       -     ( f(l1:l2,m+2,n+2)-f(l1:l2,m+2,n-2)  &
                              +f(l1:l2,m-2,n-2)-f(l1:l2,m-2,n+2)) &
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./144.)*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
                         64.*( f(l1+1:l2+1,m,n+1)-f(l1-1:l2-1,m,n+1)  &
                              +f(l1-1:l2-1,m,n-1)-f(l1+1:l2+1,m,n-1)) &
                       -     ( f(l1+2:l2+2,m,n+2)-f(l1-2:l2-2,m,n+2)  &
                              +f(l1-2:l2-2,m,n-2)-f(l1+2:l2+2,m,n-2)) &
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
            fac=(1./144)*dx_1(l1:l2)*dy_1(m)
            df=fac*( &
               8.*(( 8.*(f(l1+1:l2+1,m+1,n)-f(l1-1:l2-1,m+1,n))  &
                    -   (f(l1+2:l2+2,m+1,n)-f(l1-2:l2-2,m+1,n))) &
                  -( 8.*(f(l1+1:l2+1,m-1,n)-f(l1-1:l2-1,m-1,n))  &
                    -   (f(l1+2:l2+2,m-1,n)-f(l1-2:l2-2,m-1,n))))&
              -   (( 8.*(f(l1+1:l2+1,m+2,n)-f(l1-1:l2-1,m+2,n))  &
                    -   (f(l1+2:l2+2,m+2,n)-f(l1-2:l2-2,m+2,n))) &
                  -( 8.*(f(l1+1:l2+1,m-2,n)-f(l1-1:l2-1,m-2,n))  &
                    -   (f(l1+2:l2+2,m-2,n)-f(l1-2:l2-2,m-2,n))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif (i+j==5) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=(1./144)*dy_1(m)*dz_1(n)
            df=fac*( &
               8.*(( 8.*(f(l1:l2,m+1,n+1)-f(l1:l2,m-1,n+1))  &
                    -   (f(l1:l2,m+2,n+1)-f(l1:l2,m-2,n+1))) &
                  -( 8.*(f(l1:l2,m+1,n-1)-f(l1:l2,m-1,n-1))  &
                    -   (f(l1:l2,m+2,n-1)-f(l1:l2,m-2,n-1))))&
              -   (( 8.*(f(l1:l2,m+1,n+2)-f(l1:l2,m-1,n+2))  &
                    -   (f(l1:l2,m+2,n+2)-f(l1:l2,m-2,n+2))) &
                  -( 8.*(f(l1:l2,m+1,n-2)-f(l1:l2,m-1,n-2))  &
                    -   (f(l1:l2,m+2,n-2)-f(l1:l2,m-2,n-2))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==1).or.(i==1.and.j==3)) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./144)*dz_1(n)*dx_1(l1:l2)
            df=fac*( &
               8.*(( 8.*(f(l1+1:l2+1,m,n+1)-f(l1-1:l2-1,m,n+1))  &
                    -   (f(l1+2:l2+2,m,n+1)-f(l1-2:l2-2,m,n+1))) &
                  -( 8.*(f(l1+1:l2+1,m,n-1)-f(l1-1:l2-1,m,n-1))  &
                    -   (f(l1+2:l2+2,m,n-1)-f(l1-2:l2-2,m,n-1))))&
              -   (( 8.*(f(l1+1:l2+1,m,n+2)-f(l1-1:l2-1,m,n+2))  &
                    -   (f(l1+2:l2+2,m,n+2)-f(l1-2:l2-2,m,n+2))) &
                  -( 8.*(f(l1+1:l2+1,m,n-2)-f(l1-1:l2-1,m,n-2))  &
                    -   (f(l1+2:l2+2,m,n-2)-f(l1-2:l2-2,m,n-2))))&
                   )
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or z-direction'
          endif
        endif
!
      endif
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
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(out) :: df
      integer, intent(in) :: i, j, k
!
      call fatal_error('der5i1j','not implemented in deriv_4th')
      call keep_compiler_quiet(df)
!
    endsubroutine der5i1j
!***********************************************************************
    subroutine der4i2j(f,k,df,i,j)
!
!  Calculate 6th derivative with respect to two different directions.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx), intent(in) :: df
      integer, intent(in) :: i, j, k
!
      call fatal_error("der4i2j","not implemented in deriv_4th")
      call keep_compiler_quiet(df)
!
    endsubroutine der4i2j
!***********************************************************************
    subroutine der2i2j2k(f,k,df)
!
!  Mixed 6th derivative of der2x(der2y(der2z(f))). Worked out symbolically
!  in python. Result as spit from the python routine.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz,mfarray),intent(in) :: f
      real, dimension (nx) :: fac
      integer,intent(in) :: k
      real, dimension(nx), intent(out) :: df
!
      call fatal_error("der2i2j2k","not implemented in deriv_4th")
      call keep_compiler_quiet(df)
!
    endsubroutine der2i2j2k
!***********************************************************************
    subroutine der3i3j(f,k,df,i,j)
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx), intent(out) :: df
      real, dimension (nx) :: fac
      integer, intent(in) :: k,i,j
!
      call fatal_error("der3i3j","not implemented in deriv_4th")
      call keep_compiler_quiet(df)
!
    endsubroutine der3i3j
!***********************************************************************          
    subroutine der3i2j1k(f,ik,df,i,j,k)
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx), intent(out) :: df
      real, dimension (nx) :: fac
      integer, intent(in) :: ik,i,j,k
!
      call fatal_error("der3i2j1k","not implemented in deriv_4th")
      call keep_compiler_quiet(df)
!
    endsubroutine der3i2j1k
!***********************************************************************
    subroutine der4i1j1k(f,ik,df,i,j,k)
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx), intent(out) :: df
      real, dimension (nx) :: fac
      integer, intent(in) :: ik,i,j,k
!
      call fatal_error("der4i1j1k","not implemented in deriv_10th")
      call keep_compiler_quiet(df)
!
    endsubroutine der4i1j1k
!***********************************************************************
    subroutine der_upwind1st(f,uu,k,df,j)
!
!  First order upwind derivative of variable
!  Useful for advecting non-logarithmic variables
!
!  09-Sep-2024/PABourdin: adapted from 6th-order
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: df
      integer :: j,k,l
!
      intent(in)  :: f,uu,k,j
      intent(out) :: df
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
!  Calculate x/y/z-derivative on a yz/xz/xy-slice at gridpoint pos.
!  Uses a one-sided 4th order stencil.
!  sgn = +1 for forward difference, sgn = -1 for backwards difference.
!
!  Because of its original intended use in relation to solving
!  characteristic equations on boundaries (NSCBC), this sub should
!  return only PARTIAL derivatives, NOT COVARIANT. Applying the right
!  scaling factors and connection terms should instead be done when
!  solving the characteristic equations.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:,:) :: df
      real :: fac
      integer :: pos,k,sgn,j
!
      intent(in)  :: f,k,pos,sgn,j
      intent(out) :: df
!
      call fatal_error('der_onesided_4_slice_main','not implemented in deriv_4th')
      call keep_compiler_quiet(df)
!
    endsubroutine der_onesided_4_slice_main
!***********************************************************************
    subroutine der_onesided_4_slice_main_pt(f,sgn,k,df,lll,mmm,nnn,j)
!
!  made using der_onesided_4_slice_main. One sided derivative is calculated
!  at one point (lll,mmm,nnn)
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz,mfarray) :: f
      real  :: df
      real :: fac
      integer :: pos,lll,mmm,nnn,k,sgn,j
!
      intent(in)  :: f,k,lll,mmm,nnn,sgn,j
      intent(out) :: df
!
      call fatal_error('der_onesided_4_slice_other','not implemented in deriv_4th')
      call keep_compiler_quiet(df)
!
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
!   09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz) :: f
      real, dimension (:,:) :: df
      real :: fac
      integer :: pos,sgn,j
!
      intent(in)  :: f,pos,sgn,j
      intent(out) :: df
!
      call fatal_error('der_onesided_4_slice_other','not implemented in deriv_4th')
      call keep_compiler_quiet(df)
!
    endsubroutine der_onesided_4_slice_other
!***********************************************************************
    subroutine der_onesided_4_slice_other_pt(f,sgn,df,lll,mmm,nnn,j)
!
!  One sided derivative is calculated
!  at one point (lll,mmm,nnn).
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension (mx,my,mz) :: f
      real :: df
      real :: fac
      integer :: pos,lll,mmm,nnn,sgn,j
!
      intent(in)  :: f,lll,mmm,nnn,sgn,j
      intent(out) :: df
!
      call fatal_error('der_onesided_4_slice_main_pt','not implemented in deriv_4th')
      call keep_compiler_quiet(df)
!
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
!  09-Sep-2024/PABourdin: copied from 2nd-order
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
      call keep_compiler_quiet(df)
!
    endsubroutine deri_3d_inds
!************************************************************************
    logical function heatflux_deriv_x(f, inh, fac, topbot)
!
!   dummy routine
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension(mx,my,mz,mfarray), intent(IN):: f
      real, dimension(my,mz)           , intent(IN):: inh
      real                             , intent(IN):: fac
      integer                          , intent(IN):: topbot
!
      heatflux_deriv_x = .false.
!
    endfunction heatflux_deriv_x
!***********************************************************************
    subroutine set_ghosts_for_onesided_ders(f,topbot,j,idir,l2nd_)
!
!  Dummy.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      logical, optional :: l2nd_
!
      call fatal_error('set_ghosts_for_onesided_ders','Not implemented for deriv_4th.')
!
    endsubroutine set_ghosts_for_onesided_ders
!***********************************************************************
    subroutine bval_from_neumann_scl(f,topbot,j,idir,val)
!
!  Dummy.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real :: val
!
      call fatal_error('bval_from_neumann_scl','Not implemented for deriv_4th.')
!
    endsubroutine bval_from_neumann_scl
!***********************************************************************
    subroutine bval_from_neumann_arr(f,topbot,j,idir,val)
!
!  Dummy.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val
!
      call fatal_error('bval_from_neumann_arr','Not implemented for deriv_4th.')
!
    endsubroutine bval_from_neumann_arr
!***********************************************************************
    subroutine bval_from_3rd_scl(f,topbot,j,idir,val)
!
!  Dummy.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real :: val
!
      call fatal_error('bval_from_3rd_scl','Not implemented for deriv_4th.')
!
    endsubroutine bval_from_3rd_scl
!***********************************************************************
    subroutine bval_from_3rd_arr(f,topbot,j,idir,val,func)
!
!  Calculates the boundary value from the Neumann BC d f/d x_i = val employing
!  one-sided difference formulae. val depends on x,y.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val
      external :: func
!
      call fatal_error('bval_from_3rd_arr','not implemented for deriv_4th')
!
    endsubroutine bval_from_3rd_arr
!***********************************************************************
    subroutine bval_from_4th_scl(f,topbot,j,idir,val)
!
!  Dummy.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real :: val
!
      call fatal_error('bval_from_4th_scl','Not implemented for deriv_4th.')
!
    endsubroutine bval_from_4th_scl
!***********************************************************************
    subroutine bval_from_4th_arr(f,topbot,j,idir,val)
!
!  Calculates the boundary value from the 4th kind BC d^2 f/d x_i^2 = val*f
!  employing one-sided difference formulae. val depends on x,y.
!
!  09-Sep-2024/PABourdin: copied from 2nd-order
!
      real, dimension(mx,my,mz,*) :: f
      integer, intent(IN) :: topbot
      integer :: j,idir
      real, dimension(:,:) :: val
!
      call fatal_error('bval_from_4th_arr','Not implemented for deriv_4th.')
!
    endsubroutine bval_from_4th_arr
!***********************************************************************
 endmodule Deriv
