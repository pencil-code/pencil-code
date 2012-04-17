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
! Experimental module with explicit difference formulae for nonequidistant grids
! not yet fully tested and documented
!
! 26-Mar-12/MR: derived from deriv.f90
!
  use Messages, only: fatal_error, warning
  use Cdata
!
  implicit none
!
  private
!
  public :: initialize_deriv, finalize_deriv
  public :: deri_3d_inds
  public :: der, der2, der3, der4, der5, der6, derij, der5i1j
  public :: der6_other, der_pencil, der2_pencil
  public :: der_z,der2_z
  public :: der_upwind1st
  public :: der_onesided_4_slice
  public :: der_onesided_4_slice_other
  public :: der2_minmod
  public :: heatflux_boundcond_x
!
  real, dimension( 4,6) :: der_coef
  real, dimension(nx,6) :: r1i
  real, dimension(my,6) :: sth1i
!
  real, dimension(:,:,:), allocatable :: coeffsx, coeffsy, coeffsz, coeffsx_1, coeffsy_1, coeffsz_1
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
  interface  der_onesided_4_slice                  ! Overload the der function
    module procedure  der_onesided_4_slice_main    ! derivative of an 'mvar' variable
    module procedure  der_onesided_4_slice_main_pt
    module procedure  der_onesided_4_slice_other   ! derivative of another field
    module procedure  der_onesided_4_slice_other_pt
  endinterface

!  interface deri_3d                  !not working, why?             
!    module procedure deri_3d_ind
!    module procedure deri_3d_inds
!  endinterface  
!
  contains
!
!***********************************************************************
    subroutine initialize_deriv()
!
!  Initialize stencil coefficients
!
      integer :: i

      do i=1,6 

        if (lspherical_coords) then                         !-> grid.f90
          r1i  (:,i) = r1_mn**i
          sth1i(:,i) = sin1th**i
        endif
        if (lcylindrical_coords) r1i(:,i) = rcyl_mn1**i
        
      enddo

      der_coef(:,1) = (/  0.   ,   3./4., -3./20.,  1./60./)
      der_coef(:,3) = (/  0.   , -13./8.,  1.    , -1./8. /)
      der_coef(:,4) = (/ 28./3., -13./2.,  2.    , -1./6. /)
      der_coef(:,5) = (/  0.   ,   2.5  , -2.    ,  0.5   /)
      der_coef(:,6) = (/-20.   ,  15.   , -6.    ,  1.    /)
                                                           
      select case (der2_type)
!
      case ('standard')
        der_coef(:,2) = (/ -49./18., 3./2.,-3./20.,1./90. /)
!
      case ('tuned1')
        der_coef(:,2) = (/ -0.75, 0.34375, 0.125, -0.09375 /)  !tbc
!
      case default
        call fatal_error('initialize_deriv',"der2_type doesn't exist")
!
      endselect

      if ( lequidist(1) ) then

        if (.not.allocated(coeffsx)) allocate( coeffsx(-3:3,8,2) )
        
        do i=1,6
          coeffsx( 0:3,i,1) = der_coef(:,i)*dx_1(1)**i
          coeffsx( 0:3,i,2) = der_coef(:,i)
        enddo

      else

        if (.not.allocated(coeffsx)) allocate( coeffsx(-3:3,8,nx) )
        if (.not.allocated(coeffsx_1)) allocate( coeffsx_1(-2:2,2,2) )

        do i=l1,l2
          call calc_coeffs( x(i-nghost+1:i+nghost)-x(i-nghost:i+nghost-1), coeffsx(-3,1,i-l1+1) )   !use dx_1 ?
        enddo
!        call test( x, dx_1, nx, coeffsx )

        call calc_coeffs_1( x(l1-1:l1+2)-x(l1-2:l1+1), coeffsx_1(-2,1,1) )
        call calc_coeffs_1( x(l2-1:l2+2)-x(l2-2:l2+1), coeffsx_1(-2,1,2) )

      endif

      if ( lequidist(2) ) then

        if (.not.allocated(coeffsy)) allocate( coeffsy(-3:3,8,2) )
        do i=1,6
          coeffsy( 0:3,i,1) = der_coef(:,i)*dy_1(1)**i
          coeffsy( 0:3,i,2) = der_coef(:,i)
        enddo

      else

        if (.not.allocated(coeffsy)) allocate( coeffsy(-3:3,8,ny) )
        if (.not.allocated(coeffsy_1)) allocate( coeffsy_1(-2:2,2,2) )

        do i=m1,m2
          call calc_coeffs( y(i-nghost+1:i+nghost)-y(i-nghost:i+nghost-1), coeffsy(-3,1,i-m1+1) )
        enddo

        call calc_coeffs_1( y(m1-1:m1+2)-y(m1-2:m1+1), coeffsy_1(-2,1,1) )
        call calc_coeffs_1( y(m2-1:m2+2)-y(m2-2:m2+1), coeffsy_1(-2,1,2) )

      endif

      if ( lequidist(3) ) then

        if (.not.allocated(coeffsz)) allocate( coeffsz(-3:3,8,2) )

        do i=1,6
          coeffsz( 0:3,i,1) = der_coef(:,i)*dz_1(1)**i
          coeffsz( 0:3,i,2) = der_coef(:,i)
        enddo

      else

        if (.not.allocated(coeffsz)) allocate( coeffsz(-3:3,8,nz) )
        if (.not.allocated(coeffsz_1)) allocate( coeffsz_1(-2:2,2,2) )

        do i=n1,n2
          call calc_coeffs( z(i-nghost+1:i+nghost)-z(i-nghost:i+nghost-1), coeffsz(-3,1,i-n1+1) )
        enddo
!        call test( z, dz_1, nz, coeffsz )
     
        call calc_coeffs_1( z(n1-1:n1+2)-z(n1-2:n1+1), coeffsz_1(-2,1,1) )
        call calc_coeffs_1( z(n2-1:n2+2)-z(n2-2:n2+1), coeffsz_1(-2,1,2) )

!        call test_1( z, dz_1, nz, coeffsz_1 )
      endif
!
    endsubroutine initialize_deriv
!***********************************************************************
    subroutine finalize_deriv()

      deallocate( coeffsx, coeffsy, coeffsz )

    endsubroutine finalize_deriv
!***********************************************************************
    subroutine der_main(f,k,df,j)
!
!  calculates derivative df_k/dx_j
!  accurate to 6th order, explicit, periodic (??)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: df
      integer :: j,k
!
      intent(in)  :: f,k,j
      intent(out) :: df
!
      call deri_3d( f(1,1,1,k), df, 1, j )
!
    endsubroutine der_main
!***********************************************************************
    subroutine der_other(f,df,j)
!
!  Along one pencil in NON f variable
!  calculate derivative of a scalar, get scalar
!  accurate to 6th order, explicit, periodic
!
!  26-nov-02/tony: coded - duplicate der_main but without k subscript
!                          then overload the der interface.
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  21-feb-07/axel: added 1/r and 1/pomega factors for non-coord basis
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx) :: df
      integer :: j
!
      intent(in)  :: f,j
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(1,icount_der_other,j,1) = &
!debug                          der_call_count(1,icount_der_other,j,1) + 1
!
      call deri_3d(f,df,1,j)
!
    endsubroutine der_other
!***********************************************************************
    subroutine deri_3d(f,df,i,j,lignored,lnometric)    
!
!  Along one pencil in NON f variable
!  calculate derivative of a scalar, get scalar
!  accurate to 6th order, explicit, periodic
!
!  26-nov-02/tony: coded - duplicate der_main but without k subscript
!                          then overload the der interface.
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  21-feb-07/axel: added 1/r and 1/pomega factors for non-cart basis
!
      real, dimension (mx,my,mz)          :: f
      real, dimension (nx)                :: df
      integer                             :: i,j
      logical,                   optional :: lignored, lnometric
!
      intent(in)  :: f,i,j,lignored,lnometric
      intent(out) :: df
!
      integer :: mm, nn, is
      logical :: lmetric
!
!debug      if (loptimise_ders) der_call_count(1,icount_der_other,j,1) = &
!debug                          der_call_count(1,icount_der_other,j,1) + 1
!
      is = 1
      if ( present(lignored) ) then 
        if (lignored) is=2 
      endif

      if (present(lnometric)) then
        lmetric = .not.lnometric
      else
        lmetric = .true.
      endif

      if (j==1) then
        if (nxgrid/=1) then
          if ( lequidist(j) ) then
            if ( mod(i,2) == 0 ) then 
              df =  coeffsx(0,i,is)* f(l1  :l2  ,m,n)                   &    
                  + coeffsx(1,i,is)*(f(l1+1:l2+1,m,n)+f(l1-1:l2-1,m,n)) &
                  + coeffsx(2,i,is)*(f(l1+2:l2+2,m,n)+f(l1-2:l2-2,m,n)) &
                  + coeffsx(3,i,is)*(f(l1+3:l2+3,m,n)+f(l1-3:l2-3,m,n))
            else
              df =  coeffsx(1,i,is)*(f(l1+1:l2+1,m,n)-f(l1-1:l2-1,m,n)) &
                  + coeffsx(2,i,is)*(f(l1+2:l2+2,m,n)-f(l1-2:l2-2,m,n)) &
                  + coeffsx(3,i,is)*(f(l1+3:l2+3,m,n)-f(l1-3:l2-3,m,n))
            endif
          else 
            df =  coeffsx(-3,i,:)*f(l1-3:l2-3,m,n) & 
                + coeffsx(-2,i,:)*f(l1-2:l2-2,m,n) &
                + coeffsx(-1,i,:)*f(l1-1:l2-1,m,n) &
                + coeffsx( 0,i,:)*f(l1  :l2  ,m,n) &
                + coeffsx( 1,i,:)*f(l1+1:l2+1,m,n) &
                + coeffsx( 2,i,:)*f(l1+2:l2+2,m,n) &
                + coeffsx( 3,i,:)*f(l1+3:l2+3,m,n)
             !!lignore???
          endif
        else
          df=0.
          if (ip<=5) print*, 'deri_3d: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          if ( lequidist(j) ) then
            if ( mod(i,2) == 0 ) then
              df =  coeffsy(0,i,is)* f(l1:l2,m  ,n)                 & 
                  + coeffsy(1,i,is)*(f(l1:l2,m+1,n)+f(l1:l2,m-1,n)) &
                  + coeffsy(2,i,is)*(f(l1:l2,m+2,n)+f(l1:l2,m-2,n)) &
                  + coeffsy(3,i,is)*(f(l1:l2,m+3,n)+f(l1:l2,m-3,n))   
            else
              df =  coeffsy(1,i,is)*(f(l1:l2,m+1,n)-f(l1:l2,m-1,n)) &
                  + coeffsy(2,i,is)*(f(l1:l2,m+2,n)-f(l1:l2,m-2,n)) &
                  + coeffsy(3,i,is)*(f(l1:l2,m+3,n)-f(l1:l2,m-3,n))
            endif
          else 
            mm = m-m1+1
            df =  coeffsy(-3,i,mm)*f(l1:l2,m-3,n) & 
                + coeffsy(-2,i,mm)*f(l1:l2,m-2,n) &
                + coeffsy(-1,i,mm)*f(l1:l2,m-1,n) &
                + coeffsy( 0,i,mm)*f(l1:l2,m  ,n) &
                + coeffsy( 1,i,mm)*f(l1:l2,m+1,n) &
                + coeffsy( 2,i,mm)*f(l1:l2,m+2,n) &
                + coeffsy( 3,i,mm)*f(l1:l2,m+3,n)
             !!lignore???
          endif
          if (lmetric.and.(lspherical_coords .or. lcylindrical_coords)) df = df*r1i(:,i)
        else
          df=0.
          if (ip<=5) print*, 'deri_3d: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          if ( lequidist(j) ) then
            if ( mod(i,2) == 0 ) then
              df =  coeffsz(0,i,is)* f(l1:l2,m,n  )                 &
                  + coeffsz(1,i,is)*(f(l1:l2,m,n+1)+f(l1:l2,m,n-1)) &
                  + coeffsz(2,i,is)*(f(l1:l2,m,n+2)+f(l1:l2,m,n-2)) &
                  + coeffsz(3,i,is)*(f(l1:l2,m,n+3)+f(l1:l2,m,n-3))     
            else
              df =  coeffsz(1,i,is)*(f(l1:l2,m,n+1)-f(l1:l2,m,n-1)) &
                  + coeffsz(2,i,is)*(f(l1:l2,m,n+2)-f(l1:l2,m,n-2)) &
                  + coeffsz(3,i,is)*(f(l1:l2,m,n+3)-f(l1:l2,m,n-3))
            endif
          else
            nn = n-n1+1
            df =  coeffsz(-3,i,nn)*f(l1:l2,m,n-3) & 
                + coeffsz(-2,i,nn)*f(l1:l2,m,n-2) &
                + coeffsz(-1,i,nn)*f(l1:l2,m,n-1) &
                + coeffsz( 0,i,nn)*f(l1:l2,m,n  ) &
                + coeffsz( 1,i,nn)*f(l1:l2,m,n+1) &
                + coeffsz( 2,i,nn)*f(l1:l2,m,n+2) &
                + coeffsz( 3,i,nn)*f(l1:l2,m,n+3)
            !!lignore??
          endif
          if (lmetric.and.lspherical_coords) df=df*r1i(:,i)*sth1i(m,i)
        else
          df=0.
          if (ip<=5) print*, 'deri_3d: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine deri_3d
!***********************************************************************
    subroutine deri_3d_inds(f,df,inds,j,lnometric)
!
! only for first derivatives, i.e., inds(i) = 1,7,8 !!!
!
      real, dimension (mx,my,mz)          :: f
      real, dimension (nx)                :: df
      integer, dimension(nx)              :: inds
      integer                             :: j
      logical,                   optional :: lnometric
!
      intent(in)  :: f,j,inds,lnometric
      intent(out) :: df
!
      integer :: mm, nn, ii, li, ini
      logical :: lmetric

      if (present(lnometric)) then
        lmetric = .not.lnometric
      else
        lmetric = .true.
      endif

      if (j==1) then
        if (nxgrid/=1) then
          do ii=1,nx
            li = l1+ii-1
            ini = inds(ii)
            df =  coeffsx(-3,ini,ii)*f(li-3,m,n) & 
                + coeffsx(-2,ini,ii)*f(li-2,m,n) &
                + coeffsx(-1,ini,ii)*f(li-1,m,n) &
                + coeffsx( 0,ini,ii)*f(li,  m,n) &
                + coeffsx( 1,ini,ii)*f(li+1,m,n) &
                + coeffsx( 2,ini,ii)*f(li+2,m,n) &
                + coeffsx( 3,ini,ii)*f(li+3,m,n)
          enddo
        else
          df=0.
          if (ip<=5) print*, 'deri_3d: Degenerate case in x-direction'
        endif
      elseif (j==2) then
        if (nygrid/=1) then
          mm = m-m1+1
          df =  coeffsy(-3,inds,mm)*f(l1:l2,m-3,n) & 
              + coeffsy(-2,inds,mm)*f(l1:l2,m-2,n) &
              + coeffsy(-1,inds,mm)*f(l1:l2,m-1,n) &
              + coeffsy( 0,inds,mm)*f(l1:l2,m  ,n) &
              + coeffsy( 1,inds,mm)*f(l1:l2,m+1,n) &
              + coeffsy( 2,inds,mm)*f(l1:l2,m+2,n) &
              + coeffsy( 3,inds,mm)*f(l1:l2,m+3,n)
          if (lmetric.and.(lspherical_coords .or. lcylindrical_coords)) df = df*r1i(:,1)
        else
          df=0.
          if (ip<=5) print*, 'deri_3d: Degenerate case in y-direction'
        endif
      elseif (j==3) then
        if (nzgrid/=1) then
          nn = n-n1+1
          df =  coeffsz(-3,inds,nn)*f(l1:l2,m,n-3) & 
              + coeffsz(-2,inds,nn)*f(l1:l2,m,n-2) &
              + coeffsz(-1,inds,nn)*f(l1:l2,m,n-1) &
              + coeffsz( 0,inds,nn)*f(l1:l2,m,n  ) &
              + coeffsz( 1,inds,nn)*f(l1:l2,m,n+1) &
              + coeffsz( 2,inds,nn)*f(l1:l2,m,n+2) &
              + coeffsz( 3,inds,nn)*f(l1:l2,m,n+3)
          if (lmetric.and.lspherical_coords) df=df*r1i(:,1)*sth1i(m,1)
        else
          df=0.
          if (ip<=5) print*, 'deri_3d: Degenerate case in z-direction'
        endif
      endif
!
    endsubroutine deri_3d_inds
!***********************************************************************
    subroutine deri(f,df,i1,i2,i,coefs,lequi)
!
!  i-th derivative operating on a 1-D array
!
!   9-feb-07/axel: adapted from der_main; note that f is not the f array!
!
      real, dimension(*)       , intent(in)  :: f
      real, dimension(-3:3,8,*), intent(in)  :: coefs
      real, dimension(*)       , intent(out) :: df
      integer                  , intent(in)  :: i1,i2,i
      logical                  , intent(in)  :: lequi
!
      integer :: n1

      n1 = i2-i1+1

      if ( n1/=1 ) then
        if ( lequi ) then

          if ( mod(i,2) == 0 ) then
            df(1:n1) =  coefs(0,i,1)* f(i1  :i2  )               &
                      + coefs(1,i,1)*(f(i1+1:i2+1)+f(i1-1:i2-1)) &
                      + coefs(2,i,1)*(f(i1+2:i2+2)+f(i1-2:i2-2)) &
                      + coefs(3,i,1)*(f(i1+3:i2+3)+f(i1-3:i2-3))   
          else
            df(1:n1) =  coefs(1,i,1)*(f(i1+1:i2+1)-f(i1-1:i2-1)) &
                      + coefs(2,i,1)*(f(i1+2:i2+2)-f(i1-2:i2-2)) &
                      + coefs(3,i,1)*(f(i1+3:i2+3)-f(i1-3:i2-3))
          endif

        else
          df(1:n1) =  coefs(-3,i,1:n1)*f(i1-3:i2-3) &
                    + coefs(-2,i,1:n1)*f(i1-2:i2-2) &
                    + coefs(-1,i,1:n1)*f(i1-1:i2-1) & 
                    + coefs( 0,i,1:n1)*f(i1  :i2  ) &
                    + coefs( 1,i,1:n1)*f(i1+1:i2+1) &
                    + coefs( 2,i,1:n1)*f(i1+2:i2+2) &
                    + coefs( 3,i,1:n1)*f(i1+3:i2+3)
        endif
        
      else
        df(1:n1) = 0.
        if (ip<=5) print*, 'deri: Degenerate case!'
      endif
!
    endsubroutine deri
!***********************************************************************
    subroutine deri_2d(f,df,j,i,coefs,lequi)
!
!  i-th derivative operating on a 2-D array
!
!   9-feb-07/axel: adapted from der_main; note that f is not the f array!
!
      real, dimension (:,:)     , intent(in)  :: f
      real, dimension (-3:3,8,*), intent(in)  :: coefs
      real, dimension (*)       , intent(out) :: df
      integer                   , intent(in)  :: i,j
      logical                   , intent(in)  :: lequi
!
      integer :: n1

      n1 = size(f,1)

      if ( n1/=1 ) then

        if ( lequi ) then

          if ( mod(i,2) == 0 ) then
            df(1:n1) =  coefs(0,i,1)* f(:,j  )           &
                      + coefs(1,i,1)*(f(:,j+1)+f(:,j-1)) &
                      + coefs(2,i,1)*(f(:,j+2)+f(:,j-2)) &
                      + coefs(3,i,1)*(f(:,j+3)+f(:,j-3))   
          else
            df(1:n1) =  coefs(1,i,1)*(f(:,j+1)-f(:,j-1)) &
                      + coefs(2,i,1)*(f(:,j+2)-f(:,j-2)) &
                      + coefs(3,i,1)*(f(:,j+3)-f(:,j-3))
          endif

        else
          df(1:n1) =  coefs(-3,i,1:n1)*f(:,j-3) &
                    + coefs(-2,i,1:n1)*f(:,j-2) &
                    + coefs(-1,i,1:n1)*f(:,j-1) & 
                    + coefs( 0,i,1:n1)*f(:,j  ) &
                    + coefs( 1,i,1:n1)*f(:,j+1) &
                    + coefs( 2,i,1:n1)*f(:,j+2) &
                    + coefs( 3,i,1:n1)*f(:,j+3)
        endif
        
      else
        df(1:n1) = 0.
        if (ip<=5) print*, 'deri_2d: Degenerate case!'
      endif
!
    endsubroutine deri_2d
!***********************************************************************
    subroutine der_z(f,df)
!
!  first z derivative operating on a z-dependent 1-D array
!
!   9-feb-07/axel: adapted from der_main; note that f is not the f array!
!
      real, dimension (mz), intent(in)  :: f
      real, dimension (nz), intent(out) :: df

      call deri(f,df,n1,n2,1,coeffsz,lequidist(3))
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
      call deri(f,df2,n1,n2,2,coeffsz,lequidist(3))
!
    endsubroutine der2_z
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

        call deri(pencil,df,l1,l2,1,coeffsx,lequidist(j))

      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) then
          if (lroot) print*, 'der_pencil: pencil must be of size my for y derivative'
          call fatal_error('der_pencil','')
        endif

        call deri(pencil,df,m1,m2,1,coeffsy,lequidist(j))
        if (lspherical_coords.or.lcylindrical_coords) df = df*r1i(:,1)        !tbc

      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) then
          if (lroot) print*, 'der_pencil: pencil must be of size mz for z derivative'
          call fatal_error('der_pencil','')
        endif

        call deri(pencil,df,n1,n2,1,coeffsz,lequidist(j))
        if (lspherical_coords) df=df*r1i(:,1)*sth1i(m,1)   !tbc

      else
        if (lroot) print*, 'der_pencil: no such direction j=', j
        call fatal_error('der_pencil','')
      endif
!
    endsubroutine der_pencil
!***********************************************************************
    subroutine der2_main(f,k,df2,j)
!
!  calculate 2nd derivative d^2f_k/dx_j^2
!  accurate to 6th order, explicit, periodic
!  replace cshifts by explicit construction -> x6.5 faster!
!
!   1-oct-97/axel: coded
!   1-apr-01/axel+wolf: pencil formulation
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx)               :: df2
      integer                            :: j,k
!
      intent(in)  :: f,k,j
      intent(out) :: df2
!
!debug      if (loptimise_ders) der_call_count(k,icount_der2,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der2,j,1) + 1 !DERCOUNT
!
      call deri_3d( f(1,1,1,k), df2, 2, j )
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
      real, dimension (nx)       :: df2
      integer                    :: j
!
      intent(in)  :: f,j
      intent(out) :: df2
!
!debug      if (loptimise_ders) der_call_count(k,icount_der2,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der2,j,1) + 1 !DERCOUNT
!
!
      call deri_3d(f,df2,2,j)
!
    endsubroutine der2_other
!***********************************************************************
    subroutine der2_pencil(j,pencil,df2)
!
!  Calculate 2nd derivative of any x, y or z pencil.
!
!  01-nov-07/anders: adapted from der2
!
! non-cartesian???
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

        call deri(pencil,df2,l1,l2,2,coeffsx,lequidist(j))

      else if (j==2) then
!
!  y-derivative
!
        if (size(pencil)/=my) then
          if (lroot) &
              print*, 'der2_pencil: pencil must be of size my for y derivative'
          call fatal_error('der2_pencil','')
        endif

        call deri(pencil,df2,m1,m2,2,coeffsy,lequidist(j))
        if (lspherical_coords.or.lcylindrical_coords) df2 = df2*r1i(:,2)        !tbc

      else if (j==3) then
!
!  z-derivative
!
        if (size(pencil)/=mz) then
          if (lroot) &
              print*, 'der2_pencil: pencil must be of size mz for z derivative'
          call fatal_error('der2_pencil','')
        endif

        call deri(pencil,df2,n1,n2,2,coeffsz,lequidist(j))
        if (lspherical_coords) df2=df2*r1i(:,2)*sth1i(m,2)   !tbc

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
      real, dimension (nx)               :: df
      integer                            :: j,k
      logical, optional                  :: ignoredx
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der5,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der5,j,1) + 1 !DERCOUNT
!
      if (present(ignoredx)) then
        call deri_3d(f(1,1,1,k),df,3,j,ignoredx)
      else
        call deri_3d(f(1,1,1,k),df,3,j)
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
      real, dimension (nx)               :: df
      integer                            :: j,k
      logical, optional                  :: ignoredx,upwind
!
      intent(in)  :: f,k,j,ignoredx,upwind
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der4,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der4,j,1) + 1 !DERCOUNT
!
      if (present(upwind)) &
        call warning('der4','upwinding not implemented')
!
      if (present(ignoredx)) then
        call deri_3d(f(1,1,1,k),df,4,j,ignoredx)
      else
        call deri_3d(f(1,1,1,k),df,4,j)
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
      real, dimension (nx)               :: df
      integer                            :: j,k
      logical, optional                  :: ignoredx
!
      intent(in)  :: f,k,j,ignoredx
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der5,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der5,j,1) + 1 !DERCOUNT
!
      if (present(ignoredx)) then
        call deri_3d(f(1,1,1,k),df,5,j,ignoredx)
      else
        call deri_3d(f(1,1,1,k),df,5,j)
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
      real, dimension (nx) :: df
      integer :: j,k
      logical, optional :: ignoredx,upwind
!
      intent(in)  :: f,k,j,ignoredx,upwind
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_der6,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der6,j,1) + 1 !DERCOUNT
!
      if (present(ignoredx)) then
        if (present(upwind)) then
          call der6_other(f(1,1,1,k),df,j,ignoredx,upwind)
        else
          call der6_other(f(1,1,1,k),df,j,ignoredx)
        endif
      else
        if (present(upwind)) then
          call der6_other(f(1,1,1,k),df,j,UPWIND=upwind)
        else
          call der6_other(f(1,1,1,k),df,j)
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
            tmp1= f(i,m,n,j)-f(i-1,m,n,j)
            tmp2=(f(i+1,m,n,j)-f(i-1,m,n,j))/4.
            tmp3= f(i+1,m,n,j)-f(i,m,n,j)
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
          if (lspherical_coords.or.lcylindrical_coords) fac = fac*r1i(:,1)
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
          if (lspherical_coords) fac = fac*r1i(:,1)*sth1i(m,i)
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
      real, dimension (mx,my,mz)         :: f
      real, dimension (nx)               :: df
      integer                            :: j
      logical,                  optional :: ignoredx,upwind
!
      intent(in)  :: f,j,ignoredx,upwind
      intent(out) :: df
!
      logical :: igndx,upwnd
!
!debug      if (loptimise_ders) der_call_count(k,icount_der6,j,1) = & !DERCOUNT
!debug                          der_call_count(k,icount_der6,j,1) + 1 !DERCOUNT
!
      if (present(upwind)) then
        upwnd = upwind
        if ( upwnd .and. .not.lequidist(j) ) &
          call fatal_error('der6_other','Not to be used '// &    !tbc
          'for upwinding when grid is non-equidistant')
      else
        upwnd = .false.
        !if (.not.lcartesian_coords) &
        !     call fatal_error('der6_other','in non-cartesian coordinates '// &    !tbc
        !     'just works if upwinding is used')
      endif
!
      if (present(ignoredx)) then
        igndx = ignoredx
        upwnd = .false.   !tbc
      else
        igndx = upwnd
      endif

      call deri_3d(f,df,6,j,igndx,lnometric=upwnd)
!
      if ( upwnd .and. lequidist(j) ) then
        if (j==1) then
          df = df*(dx_1(1)/60.)           ! as dx is ignored for equidistant grids
        elseif (j==2) then
          df = df*(dy_1(1)/60.)
        elseif (j==3) then
          df = df*(dz_1(1)/60.)
        endif
      endif
!
    endsubroutine der6_other
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
      real, dimension (nx) :: df
      integer :: i,j,k
!
      intent(in) :: f,k,i,j
      intent(out) :: df
!
!debug      if (loptimise_ders) der_call_count(k,icount_derij,i,j) = & !DERCOUNT
!debug                          der_call_count(k,icount_derij,i,j) + 1 !DERCOUNT
!
    call derij_other(f(1,1,1,k),df,i,j)
!
    endsubroutine derij_main
!***********************************************************************
    subroutine derij_other(f,df,i,j)
!
!  calculate 2nd derivative with respect to two different directions, i,j
!  input: scalar, output: scalar
!  accurate to 6th order, explicit, periodic
!   8-sep-01/axel: coded
!  25-jun-04/tobi+wolf: adapted for non-equidistant grids
!  14-nov-06/wolf: implemented bidiagonal scheme
!
      real, dimension (mx,my,mz) :: f
      real, dimension (nx)       :: df
      integer                    :: i,j
!
      intent(in)  :: f,i,j
      intent(out) :: df
!
      real, dimension (nx) :: fac
      real, dimension (nx,-3:3) :: work
      integer :: ii
      logical :: lbidiag
!
!debug      if (loptimise_ders) der_call_count(k,icount_derij,i,j) = & !DERCOUNT
!debug                          der_call_count(k,icount_derij,i,j) + 1 !DERCOUNT
!
      if ( i==j ) return     !warning??

      lbidiag = lbidiagonal_derij .and. lequidist(i) .and. lequidist(j) 

      if (lbidiag) then
        !
        ! Use bidiagonal mixed-derivative operator, i.e.
        ! employ only the three neighbouring points on each of the four
        ! half-diagonals. This gives 6th-order mixed derivatives as the
        ! version below, but involves just 12 points instead of 36.
        !
        if ( i+j==3 ) then
          if (nxgrid/=1.and.nygrid/=1) then
            fac=(1./720.)*dx_1(1)*dy_1(1)
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
        elseif ( i+j==5 ) then
          if (nygrid/=1.and.nzgrid/=1) then
            fac=(1./720.)*dy_1(1)*dz_1(1)
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
        elseif ( i+j==4 ) then
          if (nzgrid/=1.and.nxgrid/=1) then
            fac=(1./720.)*dz_1(1)*dx_1(1)
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
      else                      ! not using bidiagonal mixed derivatives, necessary for non-uniform grid
        !
        ! This is the old, straight-forward scheme
        !
        if ( i+j==3 ) then
          if (nxgrid/=1.and.nygrid/=1) then

            do ii=-3,+3
              if ( ii/=0 .or. .not.lequidist(2) ) &
                call deri(f(1,m+ii,n),work(1,ii),l1,l2,1,coeffsx,lequidist(1))  ! x-derivative
            enddo

            call deri_2d(work,df,4,1,coeffsy,lequidist(2))                      ! y-derivative

          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in x- or y-direction'
          endif
        elseif ( i+j==5 ) then
          if (nygrid/=1.and.nzgrid/=1) then

            do ii=-3,+3
              if ( ii/=0 .or. .not.lequidist(3) ) &
                call deri_2d(f(l1:l2,:,n+ii),work(1,ii),m,1,coeffsy,lequidist(2))   ! y-derivative   performance?
            enddo

            call deri_2d(work,df,4,1,coeffsz,lequidist(3))                          ! z-derivative    
          else
            df=0.
            if (ip<=5) print*, 'derij: Degenerate case in y- or z-direction'
          endif
        elseif ( i+j==4 ) then
          if (nzgrid/=1.and.nxgrid/=1) then

            do ii=-3,+3
              if ( ii/=0 .or. .not.lequidist(3) ) &
                call deri(f(1,m,n+ii),work(1,ii),l1,l2,1,coeffsx,lequidist(1))  ! x-derivative
            enddo

            call deri_2d(work,df,4,1,coeffsz,lequidist(3))                      ! z-derivative 
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
        if ( i+j==3 ) df=df*r1i(:,1)                 !(minus extra terms)
        if ( i+j==4 ) df=df*r1i(:,1)*sth1i(m,1)      !(minus extra terms)
        if ( i+j==5 ) df=df*r1i(:,2)*sth1i(m,1)      !(minus extra terms)
      endif
!
      if (lcylindrical_coords) then
        if ( i+j==3 .or. i+j==5 ) df=df*r1i(:,1)
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
      real, dimension (nx) :: df
      integer :: i,j,k
!
      intent(in) :: f,k,i,j
      intent(out) :: df

      real, dimension(nx,-3:3) :: work
      integer :: ii

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
      else
        if ((i==1.and.j==2)) then                       ! 5th deriv in x, 1st in y
          if (nxgrid/=1.and.nygrid/=1) then

              do ii=-3,+3
                if ( ii/=0 .or. .not.lequidist(2) ) &
                  call deri(f(1,m+ii,n,k),work(1,ii),l1,l2,5,coeffsx,lequidist(1))      ! x-derivative
              enddo

              call deri_2d(work,df,4,1,coeffsy,lequidist(2))                            ! y-derivative
              if (lspherical_coords.or.lcylindrical_coords) df = df*r1i(:,1)            !(minus extra terms) tbc

          else
            df=0.
            if (ip<=5) print*, 'der5i1j: Degenerate case in x- or y-direction'
          endif
        elseif ((i==2.and.j==1)) then                      ! 5th deriv in y, 1st in x
          if (nygrid/=1.and.nxgrid/=1) then

              do ii=-3,+3
                if ( ii/=0 .or. .not.lequidist(2) ) &
                  call deri(f(1,m+ii,n,k),work(1,ii),l1,l2,1,coeffsx,lequidist(1))      ! x-derivative
              enddo

              call deri_2d(work,df,4,5,coeffsy,lequidist(2))                            ! y-derivative
              if (lspherical_coords.or.lcylindrical_coords) df = df*r1i(:,5)            !(minus extra terms) tbc

          else
            df=0.
            if (ip<=5) print*, 'der5i1j: Degenerate case in y- or x-direction'
          endif
        elseif ((i==1.and.j==3)) then                       ! 5th deriv in x, 1st in z
          if (nxgrid/=1.and.nzgrid/=1) then

              do ii=-3,+3
                if ( ii/=0 .or. .not.lequidist(3) ) &
                  call deri(f(1,m,n+ii,k),work(1,ii),l1,l2,5,coeffsx,lequidist(1))      ! x-derivative
              enddo

              call deri_2d(work,df,4,1,coeffsz,lequidist(3))                            ! z-derivative 
              if (lspherical_coords) df = df*r1i(:,1)*sth1i(m,1)                        !(minus extra terms) tbc

          else
            df=0.
            if (ip<=5) print*, 'der5i1j: Degenerate case in x- or z-direction'
          endif
        elseif ((i==3.and.j==1)) then                       ! 5th deriv in z, 1st in x
          if (nzgrid/=1.and.nygrid/=1) then

              do ii=-3,+3
                if ( ii/=0 .or. .not.lequidist(3) ) &
                  call deri(f(1,m,n+ii,k),work(1,ii),l1,l2,1,coeffsx,lequidist(1))      ! x-derivative
              enddo

              call deri_2d(work,df,4,5,coeffsz,lequidist(3))                            ! z-derivative 
              if (lspherical_coords) df = df*r1i(:,5)*sth1i(m,5)                        ! (minus extra terms) tbc

          else
            df=0.
            if (ip<=5) print*, 'der5i1j: Degenerate case in z- or x-direction'
          endif
        elseif ((i==2.and.j==3)) then                       ! 5th deriv in y, 1st in z
          if (nygrid/=1.and.nzgrid/=1) then

              do ii=-3,+3
                if ( ii/=0 .or. .not.lequidist(3) ) &
                  call deri_2d(f(l1:l2,:,n+ii,k),work(1,ii),m,5,coeffsy,lequidist(2))   ! y-derivative   performance?
              enddo

              call deri_2d(work,df,4,1,coeffsz,lequidist(3))                            ! z-derivative    
              if (lspherical_coords  ) df = df*r1i(:,6)*sth1i(m,1)                      ! (minus extra terms) tbc
              if (lcylindrical_coords) df = df*r1i(:,5)                                 ! (minus extra terms) tbc

          else
            df=0.
            if (ip<=5) print*, 'der5i1j: Degenerate case in y- or z-direction'
          endif
        elseif ((i==3.and.j==2)) then                       ! 5th deriv in z, 1st in y
          if (nzgrid/=1.and.nygrid/=1) then

              do ii=-3,+3
                if ( ii/=0 .or. .not.lequidist(3) ) &
                  call deri_2d(f(l1:l2,:,n+ii,k),work(1,ii),m,1,coeffsy,lequidist(2))   ! y-derivative   performance?
              enddo

              call deri_2d(work,df,4,5,coeffsz,lequidist(3))                            ! z-derivative    
              if (lspherical_coords  ) df = df*r1i(:,6)*sth1i(m,5)                      ! (minus extra terms) tbc
              if (lcylindrical_coords) df = df*r1i(:,1)                                 ! (minus extra terms) tbc

          else
            df=0.
            if (ip<=5) print*, 'der5i1j: Degenerate case in z- or y-direction'
          endif
        else
          print*, 'der5i1j: no such value for i,j=', i, j
          call fatal_error('der5i1j','')
        endif
!
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
          if (ip<=5) print*, 'der_onesided_4_slice: Degenerate case in x-direction'
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
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_4_slice_other_pt: Degenerate case in z-direction'
        endif
      endif
    endsubroutine der_onesided_4_slice_other_pt
!***********************************************************************
   subroutine calc_coeffs( grid, coeffs )

   real, dimension(-2:3), intent(in)  :: grid
   real, dimension(-3:3,8), intent(out) :: coeffs

   real :: h1, h2, h3, hm1, hm2, hm3, h12, h23, h13, hm12, hm23, hm13, &
           h1m1, h1m12, h1m13, h12m1, h13m1, h12m12, h13m12, h12m13, h13m13

   h1    = grid(1);   h2    = grid(2); h3    = grid(3)
   hm1   = grid(0);   hm2   = grid(-1);hm3   = grid(-2)
   h12   = h1 + h2;   h23   = h2+h3;   h13   = h12+h3
   hm12  = hm1 + hm2; hm23  = hm2+hm3; hm13  = hm12+hm3
   h1m1  = h1 + hm1;  h1m12 = h1+hm12;      h1m13 = h1 + hm13  
   h12m1 = h12 + hm1; h12m12 = h12 + hm12; h12m13 = h12 + hm13  
   h13m1 = h13 + hm1; h13m12 = h13 + hm12; h13m13 = h13 + hm13  

   coeffs(:,1) = (/ -h1*h12*h13*hm1*hm12/(hm3*hm13*hm23*h1m13*h12m13*h13m13), &
                     h1*h12*h13*hm1*hm13/(hm2*hm3*hm12*h1m12*h12m12*h13m12),  &
                    -h1*h12*h13*hm12*hm13/(hm1*hm2*hm23*h1m1*h12m1*h13m1),    &
                     0.,                                                      &
                     h12*h13*hm1*hm12*hm13/(h1*h2*h23*h1m1*h1m12*h1m13),      &
                    -h1*h13*hm1*hm12*hm13/(h2*h3*h12*h12m1*h12m12*h12m13),    &
                     h1*h12*hm1*hm12*hm13/(h3*h13*h23*h13m1*h13m12*h13m13)     /)

   coeffs(0,1) = -sum(coeffs(:,1))

   coeffs(:,2) = 2.*(/ (h12*h13*hm1*hm12 + h1*(h13*hm1*hm12 + h12*(hm1*hm12 - h13*(hm1 + hm12))))/       &
                       (hm3*hm13*hm23*h1m13*h12m13*h13m13),                                              &
!
                      -(h12*h13*hm1*hm13 + h1*(h13*hm1*hm13 + h12*(hm1*hm13 - h13*(hm1 + hm13))))/       &
                       (hm2*hm3*hm12*h1m12*h12m12*h13m12),                                               &
!
                       (h12*h13*hm12*hm13 + h1*(h13*hm12*hm13 + h12*(hm12*hm13 - h13*(hm12 + hm13))))/   &
                       (hm1*hm2*hm23*h1m1*h12m1*h13m1),                                                  &
!
                       0.,                                                                               &
!
                      -(h13*hm1*hm12*hm13 + h12*(hm1*hm12*hm13 - h13*(hm12*hm13 + hm1*(hm12 + hm13))))/  &
                       (h1*h2*h23*h1m1*h1m12*h1m13),                                                     &
!
                       (h13*hm1*hm12*hm13 + h1*(hm1*hm12*hm13 -  h13*(hm12*hm13 + hm1 * (hm12 + hm13))))/&
                       (h2*h3*h12*h12m1*h12m12*h12m13),                                                  &
!
                      -(h12*hm1*hm12*hm13 + h1*(hm1*hm12*hm13 - h12*(hm12*hm13 + hm1*(hm12 + hm13))))/   &
                       (h3*h13*h23*h13m1*h13m12*h13m13)                                                   /)

   coeffs(0,2) = -sum(coeffs(:,2))

   coeffs(:,3) = 6.*(/-(h13*hm1*hm12 + h12*(hm1*hm12 - h13*(hm1 + hm12)) +              &
                        h1*(hm1*hm12 - h13*(hm1 + hm12) + h12*(h13 - hm1 - hm12)))/     &
                       (hm3*hm13*hm23*h1m13*h12m13*h13m13),                             &
!
                       (h13*hm1*hm13 + h12*(hm1*hm13 - h13*(hm1 + hm13)) +              &
                        h1*(hm1*hm13 - h13*(hm1 + hm13) + h12*(h13 - hm1 - hm13)))/     &
                       (hm2*hm3*hm12*h1m12*h12m12*h13m12),                              &
!
                      -(h13*hm12*hm13 + h12*(hm12*hm13 - h13*(hm12 + hm13)) +           &
                        h1*(hm12*hm13 - h13*(hm12 + hm13) + h12*(h13 - hm12 - hm13)))/  &
                       (hm1*hm2*hm23*h1m1*h12m1*h13m1),                                 &
!
                       0.,                                                              &
!
                       (hm1*hm12*hm13 - h13*(hm12*hm13 + hm1*(hm12 + hm13)) -           &
                        h12*(hm12*hm13 + hm1*(hm12 + hm13) - h13*(hm1 + hm12 + hm13)))/ &
                       (h1*h2*h23*h1m1*h1m12*h1m13),                                    &
!
                      -(hm1*hm12*hm13 - h13*(hm12*hm13 + hm1*(hm12 + hm13)) -           &
                        h1*(hm12*hm13 + hm1*(hm12 + hm13) - h13*(hm1 + hm12 + hm13)))/  &
                       (h2*h3*h12*h12m1*h12m12*h12m13),                                 &
!
                       (hm1*hm12*hm13 - h12*(hm12*hm13 + hm1*(hm12 + hm13)) -           &
                        h1*(hm12*hm13 + hm1*(hm12 + hm13) - h12*(hm1 + hm12 + hm13)))/  &
                       (h3*h13*h23*h13m1*h13m12*h13m13)                                  /)

   coeffs(0,3) = -sum(coeffs(:,3))

   coeffs(:,4) = 24.*(/  (-h13*(hm1+hm12) + hm1*hm12  +(h1+h12)*(h13-hm1-hm12) + h1*h12)/ &
                         (hm3*hm13*hm23*h1m13*h12m13*h13m13), &
                        -(-h13*(hm1+hm13) + hm1*hm13  +(h1+h12)*(h13-hm1-hm13) + h1*h12)/ &
                         (hm2*hm3*hm12*h1m12*h12m12*h13m12),  &
                         (-h13*(hm12+hm13)+ hm12*hm13 +(h1+h12)*(h13-hm12-hm13)+ h1*h12)/ &
                         (hm1*hm2*hm23*h1m1*h12m1*h13m1),     &  
                         0.,                                  &
                         (hm1*(hm12+hm13)+ hm12*hm13 -(h12+h13)*(hm1+hm12+hm13)+ h12*h13)/ &
                         (h1*h2*h23*h1m1*h1m12*h1m13),       &
                        -(hm1*(hm12+hm13)+ hm12*hm13 -(h1+h13)*(hm1+hm12+hm13) + h1*h13 )/ &
                         (h2*h3*h12*h12m1*h12m12*h12m13),    &
                         (hm1*(hm12+hm13)+ hm12*hm13 -(h1+h12)*(hm1+hm12+hm13) + h1*h12 )/ &
                         (h3*h13*h23*h13m1*h13m12*h13m13)     /)

   coeffs(0,4) = -sum(coeffs(:,4))

   coeffs(:,5) = 120.*(/-(h1 + h12 + h13 - hm1  - hm12)/(hm3*hm13*hm23*h1m13*h12m13*h13m13), &     
                         (h1 + h12 + h13 - hm1  - hm13)/(hm2*hm3*hm12*h1m12*h12m12*h13m12),  &
                        -(h1 + h12 + h13 - hm12 - hm13)/(hm1*hm2*hm23*h1m1*h12m1*h13m1),     &
                          0.,                                                                &
                        -(h12 + h13 - hm1 - hm12 - hm13)/(h1*h2*h23*h1m1*h1m12*h1m13),       &
                         (h1  + h13 - hm1 - hm12 - hm13)/(h2*h3*h12*h12m1*h12m12*h12m13),    &
                        -(h1  + h12 - hm1 - hm12 - hm13)/(h3*h13*h23*h13m1*h13m12*h13m13)     /)                 

   coeffs(0,5) = -sum(coeffs(:,5))

   coeffs(:,6) = 720.*(/ 1./(hm3*hm13*hm23*h1m13*h12m13*h13m13),&
                        -1./(hm2*hm3*hm12*h1m12*h12m12*h13m12), &
                         1./(hm1*hm2*hm23*h1m1*h12m1*h13m1),    &
                         0.,                                    &
                         1./(h1*h2*h23*h1m1*h1m12*h1m13),       &
                        -1./(h2*h3*h12*h12m1*h12m12*h12m13),    &
                         1./(h3*h13*h23*h13m1*h13m12*h13m13)     /)

   coeffs(0,6) = -sum(coeffs(:,6))

!upwind 1st deriv, 5th order, u>0

   coeffs((/-3,-2,-1,1,2,3/),7) = (h1*h12*hm1*hm12*hm13/720.)*coeffs((/-3,-2,-1,1,2,3/),6) 
   coeffs(0,7) = -sum(coeffs(:,7))

   coeffs(:,7) = (h1*h12*hm1*hm12*hm13/720.)*coeffs(:,6)

!upwind 1st deriv, 5th order, u<0, sign already inverted!   tbc

   coeffs((/-3,-2,-1,1,2,3/),8) = (h1*h12*h13*hm1*hm12/720.)*coeffs((/-3,-2,-1,1,2,3/),6) 
   coeffs(0,8) = -sum(coeffs(:,8))

   coeffs(:,8) = (h1*h12*h13*hm1*hm12/720.)*coeffs(:,6)

   endsubroutine calc_coeffs
!*************************************************************************************************************
   subroutine calc_coeffs_1( grid, coeffs )

   real, dimension(-1:2), intent(in)  :: grid
   real, dimension(-2:2,2), intent(out) :: coeffs

   real :: h1, h2, hm1, hm2, h12, hm12, h1m1, h1m12, h12m1, h12m12

   h1    = grid(1);   h2    = grid(2)  ; hm1  = grid(0);   hm2   = grid(-1)
   h12   = h1 + h2;   hm12  = hm1 + hm2; h1m1 = h1 + hm1;  h1m12 = h1+hm12 
   h12m1 = h12 + hm1; h12m12 = h12 + hm12  

! 2nd order 1st derivative

   coeffs(:,1) = (/ 0., -h1/(hm1*h1m1), 0., hm1/(h1*h1m1), 0. /)
   coeffs(0,1) = -sum(coeffs(:,1))

! 4th order 1st dervative
 
   coeffs(:,2) = (/ h1*h12*hm1/(hm2*hm12*h1m12*h12m12), -h1*h12*hm12/(hm1*hm2*h1m1*h12m1), 0.,  &
                   +h12*hm1*hm12/(h1*h2*h1m1*h1m12),    -h1*hm1*hm12/(h2*h12*h12m1*h12m12)       /)
   coeffs(0,2) = -sum(coeffs(:,2))

   endsubroutine calc_coeffs_1
!*************************************************************************************************************
   real function derivative( func, i, coeffs, ldiff )

   integer                , intent(in)           :: i
   real, dimension(-2:*)  , intent(in)           :: func
   real, dimension(-3:3,*), intent(in)           :: coeffs
   logical,                 intent(in), optional :: ldiff

   real, dimension(-3:3) :: funcl
 
   if ( present(ldiff) ) then
     funcl = func(i-3:i+3) - func(i)
     derivative = sum(funcl*coeffs(:,i))
   else
     derivative = sum(func(i-3:i+3)*coeffs(:,i))
   endif

   endfunction derivative
!*************************************************************************************************************
   real function derivative_1( func, i, coeffs, ldiff )
                          
   integer              , intent(in)           :: i
   real, dimension(-2:*), intent(in)           :: func
   real, dimension(-2:2), intent(in)           :: coeffs
   logical,               intent(in), optional :: ldiff

   real, dimension(-2:2) :: funcl
 
   if ( present(ldiff) ) then
     funcl = func(i-2:i+2) - func(i)
     derivative_1 = sum(funcl*coeffs)
   else
     derivative_1 = sum(func(i-2:i+2)*coeffs)
   endif

   endfunction derivative_1
!*************************************************************************************************************
   subroutine heatflux_boundcond_x( f, inh, fac, im, topbot )

     real, dimension(mx,my,mz,mfarray), intent(INOUT):: f
     real, dimension(my,mz)           , intent(IN)   :: inh
     real                             , intent(IN)   :: fac
     integer                          , intent(IN)   :: im,topbot

     real, dimension(my,mz) :: work_yz
     integer                :: i,ic

     work_yz = inh

     select case (im)
     case (1,2)
       do i=-im,im                                                    ! dlnrho, 1st+2nd order, + Fbot/(K*cs2)
         work_yz = work_yz + coeffsx_1(i,im,topbot)*f(l1+i,:,:,ilnrho)
       enddo
     case (3)
       if (topbot==1) then
         ic=1
       else
         ic=nx
       endif
       do i=-im,im                                                    ! dlnrho, 3rd order, + Fbot/(K*cs2)
         work_yz = work_yz + coeffsx(i,1,ic)*f(l1+i,:,:,ilnrho)
       enddo
     end select

     work_yz = fac*work_yz

     select case (im)
     case (1,2)
       do i=-im+1,im                                                    ! dlnrho, 1st+2nd order, + Fbot/(K*cs2)
         work_yz = work_yz + coeffsx_1(i,im,topbot)*f(l1+i,:,:,iss)
       enddo

       f(l1-im,:,:,iss) = -work_yz/coeffsx_1(-im,im,topbot)

     case (3)
       do i=-im+1,im                                                    ! dlnrho, 3rd order, + Fbot/(K*cs2)
         work_yz = work_yz + coeffsx(i,1,ic)*f(l1+i,:,:,iss)
       enddo

       f(l1-im,:,:,iss) = -work_yz/coeffsx(-im,1,ic)

     end select

   endsubroutine heatflux_boundcond_x
!*************************************************************************************************************
   subroutine test( coord, dc_1, nc, coefs )

   integer :: mc, nc, nc1, nc2
   real, dimension(nc+2*nghost) :: coord, dc_1
   real, dimension(-3:3,8,nc) :: coefs
   real, dimension(-2:2,2,2) :: coefs_1

   real, parameter, dimension(0:6) :: c=(/1.d0,1.1d0,1.2d0,0.d0,0.d0,0.d0,0.d0/)
   real, dimension(nc+2*nghost) :: func

   integer :: i,j
   logical :: t1

   t1=.false.
   goto 1
   entry test_1( coord, dc_1, nc, coefs_1 )

   t1=.true.
   
 1 mc = nc+2*nghost

   do i=1,mc
     func(i) = c(6)
     do j=5,0,-1
       func(i) = func(i)*coord(i) + c(j)
     enddo
   enddo

   nc = mc-2*nghost; nc1 = nghost+1; nc2 = nc1+nc-1

   if (t1) then
     print*, derivative_1( func, 1, coefs_1(-2,1,1) )
     print*, c(1) + 2.e0*c(2)*coord(nc1) + 3.e0*c(3)*coord(nc1)**2 + 4.e0*c(4)*coord(nc1)**3 &
                  + 5.e0*c(5)*coord(nc1)**4 + 6.e0*c(6)*coord(nc1)**5

     print*, derivative_1( func, nc, coefs_1(:,1,2) )
     print*, c(1) + 2.e0*c(2)*coord(nc2) + 3.e0*c(3)*coord(nc2)**2 + 4.e0*c(4)*coord(nc2)**3 &
                  + 5.e0*c(5)*coord(nc2)**4 + 6.e0*c(6)*coord(nc2)**5
     stop
   endif

   !print*, (derivative( func, i, coefs(:,3,:) ), i=1,nc)
   !print*, '==='
   !print*, 6.e0*c(3) + 24.e0*c(4)*coord(nc1:nc2) + 60.e0*c(5)*coord(nc1:nc2)**2 + 120.e0*c(6)*coord(nc1:nc2)**3
   !print*, '==='
   !print*, (derivative( func, i, coefs(:,2,:) ), i=1,nc)
   !print*, '==='
   !print*, 2.e0*c(2) + 6.e0*c(3)*coord(nc1:nc2) + 12.e0*c(4)*coord(nc1:nc2)**2 + 20.e0*c(5)*coord(nc1:nc2)**3 + 30.e0*c(6)*coord(nc1:nc2)**4
   !print*, '==='
   print*, (derivative( func, i, coefs(:,1,:) ), i=1,nc)
   print*, '==='
   print*, c(1) + 2.e0*c(2)*coord(nc1:nc2) + 3.e0*c(3)*coord(nc1:nc2)**2 + 4.e0*c(4)*coord(nc1:nc2)**3 &
                + 5.e0*c(5)*coord(nc1:nc2)**4 + 6.e0*c(6)*coord(nc1:nc2)**5
   print*, '==='
   print*, ( dc_1(i)*(45.0*(func(i+1)-func(i-1))-9.0*(func(i+2)-func(i-2))+func(i+3)-func(i-3))/60., i=nc1,nc2 )


   stop
   endsubroutine test
!*************************************************************************************************************
   endmodule Deriv
