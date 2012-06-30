! $Id$
!  This module provide a way for users to specify custom 
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc. 
!
!  The module provides a set of standard hooks into the Pencil-Code and 
!  currently allows the following customizations:                                        
!
!   Description                                     | Relevant function call 
!  ---------------------------------------------------------------------------
!   Special variable registration                   | register_special 
!     (pre parameter read)                          |
!   Special variable initialization                 | initialize_special 
!     (post parameter read)                         |
!                                                   |
!   Special initial condition                       | init_special
!    this is called last so may be used to modify   |
!    the mvar variables declared by this module     |
!    or optionally modify any of the other f array  |
!    variables.  The latter, however, should be     |
!    avoided where ever possible.                   |
!                                                   |
!   Special term in the mass (density) equation     | special_calc_density
!   Special term in the momentum (hydro) equation   | special_calc_hydro
!   Special term in the entropy equation            | special_calc_entropy
!   Special term in the induction (magnetic)        | special_calc_magnetic 
!      equation                                     |
!                                                   |
!   Special equation                                | dspecial_dt
!     NOT IMPLEMENTED FULLY YET - HOOKS NOT PLACED INTO THE PENCIL-CODE 
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 0
!
!***************************************************************

!-------------------------------------------------------------------
!
! HOW TO USE THIS FILE
! --------------------
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the 
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code CVS repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or 
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module 
! selections to say something like:
!
!    SPECIAL=special/geo_kws
!
! Where geo_kws it replaced by the filename of your new module
! upto and not including the .f90
!
!--------------------------------------------------------------------

module Special

  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages

  implicit none

  include '../special.h'
  
  type :: line_param
    real :: x0          ! x position
    real :: y0          ! y position
    real :: amp         ! amplitude of a disturbance of the vortex line
    real :: ll          ! wavelength of the above disturbance
    real :: sgn         ! sign of the argument of the line
  end type line_param

  type :: ring_param
    real :: x0          ! x position
    real :: y0          ! y position
    real :: r0          ! radius
    real :: dir         ! Propagation direction (+/-1)
  end type ring_param

  integer, parameter :: iboundary_point_SYMMETRIC = 1
  integer, parameter :: iboundary_point_ZERO      = 2
  integer, parameter :: iboundary_point_LINEAR    = 3

  type :: boundary_point
    integer :: ix          ! x position
    integer :: iy          ! y position
    integer :: iz          ! z position
    integer :: pnt_type 
    real :: bdry_value  
    real :: rr    ! Distance from the boundary
    real, dimension(3) :: xxp   ! Image point location
    integer, dimension(3) :: inear  ! Image point nearest grid point
  end type boundary_point

  integer, parameter :: nboundary_points = 3000
 
  type (boundary_point), dimension(nboundary_points) :: boundary_pnts
  integer :: nboundary_pnts = 0

  logical :: ltest_sphere = .false.
  logical :: limag_time = .false.
  real :: diff_boundary = 0.000
  real :: vortex_spacing = 12.000
  real :: ampl = 0.000
  real :: frame_Ux = 0.
  real :: test_sphere_radius = 0.
  character(len=50) :: initgpe = 'constant'

! input parameters
  namelist /gpe_init_pars/ initgpe, vortex_spacing, ampl, &
                          test_sphere_radius

! run parameters
  namelist /gpe_run_pars/ diff_boundary, &
                          limag_time, &
                          frame_Ux, test_sphere_radius

!!
!! Declare any index variables necessary for main or 
!! 
   integer :: ipsi_real=0
   integer :: ipsi_imag=0
!!  
!! other variables (needs to be consistent with reset list below)
!!
   integer :: i_modpsim=0
!!

  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables 
!  which should be know to be able to evaluate
! 
!  6-oct-03/tony: coded
!
      use FArrayManager
! 
! Set any required f-array indexes to the next available slot 
!  
      call farray_register_pde('psi_real',ipsi_real)
      call farray_register_pde('psi_imag',ipsi_imag)
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
      real :: rr, r2, bdry_depth, inner_radius, proximity
      real, dimension(3) :: xxp, dr
      integer :: l
!!
!!  Initialize any module variables which are parameter dependent  
!!
      if (test_sphere_radius > 0.) ltest_sphere=.true.
!
      if (ltest_sphere) then
        bdry_depth=nghost*dxmax
        open(82,file='data/xxp_in.dat')
        open(83,file='data/xxp_out.dat')
        open(84,file='data/xxp_out_near.dat')
        inner_radius=max(test_sphere_radius-bdry_depth,0.)
        print*,'initialize_special (gpe): (bdry_depth, inner_radius, test_sphere_radius) = ', &
           bdry_depth, inner_radius, test_sphere_radius
        do n=n1,n2; do m=m1,m2; do l=l1,l2
          r2 = x(l)**2 + y(m)**2 + z(n)**2
          rr=sqrt(r2)
          if ((rr < test_sphere_radius)) then
            if (rr .lt. inner_radius) cycle

            nboundary_pnts=nboundary_pnts+1

            if (nboundary_pnts > nboundary_points) then
              close(82)
              close(83)
              close(84)
              call fatal_error('initialize_special (gpe)', &
                               'Too many boundary points')
            endif
            boundary_pnts(nboundary_pnts)%ix=l
            boundary_pnts(nboundary_pnts)%iy=m
            boundary_pnts(nboundary_pnts)%iz=n
            xxp=(/ x(l), y(m), z(n) /)
            write(82,'(3e17.8)') xxp
            dr=(xxp/rr)*max(test_sphere_radius-rr,0.)
            proximity=(test_sphere_radius-rr)/sqrt(dx**2+dy**2+dz**2)
            boundary_pnts(nboundary_pnts)%rr=rr

            if (proximity <= 0.5) then
              boundary_pnts(nboundary_pnts)%pnt_type=iboundary_point_ZERO
            elseif (proximity <= 1.) then
              boundary_pnts(nboundary_pnts)%pnt_type=iboundary_point_LINEAR
              boundary_pnts(nboundary_pnts)%xxp=xxp+3.*dr
            else
              boundary_pnts(nboundary_pnts)%pnt_type=iboundary_point_SYMMETRIC
              boundary_pnts(nboundary_pnts)%xxp=xxp+2.*dr
            endif 
            boundary_pnts(nboundary_pnts)%inear = &
                         int(((xxp - (/ x(l1), y(m1), z(n1) /)) / (/ dx, dy, dz/))+0.5)
          
            write(83,'(3e17.8)') boundary_pnts(nboundary_pnts)%xxp
            write(84,'(3e17.8)') (/ x(boundary_pnts(nboundary_pnts)%inear(1)), &
                                    y(boundary_pnts(nboundary_pnts)%inear(2)), &
                                    z(boundary_pnts(nboundary_pnts)%inear(3)) /)
          endif
        enddo; enddo; enddo
        close(82)
        close(83)
        close(84)
        if (lroot) print*,'Found ', nboundary_pnts, ' boundary points'
      endif
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      intent(inout) :: f

      type (line_param), parameter :: vl0 = line_param( 0.0, 0.0,0.,33.0, 1.0)
      type (line_param) :: vl1
     ! type (line_param), parameter :: vl1 = line_param( 0.0, 1.1, 0.1,14.6, 1.0)
     ! type (line_param), parameter :: vl2 = line_param(-3.0, 3.0,-0.0,33.0,-1.0)
      type (line_param) :: vl2
      type (line_param) :: vl3
      type (line_param) :: vl4
     ! type (line_param), parameter :: vl3 = line_param( 0.0,-1.1,-0.1,14.6,-1.0)
     ! type (line_param), parameter :: vl4 = line_param(-3.0,-3.0,-0.0,33.0, 1.0)
      type (ring_param) :: vr1
     ! type (line_param) :: vl1
      vl1 = line_param( 0.0, vortex_spacing*0.5,ampl,33.0, 1.0)
      vl3 = line_param( 0.0,-vortex_spacing*0.5,-ampl,33.0,-1.0)
      vl2 = line_param( -20.0, vortex_spacing*0.5,ampl,33.0, 1.0)
      vl4 = line_param( -20.0,-vortex_spacing*0.5,-ampl,33.0,-1.0)
      vr1 = ring_param(-20.0, 0.0, 15.0, -1.0)

!!
!!  SAMPLE IMPLEMENTATION
!!
      select case (initgpe)
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('constant', '0'); 
          f(:,:,:,ipsi_real) = 1.
          f(:,:,:,ipsi_imag) = 1.
        case ('vortex-line'); 
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ipsi_real:ipsi_imag) = vortex_line(vl0)
          enddo; enddo
        case ('vortex-pair'); 
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ipsi_real:ipsi_imag) = complex_mult(vortex_line(vl1), &
                                               vortex_line(vl3))
          enddo; enddo
        case ('sphere'); 
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ipsi_real) = imaged_sphere(0.,0.,0.,1)
            !f(l1:l2,m,n,ipsi_real) = sphere(0.,0.,0.)
            f(l1:l2,m,n,ipsi_imag) = 0.
          enddo; enddo
        case ('add-vortex-ring'); 
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ipsi_real:ipsi_imag) = &
              f(l1:l2,m,n,ipsi_real:ipsi_imag) &
               * vortex_ring(vr1)
          enddo; enddo
        case ('add-vortex-pair'); 
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ipsi_real:ipsi_imag) = &
              f(l1:l2,m,n,ipsi_real:ipsi_imag) &
               * vortex_line(vl2) * vortex_line(vl4) 
          enddo; enddo
        case ('vortex-ring'); 
          do n=n1,n2; do m=m1,m2
            f(l1:l2,m,n,ipsi_real:ipsi_imag) = vortex_ring(vr1)
          enddo; enddo
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_special: No such value for initgpe: ', trim(initgpe)
          call stop_it("")
      endselect
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
! 
!  All pencils that this special module depends on are specified here.
! 
!  18-07-06/tony: coded
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f       
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!     
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!   06-oct-03/tony: coded
!
      use Mpicomm
      use Sub
      use Deriv
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df     
      real, dimension (nx) :: pimag, preal, diss, psi2
      real, dimension (nx) :: del2real, del2imag
      real, dimension (nx) :: drealdx, dimagdx
      real, dimension (nx) :: boundaries
      real :: a, b, c
      type (pencil_case) :: p

!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) &
        print*,'dspecial_dt: SOLVE dpsi_dt (Gross-Pitaevskii Equation)'
!!      if (headtt) call identify_bcs('ss',iss)
!
    preal = f(l1:l2,m,n,ipsi_real)
    pimag = f(l1:l2,m,n,ipsi_imag)
!
!  calculate dpsi/dx
!
     if (frame_Ux /= 0.) then
        call der(f,ipsi_real,drealdx,1)
        call der(f,ipsi_imag,dimagdx,1)
     endif
!
!  calculate the position 3 mesh points away from the boundaries
!  in the positive x, y, and z directions
!
     a = 0.5*Lxyz(1)-3*dxmax
     b = 0.5*Lxyz(2)-3*dxmax
     c = 0.5*Lxyz(3)-3*dxmax
!
!  calculate mask for damping term
!
     if (diff_boundary /= 0.) then
       diss = diff_boundary *((1.0+tanh(x(l1:l2)-a)*tanh(x(l1:l2)+a)) + &
                   (1.0+tanh(y(m)-b)*tanh(y(m)+b)) + &
                   (1.0+tanh(z(n)-c)*tanh(z(n)+c)))
     endif
!
!  calculate del2(psi)
!
!    call der(f, ipsi_real, dpsi)
    !call deriv_z(in_var, dz)
    call del2(f,ipsi_real,del2real)   
    call del2(f,ipsi_imag,del2imag)   

    psi2 = preal**2 + pimag**2

!    if (ltest_sphere) boundaries = sphere_sharp(0.,0.,0.)

    if (limag_time) then
      df(l1:l2,m,n,ipsi_real) = df(l1:l2,m,n,ipsi_real) + &
         (0.5 * ((del2real + diss * del2imag) &
           + (1. - psi2) * (preal + diss * pimag))) !* boundaries
!
      df(l1:l2,m,n,ipsi_imag) = df(l1:l2,m,n,ipsi_imag) + &
         (0.5 * ((del2imag - diss * del2real) &
           + (psi2 - 1.) * (diss * preal - pimag))) !* boundaries
!
      if (frame_Ux /= 0.) then
        df(l1:l2,m,n,ipsi_real) = df(l1:l2,m,n,ipsi_real) + &
           frame_Ux * dimagdx !* boundaries

        df(l1:l2,m,n,ipsi_imag) = df(l1:l2,m,n,ipsi_imag) - &
           frame_Ux * drealdx !* boundaries
      endif
    else
!
!  dpsi/dt = hbar/(2m) * [ diss*del2(psi) + i*del2(psi) ] + (1-|psi|^2)*psi
!  (but use hbar=m=1)
!
      df(l1:l2,m,n,ipsi_real) = df(l1:l2,m,n,ipsi_real) + &
         (0.5 * ((diss * del2real - del2imag) &
           + (1. - psi2) * (diss * preal - pimag))) !*boundaries

      df(l1:l2,m,n,ipsi_imag) = df(l1:l2,m,n,ipsi_imag) + &
         (0.5 * ((del2real + diss * del2imag) &
           + (1. - psi2) * (preal + diss * pimag))) !* boundaries
!
!  dpsi/dt = ... +Ux*dpsi/dx
!
      if (frame_Ux /= 0.) then
        df(l1:l2,m,n,ipsi_real) = df(l1:l2,m,n,ipsi_real) + &
           frame_Ux * drealdx !* boundaries

        df(l1:l2,m,n,ipsi_imag) = df(l1:l2,m,n,ipsi_imag) + &
           frame_Ux * dimagdx !* boundaries 
      endif
    endif
 
!    rhs = 0.5*(eye+diss) * ( laplacian(in_var) + &
!                    (1.0-abs(in_var(:,jsta:jend,ksta:kend))**2)*&
!                             in_var(:,jsta:jend,ksta:kend) ) + &
 !                    Urhs*dpsidx
                     
!    rhs = eye * ( laplacian(in_var) - &
!                    (abs(in_var(:,jsta:jend,ksta:kend))**2)*&
!                         in_var(:,jsta:jend,ksta:kend) ) + &
!                     Urhs*dpsidx


!!
!! SAMPLE DIAGNOSTIC IMPLEMENTATION
!!
      if (ldiagnos) then
        if (i_modpsim/=0) then
          call sum_mn_name(sqrt(psi2),i_modpsim)
! see also integrate_mn_name
        endif
      endif

! Keep compiler quiet by ensuring every parameter is used
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)

    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=gpe_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=gpe_init_pars,ERR=99)
      endif

99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=gpe_init_pars)

    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
    
      if (present(iostat)) then
        read(unit,NML=gpe_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=gpe_run_pars,ERR=99)
      endif

99    return
endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      write(unit,NML=gpe_run_pars)

    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Sub
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite

!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
        i_modpsim=0
      endif
!!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'modpsim',i_modpsim)
      enddo
!!
!!  write column where which magnetic variable is stored
      if (lwr) then
        write(3,*) 'i_modpsim=',i_modpsim
        write(3,*) 'ipsi_real=',ipsi_real
        write(3,*) 'ipsi_imag=',ipsi_imag
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of gross_pitaevskii variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (slice_data) :: slices
!
      integer :: inamev
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  psi2 - Absolute value of the wave function squared
!
        case ('psi2')
          slices%yz=f(ix_loc,m1:m2,n1:n2,ipsi_real)**2 &
                  + f(ix_loc,m1:m2,n1:n2,ipsi_imag)**2 
          slices%xz=f(l1:l2,iy_loc,n1:n2,ipsi_real)**2 &
                  + f(l1:l2,iy_loc,n1:n2,ipsi_imag)**2 
          slices%xy=f(l1:l2,m1:m2,iz_loc,ipsi_real)**2 &
                  + f(l1:l2,m1:m2,iz_loc,ipsi_imag)**2
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,ipsi_real)**2 &
                  + f(l1:l2,m1:m2,iz2_loc,ipsi_imag)**2 
          slices%ready = .true.
!
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine calc_lspecial_pars(f)
!
!  dummy routine
!
!  15-jan-08/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lspecial_pars
!***********************************************************************
    subroutine special_calc_density(df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
!!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(df,p)
!
!   calculate a additional 'special' term on the right hand side of the 
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p

!!
!!  SAMPLE IMPLEMENTATION
!!     (remember one must ALWAYS add to df)
!!  
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
!!
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
!***********************************************************************
    function vortex_line(vl)
! Vortex line initial condition
      real, dimension(nx,2) :: vortex_line
      type (line_param), intent(in) :: vl
      real, dimension(nx) :: vort_r, vort_theta
!  
      call get_r(vl%x0, vl%y0, vl%amp, vl%ll, vort_r)
      call get_theta(vl%x0, vl%y0, vl%amp, vl%ll, vl%sgn, vort_theta)
!  
      vortex_line(:,1) = amp(vort_r) * cos(vort_theta)
      vortex_line(:,2) = amp(vort_r) * sin(vort_theta)
!      
    endfunction vortex_line
!***********************************************************************
    function vortex_ring(vr)
      ! Vortex ring initial condition

      real, dimension(nx,2) :: vortex_ring
      type (ring_param), intent(in) :: vr
      real :: s
      real,parameter :: scal=1.
      real, dimension(nx) :: rr1, rr2, d1, d2
      integer :: i, j, k
  
      call get_s(s, vr%y0)
      
      d1 = sqrt( (scal*(x(l1:l2)-vr%x0))**2 + (s+vr%r0)**2 )
      d2 = sqrt( (scal*(x(l1:l2)-vr%x0))**2 + (s-vr%r0)**2 )
      
      call get_rr(d1,rr1)
      call get_rr(d2,rr2)
      
      rr1 = sqrt( ((0.3437+0.0286*d1**2)) / &
                          (1.0+(0.3333*d1**2)+(0.0286*d1**4)) )
      rr2 = sqrt( ((0.3437+0.0286*d2**2)) / &
                          (1.0+(0.3333*d2**2)+(0.0286*d2**4)) )
  
      vortex_ring(:,1) = rr1*rr2*scal**2*((x(l1:l2)-vr%x0)**2 + &
                          (vr%dir**2*(s-vr%r0)*(s+vr%r0))) 
      vortex_ring(:,2) = rr1*rr2*scal**2*((x(l1:l2)-vr%x0) * &
                          (vr%dir*2.*vr%r0))
  
    endfunction vortex_ring
!***********************************************************************
!    function vortex_ring2(x0, y0, r0, dir)
!      ! Vortex ring initial condition
!      use parameters
!      implicit none
!  
!      complex, dimension(0:nx1,jsta:jend,ksta:kend) :: vortex_ring2
!      real,    intent(in)                           :: x0, y0, r0, dir
!      real,    dimension(0:nx1,ksta:kend)           :: s
!      real,    dimension(0:nx1,jsta:jend,ksta:kend) :: rr1, rr2, d1, d2
!      integer                                       :: i, j, k
!  
!      do k=ksta,kend
!        do i=0,nx1
!          s(i,k) = sqrt((x(i)-x0)**2 + z(k)**2)
!        end do
!      end do
!      
!      do k=ksta,kend
!        do j=jsta,jend
!          do i=0,nx1
!            d1(i,j,k) = sqrt( (y(j)-y0)**2 + (s(i,k)+r0)**2 )
!            d2(i,j,k) = sqrt( (y(j)-y0)**2 + (s(i,k)-r0)**2 )
!          end do
!        end do
!      end do
!      
!      rr1 = sqrt( ((0.3437+0.0286*d1**2)) / &
!                          (1.0+(0.3333*d1**2)+(0.0286*d1**4)) )
!      rr2 = sqrt( ((0.3437+0.0286*d2**2)) / &
!                          (1.0+(0.3333*d2**2)+(0.0286*d2**4)) )
!  
!      do k=ksta,kend
!        do j=jsta,jend
!          do i=0,nx1
!            vortex_ring2(i,j,k) = rr1(i,j,k)*((y(j)-y0)+dir*eye*(s(i,k)+r0)) * &
!                                  rr2(i,j,k)*((y(j)-y0)-dir*eye*(s(i,k)-r0))
!          end do
!        end do
!      end do
!  
!      return
!    endfunction vortex_ring2
!***********************************************************************
    subroutine get_r(vort_x0, vort_y0, vort_a, vort_ll, vort_r)
! Get the cylindrical-polar radius r**2=x**2+y**2
      real, intent(in)  :: vort_x0, vort_y0, vort_a, vort_ll
      real, dimension(nx), intent(out) :: vort_r
!  
       vort_r = sqrt((x(l1:l2)-vort_x0)**2 +  &
            spread((y(m)-vort_y0-vort_a*cos(2.0*pi*z(n)/vort_ll))**2,1,nx))
!  
    endsubroutine get_r
!***********************************************************************
    subroutine get_s(s, sy0)
! Another radial variable
      real, intent(in)  :: sy0
      real, intent(out) :: s
!  
      s = sqrt((y(m)-sy0)**2 + z(n)**2)
!
    endsubroutine get_s
!***********************************************************************
    subroutine get_theta(vort_x0, vort_y0, vort_a, vort_ll, vort_sgn, vort_theta)
! Get the argument theta=arctan(y/x)
      real, intent(in) :: vort_x0, vort_y0, vort_a, vort_ll, vort_sgn
      real, dimension(nx), intent(out) :: vort_theta
!  
      vort_theta = vort_sgn * atan2( &
            spread(y(m)-vort_y0-vort_a*cos(2.0*pi*z(n)/vort_ll),1,nx), &
                         x(l1:l2)-vort_x0)
!  
    endsubroutine get_theta
!***********************************************************************
    subroutine get_rr(r,rr)
! R in psi=R(r)exp(i*theta)
      real, dimension(nx), intent(in)  :: r
      real, dimension(nx), intent(out) :: rr
!      
      rr = sqrt( ((0.3437+0.0286*r**2)) / &
                  (1.0+(0.3333*r**2)+(0.0286*r**4)) )
!  
    endsubroutine get_rr
!***********************************************************************
  function complex_mult(a,b)
    ! Amplitude of a vortex line

    real, dimension(nx,2), intent(in) :: a, b
    real, dimension(nx,2) :: complex_mult

    complex_mult(:,1) = a(:,1)*b(:,1)-a(:,2)*b(:,2)
    complex_mult(:,2) = a(:,2)*b(:,1)+a(:,1)*b(:,2)

  endfunction complex_mult
!***********************************************************************
  function amp(vort_r)
    ! Amplitude of a vortex line

    real, dimension(nx) :: amp
    real, dimension(nx), intent(in) :: vort_r
    real, parameter :: c1 = -0.7
    real, parameter :: c2 = 1.15

    amp = 1.0 - exp(c1*vort_r**c2)

  endfunction amp
!***********************************************************************
  function imaged_sphere(sx,sy,sz,level)
    implicit none

    real, intent(in) :: sx,sy,sz
    real, dimension(nx) :: imaged_sphere, temp
    real :: local_sx, local_sy, local_sz
    integer :: level
    integer :: i,j,k

    temp=1.
    do k=-level,level
      local_sz=sz-k*Lxyz(3) 
      do j=-level,level
        local_sy=sy-j*Lxyz(2) 
        do i=-level,level
          local_sx=sx-i*Lxyz(1) 
          temp=temp*sphere(local_sx,local_sy,local_sz)
        enddo
      enddo
    enddo

    imaged_sphere=temp
    !sphere(i,j,k) = max(0.5*(1.0+&
    !                         tanh(sqrt(x(i)**2+y(j)**2+z(k)**2)-rad)-&
    !                         eps),0.0)

  endfunction imaged_sphere
!***********************************************************************
  function sphere(sx,sy,sz)
    implicit none

    real, intent(in) :: sx,sy,sz
    real, dimension(nx) :: sphere
    real, parameter :: eps = 2.0

!    sphere = 0.5*(1.0+tanh(sqrt((x(l1:l2)-sx)**2+(y(m)-sy)**2+(z(n)-sz)**2)-test_sphere_radius-eps))
    sphere = tanh(sqrt((x(l1:l2)-sx)**2+(y(m)-sy)**2+(z(n)-sz)**2)-test_sphere_radius)

  endfunction sphere
!***********************************************************************
  function sphere_sharp(sx,sy,sz)
    implicit none

    real, intent(in) :: sx,sy,sz
    real, dimension(nx) :: sphere_sharp
    real, parameter :: rad = 10.0
    real, parameter :: eps = 2.0

    sphere_sharp = max(tanh(sqrt((x(l1:l2)-sx)**2+(y(m)-sy)**2+(z(n)-sz)**2)-rad-eps),0.)

  endfunction sphere_sharp
!***********************************************************************
    subroutine gpe_interpolate_quadratic_spline(f,ivar1,ivar2,xxp,gp,inear)
!
!  Quadratic spline interpolation of the function g to the point xxp=(xp,yp,zp).
!
!  10-jun-06/anders: coded
!
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (3) :: xxp
      integer, dimension (3) :: inear
      integer :: ivar1, ivar2
      real, dimension (ivar2-ivar1+1) :: gp
!
      real :: fac_x_m1, fac_x_00, fac_x_p1
      real :: fac_y_m1, fac_y_00, fac_y_p1
      real :: fac_z_m1, fac_z_00, fac_z_p1
      real :: dxp0, dyp0, dzp0
      integer :: i, ix0, iy0, iz0
      logical :: lfirstcall=.true.
!
      intent(in)  :: f, xxp, ivar1
      intent(out) :: gp
!
!  Redefine the interpolation point in coordinates relative to nearest grid
!  point and normalize with the cell size.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      dxp0=(xxp(1)-x(ix0))*dx_1(ix0)
      dyp0=(xxp(2)-y(iy0))*dy_1(iy0)
      dzp0=(xxp(3)-z(iz0))*dz_1(iz0)
!
!  Interpolation formulae.
!
      if (dimensionality==0) then
        gp=f(ix0,iy0,iz0,ivar1:ivar2)
      elseif (dimensionality==1) then
        if (nxgrid/=1) then
          gp = 0.5*(0.5-dxp0)**2*f(ix0-1,iy0,iz0,ivar1:ivar2) + &
                  (0.75-dxp0**2)*f(ix0  ,iy0,iz0,ivar1:ivar2) + &
               0.5*(0.5+dxp0)**2*f(ix0+1,iy0,iz0,ivar1:ivar2)
        endif
        if (nygrid/=1) then
          gp = 0.5*(0.5-dyp0)**2*f(ix0,iy0-1,iz0,ivar1:ivar2) + &
                  (0.75-dyp0**2)*f(ix0,iy0  ,iz0,ivar1:ivar2) + &
               0.5*(0.5+dyp0)**2*f(ix0,iy0+1,iz0,ivar1:ivar2)
        endif
        if (nzgrid/=1) then
          gp = 0.5*(0.5-dzp0)**2*f(ix0,iy0,iz0-1,ivar1:ivar2) + &
                  (0.75-dzp0**2)*f(ix0,iy0,iz0  ,ivar1:ivar2) + &
               0.5*(0.5+dzp0)**2*f(ix0,iy0,iz0+1,ivar1:ivar2)
        endif
      elseif (dimensionality==2) then
        if (nxgrid==1) then
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_y_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_y_00*( f(ix0,iy0  ,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0  ,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_z_00*( f(ix0,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
              fac_y_p1*( f(ix0,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_y_m1*( f(ix0,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 )
        elseif (nygrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_x_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_x_00*( f(ix0  ,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0  ,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_z_00*( f(ix0+1,iy0,iz0  ,ivar1:ivar2)*fac_x_p1 + &
                         f(ix0-1,iy0,iz0  ,ivar1:ivar2)*fac_x_m1 ) + &
              fac_x_p1*( f(ix0+1,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0+1,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_x_m1*( f(ix0-1,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0-1,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 )
        elseif (nzgrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
!
          gp= fac_x_00*fac_y_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_x_00*( f(ix0  ,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0  ,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 ) + &
              fac_y_00*( f(ix0+1,iy0  ,iz0,ivar1:ivar2)*fac_x_p1 + &
                         f(ix0-1,iy0  ,iz0,ivar1:ivar2)*fac_x_m1 ) + &
              fac_x_p1*( f(ix0+1,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0+1,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 ) + &
              fac_x_m1*( f(ix0-1,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0-1,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 )
        endif
      elseif (dimensionality==3) then
        fac_x_m1 = 0.5*(0.5-dxp0)**2
        fac_x_00 = 0.75-dxp0**2
        fac_x_p1 = 0.5*(0.5+dxp0)**2
        fac_y_m1 = 0.5*(0.5-dyp0)**2
        fac_y_00 = 0.75-dyp0**2
        fac_y_p1 = 0.5*(0.5+dyp0)**2
        fac_z_m1 = 0.5*(0.5-dzp0)**2
        fac_z_00 = 0.75-dzp0**2
        fac_z_p1 = 0.5*(0.5+dzp0)**2
!
        gp= fac_x_00*fac_y_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
            fac_x_00*fac_y_00*( f(ix0  ,iy0  ,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0  ,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_z_00*( f(ix0  ,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0  ,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
            fac_y_00*fac_z_00*( f(ix0+1,iy0  ,iz0  ,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0  ,ivar1:ivar2)*fac_x_m1 ) + &
            fac_x_p1*fac_y_p1*( f(ix0+1,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0+1,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_p1*fac_y_m1*( f(ix0+1,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0+1,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_m1*fac_y_p1*( f(ix0-1,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0-1,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_m1*fac_y_m1*( f(ix0-1,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0-1,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_y_p1*( f(ix0  ,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_y_m1*( f(ix0  ,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_y_00*fac_z_p1*( f(ix0+1,iy0  ,iz0+1,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0+1,ivar1:ivar2)*fac_x_m1 ) + &
            fac_y_00*fac_z_m1*( f(ix0+1,iy0  ,iz0-1,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0-1,ivar1:ivar2)*fac_x_m1 ) + &
            fac_z_00*fac_x_p1*( f(ix0+1,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0+1,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
            fac_z_00*fac_x_m1*( f(ix0-1,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0-1,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 )
      endif
!
    endsubroutine gpe_interpolate_quadratic_spline
!***********************************************************************
    subroutine special_before_boundary(f)
!
!   Possibility to modify the f array before the boundaries are 
!   communicated.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
      real, dimension (2)  :: bdry_value
      integer :: i

      if (ltest_sphere) then
        do i=1,nboundary_pnts
          if (boundary_pnts(i)%pnt_type /= iboundary_point_ZERO) then
              call gpe_interpolate_quadratic_spline(f, &
                                 ipsi_real, ipsi_imag, &
                                 boundary_pnts(i)%xxp, &
                                 boundary_pnts(i)%bdry_value, &
                                 boundary_pnts(i)%inear)
          endif
        enddo  
        do i=1,nboundary_pnts
          select case (boundary_pnts(i)%pnt_type)
            case (iboundary_point_SYMMETRIC)
              f( boundary_pnts(i)%ix, &
                 boundary_pnts(i)%iy, &
                 boundary_pnts(i)%iz, ipsi_real:ipsi_imag ) = - boundary_pnts(i)%bdry_value
            case (iboundary_point_LINEAR)
              f( boundary_pnts(i)%ix, &
                 boundary_pnts(i)%iy, &
                 boundary_pnts(i)%iz, ipsi_real:ipsi_imag ) = - 0.5 * boundary_pnts(i)%bdry_value
            case (iboundary_point_ZERO)
              f( boundary_pnts(i)%ix, &
                 boundary_pnts(i)%iy, &
                 boundary_pnts(i)%iz, ipsi_real:ipsi_imag ) = 0.
          endselect
        enddo  
      endif 

    endsubroutine special_before_boundary
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special

