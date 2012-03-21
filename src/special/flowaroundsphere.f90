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
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special

  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include '../special.h'

!
! Sphere geometry
!
  real :: sphere_radius = 0.0,cylinder_temp=293.0
  real :: sph_dia, sph_rad
  integer :: sph_center_x
  integer :: sph_nx, sph_ny
  integer :: sph_l1, sph_l2, sph_n1, sph_n2
!
! Grid geometry, northern and southern hemisphere
!
  logical :: northern = .true., southern = .true.
  integer :: nequator = n1 + nz/2
!
! Div.
!
  character(len=24) :: initspecial='nothing', sphere_type='none',&
      special_inituu='nothing'
  real :: special_infuu=-1., skin_depth=0.0
  logical :: ldisturbance = .false.
  logical :: linlet_northern = .true.
  logical :: ltest = .false.
!
! Input parameters
!
  namelist /flowaroundsphere_init_pars/ &
          initspecial,special_inituu,special_infuu,&
          sphere_radius,sphere_type,&
          sph_l1,sph_l2,sph_n1,sph_n2,&
          sph_center_x,sph_nx,sph_ny,&
          sph_dia,sph_rad,&
          northern,southern,nequator,ltest,&
          skin_depth,cylinder_temp
!
! Run parameters
!
  namelist /flowaroundsphere_run_pars/ &
          special_infuu,&
          sphere_radius,sphere_type,&
          sph_l1,sph_l2,sph_n1,sph_n2,&
          sph_center_x,sph_nx,sph_ny,&
          sph_dia,sph_rad,&
          northern,southern,nequator,&
          ldisturbance,linlet_northern,ltest
!
! Other variables (needs to be consistent with reset list below)
!

  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!
!  6-oct-03/tony: coded
!
      use Cdata
      use FArrayManager, only: farray_register_auxiliary
!
! Write auxiliary variables to var.dat
!
!      lwrite_aux = .true.
!
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      real phi1,phi2
      integer i
!
!  Set up variables for borders between northern and southern hemisphere.
!  The equator is at z=pi, northern hemisphere 0<z<pi, southern pi<z<2pi.
!  The variables will be set up s.t. northern = T if any part of the
!  domain of this processor is on the northern hemisphere, and likewise for
!  southern. n1:nequator-1 will be the part on the northern hemisphere, and
!  nequator:n2 will be the southern part.
!
      phi1=z(n1)
      phi2=z(n2)
      northern = (phi1 < pi)
      southern = (phi2 > pi)
      if (northern) then
        if (southern) then
          nequator = binary_search(pi,z,n1,n2)
        else
          nequator = n2+1
        endif
      else
        nequator = n1
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
      use Cdata
      use Mpicomm
      use Sub
      use Initcond,    only:set_thermodynamical_quantities,gaussian3d
      use EquationOfState, only: cs20
      use FArrayManager, only: farray_use_global

      real, dimension (mx,my,mz,mfarray) :: f

      intent(inout) :: f
      
      integer, pointer :: iglobal_cs2,iglobal_glnTT
      real :: a2,rr2,pphi,wall_smoothing,rr2_low,rr2_high,shiftx,shifty
      real :: wall_smoothing_temp
      integer i,j,k,cyl

      select case (initspecial)
        case ('nothing')
          if (lroot) print*,'init_special: nothing'
          sphere_type = 'none'
          sph_dia = 2.*sphere_radius
        case ('sphere')
          print*,'init_special: sphere'
          sph_center_x = l1
          sph_rad = sphere_radius
          sph_dia = 2.*sphere_radius
          sph_l1 = sph_center_x
          sph_l2 = binary_search(sph_rad,x,l1,l2)
          sph_nx = sph_l2-sph_l1
          sph_n1 = n1
          sph_n2 = n2

          print*,'sphere center: x=',x(sph_center_x)
          print*,'spere edge: x=',x(sph_l2)
          print*,'grid resolution at x=',x(sph_l1),':  dx_1=',dx_1(sph_l1)
          print*,'grid resolution at x=',x(sph_l2),':  dx_1=',dx_1(sph_l2)
          print*,'grid resolution at x=',x(10*sph_l2),':  dx_1=',dx_1(10*sph_l2)
          print*,'grid resolution at x=',x(l2),':  dx_1=',dx_1(l2)
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_special: No such value for initspecial: ',&
              trim(initspecial)
          call stop_it("")
      endselect

      select case (special_inituu)
!
!   This overrides any initial conditions set in the Hydro module.
!
        case ('nothing')
          if (lroot) print*,'special_inituu: nothing'
        case ('cylinderstream')
!   Stream functions for flow around a cylinder as initial condition. 
          a2 = sph_rad**2
          do i=l1,l2
            if (x(i) < sph_rad) cycle
            rr2 = x(i)**2
            do j=n1,n2
              pphi = z(j)
              f(i,m1:m2,j,iux) = special_infuu*sin(pphi)&
                                      *(1. - a2/rr2)
              f(i,m1:m2,j,iuz) = special_infuu*cos(pphi)&
                                      *(1. + a2/rr2)
            end do
          end do
  case ('cylinderstream_nils_x')
!   Stream functions for flow around a cylinder as initial condition. 
          a2 = sph_rad**2
          f(:,:,:,iux:iuz)=0
          shiftx=0
          do i=l1,l2
            do j=m1,m2
              do k=n1,n2
                rr2 = x(i)**2+y(j)**2
                if (rr2 > a2) then
                  do cyl=0,100
                    if (cyl==0) then
                      wall_smoothing=1-exp(-(rr2-a2)/skin_depth**2)
                      f(i,j,k,iuy) = -special_infuu*&
                           2*x(i)*y(j)*a2/rr2**2*wall_smoothing
                      f(i,j,k,iux) = special_infuu*&
                           (1. - a2/rr2 + 2*y(j)**2*a2/rr2**2)&
                           *wall_smoothing
                      if (ilnTT .ne. 0) then
                        wall_smoothing_temp=1-exp(-(rr2-a2)/(sqrt(a2))**2)
                        f(i,j,k,ilnTT) = wall_smoothing_temp*f(i,j,k,ilnTT)&
                             +cylinder_temp*(1-wall_smoothing_temp)
                        f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)*f(l2,m2,n2,ilnTT)/f(i,j,k,ilnTT)
                      endif
                    else
                      shifty=cyl*Lxyz(2)
                      rr2_low =(x(i)+shiftx)**2+(y(j)+shifty)**2
                      rr2_high=(x(i)-shiftx)**2+(y(j)-shifty)**2
                      f(i,j,k,iux) = f(i,j,k,iux)+special_infuu*( &
                           +2*(y(j)-shifty)**2*a2/rr2_high**2-a2/rr2_high&
                           +2*(y(j)+shifty)**2*a2/rr2_low**2 -a2/rr2_low)
                      f(i,j,k,iuy) = f(i,j,k,iuy)-special_infuu*( &
                           +2*(x(i)-shiftx)*(y(j)-shifty)&
                           *a2/rr2_high**2&
                           +2*(x(i)+shiftx)*(y(j)+shifty)&
                           *a2/rr2_low**2)
                    endif
                  enddo
                else
                  if (ilnTT .ne. 0) then
                    f(i,j,k,ilnTT) = cylinder_temp
                    f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)*f(l2,m2,n2,ilnTT)/cylinder_temp
                  endif
                end if
              end do
            end do
          end do
  case ('cylinderstream_nils_y')
!   Stream functions for flow around a cylinder as initial condition. 
          a2 = sph_rad**2
          f(:,:,:,iux:iuz)=0
          shifty=0
          do i=l1,l2
            do j=m1,m2
              do k=n1,n2
                rr2 = x(i)**2+y(j)**2
                if (rr2 > a2) then
                  do cyl=0,100
                    if (cyl==0) then
                      wall_smoothing=1-exp(-(rr2-a2)/skin_depth**2)
                      f(i,j,k,iux) = -special_infuu*&
                           2*x(i)*y(j)*a2/rr2**2*wall_smoothing
                      f(i,j,k,iuy) = special_infuu*&
                           (1. - a2/rr2 + 2*x(i)**2*a2/rr2**2)&
                           *wall_smoothing
                      if (ilnTT .ne. 0) then
                        wall_smoothing_temp=1-exp(-(rr2-a2)/(sqrt(a2))**2)
                        f(i,j,k,ilnTT) = wall_smoothing_temp*f(i,j,k,ilnTT)&
                             +cylinder_temp*(1-wall_smoothing_temp)
                        f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)*f(l2,m2,n2,ilnTT)/f(i,j,k,ilnTT)
                      endif

                    else
                      shiftx=cyl*Lxyz(1)
                      rr2_low =(x(i)+shiftx)**2+(y(j)+shifty)**2
                      rr2_high=(x(i)-shiftx)**2+(y(j)-shifty)**2
                      f(i,j,k,iuy) = f(i,j,k,iuy)+special_infuu*( &
                           +2*(x(i)-shiftx)**2*a2/rr2_high**2-a2/rr2_high&
                           +2*(x(i)+shiftx)**2*a2/rr2_low**2 -a2/rr2_low)
                      f(i,j,k,iux) = f(i,j,k,iux)-special_infuu*( &
                           +2*(x(i)-shiftx)*(y(j)-shifty)&
                           *a2/rr2_high**2&
                           +2*(x(i)+shiftx)*(y(j)+shifty)&
                           *a2/rr2_low**2)
                    endif
                  enddo  
                else
                  if (ilnTT .ne. 0) then
                    f(i,j,k,ilnTT) = cylinder_temp
                    f(i,j,k,ilnrho)=f(l2,m2,n2,ilnrho)*f(l2,m2,n2,ilnTT)/cylinder_temp
                  endif
                end if
              end do
            end do
          end do
        case ('const-x')
!
!   ux = const., uy = 0 in domain.
!
          do i=l1,l2
            do j=1,my
              f(i,j,n1:n2,iux) = &
                  cos(z(n1:n2))*special_infuu
              f(i,j,n1:n2,iuz) = &
                  -sin(z(n1:n2))*special_infuu
            end do
          enddo
        case ('const-y')
!
!   ux = 0, uy = const. in domain.
!
          do i=l1,l2
            do j=1,my
              f(i,j,n1:n2,iux) = &
                  sin(z(n1:n2))*special_infuu
              f(i,j,n1:n2,iuz) = &
                  cos(z(n1:n2))*special_infuu
            end do
          enddo
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_special: No such value for special_inituu: ', &
              trim(special_inituu)
          call stop_it("")
      endselect

      select case (sphere_type)
        case ('none')
          if (lroot) print*,'sphere_type: none'
        case ('freeze','nscbc_wall')
          call sph_vel(f,0.)
        case ('ghosts')
          call sph_vel(f,0.)
          call sph_ghosts(f)
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_special: No such value for sphere_type: ', &
              trim(sphere_type)
          call stop_it("")
      endselect
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
!   Manipulate Hydro pencils.
!   Most basic pencils should come first, as others may depend on them.
!
      real, dimension (mx,my,mz,mfarray) :: f
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
      use Cdata
      use Mpicomm
      use Sub

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p

      intent(inout) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!!      if (headtt) call identify_bcs('ss',iss)
!

    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=flowaroundsphere_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=flowaroundsphere_init_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=flowaroundsphere_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=flowaroundsphere_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=flowaroundsphere_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=flowaroundsphere_run_pars)
!
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
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
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
    subroutine special_calc_density(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
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
      call keep_compiler_quiet(f,df)
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
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
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!!***********************************************************************
    subroutine special_calc_entropy(f,df,p)
!
!   calculate a additional 'special' term on the right hand side of the
!   entropy equation.
!
!   Some precalculated pencils of data are passed in for efficiency
!   others may be calculated directly from the f array
!
!   06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
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
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_entropy
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_boundtreat(f,df)
!
!   Special boundary treatment.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      integer i,j

      intent(inout) :: f,df

      if (sphere_type == 'freeze') then
!   Set df to zero inside sphere
        df(sph_l1:sph_l2-1,:,sph_n1:sph_n2,:) = 0.
      end if

      if (sphere_type == 'ghosts') then
!   Set up ghosts points: Antisymmetric for velocities, symmetric for density.
        call sph_ghosts(f)
!   Set df to zero inside sphere.
        df(sph_l1:sph_l2-1,:,sph_n1:sph_n2,:) = 0.
      endif

      if (sphere_type == 'nscbc_wall') then
!   Do characteristic wall treatment.
        call bc_nscbc_wall(f,df,sph_l2)
!   Set up ghost points near wall inside sphere so that grid points near wall
!   outside sphere gets 1st r-derivative as if calculated by semi-onesided stencils.
!
!   NB: The 'one-sided' condition has shown to be unstable,
!   and make the simulation crash, so I am not sure if we should use it...
!   This is also true for the domain boundaries (using the '1s' condition).
!   Must have equidistant grid?
        !call sph_ghosts_1s(f)
!   Alternatively: Antisymmetric for velocities, symmetric for density.
        call sph_ghosts(f)
!   Set df to zero inside sphere.
        df(sph_l1:sph_l2-1,:,sph_n1:sph_n2,:) = 0.
      endif
!
!   Non-reflecting boundaries
!
      !call bc_nscbc_prf_r(f,df,l2-2)
      !call bc_nscbc_prf_r(f,df,l2-1)
      call bc_nscbc_prf_r(f,df,l2)
!
!   Make a small disturbance to check unstability.
!
      if (ldisturbance .and. northern .and. t > 0.200 .and. t < 0.201) then
        do i=1,my
          do j=0,nghost
              f(l2+j,i,n1:(nequator-1)/2,iux) = &
                  sin(z(n1:(nequator-1)/2))*special_infuu*1.01
              f(l2+j,i,n1:(nequator-1)/2,iuz) = &
                  cos(z(n1:(nequator-1)/2))*special_infuu*1.01
         enddo
        enddo
      endif
    endsubroutine
!***********************************************************************
!
!   Misc. functions and subroutines
!
!***********************************************************************

!***********************************************************************
    subroutine bc_nscbc_prf_r(f,df,lll)
!
!   Set du_r, du_phi and dlnrho at a partially reflecting outlet normal
!   to r-direction acc. to characteristic boundary treatment.
!   Uses a one-sided finite diff. stencil.
!
!   7-jul-08/arne: coded.
!
      use MpiComm, only: stop_it
      use EquationOfState, only: cs0, cs20

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(ny,nz) :: dur_dr,duphi_dr,dlnrho_dr,&
                                dur_dphi,duphi_dphi,dlnrho_dphi,rho0,&
                                r_1,u_r,u_phi,&
                                L_1,L_5,L_4
      integer lll,i,nnn,mmm
      real, parameter :: sigma = 1.
      real :: u_T

      intent(in) :: f
      intent(out) :: df
  
      !lll = l2
      call der_onesided_yzslice(f,-1,ilnrho,dlnrho_dr,lll,1)
      call der_onesided_yzslice(f,-1,iux,dur_dr,lll,1)
      call der_onesided_yzslice(f,-1,iuz,duphi_dr,lll,1)
      call der_onesided_yzslice(f,-1,ilnrho,dlnrho_dphi,lll,3)
      call der_onesided_yzslice(f,-1,iux,dur_dphi,lll,3)
      call der_onesided_yzslice(f,-1,iuz,duphi_dphi,lll,3)

      ! ``dp = cs20*rho0*dlnrho''
      r_1 = 1./x(lll)
      u_r = f(lll,m1:m2,n1:n2,iux)
      u_phi = f(lll,m1:m2,n1:n2,iuz)
      rho0 = exp(f(lll,m1:m2,n1:n2,ilnrho))
      u_T = special_infuu

      ! It would be better to insert the expression for L_1 into the eqs. for
      ! df below, but it is more convenient to set L_1 explicitly for testing
      ! purposes.
      L_5 = (u_r + cs0)*(cs20*rho0*dlnrho_dr + rho0*cs0*dur_dr)
      !L_1 = 0 ! This one causes drift...
      !L_1 = -cs20*rho0*(2*u_r*r_1) ! This one reflects!!
      L_1 = -cs20*rho0*(2*u_r*r_1 + duphi_dphi*r_1) ! Chrashes around phi=0 and pi 
                                                     ! after some time. Best sofar.
      !L_1 = -cs0*rho0*r_1*u_phi**2 - cs20*rho0*r_1*(2*u_r + duphi_dphi) ! Drift
      !L_1 = -cs0*rho0*r_1*u_phi**2 - cs20*rho0*r_1*(2*u_r) ! Drift and reflection!
      !L_1 = cs0*rho0*r_1*(u_phi*dur_dphi - u_phi**2) &
      !      - cs20*rho0*r_1*(2*u_r + duphi_dphi) ! Nix
      !L_1 = cs0*rho0*r_1*(u_phi*dur_dphi - u_phi**2) &
      !      -u_phi*rho0*r_1*dlnrho_dphi- cs20*rho0*r_1*(2*u_r + duphi_dphi) ! Nix
      
      
      where (u_r < 0)
        !L_4 = -u_phi*r_1*duphi_dphi
        L_4 = 0
        !L_4 = u_r*duphi_dr
      elsewhere
        L_4 = u_r*duphi_dr
      endwhere

      if (linlet_northern .and. northern) then
        do mmm=1,ny
          L_1(mmm,1:nequator-1-nghost) = L_1(mmm,1:nequator-1-nghost)&
            -sigma*cs20*rho0(mmm,1:nequator-1-nghost)&
            *(u_r(mmm,1:nequator-1-nghost)-sin(z(n1:nequator-1))*(u_T))

          L_4(mmm,1:nequator-1-nghost) = L_4(mmm,1:nequator-1-nghost)&
            +sigma*cs20*rho0(mmm,1:nequator-1-nghost)&
            *(u_phi(mmm,1:nequator-1-nghost)-cos(z(n1:nequator-1))*(u_T))
        enddo
      endif

      df(lll,m1:m2,n1:n2,ilnrho) =    & 
          itsub*(                     &
          -1./(2.*cs20*rho0)*(L_5+L_1)&
          - u_phi*r_1*dlnrho_dphi     &
          - 2*u_r*r_1                 &
          - duphi_dphi*r_1            &
          )
      df(lll,m1:m2,n1:n2,iux) =       &
          itsub*(                     &
          -1./(2.*cs0*rho0)*(L_5-L_1) &
          - u_phi*r_1*dur_dphi        &
          + u_phi**2*r_1              &
          )
      df(lll,m1:m2,n1:n2,iuz) =       &
          itsub*(                     &
          -L_4                        &
          - u_phi*r_1*duphi_dphi      &
          - cs20*r_1*dlnrho_dphi      &
          - u_r*u_phi*r_1             &
          )
    endsubroutine
!***********************************************************************
    subroutine bc_nscbc_wall(f,df,pos)
!
!   Set dlnrho at a wall normal to
!   to r-direction acc. to characteristic boundary treatment.
!   Uses a one-sided finite diff. stencil.
!
!   7-jul-08/arne: coded.
!
      use MpiComm, only: stop_it
      use EquationOfState, only: cs0, cs20

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(ny,nz) :: dur_dr,dlnrho_dr,&
                                duphi_dphi,dlnrho_dphi,&
                                r_1,u_r,u_phi,rho0,&
                                L_1
      integer pos,lll,i

      intent(in) :: f,pos
      intent(inout) :: df
  
      lll = pos
      call der_onesided_yzslice(f,+1,ilnrho,dlnrho_dr,lll,1)
      call der_onesided_yzslice(f,+1,iux,dur_dr,lll,1)
      call der_onesided_yzslice(f,-1,iuz,duphi_dphi,lll,3)
      call der_onesided_yzslice(f,-1,iuz,duphi_dphi,lll,3)

      ! ``dp = cs20*rho0*dlnrho''
      r_1 = 1./x(lll)
      u_r = f(lll,m1:m2,n1:n2,iux)
      u_phi = f(lll,m1:m2,n1:n2,iuz)
      rho0 = exp(f(lll,m1:m2,n1:n2,ilnrho))

      L_1 = (u_r - cs0)*(cs20*rho0*dlnrho_dr - rho0*cs0*dur_dr)
      ! L_5 = L_1

      df(lll,m1:m2,n1:n2,ilnrho) = &
          -1./(cs20*rho0)*L_1      &
          !- u_phi*r_1*dlnrho_dphi  &
          !- 2*u_r*r_1              &
          - duphi_dphi*r_1
      df(lll,m1:m2,n1:n2,iux) = 0
      df(lll,m1:m2,n1:n2,iuy) = 0
      df(lll,m1:m2,n1:n2,iuz) = 0
    endsubroutine
!***********************************************************************
    subroutine sph_vel(f,val)
!
!   Set velocity inside sphere.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: val

      call sph_set(f,val,iux)
      call sph_set(f,val,iuy)
      call sph_set(f,val,iuz)
    endsubroutine
!***********************************************************************
    subroutine sph_set(f,val,j)
!
!   Set f-array variable inside sphere. If j is not specified, all
!   variables are set to val.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer, optional :: j
      real, intent(in) :: val

      if (present(j)) then
        f(sph_l1:sph_l2,:,sph_n1:sph_n2,j) = val
      else
        f(sph_l1:sph_l2,:,sph_n1:sph_n2,:) = val
      endif
    endsubroutine
!***********************************************************************
    subroutine sph_ghosts(f)
!
!   Set up ghost sonez inside sphere.
!   Antisymmetric velocities, and symmetric density.
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer sgn,i

      do i=1,nghost
        f(sph_l2-i,:,sph_n1:sph_n2,iux:iuz) = &
            -f(sph_l2+i,:,sph_n1:sph_n2,iux:iuz)
        f(sph_l2-i,:,sph_n1:sph_n2,ilnrho) = &
            f(sph_l2+i,:,sph_n1:sph_n2,ilnrho)
      end do
    endsubroutine
!***********************************************************************
    subroutine sph_ghosts_1s(f)
!
!   Set up ghost zones inside sphere s.t. the two grid points outside
!   sphere gets x-derivative as if calculated by ``semi-onesided''
!   stencils.
!   These expressions result from combining Eqs(207)-(209), astro-ph/0109497,
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer sgn,i,k

      do i=1,2
        k = sph_l2-i
        f(k,:,sph_n1:sph_n2,:) = &
            +  7*f(k+1,:,sph_n1:sph_n2,:)&
            - 21*f(k+2,:,sph_n1:sph_n2,:)&
            + 35*f(k+3,:,sph_n1:sph_n2,:)&
            - 35*f(k+4,:,sph_n1:sph_n2,:)&
            + 21*f(k+5,:,sph_n1:sph_n2,:)&
            -  7*f(k+6,:,sph_n1:sph_n2,:)&
            +  1*f(k+7,:,sph_n1:sph_n2,:)
      enddo

    endsubroutine
!***********************************************************************
    subroutine sph_ghosts_s1s(f)
!
!   Set up ghost zones inside sphere s.t. the two grid points outside
!   sphere gets 1st x-derivative as if calculated by ``semi-onesided''
!   stencils.
!   These expressions result from combining Eqs(207)-(209), astro-ph/0109497,
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      integer sgn,i

      f(sph_l2-1,:,sph_n1:sph_n2,:) = &
          +  7*f(sph_l2  ,:,sph_n1:sph_n2,:)&
          - 21*f(sph_l2+1,:,sph_n1:sph_n2,:)&
          + 35*f(sph_l2+2,:,sph_n1:sph_n2,:)&
          - 35*f(sph_l2+3,:,sph_n1:sph_n2,:)&
          + 21*f(sph_l2+4,:,sph_n1:sph_n2,:)&
          -  7*f(sph_l2+5,:,sph_n1:sph_n2,:)&
          +  1*f(sph_l2+6,:,sph_n1:sph_n2,:)

      f(sph_l2-2,:,sph_n1:sph_n2,:) = &
          + 28*f(sph_l2  ,:,sph_n1:sph_n2,:)&
          -112*f(sph_l2+1,:,sph_n1:sph_n2,:)&
          +210*f(sph_l2+2,:,sph_n1:sph_n2,:)&
          -224*f(sph_l2+3,:,sph_n1:sph_n2,:)&
          +144*f(sph_l2+4,:,sph_n1:sph_n2,:)&
          - 48*f(sph_l2+5,:,sph_n1:sph_n2,:)&
          +  7*f(sph_l2+6,:,sph_n1:sph_n2,:)
    endsubroutine
!***********************************************************************
    integer function binary_search(value,vector,lowstart,highstart)
!
! Returns the element index of element nearest to value in vector.
!
      real :: value
      real, dimension (:) :: vector
      integer :: lowstart,highstart
      integer :: low,high,mid

      low = lowstart
      high = highstart

      do while(high>low)
        mid = (low+high)/2
        if (value > vector(mid)) then
          low = mid+1
        else if (value < vector(mid)) then
          high = mid
        else
          exit
        end if
      end do

      binary_search = low
    endfunction
!***********************************************************************
    subroutine der_onesided_yzslice(f,sgn,k,df,pos,j)
      use Cdata
!
!   Uses a one-sided 4th order stencil for x-derivative, 6th order
!   central for z-derivative...
!   sgn = +1 for forward difference, sgn = -1 for backwards difference.
!
!   NB! Because of its original intended use in relation to solving
!   characteristic equations on boundaries, this sub returns PARTIAL
!   derivatives, NOT COVARIANT.
!
!   7-jul-08/arne: coded.
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:,:) :: df
      real :: fac,r_1
      integer :: pos,k,sgn,j
      integer :: nnn

      intent(in)  :: f,k,pos,sgn,j
      intent(out) :: df

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
          if (ip<=5) print*, 'der_onesided_yzslice: Degenerate case in x-direction'
        endif
      elseif (j==2) then
      elseif (j==3) then
        if (nzgrid/=1) then
          do nnn=n1,n2
            fac=(1./60)*dz_1(nnn)
            df(:,nnn-nghost)=&
                  fac*(+ 45.0*(f(pos,m1:m2,nnn+1,k)-f(pos,m1:m2,nnn-1,k)) &
                       -  9.0*(f(pos,m1:m2,nnn+2,k)-f(pos,m1:m2,nnn-2,k)) &
                       +      (f(pos,m1:m2,nnn+3,k)-f(pos,m1:m2,nnn-3,k)))
          enddo
        else
          df=0.
          if (ip<=5) print*, 'der_onesided_yzslice: Degenerate case in z-direction'
        endif
      endif
    endsubroutine
!***********************************************************************

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

