! $Id$

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lgrav = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED gg(3)
!
!***************************************************************
module Gravity
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'gravity.h'
!
  interface potential
    module procedure potential_global
    module procedure potential_penc
    module procedure potential_point
  endinterface
!
  interface acceleration
    module procedure acceleration_penc
    module procedure acceleration_penc_1D
    module procedure acceleration_point
  endinterface
!
  real :: z1,z2,zref,zgrav,gravz,zinfty,nu_epicycle=1.
  real :: lnrho_bot,lnrho_top,ss_bot,ss_top
  real :: gravz_const=1.,reduced_top=1.
  real :: g0=0.,r0_pot=0.,qgshear=1.5
  integer :: n_pot=10
  character (len=labellen) :: gravz_profile='const'  !(used by Density)
  logical :: lnumerical_equilibrium=.false.
!
  contains
!***********************************************************************
    subroutine register_gravity()
!
!  Initialise gravity flags.
!
!  9-jan-02/wolf: coded
! 28-mar-02/axel: adapted from grav_z
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
!
      lgravz = .false.
      lgravr = .false.
!
    endsubroutine register_gravity
!***********************************************************************
    subroutine initialize_gravity(f,lstarting)
!
!  Set up some variables for gravity; do nothing in nograv
!  16-jul-02/wolf: coded
!  22-nov-02/tony: renamed from setup_grav
!
      real, dimension(mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_gravity
!***********************************************************************
    subroutine read_gravity_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_gravity_init_pars
!***********************************************************************
    subroutine write_gravity_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_gravity_init_pars
!***********************************************************************
    subroutine read_gravity_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_gravity_run_pars
!***********************************************************************
    subroutine write_gravity_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_gravity_run_pars
!***********************************************************************
    subroutine init_gg(f)
!
!  initialise gravity; called from start.f90
!   9-jan-02/wolf: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
! Not doing anything (this might change if we decide to store gg)
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_gg
!***********************************************************************
    subroutine pencil_criteria_gravity()
!
!  All pencils that the Gravity module depends on are specified here.
!
!  20-11-04/anders: coded
!
    endsubroutine pencil_criteria_gravity
!***********************************************************************
    subroutine pencil_interdep_gravity(lpencil_in)
!
!  Interdependency among pencils from the Gravity module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_gravity
!***********************************************************************
    subroutine calc_pencils_gravity(f,p)
!
!  Calculate Gravity pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  12-nov-04/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      if (lpencil(i_gg)) p%gg=0.
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_gravity
!***********************************************************************
    subroutine duu_dt_grav(f,df,p)
!
!  Add nothing to duu/dt.
!
!  28-mar-02/axel: adapted from grav_z
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential_global(pot,pot0)
!
!  Gravity potential.
!
!  28-mar-02/axel: adapted from grav_z
!
      real, dimension (mx,my,mz) :: pot
      real, optional :: pot0
!
      intent(out) :: pot,pot0
!
      if (lroot) print*,'potential: note, GRAV=nograv is not OK'
      pot = 0.
      if (present(pot0)) pot0 = 0.
!
    endsubroutine potential_global
!***********************************************************************
    subroutine potential_penc(xmn,ymn,zmn,pot,pot0,grav,rmn)
!
!  Gravity potential.
!
!  28-mar-02/axel: adapted from grav_z
!
      real, dimension (nx) :: pot
      real, optional :: ymn,zmn,pot0
      real, optional, dimension (nx) :: xmn,rmn
      real, optional, dimension (nx,3) :: grav
!
      intent(in) :: xmn,ymn,zmn,rmn
      intent(out) :: pot,pot0
!
      if (lroot) print*,'potential: note, GRAV=nograv is not OK'
      pot = 0.
      if (present(pot0)) pot0 = 0.
!
      call keep_compiler_quiet(present(xmn))
      call keep_compiler_quiet(present(ymn))
      call keep_compiler_quiet(present(zmn))
      call keep_compiler_quiet(present(rmn))
      call keep_compiler_quiet(present(grav))
!
    endsubroutine potential_penc
!***********************************************************************
    subroutine potential_point(x,y,z,r, pot,pot0, grav)
!
!  Gravity potential in one point.
!
!  20-dec-03/wolf: coded
!
      real :: pot
      real, optional :: x,y,z,r
      real, optional :: pot0,grav
!
      intent(in)  :: x,y,z,r
      intent(out) :: pot
!
      call fatal_error("nograv","potential_point not implemented")
!
      pot = 0.
!
      call keep_compiler_quiet(present(x))
      call keep_compiler_quiet(present(y))
      call keep_compiler_quiet(present(z))
      call keep_compiler_quiet(present(r))
      call keep_compiler_quiet(present(pot0))
      call keep_compiler_quiet(present(grav))
!
    endsubroutine potential_point
!***********************************************************************
    subroutine acceleration_penc(gg)
!
!  Calculates gravitational acceleration on a pencil.
!
!  21-apr-07/tobi: adapted from potential_penc
!
      real, dimension (:,:), intent (out) :: gg
!
!  Calculate acceleration from master pencils defined in initialize_gravity.
!
      call fatal_error("acceleration_penc","Not implemented")
!
      call keep_compiler_quiet(gg)
!
    endsubroutine acceleration_penc
!***********************************************************************
    subroutine acceleration_penc_1D(gr)
!
!  Calculates gravitational acceleration on a pencil.
!
!  21-apr-07/tobi: adapted from potential_penc
!
      real, dimension (nx), intent (out) :: gr
!
!  Calculate acceleration from master pencils defined in initialize_gravity.
!
      call fatal_error("acceleration_penc_1D","Not implemented")
!
      call keep_compiler_quiet(gr)
!
    endsubroutine acceleration_penc_1D
!***********************************************************************
    subroutine acceleration_point(x,y,z,r,g_r)
!
!  Gravity in one point.
!
!  18-nov-08/wlad: coded
!
      real :: g_r
      real, optional :: x,y,z,r
!
      intent(in)  :: x,y,z,r
      intent(out) :: g_r
!
      call fatal_error('acceleration_penc_1D','Not implemented')
!
      g_r = 0.0
!
      call keep_compiler_quiet(present(x))
      call keep_compiler_quiet(present(y))
      call keep_compiler_quiet(present(z))
      call keep_compiler_quiet(present(r))
!
    endsubroutine acceleration_point
!***********************************************************************
    subroutine rprint_gravity(lreset,lwrite)
!
!  Reads and registers print parameters relevant for gravity advance.
!
!  26-apr-03/axel: dummy
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Write column, idiag_XYZ, where our variable XYZ is stored.
!  IDL needs this even if everything is zero.
!
      if (lwr) then
        write(3,*) 'igg=',igg
        write(3,*) 'igx=',igx
        write(3,*) 'igy=',igy
        write(3,*) 'igz=',igz
      endif
!
      call keep_compiler_quiet(lreset)
!
    endsubroutine rprint_gravity
!***********************************************************************
    subroutine compute_gravity_star(f, wheat, luminosity, star_cte)
!
!  5-jan-10/boris: coded
!
    real, dimension (mx,my,mz,mfarray) :: f
    real :: wheat, luminosity, star_cte
!
    call keep_compiler_quiet(f)
    call keep_compiler_quiet(wheat,luminosity,star_cte)
!
    endsubroutine compute_gravity_star
!***********************************************************************
    subroutine get_xgravity(xgrav)
!
!  Dummy routine for getting the gravity profile into the
!  intial conditions
!
!  04-oct-10/bing: coded
!
      real, dimension(nx) :: xgrav
!
      call keep_compiler_quiet(xgrav)
!
    endsubroutine get_xgravity
!***********************************************************************
    subroutine set_consistent_gravity(ginput,gtype,gprofile,lsuccess)
      real :: ginput
      character (len=labellen) :: gtype,gprofile
      character (len=labellen) :: gprof
      logical :: lsuccess
      logical :: lconsistent=.true.
!
! This routine should never be called in the way it is written now.
!
      lsuccess=.false.
!
! gravity parameters set consistently.
!
    endsubroutine set_consistent_gravity
!***********************************************************************
endmodule Gravity
