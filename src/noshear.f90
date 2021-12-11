! $Id$
!
!  This modules deals with all aspects of shear; if no
!  shear is invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  shear relevant subroutines listed in here.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lshear = .false.
!
!***************************************************************
module Shear
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'shear.h'
!
  contains
!***********************************************************************
    subroutine register_shear
!
!  Initialise variables.
!
!  2-july-02/nils: coded
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_shear
!***********************************************************************
    subroutine initialize_shear
!
!  21-nov-02/tony: coded
!  17-jul-04/axel: Sshear=0 is needed for forcing_hel to work correctly.
!
      Sshear=0.0
!
    endsubroutine initialize_shear
!***********************************************************************
    subroutine read_shear_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_shear_init_pars
!***********************************************************************
    subroutine write_shear_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_shear_init_pars
!***********************************************************************
    subroutine read_shear_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_shear_run_pars
!***********************************************************************
    subroutine write_shear_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_shear_run_pars
!***********************************************************************
    subroutine shear_before_boundary(f)
!
!  Actions to take before boundary conditions are set.
!
!   1-may-08/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine shear_before_boundary
!***********************************************************************
    subroutine pencil_criteria_shear
!
!  All pencils that the Shear module depends on are specified here.
!
!  01-may-09/wlad: coded
!
    endsubroutine pencil_criteria_shear
!***********************************************************************
    subroutine pencil_interdep_shear(lpencil_in)
!
!  Interdependency among pencils from the Shear module is specified here.
!
!  01-may-09/wlad: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_shear
!***********************************************************************
    subroutine calc_pencils_shear(f,p)
!
!  Calculate Shear pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  01-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_shear
!***********************************************************************
    subroutine shearing(f,df,p)
!
!  Calculates the actual shear terms
!
!  2-july-02/nils: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(p)
!
    endsubroutine shearing
!***********************************************************************
    subroutine advance_shear(f,df,dt_shear)
!
!  Dummy routine: deltay remains unchanged.
!
!  18-aug-02/axel: incorporated from nompicomm.f90
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real :: dt_shear
!
!  Print identifier.
!
      if (headtt.or.ldebug) print*,'advance_shear: deltay=const=',deltay
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dt_shear)
!
    endsubroutine advance_shear
!***********************************************************************
    subroutine sheared_advection_fft(a, comp_start, comp_end, dt_shear)
!
!  dummy
!
      integer, intent(in) :: comp_start, comp_end
      real, dimension(:,:,:,:), intent(inout) :: a
      real, intent(in) :: dt_shear
!
      call keep_compiler_quiet(a)
      call keep_compiler_quiet(comp_start,comp_end)
      call keep_compiler_quiet(dt_shear)
!
    endsubroutine sheared_advection_fft
!***********************************************************************
    subroutine boundcond_shear(f,ivar1,ivar2)
!
!  Dummy routine.
!
!  01-oct-07/anders: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(ivar1)
      call keep_compiler_quiet(ivar2)
!
    endsubroutine boundcond_shear
!***********************************************************************
    subroutine shear_variables(f,df,nvars,jstart,jstep,shear1)
!
!  Dummy routine.
!
!  28-apr-11/wlad: coded
!  02-aug-11/MR: added optionals jstep,shear1
!
      real, dimension(mx,my,mz,mfarray), intent(in)           :: f
      real, dimension(mx,my,mz,mvar)   , intent(inout)        :: df
      integer,                           intent(in)           :: nvars, jstart
      integer,                           intent(in), optional :: jstep
      logical,                           intent(in), optional :: shear1
!
      df = df+0.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(nvars,jstart)
      if ( present(jstep)  ) call keep_compiler_quiet(jstep)
      if ( present(shear1) ) call keep_compiler_quiet(shear1)
!
    endsubroutine shear_variables
!***********************************************************************
    subroutine rprint_shear(lreset,lwrite)
!
!  Dummy routine.
!
!  02-jul-04/tobi: coded
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_shear
!***********************************************************************
    subroutine get_uy0_shear(uy0_shear, x)
!
!  Gets the shear velocity.
!
!  08-oct-13/ccyang: coded
!
      real, dimension(:), intent(out) :: uy0_shear
      real, dimension(:), intent(in), optional :: x
!
      uy0_shear = 0.0
      if (present(x)) call keep_compiler_quiet(x)
!
    endsubroutine
!***********************************************************************
    subroutine get_hyper3x_mesh(lhyper3x_mesh_out, diff_hyper3x_mesh_out)
!
!  Gets module variables lhyper3x_mesh and diff_hyper3x_mesh.
!
!  04-jun-14/ccyang: coded
!
      logical, intent(out) :: lhyper3x_mesh_out
      real, intent(out) :: diff_hyper3x_mesh_out
!
      lhyper3x_mesh_out = .false.
      diff_hyper3x_mesh_out = 0.0
!
    endsubroutine get_hyper3x_mesh
!***********************************************************************
    subroutine shear_frame_transform(a,tshift)
!
!  Transforms a variable a from lab frame to shear frame in the x space
!
!  31-jun-21/hongzhe: coded
!
      real, dimension(nx,ny,nz), intent(inout) :: a
      rreal, optional :: tshift
!
      call keep_compiler_quiet(a)
!
    endsubroutine shear_frame_transform
!***********************************************************************
endmodule Shear
