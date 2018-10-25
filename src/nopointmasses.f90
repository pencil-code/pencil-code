! $Id$
!
!  This module takes care of direct N-body gravity between point masses.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lpointmasses=.false.
!
!***************************************************************
module PointMasses
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id
!
  implicit none
!
  include 'pointmasses.h'
!
  contains
!***********************************************************************
    subroutine register_pointmasses
!
!  Set up indices for access to the f and fq.
!
!  27-aug-06/wlad: adapted
!
      if (lroot) call svn_id( &
          "$Id$")
!
    endsubroutine register_pointmasses
!***********************************************************************
    subroutine initialize_pointmasses(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_pointmasses
!***********************************************************************
    subroutine pencil_criteria_pointmasses
!
    endsubroutine pencil_criteria_pointmasses
!***********************************************************************
    subroutine pencil_interdep_pointmasses(lpencil_in)
!
!  Interdependency among pencils provided by the pointmasses module
!  is specified here.
!
!  22-sep-06/wlad: adapted
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_pointmasses
!***********************************************************************
    subroutine calc_pencils_pointmasses(f,p)
!
!  Calculate point mass particle pencils
!
!  22-sep-06/wlad: adapted
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_pointmasses
!***********************************************************************
    subroutine init_pointmasses(f)
!
!  Initial positions and velocities of point mass particles.
!  Overwrite the position asserted by the dust module
!
!  17-nov-05/anders+wlad: adapted
!
      use General, only: random_number_wrapper
      use Sub
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent (in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_pointmasses
!***********************************************************************
    subroutine pointmasses_pde_pencil(f,df,p)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine pointmasses_pde_pencil
!***********************************************************************         
    subroutine pointmasses_pde(f,df)
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      call keep_compiler_quiet(f,df)
!
    endsubroutine  pointmasses_pde
!***********************************************************************
    subroutine read_pointmasses_init_pars(iostat)
!
      integer, intent(out) :: iostat
      iostat=0
!
    endsubroutine read_pointmasses_init_pars
!***********************************************************************
    subroutine write_pointmasses_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_pointmasses_init_pars
!***********************************************************************
    subroutine read_pointmasses_run_pars(iostat)
!
      integer, intent(out) :: iostat
      iostat=0
!
    endsubroutine read_pointmasses_run_pars
!***********************************************************************
    subroutine write_pointmasses_run_pars(unit)
!
      integer, intent(in) :: unit
!      
      call keep_compiler_quiet(unit)
!
    endsubroutine write_pointmasses_run_pars
!***********************************************************************
    subroutine get_totalmass(tmass)
!
      real :: tmass
      call keep_compiler_quiet(tmass)
!
    endsubroutine get_totalmass
!***********************************************************************
    subroutine pointmasses_read_snapshot(filename)
!
      character (len=*) :: filename
      call keep_compiler_quiet(filename)
!
    endsubroutine pointmasses_read_snapshot
!***********************************************************************
    subroutine pointmasses_write_snapshot(file,enum,flist)
!
      character (len=*) :: file,flist
      logical :: enum
      optional :: flist
!
      call keep_compiler_quiet(file,flist)
      call keep_compiler_quiet(enum)
!
    endsubroutine pointmasses_write_snapshot
!***********************************************************************
    subroutine pointmasses_write_qdim(filename)
!
      character (len=*) :: filename
      call keep_compiler_quiet(filename)
!
    endsubroutine pointmasses_write_qdim
!***********************************************************************
    subroutine rprint_pointmasses(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset,lwrite)
!
    endsubroutine rprint_pointmasses
!***********************************************************************
    subroutine boundconds_pointmasses
!
    endsubroutine boundconds_pointmasses
!***********************************************************************
    subroutine pointmasses_timestep_first(f)
!    
      real, dimension (mx,my,mz,mfarray) :: f
      call keep_compiler_quiet(f)
!
   endsubroutine pointmasses_timestep_first
!***********************************************************************
    subroutine pointmasses_timestep_second(f)
!
      real, dimension (mx,my,mz,mfarray) :: f
      call keep_compiler_quiet(f)
!
    endsubroutine pointmasses_timestep_second
!***********************************************************************
  endmodule PointMasses
