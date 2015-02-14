! $Id$
!
!  This module takes care of drag forces between particles and gas.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_drag = .true.
!
!***************************************************************
module Particles_drag
!
  use Cdata
  use Cparam
  use Messages
  use Particles_cdata
!
  implicit none
!
  include 'particles_drag.h'
!
!  Initialization parameters
!
  real :: taus = 0.0
!
  namelist /particles_drag_init_pars/ taus
!
!  Runtime parameters
!
  namelist /particles_drag_run_pars/ taus
!
  contains
!***********************************************************************
    subroutine register_particles_drag()
!
!  Register this module.
!
!  14-dec-14/ccyang: coded.
!
      if (lroot) call svn_id("$Id$")
!
    endsubroutine register_particles_drag
!***********************************************************************
    subroutine initialize_particles_drag()
!
!  Perform any post-parameter-read initialization, i.e., calculate
!  derived parameters.
!
!  14-dec-14/ccyang: coded.
!
    endsubroutine initialize_particles_drag
!***********************************************************************
    subroutine read_particles_drag_init_pars(unit, iostat)
!
!  Read initialization parameters from namelist particles_drag_init_pars.
!
!  14-feb-15/ccyang: coded.
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      integer :: stat
!
      read(unit, NML=particles_drag_init_pars, IOSTAT=stat)
      if (present(iostat)) then
        iostat = stat
      else if (stat /= 0) then
        call fatal_error('read_particles_drag_init_pars', 'cannot read particles_drag_init_pars. ')
      endif
!
    endsubroutine read_particles_drag_init_pars
!***********************************************************************
    subroutine write_particles_drag_init_pars(unit)
!
!  Write initialization parameters to namelist particles_drag_init_pars.
!
!  14-feb-15/ccyang: coded.
!
      integer, intent(in) :: unit
!
      integer :: stat
!
      write(unit, NML=particles_drag_init_pars, IOSTAT=stat)
      if (stat /= 0) call fatal_error('write_particles_drag_init_pars', 'cannot write particles_drag_init_pars. ')
!
    endsubroutine write_particles_drag_init_pars
!***********************************************************************
    subroutine read_particles_drag_run_pars(unit, iostat)
!
!  Read runtime parameters from namelist particles_drag_run_pars.
!
!  14-dec-14/ccyang: coded.
!
      include 'unit.h'
      integer, intent(inout), optional :: iostat
!
      integer :: stat
!
      read(unit, NML=particles_drag_run_pars, IOSTAT=stat)
      if (present(iostat)) then
        iostat = stat
      else if (stat /= 0) then
        call fatal_error('read_particles_drag_run_pars', 'cannot read particles_drag_run_pars. ')
      endif
!
    endsubroutine read_particles_drag_run_pars
!***********************************************************************
    subroutine write_particles_drag_run_pars(unit)
!
!  Write runtime parameters to namelist particles_drag_run_pars.
!
!  14-dec-14/ccyang: coded.
!
      integer, intent(in) :: unit
!
      integer :: stat
!
      write(unit, NML=particles_drag_run_pars, IOSTAT=stat)
      if (stat /= 0) call fatal_error('write_particles_drag_run_pars', 'cannot write particles_drag_run_pars. ')
!
    endsubroutine write_particles_drag_run_pars
!***********************************************************************
    subroutine integrate_drag(f, fp)
!
!  Wrapper for the integration of the drag force between particles and
!  gas.
!
!  05-feb-15/ccyang: coded.
!
      use General, only: keep_compiler_quiet
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine integrate_drag
!***********************************************************************
endmodule Particles_drag
