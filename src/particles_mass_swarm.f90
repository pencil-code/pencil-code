! $Id$
!
! This module handles the mass of super-particles.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
!
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MPVAR CONTRIBUTION 1
!
! CPARAM logical, parameter :: lparticles_mass = .true.
!
!***************************************************************
module Particles_mass
!
  use Cdata
  use Cparam
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
  use Particles_sub
!
  implicit none
!
  include 'particles_mass.h'
!
! Diagnostic variables
!
  integer :: idiag_mpm = 0    ! DIAG_DOC: $\overline{m_p}$
  integer :: idiag_mpmin = 0  ! DIAG_DOC: $\min_j m_{p,j}$
  integer :: idiag_mpmax = 0  ! DIAG_DOC: $\max_j m_{p,j}$
!
  contains
!***********************************************************************
    subroutine register_particles_mass()
!
! Set up indices for access to the fp and dfp arrays.
!
! 18-jun-17/ccyang: coded
!
      if (lroot) call svn_id("$Id$")
!
! Index for particle mass.
!
      call append_npvar('imp',imp)
!
    endsubroutine register_particles_mass
!***********************************************************************
    subroutine initialize_particles_mass(f)
!
! Perform any post-parameter-read initialization i.e. calculate derived
! parameters.
!
! 18-jun-17/ccyang: dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_mass
!***********************************************************************
    subroutine init_particles_mass(f, fp)
!
! Initialize particle mass.
!
! 18-jun-17/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mpar_loc,mparray), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
! Assign mp_swarm to all particles.
!
      fp(1:npar_loc,imp) = mp_swarm
!
    endsubroutine init_particles_mass
!***********************************************************************
    subroutine pencil_criteria_par_mass()
!
! All pencils that the Particles_mass module depends on are specified
! here.
!
! 18-jun-17/ccyang: dummy
!
    endsubroutine pencil_criteria_par_mass
!***********************************************************************
    subroutine dpmass_dt(f, df, fp, dfp, ineargrid)
!
! Evolution of particle mass.
!
! 04-jul-17/ccyang: coded
!
      use Particles_sub, only: sum_par_name, max_par_name
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(in) :: df
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mpar_loc,mpvar), intent(in) :: dfp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
! Diagnostic output
!
      diag: if (ldiagnos) then
        if (idiag_mpm /= 0) call sum_par_name(fp(1:npar_loc,imp), idiag_mpm)
        if (idiag_mpmin /= 0) call max_par_name(-fp(1:npar_loc,imp), idiag_mpmin, lneg=.true.)
        if (idiag_mpmax /= 0) call max_par_name(fp(1:npar_loc,imp), idiag_mpmax)
      endif diag
!
    endsubroutine dpmass_dt
!***********************************************************************
    subroutine dpmass_dt_pencil(f, df, fp, dfp, p, ineargrid)
!
! Evolution of particle mass in pencils.
!
! 04-jul-17/ccyang: dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(mx,my,mz,mvar), intent(in) :: df
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      real, dimension(mpar_loc,mpvar), intent(in) :: dfp
      type(pencil_case), intent(in) :: p
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine dpmass_dt_pencil
!***********************************************************************
    subroutine read_particles_mass_init_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_mass_init_pars
!***********************************************************************
    subroutine write_particles_mass_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_mass_init_pars
!***********************************************************************
    subroutine read_particles_mass_run_pars(iostat)
!
      integer, intent(out) :: iostat
!
      iostat = 0
!
    endsubroutine read_particles_mass_run_pars
!***********************************************************************
    subroutine write_particles_mass_run_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_mass_run_pars
!***********************************************************************
    subroutine rprint_particles_mass(lreset, lwrite)
!
! Read and register print parameters relevant for particle mass.
!
! 04-jul-17/ccyang: coded
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
!
      logical, intent(in) :: lreset
      logical, intent(in), optional :: lwrite
!
      logical :: lwr
      integer :: iname
!
! Write information to index.pro.
!
      lwr = .false.
      if (present(lwrite)) lwr = lwrite
      if (lwr) call farray_index_append('imp', imp)
!
! Reset diagnostic variables if requested.
!
      reset: if (lreset) then
        idiag_mpm = 0
        idiag_mpmin = 0
        idiag_mpmax = 0
      endif reset
!
! Parse the names in print.in.
!
      parse: do iname = 1, nname
        call parse_name(iname, cname(iname), cform(iname), 'mpm', idiag_mpm)
        call parse_name(iname, cname(iname), cform(iname), 'mpmin', idiag_mpmin)
        call parse_name(iname, cname(iname), cform(iname), 'mpmax', idiag_mpmax)
      enddo parse
!
    endsubroutine rprint_particles_mass
!***********************************************************************
endmodule Particles_mass
