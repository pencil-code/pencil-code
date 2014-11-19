! $Id: particles_coagulation.f90 19828 2012-11-27 09:58:06Z kalle.jansson.89 $
!
!  This modules takes care of adapting the number of particles in a grid cell
!  to a desired value. This module is based on an original idea by Jacob Trier
!  Frederiksen.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_adaptation= .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
module Particles_adaptation
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
  include 'particles_adaptation.h'
!
  contains
!***********************************************************************
    subroutine initialize_particles_adaptation(f)
!
!  03-apr-13/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_particles_adaptation
!***********************************************************************
    subroutine particles_adaptation_pencils(f,fp,dfp,ipar,ineargrid)
!
!  03-apr-13/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      real, dimension (mpar_loc,mpvar) :: dfp
      integer, dimension (mpar_loc) :: ipar
      integer, dimension (mpar_loc,3) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ipar)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine particles_adaptation_pencils
!***********************************************************************
    subroutine read_particles_adapt_run_pars(unit,iostat)
!
!  03-apr-13/anders: dummy
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      if (present(iostat)) call keep_compiler_quiet(iostat)
!
    endsubroutine read_particles_adapt_run_pars
!***********************************************************************
    subroutine write_particles_adapt_run_pars(unit)
!
!  03-apr-13/anders: dummy
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_particles_adapt_run_pars
!*******************************************************************
    subroutine rprint_particles_adaptation(lreset,lwrite)
!
!  03-apr-13/anders: dummy
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      if (present(lwrite)) call keep_compiler_quiet(lwrite)
!
    endsubroutine rprint_particles_adaptation
!***********************************************************************
endmodule Particles_adaptation
