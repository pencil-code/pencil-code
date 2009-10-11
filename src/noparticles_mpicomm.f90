! $Id: particles_sub.f90 11877 2009-10-11 13:20:36Z ajohan@strw.leidenuniv.nl $
!
!  This module deals with communication of particles between processors.
!
module Particles_mpicomm
!
  use Cdata
  use Messages
  use Particles_cdata
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  private
!
  public :: redist_particles_procs
!
  contains
!***********************************************************************
    subroutine redist_particles_procs(fp,npar_loc,ipar,dfp,linsert)
!
!  11-oct-09/anders: dummy
!
      use Mpicomm
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mpar_loc,mpvar), optional :: dfp
      logical, optional :: linsert
      integer, dimension (mpar_loc) :: ipar
      integer :: npar_loc
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(npar_loc)
      call keep_compiler_quiet(ipar)
      if (present(dfp)) call keep_compiler_quiet(dfp)
      if (present(linsert)) call keep_compiler_quiet(linsert)
!
    endsubroutine redist_particles_procs
!***********************************************************************
endmodule Particles_mpicomm
