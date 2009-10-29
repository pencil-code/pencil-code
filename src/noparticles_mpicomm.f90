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
  include 'particles_mpicomm.h'
!
  contains
!***********************************************************************
    subroutine initialize_particles_mpicomm(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      integer :: iblock
!
      intent (in) :: f, lstarting
!
!  Distribute particles evenly among processors to begin with.
!
      if (lstarting) call dist_particles_evenly_procs(ipar)
!
    endsubroutine initialize_particles_mpicomm
!***********************************************************************
    subroutine dist_particles_evenly_procs(ipar)
!
!  Distribute particles evenly among processors - single CPU version.
!
!  05-jan-05/anders: coded
!
      integer, dimension (mpar_loc) :: ipar
!
      integer :: i, k, jspec, npar_per_species
!
      intent (out) :: ipar
!
!  Set index interval of particles.
!
      if (lparticles_nbody.and.(npar==nspar)) then
        npar_loc=npar
        do k=1,nspar
          ipar(k)=k
          ipar_nbody(k)=k
        enddo
      else if (linsert_particles_continuously) then
        npar_loc=0
        ipar=0
      else
        npar_loc=npar
        do k=1,npar_loc
          ipar(k)=k
        enddo
      endif
!
!  Must have same number of particles in each species.
!
      if (mod(npar,npar_species)/=0) then
        print*, 'dist_particles_evenly_procs: npar_species '// &
            'must be a whole multiple of npar!'
        print*, 'npar_species, npar=', npar_species, npar
        call fatal_error('dist_particles_evenly_procs','')
      endif
!
      npar_per_species=npar/npar_species
!
!  Calculate right index for each species.
!
      ipar_fence_species(1)=npar_per_species
      ipar_fence_species(npar_species)=npar
      ipar_fence_species(1)=npar_per_species
      ipar_fence_species(npar_species)=npar
!
      do jspec=2,npar_species-1
        ipar_fence_species(jspec)= &
            ipar_fence_species(jspec-1)+npar_per_species
      enddo
!
      print*, 'dist_particles_evenly_procs: npar_per_species   =', &
          npar_per_species
      print*, 'dist_particles_evenly_procs: ipar_fence_species =', &
          ipar_fence_species
!
    endsubroutine dist_particles_evenly_procs
!***********************************************************************
    subroutine redist_particles_procs(fp,ipar,dfp,linsert)
!
!  11-oct-09/anders: dummy
!
      use Mpicomm
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
      logical, optional :: linsert
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ipar)
      if (present(dfp)) call keep_compiler_quiet(dfp)
      if (present(linsert)) call keep_compiler_quiet(linsert)
!
    endsubroutine redist_particles_procs
!***********************************************************************
endmodule Particles_mpicomm
