! $Id$
!
!  This module deals with communication of particles between processors.
!
module Particles_mpicomm
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use Particles_cdata
!
  implicit none
!
  include 'particles_mpicomm.h'
!
  integer, parameter :: nxb=1, nyb=1, nzb=1, nbx=1, nby=1, nbz=1
  integer, parameter :: nbricks=0, nghostb=0, mxb=1,myb=1, mzb=1
  integer, parameter :: l1b=1, l2b=1, m1b=1, m2b=1, n1b=1, n2b=1
  integer, dimension (0:0) :: k1_iblock=0, k2_iblock=0, npar_iblock=0
  integer, dimension (1) :: inearblock
  integer, dimension (0:0) :: ibrick_parent_block, iproc_parent_block
  integer, dimension (0:0) :: iproc_foster_brick
  integer, dimension (1) :: iproc_parent_list, iproc_foster_list
  integer :: nbrick_foster=0, nproc_parent=0, nproc_foster=0, nblock_loc=0
  real, dimension (1,0:0) :: xbrick=0, ybrick=0, zbrick=0, xb=0, yb=0, zb=0
  real, dimension (1,0:0) :: dx1brick=0, dy1brick=0, dz1brick=0
  real, dimension (1,0:0) :: dx1b=0, dy1b=0, dz1b=0
  real, dimension (1,0:0) :: dVol1xbrick=0, dVol1ybrick=0, dVol1zbrick=0
  real, dimension (1,0:0) :: dVol1xb=0, dVol1yb=0, dVol1zb=0
  real, dimension (1,1,1,1,0:0) :: fb, dfb
  integer :: it1_loadbalance=1
  logical :: lfill_blocks_density=.false., lfill_blocks_velocity=.false.
  logical :: lfill_blocks_gpotself=.false., lfill_bricks_velocity=.false.
  logical :: lreblock_particles_run=.false., lbrick_partition=.false.
  logical :: ladopt_own_light_bricks=.false.
!
  contains
!***********************************************************************
    subroutine initialize_particles_mpicomm(f)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
!
!  Distribute particles evenly among processors to begin with.
!
      if (lstart) call dist_particles_evenly_procs(ipar)
!
!  Set neighbor process ranks for particles.
!
      neighbors_par = root
      if (bcpx /= 'p' .or. nxgrid == 1 .or. lshear .and. nygrid > 1) neighbors_par((/-1,1/),:,:) = -1
      if (bcpy /= 'p' .or. nygrid == 1) neighbors_par(:,(/-1,1/),:) = -1
      if (bcpz /= 'p' .or. nzgrid == 1) neighbors_par(:,:,(/-1,1/)) = -1
!
      call keep_compiler_quiet(f)
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
      integer :: k, jspec, npar_per_species
!
      intent (out) :: ipar
!
!  Set index interval of particles.
!
      if (linsert_particles_continuously) then
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
    subroutine migrate_particles(fp,ipar,dfp,linsert)
!
!  11-oct-09/anders: dummy
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
      logical, optional :: linsert
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ipar)
      if (present(dfp)) call keep_compiler_quiet(dfp)
      if (present(linsert)) call keep_compiler_quiet(linsert)
!
    endsubroutine migrate_particles
!***********************************************************************
    subroutine load_balance_particles(f,fp,ipar)
!
!  This subroutine counts particles in the bricks at the local processor
!  and distributes the bricks in such a away that there is approximately
!  equal number of particles per processor.
!
!  16-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc) :: ipar
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ipar)
!
    endsubroutine load_balance_particles
!***********************************************************************
    subroutine output_blocks(filename)
!
!  Write block domain decomposition to file.
!
!  16-nov-09/anders: dummy
!
      character(len=*) :: filename
!
      call keep_compiler_quiet(filename)
!
    endsubroutine output_blocks
!***********************************************************************
    subroutine input_blocks(filename)
!
!  Read block domain decomposition from file.
!
!  16-nov-09/anders: dummy
!
      character(len=*) :: filename
!
      call keep_compiler_quiet(filename)
!
    endsubroutine input_blocks
!***********************************************************************
    subroutine sort_blocks
!
!  Sort the blocks by parent processor and by parent brick index.
!
!  18-nov-09/anders: dummy
!
    endsubroutine sort_blocks
!***********************************************************************
    subroutine get_brick_index(xxp, iproc, ibrick)
!
!  10-jan-12/ccyang: dummy
!
      real, dimension(3), intent(in) :: xxp
      integer, intent(out) :: iproc, ibrick
!
      call keep_compiler_quiet(xxp)
      call keep_compiler_quiet(iproc)
      call keep_compiler_quiet(ibrick)
!
    endsubroutine get_brick_index
!***********************************************************************
    subroutine particles_neighbour_proc
!
! There is no neighbouring processors    
!      

    endsubroutine particles_neighbour_proc
!***********************************************************************
    subroutine communicate_fpbuf(to_neigh,from_neigh,her_npbuf,my_npbuf)
! dummy

      integer :: to_neigh,from_neigh,her_npbuf,my_npbuf
      call keep_compiler_quiet(to_neigh,from_neigh,her_npbuf,my_npbuf)
!
    endsubroutine communicate_fpbuf
!***********************************************************************
endmodule Particles_mpicomm
