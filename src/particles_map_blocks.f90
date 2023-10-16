! $Id$
!
!  This module contains subroutines useful for mapping particles on the mesh.
!
!  This version is for block domain decomposition of particles.
!
!  In block domain decomposition the main mesh is divided among the processors
!  in the usual way. The local mesh is then divided into so-called "bricks",
!  small volumes of grid points. Particles are counted in each of those bricks,
!  and then the bricks are distributed among the processors so that each
!  processor has approximately the same number of particles. A brick fostered
!  by a processor is called a "block".
!
!  In each time-step the relevant dynamical variables must be transferred from
!  bricks at the parent processors to blocks at the foster processors. This
!  can be e.g. gas velocity field (for drag force) or gravitational potential
!  (for self-gravity).
!
!  A processor can open up new bricks in its own domain, if a particle moves
!  into an empty brick. Full load balancing is performed at regular intervals.
!  Here each processor count particles in their blocks and sends the
!  information to the parent processors. The parent processors then decide on a
!  new division of particle blocks.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lparticles_blocks = .true.
!
!***************************************************************
module Particles_map
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Mpicomm
  use Messages
  use Particles_cdata
  use Particles_mpicomm
!
  implicit none
!
  include 'particles_map.h'
!
  contains
!***********************************************************************
    subroutine initialize_particles_map()
!
!  Perform any post-parameter-read initialization.
!
!  29-mar-16/ccyang: coded.
!
!  Note: Currently, this subroutine is called after modules
!    Particles_mpicomm and Particles.
!
!  Check the particle-mesh interpolation method.
!
      pm: select case (particle_mesh)
!
      case ('ngp', 'NGP') pm
!       Nearest-Grid-Point
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .false.
        if (lroot) print *, 'initialize_particles_map: selected nearest-grid-point for particle-mesh method. '
!
      case ('cic', 'CIC') pm
!       Cloud-In-Cell
        lparticlemesh_cic = .true.
        lparticlemesh_tsc = .false.
        if (lroot) print *, 'initialize_particles_map: selected cloud-in-cell for particle-mesh method. '
!
      case ('tsc', 'TSC') pm
!       Triangular-Shaped-Cloud
        lparticlemesh_cic = .false.
        lparticlemesh_tsc = .true.
        if (lroot) print *, 'initialize_particles_map: selected triangular-shaped-cloud for particle-mesh method. '
!
      case ('') pm
!       Let the logical switches decide.
!       TSC assignment/interpolation overwrites CIC in case they are both set.
        switch: if (lparticlemesh_tsc) then
          lparticlemesh_cic = .false.
          particle_mesh = 'tsc'
        elseif (lparticlemesh_cic) then switch
          particle_mesh = 'cic'
        else switch
          particle_mesh = 'ngp'
        endif switch
        if (lroot) print *, 'initialize_particles_map: particle_mesh = ' // trim(particle_mesh)
!
      case default pm
        call fatal_error('initialize_particles_map', 'unknown particle-mesh type ' // trim(particle_mesh))
!
      endselect pm
!
    endsubroutine initialize_particles_map
!***********************************************************************
    subroutine map_nearest_grid(fp,ineargrid,k1_opt,k2_opt)
!
!  Find processor, brick, and grid point index of all or some of the
!  particles.
!
!  01-nov-09/anders: coded
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, optional :: k1_opt, k2_opt
!
      integer, dimension (0:nblockmax-1) :: ibrick_global_arr
      integer :: k1, k2, k, status
      integer :: iblockl, iblocku, iblockm, ibrick_global_par
      integer :: iproc_parent_par, ibrick_parent_par
      integer :: iproc_parent_par_previous, ibrick_parent_par_previous
      logical :: lbinary_search
!
      intent(in)  :: fp
      intent(out) :: ineargrid
!
      if (present(k1_opt)) then
        k1=k1_opt
      else
        k1=1
      endif
!
      if (present(k2_opt)) then
        k2=k2_opt
      else
        k2=npar_loc
      endif
!
!  Sort blocks by parent processor and by parent brick and create global
!  brick array.
!
      ibrick_global_arr(0:nblock_loc-1)= &
          iproc_parent_block(0:nblock_loc-1)*nbricks+ &
          ibrick_parent_block(0:nblock_loc-1)
!
      do k=k1,k2
!
!  Calculate processor, brick, and grid point index of particle.
!
        call get_brick_index(fp(k,(/ixp,iyp,izp/)), iproc_parent_par, &
            ibrick_parent_par, ineargrid(k,:), status=status)
        if (status < 0) then
          print*, 'map_nearest_grid: error in finding the grid index of particle'
          print*, 'map_nearest_grid: it, itsub, iproc, k, ipar=', it, itsub, iproc, k, ipar(k)
          call fatal_error_local('map_nearest_grid','')
        endif
!
!  Check if nearest block is the same as for previous particle.
!
        lbinary_search=.true.
        if (k>=k1+1) then
          if (iproc_parent_par==iproc_parent_par_previous .and. &
              ibrick_parent_par==ibrick_parent_par_previous) then
            inearblock(k)=inearblock(k-1)
            lbinary_search=.false.
          endif
        endif
!
!  Find nearest block by binary search.
!
        if (lbinary_search) then
          ibrick_global_par=iproc_parent_par*nbricks+ibrick_parent_par
          iblockl=0; iblocku=nblock_loc-1
          do while (abs(iblocku-iblockl)>1)
            iblockm=(iblockl+iblocku)/2
            if (ibrick_global_par>ibrick_global_arr(iblockm)) then
              iblockl=iblockm
            else
              iblocku=iblockm
            endif
          enddo
          if (ibrick_global_arr(iblockl)==ibrick_global_par) then
            inearblock(k)=iblockl
          elseif (ibrick_global_arr(iblocku)==ibrick_global_par) then
            inearblock(k)=iblocku
          else
            print*, 'map_nearest_grid: particle does not belong to any '// &
                'adopted block'
            print*, 'map_nearest_grid: it, itsub, iproc, k, ipar=', &
                it, itsub, iproc, k, ipar(k)
            print*, 'map_nearest_grid: ipx , ipy , ipz , iproc =', &
                ipx, ipy, ipz, iproc
            print*, 'map_nearest_grid: ibrick_global_par, iblockl, iblocku =', &
                ibrick_global_par, iblockl, iblocku
            print*, 'map_nearest_grid: ibrick_global_arr =', &
                ibrick_global_arr(0:nblock_loc-1)
            print*, 'map_nearest_grid: fp=', fp(k,:)
            print*, 'map_nearest_grid: xl=', xb(:,iblockl)
            print*, 'map_nearest_grid: yl=', yb(:,iblockl)
            print*, 'map_nearest_grid: zl=', zb(:,iblockl)
            print*, 'map_nearest_grid: xu=', xb(:,iblocku)
            print*, 'map_nearest_grid: yu=', yb(:,iblocku)
            print*, 'map_nearest_grid: zu=', zb(:,iblocku)
            call fatal_error_local('map_nearest_grid','')
          endif
        endif
        iproc_parent_par_previous=iproc_parent_par
        ibrick_parent_par_previous=ibrick_parent_par
      enddo
!
    endsubroutine map_nearest_grid
!***********************************************************************
    subroutine map_xxp_grid(f,fp,ineargrid,lmapsink_opt)
!
!  Map the particles as a continuous density field on the grid.
!
!  01-nov-09/anders: coded
!
      use GhostFold,     only: fold_f
      use Particles_sub, only: get_rhopswarm, weigh_particle
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
      logical, intent(in), optional :: lmapsink_opt
!
      real :: weight0, weight, weight_x, weight_y, weight_z
      real :: rhop_swarm_par
      integer :: k, ix0, iy0, iz0, ixx, iyy, izz, ib
      integer :: ixx0, ixx1, iyy0, iyy1, izz0, izz1, irhopm
      integer :: lb, mb, nb
      logical :: lsink, lmapsink
!
!  Possible to map sink particles by temporarily switching irhop to irhops.
!
      if (present(lmapsink_opt)) then
        lmapsink=lmapsink_opt
        if (lmapsink) then
          irhopm=irhop
          irhop=ipotself
        endif
      else
        lmapsink=.false.
      endif
!
!  Calculate the number of particles in each grid cell.
!
      if (inp/=0 .and. (mod(it,it1_loadbalance)==0.or.(.not.lnocalc_np)) &
           .and. (.not. lmapsink)) then
        f(:,:,:,inp)=0.0
        fb(:,:,:,inp,0:nblock_loc-1)=0.0
        do ib=0,nblock_loc-1
          if (npar_iblock(ib)/=0) then
            do k=k1_iblock(ib),k2_iblock(ib)
              ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
              fb(ix0,iy0,iz0,inp,ib)=fb(ix0,iy0,iz0,inp,ib)+1.0
            enddo
          endif
        enddo
      endif
!
!  Calculate the smooth number of particles in each grid cell. Three methods are
!  implemented for assigning a particle to the mesh (see Hockney & Eastwood):
!
!    0. NGP (Nearest Grid Point)
!       The entire effect of the particle goes to the nearest grid point.
!    1. CIC (Cloud In Cell)
!       The particle has a region of influence with the size of a grid cell.
!       This is equivalent to a first order (spline) interpolation scheme.
!    2. TSC (Triangular Shaped Cloud)
!       The particle is spread over a length of two grid cells, but with
!       a density that falls linearly outwards.
!       This is equivalent to a second order spline interpolation scheme.
!
      if (irhop/=0 .and. (.not.lnocalc_rhop)) then
        f(:,:,:,irhop)=0.0
        fb(:,:,:,irhop,0:nblock_loc-1)=0.0
        if (lparticlemesh_cic) then
!
!  Cloud In Cell (CIC) scheme.
!
          do ib=0,nblock_loc-1
            if (npar_iblock(ib)/=0) then
              do k=k1_iblock(ib),k2_iblock(ib)
                lsink=.false.
                if (lparticles_sink) then
                  if (fp(k,iaps)>0.0) lsink=.true.
                endif
                if (lmapsink) lsink=.not.lsink
                if (.not.lsink) then
                  ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
                  ixx0=ix0; iyy0=iy0; izz0=iz0
                  ixx1=ix0; iyy1=iy0; izz1=iz0
                  if ( (xb(ix0,ib)>fp(k,ixp)) .and. nxgrid/=1) ixx0=ixx0-1
                  if ( (yb(iy0,ib)>fp(k,iyp)) .and. nygrid/=1) iyy0=iyy0-1
                  if ( (zb(iz0,ib)>fp(k,izp)) .and. nzgrid/=1) izz0=izz0-1
                  if (nxgrid/=1) ixx1=ixx0+1
                  if (nygrid/=1) iyy1=iyy0+1
                  if (nzgrid/=1) izz1=izz0+1
!
!  Calculate mass density per superparticle.
!
                  if (lparticles_density) then
                    weight0=fp(k,irhopswarm)
                  elseif (lparticles_radius.and.lparticles_number) then
                    weight0=four_pi_rhopmat_over_three*fp(k,iap)**3* &
                        fp(k,inpswarm)
                  elseif (lparticles_radius) then
                    weight0=four_pi_rhopmat_over_three*fp(k,iap)**3*np_swarm
                  elseif (lparticles_number) then
                    weight0=mpmat*fp(k,inpswarm)
                  else
                    weight0=1.0
                  endif
!
                  do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
!
                    weight=weight0
!
                    if (nxgrid/=1) weight=weight* &
                        ( 1.0-abs(fp(k,ixp)-xb(ixx,ib))*dx1b(ixx,ib) )
                    if (nygrid/=1) weight=weight* &
                        ( 1.0-abs(fp(k,iyp)-yb(iyy,ib))*dy1b(iyy,ib) )
                    if (nzgrid/=1) weight=weight* &
                        ( 1.0-abs(fp(k,izp)-zb(izz,ib))*dz1b(izz,ib) )
!
                    fb(ixx,iyy,izz,irhop,ib)=fb(ixx,iyy,izz,irhop,ib) + weight
!
!  For debugging:
!
!                    if (weight<0.0 .or. weight>1.0) then
!                      print*, 'map_xxp_grid: weight is wrong'
!                      print*, 'map_xxp_grid: iproc, it, itsub, ipar=', &
!                          iproc, it, itsub, ipar(k)
!                      print*, 'map_xxp_grid: iblock, inearblock=', &
!                          ib, inearblock(k)
!                      print*, 'map_xxp_grid: weight=', weight
!                      print*, 'map_xxp_grid: xxp=', fp(k,ixp:izp)
!                      print*, 'map_xxp_grid: xb=', xb(:,ib)
!                      print*, 'map_xxp_grid: yb=', yb(:,ib)
!                      print*, 'map_xxp_grid: zb=', zb(:,ib)
!                    endif
                  enddo; enddo; enddo
                endif
              enddo
            endif
          enddo
        elseif (lparticlemesh_tsc) then
!
!  Triangular Shaped Cloud (TSC) scheme.
!
          do ib=0,nblock_loc-1
            if (npar_iblock(ib)/=0) then
              do k=k1_iblock(ib),k2_iblock(ib)
                lsink=.false.
                if (lparticles_sink) then
                  if (fp(k,iaps)>0.0) lsink=.true.
                endif
                if (lmapsink) lsink=.not.lsink
                if (.not.lsink) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
                  ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
                  call tsc_index_range(ix0, nxgrid, ixx0, ixx1)
                  call tsc_index_range(iy0, nygrid, iyy0, iyy1)
                  call tsc_index_range(iz0, nzgrid, izz0, izz1)
!
!  Calculate mass density per superparticle.
!
                  if (lparticles_density) then
                    weight0=fp(k,irhopswarm)
                  elseif (lparticles_radius.and.lparticles_number) then
                    weight0=four_pi_rhopmat_over_three*fp(k,iap)**3* &
                        fp(k,inpswarm)
                  elseif (lparticles_radius) then
                    weight0=four_pi_rhopmat_over_three*fp(k,iap)**3*np_swarm
                  elseif (lparticles_number) then
                    weight0=mpmat*fp(k,inpswarm)
                  else
                    weight0=1.0
                  endif
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
                  do izz=izz0,izz1 ; do iyy=iyy0,iyy1 ; do ixx=ixx0,ixx1
                    weight = weight0 * weigh_particle(abs(fp(k,ixp) - xb(ixx,ib)) * dx1b(ixx,ib), &
                                                      abs(fp(k,iyp) - yb(iyy,ib)) * dy1b(iyy,ib), &
                                                      abs(fp(k,izp) - zb(izz,ib)) * dz1b(izz,ib))
                    fb(ixx,iyy,izz,irhop,ib)=fb(ixx,iyy,izz,irhop,ib) + weight
!
!  For debugging:
!
!                    if (weight<0.0 .or. weight>1.0) then
!                      print*, 'map_xxp_grid: weight is wrong'
!                      print*, 'map_xxp_grid: iproc, it, itsub, ipar=', &
!                          iproc, it, itsub, ipar(k)
!                      print*, 'map_xxp_grid: iblock, inearblock=', &
!                          ib, inearblock(k)
!                      print*, 'map_xxp_grid: weight=', weight
!                      print*, 'map_xxp_grid: xxp=', fp(k,ixp:izp)
!                      print*, 'map_xxp_grid: xb=', xb(:,ib)
!                      print*, 'map_xxp_grid: yb=', yb(:,ib)
!                      print*, 'map_xxp_grid: zb=', zb(:,ib)
!                    endif
                  enddo; enddo; enddo
                endif
              enddo
            endif
          enddo
        else
!
!  Nearest Grid Point (NGP) method.
!
          if (lparticles_radius.or.lparticles_number.or.lparticles_density) then
            if (npar_iblock(ib)/=0) then
              do k=k1_iblock(ib),k2_iblock(ib)
                lsink=.false.
                if (lparticles_sink) then
                  if (fp(k,iaps)>0.0) lsink=.true.
                endif
                if (lmapsink) lsink=.not.lsink
                if (.not.lsink) then
                  ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
!
!  Calculate mass density per superparticle.
!
                  if (lparticles_density) then
                    weight0=fp(k,irhopswarm)
                  elseif (lparticles_radius.and.lparticles_number) then
                    weight0=four_pi_rhopmat_over_three*fp(k,iap)**3* &
                        fp(k,inpswarm)
                  elseif (lparticles_radius) then
                    weight0=four_pi_rhopmat_over_three*fp(k,iap)**3*np_swarm
                  elseif (lparticles_number) then
                    weight0=mpmat*fp(k,inpswarm)
                  else
                    weight0=1.0
                  endif
!
                  f(ix0,iy0,iz0,irhop)=f(ix0,iy0,iz0,irhop) + weight0
!
                endif
              enddo
            endif
          else
            fb(:,:,:,irhop,0:nblock_loc-1)=fb(:,:,:,inp,0:nblock_loc-1)
          endif
        endif
!
!  Multiply assigned particle number density by the mass density per particle.
!
        if (.not.(lparticles_radius.or.lparticles_number.or. &
            lparticles_density)) then
          if (lcartesian_coords.and.(all(lequidist))) then
            fb(:,:,:,irhop,0:nblock_loc-1)= &
                 rhop_swarm*fb(:,:,:,irhop,0:nblock_loc-1)
          else
!
            do ib=0,nblock_loc-1
              do nb=1,mzb ; do mb=1,myb ; do lb=1,mxb
                rhop_swarm_par  = mp_swarm*&
                     dVol1xb(lb,ib)*dVol1yb(mb,ib)*dVol1zb(nb,ib)
                fb(lb,mb,nb,irhop,ib) = rhop_swarm_par*fb(lb,mb,nb,irhop,ib)
              enddo;enddo;enddo
            enddo
!
          endif
        endif
!
!  Fill the bricks on each processor with particle density assigned on the
!  blocks.
!
        if (inp/=0 .and. (mod(it,it1_loadbalance)==0.or.(.not.lnocalc_np))) &
            call fill_bricks_with_blocks(f,fb,mfarray,inp)
        if (irhop/=0 .and. (.not.lnocalc_rhop)) &
            call fill_bricks_with_blocks(f,fb,mfarray,irhop)
!
!  Fold first ghost zone of f.
!
        if (lparticlemesh_cic.or.lparticlemesh_tsc) call fold_f(f,irhop,irhop)
!
!  Give folded rhop back to the blocks on the foster parents.
!
        if (lparticlemesh_cic.or.lparticlemesh_tsc) &
             call fill_blocks_with_bricks(f,fb,mfarray,irhop)
!
      else
!
!  Only particle number density.
!
        call fill_bricks_with_blocks(f,fb,mfarray,inp)
      endif
!
!  Restore normal particle index if mapping sink
!
      if (lmapsink) irhop=irhopm
!
    endsubroutine map_xxp_grid
!***********************************************************************
    subroutine map_vvp_grid(f,fp,ineargrid)
!
!  Map the particle velocities as vector field on the grid.
!
!  16-nov-09/anders: dummy
!  17-may-23/ccyang: under construction
!
      use GhostFold, only: fold_f
      use Particles_sub, only: weigh_particle
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mpar_loc,mparray), intent(in) :: fp
      integer, dimension(mpar_loc,3), intent(in) :: ineargrid
!
      integer :: ix0, ixx0, ixx1
      integer :: iy0, iyy0, iyy1
      integer :: iz0, izz0, izz1
      integer :: ixx, iyy, izz
      integer :: ib, ivp, k
      real :: weight, dxi1, dxi2, dxi3
!
      uup: if (iuup /= 0) then
        f(:,:,:,iupx:iupz) = 0.0
        fb(:,:,:,iupx:iupz,0:nblock_loc-1) = 0.0
!
!  Sanity check.
!
        if (lparticles_density .or. lparticles_radius .or. lparticles_number) &
            call fatal_error("map_vvp_grid", "variable particle mass is not implemented")
!
        blocks: do ib = 0, nblock_loc - 1
          tsc: if (npar_iblock(ib) /= 0) then
            par: do k = k1_iblock(ib), k2_iblock(ib)
!
!  Find neighboring grid points of a particle.
!
              ix0 = ineargrid(k,1)
              iy0 = ineargrid(k,2)
              iz0 = ineargrid(k,3)
              pm: if (lparticlemesh_tsc) then
                call tsc_index_range(ix0, nxgrid, ixx0, ixx1)
                call tsc_index_range(iy0, nygrid, iyy0, iyy1)
                call tsc_index_range(iz0, nzgrid, izz0, izz1)
              else pm
                call fatal_error("map_vvp_grid", "not implemented for non-TSC scheme")
              endif pm
!
!  Add the particle momentum to neighboring grid points.
!
              zscan: do izz = izz0, izz1
                dxi3 = abs(fp(k,izp) - zb(izz,ib)) * dz1b(izz,ib)
                yscan: do iyy = iyy0, iyy1
                  dxi2 = abs(fp(k,iyp) - yb(iyy,ib)) * dy1b(iyy,ib)
                  xscan: do ixx = ixx0, ixx1
                    dxi1 = abs(fp(k,ixp) - xb(ixx,ib)) * dx1b(ixx,ib)
                    weight = mp_swarm * weigh_particle(dxi1, dxi2, dxi3)
                    fb(ixx,iyy,izz,iupx:iupz,ib) = fb(ixx,iyy,izz,iupx:iupz,ib) + weight * fp(k,ivpx:ivpz)
                  enddo xscan
                enddo yscan
              enddo zscan
!
            enddo par
          endif tsc
        enddo blocks
!
!  Communicate the particle momentum.
!
        do ivp = iupx, iupz
          call fill_bricks_with_blocks(f, fb, mfarray, ivp)
        enddo
!
        comm: if (particle_mesh /= "NGP" .and. particle_mesh /= "ngp") then
          call fold_f(f, iupx, iupz)
          uup1: do ivp = iupx, iupz
            where (f(l1:l2,m1:m2,n1:n2,irhop) > 0.0)
              f(l1:l2,m1:m2,n1:n2,ivp) = f(l1:l2,m1:m2,n1:n2,ivp) / f(l1:l2,m1:m2,n1:n2,irhop)
            elsewhere
              f(l1:l2,m1:m2,n1:n2,ivp) = 0.0
            endwhere
            do izz = n1, n2
              do iyy = m1, m2
                f(l1:l2,iyy,izz,ivp) = f(l1:l2,iyy,izz,ivp) / (dVol_x(l1:l2) * dVol_y(iyy) * dVol_z(izz))
              enddo
            enddo
            call fill_blocks_with_bricks(f, fb, mfarray, ivp)
          enddo uup1
        endif comm
      endif uup
!
    endsubroutine map_vvp_grid
!***********************************************************************
    subroutine sort_particles_iblock(fp,ineargrid,ipar,dfp)
!
!  Sort the particles so that they appear in order of the global brick index.
!  That is, sorted first by processor number and then by local brick index.
!
!  12-oct-09/anders: coded
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
!
      integer, dimension (mpar_loc) :: ipark_sorted
      integer, dimension (0:nblockmax-1) :: kk
      integer :: k
!
      intent(inout) :: fp, ineargrid, dfp, ipar
!
!  Determine beginning and ending index of particles from each block.
!
      call particle_block_index()
!
!  Sort particles by blocks (counting sort).
!
      kk=k1_iblock
      do k=1,npar_loc
        ipark_sorted(kk(inearblock(k)))=k
        kk(inearblock(k))=kk(inearblock(k))+1
      enddo
!
      if (npar_loc>0) then
        ineargrid(1:npar_loc,:)=ineargrid(ipark_sorted(1:npar_loc),:)
        inearblock(1:npar_loc)=inearblock(ipark_sorted(1:npar_loc))
        ipar(1:npar_loc)=ipar(ipark_sorted(1:npar_loc))
        fp(1:npar_loc,:)=fp(ipark_sorted(1:npar_loc),:)
        if (present(dfp)) dfp(1:npar_loc,:)=dfp(ipark_sorted(1:npar_loc),:)
      endif
!
!  Possible to randomize particles inside each block. This screws with the
!  pencil consistency check, so we turn it off when the test is running.
!
      if (lrandom_particle_blocks .and. (.not.lpencil_check_at_work)) then
        if (present(dfp)) then
          call random_particle_blocks(fp,ineargrid,inearblock,ipar,dfp)
        else
          call random_particle_blocks(fp,ineargrid,inearblock,ipar)
        endif
      endif
!
    endsubroutine sort_particles_iblock
!***********************************************************************
    subroutine random_particle_blocks(fp,ineargrid,inearblock,ipar,dfp)
!
!  Randomize particles within each block to avoid low index particles
!  always being considered first.
!
!  Slows down simulation by around 10%.
!
!  31-jan-10/anders: coded
!
      use General, only: random_number_wrapper
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (mpar_loc) :: inearblock, ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
!
      real, dimension (mparray) :: fp_swap
      real, dimension (mpvar) :: dfp_swap
      real :: r
      integer, dimension (3) :: ineargrid_swap
      integer :: inearblock_swap, ipar_swap, iblock, k, kswap
!
      intent (out) :: fp, ineargrid, ipar, dfp
!
      do iblock=0,nblock_loc-1
        if (npar_iblock(iblock)>=2) then
          do k=k1_iblock(iblock),k2_iblock(iblock)
            call random_number_wrapper(r)
            kswap=k1_iblock(iblock)+floor(r*npar_iblock(iblock))
            if (kswap/=k) then
              fp_swap=fp(kswap,:)
              ineargrid_swap=ineargrid(kswap,:)
              inearblock_swap=inearblock(kswap)
              ipar_swap=ipar(kswap)
              if (present(dfp)) dfp_swap=dfp(kswap,:)
              fp(kswap,:)=fp(k,:)
              ineargrid(kswap,:)=ineargrid(k,:)
              inearblock(kswap)=inearblock(k)
              ipar(kswap)=ipar(k)
              if (present(dfp)) dfp(kswap,:)=dfp(k,:)
              fp(k,:)=fp_swap
              ineargrid(k,:)=ineargrid_swap
              inearblock(k)=inearblock_swap
              ipar(k)=ipar_swap
              if (present(dfp)) dfp(k,:)=dfp_swap
            endif
          enddo
        endif
      enddo
!
    endsubroutine random_particle_blocks
!***********************************************************************
    subroutine particle_block_index()
!
!  Calculate the beginning and ending index of particles in each block.
!
!  18-nov-09/anders: coded
!
      integer :: k, iblock
!
      npar_iblock=0
!
!  Calculate the number of particles adopted from each block.
!
      do k=1,npar_loc
        npar_iblock(inearblock(k))=npar_iblock(inearblock(k))+1
      enddo
!
!  Calculate beginning and ending particle index for each block.
!
      k=0
      do iblock=0,nblock_loc-1
        if (npar_iblock(iblock)/=0) then
          k1_iblock(iblock)=k+1
          k2_iblock(iblock)=k1_iblock(iblock)+npar_iblock(iblock)-1
          k=k+npar_iblock(iblock)
        else
          k1_iblock(iblock)=0
          k2_iblock(iblock)=0
        endif
      enddo
!
    endsubroutine particle_block_index
!***********************************************************************
    subroutine fill_blocks_with_bricks(a,ab,marray,ivar)
!
!  Fill adopted blocks with bricks from the f-array.
!
!  04-nov-09/anders: coded
!
      use Mpicomm, only: mpirecv_nonblock_real,mpisend_nonblock_real,mpiwait
!
      integer :: marray, ivar
      real, dimension (mx,my,mz,marray) :: a
      real, dimension (mxb,myb,mzb,marray,0:nblockmax-1) :: ab
!
      integer, dimension (2*nblockmax) :: ireq_array
      real, dimension (mxb,myb,mzb,0:nbricks-1) :: ab_send
      real, dimension (mxb,myb,mzb,0:nblockmax-1) :: ab_recv
      integer :: ibx, iby, ibz, ix1, ix2, iy1, iy2, iz1, iz2
      integer :: ireq, nreq, iblock, iblock1, iblock2
      integer :: ibrick, ibrick1, ibrick2
      integer :: ibrick_send, ibrick1_send, ibrick2_send
      integer :: iproc_recv, iproc_send, nblock_recv, nbrick_send, tag_id
!
      intent (in) :: a, ivar
      intent (out) :: ab
!
      tag_id=100
      nreq=0
!
!  Fill up blocks with bricks from local processor.
!
      do iblock=0,nblock_loc-1
        if (iproc_parent_block(iblock)==iproc) then
          ibrick=ibrick_parent_block(iblock)
          ibx=modulo(ibrick,nbx)
          iby=modulo(ibrick/nbx,nby)
          ibz=ibrick/(nbx*nby)
          ix1=l1-1+ibx*nxb; ix2=l1+(ibx+1)*nxb
          iy1=m1-1+iby*nyb; iy2=m1+(iby+1)*nyb
          iz1=n1-1+ibz*nzb; iz2=n1+(ibz+1)*nzb
          ab(:,:,:,ivar,iblock)=a(ix1:ix2,iy1:iy2,iz1:iz2,ivar)
        endif
      enddo
!
!  Fill up buffer array with bricks that need to be send to other processors.
!
      ibrick_send=0
      do ibrick=0,nbricks-1
        if (iproc_foster_brick(ibrick)/=iproc .and. &
            iproc_foster_brick(ibrick)/=-1) then
          ibx=modulo(ibrick,nbx)
          iby=modulo(ibrick/nbx,nby)
          ibz=ibrick/(nbx*nby)
          ix1=l1-1+ibx*nxb; ix2=l1+(ibx+1)*nxb
          iy1=m1-1+iby*nyb; iy2=m1+(iby+1)*nyb
          iz1=n1-1+ibz*nzb; iz2=n1+(ibz+1)*nzb
          ab_send(:,:,:,ibrick_send)=a(ix1:ix2,iy1:iy2,iz1:iz2,ivar)
          ibrick_send=ibrick_send+1
        endif
      enddo
!
!  Receive blocks with non-blocking MPI.
!
      iblock=0
      do while (iblock<nblock_loc)
        iproc_recv=iproc_parent_block(iblock)
        iblock1=iblock
        iblock2=iblock
        do while (iblock2<nblock_loc-1)
          if (iproc_parent_block(iblock2+1)==iproc_recv) then
            iblock2=iblock2+1
          else
            exit
          endif
        enddo
        if (iproc_parent_block(iblock)/=iproc) then
          nblock_recv=iblock2-iblock1+1
          call mpirecv_nonblock_real(ab_recv(:,:,:,iblock1:iblock2),&
               (/mxb,myb,mzb,nblock_recv/),&
               iproc_recv,&
               tag_id+iproc_recv,&
               ireq)
          nreq=nreq+1
          ireq_array(nreq)=ireq
        endif
        iblock=iblock2+1
      enddo
!
!  Send bricks with non-blocking MPI.
!
      ibrick=0
      ibrick1_send=0
      ibrick2_send=0
      do while (ibrick<nbricks)
        iproc_send=iproc_foster_brick(ibrick)
        ibrick1=ibrick
        ibrick2=ibrick
        do while (ibrick2<nbricks-1)
          if (iproc_foster_brick(ibrick2+1)==iproc_send) then
            ibrick2=ibrick2+1
            if (iproc_foster_brick(ibrick1)/=iproc .and. &
                iproc_foster_brick(ibrick1)/=-1) ibrick2_send=ibrick2_send+1
          else
            if (iproc_foster_brick(ibrick2+1)==iproc .or. &
                iproc_foster_brick(ibrick2+1)==-1) then
              ibrick2=ibrick2+1
            else
              exit
            endif
          endif
        enddo
        if (iproc_foster_brick(ibrick)/=iproc .and. &
            iproc_foster_brick(ibrick)/=-1) then
          nbrick_send=ibrick2_send-ibrick1_send+1
!          
          call mpisend_nonblock_real(ab_send(:,:,:,ibrick1_send:ibrick2_send),&
               (/mxb,myb,mzb,nbrick_send/),&
               iproc_send,&
               tag_id+iproc,&
               ireq)
          nreq=nreq+1
          ireq_array(nreq)=ireq
          ibrick1_send=ibrick2_send+1
          ibrick2_send=ibrick2_send+1
        endif
        ibrick=ibrick2+1
      enddo
!
!  Wait for non-blocking MPI calls to finish.
!
      if (nreq>0) then
        do ireq=1,nreq
          call mpiwait(ireq_array(ireq))
        enddo
      endif
!
!  Fill up blocks with bricks from receive buffer.
!
      iblock=0
      do while (iblock<nblock_loc)
        iproc_recv=iproc_parent_block(iblock)
        iblock1=iblock
        iblock2=iblock
        do while (iblock2<nblock_loc-1)
          if (iproc_parent_block(iblock2+1)==iproc_recv) then
            iblock2=iblock2+1
          else
            exit
          endif
        enddo
        if (iproc_parent_block(iblock)/=iproc) then
          ab(:,:,:,ivar,iblock1:iblock2)=ab_recv(:,:,:,iblock1:iblock2)
        endif
        iblock=iblock2+1
      enddo
!
    endsubroutine fill_blocks_with_bricks
!***********************************************************************
    subroutine fill_bricks_with_blocks(a,ab,marray,ivar,nosum_opt)
!
!  Fill bricks (i.e. the f-array) with blocks adopted by other processors.
!
!  04-nov-09/anders: coded
!
      use Mpicomm, only: mpirecv_nonblock_real,mpisend_nonblock_real,mpiwait
!
      integer :: marray, ivar
      real, dimension (mx,my,mz,marray) :: a
      real, dimension (mxb,myb,mzb,marray,0:nblockmax-1) :: ab
      logical, optional :: nosum_opt
!
      integer, dimension (2*nblockmax) :: ireq_array
      real, dimension (mxb,myb,mzb,0:nblockmax-1) :: ab_send
      real, dimension (mxb,myb,mzb,0:nbricks-1) :: ab_recv
      integer :: ibx, iby, ibz, ix1, ix2, iy1, iy2, iz1, iz2
      integer :: ireq, nreq, iblock, iblock1, iblock2
      integer :: ibrick, ibrick1, ibrick2
      integer :: ibrick_recv, ibrick1_recv, ibrick2_recv
      integer :: iblock_send, iblock1_send, iblock2_send
      integer :: iproc_recv, iproc_send, nbrick_recv, nblock_send, tag_id
      logical :: nosum
!
      intent (in) :: ab, ivar
      intent (out) :: a
!
      if (present(nosum_opt)) then
        nosum=nosum_opt
      else
        nosum=.false.
      endif
!
      tag_id=1000
      nreq=0
!
!  Fill up bricks with blocks from local processor.
!
      do iblock=0,nblock_loc-1
        if (iproc_parent_block(iblock)==iproc) then
          ibrick=ibrick_parent_block(iblock)
          ibx=modulo(ibrick,nbx)
          iby=modulo(ibrick/nbx,nby)
          ibz=ibrick/(nbx*nby)
          ix1=l1-1+ibx*nxb; ix2=l1+(ibx+1)*nxb
          iy1=m1-1+iby*nyb; iy2=m1+(iby+1)*nyb
          iz1=n1-1+ibz*nzb; iz2=n1+(ibz+1)*nzb
          if (nosum) then
            a(ix1:ix2,iy1:iy2,iz1:iz2,ivar)=ab(:,:,:,ivar,iblock)
          else
            a(ix1:ix2,iy1:iy2,iz1:iz2,ivar)= &
                a(ix1:ix2,iy1:iy2,iz1:iz2,ivar)+ab(:,:,:,ivar,iblock)
          endif
        endif
      enddo
!
!  Fill up buffer array with blocks that need to be send to other processors.
!
      iblock=0
      iblock_send=0
      do while (iblock<nblock_loc)
        if (iproc_parent_block(iblock)/=iproc .and. &
            iproc_parent_block(iblock)/=-1) then
          ab_send(:,:,:,iblock_send)=ab(:,:,:,ivar,iblock)
          iblock_send=iblock_send+1
        endif
        iblock=iblock+1
      enddo
!
!  Receive bricks with non-blocking MPI.
!
      ibrick=0
      ibrick1_recv=0
      ibrick2_recv=0
      do while (ibrick<nbricks)
        iproc_recv=iproc_foster_brick(ibrick)
        ibrick1=ibrick
        ibrick2=ibrick
        do while (ibrick2<nbricks-1)
          if (iproc_foster_brick(ibrick2+1)==iproc_recv) then
            ibrick2=ibrick2+1
            if (iproc_foster_brick(ibrick1)/=iproc .and. &
                iproc_foster_brick(ibrick1)/=-1) ibrick2_recv=ibrick2_recv+1
          else
            if (iproc_foster_brick(ibrick2+1)==-1 .or. &
                iproc_foster_brick(ibrick2+1)==iproc) then
              ibrick2=ibrick2+1
            else
              exit
            endif
          endif
        enddo
        if (iproc_foster_brick(ibrick)/=iproc .and. &
            iproc_foster_brick(ibrick)/=-1) then
          nbrick_recv=ibrick2_recv-ibrick1_recv+1
          call mpirecv_nonblock_real(ab_recv(:,:,:,ibrick1_recv:ibrick2_recv),&
               (/mxb,myb,mzb,nbrick_recv/),&
               iproc_recv,&
               tag_id+iproc_recv,&
               ireq)
          nreq=nreq+1
          ireq_array(nreq)=ireq
          ibrick1_recv=ibrick2_recv+1
          ibrick2_recv=ibrick2_recv+1
        endif
        ibrick=ibrick2+1
      enddo
!
!  Send blocks with non-blocking MPI.
!
      iblock=0
      iblock1_send=0
      iblock2_send=0
      do while (iblock<nblock_loc)
        iproc_send=iproc_parent_block(iblock)
        iblock1=iblock
        iblock2=iblock
        do while (iblock2<nblock_loc-1)
          if (iproc_parent_block(iblock2+1)==iproc_send) then
            iblock2=iblock2+1
            if (iproc_parent_block(iblock1)/=iproc .and. &
                iproc_parent_block(iblock1)/=-1) iblock2_send=iblock2_send+1
          else
            if (iproc_parent_block(iblock2+1)==-1) then
              iblock2=iblock2+1
            else
              exit
            endif
          endif
        enddo
        if (iproc_parent_block(iblock)/=iproc .and. &
            iproc_parent_block(iblock)/=-1) then
          nblock_send=iblock2_send-iblock1_send+1
          call mpisend_nonblock_real(ab_send(:,:,:,iblock1_send:iblock2_send),&
               (/mxb,myb,mzb,nblock_send/),&
               iproc_send,&
               tag_id+iproc,&
               ireq)
          nreq=nreq+1
          ireq_array(nreq)=ireq
          iblock1_send=iblock2_send+1
          iblock2_send=iblock2_send+1
        endif
        iblock=iblock2+1
      enddo
!
!  Wait for non-blocking MPI calls to finish.
!
      if (nreq>0) then
        do ireq=1,nreq
          call mpiwait(ireq_array(ireq))
        enddo
      endif
!
!  Fill up bricks with blocks from receive buffer.
!
      ibrick_recv=0
      do ibrick=0,nbricks-1
        if (iproc_foster_brick(ibrick)/=iproc .and. &
            iproc_foster_brick(ibrick)/=-1) then
          ibx=modulo(ibrick,nbx)
          iby=modulo(ibrick/nbx,nby)
          ibz=ibrick/(nbx*nby)
          ix1=l1-1+ibx*nxb; ix2=l1+(ibx+1)*nxb
          iy1=m1-1+iby*nyb; iy2=m1+(iby+1)*nyb
          iz1=n1-1+ibz*nzb; iz2=n1+(ibz+1)*nzb
          if (nosum) then
            a(ix1:ix2,iy1:iy2,iz1:iz2,ivar)=ab_recv(:,:,:,ibrick_recv)
          else
            a(ix1:ix2,iy1:iy2,iz1:iz2,ivar)= &
                a(ix1:ix2,iy1:iy2,iz1:iz2,ivar)+ab_recv(:,:,:,ibrick_recv)
          endif
          ibrick_recv=ibrick_recv+1
        endif
      enddo
!
    endsubroutine fill_bricks_with_blocks
!***********************************************************************
    subroutine interpolate_linear(f,ivar1,ivar2,xxp,gp,inear,iblock,ipar)
!
!  Interpolate the value of g to arbitrary (xp, yp, zp) coordinate
!  using the linear interpolation formula
!
!    g(x,y,z) = A*x*y*z + B*x*y + C*x*z + D*y*z + E*x + F*y + G*z + H .
!
!  The coefficients are determined by the 8 grid points surrounding the
!  interpolation point.
!
!  30-dec-04/anders: coded
!
      use Solid_Cells
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (3) :: inear
      integer :: iblock, ipar
!
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8
      real :: xp0, yp0, zp0
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      integer :: i, ix0, iy0, iz0, ib
      logical :: lfirstcall=.true.
!
      intent(in)  :: xxp, ivar1, inear
      intent(out)  :: gp
!
      call keep_compiler_quiet(f)
!
!  Abbreviations.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      ib=iblock
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      if ( (xb(ix0,iblock)>xxp(1)) .and. nxgrid/=1) ix0=ix0-1
      if ( (yb(iy0,iblock)>xxp(2)) .and. nygrid/=1) iy0=iy0-1
      if ( (zb(iz0,iblock)>xxp(3)) .and. nzgrid/=1) iz0=iz0-1
!
!  Check if the grid point interval is really correct.
!
      if ((xb(ix0,ib)<=xxp(1) .and. xb(ix0+1,ib)>=xxp(1) .or. nxgrid==1) .and. &
          (yb(iy0,ib)<=xxp(2) .and. yb(iy0+1,ib)>=xxp(2) .or. nygrid==1) .and. &
          (zb(iz0,ib)<=xxp(3) .and. zb(iz0+1,ib)>=xxp(3) .or. nzgrid==1)) then
        ! Everything okay
      else
        print*, 'interpolate_linear: Interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'ipar = ', ipar
        print*, 'mxb, xb(1), xb(mx) = ', mxb, xb(1,ib), xb(mxb,ib)
        print*, 'myb, yb(1), yb(my) = ', myb, yb(1,ib), yb(myb,ib)
        print*, 'mzb, zb(1), zb(mz) = ', mzb, zb(1,ib), zb(mzb,ib)
        print*, 'iblock, ix0, iy0, iz0 = ', iblock, ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), xb(ix0,ib), xb(ix0+1,ib)
        print*, 'yp, yp0, yp1 = ', xxp(2), yb(iy0,ib), yb(iy0+1,ib)
        print*, 'zp, zp0, zp1 = ', xxp(3), zb(iz0,ib), zb(iz0+1,ib)
        call fatal_error('interpolate_linear','')
      endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!  Set it equal to 0 for dimensions having 1 grid points; this will make sure
!  that the interpolation is bilinear for 2D grids.
!
      xp0=0; yp0=0; zp0=0
      if (nxgrid/=1) xp0=xxp(1)-xb(ix0,ib)
      if (nygrid/=1) yp0=xxp(2)-yb(iy0,ib)
      if (nzgrid/=1) zp0=xxp(3)-zb(iz0,ib)
!
!  Calculate derived grid spacing parameters needed for interpolation.
!  For an equidistant grid we only need to do this at the first call.
!
      if (lequidist(1)) then
        if (lfirstcall) dx1=dx1b(ix0,ib) !1/dx
      else
        dx1=dx1b(ix0,ib)
      endif
!
      if (lequidist(2)) then
        if (lfirstcall) dy1=dy1b(iy0,ib)
      else
        dy1=dy1b(iy0,ib)
      endif
!
      if (lequidist(3)) then
        if (lfirstcall) dz1=dz1b(iz0,ib)
      else
        dz1=dz1b(iz0,ib)
      endif
!
      if ( (.not. all(lequidist)) .or. lfirstcall) then
        dxdy1=dx1*dy1; dxdz1=dx1*dz1; dydz1=dy1*dz1
        dxdydz1=dx1*dy1*dz1
      endif
!
!  Function values at all corners.
!
      g1=fb(ix0  ,iy0  ,iz0  ,ivar1:ivar2,ib)
      g2=fb(ix0+1,iy0  ,iz0  ,ivar1:ivar2,ib)
      g3=fb(ix0  ,iy0+1,iz0  ,ivar1:ivar2,ib)
      g4=fb(ix0+1,iy0+1,iz0  ,ivar1:ivar2,ib)
      g5=fb(ix0  ,iy0  ,iz0+1,ivar1:ivar2,ib)
      g6=fb(ix0+1,iy0  ,iz0+1,ivar1:ivar2,ib)
      g7=fb(ix0  ,iy0+1,iz0+1,ivar1:ivar2,ib)
      g8=fb(ix0+1,iy0+1,iz0+1,ivar1:ivar2,ib)
!
!  Interpolation formula.
!
      gp = g1 + xp0*dx1*(-g1+g2) + yp0*dy1*(-g1+g3) + zp0*dz1*(-g1+g5) + &
          xp0*yp0*dxdy1*(g1-g2-g3+g4) + xp0*zp0*dxdz1*(g1-g2-g5+g6) + &
          yp0*zp0*dydz1*(g1-g3-g5+g7) + &
          xp0*yp0*zp0*dxdydz1*(-g1+g2+g3-g4+g5-g6-g7+g8)
!
!  Do a reality check on the interpolation scheme.
!
      if (linterp_reality_check) then
        do i=1,ivar2-ivar1+1
          if (gp(i)>max(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
            print*, 'interpolate_linear: interpolated value is LARGER than'
            print*, 'interpolate_linear: a values at the corner points!'
            print*, 'interpolate_linear: ipar, xxp=', ipar, xxp
            print*, 'interpolate_linear: x0, y0, z0=', &
                xb(ix0,ib), yb(iy0,ib), zb(iz0,ib)
            print*, 'interpolate_linear: i, gp(i)=', i, gp(i)
            print*, 'interpolate_linear: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
          if (gp(i)<min(g1(i),g2(i),g3(i),g4(i),g5(i),g6(i),g7(i),g8(i))) then
            print*, 'interpolate_linear: interpolated value is smaller than'
            print*, 'interpolate_linear: a values at the corner points!'
            print*, 'interpolate_linear: xxp=', xxp
            print*, 'interpolate_linear: x0, y0, z0=', &
                xb(ix0,ib), yb(iy0,ib), zb(iz0,ib)
            print*, 'interpolate_linear: i, gp(i)=', i, gp(i)
            print*, 'interpolate_linear: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
        enddo
      endif
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine interpolate_linear
!***********************************************************************
    subroutine interpolate_quadratic(f,ivar1,ivar2,xxp,gp,inear,iblock,ipar)
!
!  Quadratic interpolation of g to arbitrary (xp, yp, zp) coordinate
!  using the biquadratic interpolation function
!
!    g(x,y,z) = (1+x+x^2)*(1+y+y^2)*(1+z+z^2)
!
!  The coefficients (9, one for each unique term) are determined by the 9
!  grid points surrounding the interpolation point.
!
!  The interpolation matrix M is defined through the relation
!    M#c = g
!  Here c are the coefficients and g is the value of the function at the grid
!  points. An equidistant grid has the following value of M:
!
!    invmat(:,1)=(/ 0.00, 0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00, 0.00/)
!    invmat(:,2)=(/ 0.00, 0.00, 0.00,-0.50, 0.00, 0.50, 0.00, 0.00, 0.00/)
!    invmat(:,3)=(/ 0.00, 0.00, 0.00, 0.50,-1.00, 0.50, 0.00, 0.00, 0.00/)
!    invmat(:,4)=(/ 0.00,-0.50, 0.00, 0.00, 0.00, 0.00, 0.00, 0.50, 0.00/)
!    invmat(:,5)=(/ 0.00, 0.50, 0.00, 0.00,-1.00, 0.00, 0.00, 0.50, 0.00/)
!    invmat(:,6)=(/ 0.25, 0.00,-0.25, 0.00, 0.00, 0.00,-0.25, 0.00, 0.25/)
!    invmat(:,7)=(/-0.25, 0.50,-0.25, 0.00, 0.00, 0.00, 0.25,-0.50, 0.25/)
!    invmat(:,8)=(/-0.25, 0.00, 0.25, 0.50, 0.00,-0.50,-0.25, 0.00, 0.25/)
!    invmat(:,9)=(/ 0.25,-0.50, 0.25,-0.50, 1.00,-0.50, 0.25,-0.50, 0.25/)
!
!    invmat(:,1)=invmat(:,1)
!    invmat(:,2)=invmat(:,2)/dx
!    invmat(:,3)=invmat(:,3)/dx**2
!    invmat(:,4)=invmat(:,4)/dz
!    invmat(:,5)=invmat(:,5)/dz**2
!    invmat(:,6)=invmat(:,6)/(dx*dz)
!    invmat(:,7)=invmat(:,7)/(dx**2*dz)
!    invmat(:,8)=invmat(:,8)/(dx*dz**2)
!    invmat(:,9)=invmat(:,9)/(dx**2*dz**2)
!
!  Space coordinates are defined such that the nearest grid point is at (0,0).
!  The grid points are counted from lower left:
!
!    7  8  9
!    4  5  6
!    1  2  3
!
!  The nearest grid point has index number 5.
!
!  09-jun-06/anders: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (3) :: inear
      integer :: iblock, ipar
!
      real, dimension (9,ivar2-ivar1+1) :: cc
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8, g9
      real :: dxp, dzp
      real, save :: dx1, dx2, dz1, dz2
      real, save :: dx1dz1, dx2dz1, dx1dz2, dx2dz2
      integer :: ix0, iy0, iz0, ib
      logical, save :: lfirstcall=.true.
!
      intent(in)  :: xxp, ivar1, inear
      intent(out) :: gp
!
      call keep_compiler_quiet(f)
!
!  Abbreviations.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      ib=iblock
!
!  Not implemented in y-direction yet (but is easy to generalise).
!
      if (nygrid/=1) then
        if (lroot) print*, 'interpolate_quadratic: not implemented in y'
        call fatal_error('interpolate_quadratic','')
      endif
!
!  A few values that only need to be calculated once for equidistant grids.
!
      if (lequidist(1)) then
        if (lfirstcall) then
          dx1=1/dx; dx2=1/dx**2
        endif
      else
        dx1=dx1b(ix0,ib); dx2=dx1**2
      endif
!
      if (lequidist(3)) then
        if (lfirstcall) then
          dz1=1/dz; dz2=1/dz**2
        endif
      else
        dz1=dz1b(iz0,ib); dz2=dz1**2
      endif
!
      if (lequidist(1).and.lequidist(3)) then
        if (lfirstcall) then
          dx1dz1=1/(dx*dz)
          dx2dz1=1/(dx**2*dz); dx1dz2=1/(dx*dz**2); dx2dz2=1/(dx**2*dz**2)
        endif
      else
        dx1dz1=dx1*dz1
        dx2dz1=dx2*dz1; dx1dz2=dx1*dz1; dx2dz2=dx2*dz2
      endif
!
!  Define function values at the grid points.
!
      g1=fb(ix0-1,iy0,iz0-1,ivar1:ivar2,ib)
      g2=fb(ix0  ,iy0,iz0-1,ivar1:ivar2,ib)
      g3=fb(ix0+1,iy0,iz0-1,ivar1:ivar2,ib)
      g4=fb(ix0-1,iy0,iz0  ,ivar1:ivar2,ib)
      g5=fb(ix0  ,iy0,iz0  ,ivar1:ivar2,ib)
      g6=fb(ix0+1,iy0,iz0  ,ivar1:ivar2,ib)
      g7=fb(ix0-1,iy0,iz0+1,ivar1:ivar2,ib)
      g8=fb(ix0  ,iy0,iz0+1,ivar1:ivar2,ib)
      g9=fb(ix0+1,iy0,iz0+1,ivar1:ivar2,ib)
!
!  Calculate the coefficients of the interpolation formula (see introduction).
!
      cc(1,:)=                                g5
      cc(2,:)=dx1   *0.5 *(             -g4     +  g6           )
      cc(3,:)=dx2   *0.5 *(              g4-2*g5+  g6           )
      cc(4,:)=dz1   *0.5 *(     -g2                     +  g8   )
      cc(5,:)=dz2   *0.5 *(      g2        -2*g5        +  g8   )
      cc(6,:)=dx1dz1*0.25*( g1     -g3               -g7     +g9)
      cc(7,:)=dx2dz1*0.25*(-g1+2*g2-g3               +g7-2*g8+g9)
      cc(8,:)=dx1dz2*0.25*(-g1     +g3+2*g4     -2*g6-g7     +g9)
      cc(9,:)=dx2dz2*0.25*( g1-2*g2+g3-2*g4+4*g5-2*g6+g7-2*g8+g9)
!
!  Calculate the value of the interpolation function at the point (dxp,dzp).
!
      dxp=xxp(1)-xb(ix0,ib)
      dzp=xxp(3)-zb(iz0,ib)
!
      gp = cc(1,:)            + cc(2,:)*dxp        + cc(3,:)*dxp**2        + &
           cc(4,:)*dzp        + cc(5,:)*dzp**2     + cc(6,:)*dxp*dzp       + &
           cc(7,:)*dxp**2*dzp + cc(8,:)*dxp*dzp**2 + cc(9,:)*dxp**2*dzp**2
!
      call keep_compiler_quiet(ipar)
!
      if (lfirstcall) lfirstcall=.false.
!
    endsubroutine interpolate_quadratic
!***********************************************************************
    subroutine interpolate_quadratic_spline(f,ivar1,ivar2,xxp,gp,inear,iblock,ipar)
!
!  Quadratic spline interpolation of the function g to the point xxp=(xp,yp,zp).
!
!  10-jun-06/anders: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: ivar1, ivar2
      real, dimension (3) :: xxp
      real, dimension (ivar2-ivar1+1) :: gp
      integer, dimension (3) :: inear
      integer :: iblock, ipar
!
      real :: fac_x_m1, fac_x_00, fac_x_p1
      real :: fac_y_m1, fac_y_00, fac_y_p1
      real :: fac_z_m1, fac_z_00, fac_z_p1
      real :: dxp0, dyp0, dzp0
      integer :: ix0, iy0, iz0, ib
!
      intent(in)  :: xxp, ivar1, inear
      intent(out) :: gp
!
      call keep_compiler_quiet(f)
!
!  Abbreviations.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      ib=iblock
!
!  Redefine the interpolation point in coordinates relative to nearest grid
!  point and normalize with the cell size.
!
      dxp0=(xxp(1)-xb(ix0,ib))*dx1b(ix0,ib)
      dyp0=(xxp(2)-yb(iy0,ib))*dy1b(iy0,ib)
      dzp0=(xxp(3)-zb(iz0,ib))*dz1b(iz0,ib)
!
!  Interpolation formulae.
!
      if (dimensionality==0) then
        gp=fb(ix0,iy0,iz0,ivar1:ivar2,ib)
      elseif (dimensionality==1) then
        if (nxgrid/=1) then
          gp = 0.5*(0.5-dxp0)**2*fb(ix0-1,iy0,iz0,ivar1:ivar2,ib) + &
                  (0.75-dxp0**2)*fb(ix0  ,iy0,iz0,ivar1:ivar2,ib) + &
               0.5*(0.5+dxp0)**2*fb(ix0+1,iy0,iz0,ivar1:ivar2,ib)
        endif
        if (nygrid/=1) then
          gp = 0.5*(0.5-dyp0)**2*fb(ix0,iy0-1,iz0,ivar1:ivar2,ib) + &
                  (0.75-dyp0**2)*fb(ix0,iy0  ,iz0,ivar1:ivar2,ib) + &
               0.5*(0.5+dyp0)**2*fb(ix0,iy0+1,iz0,ivar1:ivar2,ib)
        endif
        if (nzgrid/=1) then
          gp = 0.5*(0.5-dzp0)**2*fb(ix0,iy0,iz0-1,ivar1:ivar2,ib) + &
                  (0.75-dzp0**2)*fb(ix0,iy0,iz0  ,ivar1:ivar2,ib) + &
               0.5*(0.5+dzp0)**2*fb(ix0,iy0,iz0+1,ivar1:ivar2,ib)
        endif
      elseif (dimensionality==2) then
        if (nxgrid==1) then
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_y_00*fac_z_00*fb(ix0,iy0,iz0,ivar1:ivar2,ib) + &
              fac_y_00*( fb(ix0,iy0  ,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                         fb(ix0,iy0  ,iz0-1,ivar1:ivar2,ib)*fac_z_m1 ) + &
              fac_z_00*( fb(ix0,iy0+1,iz0  ,ivar1:ivar2,ib)*fac_y_p1 + &
                         fb(ix0,iy0-1,iz0  ,ivar1:ivar2,ib)*fac_y_m1 ) + &
              fac_y_p1*( fb(ix0,iy0+1,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                         fb(ix0,iy0+1,iz0-1,ivar1:ivar2,ib)*fac_z_m1 ) + &
              fac_y_m1*( fb(ix0,iy0-1,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                         fb(ix0,iy0-1,iz0-1,ivar1:ivar2,ib)*fac_z_m1 )
        elseif (nygrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_x_00*fac_z_00*fb(ix0,iy0,iz0,ivar1:ivar2,ib) + &
              fac_x_00*( fb(ix0  ,iy0,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                         fb(ix0  ,iy0,iz0-1,ivar1:ivar2,ib)*fac_z_m1 ) + &
              fac_z_00*( fb(ix0+1,iy0,iz0  ,ivar1:ivar2,ib)*fac_x_p1 + &
                         fb(ix0-1,iy0,iz0  ,ivar1:ivar2,ib)*fac_x_m1 ) + &
              fac_x_p1*( fb(ix0+1,iy0,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                         fb(ix0+1,iy0,iz0-1,ivar1:ivar2,ib)*fac_z_m1 ) + &
              fac_x_m1*( fb(ix0-1,iy0,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                         fb(ix0-1,iy0,iz0-1,ivar1:ivar2,ib)*fac_z_m1 )
        elseif (nzgrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
!
          gp= fac_x_00*fac_y_00*fb(ix0,iy0,iz0,ivar1:ivar2,ib) + &
              fac_x_00*( fb(ix0  ,iy0+1,iz0,ivar1:ivar2,ib)*fac_y_p1 + &
                         fb(ix0  ,iy0-1,iz0,ivar1:ivar2,ib)*fac_y_m1 ) + &
              fac_y_00*( fb(ix0+1,iy0  ,iz0,ivar1:ivar2,ib)*fac_x_p1 + &
                         fb(ix0-1,iy0  ,iz0,ivar1:ivar2,ib)*fac_x_m1 ) + &
              fac_x_p1*( fb(ix0+1,iy0+1,iz0,ivar1:ivar2,ib)*fac_y_p1 + &
                         fb(ix0+1,iy0-1,iz0,ivar1:ivar2,ib)*fac_y_m1 ) + &
              fac_x_m1*( fb(ix0-1,iy0+1,iz0,ivar1:ivar2,ib)*fac_y_p1 + &
                         fb(ix0-1,iy0-1,iz0,ivar1:ivar2,ib)*fac_y_m1 )
        endif
      elseif (dimensionality==3) then
        fac_x_m1 = 0.5*(0.5-dxp0)**2
        fac_x_00 = 0.75-dxp0**2
        fac_x_p1 = 0.5*(0.5+dxp0)**2
        fac_y_m1 = 0.5*(0.5-dyp0)**2
        fac_y_00 = 0.75-dyp0**2
        fac_y_p1 = 0.5*(0.5+dyp0)**2
        fac_z_m1 = 0.5*(0.5-dzp0)**2
        fac_z_00 = 0.75-dzp0**2
        fac_z_p1 = 0.5*(0.5+dzp0)**2
!
        gp= fac_x_00*fac_y_00*fac_z_00*fb(ix0,iy0,iz0,ivar1:ivar2,ib) + &
            fac_x_00*fac_y_00*(fb(ix0  ,iy0  ,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                               fb(ix0  ,iy0  ,iz0-1,ivar1:ivar2,ib)*fac_z_m1)+ &
            fac_x_00*fac_z_00*(fb(ix0  ,iy0+1,iz0  ,ivar1:ivar2,ib)*fac_y_p1 + &
                               fb(ix0  ,iy0-1,iz0  ,ivar1:ivar2,ib)*fac_y_m1)+ &
            fac_y_00*fac_z_00*(fb(ix0+1,iy0  ,iz0  ,ivar1:ivar2,ib)*fac_x_p1 + &
                               fb(ix0-1,iy0  ,iz0  ,ivar1:ivar2,ib)*fac_x_m1)+ &
            fac_x_p1*fac_y_p1*(fb(ix0+1,iy0+1,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                               fb(ix0+1,iy0+1,iz0-1,ivar1:ivar2,ib)*fac_z_m1)+ &
            fac_x_p1*fac_y_m1*(fb(ix0+1,iy0-1,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                               fb(ix0+1,iy0-1,iz0-1,ivar1:ivar2,ib)*fac_z_m1)+ &
            fac_x_m1*fac_y_p1*(fb(ix0-1,iy0+1,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                               fb(ix0-1,iy0+1,iz0-1,ivar1:ivar2,ib)*fac_z_m1)+ &
            fac_x_m1*fac_y_m1*(fb(ix0-1,iy0-1,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                               fb(ix0-1,iy0-1,iz0-1,ivar1:ivar2,ib)*fac_z_m1)+ &
            fac_x_00*fac_y_p1*(fb(ix0  ,iy0+1,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                               fb(ix0  ,iy0+1,iz0-1,ivar1:ivar2,ib)*fac_z_m1)+ &
            fac_x_00*fac_y_m1*(fb(ix0  ,iy0-1,iz0+1,ivar1:ivar2,ib)*fac_z_p1 + &
                               fb(ix0  ,iy0-1,iz0-1,ivar1:ivar2,ib)*fac_z_m1)+ &
            fac_y_00*fac_z_p1*(fb(ix0+1,iy0  ,iz0+1,ivar1:ivar2,ib)*fac_x_p1 + &
                               fb(ix0-1,iy0  ,iz0+1,ivar1:ivar2,ib)*fac_x_m1)+ &
            fac_y_00*fac_z_m1*(fb(ix0+1,iy0  ,iz0-1,ivar1:ivar2,ib)*fac_x_p1 + &
                               fb(ix0-1,iy0  ,iz0-1,ivar1:ivar2,ib)*fac_x_m1)+ &
            fac_z_00*fac_x_p1*(fb(ix0+1,iy0+1,iz0  ,ivar1:ivar2,ib)*fac_y_p1 + &
                               fb(ix0+1,iy0-1,iz0  ,ivar1:ivar2,ib)*fac_y_m1)+ &
            fac_z_00*fac_x_m1*(fb(ix0-1,iy0+1,iz0  ,ivar1:ivar2,ib)*fac_y_p1 + &
                               fb(ix0-1,iy0-1,iz0  ,ivar1:ivar2,ib)*fac_y_m1 )
      endif
!
      call keep_compiler_quiet(ipar)
!
    endsubroutine interpolate_quadratic_spline
!***********************************************************************
    subroutine sort_particles_imn(fp,ineargrid,ipar,dfp,f)
!
!  Sort the particles so that they appear in the same order as the (m,n) loop.
!
!  16-nov-09/anders: dummy
!
      real, dimension (mx,my,mz,mfarray),optional :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
!
      call fatal_error('sort_particles_imn', &
          'not implemented for block domain decomposition')
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(ipar)
      call keep_compiler_quiet(dfp)
!
    endsubroutine sort_particles_imn
!***********************************************************************
    subroutine boundcond_neighbour_list
!
! Copy the number of neighbours to the boundary points of the 
! neighbour list
!
! Dummy so far.
!
!  12-qpr-15/MR: added
!
    endsubroutine boundcond_neighbour_list
!***********************************************************************
    subroutine shepherd_neighbour_pencil(fp,ineargrid,kshepherd,kneighbour)
!
!  Create a shepherd/neighbour list of particles in the pencil.
!
!  16-nov-09/anders: dummy
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (nx) :: kshepherd
      integer, dimension (:) :: kneighbour
!
      intent (in) :: fp, ineargrid
      intent (out) :: kshepherd, kneighbour
!
      call fatal_error('shepherd_neighbour_pencil', &
          'not implemented for block domain decomposition')
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(kshepherd)
      call keep_compiler_quiet(kneighbour)
!
    endsubroutine shepherd_neighbour_pencil
!***********************************************************************
    subroutine shepherd_neighbour_block(fp,ineargrid,kshepherd,kneighbour, &
        iblock)
!
!  Create a shepherd/neighbour list of particles in the block.
!
!  17-nov-09/anders: coded
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (nxb,nyb,nzb) :: kshepherd
      integer, dimension (:) :: kneighbour
      integer :: iblock
!
      integer :: k, ix0, iy0, iz0
!
      intent (in) :: fp, ineargrid
      intent (out) :: kshepherd, kneighbour
!
      call keep_compiler_quiet(fp)
!
      kshepherd=0
      if (iblock==0) kneighbour=0
!
      if (npar_iblock(iblock)/=0) then
        do k=k1_iblock(iblock),k2_iblock(iblock)
          ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
          kneighbour(k)=kshepherd(ix0-nghostb,iy0-nghostb,iz0-nghostb)
          kshepherd(ix0-nghostb,iy0-nghostb,iz0-nghostb)=k
        enddo
      endif
!
    endsubroutine shepherd_neighbour_block
!***********************************************************************
    subroutine shepherd_neighbour_pencil3d(fp,ineargrid,kshepherd,kneighbour)
!
!  17-dec-2011: ought to be coded by AlexHubbard
!
!  Create a shepherd/neighbour list of particles in the pencil.
!  On collisional grid
!  Adapted from particles_map
!
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (:,:,:) :: kshepherd
      integer, dimension (mpar_loc) :: kneighbour
!
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(kshepherd)
      call keep_compiler_quiet(kneighbour)
!
      call not_implemented("shepherd_neighbour_pencil3d", &
           "This is just a dummy coded to make the code compile")
!
    endsubroutine shepherd_neighbour_pencil3d
!***********************************************************************
    subroutine interpolation_consistency_check()
!
!  Check that all interpolation requirements are satisfied:
!
!  16-nov-09/anders: dummy
!
      if (interp%luu .or. interp%loo .or. interp%lTT) &
          call fatal_error('interpolation_consistency_check', &
          'not implemented for block domain decomposition')
!
    endsubroutine interpolation_consistency_check
!***********************************************************************
    subroutine interpolate_quantities(f,fp,p,ineargrid)
!
!  Interpolate the needed sub-grid quantities according to preselected
!  interpolation policies.
!
!  16-nov-09/anders: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
      type (pencil_case) :: p
!
      if (interp%luu .or. interp%loo .or. interp%lTT .or. interp%lgradTT) &
          call fatal_error('interpolate_quantities', &
          'not implemented for block domain decomposition')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine interpolate_quantities
!***********************************************************************
    subroutine cleanup_interpolated_quantities
!
!  Deallocate memory from particle pencil interpolation variables
!
!  16-nov-09/anders: dummy
!
      if (interp%luu .or. interp%loo .or. interp%lTT) &
          call fatal_error('cleanup_interpolated_quantities', &
          'not implemented for block domain decomposition')
!
    endsubroutine cleanup_interpolated_quantities
!***********************************************************************
    subroutine interp_field_pencil_0(f,i1,i2,fp,ineargrid,vec,policy)
!
!  Overloaded interpolation wrapper for scalar fields.
!
!  16-nov-09/anders: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i1,i2
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
      real, dimension(:) :: vec
      integer :: policy
!
      if (interp%luu .or. interp%loo .or. interp%lTT) &
          call fatal_error('interp_field_pencil_0', &
          'not implemented for block domain decomposition')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(i1)
      call keep_compiler_quiet(i2)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(vec)
      call keep_compiler_quiet(policy)
!
    endsubroutine interp_field_pencil_0
!***********************************************************************
    subroutine interp_field_pencil_1(f,i1,i2,fp,ineargrid,vec,policy)
!
!  Overloaded interpolation wrapper for vector fields.
!
!  16-nov-09/anders: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i1,i2
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
      real, dimension(:,:) :: vec
      integer :: policy
!
      if (interp%luu .or. interp%loo .or. interp%lTT) &
          call fatal_error('interp_field_pencil_1', &
          'not implemented for block domain decomposition')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(i1)
      call keep_compiler_quiet(i2)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(vec)
      call keep_compiler_quiet(policy)
!
    endsubroutine interp_field_pencil_1
!***********************************************************************
    subroutine interp_field_pencil(f,i1,i2,fp,ineargrid,uvec2,vec,policy)
!
!  Interpolate stream field to all sub grid particle positions in the
!  current pencil.
!
!  16-nov-09/anders: dummy
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: i1,i2
      real, dimension (mpar_loc,mparray) :: fp
      integer, dimension(mpar_loc,3) :: ineargrid
      integer :: uvec2,policy
      real, dimension(k1_imn(imn):k2_imn(imn),uvec2) :: vec
!
      if (interp%luu .or. interp%loo .or. interp%lTT) &
          call fatal_error('interp_field_pencil', &
          'not implemented for block domain decomposition')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(i1)
      call keep_compiler_quiet(i2)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
      call keep_compiler_quiet(vec)
      call keep_compiler_quiet(policy)
!
    endsubroutine interp_field_pencil
!***********************************************************************
!***********************************************************************
! LOCAL SUBROUTINES OR FUNCTIONS
!***********************************************************************
!***********************************************************************
    elemental subroutine tsc_index_range(ix0, nxgrid, ixx0, ixx1)
!
!  Get the index range surrounding a particle for TSC.
!
!  17-may-23/ccyang: coded
!
      integer, intent(in) :: ix0, nxgrid
      integer, intent(out) :: ixx0, ixx1
!
      range: if (nxgrid > 1) then
        ixx0 = ix0 - 1
        ixx1 = ix0 + 1
      else range
        ixx0 = ix0
        ixx1 = ix0
      endif range
!
    endsubroutine tsc_index_range
!***********************************************************************
endmodule Particles_map
