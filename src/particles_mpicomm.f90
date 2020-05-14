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
  include 'mpif.h'
!
  integer, parameter :: nxb=1, nyb=1, nzb=1, nbx=1, nby=1, nbz=1
  integer, parameter :: nbricks=0, nghostb=0, mxb=1, myb=1, mzb=1
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
  real :: xref_par=0.0, yref_par=0.0, zref_par=0.0
  integer :: it1_loadbalance=-1
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
      call keep_compiler_quiet(f)
!
!  Distribute particles evenly among processors to begin with.
!  DM: for some initial conditions it may be better not to distribute
!  particles evenly so the logical variable ldist_particles_evenly has been
!  introduced. This variable is true by default.
!
!  WL: This logical has been removed, right? This comment does not seem to 
!  be relevant anymore. 
!
      if (lstart) call dist_particles_evenly_procs(ipar)
!
!  Set neighbor process ranks for particles.
!
      neighbors_par = neighbors
!
      xdir: if (bcpx /= 'p' .or. nxgrid == 1 .or. lshear .and. nygrid > 1) then
        if (lfirst_proc_x) neighbors_par(-1,:,:) = -1
        if (llast_proc_x) neighbors_par(+1,:,:) = -1
      endif xdir
!
      ydir: if (bcpy /= 'p' .or. nygrid == 1) then
        if (lfirst_proc_y) neighbors_par(:,-1,:) = -1
        if (llast_proc_y) neighbors_par(:,+1,:) = -1
      endif ydir
!
      zdir: if (bcpz /= 'p' .or. nzgrid == 1) then
        if (lfirst_proc_z) neighbors_par(:,:,-1) = -1
        if (llast_proc_z) neighbors_par(:,:,+1) = -1
      endif zdir
!
    endsubroutine initialize_particles_mpicomm
!***********************************************************************
    subroutine dist_particles_evenly_procs(ipar)
!
!  Distribute particles evenly among processors.
!
!  05-jan-05/anders: coded
!
      use Mpicomm,  only: mpibcast_int
!
      integer, dimension (mpar_loc) :: ipar
!
      integer :: i, k, jspec, npar_per_species, npar_rest, icycle, iproc_rec
      integer :: npar_per_species_missed
      real :: ndim,nfracx,nfracy,nfracz
      integer, dimension (0:ncpus-1) :: ipar1, ipar2
      integer, dimension (0:ncpus-1) :: npar_loc_array, npar_rest_array
!
      intent (inout) :: ipar
!
!  Set index interval of particles that belong to the local processor.
!
      if (linsert_particles_continuously) then
        npar_loc=0
        ipar=0
      else
!
!  Place particles evenly on all processors. Some processors may get an extra
!  particle if the particle number is not divisible by the number of processors.
!
        if (lcartesian_coords) then 
          npar_loc =npar/ncpus
          npar_rest=npar-npar_loc*ncpus
          do i=0,ncpus-1
            if (i<npar_rest) then
              npar_loc_array(i)=npar_loc+1
            else
              npar_loc_array(i)=npar_loc
            endif
          enddo
        else
!
!  Take into account that in cylindrical and spherical coordinates the volume occupied by
!  different processors in the radial direction is different. The line is correct also in 
!  the case of nprocx=1, in which case the stuff in parentheses reduces to 1/ncpus.
!
          if (lcylindrical_coords) then 
            ndim=2-dustdensity_powerlaw
          elseif (lspherical_coords) then 
            ndim=3-dustdensity_powerlaw
          else
            call fatal_error("dist_particles_evenly","The world is flat, and we never got here.")
          endif
!
!  Sanity check
!
          if (ndim==0) then
            if (lroot) then
              print*,''
              print*,'For cylindrical coordinates dustdensity_powerlaw cannot be exactly 2.'
              print*,'For spherical coordinates dustdensity_powerlaw cannot be exactly 3.'
              print*,'Because of the dimensionality of the Jacobian, it becomes a mathematical'
              print*,'pole for calculating npar_loc and the dust density distribution.'
            endif
            call fatal_error("dist_particles_evenly","")
          endif
!
          do i=0,ncpus-1 
            nfracx = (xyz1_loc(1)**ndim - xyz0_loc(1)**ndim)/(xyz1(1)**ndim - xyz0(1)**ndim)
!
            nfracy = (xyz1_loc(2) - xyz0_loc(2))/(xyz1(2) - xyz0(2))
            nfracz = (xyz1_loc(3) - xyz0_loc(3))/(xyz1(3) - xyz0(3))
            npar_loc_array(i)=npar * nfracx * nfracy * nfracz
          enddo
!
        endif
        npar_loc=npar_loc_array(iproc)
!
        if (lroot) print*, 'dist_particles_evenly_procs: npar_loc_array     =',&
            npar_loc_array
!
!  If there are zero particles on a processor, set ipar1 and ipar2 to zero.
!
        if (npar_species==1) then
          do i=0,ncpus-1
            if (npar_loc_array(i)==0) then
              ipar1(i)=0
              ipar2(i)=0
            else
              if (i==0) then
                ipar1(i)=1
              else
                ipar1(i)=maxval(ipar2(0:i-1))+1
              endif
              ipar2(i)=ipar1(i) + npar_loc_array(i) - 1
            endif
          enddo
!
          if (lroot) then
            print*, 'dist_particles_evenly_procs: ipar1 =', ipar1
            print*, 'dist_particles_evenly_procs: ipar2 =', ipar2
          endif
!
!  Fill in particle index between ipar1 and ipar2.
!
          do k=1,npar_loc
            ipar(k)=k-1+ipar1(iproc)
          enddo
!
        else !start multiple particle species
!
!  Must have same number of particles in each species.
!
          if (mod(npar,npar_species)/=0) then
            if (lroot) then
              print*, 'dist_particles_evenly_procs: npar_species '// &
                  'must be a whole multiple of npar!'
              print*, 'npar_species, npar=', npar_species, npar
            endif
            call fatal_error('dist_particles_evenly_procs','')
          endif
!
!  Distribute particle species evenly among processors.
!
          npar_per_species=npar/npar_species
!
          do jspec=1,npar_species
            do k=1,npar_per_species/ncpus
              ipar(k+npar_per_species/ncpus*(jspec-1))= &
                 k+npar_per_species*(jspec-1)+iproc*(npar_per_species/ncpus)
            enddo
          enddo
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
          if (lroot) then
            print*, 'dist_particles_evenly_procs: npar_per_species   =', &
                npar_per_species
            print*, 'dist_particles_evenly_procs: ipar_fence_species =', &
                ipar_fence_species
          endif
!
!  It is not always possible to have the same number of particles of each
!  species at each processor. In that case we place the remaining particles
!  by hand in order of increasing processor number.
!
          npar_rest_array=npar_loc_array-npar_per_species/ncpus*npar_species
!
          if (lroot .and. (minval(npar_rest_array)/=0)) &
              print*, 'dist_particles_evenly_procs: npar_rest_array    =', &
              npar_rest_array
!
          npar_per_species_missed = &
              npar_per_species-ncpus*(npar_per_species/ncpus)
          icycle   =1
          iproc_rec=0
          do k=1,npar_per_species_missed*npar_species
            jspec=(k-1)/npar_per_species_missed+1
            if (iproc==iproc_rec) &
                ipar(npar_loc-npar_rest_array(iproc_rec)+icycle)= &
                ipar_fence_species(jspec)-jspec*npar_per_species_missed+k
            if (lroot) print*, 'dist_particles_evenly_procs: placed ', &
                'particle', ipar_fence_species(jspec)- &
                jspec*npar_per_species_missed+k, &
                ' on proc', iproc_rec, ' in species', jspec
            iproc_rec=iproc_rec+1
            if (iproc_rec==ncpus) then
              iproc_rec=0
              icycle=icycle+1
            endif
          enddo
        endif
!
      endif
!      write(*,*)'DM','proc no=',iproc,'npar_loc=',npar_loc
!
    endsubroutine dist_particles_evenly_procs
!***********************************************************************
    subroutine migrate_particles(fp,ipar,dfp,linsert)
!
!  Redistribute particles among processors based on the local yz-interval
!  of each processor.
!
!  01-jan-05/anders: coded
!
!  TODO:
!    - For ncpus>>1000 this subroutine possibly uses too much memory.
!    - Optimize lmigration_redo to not reanalyse particles that have already
!      been considered.
!
      use Mpicomm
      use Diagnostics, only: max_name
!
      real, dimension (mpar_loc,mpcom) :: fp
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
      logical, optional :: linsert
!
      real, dimension (npar_mig,mpcom) :: fp_mig
      real, dimension (npar_mig,mpvar) :: dfp_mig
      integer, dimension (npar_mig) :: ipar_mig, iproc_rec_array
      integer, dimension (npar_mig) :: isort_array
      integer, dimension (0:ncpus-1) :: nmig_leave, nmig_enter
      integer, dimension (0:ncpus-1) :: ileave_low, ileave_high
      integer, dimension (0:ncpus-1) :: iproc_rec_count
      integer, save :: nmig_max = 0
      integer :: i, j, k, iproc_rec, ipx_rec, ipy_rec, ipz_rec
      integer :: nmig_leave_total, nmig_left, ileave_high_max, buffer_max
      real, dimension (:,:), allocatable :: buffer
      logical :: lredo, lredo_all
      integer :: itag_nmig=500, itag_ipar=510, itag_fp=520, itag_dfp=530
!
      intent (inout) :: fp, ipar, dfp
!
!  Possible to iterate until all particles have migrated.
!
      nmig_left = 0
      lredo=.false.; lredo_all=.true.
      do while (lredo_all)
        lredo=.false.
!
!  Find out which particles are not in the local processor's xyz-interval.
!  Need to use special definition of processor boundaries for migration since
!  the physical processor boundaries may be closer to a ghost point than to a
!  physical grid point. Thus we need to define the boundary as the average
!  coordinate value of the last grid point and the first ghost point.
!
        nmig_leave=0
        nmig_leave_total=0
        nmig_enter=0
!
        do k=npar_loc,1,-1
!  Find x index of receiving processor.
          ipx_rec=ipx
          if (fp(k,ixp)>=procx_bounds(ipx+1).and.ipx<nprocx-1) then
            do j=ipx+1,nprocx-1
              if (fp(k,ixp)<procx_bounds(j+1)) then
                ipx_rec=j
                exit
              endif
            enddo
          else if (fp(k,ixp)<procx_bounds(ipx).and.ipx>0) then
            do j=ipx-1,0,-1
              if (fp(k,ixp)>procx_bounds(j)) then
                ipx_rec=j
                exit
              endif
            enddo
          endif
!  Find y index of receiving processor.
          ipy_rec=ipy
          if (fp(k,iyp)>=procy_bounds(ipy+1).and.ipy<nprocy-1) then
            do j=ipy+1,nprocy-1
              if (fp(k,iyp)<procy_bounds(j+1)) then
                ipy_rec=j
                exit
              endif
            enddo
          else if (fp(k,iyp)<procy_bounds(ipy).and.ipy>0) then
            do j=ipy-1,0,-1
              if (fp(k,iyp)>procy_bounds(j)) then
                ipy_rec=j
                exit
              endif
            enddo
          endif
!  Find z index of receiving processor.
          ipz_rec=ipz
          if (fp(k,izp)>=procz_bounds(ipz+1).and.ipz<nprocz-1) then
            do j=ipz+1,nprocz-1
              if (fp(k,izp)<procz_bounds(j+1)) then
                ipz_rec=j
                exit
              endif
            enddo
          else if (fp(k,izp)<procz_bounds(ipz).and.ipz>0) then
            do j=ipz-1,0,-1
              if (fp(k,izp)>procz_bounds(j)) then
                ipz_rec=j
                exit
              endif
            enddo
          endif
!
!  Fix for error that might occur when a particle lands exactly at the
!  box boundary and is assigned to a non-existing processor.
!
          if (lcheck_exact_frontier) then
            if (nprocx/=1) then
!
!  Check if the particle is really closer to a grid cell than to a
!  ghost one. otherwise, this is a more serious particle position
!  problem, that should be allowed to lead to a crash.
!
              if (ipx_rec==-1) then
                if (xyz0(1)-fp(k,ixp)<=(x(l1)-x(l1-1))/2) ipx_rec=0
              endif
              if (ipx_rec==nprocx) then
                if (fp(k,ixp)-xyz1(1)<=(x(l2+1)-x(l2))/2) ipx_rec=nprocx-1
              endif
            endif
            if (nprocy/=1) then
              if (ipy_rec==-1) then
                if (xyz0(2)-fp(k,iyp)<=(y(m1)-y(m1-1))/2) ipy_rec=0
              endif
              if (ipy_rec==nprocy) then
                if (fp(k,iyp)-xyz1(2)<=(y(m2+1)-y(m2))/2) ipy_rec=nprocy-1
              endif
            endif
            if (nprocz/=1) then
              if (ipz_rec==-1) then
                if (xyz0(3)-fp(k,izp)<=(z(n1)-z(n1-1))/2) ipz_rec=0
              endif
              if (ipz_rec==nprocz) then
                if (fp(k,izp)-xyz1(3)<=(z(n2+1)-z(n2))/2) ipz_rec=nprocz-1
              endif
            endif
          endif
!
!  Calculate serial index of receiving processor.
!
          iproc_rec=ipx_rec+nprocx*ipy_rec+nprocx*nprocy*ipz_rec
!
!  Migrate particle if it is no longer at the current processor.
!
          if (iproc_rec/=iproc) then
            if (ip<=7) print '(a,i7,a,i3,a,i3)', &
                'migrate_particles: Particle ', ipar(k), &
                ' moves out of proc ', iproc, &
                ' and into proc ', iproc_rec
!
!  Check that particle wants to migrate to neighbouring processor.
!
            if (.not. (lstart .or. present(linsert)) .and. &
                (.not.any(iproc_rec==iproc_comm(1:nproc_comm))) ) then
              print*, 'migrate_particles: particle ', ipar(k), ' wants to'
              print*, '    migrate to a processor that is not a neighbour!'
              print*, 'migrate_particles: iproc, iproc_rec=', &
                  iproc, iproc_rec
              print*, 'migrate_particles: iproc_comm=', &
                  iproc_comm(1:nproc_comm)
              print*, 'migrate_particles: xxp=', fp(k,ixp:izp)
              print*, 'migrate_particles: deltay=', deltay
              call fatal_error_local("","")
            endif
!
!  Check that particle wants to migrate to existing processor.
!
            if (iproc_rec>=ncpus .or. iproc_rec<0) then
              call warning('migrate_particles','',iproc)
              print*, 'migrate_particles: receiving proc does not exist'
              print*, 'migrate_particles: iproc, iproc_rec=', &
                  iproc, iproc_rec
              print*, 'migrate_particles: ipar(k), xxp=', &
                  ipar(k), fp(k,ixp:izp)
              print*, 'migrate_particles: x0_mig, x1_mig=', &
                  procx_bounds(ipx), procx_bounds(ipx+1)
              print*, 'migrate_particles: y0_mig, y1_mig=', &
                  procy_bounds(ipy), procy_bounds(ipy+1)
              print*, 'migrate_particles: z0_mig, z1_mig=', &
                  procz_bounds(ipz), procz_bounds(ipz+1)
              call fatal_error_local("","")
            endif
!
!  Copy migrating particle to the end of the fp array.
!
            nmig_leave(iproc_rec)=nmig_leave(iproc_rec)+1
            nmig_leave_total     =nmig_leave_total     +1
            if (sum(nmig_leave)>npar_mig) then
              if (lrun) then
                print '(a,i3,a,i3,a)', &
                    'migrate_particles: too many particles migrating '// &
                    'from proc ', iproc, ' to proc ', iproc_rec
                print*, '                       (npar_mig=', npar_mig, 'nmig=',sum(nmig_leave),')'
              endif
              if (lstart.or.lmigration_redo) then
                if (lrun) then
                  print*, '                       Going to do one more '// &
                      'migration iteration!'
                  print*, '                       (this is time consuming - '//&
                      'consider setting npar_mig'
                  print*, '                        higher in cparam.local)'
                endif
                nmig_leave(iproc_rec)=nmig_leave(iproc_rec)-1
                nmig_leave_total     =nmig_leave_total     -1
                lredo=.true.
                exit
              else
                call fatal_error('migrate_particles','too many particles migrating',.true.)
              endif
            endif
            fp_mig(nmig_leave_total,:)=fp(k,:)
            if (present(dfp)) dfp_mig(nmig_leave_total,:)=dfp(k,:)
            ipar_mig(nmig_leave_total)=ipar(k)
            iproc_rec_array(nmig_leave_total)=iproc_rec
!
!  Move the particle with the highest index number to the empty spot left by
!  the migrating particle.
!
            fp(k,:)=fp(npar_loc,:)
            if (present(dfp)) dfp(k,:)=dfp(npar_loc,:)
            ipar(k)=ipar(npar_loc)
!
!  Reduce the number of particles by one.
!
            npar_loc=npar_loc-1
          endif
        enddo
!
!  Print out information about number of migrating particles.
!
        if (ip<=8) print*, 'migrate_particles: iproc, nmigrate = ', &
            iproc, sum(nmig_leave)
!
!  Diagnostic about number of migrating particles.
!  WARNING: in time-steps where snapshots are written, this diagnostic
!  parameter will be zero (quite confusing)!
!
        if (ldiagnos.and.(idiag_nmigmax/=0)) &
            call max_name(sum(nmig_leave),idiag_nmigmax)
!
!  Share information about number of migrating particles. For the initial
!  condition we allow all processors to communicate particles. However, this
!  is extremely slow when the processor number is high (>>100). Thus we only
!  allow communication with surrounding processors during the run.
!
        if (lstart .or. present(linsert)) then
          do i=0,ncpus-1
            if (iproc/=i) then
              call mpirecv_int(nmig_enter(i),i,itag_nmig)
            else
              do j=0,ncpus-1
                if (iproc/=j) call mpisend_int(nmig_leave(j),j,itag_nmig)
              enddo
            endif
          enddo
        else
          do i=1,nproc_comm
            call mpisend_int(nmig_leave(iproc_comm(i)),iproc_comm(i),itag_nmig+iproc)
          enddo
          do i=1,nproc_comm
            call mpirecv_int(nmig_enter(iproc_comm(i)),iproc_comm(i),itag_nmig+iproc_comm(i))
          enddo
        endif
!
!  Check that there is room for the new particles at each processor.
!
        if (npar_loc+sum(nmig_enter)>mpar_loc) then
          print*, 'migrate_particles: Too many particles want to be at proc', iproc
          print*, 'migrate_particles: npar_loc, mpar_loc, nmig=', &
              npar_loc, mpar_loc, sum(nmig_enter)
          call fatal_error_local('migrate_particles','')
        endif
        call fatal_error_local_collect()
!
!  Sort array of migrating particles in order of receiving processor.
!
        if (nmig_leave_total>=1) then
          ileave_high_max=0
          do i=0,ncpus-1
            if (nmig_leave(i)>=1) then
              ileave_low(i)  =ileave_high_max+1
              ileave_high(i) =ileave_low(i)+nmig_leave(i)-1
              ileave_high_max=ileave_high_max+nmig_leave(i)
            else
              ileave_low(i) =0
              ileave_high(i)=0
            endif
          enddo
          iproc_rec_count=0
          do k=1,nmig_leave_total
            isort_array(ileave_low(iproc_rec_array(k))+iproc_rec_count(iproc_rec_array(k)))=k
            iproc_rec_count(iproc_rec_array(k))= &
                iproc_rec_count(iproc_rec_array(k))+1
          enddo
          ipar_mig(1:nmig_leave_total)= &
              ipar_mig(isort_array(1:nmig_leave_total))
          fp_mig(1:nmig_leave_total,:)= &
              fp_mig(isort_array(1:nmig_leave_total),:)
          if (present(dfp)) dfp_mig(1:nmig_leave_total,:)= &
              dfp_mig(isort_array(1:nmig_leave_total),:)
        endif
!
!  Set communication buffers.
!
        buffer_max = 0
        do i=0,ncpus-1
          buffer_max = max(buffer_max, nmig_leave(i))
          buffer_max = max(buffer_max, nmig_enter(i))
        enddo
        call mpireduce_max_int(max (0,buffer_max), buffer_max)
        call mpibcast(buffer_max)
        if (buffer_max > 0) allocate (buffer(buffer_max,mpcom))
!
!  Set to receive.
!
        do i=0,ncpus-1
          if (iproc/=i .and. nmig_enter(i)/=0) then
            call mpirecv_real(buffer,(/buffer_max,mpcom/),i,itag_fp+2)
            fp(npar_loc+1:npar_loc+nmig_enter(i),:) = buffer(1:nmig_enter(i),:)
            call mpirecv_int(ipar(npar_loc+1:npar_loc+nmig_enter(i)),nmig_enter(i),i,itag_ipar)
            if (present(dfp)) then
!              call mpirecv_real(dfp(npar_loc+1:npar_loc+nmig_enter(i),:),(/nmig_enter(i),mpvar/),i,itag_dfp)
              call mpirecv_real(buffer(:,1:mpvar),(/buffer_max,mpvar/),i,itag_dfp)
              dfp(npar_loc+1:npar_loc+nmig_enter(i),:) = buffer(1:nmig_enter(i),1:mpvar)
            endif
            if (ip<=6) then
              print*, 'migrate_particles: iproc, iproc_send=', iproc, i
              print*, 'migrate_particles: received fp=', &
                  fp(npar_loc+1:npar_loc+nmig_enter(i),:)
              print*, 'migrate_particles: received ipar=', &
                  ipar(npar_loc+1:npar_loc+nmig_enter(i))
              if (present(dfp)) &
                  print*, 'migrate_particles: received dfp=',&
                  dfp(npar_loc+1:npar_loc+nmig_enter(i),:)
            endif
!
!  Check that received particles are really at the right processor.
!
            if (lmigration_real_check) then
              do k=npar_loc+1,npar_loc+nmig_enter(i)
                if (nxgrid/=1) then
                  if (fp(k,ixp)<procx_bounds(ipx) .or. &
                      fp(k,ixp)>=procx_bounds(ipx+1)) then
                    print*, 'migrate_particles: received particle '// &
                        'closer to ghost point than to physical grid point!'
                    print*, 'migrate_particles: ipar, xxp=', &
                        ipar(k), fp(k,ixp:izp)
                    print*, 'migrate_particles: x0_mig, x1_mig=', &
                        procx_bounds(ipx), procx_bounds(ipx+1)
                  endif
                endif
                if (nygrid/=1) then
                  if (fp(k,iyp)<procy_bounds(ipy) .or. &
                      fp(k,iyp)>=procy_bounds(ipy+1)) then
                    print*, 'migrate_particles: received particle '// &
                        'closer to ghost point than to physical grid point!'
                    print*, 'migrate_particles: ipar, xxp=', &
                        ipar(k), fp(k,ixp:izp)
                    print*, 'migrate_particles: y0_mig, y1_mig=', &
                        procy_bounds(ipy), procy_bounds(ipy+1)
                  endif
                endif
                if (nzgrid/=1) then
                  if (fp(k,izp)<procz_bounds(ipz) .or. &
                      fp(k,izp)>=procz_bounds(ipz+1)) then
                    print*, 'migrate_particles: received particle '// &
                        'closer to ghost point than to physical grid point!'
                    print*, 'migrate_particles: ipar, xxp=', &
                        ipar(k), fp(k,ixp:izp)
                    print*, 'migrate_particles: z0_mig, z1_mig=', &
                        procz_bounds(ipz), procz_bounds(ipz+1)
                  endif
                endif
              enddo
            endif
            npar_loc=npar_loc+nmig_enter(i)
            if (npar_loc>mpar_loc) then
              print*, 'migrate_particles: Too many particles at proc', iproc
              print*, 'migrate_particles: npar_loc, mpar_loc=', &
                  npar_loc, mpar_loc
              call fatal_error('migrate_particles','')
            endif
          endif
!
!  Directed send.
!
          if (iproc==i) then
            do j=0,ncpus-1
              if (iproc/=j .and. nmig_leave(j)/=0) then
                buffer(1:nmig_leave(j),:) = fp_mig(ileave_low(j):ileave_high(j),:)
                call mpisend_real(buffer,(/buffer_max,mpcom/),j,itag_fp+2)
                call mpisend_int(ipar_mig(ileave_low(j):ileave_high(j)),nmig_leave(j),j,itag_ipar)
                if (present(dfp)) then
!                  call mpisend_real(dfp_mig(ileave_low(j):ileave_high(j),:),(/nmig_leave(j),mpvar/),j,itag_dfp)
                  buffer(1:nmig_leave(j),1:mpvar) = dfp_mig(ileave_low(j):ileave_high(j),:)
                  call mpisend_real(buffer(:,1:mpvar),(/buffer_max,mpvar/),j,itag_dfp)
                endif
                if (ip<=6) then
                  print*, 'migrate_particles: iproc, iproc_rec=', iproc, j
                  print*, 'migrate_particles: sent fp=', &
                      fp_mig(ileave_low(j):ileave_high(j),:)
                  print*, 'migrate_particles: sent ipar=', &
                      ipar_mig(ileave_low(j):ileave_high(j))
                  if (present(dfp)) &
                      print*, 'migrate_particles: sent dfp=', &
                      dfp_mig(ileave_low(j):ileave_high(j),:)
                endif
              endif
            enddo
          endif
        enddo
        if (buffer_max > 0) deallocate (buffer)
        nmig_left = nmig_left + nmig_leave_total
!
!  Sum up processors that have not had place to let all migrating particles go.
!
        if (lstart.or.lmigration_redo) then   !  5-10% slowdown of code
          call mpireduce_or(lredo, lredo_all)
          call mpibcast_logical(lredo_all)
        else
          lredo_all=.false.
        endif
!
!  If sum is not zero, then the while loop will be executed once more.
!
      enddo
!
!  Report maximum number of particles migrated per time step since the last report.
!
      nmig_max = max(nmig_max, nmig_left)
      nmigmmax: if (ldiagnos .and. idiag_nmigmmax /= 0) then
        call max_name(nmig_max, idiag_nmigmmax)
        nmig_max = 0
      endif nmigmmax
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
      real, dimension (mpar_loc,mpcom) :: fp
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
    subroutine sort_blocks()
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
    subroutine communicate_fpbuf(to_neigh,from_neigh,her_npbuf,my_npbuf)

      use Mpicomm
      integer, intent(in) :: to_neigh,from_neigh
      integer, intent(in) :: her_npbuf,my_npbuf
      integer :: comm,mpierr
      integer :: tolow=-1,toup=1
      integer :: irecv_rq_fromlow,irecv_rq_fromupp,isend_rq_tolow,isend_rq_toupp
      integer, dimension (MPI_STATUS_SIZE) :: irecv_stat_fl,irecv_stat_fu,&
           isend_stat_tl,isend_stat_tu
!---------
      comm=MPI_COMM_WORLD
!
! We assume below that send and receive is non-blocking.
!
      call MPI_IRECV(her_npbuf,1,MPI_INTEGER, &
           from_neigh,toup,comm,irecv_rq_fromlow,mpierr)
      call MPI_ISEND(my_npbuf,1,MPI_INTEGER, &
           to_neigh,toup,comm,isend_rq_toupp,mpierr)
      call MPI_WAIT(irecv_rq_fromlow,irecv_stat_fl,mpierr)
      call MPI_WAIT(isend_rq_toupp,isend_stat_tu,mpierr)
!
      call MPI_IRECV(fp_buffer_in,mparray*nslab,MPI_REAL, &
           from_neigh,toup,comm,irecv_rq_fromlow,mpierr)
      call MPI_ISEND(fp_buffer_out,mparray*nslab,MPI_REAL, &
           to_neigh,toup,comm,isend_rq_toupp,mpierr)
      call MPI_WAIT(irecv_rq_fromlow,irecv_stat_fl,mpierr)
      call MPI_WAIT(isend_rq_toupp,isend_stat_tu,mpierr)
!
    endsubroutine communicate_fpbuf
!***********************************************************************
endmodule Particles_mpicomm
