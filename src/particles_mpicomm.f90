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
!  Redistribute particles among processors based on the local yz-interval
!  of each processor.
!
!  01-jan-05/anders: coded
!
!  TODO: For ncpus>>1000 this subroutine possibly uses too much memory.
!
      use Mpicomm
      use Diagnostics, only: max_name
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mpar_loc,mpvar), optional :: dfp
      logical, optional :: linsert
      integer, dimension (mpar_loc) :: ipar
      integer :: npar_loc
!
      real, dimension (npar_mig,mpvar) :: fp_mig, dfp_mig
      integer, dimension (npar_mig) :: ipar_mig, iproc_rec_array
      integer, dimension (npar_mig) :: isort_array
      integer, dimension (0:ncpus-1) :: nmig_leave, nmig_enter
      integer, dimension (0:ncpus-1) :: ileave_low, ileave_high
      integer, dimension (0:ncpus-1) :: iproc_rec_count
      integer, dimension (26), save :: iproc_comm=-1
      integer, save :: nproc_comm=0
      integer :: dipx, dipy, dipz
      integer :: i, j, k, iproc_rec, ipx_rec, ipy_rec, ipz_rec
      integer :: nmig_leave_total, ileave_high_max
      logical :: lredo, lredo_all
      integer :: itag_nmig=500, itag_ipar=510, itag_fp=520, itag_dfp=530
!
      intent (inout) :: fp, npar_loc, ipar, dfp
!
!  Create list of processors that we allow migration to and from.
!
      if (lshear .or. (it==1 .and. itsub==1)) then
        iproc_comm=-1; nproc_comm=0
        do dipx=-1,1; do dipy=-1,1; do dipz=-1,1
          ipx_rec=ipx+dipx
          ipy_rec=ipy+dipy
          if (lshear) then
            if (ipx_rec<0)        ipy_rec=ipy_rec-ceiling(deltay/Lxyz_loc(2))
            if (ipx_rec>nprocx-1) ipy_rec=ipy_rec+ceiling(deltay/Lxyz_loc(2))
          endif
          ipz_rec=ipz+dipz
          do while (ipx_rec<0);        ipx_rec=ipx_rec+nprocx; enddo
          do while (ipx_rec>nprocx-1); ipx_rec=ipx_rec-nprocx; enddo
          do while (ipy_rec<0);        ipy_rec=ipy_rec+nprocy; enddo
          do while (ipy_rec>nprocy-1); ipy_rec=ipy_rec-nprocy; enddo
          do while (ipz_rec<0);        ipz_rec=ipz_rec+nprocz; enddo
          do while (ipz_rec>nprocz-1); ipz_rec=ipz_rec-nprocz; enddo
          iproc_rec=ipx_rec+ipy_rec*nprocx+ipz_rec*nprocx*nprocy
          if (iproc_rec/=iproc) then
            if (nproc_comm==0) then
              nproc_comm=nproc_comm+1
              iproc_comm(nproc_comm)=iproc_rec
            elseif ( (.not.any(iproc_rec==iproc_comm(1:nproc_comm))) .and. &
                 (iproc_rec/=iproc) ) then
              nproc_comm=nproc_comm+1
              iproc_comm(nproc_comm)=iproc_rec
            endif
          endif
        enddo; enddo; enddo
      endif
!
!  Possible to iterate until all particles have migrated.
!
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
                'redist_particles_procs: Particle ', ipar(k), &
                ' moves out of proc ', iproc, &
                ' and into proc ', iproc_rec
            if (iproc_rec>=ncpus .or. iproc_rec<0) then
              call warning('redist_particles_procs','',iproc)
              print*, 'redist_particles_procs: receiving proc does not exist'
              print*, 'redist_particles_procs: iproc, iproc_rec=', &
                  iproc, iproc_rec
              print*, 'redist_particles_procs: ipar(k), xxp=', &
                  ipar(k), fp(k,ixp:izp)
              print*, 'redist_particles_procs: x0_mig, x1_mig=', &
                  procx_bounds(ipx), procx_bounds(ipx+1)
              print*, 'redist_particles_procs: y0_mig, y1_mig=', &
                  procy_bounds(ipy), procy_bounds(ipy+1)
              print*, 'redist_particles_procs: z0_mig, z1_mig=', &
                  procz_bounds(ipz), procz_bounds(ipz+1)
              call fatal_error_local("","")
            endif
!
!  Copy migrating particle to the end of the fp array.
!
            nmig_leave(iproc_rec)=nmig_leave(iproc_rec)+1
            nmig_leave_total     =nmig_leave_total     +1
            if (sum(nmig_leave)>npar_mig) then
              if (.not. lstart) then
                print '(a,i3,a,i3,a)', &
                    'redist_particles_procs: too many particles migrating '// &
                    'from proc ', iproc, ' to proc ', iproc_rec
                print*, '                       (npar_mig=', npar_mig, 'nmig=',sum(nmig_leave),')'
              endif
              if (lstart.or.lmigration_redo) then
                if (.not. lstart) then
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
                call fatal_error('redist_particles_procs','')
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
        if (ip<=8) print*, 'redist_particles_procs: iproc, nmigrate = ', &
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
              call mpirecv_int(nmig_enter(i),1,i,itag_nmig)
            else
              do j=0,ncpus-1
                if (iproc/=j) call mpisend_int(nmig_leave(j),1,j,itag_nmig)
              enddo
            endif
          enddo
        else
          do i=1,nproc_comm
            call mpisend_int(nmig_leave(iproc_comm(i)),1,iproc_comm(i),itag_nmig+iproc)
          enddo
          do i=1,nproc_comm
            call mpirecv_int(nmig_enter(iproc_comm(i)),1,iproc_comm(i),itag_nmig+iproc_comm(i))
          enddo
        endif
!
!  Check that there is room for the new particles at each processor.
!
        if (npar_loc+sum(nmig_enter)>mpar_loc) then
          print*, 'redist_particles_proc: Too many particles want to be at proc', iproc
          print*, 'redist_particles_proc: npar_loc, mpar_loc, nmig=', &
              npar_loc, mpar_loc, sum(nmig_enter)
          call fatal_error_local('redist_particles_proc','')
        endif
        call fatal_error_local_collect()
!
!  Sort array of migrating particle in order of receiving processor.
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
!  Set to receive.
!
        do i=0,ncpus-1
          if (iproc/=i .and. nmig_enter(i)/=0) then
            call mpirecv_real(fp(npar_loc+1:npar_loc+nmig_enter(i),:), &
                (/nmig_enter(i),mpvar/),i,itag_fp)
            call mpirecv_int(ipar(npar_loc+1:npar_loc+nmig_enter(i)), &
                nmig_enter(i),i,itag_ipar)
            if (present(dfp)) &
                call mpirecv_real(dfp(npar_loc+1:npar_loc+nmig_enter(i),:), &
                (/nmig_enter(i),mpvar/),i,itag_dfp)
            if (ip<=6) then
              print*, 'redist_particles_procs: iproc, iproc_send=', iproc, i
              print*, 'redist_particles_procs: received fp=', &
                  fp(npar_loc+1:npar_loc+nmig_enter(i),:)
              print*, 'redist_particles_procs: received ipar=', &
                  ipar(npar_loc+1:npar_loc+nmig_enter(i))
              if (present(dfp)) &
                  print*, 'redist_particles_procs: received dfp=',&
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
                    print*, 'redist_particles_procs: received particle '// &
                        'closer to ghost point than to physical grid point!'
                    print*, 'redist_particles_procs: ipar, xxp=', &
                        ipar(k), fp(k,ixp:izp)
                    print*, 'redist_particles_procs: x0_mig, x1_mig=', &
                        procx_bounds(ipx), procx_bounds(ipx+1)
                  endif
                endif
                if (nygrid/=1) then
                  if (fp(k,iyp)<procy_bounds(ipy) .or. &
                      fp(k,iyp)>=procy_bounds(ipy+1)) then
                    print*, 'redist_particles_procs: received particle '// &
                        'closer to ghost point than to physical grid point!'
                    print*, 'redist_particles_procs: ipar, xxp=', &
                        ipar(k), fp(k,ixp:izp)
                    print*, 'redist_particles_procs: y0_mig, y1_mig=', &
                        procy_bounds(ipy), procy_bounds(ipy+1)
                  endif
                endif
                if (nzgrid/=1) then
                  if (fp(k,izp)<procz_bounds(ipz) .or. &
                      fp(k,izp)>=procz_bounds(ipz+1)) then
                    print*, 'redist_particles_procs: received particle '// &
                        'closer to ghost point than to physical grid point!'
                    print*, 'redist_particles_procs: ipar, xxp=', &
                        ipar(k), fp(k,ixp:izp)
                    print*, 'redist_particles_procs: z0_mig, z1_mig=', &
                        procz_bounds(ipz), procz_bounds(ipz+1)
                  endif
                endif
              enddo
            endif
            npar_loc=npar_loc+nmig_enter(i)
            if (npar_loc>mpar_loc) then
              print*, 'redist_particles_proc: Too many particles at proc', iproc
              print*, 'redist_particles_proc: npar_loc, mpar_loc=', &
                  npar_loc, mpar_loc
              call fatal_error('redist_particles_proc','')
            endif
          endif
!
!  Directed send.
!
          if (iproc==i) then
            do j=0,ncpus-1
              if (iproc/=j .and. nmig_leave(j)/=0) then
                call mpisend_real(fp_mig(ileave_low(j):ileave_high(j),:), &
                    (/nmig_leave(j),mpvar/),j,itag_fp)
                call mpisend_int(ipar_mig(ileave_low(j):ileave_high(j)), &
                    nmig_leave(j),j,itag_ipar)
                if (present(dfp)) &
                    call mpisend_real(dfp_mig(ileave_low(j):ileave_high(j),:), &
                    (/nmig_leave(j),mpvar/),j,itag_dfp)
                if (ip<=6) then
                  print*, 'redist_particles_proc: iproc, iproc_rec=', iproc, j
                  print*, 'redist_particles_proc: sent fp=', &
                      fp_mig(ileave_low(j):ileave_high(j),:)
                  print*, 'redist_particles_proc: sent ipar=', &
                      ipar_mig(ileave_low(j):ileave_high(j))
                  if (present(dfp)) &
                      print*, 'redist_particles_proc: sent dfp=', &
                      dfp_mig(ileave_low(j):ileave_high(j),:)
                endif
              endif
            enddo
          endif
        enddo
!
!  Sum up processors that have not had place to let all migrating particles go.
!
        if (lstart.or.lmigration_redo) then   !  5-10% slowdown of code
          call mpireduce_or(lredo, lredo_all)
          call mpibcast_logical(lredo_all, 1)
        else
          lredo_all=.false.
        endif
!
!  If sum is not zero, then the while loop while be executed once more.
!
      enddo
!
    endsubroutine redist_particles_procs
!***********************************************************************
endmodule Particles_mpicomm
