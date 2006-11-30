! $Id: particles_sub.f90,v 1.93 2006-11-30 09:03:36 dobler Exp $
!
!  This module contains subroutines useful for the Particle module.
!
module Particles_sub

  use Cdata
  use Particles_cdata

  implicit none

  private

  public :: input_particles, output_particles
  public :: wsnap_particles, boundconds_particles
  public :: redist_particles_procs, dist_particles_evenly_procs
  public :: sum_par_name, max_par_name, sum_par_name_nw, integrate_par_name
  public :: interpolate_linear
  public :: interpolate_quadratic, interpolate_quadratic_spline
  public :: map_nearest_grid, map_xxp_grid, sort_particles_imn
  public :: particle_pencil_index, find_closest_gridpoint
  public :: shepherd_neighbour

  contains

!***********************************************************************
    subroutine input_particles(filename,fp,npar_loc,ipar)
!
!  Read snapshot file with particle data.
!
!  29-dec-04/anders: adapted from input
!
      real, dimension (mpar_loc,mpvar) :: fp
      character (len=*) :: filename
      integer, dimension (mpar_loc) :: ipar
      integer :: npar_loc
!
      intent (in) :: filename
      intent (out) :: fp,npar_loc,ipar
!
      open(1,FILE=filename,FORM='unformatted')
!
!  First read the number of particles present at the processor and the index
!  numbers of the particles.
!
        read(1) npar_loc
        if (npar_loc/=0) read(1) ipar(1:npar_loc)
!
!  Then read particle data.
!
        if (npar_loc/=0) read(1) fp(1:npar_loc,:)
!
!  Read snapshot time.
!
!        read(1) t
!
        if (ip<=8) print*, 'input_particles: read ', filename
!
      close(1)
!
    endsubroutine input_particles
!***********************************************************************
    subroutine output_particles(filename,fp,npar_loc,ipar)
!
!  Write snapshot file with particle data.
!
!  29-dec-04/anders: adapted from output
!
      use IO, only: lun_output
!
      character(len=*) :: filename
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension(mpar_loc) :: ipar
      integer :: npar_loc
!
      intent (in) :: filename, npar_loc, ipar
!
      if (ip<=8.and.lroot) print*,'output_particles: writing snapshot file '// &
          filename
!
      open(lun_output,FILE=filename,FORM='unformatted')
!
!  First write the number of particles present at the processor and the index
!  numbers of the particles.
!
        write(lun_output) npar_loc
        if (npar_loc/=0) write(lun_output) ipar(1:npar_loc)
!
!  Then write particle data.
!
        if (npar_loc/=0) write(lun_output) fp(1:npar_loc,:)
!
!  Write time and grid parameters.
!
        write(lun_output) t, x, y, z, dx, dy, dz
!
      close(lun_output)
!
    endsubroutine output_particles
!***********************************************************************
    subroutine wsnap_particles(snapbase,fp,enum,lsnap,dsnap_par_minor, &
        npar_loc,ipar,flist)
!
!  Write particle snapshot file, labelled consecutively if enum==.true.
!  Otherwise just write a snapshot without label (used e.g. for pvar.dat)
!
!  29-dec-04/anders: adapted from wsnap
!
      use General
      use IO
      use Sub
!
      real, dimension (mpar_loc,mpvar) :: fp
      real :: dsnap_par_minor
      integer, dimension (mpar_loc) :: ipar
      integer :: npar_loc
      logical :: enum, lsnap
      character (len=*) :: snapbase, flist
!
      integer, save :: ifirst=0, nsnap, nsnap_minor
      real, save :: tsnap,tsnap_minor
      logical :: lsnap_minor=.false.
      character (len=fnlen) :: snapname, filename_diag
      character (len=4) :: nsnap_ch,nsnap_minor_ch,nsnap_ch_last
!
      optional :: flist
!
!  Output snapshot with label in 'tsnap' time intervals
!  file keeps the information about number and time of last snapshot
!
      if (enum) then
        call safe_character_assign(filename_diag,trim(datadir)//'/tsnap.dat')
!
!  at first call, need to initialize tsnap
!  tsnap calculated in read_snaptime, but only available to root processor
!
        if (ifirst==0) then
          call read_snaptime(filename_diag,tsnap,nsnap,dsnap,t)
          ifirst=1
          tsnap_minor=dsnap_par_minor
          nsnap_minor=1
        endif
!
!  Possible to output minor particle snapshots (e.g. for a movie).
!
        if (dsnap_par_minor/=0.) &
            call update_snaptime(filename_diag,tsnap_minor,nsnap_minor, &
            dsnap_par_minor,t,lsnap_minor,nsnap_minor_ch,ENUM=.true.)
        if (lsnap_minor) then
          call chn(nsnap-1,nsnap_ch_last,'')
          snapname=snapbase//trim(nsnap_ch_last)//'.'//trim(nsnap_minor_ch)
          call boundconds_particles(fp,npar_loc,ipar)
          call output_particles(snapname,fp,npar_loc,ipar)
          if(ip<=10 .and. lroot) &
              print*,'wsnap_particles: written snapshot ', snapname
          if (present(flist)) call log_filename_to_file(snapname,flist)
        endif
!
!  Regular data snapshots must come synchronized with the fluid snapshots.
!
        call update_snaptime(filename_diag,tsnap,nsnap,dsnap,t,lsnap,nsnap_ch, &
            ENUM=.true.)
        if (lsnap) then
          snapname=snapbase//nsnap_ch
          call boundconds_particles(fp,npar_loc,ipar)
          call output_particles(snapname,fp,npar_loc,ipar)
          if(ip<=10 .and. lroot) &
              print*,'wsnap_particles: written snapshot ', snapname
          if (present(flist)) call log_filename_to_file(snapname,flist)
          nsnap_minor=1
        endif
!
      else
!
!  Write snapshot without label
!
        snapname=snapbase
        call boundconds_particles(fp,npar_loc,ipar)
        call output_particles(snapname,fp,npar_loc,ipar)
        if(ip<=10 .and. lroot) &
             print*,'wsnap_particles: written snapshot ', snapname
        if (present(flist)) call log_filename_to_file(snapname,flist)
      endif
!
    endsubroutine wsnap_particles
!***********************************************************************
    subroutine boundconds_particles(fp,npar_loc,ipar,dfp)
!
!  Global boundary conditions for particles.
!
!  30-dec-04/anders: coded
!
      use Messages, only: fatal_error_local
      use Mpicomm
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mpar_loc,mpvar), optional :: dfp
      integer, dimension (mpar_loc) :: ipar
      integer :: npar_loc
!
      real :: xold,yold,rad,r1old
!
      integer :: k
!
      intent (inout) :: fp, npar_loc, ipar, dfp

      if (.not.lcartesian_mig) then
!     radial boundary condition
         do k=1,npar_loc
            xold=fp(k,ixp) ; yold=fp(k,iyp)
            rad = sqrt(fp(k,ixp)**2 + fp(k,iyp)**2)
            r1old = 1./max(rad,tini)
!      rp < r_int
            if (rad .lt. r_int) then
               rad=rad+(r_ext-r_int)
!
!  Particle position must never need more than one addition of Lr to get back
!  in the box. Often a NaN or Inf in the particle position will show up as a
!  problem here.
!
               if (rad .lt. r_int) then
                  print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                       ' was further than r_ext outside the simulation box!'
                  print*, 'This must never happen.'
                  print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                  call fatal_error_local('boundconds_particles','')
               endif
            endif
!      rp > r_ext
            if (rad >= r_ext) then
               rad=rad-(r_ext-r_int)
               if (rad>=r_ext) then
                  print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                       ' was further than r_ext outside the simulation box!'
                  print*, 'This must never happen.'
                  print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                  call fatal_error_local('boundconds_particles','')
               endif
            endif
!
            fp(k,ixp) = rad *xold*r1old !r*cos(theta)
            fp(k,iyp) = rad *yold*r1old !r*sin(theta)
!
         enddo
      else
!
!  Boundary condition in the x-direction.
!
      if (nxgrid/=1) then
        if (bcpx=='p') then
          do k=1,npar_loc
!  xp < x0
            if (fp(k,ixp)< xyz0(1)) then
              fp(k,ixp)=fp(k,ixp)+Lxyz(1)
              if (lshear.and.nygrid/=1) fp(k,iyp)=fp(k,iyp)-deltay
!  Particle position must never need more than one addition of Lx to get back
!  in the box. Often a NaN or Inf in the particle position will show up as a
!  problem here.
              if (fp(k,ixp)< xyz0(1)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                    ' was further than Lx outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
!  xp > x1
            if (fp(k,ixp)>=xyz1(1)) then
              fp(k,ixp)=fp(k,ixp)-Lxyz(1)
              if (lshear.and.nygrid/=1) fp(k,iyp)=fp(k,iyp)+deltay
              if (fp(k,ixp)>=xyz1(1)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                    ' was further than Lx outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
          enddo
        else
          print*, 'boundconds_particles: No such boundary condition bcpx=', bcpx
          call stop_it('boundconds_particles')
        endif
      endif
!
!  Boundary condition in the y-direction.
!
      if (nygrid/=1) then
        if (bcpy=='p') then
!  yp < y0
          do k=1,npar_loc
            if (fp(k,iyp)< xyz0(2)) then
              fp(k,iyp)=fp(k,iyp)+Lxyz(2)
              if (fp(k,iyp)< xyz0(2)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                    ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
!  yp > y1
            if (fp(k,iyp)>=xyz1(2)) then
              fp(k,iyp)=fp(k,iyp)-Lxyz(2)
              if (fp(k,iyp)>=xyz1(2)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                    ' was further than Ly outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
          enddo
        else
          print*, 'boundconds_particles: No such boundary condition bcpy=', bcpy
          call stop_it('boundconds_particles')
        endif
      endif
!
!  Boundary condition in the z-direction.
!
      if (nzgrid/=1) then
        if (bcpz=='p') then
          do k=1,npar_loc
!  zp < z0
            if (fp(k,izp)< xyz0(3)) then
              fp(k,izp)=fp(k,izp)+Lxyz(3)
              if (fp(k,izp)< xyz0(3)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                    ' was further than Lz outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
!  zp > z1
            if (fp(k,izp)>=xyz1(3)) then
              fp(k,izp)=fp(k,izp)-Lxyz(3)
              if (fp(k,izp)>=xyz1(3)) then
                print*, 'boundconds_particles: ERROR - particle ', ipar(k), &
                    ' was further than Lz outside the simulation box!'
                print*, 'This must never happen.'
                print*, 'iproc, ipar, xxp=', iproc, ipar(k), fp(k,ixp:izp)
                call fatal_error_local('boundconds_particles','')
              endif
            endif
          enddo
        else
          print*, 'boundconds_particles: No such boundary condition bcpz=', bcpz
          call stop_it('boundconds_particles')
        endif
      endif

   endif

!
!  Redistribute particles among processors (internal boundary conditions).
!
      if (lmpicomm) then
         if (present(dfp)) then
            call redist_particles_procs(fp,npar_loc,ipar,dfp)
         else
            call redist_particles_procs(fp,npar_loc,ipar)
         endif
      endif
!
    endsubroutine boundconds_particles
!***********************************************************************
    subroutine redist_particles_procs(fp,npar_loc,ipar,dfp)
!
!  Redistribute particles among processors based on the local yz-interval
!  of each processor.
!
!  01-jan-05/anders: coded
!
      use Messages, only: fatal_error, warning
      use Mpicomm
      use Sub, only: max_name
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mpar_loc,mpvar), optional :: dfp
      integer, dimension (mpar_loc) :: ipar
      integer :: npar_loc
!
      real, dimension (0:ncpus-1,npar_mig,mpvar) :: fp_mig, dfp_mig
      integer, dimension (0:ncpus-1,npar_mig) :: ipar_mig
      integer, dimension (0:ncpus-1,0:ncpus-1) :: nmig
      integer :: i, j, k, iproc_rec, ipy_rec, ipz_rec
      logical :: lredo, lredo_all
!
      intent (inout) :: fp, npar_loc, ipar, dfp
!
!  Possible to iterate untill all particles have migrated.
!
      lredo=.false.; lredo_all=.true.
      do while (lredo_all)
        lredo=.false.
!
!  Find out which particles are not in the local processor's yz-interval.
!
        nmig=0
        do k=npar_loc,1,-1
!  Find y index of receiving processor.
          ipy_rec=ipy
          if (fp(k,iyp)>=xyz1_loc(2)) then
            do while ( (fp(k,iyp)-(ipy_rec+1)*Lxyz_loc(2)) >= xyz0(2) )
              ipy_rec=ipy_rec+1
            enddo
          else if (fp(k,iyp)< xyz0_loc(2)) then
            do while ( (fp(k,iyp)+(nprocy-ipy_rec)*Lxyz_loc(2)) <  xyz1(2) )
              ipy_rec=ipy_rec-1
            enddo
          endif
!  Find z index of receiving processor.
          ipz_rec=ipz
          if (fp(k,izp)>=xyz1_loc(3)) then
            do while ( (fp(k,izp)-(ipz_rec+1)*Lxyz_loc(3)) >= xyz0(3) )
              ipz_rec=ipz_rec+1
            enddo
          else if (fp(k,izp)< xyz0_loc(3)) then
            do while ( (fp(k,izp)+(nprocz-ipz_rec)*Lxyz_loc(3)) <  xyz1(3) )
              ipz_rec=ipz_rec-1
            enddo
          endif
!  Calculate serial index of receiving processor.
          iproc_rec=ipy_rec+nprocy*ipz_rec
!  Migrate particle if it is no longer at the current processor.
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
            endif
!  Copy migrating particle to the end of the fp array.
            nmig(iproc,iproc_rec)=nmig(iproc,iproc_rec)+1
            if (nmig(iproc,iproc_rec)>npar_mig) then
              print '(a,i3,a,i3,a)', &
                  'redist_particles_procs: too many particles migrating '// &
                  'from proc ', iproc, ' to proc ', iproc_rec
              print*, '                       (npar_mig=', npar_mig, ')'
              if (ip<=8) then
                print*, 'redist_particles_procs: iproc, npar_mig=', &
                    iproc, npar_mig
                print*, 'redist_particles_procs: nmig=', nmig(iproc,:)
              endif
              if (lmigration_redo) then
                print*, '                       Going to do one more '// &
                    'migration iteration!'
                print*, '                       (this is time consuming - '// &
                    'consider setting npar_mig'
                print*, '                        higher in cparam.local)'
                nmig(iproc,iproc_rec)=nmig(iproc,iproc_rec)-1
                lredo=.true.
                exit
              else
                call fatal_error('redist_particles_procs','')
              endif
            endif
            fp_mig(iproc_rec,nmig(iproc,iproc_rec),:)=fp(k,:)
            if (present(dfp)) &
                dfp_mig(iproc_rec,nmig(iproc,iproc_rec),:)=dfp(k,:)
            ipar_mig(iproc_rec,nmig(iproc,iproc_rec))=ipar(k)
!  Move the particle with the highest index number to the empty spot left by
!  the migrating particle.
            fp(k,:)=fp(npar_loc,:)
            if (present(dfp)) dfp(k,:)=dfp(npar_loc,:)
            ipar(k)=ipar(npar_loc)
!  Reduce the number of particles by one.
            npar_loc=npar_loc-1
          endif
        enddo
!
!  Print out information about number of migrating particles.
!
        if (ip<=8) print*, 'redist_particles_procs: iproc, nmigrate = ', &
            iproc, sum(nmig(iproc,:))
!
!  Diagnostic about number of migrating particles.
!  WARNING: in time-steps where snapshots are written, this diagnostic
!  parameter will be zero (quite confusing)!
!
        if (ldiagnos.and.(idiag_nmigmax/=0)) &
            call max_name(sum(nmig(iproc,:)),idiag_nmigmax)
!
!  Share information about number of migrating particles.
!
        do i=0,ncpus-1
          if (iproc/=i) then
            call mpirecv_int(nmig(i,iproc), 1, i, 111)
          else
            do j=0,ncpus-1
              if (iproc/=j) call mpisend_int(nmig(iproc,j), 1, j, 111)
            enddo
          endif
        enddo
!
!  Set to receive.
!
        do i=0,ncpus-1
          if (iproc/=i .and. nmig(i,iproc)/=0) then
            call mpirecv_real(fp(npar_loc+1:npar_loc+nmig(i,iproc),:), &
                (/nmig(i,iproc),mpvar/), i, 222)
            call mpirecv_int(ipar(npar_loc+1:npar_loc+nmig(i,iproc)), &
                nmig(i,iproc), i, 223)
            if (present(dfp)) &
                call mpirecv_real(dfp(npar_loc+1:npar_loc+nmig(i,iproc),:), &
                (/nmig(i,iproc),mpvar/), i, 333)
            if (ip<=6) then
              print*, 'redist_particles_proc: iproc, iproc_send=', iproc, i
              print*, 'redist_particles_proc: received fp=', &
                  fp(npar_loc+1:npar_loc+nmig(i,iproc),:)
              print*, 'redist_particles_proc: received ipar=', &
                  ipar(npar_loc+1:npar_loc+nmig(i,iproc))
              if (present(dfp)) print*, 'redist_particles_proc: received dfp=',&
                  dfp(npar_loc+1:npar_loc+nmig(i,iproc),:)
            endif
            npar_loc=npar_loc+nmig(i,iproc)
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
              if (iproc/=j .and. nmig(iproc,j)/=0) then
                call mpisend_real(fp_mig(j,1:nmig(iproc,j),:), &
                    (/nmig(iproc,j),mpvar/), j, 222)
                call mpisend_int(ipar_mig(j,1:nmig(iproc,j)), &
                    nmig(iproc,j), j, 223)
                if (present(dfp)) &
                    call mpisend_real(dfp_mig(j,1:nmig(iproc,j),:), &
                    (/nmig(iproc,j),mpvar/), j, 333)
                if (ip<=6) then
                  print*, 'redist_particles_proc: iproc, iproc_rec=', iproc, j
                  print*, 'redist_particles_proc: sent fp=', &
                      fp_mig(j,1:nmig(iproc,j),:)
                  print*, 'redist_particles_proc: sent ipar=', &
                      ipar_mig(j,1:nmig(iproc,j))
                  if (present(dfp)) &
                      print*, 'redist_particles_proc: sent dfp=', &
                      dfp_mig(j,1:nmig(iproc,j),:)
                endif
              endif
            enddo
          endif
        enddo
!
!  Sum up processors that have not had place to let all migrating particles go.
!
        if (lmigration_redo) then   !  5-10% slowdown of code
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
    subroutine dist_particles_evenly_procs(npar_loc,ipar)
!
!  Distribute particles evenly among processors.
!
!  05-jan-05/anders: coded
!
      integer :: npar_loc
      integer, dimension (mpar_loc) :: ipar
!
      integer :: i,k
      integer, dimension (0:ncpus-1) :: ipar1, ipar2
!
      intent (out) :: npar_loc,ipar
!
      do i=0,ncpus-1
        if (i==0) then
          ipar1(i)=1
        else
          ipar1(i)=ipar2(i-1)+1
        endif
        ipar2(i)=ipar1(i) + (npar/ncpus-1) + (ncpus-(i+1)+mod(npar,ncpus))/ncpus
      enddo
!
      if (lroot) then
        print*, 'dist_particles_evenly_procs: ipar1=', ipar1
        print*, 'dist_particles_evenly_procs: ipar2=', ipar2
      endif
!
!  Set index interval of particles that belong to the local processor.
!
      npar_loc=ipar2(iproc)-ipar1(iproc)+1
      do k=1,npar_loc
        ipar(k)=k-1+ipar1(iproc)
      enddo
!
    endsubroutine dist_particles_evenly_procs
!***********************************************************************
    subroutine sum_par_name(a,iname,lsqrt)
!
!  Successively calculate sum of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  02-jan-05/anders: adapted from sum_mn_name
!
      use Cdata
      use Messages, only: fatal_error
!
      real, dimension (:) :: a
      integer :: iname
      logical, optional :: lsqrt
!
      integer, save :: icount=0
!
      if (iname/=0) then
!
        if (icount==0) fname(iname)=0
!
        fname(iname)=fname(iname)+sum(a)
!
!  Set corresponding entry in itype_name
!
        if (present(lsqrt)) then
          itype_name(iname)=ilabel_sum_sqrt_par
        else
          itype_name(iname)=ilabel_sum_par
        endif
!
!  Reset sum when npar_loc particles have been considered.
!
        icount=icount+size(a)
        if (icount==npar_loc) then
          icount=0
        elseif (icount>=npar_loc) then
          print*, 'sum_par_name: Too many particles entered this sub.'
          print*, 'sum_par_name: Can only do statistics on npar_loc particles!'
          call fatal_error('sum_par_name','')
        endif
!
      endif
!
    endsubroutine sum_par_name
!***********************************************************************
    subroutine sum_par_name_nw(a,iname,lsqrt)
!
!  successively calculate sum of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  22-aug-05/anders: adapted from sum_par_name
!
      use Cdata
      use Messages, only: fatal_error
!
      real, dimension (:) :: a
      integer :: iname
      logical, optional :: lsqrt
!
      integer, save :: icount=0
!
      if (iname/=0) then
!
        if (icount==0) fname(iname)=0
!
        fname(iname)=fname(iname)+sum(a)
!
!  Set corresponding entry in itype_name
!
        if (present(lsqrt)) then
          itype_name(iname)=ilabel_sum_sqrt
        else
          itype_name(iname)=ilabel_sum
        endif
!
        icount=icount+size(a)
        if (icount==nw) then
          icount=0
        elseif (icount>=nw) then
          print*, 'sum_par_name_nw: Too many grid points entered this sub.'
          print*, 'sum_par_name_nw: Can only do statistics on nw grid points!'
          call fatal_error('sum_par_name_nw','')
        endif
!
      endif
!
    endsubroutine sum_par_name_nw
!***********************************************************************
    subroutine max_par_name(a,iname,lneg)
!
!  Successively calculate maximum of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  28-nov-05/anders: adapted from max_par_name
!
      use Cdata
!
      real, dimension (:) :: a
      integer :: iname
      logical, optional :: lneg
!
      if (iname/=0) then
!
        fname(iname)=0.
        fname(iname)=fname(iname)+maxval(a)
!
!  set corresponding entry in itype_name
!
        if (present(lneg)) then
          itype_name(iname)=ilabel_max_neg
        else
          itype_name(iname)=ilabel_max
        endif
!
      endif
!
    endsubroutine max_par_name
!***********************************************************************
    subroutine integrate_par_name(a,iname)
!
!  Calculate integral of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  29-nov-05/anders: adapted from sum_par_name
!
      use Cdata
      use Messages, only: fatal_error
!
      real, dimension (:) :: a
      integer :: iname
!
      integer, save :: icount=0
!
      if (iname/=0) then
!
        if (icount==0) fname(iname)=0
!
        fname(iname)=fname(iname)+sum(a)
!
!  Set corresponding entry in itype_name
!
        itype_name(iname)=ilabel_integrate
!
!  Reset sum when npar_loc particles have been considered.
!
        icount=icount+size(a)
        if (icount==npar_loc) then
          icount=0
        elseif (icount>=npar_loc) then
          print*, 'integral_par_name: Too many particles entered this sub.'
          print*, 'integral_par_name: Can only do statistics on npar_loc particles!'
          call fatal_error('integral_par_name','')
        endif
!
      endif
!
    endsubroutine integrate_par_name
!***********************************************************************
    subroutine interpolate_linear(f,ivar1,ivar2,xxp,gp,inear,ipar)
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
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (3) :: xxp
      integer, dimension (3) :: inear
      integer :: ivar1, ivar2
      real, dimension (ivar2-ivar1+1) :: gp
      integer, optional :: ipar
!
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8
      real :: xp0, yp0, zp0
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      integer :: i, ix0, iy0, iz0
      logical :: lfirstcall=.true.
!
      intent(in)  :: f, xxp, ivar1
      intent(out) :: gp
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      if ( (x(ix0)>xxp(1)) .and. nxgrid/=1) ix0=ix0-1
      if ( (y(iy0)>xxp(2)) .and. nygrid/=1) iy0=iy0-1
      if ( (z(iz0)>xxp(3)) .and. nzgrid/=1) iz0=iz0-1
!
!  Check if the grid point interval is really correct.
!
      if ((x(ix0)<=xxp(1) .and. x(ix0+1)>=xxp(1) .or. nxgrid==1) .and. &
          (y(iy0)<=xxp(2) .and. y(iy0+1)>=xxp(2) .or. nygrid==1) .and. &
          (z(iz0)<=xxp(3) .and. z(iz0+1)>=xxp(3) .or. nzgrid==1)) then
        ! Everything okay
      else
        print*, 'interpolate_linear: Interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        if (present(ipar)) print*, 'ipar= ', ipar
        print*, 'mx, x(1), x(mx) = ', mx, x(1), x(mx)
        print*, 'my, y(1), y(my) = ', my, y(1), y(my)
        print*, 'mz, z(1), z(mz) = ', mz, z(1), z(mz)
        print*, 'ix0, iy0, iz0 = ', ix0, iy0, iz0
        print*, 'xp, xp0, xp1 = ', xxp(1), x(ix0), x(ix0+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), y(iy0), y(iy0+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), z(iz0), z(iz0+1)
        call stop_it('interpolate_linear')
      endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!
      xp0=xxp(1)-x(ix0)
      yp0=xxp(2)-y(iy0)
      zp0=xxp(3)-z(iz0)
!
!  Calculate derived grid spacing variables in the first call to this sub.
!
      if (lfirstcall) then
        dx1=1/dx;  dy1=1/dy;  dz1=1/dz
        dxdy1=1/(dx*dy);  dxdz1=1/(dx*dz);  dydz1=1/(dy*dz)
        dxdydz1=1/(dx*dy*dz)
        lfirstcall=.false.
      endif
!
!  Function values at all corners.
!
      g1=f(ix0  ,iy0  ,iz0  ,ivar1:ivar2)
      g2=f(ix0+1,iy0  ,iz0  ,ivar1:ivar2)
      g3=f(ix0  ,iy0+1,iz0  ,ivar1:ivar2)
      g4=f(ix0+1,iy0+1,iz0  ,ivar1:ivar2)
      g5=f(ix0  ,iy0  ,iz0+1,ivar1:ivar2)
      g6=f(ix0+1,iy0  ,iz0+1,ivar1:ivar2)
      g7=f(ix0  ,iy0+1,iz0+1,ivar1:ivar2)
      g8=f(ix0+1,iy0+1,iz0+1,ivar1:ivar2)
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
                x(ix0), y(iy0), z(iz0)
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
                x(ix0), y(iy0), z(iz0)
            print*, 'interpolate_linear: i, gp(i)=', i, gp(i)
            print*, 'interpolate_linear: g1...g8=', &
                g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i)
            print*, '------------------'
          endif
        enddo
      endif
!
    endsubroutine interpolate_linear
!***********************************************************************
    subroutine interpolate_quadratic(f,ivar1,ivar2,xxp,gp,inear,ipar)
!
!  Quadratic interpolation of g to arbitrary (xp, yp, zp) coordinate
!  using the biquadratic interpolation function
!
!    g(x,y,z) = (1+x+x^2)*(1+z+z^2)
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
      use Cdata
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (3) :: xxp
      integer, dimension (3) :: inear
      integer :: ivar1, ivar2
      real, dimension (ivar2-ivar1+1) :: gp
      integer, optional :: ipar
!
      real, dimension (9,ivar2-ivar1+1) :: cc
      real, dimension (ivar2-ivar1+1) :: g1, g2, g3, g4, g5, g6, g7, g8, g9
      real :: dxp, dzp
      real, save :: dx1, dx2, dz1, dz2, dx1dz1, dx2dz1, dx1dz2, dx2dz2
      integer :: ix0, iy0, iz0
      logical, save :: lfirstcall=.true.
!
      intent(in)  :: f, xxp, ivar1
      intent(out) :: gp
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
!
!  Not implemented in y-direction yet (but is easy to generalise).
!
      if (nygrid/=1) then
        if (lroot) print*, 'interpolate_quadratic: not implemented in y'
        call fatal_error('interpolate_quadratic','')
      endif
!
!  A few values that only need to be calcultad once.
!
      if (lfirstcall) then
        dx1=1/dx; dx2=1/dx**2
        dz1=1/dz; dz2=1/dz**2
        dx1dz1=1/(dx*dz)
        dx2dz1=1/(dx**2*dz); dx1dz2=1/(dx*dz**2); dx2dz2=1/(dx**2*dz**2)
        lfirstcall=.false.
      endif
!
!  Define function values at the grid points.
!
      g1=f(ix0-1,iy0,iz0-1,ivar1:ivar2)
      g2=f(ix0  ,iy0,iz0-1,ivar1:ivar2)
      g3=f(ix0+1,iy0,iz0-1,ivar1:ivar2)
      g4=f(ix0-1,iy0,iz0  ,ivar1:ivar2)
      g5=f(ix0  ,iy0,iz0  ,ivar1:ivar2)
      g6=f(ix0+1,iy0,iz0  ,ivar1:ivar2)
      g7=f(ix0-1,iy0,iz0+1,ivar1:ivar2)
      g8=f(ix0  ,iy0,iz0+1,ivar1:ivar2)
      g9=f(ix0+1,iy0,iz0+1,ivar1:ivar2)
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
      dxp=xxp(1)-x(ix0)
      dzp=xxp(3)-z(iz0)
!
      gp = cc(1,:)            + cc(2,:)*dxp        + cc(3,:)*dxp**2        + &
           cc(4,:)*dzp        + cc(5,:)*dzp**2     + cc(6,:)*dxp*dzp       + &
           cc(7,:)*dxp**2*dzp + cc(8,:)*dxp*dzp**2 + cc(9,:)*dxp**2*dzp**2
!
      if (NO_WARN) print*, ipar
!
    endsubroutine interpolate_quadratic
!***********************************************************************
    subroutine interpolate_quadratic_spline(f,ivar1,ivar2,xxp,gp,inear,ipar)
!
!  Quadratic spline interpolation of the function g to the point xxp=(xp,yp,zp).
!
!  10-jun-06/anders: coded
!
      use Cdata
      use Messages, only: fatal_error
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (3) :: xxp
      integer, dimension (3) :: inear
      integer :: ivar1, ivar2
      real, dimension (ivar2-ivar1+1) :: gp
      integer, optional :: ipar
!
      real :: fac_x_m1, fac_x_00, fac_x_p1
      real :: fac_y_m1, fac_y_00, fac_y_p1
      real :: fac_z_m1, fac_z_00, fac_z_p1
      real :: dxp0, dyp0, dzp0
      integer :: ix0, iy0, iz0
!
      intent(in)  :: f, xxp, ivar1
      intent(out) :: gp
!
!  Redefine the interpolation point in coordinates relative to nearest grid
!  point and normalize with the cell size.
!
      ix0=inear(1); iy0=inear(2); iz0=inear(3)
      dxp0=(xxp(1)-x(ix0))*dx_1(ix0)
      dyp0=(xxp(2)-y(iy0))*dy_1(iy0)
      dzp0=(xxp(3)-z(iz0))*dz_1(iz0)
!
!  Interpolation formulae.
!
      if (dimensionality==0) then
        gp=f(ix0,iy0,iz0,ivar1:ivar2)
      elseif (dimensionality==1) then
        if (nxgrid/=1) then
          gp = 0.5*(0.5-dxp0)**2*f(ix0-1,iy0,iz0,ivar1:ivar2) + &
                  (0.75-dxp0**2)*f(ix0  ,iy0,iz0,ivar1:ivar2) + &
               0.5*(0.5+dxp0)**2*f(ix0+1,iy0,iz0,ivar1:ivar2)
        endif
        if (nygrid/=1) then
          gp = 0.5*(0.5-dyp0)**2*f(ix0,iy0-1,iz0,ivar1:ivar2) + &
                  (0.75-dyp0**2)*f(ix0,iy0  ,iz0,ivar1:ivar2) + &
               0.5*(0.5+dyp0)**2*f(ix0,iy0+1,iz0,ivar1:ivar2)
        endif
        if (nzgrid/=1) then
          gp = 0.5*(0.5-dzp0)**2*f(ix0,iy0,iz0-1,ivar1:ivar2) + &
                  (0.75-dzp0**2)*f(ix0,iy0,iz0  ,ivar1:ivar2) + &
               0.5*(0.5+dzp0)**2*f(ix0,iy0,iz0+1,ivar1:ivar2)
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
          gp= fac_y_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_y_00*( f(ix0,iy0  ,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0  ,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_z_00*( f(ix0,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
              fac_y_p1*( f(ix0,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_y_m1*( f(ix0,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 )
        elseif (nygrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_z_m1 = 0.5*(0.5-dzp0)**2
          fac_z_00 = 0.75-dzp0**2
          fac_z_p1 = 0.5*(0.5+dzp0)**2
!
          gp= fac_x_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_x_00*( f(ix0  ,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0  ,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_z_00*( f(ix0+1,iy0,iz0  ,ivar1:ivar2)*fac_x_p1 + &
                         f(ix0-1,iy0,iz0  ,ivar1:ivar2)*fac_x_m1 ) + &
              fac_x_p1*( f(ix0+1,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0+1,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
              fac_x_m1*( f(ix0-1,iy0,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                         f(ix0-1,iy0,iz0-1,ivar1:ivar2)*fac_z_m1 )
        elseif (nzgrid==1) then
          fac_x_m1 = 0.5*(0.5-dxp0)**2
          fac_x_00 = 0.75-dxp0**2
          fac_x_p1 = 0.5*(0.5+dxp0)**2
          fac_y_m1 = 0.5*(0.5-dyp0)**2
          fac_y_00 = 0.75-dyp0**2
          fac_y_p1 = 0.5*(0.5+dyp0)**2
!
          gp= fac_x_00*fac_y_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
              fac_x_00*( f(ix0  ,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0  ,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 ) + &
              fac_y_00*( f(ix0+1,iy0  ,iz0,ivar1:ivar2)*fac_x_p1 + &
                         f(ix0-1,iy0  ,iz0,ivar1:ivar2)*fac_x_m1 ) + &
              fac_x_p1*( f(ix0+1,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0+1,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 ) + &
              fac_x_m1*( f(ix0-1,iy0+1,iz0,ivar1:ivar2)*fac_y_p1 + &
                         f(ix0-1,iy0-1,iz0,ivar1:ivar2)*fac_y_m1 )
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
        gp= fac_x_00*fac_y_00*fac_z_00*f(ix0,iy0,iz0,ivar1:ivar2) + &
            fac_x_00*fac_y_00*( f(ix0  ,iy0  ,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0  ,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_z_00*( f(ix0  ,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0  ,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
            fac_y_00*fac_z_00*( f(ix0+1,iy0  ,iz0  ,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0  ,ivar1:ivar2)*fac_x_m1 ) + &
            fac_x_p1*fac_y_p1*( f(ix0+1,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0+1,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_p1*fac_y_m1*( f(ix0+1,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0+1,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_m1*fac_y_p1*( f(ix0-1,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0-1,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_m1*fac_y_m1*( f(ix0-1,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0-1,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_y_p1*( f(ix0  ,iy0+1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0+1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_x_00*fac_y_m1*( f(ix0  ,iy0-1,iz0+1,ivar1:ivar2)*fac_z_p1 + &
                                f(ix0  ,iy0-1,iz0-1,ivar1:ivar2)*fac_z_m1 ) + &
            fac_y_00*fac_z_p1*( f(ix0+1,iy0  ,iz0+1,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0+1,ivar1:ivar2)*fac_x_m1 ) + &
            fac_y_00*fac_z_m1*( f(ix0+1,iy0  ,iz0-1,ivar1:ivar2)*fac_x_p1 + &
                                f(ix0-1,iy0  ,iz0-1,ivar1:ivar2)*fac_x_m1 ) + &
            fac_z_00*fac_x_p1*( f(ix0+1,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0+1,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 ) + &
            fac_z_00*fac_x_m1*( f(ix0-1,iy0+1,iz0  ,ivar1:ivar2)*fac_y_p1 + &
                                f(ix0-1,iy0-1,iz0  ,ivar1:ivar2)*fac_y_m1 )
      endif
!
      if (NO_WARN) print*, ipar
!
    endsubroutine interpolate_quadratic_spline
!***********************************************************************
    subroutine map_nearest_grid(fp,ineargrid)
!
!  Find index (ix0, iy0, iz0) of nearest grid point of all the particles.
!
!  23-jan-05/anders: coded
!
      use Cdata
      use Mpicomm, only: stop_it
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real, save :: dx1, dy1, dz1
      integer :: k, ix0, iy0, iz0
      logical, save :: lfirstcall=.true.
!
      intent(in)  :: fp
      intent(out) :: ineargrid
!
!  Default values in case of missing directions.
!
      ix0=nghost+1; iy0=nghost+1; iz0=nghost+1
!
      if (lfirstcall) then
        if (.not. all(lequidist)) then
          print*, 'map_xx_vv_grid: only works for equidistant grid!'
          call stop_it('map_xx_vv_grid')
        endif
        dx1=1/dx; dy1=1/dy; dz1=1/dz
        lfirstcall=.false.
      endif
!
      do k=1,npar_loc
        if (nxgrid/=1) ix0 = nint((fp(k,ixp)-x(1))*dx1) + 1
        if (nygrid/=1) iy0 = nint((fp(k,iyp)-y(1))*dy1) + 1
        if (nzgrid/=1) iz0 = nint((fp(k,izp)-z(1))*dz1) + 1
        ineargrid(k,1)=ix0; ineargrid(k,2)=iy0; ineargrid(k,3)=iz0
      enddo
!
    endsubroutine map_nearest_grid
!***********************************************************************
    subroutine sort_particles_imn(fp,ineargrid,ipar,dfp)
!
!  Sort the particles so that they appear in the same order as the (m,n) loop.
!
!  20-apr-06/anders: coded
!
      use Cdata
      use General, only: safe_character_assign
      use Messages, only: fatal_error
!
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (mpar_loc) :: ipar
      real, dimension (mpar_loc,mpvar), optional :: dfp
!
      real, dimension (mpvar) :: fp_tmp, dfp_tmp
      integer, dimension (3) :: ineargrid_tmp
      integer, dimension (ny*nz) :: kk
      integer :: ilmn_par_tmp, ipark_sorted_tmp, ipar_tmp
      integer, dimension (mpar_loc) :: ilmn_par, ipark_sorted
      integer :: j, k, ix0, iy0, iz0, ih, lun
      integer, save :: ncount=-1, isorttype=3
      integer, parameter :: nshellsort=21
      integer, dimension (nshellsort), parameter :: &
          hshellsort=(/ 14057, 9371, 6247, 4177, 2777, 1861, 1237, 823, 557, &
                        367, 251, 163, 109, 73, 37, 19, 11, 7, 5, 3, 1 /)
      logical, save :: lrunningsort=.false.
      character (len=fnlen) :: filename
!
      intent(inout)  :: fp, ineargrid, ipar, dfp
!
!  Determine beginning and ending index of particles in pencil (m,n).
!
      call particle_pencil_index(ineargrid)
!
!  Choose sort algorithm.
!  WARNING - choosing the wrong one might make the code unnecessarily slow.
!
      if ( lstart .or. &
          (lshear .and. Omega/=0.0) .and. (nxgrid>1 .and. nygrid>1) ) then
        isorttype=3
        lrunningsort=.false.
      else
        isorttype=1
        lrunningsort=.true.
      endif
      ncount=0
!
!  Calculate integer value to sort after.
!
      do k=1,npar_loc
        ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
        ilmn_par(k)=imn_array(iy0,iz0)!-1)*ny*nz+ix0
        ipark_sorted(k)=k
      enddo
!
!  Sort using either straight insertion (1), shell sorting (2) or counting
!  sort (3).
!
      select case (isorttype)
!  Straight insertion.
      case (1)
        do k=2,npar_loc

          j=k

          do while ( ilmn_par(k)<ilmn_par(j-1) )
            j=j-1
            if (j==1) exit
          enddo

          if (j/=k) then
            ncount=ncount+k-j
!
            ilmn_par_tmp=ilmn_par(k)
            ilmn_par(j+1:k)=ilmn_par(j:k-1)
            ilmn_par(j)=ilmn_par_tmp
            ipark_sorted_tmp=ipark_sorted(k)
            ipark_sorted(j+1:k)=ipark_sorted(j:k-1)
            ipark_sorted(j)=ipark_sorted_tmp
!  Sort particle data on the fly (practical for almost ordered arrays).
            if (lrunningsort) then
              fp_tmp=fp(k,:)
              fp(j+1:k,:)=fp(j:k-1,:)
              fp(j,:)=fp_tmp
!
              if (present(dfp)) then
                dfp_tmp=dfp(k,:)
                dfp(j+1:k,:)=dfp(j:k-1,:)
                dfp(j,:)=dfp_tmp
              endif
!
              ineargrid_tmp=ineargrid(k,:)
              ineargrid(j+1:k,:)=ineargrid(j:k-1,:)
              ineargrid(j,:)=ineargrid_tmp
!
              ipar_tmp=ipar(k)
              ipar(j+1:k)=ipar(j:k-1)
              ipar(j)=ipar_tmp
!
            endif
          endif
        enddo
!  Shell sort.
      case (2)

        do ih=1,21
          do k=1+hshellsort(ih),npar_loc
            ilmn_par_tmp=ilmn_par(k)
            ipark_sorted_tmp=ipark_sorted(k)
            j=k
            do while (ilmn_par(j-hshellsort(ih)) > ilmn_par_tmp)
              ncount=ncount+1
              ilmn_par(j)=ilmn_par(j-hshellsort(ih))
              ilmn_par(j-hshellsort(ih))=ilmn_par_tmp
              ipark_sorted(j)=ipark_sorted(j-hshellsort(ih))
              ipark_sorted(j-hshellsort(ih))=ipark_sorted_tmp
              j=j-hshellsort(ih)
              if (j-hshellsort(ih)<1) exit
            enddo
          enddo
        enddo
!  Counting sort.
      case(3)

        kk=k1_imn
        do k=1,npar_loc
          ipark_sorted(kk(ilmn_par(k)))=k
          kk(ilmn_par(k))=kk(ilmn_par(k))+1
        enddo
        ncount=npar_loc

      endselect
!
!  Sort particle data according to sorting index.
!
      if (lrunningsort .and. isorttype/=1) then
        if (lroot) print*, 'sort_particles_imn: lrunningsort is only '// &
            'allowed for straight insertion sort.'
        call fatal_error('sort_particles_imn','')
      endif
!
      if ( (.not. lrunningsort) .and. (ncount/=0) ) then
         fp(1:npar_loc,:)=fp(ipark_sorted(1:npar_loc),:)
         if (present(dfp)) dfp(1:npar_loc,:)=dfp(ipark_sorted(1:npar_loc),:)
         ineargrid(1:npar_loc,:)=ineargrid(ipark_sorted(1:npar_loc),:)
         ipar(1:npar_loc)=ipar(ipark_sorted(1:npar_loc))
      endif
!
      if (lroot.and.ldiagnos) then
        call safe_character_assign(filename,trim(datadir)//'/sort_particles.dat')
        lun=1
        open(lun,file=trim(filename),action='write',position='append')
        write(lun,'(A15,f7.3)') '------------ t=', t
        write(lun,'(A40,3i9,l9)')  'iproc, ncount, isorttype, lrunningsort=', &
            iproc, ncount, isorttype, lrunningsort
        close (lun)
      endif
!
      if (ip<=8) print '(A,i4,i8,i4,l4)', &
          'sort_particles_imn: iproc, ncount, isorttype, lrunningsort=', &
          iproc, ncount, isorttype, lrunningsort
!
    endsubroutine sort_particles_imn
!***********************************************************************
    subroutine map_xxp_grid(f,fp,ineargrid)
!
!  Calculate the number of particles in each grid cell.
!
!  27-nov-05/anders: coded
!
      use Cdata
      use GhostFold, only: fold_f
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
!
      real :: weight, weight_x, weight_y, weight_z
      integer :: k, ix0, iy0, iz0, ixx, iyy, izz
      integer :: ixx0, ixx1, iyy0, iyy1, izz0, izz1
!
      intent(in)  :: fp, ineargrid
      intent(out) :: f
!
!  Calculate the number of particles in each grid cell.
!
      if (inp/=0) then
        f(:,:,:,inp)=0.0
        do k=1,npar_loc
          ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
          f(ix0,iy0,iz0,inp) = f(ix0,iy0,iz0,inp) + 1.0
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
      if (irhop/=0) then
        f(:,:,:,irhop)=0.0
        if (lparticlemesh_cic) then
!
!  Cloud In Cell (CIC) scheme.
!
          do k=1,npar_loc
            ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
            ixx0=ix0; iyy0=iy0; izz0=iz0
            ixx1=ix0; iyy1=iy0; izz1=iz0
            if ( (x(ix0)>fp(k,ixp)) .and. nxgrid/=1) ixx0=ixx0-1
            if ( (y(iy0)>fp(k,iyp)) .and. nygrid/=1) iyy0=iyy0-1
            if ( (z(iz0)>fp(k,izp)) .and. nzgrid/=1) izz0=izz0-1
            if (nxgrid/=1) ixx1=ixx0+1
            if (nygrid/=1) iyy1=iyy0+1
            if (nzgrid/=1) izz1=izz0+1
            do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
              weight=1.0
              if (nxgrid/=1) &
                  weight=weight*( 1.0-abs(fp(k,ixp)-x(ixx))*dx_1(ixx) )
              if (nygrid/=1) &
                  weight=weight*( 1.0-abs(fp(k,iyp)-y(iyy))*dy_1(iyy) )
              if (nzgrid/=1) &
                  weight=weight*( 1.0-abs(fp(k,izp)-z(izz))*dz_1(izz) )
              f(ixx,iyy,izz,irhop)=f(ixx,iyy,izz,irhop) + weight
            enddo; enddo; enddo
          enddo
!
!  Triangular Shaped Cloud (TSC) scheme.
!
        elseif (lparticlemesh_tsc) then
!
!  Particle influences the 27 surrounding grid points, but has a density that
!  decreases with the distance from the particle centre.
!
          do k=1,npar_loc
            ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
            if (nxgrid/=1) then
              ixx0=ix0-1; ixx1=ix0+1
            else
              ixx0=ix0  ; ixx1=ix0
            endif
            if (nygrid/=1) then
              iyy0=iy0-1; iyy1=iy0+1
            else
              iyy0=iy0  ; iyy1=iy0
            endif
            if (nzgrid/=1) then
              izz0=iz0-1; izz1=iz0+1
            else
              izz0=iz0  ; izz1=iz0
            endif
!
!  The nearest grid point is influenced differently than the left and right
!  neighbours are. A particle that is situated exactly on a grid point gives
!  3/4 contribution to that grid point and 1/8 to each of the neighbours.
!
            do izz=izz0,izz1; do iyy=iyy0,iyy1; do ixx=ixx0,ixx1
              if ( ((ixx-ix0)==-1) .or. ((ixx-ix0)==+1) ) then
                weight_x = 1.125 - 1.5* abs(fp(k,ixp)-x(ixx))*dx_1(ixx) + &
                                   0.5*(abs(fp(k,ixp)-x(ixx))*dx_1(ixx))**2
              else
                if (nxgrid/=1) &
                weight_x = 0.75  -       ((fp(k,ixp)-x(ixx))*dx_1(ixx))**2
              endif
              if ( ((iyy-iy0)==-1) .or. ((iyy-iy0)==+1) ) then
                weight_y = 1.125 - 1.5* abs(fp(k,iyp)-y(iyy))*dy_1(iyy) + &
                                   0.5*(abs(fp(k,iyp)-y(iyy))*dy_1(iyy))**2
              else
                if (nygrid/=1) &
                weight_y = 0.75  -       ((fp(k,iyp)-y(iyy))*dy_1(iyy))**2
              endif
              if ( ((izz-iz0)==-1) .or. ((izz-iz0)==+1) ) then
                weight_z = 1.125 - 1.5* abs(fp(k,izp)-z(izz))*dz_1(izz) + &
                                   0.5*(abs(fp(k,izp)-z(izz))*dz_1(izz))**2
              else
                if (nzgrid/=1) &
                weight_z = 0.75  -       ((fp(k,izp)-z(izz))*dz_1(izz))**2
              endif

              weight=1.0

              if (nxgrid/=1) weight=weight*weight_x
              if (nygrid/=1) weight=weight*weight_y
              if (nzgrid/=1) weight=weight*weight_z
              f(ixx,iyy,izz,irhop)=f(ixx,iyy,izz,irhop) + weight
            enddo; enddo; enddo
          enddo
!
!  Nearest Grid Point (NGP) method.
!
        else
          f(l1:l2,m1:m2,n1:n2,irhop)=f(l1:l2,m1:m2,n1:n2,inp)
        endif
!
!  Fold first ghost zone of f.
!
        if (lparticlemesh_cic.or.lparticlemesh_tsc) call fold_f(f,irhop,irhop)
        f(l1:l2,m1:m2,n1:n2,irhop)=rhop_tilde*f(l1:l2,m1:m2,n1:n2,irhop)
!        call sharpen_tsc_density(f)
      endif
!
    endsubroutine map_xxp_grid
!***********************************************************************
    subroutine sharpen_tsc_density(f)
!
!  Sharpen density amplitudes (experimental).
!
!   9-nov-06/anders: coded
!
      use Fourier
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (nx,ny,nz) :: a1, b1
      integer :: ikx, iky, ikz
      real :: k2, k4
!
      a1=f(l1:l2,m1:m2,n1:n2,irhop)
      b1=0.0
      if (lshear) then
        call fourier_transform_shear(a1,b1)
      else
        call fourier_transform(a1,b1)
      endif
      do ikz=1,nz; do iky=1,ny; do ikx=1,nx
        k2 = kx_fft(ikx)**2 + ky_fft(iky+ipy*ny)**2 + kz_fft(ikz+ipz*nz)**2
        k4 = kx_fft(ikx)**4 + ky_fft(iky+ipy*ny)**4 + kz_fft(ikz+ipz*nz)**4
        a1(ikx,iky,ikz)=a1(ikx,iky,ikz)/(1-dx**2*k2/8+13*dx**4*k4/1920)
        b1(ikx,iky,ikz)=b1(ikx,iky,ikz)/(1-dx**2*k2/8+13*dx**4*k4/1920)
      enddo; enddo; enddo
      if (lshear) then
        call fourier_transform_shear(a1,b1,linv=.true.)
      else
        call fourier_transform(a1,b1,linv=.true.)
      endif
      f(l1:l2,m1:m2,n1:n2,irhop)=a1
!
    endsubroutine sharpen_tsc_density
!***********************************************************************
    subroutine particle_pencil_index(ineargrid)
!
!  Calculate the beginning and ending index of particles in a pencil.
!
!  24-apr-06/anders: coded
!
      use Cdata
      use Mpicomm, only: stop_it
!
      integer, dimension (mpar_loc,3) :: ineargrid
!
      integer :: k, ix0, iy0, iz0
!
      intent(in)  :: ineargrid
!
      npar_imn=0
!
!  Calculate the number of particles in each pencil.
!
      do k=1,npar_loc
        ix0=ineargrid(k,1); iy0=ineargrid(k,2); iz0=ineargrid(k,3)
        npar_imn(imn_array(iy0,iz0))=npar_imn(imn_array(iy0,iz0))+1
      enddo
!
!  Calculate beginning and ending particle index for each pencil.
!
      k=0
      do imn=1,ny*nz
        if (npar_imn(imn)/=0) then
          k1_imn(imn)=k + 1
          k2_imn(imn)=k1_imn(imn) + npar_imn(imn) - 1
          k=k+npar_imn(imn)
        else
          k1_imn(imn)=0
          k2_imn(imn)=0
        endif
      enddo
!
    endsubroutine particle_pencil_index
!***********************************************************************
    subroutine find_closest_gridpoint(xxp,ix0,iy0,iz0)
!
!  Find index (ix0, iy0, iz0) of closest grid point to the point xxp.
!
!  23-aug-05/anders: coded
!
      use Cdata
      use Messages, only: fatal_error
!
      real, dimension(3) :: xxp
      integer :: ix0, iy0, iz0
!
      real, save :: dx1, dy1, dz1
      logical, save :: lfirstcall=.true.
!
      intent(in)  :: xxp
      intent(out) :: ix0, iy0, iz0
!
      ix0=nghost+1; iy0=nghost+1; iz0=nghost+1
!
      if (lfirstcall) then
        if (.not. all(lequidist)) then
          print*, 'find_closeset_gridpoint: only works for equidistant grid!'
          call fatal_error('find_closest_gridpoint','')
        endif
        dx1=1/dx; dy1=1/dy; dz1=1/dz
        lfirstcall=.false.
      endif
!
      if (nxgrid/=1) ix0 = nint((xxp(1)-x(1))*dx1) + 1
      if (nygrid/=1) iy0 = nint((xxp(2)-y(1))*dy1) + 1
      if (nzgrid/=1) iz0 = nint((xxp(3)-z(1))*dz1) + 1
!
    endsubroutine find_closest_gridpoint
!***********************************************************************
    subroutine shepherd_neighbour(f,fp,ineargrid,kshepherd,kneighbour)
!
!  Create a shepherd/neighbour list of particles in the pencil.
!
!  24-oct-05/anders: coded
!
      use Cdata, only: nghost
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mpar_loc,mpvar) :: fp
      integer, dimension (mpar_loc,3) :: ineargrid
      integer, dimension (nx) :: kshepherd
      integer, dimension (:) :: kneighbour
!
      integer :: k, ix0
!
      kshepherd=0
      if (imn==1) kneighbour=0
!
      if (npar_imn(imn)/=0) then
        do k=k1_imn(imn),k2_imn(imn)
          ix0=ineargrid(k,1)
          kneighbour(k)=kshepherd(ix0-nghost)
          kshepherd(ix0-nghost)=k
        enddo
      endif
!
    endsubroutine shepherd_neighbour
!***********************************************************************
endmodule Particles_sub
