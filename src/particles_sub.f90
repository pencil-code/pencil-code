! $Id: particles_sub.f90,v 1.3 2005-07-01 17:45:29 ajohan Exp $
!
!  This module contains subroutines useful for the Particle module.
!
module Particles_sub

  use Cdata

  implicit none

  private
 
  public :: input_particles, output_particles
  public :: wsnap_particles, boundconds_particles
  public :: redist_particles_procs, dist_particles_evenly_procs
  public :: sum_par_name, find_lowest_cornerpoint, interpolate_3d_1st
  public :: map_xx_vv_grid

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
      integer :: k
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
        read(1) ipar(1:npar_loc)
!
!  Then read particle data.
!        
        if (npar_loc/=0) read(1) fp(1:npar_loc,:)
!
!  Read snapshot time.
!
        read(1) t
!        
        if (ip<=8) print*,'input_particles: read ',filename
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
      integer :: k
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
        write(lun_output) ipar
!
!  Then write particle data.
!
        if (npar_loc/=0) write(lun_output) fp(1:npar_loc,:)
!
!  Write time
!
        write(lun_output) t
!
      close(lun_output)
!
    endsubroutine output_particles
!***********************************************************************
    subroutine wsnap_particles(snapbase,fp,msnap,enum,lsnap,dsnap_par_minor, &
        npar_loc,ipar,flist)
!
!  Write particle snapshot file, labelled consecutively if enum==.true.
!  Otherwise just write a snapshot without label (used e.g. for pvar.dat)
!
!  29-dec-04/anders: adapted from wsnap
!
      use General
      use Io
      use Sub
!
      real, dimension (mpar_loc,mpvar) :: fp
      real :: dsnap_par_minor
      integer, dimension (mpar_loc) :: ipar
      integer :: msnap, npar_loc
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
      use Mpicomm
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mpar_loc,mpvar), optional :: dfp
      integer, dimension (mpar_loc) :: ipar
      integer :: npar_loc
!
      integer :: k
!
      intent (inout) :: fp, npar_loc, ipar, dfp
!
!  Boundary condition in the x-direction.
!
      if (bcpx=='p') then
        do k=1,npar_loc
          do while (fp(k,ixp)< xyz0(1))
            fp(k,ixp)=fp(k,ixp)+Lxyz(1)
            if (lshear) fp(k,iyp)=fp(k,iyp)-deltay
          enddo
          do while (fp(k,ixp)>=xyz1(1))
            fp(k,ixp)=fp(k,ixp)-Lxyz(1)
            if (lshear) fp(k,iyp)=fp(k,iyp)+deltay
          enddo
        enddo
      else
        print*, 'boundconds_particles: No such boundary condition bcpx=', bcpx
        call stop_it('boundconds_particles')
      endif
!
!  Boundary condition in the y-direction.
!
      if (bcpy=='p') then
        do k=1,npar_loc
          do while (fp(k,iyp)< xyz0(2))
            fp(k,iyp)=fp(k,iyp)+Lxyz(2)
          enddo
          do while (fp(k,iyp)>=xyz1(2))
            fp(k,iyp)=fp(k,iyp)-Lxyz(2)
          enddo
        enddo
      else
        print*, 'boundconds_particles: No such boundary condition bcpy=', bcpy
        call stop_it('boundconds_particles')
      endif
!
!  Boundary condition in the z-direction.
!
      if (bcpz=='p') then
        do k=1,npar_loc
          do while (fp(k,izp)< xyz0(3))
            fp(k,izp)=fp(k,izp)+Lxyz(3)
          enddo
          do while (fp(k,izp)>=xyz1(3))
            fp(k,izp)=fp(k,izp)-Lxyz(3)
          enddo
        enddo
      else
        print*, 'boundconds_particles: No such boundary condition bcpz=', bcpz
        call stop_it('boundconds_particles')
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
      use Messages, only: fatal_error
!
      include 'mpif.h'
!
      real, dimension (mpar_loc,mpvar) :: fp
      real, dimension (mpar_loc,mpvar), optional :: dfp
      integer, dimension (mpar_loc) :: ipar
      integer :: npar_loc, ipy_rec, ipz_rec
!
      integer, dimension (0:ncpus-1,0:ncpus-1) :: nmig
      integer, dimension (0:ncpus-1) :: k_move, k0_move, k0_move1
      integer :: ierr=0, status=0
      integer :: i, j, k, iproc_rec
      logical, save :: lfirstcall=.true.
!
      intent (out) :: fp, npar_loc, ipar, dfp
!
!  Initialise start index of migrating particles..
!
      if (lfirstcall) then
        do i=0,ncpus-1
          k0_move(i)=mpar_loc-npar_mig*(i+1)
        enddo
        k0_move1=k0_move-1
        lfirstcall=.false.
      endif
      k_move=k0_move1
!
!  Find out which particles are not in the local processor's yz-interval.
!
      nmig=0
      do k=npar_loc,1,-1
        if (fp(k,iyp)< xyz0_loc(2) .or. fp(k,iyp)>=xyz1_loc(2) .or. &
            fp(k,izp)< xyz0_loc(3) .or. fp(k,izp)>=xyz1_loc(3)) then
!  Calculate processor index of receiving processor.            
          iproc_rec = int((fp(k,iyp)-xyz0(2))/Lxyz_loc(2)) + &
              nprocy*int((fp(k,izp)-xyz0(3))/Lxyz_loc(3))
          if (ip<=8) print '(a,i6,a,i3,a,i3)', &
              'redist_particles_procs: Particle ', ipar(k), &
              ' moves out of proc ', iproc, &
              ' and into proc ', iproc_rec
!  Copy migrating particle to the end of the fp array.
          nmig(iproc,iproc_rec)=nmig(iproc,iproc_rec)+1
          k_move(iproc_rec)=k_move(iproc_rec)+1
          if (nmig(iproc,iproc_rec)>npar_mig) then
            print '(a,i3,a,i3,a)', &
                'redist_particles_procs: too many particles migrating '// &
                'from proc ', iproc, ' to proc ', iproc_rec
            print*, 'redist_particles_procs: set npar_mig higher '// &
                'in cparam.local'
            call fatal_error('redist_particles_procs','')
          endif
          fp(k_move(iproc_rec),:)=fp(k,:)
          if (present(dfp)) dfp(k_move(iproc_rec),:)=dfp(k,:)
          ipar(k_move(iproc_rec))=ipar(k)
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
      if (ip<=7) print*, 'redist_particles_procs: iproc, nmigrate = ', &
          iproc, sum(nmig(iproc,:))
!
!  Share information about number of migrating particles.
!
      do i=0,ncpus-1
        if (iproc/=i) call MPI_RECV(nmig(i,iproc), 1, MPI_INTEGER, i, &
            111, MPI_COMM_WORLD, status, ierr)
        if (iproc==i) then
          do j=0,ncpus-1
            call MPI_SEND(nmig(iproc,j), 1, MPI_INTEGER, j, &
                111, MPI_COMM_WORLD, ierr)
          enddo
        endif
      enddo
!
!  Set to receive.
!      
      do i=0,ncpus-1
        if (iproc/=i .and. nmig(i,iproc)/=0) then
          call MPI_RECV(fp(npar_loc+1:npar_loc+nmig(i,iproc),:), &
              mpvar*nmig(i,iproc), MPI_REAL, i, &
              222, MPI_COMM_WORLD, status,ierr)
          call MPI_RECV(ipar(npar_loc+1:npar_loc+nmig(i,iproc)), &
              nmig(i,iproc), MPI_INTEGER, i, &
              223, MPI_COMM_WORLD, status, ierr)
          if (present(dfp)) &
              call MPI_RECV(dfp(npar_loc+1:npar_loc+nmig(i,iproc),:), &
              mpvar*nmig(i,iproc), MPI_REAL, i, &
              333, MPI_COMM_WORLD, status, ierr)
          npar_loc=npar_loc+nmig(i,iproc)
        endif
!
!  Directed send.
!        
        if (iproc==i) then
          do j=0,ncpus-1
            if (nmig(iproc,j)/=0) then
              call MPI_SEND(fp(k0_move(j):k0_move(j)+nmig(iproc,j)-1,:), &
                  mpvar*nmig(iproc,j),  MPI_REAL, j, &
                  222, MPI_COMM_WORLD, ierr)
              call MPI_SEND(ipar(k0_move(j):k0_move(j)+nmig(iproc,j)-1), &
                  nmig(iproc,j),  MPI_INTEGER, j, &
                  223, MPI_COMM_WORLD, ierr)
              if (present(dfp)) &
                call MPI_SEND(dfp(k0_move(j):k0_move(j)+nmig(iproc,j)-1,:), &
                    mpvar*nmig(iproc,j), MPI_REAL, j, &
                    333, MPI_COMM_WORLD, ierr)
            endif
          enddo
        endif
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
!  successively calculate sum of a, which is supplied at each call.
!  Works for particle diagnostics.
!
!  02-jan-05/anders: adapted from sum_mn_name
!
      use Cdata
!
      real, dimension (:) :: a
      integer :: iname
      logical, optional :: lsqrt
!
      if (iname/=0) then
!
        fname(iname)=0.
        fname(iname)=fname(iname)+sum(a)
!
!  set corresponding entry in itype_name
!
        if (present(lsqrt)) then
          itype_name(iname)=ilabel_sum_sqrt_par
        else
          itype_name(iname)=ilabel_sum_par
        endif
!
      endif
!
    endsubroutine sum_par_name
!***********************************************************************
    subroutine find_lowest_cornerpoint(xxp,ix0,iy0,iz0)
!
!  Find index (ix0, iy0, iz0) of lowest corner point in the grid box
!  surrounding the point xxp.
!
!  23-jan-05/anders: coded
!
      use Cdata
      use Global
      use Mpicomm, only: stop_it
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
      ix0=4; iy0=4; iz0=4
!
      if (lfirstcall) then
        if (.not. all(lequidist)) then
          print*, 'find_lowest_cornerpoint: only works for equidistant grid!'
          call stop_it('find_lowest_cornerpoint')
        endif
        dx1=1/dx; dy1=1/dy; dz1=1/dz
        lfirstcall=.false.
      endif
!      
      if (nxgrid/=1) ix0 = nint((xxp(1)-x(1))*dx1) + 1
      if (nygrid/=1) iy0 = nint((xxp(2)-y(1))*dy1) + 1
      if (nzgrid/=1) iz0 = nint((xxp(3)-z(1))*dz1) + 1
!
!      call set_global(1.0,iy0,iz0,'nd',ix0)
!
      if ( (x(ix0)>xxp(1)) .and. nxgrid/=1) ix0=ix0-1
      if ( (y(iy0)>xxp(2)) .and. nygrid/=1) iy0=iy0-1
      if ( (z(iz0)>xxp(3)) .and. nzgrid/=1) iz0=iz0-1
!
    endsubroutine find_lowest_cornerpoint
!***********************************************************************
    subroutine interpolate_3d_1st(f,ii0,xxp,gp)
!
!  Interpolate the value of g to arbitrary (xp, yp, zp) coordinate
!  using first order formula 
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (3) :: xxp, gp
      integer :: ii0
!
      real, dimension (3) :: g1, g2, g3, g4, g5, g6, g7, g8
      real :: xp0, yp0, zp0
      real, save :: dxdydz1, dxdy1, dxdz1, dydz1, dx1, dy1, dz1
      integer, dimension (3) :: ixx0
      integer :: ii1
      logical :: lfirstcall=.true.
!
      intent(in)  :: f, xxp, ii0
      intent(out) :: gp
!
!  Determine index value of lowest lying corner point of grid box surrounding
!  the interpolation point.
!
      call find_lowest_cornerpoint(xxp,ixx0(1),ixx0(2),ixx0(3))
!
!  Check if the grid point interval is really correct.
!
      if ((x(ixx0(1))<=xxp(1) .and. x(ixx0(1)+1)>=xxp(1) .or. nxgrid==1) .and. &
          (y(ixx0(2))<=xxp(2) .and. y(ixx0(2)+1)>=xxp(2) .or. nygrid==1) .and. &
          (z(ixx0(3))<=xxp(3) .and. z(ixx0(3)+1)>=xxp(3) .or. nzgrid==1)) then
        ! Everything okay
      else
        print*, 'interpolate_3d_1st: Interpolation point does not ' // &
            'lie within the calculated grid point interval.'
        print*, 'iproc = ', iproc
        print*, 'mx, x(1), x(mx) = ', mx, x(1), x(mx)
        print*, 'my, y(1), y(my) = ', my, y(1), y(my)
        print*, 'mz, z(1), z(mz) = ', mz, z(1), z(mz)
        print*, 'ixx0(1), ixx0(2), ixx0(3) = ', ixx0(1), ixx0(2), ixx0(3)
        print*, 'xp, xp0, xp1 = ', xxp(1), x(ixx0(1)), x(ixx0(1)+1)
        print*, 'yp, yp0, yp1 = ', xxp(2), y(ixx0(2)), y(ixx0(2)+1)
        print*, 'zp, zp0, zp1 = ', xxp(3), z(ixx0(3)), z(ixx0(3)+1)
        call stop_it('interpolate_3d_1st')
      endif
!
!  Redefine the interpolation point in coordinates relative to lowest corner.
!
      xp0=xxp(1)-x(ixx0(1))
      yp0=xxp(2)-y(ixx0(2))
      zp0=xxp(3)-z(ixx0(3))
!
!  Calculate derived grid spacing variables in the first call to this sub.
!      
      if (lfirstcall) then
        dxdydz1=1/(dx*dy*dz)
        dxdy1=1/(dx*dy)
        dxdz1=1/(dx*dz)
        dydz1=1/(dy*dz)
        dx1=1/dx
        dy1=1/dy
        dz1=1/dz
        lfirstcall=.false.
      endif
!
!  Function values at all corners.
!
      ii1=ii0+2
      g1=f(ixx0(1)+1,ixx0(2)+1,ixx0(3)+1,ii0:ii1)
      g2=f(ixx0(1)  ,ixx0(2)+1,ixx0(3)+1,ii0:ii1)
      g3=f(ixx0(1)+1,ixx0(2)  ,ixx0(3)+1,ii0:ii1)
      g4=f(ixx0(1)+1,ixx0(2)+1,ixx0(3)  ,ii0:ii1)
      g5=f(ixx0(1)  ,ixx0(2)  ,ixx0(3)+1,ii0:ii1)
      g6=f(ixx0(1)  ,ixx0(2)+1,ixx0(3)  ,ii0:ii1)
      g7=f(ixx0(1)+1,ixx0(2)  ,ixx0(3)  ,ii0:ii1)
      g8=f(ixx0(1)  ,ixx0(2)  ,ixx0(3)  ,ii0:ii1)
!
!  Interpolation formula.
!
      gp = xp0*yp0*zp0*dxdydz1*(g1-g2-g3-g4+g5+g6+g7-g8) &
          + xp0*yp0*dxdy1*(g4-g6-g7+g8) &
          + xp0*zp0*dxdz1*(g3-g5-g7+g8) &
          + yp0*zp0*dydz1*(g2-g5-g6+g8) &
          + xp0*dx1*(g7-g8) + yp0*dy1*(g6-g8) + zp0*dz1*(g5-g8) + g8
!
    endsubroutine interpolate_3d_1st
!***********************************************************************
    subroutine map_xx_vv_grid(xxp,vvp)
!
!  Find index (ix0, iy0, iz0) of lowest corner point in the grid box
!  surrounding the point xxp.
!
!  23-jan-05/anders: coded
!
      use Cdata
      use Global
      use Mpicomm, only: stop_it
!
      real, dimension(3) :: xxp, vvp
      integer :: ix0, iy0, iz0
!
      real, save :: dx1, dy1, dz1
      logical, save :: lfirstcall=.true.
!
      intent(in)  :: xxp, vvp
!
      ix0=4; iy0=4; iz0=4
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
      if (nxgrid/=1) ix0 = nint((xxp(1)-x(1))*dx1) + 1
      if (nygrid/=1) iy0 = nint((xxp(2)-y(1))*dy1) + 1
      if (nzgrid/=1) iz0 = nint((xxp(3)-z(1))*dz1) + 1
!
      call set_global_point(1.0,ix0,iy0,iz0,'np')
      call set_global_point(vvp,ix0,iy0,iz0,'uupsum')
!
    endsubroutine map_xx_vv_grid
!***********************************************************************

endmodule Particles_sub
