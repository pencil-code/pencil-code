! $Id$
!
!  DOCUMENT ME OR MOVE ME TO EXPERIMENTAL
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
module Fixed_point
!
  use Cdata
  use Cparam
  use Mpicomm
  use Messages
  use Streamlines
!
  implicit none
!
! a few constants
  integer :: MERGE_FIXED = 97
  integer :: FINISHED_FIXED = 98
! the arrays with the values for x, y, and z for all cores (tracer xyz)
  real, pointer, dimension (:) :: xt, yt, zt
! fixed points array for this core
  real, dimension (1000,3) :: fixed_points
! fixed points array for all fixed points  
!   real, pointer, dimension (:,:) :: fixed_points_all
  real, dimension (1000,3) :: fixed_points_all
! temporary array, which stores the transposed fixed_points array
! this is needed for the MPI communication and to keep the order of the indices
  real, dimension (3,1000) :: buffer_tmp
! total number of fixed points of one core
  integer :: fidx, fidx_other
! global number of fixed points
  integer :: fidx_all
  real :: fidx_read
!
  real, public :: tfixed_points
  integer, public :: nfixed_points
! Time value to be written together with the fixed points.
  real :: tfixed_points_write
!
  contains
!
!***********************************************************************
  subroutine fixed_points_prepare
!
!  Prepare lfixed_points for writing the fixed points into the fixed points file.
!
!  14-mar-12/simon: coded
!
    use Sub, only: read_snaptime, update_snaptime
!
    integer, save :: ifirst=0
!
    character (len=fnlen) :: file
!   fixed points file name
    character(len=1024) :: filename, str_tmp
    real, pointer, dimension (:,:) :: fixed_tmp
!   integer which checks if end of file has been reached
    integer :: IOstatus, j
!   increment variable for the processor index
    integer :: proc_idx, ierr, x_proc, y_proc
!
!  Output fixed points-data in 'tfixed_points' time intervals
!
    file = trim(datadir)//'/tfixed_points.dat'
!
    if (ifirst==0) then
      do proc_idx=0,(nprocx*nprocy*nprocz-1)
!       artificially serialize this routine to avoid accessing the same file by different cores
        call MPI_BARRIER(MPI_comm_world, ierr)
        if (proc_idx == iproc) then
          call read_snaptime(file,tfixed_points,ntracers,dfixed_points,t)
!         Read the previous fixed points from the file.
!           open(unit = 1, file = adjustl(trim('data/fixed_points.dat')), form = "unformatted")
          open(unit = 1, file = 'data/fixed_points.dat', form = "unformatted")
!         loop until we find the last entry
          IOstatus = 0
!
          read(1,iostat = IOstatus) tfixed_points_write
          if (IOstatus == 0) then
            read(1) fidx_read
            fidx_all = int(fidx_read)
!             write(*,*) "reading ", fidx_all, " many fixed points"
!             allocate(fixed_points_all(fidx,3))
            do j=1,fidx_all
              read(1) fixed_points_all(j,:)
            enddo
            close(1)
          endif
          ifirst=1
        endif
      enddo
    endif         
!
!   find out which fixed points are destinated for this core
!
    fidx = 0
    do j = 1,fidx_all
!     find the corresponding core
!     NB: the whole thing only works for nprocz = 1
      x_proc = floor((fixed_points_all(j,1)-x0)*real(nprocx)/Lx)
      y_proc = floor((fixed_points_all(j,2)-y0)*real(nprocy)/Ly)
      if ((x_proc == ipx) .and. (y_proc == ipy)) then
        fidx = fidx + 1
        fixed_points(fidx,:) = fixed_points_all(j,:)
      endif      
    enddo
!     deallocate(fixed_points_all)
!
!     ==== old code ===
!
!     if (ifirst==0) then
!       call read_snaptime(file,tfixed_points,ntracers,dtracers,t)
! !     Read the previous fixed points from the file.
!       write(str_tmp, "(I10.1,A)") iproc, '/fixed_points.dat'
!       write(filename, *) 'data/proc', adjustl(trim(str_tmp))
!       open(unit = 1, file = adjustl(trim(filename)), form = "unformatted")
! !     loop until we find the last entry
!       IOstatus = 0
! !
!       read(1,iostat = IOstatus) tfixed_points_write
!       if (IOstatus == 0) then
!         read(1) fidx_read
!         fidx = int(fidx_read)
!         allocate(fixed_tmp(fidx,3))
!         do j=1,fidx
!           read(1) fixed_tmp(j,:)
!         enddo
!         fixed_points(1:fidx,:) = fixed_tmp(:,:)
!         deallocate(fixed_tmp)
!         close(1)
!       endif
!       ifirst=1
!     endif
!
!     ==== old code end ===
!
!  This routine sets lfixed_points=T whenever its time to write the fixed points
!
    call update_snaptime(file,tfixed_points,nfixed_points,dfixed_points,t,lfixed_points)
!
!  Save current time so that the time that is written out is not
!  from the next time step
!
    if (lfixed_points) tfixed_points_write = t
!
  endsubroutine fixed_points_prepare
!***********************************************************************
  subroutine get_fixed_point(f, point, fixed_point, q, vv)
!
!   Finds the fixed point near 'point'.
!   Returns the position of the fixed point together with the integrated value q.
!
!   01-mar-12/simon: coded
!   09-aug-12/anthony: modified to use Newton's method
!
    real, dimension (mx,my,mz,mfarray) :: f
    real :: point(2), fixed_point(2)
!   the integrated quantity along the field line
    real :: q
    real, pointer, dimension (:,:,:,:) :: vv
!   the points for the extrapolation
    real, pointer, dimension (:,:) :: tracers
    integer :: iter
    real :: dl
    real :: fjac(2,2), fjin(2,2)
    real :: det, ff(2), dpoint(2)
!
    intent(out) :: fixed_point, q
!
    allocate(tracers(5,7))
!
!   step-size for calculating Jacobian by finite differences:
    dl = min(dx,dy)/100.
!
    iter = 0
    do
!
!   trace field lines at original point and for Jacobian:
!   (second order seems to be enough)
!
      tracers(1,:) = (/point(1),point(2),point(1),point(2),z(1+nghost)-ipz*nz*dz+dz,0.,1./)
      tracers(2,:) = (/point(1)-1*dl,point(2),point(1)-1*dl,point(2),z(1+nghost)-ipz*nz*dz+dz,0.,2./)
      tracers(3,:) = (/point(1)+1*dl,point(2),point(1)+1*dl,point(2),z(1+nghost)-ipz*nz*dz+dz,0.,2./)
      tracers(4,:) = (/point(1),point(2)-1*dl,point(1),point(2)-1*dl,z(1+nghost)-ipz*nz*dz+dz,0.,4./)
      tracers(5,:) = (/point(1),point(2)+1*dl,point(1),point(2)+1*dl,z(1+nghost)-ipz*nz*dz+dz,0.,4./)
      call trace_streamlines(f,tracers,5,vv)
!
!   compute Jacobian at x:
!
      fjac(1,1) = (tracers(3,3) - tracers(3,1) - &
          (tracers(2,3) - tracers(2,1)))/2.0/dl
      fjac(1,2) = (tracers(5,3) - tracers(5,1) - &
          (tracers(4,3) - tracers(4,1)))/2.0/dl
      fjac(2,1) = (tracers(3,4) - tracers(3,2) - &
          (tracers(2,4) - tracers(2,2)))/2.0/dl
      fjac(2,2) = (tracers(5,4) - tracers(5,2) - &
          (tracers(4,4) - tracers(4,2)))/2.0/dl
!
!   check function convergence:
!
      ff(1) = tracers(1,3) - tracers(1,1)
      ff(2) = tracers(1,4) - tracers(1,2)
      if (sum(abs(ff)) <= 1.0e-4) then
        q = tracers(1,7)
        fixed_point = point
        exit
      end if
!
!   invert Jacobian and do Newton step:
!
      det = fjac(1,1)*fjac(2,2) - fjac(1,2)*fjac(2,1)
      fjin(1,1) = fjac(2,2)
      fjin(2,2) = fjac(1,1)
      fjin(1,2) = -fjac(1,2)
      fjin(2,1) = -fjac(2,1)
      fjin = fjin/det
      dpoint(1) = -fjin(1,1)*ff(1) - fjin(1,2)*ff(2)
      dpoint(2) = -fjin(2,1)*ff(1) - fjin(2,2)*ff(2)
      point = point + dpoint
!
!   check root convergence:
!
      if (sum(abs(dpoint)) <= 1.0e-4) then
        q = tracers(1,7)
        fixed_point = point
        exit
      end if
!
      if (iter > 20) then
        q = tracers(1,7)
        fixed_point = point
!         fixed_point = (/x0-0.1, y0-0.1/)
        write(*,*) 'Warning: Newton not converged'
        exit
      endif
      iter = iter + 1
    enddo
!
!     write(*,*) fixed_point
    deallocate(tracers)
  end subroutine get_fixed_point
!***********************************************************************
  recursive function edge(f, sx, sy, diff1, diff2, phi_min, vv, rec) result(dtot)
!
! Computes rotation along one edge (recursively until phi_min is reached).
!
! 01-mar-12/simon: coded
!
    real, dimension (mx,my,mz,mfarray) :: f
!   the two corners
    real :: sx(2), sy(2)
!   F1(x0,y0) - (x0,y0) at the corners
    real diff1(2), diff2(2)
!   tolerated change in angle
    real :: phi_min
!   the tracer field
    real, pointer, dimension (:,:,:,:) :: vv
!   number of recursions
    integer :: rec
!   total rotation along this edge
    real :: dtot
!   tracer for any eventual field line on the edge
    real, pointer, dimension (:,:) :: tracer
!   intermediate point on the edge
    real :: xm, ym, diffm(2)
!
    allocate(tracer(1,7))
!
    dtot = atan2(diff1(1)*diff2(2) - diff2(1)*diff1(2), diff1(1)*diff2(1) + diff1(2)*diff2(2))
!
    if ((abs(dtot) > phi_min) .and. (rec < 5)) then
      xm = 0.5*(sx(1)+sx(2))
      ym = 0.5*(sy(1)+sy(2))
!     trace intermediate field line
      tracer(1,:) = (/xm,ym,xm,ym,z(1+nghost)-ipz*nz*dz+dz,0.,0./)
      call trace_streamlines(f,tracer,1,vv)
      if ((tracer(1,6) >= l_max) .or. (tracer(1,5) < zt(nzgrid)-dz)) then
!       discard any streamline which does not converge or hits the boundary
        dtot = 0.
      else
        diffm = (/tracer(1,3)-tracer(1,1), tracer(1,4)-tracer(1,2)/)
        if (diffm(1)**2 + diffm(2)**2 /= 0) &
            diffm = diffm / sqrt(diffm(1)**2 + diffm(2)**2)
        dtot = edge(f,(/sx(1),xm/), (/sy(1),ym/), diff1, diffm, phi_min, vv, rec+1) + &
            edge(f,(/xm,sx(2)/), (/ym,sy(2)/), diffm, diff2, phi_min, vv, rec+1)
      endif
    else
      dtot = dtot
    endif
!
    deallocate(tracer)
!
  end function edge
!***********************************************************************
  subroutine pindex(f, sx, sy, diff, phi_min, vv, poincare)
!
! Finds the Poincare index of this grid cell.
!
! 01-mar-12/simon: coded
!
    real, dimension (mx,my,mz,mfarray) :: f
!   corners of the grid square
    real :: sx(2), sy(2)
!   F1(x0,y0) - (x0,y0) at the corners
    real diff(4,2)
!   Poincare index for this grid cell
    real :: poincare
!   tolerated change in angle
    real :: phi_min
    real, pointer, dimension (:,:,:,:) :: vv
!
    poincare = 0
    poincare = poincare + edge(f, (/sx(1),sx(2)/), (/sy(1),sy(1)/), diff(1,:), diff(2,:), phi_min, vv, 1)
    poincare = poincare + edge(f, (/sx(2),sx(2)/), (/sy(1),sy(2)/), diff(2,:), diff(3,:), phi_min, vv, 1)
    poincare = poincare + edge(f, (/sx(2),sx(1)/), (/sy(2),sy(2)/), diff(3,:), diff(4,:), phi_min, vv, 1)
    poincare = poincare + edge(f, (/sx(1),sx(1)/), (/sy(2),sy(1)/), diff(4,:), diff(1,:), phi_min, vv, 1)
!
  end subroutine pindex
!***********************************************************************
  subroutine get_fixed_points(f, tracers, vv)
!
!  trace stream lines of the vector field stored in f(:,:,:,iaa)
!
!   13-feb-12/simon: coded
!   09-aug-12/anthony: modified to subsample identified cells
!
    use Sub
!
    real, dimension (mx,my,mz,mfarray) :: f
    real, pointer, dimension (:,:) :: tracers, tracers2, tracer_tmp
    real, pointer, dimension (:,:,:,:) :: vv
!   filename for the fixed point output
    real :: poincare, diff(4,2), phi_min
    integer :: j, l, addx, addy, proc_idx, ierr, flag
    integer, dimension (MPI_STATUS_SIZE) :: status
!   array with all finished cores
    integer :: finished_rooting(nprocx*nprocy*nprocz)
!   variables for the final non-blocking mpi communication
    integer :: request_finished_send(nprocx*nprocy*nprocz)
    integer :: request_finished_rcv(nprocx*nprocy*nprocz)
    real, pointer, dimension(:,:) :: tracers3
    real :: x3min, x3max, y3min, y3max, minx, miny, min3, x3, y3, diff3
    integer :: i1, j1, k1, nt3
!
    addx = 0; addy = 0
    if (ipx /= nprocx-1) addx = 1
    if (ipy /= nprocy-1) addy = 1
!
    allocate(xt(nxgrid*trace_sub))
    allocate(yt(nygrid*trace_sub))
    allocate(zt(nz))
!
    phi_min = pi/8.
!     phi_min = pi/8.-2.**(-15)
!
!   compute the array with the global xyz values
    do j=1,(nxgrid*trace_sub)
      xt(j) = x(1+nghost) - ipx*nx*dx + (j-1)*dx/trace_sub
    enddo
    do j=1,(nygrid*trace_sub)
      yt(j) = y(1+nghost) - ipy*ny*dy + (j-1)*dy/trace_sub
    enddo
    do j=1,nzgrid
      zt(j) = z(1+nghost) - ipz*nz*dz + (j-1)*dz
    enddo
!
!   Make sure the core boundaries are considered.
    allocate(tracer_tmp(1,7))
    if (addx == 1 .or. addy == 1) then
      allocate(tracers2(nx*ny*trace_sub**2+trace_sub*(ny*addx+nx*addy)+addx*addy,7))
!     Assign the values without the core boundaries.
      do l=1,ny*trace_sub
        do j=1,nx*trace_sub
          tracers2(j+(l-1)*(nx*trace_sub+addx),:) = tracers(j+(l-1)*nx*trace_sub,:)
        enddo
      enddo
!     Recompute values at the x-boundaries of the core.
      if (addx == 1) then
        do l=1,ny*trace_sub
          tracer_tmp(1,:) = (/x(nghost+nx+1),yt(l)+ipy*ny*dy,x(nghost+nx+1),yt(l)+ipy*ny*dy,0.,0.,0./)
          call trace_streamlines(f,tracer_tmp,1,vv)
          tracers2(l*(nx*trace_sub+addx),:) = tracer_tmp(1,:)
        enddo
      endif
!     Recompute values at the y-boundaries of the core.
      if (addy == 1) then
        do j=1,nx*trace_sub
          tracer_tmp(1,:) = (/xt(j)+ipx*nx*dx,y(nghost+ny+1),xt(j)+ipx*nx*dx,y(nghost+ny+1),0.,0.,0./)
          call trace_streamlines(f,tracer_tmp,1,vv)
          tracers2(j+(nx*trace_sub+addx)*ny*trace_sub,:) = tracer_tmp(1,:)
        enddo
      endif
!     Compute the value at the corner of the core.
      if (addx*addy == 1) then
        tracer_tmp(1,:) = &
            (/x(nghost+nx+1),y(nghost+ny+1),x(nghost+nx+1),y(nghost+ny+1),0.,0.,0./)
        call trace_streamlines(f,tracer_tmp,1,vv)
        tracers2(nx*ny*trace_sub**2+trace_sub*(ny+nx)+1,:) = tracer_tmp(1,:)
      endif
    else
      allocate(tracers2(nx*ny*trace_sub**2,7))
      tracers2 = tracers
    endif
!
!   Tell every other core that we have finished.
    finished_rooting(:) = 0
    finished_rooting(iproc+1) = 1
    do proc_idx=0,(nprocx*nprocy*nprocz-1)
      if (proc_idx /= iproc) then
        call MPI_ISEND(finished_rooting(iproc+1), 1, MPI_integer, proc_idx, FINISHED_FIXED, &
            MPI_comm_world, request_finished_send(proc_idx+1), ierr)
        if (ierr /= MPI_SUCCESS) &
            call fatal_error("streamlines", "MPI_ISEND could not send")
        call MPI_IRECV(finished_rooting(proc_idx+1), 1, MPI_integer, proc_idx, FINISHED_FIXED, &
            MPI_comm_world, request_finished_rcv(proc_idx+1), ierr)
        if (ierr /= MPI_SUCCESS) &
            call fatal_error("streamlines", "MPI_IRECV could not create a receive request")
      endif
    enddo
!
!   Trace a dummy field line to start a field receive request.
    tracer_tmp(1,:) = (/(x(1+nghost)+x(nx+nghost))/2.,(y(1+nghost)+y(ny+nghost))/2., &
        (x(1+nghost)+x(nx+nghost))/2.,(y(1+nghost)+y(ny+nghost))/2.,0.,0.,0./)
    call trace_streamlines(f,tracer_tmp,1,vv)
!
!   Wait for all cores to compute their missing stream lines.
    do
!     Check if a core has finished and update finished_rooting array.
      do proc_idx=0,(nprocx*nprocy*nprocz-1)
        if ((proc_idx /= iproc) .and. (finished_rooting(proc_idx+1) == 0)) then
          flag = 0
          call MPI_TEST(request_finished_rcv(proc_idx+1),flag,status,ierr)
          if (ierr /= MPI_SUCCESS) &
              call fatal_error("fixed_points", "MPI_TEST failed")
          if (flag == 1) then
            finished_rooting(proc_idx+1) = 1
          endif
        endif
      enddo
      if (sum(finished_rooting) == nprocx*nprocy*nprocz) exit
!     Trace a dummy field line to start a field receive request.
      tracer_tmp(1,:) = (/(x(1+nghost)+x(nx+nghost))/2.,(y(1+nghost)+y(ny+nghost))/2., &
          (x(1+nghost)+x(nx+nghost))/2.,(y(1+nghost)+y(ny+nghost))/2.,0.,0.,0./)
      call trace_streamlines(f,tracer_tmp,1,vv)
    enddo
!
    call MPI_BARRIER(MPI_comm_world, ierr)
!
!   Find possible fixed points in each grid cell.
!
!   index of the fixed point
    fidx = 1
    do j=1,(nx*trace_sub+addx-1)
      do l=1,(ny*trace_sub+addy-1)
        diff(1,:) = (/(tracers2(j+(l-1)*(nx*trace_sub+addx),3)-tracers2(j+(l-1)*(nx*trace_sub+addx),1)) , &
            (tracers2(j+(l-1)*(nx*trace_sub+addx),4)-tracers2(j+(l-1)*(nx*trace_sub+addx),2))/)            
        if (diff(1,1)**2+diff(1,2)**2 /= 0) &
            diff(1,:) = diff(1,:) / sqrt(diff(1,1)**2+diff(1,2)**2)            
        diff(2,:) = (/(tracers2(j+1+(l-1)*(nx*trace_sub+addx),3)-tracers2(j+1+(l-1)*(nx*trace_sub+addx),1)) , &
            (tracers2(j+1+(l-1)*(nx*trace_sub+addx),4)-tracers2(j+1+(l-1)*(nx*trace_sub+addx),2))/)
        if (diff(2,1)**2+diff(2,2)**2 /= 0) &
            diff(2,:) = diff(2,:) / sqrt(diff(2,1)**2+diff(2,2)**2)
        diff(3,:) = (/(tracers2(j+1+l*(nx*trace_sub+addx),3)-tracers2(j+1+l*(nx*trace_sub+addx),1)) , &
            (tracers2(j+1+l*(nx*trace_sub+addx),4)-tracers2(j+1+l*(nx*trace_sub+addx),2))/)
        if (diff(3,1)**2+diff(3,2)**2 /= 0) &
            diff(3,:) = diff(3,:) / sqrt(diff(3,1)**2+diff(3,2)**2)
        diff(4,:) = (/(tracers2(j+l*(nx*trace_sub+addx),3)-tracers2(j+l*(nx*trace_sub+addx),1)) , &
            (tracers2(j+l*(nx*trace_sub+addx),4)-tracers2(j+l*(nx*trace_sub+addx),2))/)
        if (diff(4,1)**2+diff(4,2)**2 /= 0) &
            diff(4,:) = diff(4,:) / sqrt(diff(4,1)**2+diff(4,2)**2)
!       Get the Poincare index for this grid cell
        call pindex(f, xt(j:j+1)+ipx*nx*dx, yt(l:l+1)+ipy*ny*dy, diff, phi_min, vv, poincare)
!       find the fixed point in this cell
        if (abs(poincare) >= 5.) then
!
!       subsample to get starting point for iteration:
!
          nt3 = 4
          x3min = tracers2(j+(l-1)*(nx*trace_sub+addx),1)
          y3min = tracers2(j+(l-1)*(nx*trace_sub+addx),2)
          x3max = tracers2(j+1+(l-1)*(nx*trace_sub+addx),1)
          y3max = tracers2(j+l*(nx*trace_sub+addx),2)
          allocate(tracers3(nt3*nt3,7))
          i1=1
          do j1=1,nt3
             do k1=1,nt3
               x3 = x3min + (j1-1.)/(real(nt3)-1)*(x3max - x3min)
               y3 = y3min + (k1-1.)/(real(nt3)-1)*(y3max - y3min)
               tracers3(i1,:) = (/x3, y3, x3, y3, z(1+nghost)-ipz*nz*dz+dz, &
                    0., 4./)
               i1 = i1 + 1
             end do
          end do
          call trace_streamlines(f,tracers3,nt3*nt3,vv)
          min3 = 1.0e6
          minx = x3min
          miny = y3min
          i1=1
          do j1=1,nt3
             do k1=1,nt3
               diff3 = (tracers3(i1,3) - tracers3(i1,1))**2 + &
                    (tracers3(i1,4) - tracers3(i1,2))**2
               if (diff3 < min3) then
                 min3 = diff3
                 minx = x3min + (j1-1.)/(real(nt3)-1)*(x3max - x3min)
                 miny = y3min + (k1-1.)/(real(nt3)-1)*(y3max - y3min)
               end if
               i1 = i1 + 1
             end do
          end do
          deallocate(tracers3)
!
!       iterate to find the fixed point from this starting point:        
!
          call get_fixed_point(f,(/minx, miny/), &
              fixed_points(fidx,1:2), fixed_points(fidx,3), vv)
!
!       check that fixed point lies inside the cell:
!
          if ((fixed_points(fidx,1) < tracers2(j+(l-1)*(nx*trace_sub+addx),1)) .or. &
              (fixed_points(fidx,1) > tracers2(j+1+(l-1)*(nx*trace_sub+addx),1)) .or. &
              (fixed_points(fidx,2) < tracers2(j+(l-1)*(nx*trace_sub+addx),2)) .or. &
              (fixed_points(fidx,2) > tracers2(j+l*(nx*trace_sub+addx),2))) then
            write(*,*) iproc, "warning: fixed point lies outside the cell"
          else
            fidx = fidx+1
          endif
        endif
      enddo
    enddo
    fidx = fidx - 1
!
!   Tell every other core that we have finished.
    finished_rooting(:) = 0
    finished_rooting(iproc+1) = 1
    do proc_idx=0,(nprocx*nprocy*nprocz-1)
      if (proc_idx /= iproc) then
        call MPI_ISEND(finished_rooting(iproc+1), 1, MPI_integer, proc_idx, FINISHED_FIXED, &
            MPI_comm_world, request_finished_send(proc_idx+1), ierr)
        if (ierr /= MPI_SUCCESS) &
            call fatal_error("streamlines", "MPI_ISEND could not send")
        call MPI_IRECV(finished_rooting(proc_idx+1), 1, MPI_integer, proc_idx, FINISHED_FIXED, &
            MPI_comm_world, request_finished_rcv(proc_idx+1), ierr)
        if (ierr /= MPI_SUCCESS) &
            call fatal_error("streamlines", "MPI_IRECV could not create a receive request")
      endif
    enddo
!
!   Wait for all cores to compute their missing stream lines.
    do
!     Trace a dummy field line to start a field receive request.
      tracer_tmp(1,:) = (/(x(1+nghost)+x(nx+nghost))/2.,(y(1+nghost)+y(ny+nghost))/2., &
          (x(1+nghost)+x(nx+nghost))/2.,(y(1+nghost)+y(ny+nghost))/2.,0.,0.,0./)
      call trace_streamlines(f,tracer_tmp,1,vv)
!     Check if a core has finished and update finished_rooting array.
      do proc_idx=0,(nprocx*nprocy*nprocz-1)
        if (proc_idx /= iproc) then
          flag = 0
          call MPI_TEST(request_finished_rcv(proc_idx+1),flag,status,ierr)
          if (ierr /= MPI_SUCCESS) &
              call fatal_error("fixed_points", "MPI_TEST failed")
          if (flag == 1) then
            finished_rooting(proc_idx+1) = 1
          endif
        endif
      enddo
      if (sum(finished_rooting) == nprocx*nprocy*nprocz) then
        exit
      endif
    enddo
!
    close(1)
!
    deallocate(xt)
    deallocate(yt)
    deallocate(zt)
    deallocate(tracers2)
    deallocate(tracer_tmp)
!
  endsubroutine get_fixed_points
!***********************************************************************
  subroutine wfixed_points(f,path)
!
!   Write the tracers values to tracer.dat.
!   This should be called during runtime.
!
!   14-mar-12/simon: coded
!
    use Sub
    use General, only: keep_compiler_quiet
!
    real, dimension (mx,my,mz,mfarray) :: f
    character(len=*) :: path
!   the traced field
    real, pointer, dimension (:,:,:,:) :: vv
!   filename for the tracer output
    character(len=1024) :: filename, str_tmp
    integer :: i, j, flag
    integer, dimension (MPI_STATUS_SIZE) :: status
    real :: point(2)
!   array with all finished cores
    integer :: finished_rooting(nprocx*nprocy*nprocz)
!   variables for the final non-blocking mpi communication
    integer :: request_finished_send(nprocx*nprocy*nprocz)
    integer :: request_finished_rcv(nprocx*nprocy*nprocz)
    real, pointer, dimension (:,:) :: tracer_tmp
    integer :: ierr, proc_idx
    character (len=labellen) :: trace_field='bb'
!   array with indices of fixed points to discard (double and too close ones)
    integer :: discard(1000)
!   increment variable for the processor index
    integer :: x_proc, y_proc
!
!   allocate memory for the traced field
    allocate(vv(nx,ny,nz,3))
!
    allocate(tracer_tmp(1,7))
!
    call keep_compiler_quiet(path)
!
!   TODO: include other fields as well
    if (trace_field == 'bb' .and. ipz == 0) then
!     convert the magnetic vector potential into the magnetic field
      do m=m1,m2
        do n=n1,n2
          call curl(f,iaa,vv(:,m-nghost,n-nghost,:))
        enddo
      enddo
    endif
!
    do j=1,fidx
      point = fixed_points(j,1:2)
      call get_fixed_point(f, point, fixed_points(j,1:2), fixed_points(j,3), vv)
    enddo
!
!   Make sure there is a delay before sending the finished signal, so that we
!   don't end up with the last process to finish stuck in trace_streamlines.
!
    if (nprocx*nprocy*nprocz > 1) &
        call sleep(1)
!
!   Tell every other core that we have finished.
    finished_rooting(:) = 0
    finished_rooting(iproc+1) = 1
    do proc_idx=0,(nprocx*nprocy*nprocz-1)
      if (proc_idx /= iproc) then
        call MPI_ISEND(finished_rooting(iproc+1), 1, MPI_integer, proc_idx, FINISHED_FIXED, &
            MPI_comm_world, request_finished_send(proc_idx+1), ierr)
        if (ierr /= MPI_SUCCESS) &
            call fatal_error("streamlines", "MPI_ISEND could not send")
        call MPI_IRECV(finished_rooting(proc_idx+1), 1, MPI_integer, proc_idx, FINISHED_FIXED, &
            MPI_comm_world, request_finished_rcv(proc_idx+1), ierr)
        if (ierr /= MPI_SUCCESS) &
            call fatal_error("streamlines", "MPI_IRECV could not create a receive request")
      endif
    enddo
!
!   Wait for all cores to compute their missing stream lines.
    do
!     Trace a dummy field line to start a field receive request.
      tracer_tmp(1,:) = (/(x(1+nghost)+x(nx+nghost))/2.,(y(1+nghost)+y(ny+nghost))/2., &
          (x(1+nghost)+x(nx+nghost))/2.,(y(1+nghost)+y(ny+nghost))/2.,0.,0.,0./)
      call trace_streamlines(f,tracer_tmp,1,vv)
!     Check if a core has finished and update finished_rooting array.
      do proc_idx=0,(nprocx*nprocy*nprocz-1)
        if (proc_idx /= iproc) then
          flag = 0
          call MPI_TEST(request_finished_rcv(proc_idx+1),flag,status,ierr)
          if (ierr /= MPI_SUCCESS) &
              call fatal_error("fixed_points", "MPI_TEST failed")
          if (flag == 1) then
            finished_rooting(proc_idx+1) = 1
          endif
        endif
      enddo
      if (sum(finished_rooting) == nprocx*nprocy*nprocz) then
        exit
      endif
    enddo
!
!   Wait for other cores. This ensures that uncomplete fixed points wont get written out.
    call MPI_BARRIER(MPI_comm_world, ierr)
!
!   communicate the fixed points to proc0
    if (iproc == 0) then
!       allocate(fixed_points_all(1000,3))
      fixed_points_all(1:fidx,:) = fixed_points(1:fidx,:)
!     receive the fixed_points from the other cores
      fidx_all = fidx
      do proc_idx=1,(nprocx*nprocy*nprocz-1)
        fidx_other = 0
!       receive the number of fixed points of that proc
        call MPI_RECV(fidx_other, 1, MPI_integer, proc_idx, MERGE_FIXED, MPI_comm_world, status, ierr)
        if (ierr /= MPI_SUCCESS) &
            call fatal_error("streamlines", "MPI_RECV could not receive")
!       receive the fixed points form that proc
        if (fidx_other > 0) then
          call MPI_RECV(buffer_tmp, fidx_other*3, MPI_real, proc_idx, MERGE_FIXED, MPI_comm_world, status, ierr)
          fixed_points_all(fidx_all+1:fidx_all+fidx_other,:) = transpose(buffer_tmp(:,1:fidx_other))
          if (ierr /= MPI_SUCCESS) &
              call fatal_error("streamlines", "MPI_RECV could not receive")
!               write(*,*) "fidx = ", fidx, " fidx_other = ", fidx_other
          fidx_all = fidx_all + fidx_other
        endif
      enddo
!
!     Check whether fixed points are too close or out of the domain.
!
      discard(:) = 0
      do j=1,fidx_all
        if ((fixed_points_all(j,1) < x0) .or. (fixed_points_all(j,1) > (x0+Lx)) .or. &
           (fixed_points_all(j,2) < y0) .or. (fixed_points_all(j,2) > (y0+Ly))) then
          discard(j) = 1
        else
          do i=j+1,fidx_all
            if ((abs(fixed_points_all(i,1) - fixed_points_all(j,1)) < dx/2) .and. &
                (abs(fixed_points_all(i,2) - fixed_points_all(j,2)) < dy/2)) then
              discard(i) = 1
            endif
          enddo
        endif
      enddo
!
      open(unit = 1, file = 'data/fixed_points.dat', form = "unformatted", position = "append")
      write(1) tfixed_points_write
      write(1) float(fidx_all-sum(discard))
      do j=1,fidx_all
        if (discard(j) == 0) then
          write(1) fixed_points_all(j,:)
!           write(*,*) "writing ", fidx_all-sum(discard), " many fixed points"
!           write(*,*) "writing ", fixed_points_all(j,:)
        else
!           write(*,*) "discarding ", fixed_points_all(j,:)
        endif
      enddo
      close(1)
!
!       deallocate(fixed_points_all)
    else
!       write(*,*) "sending fidx = ", fidx
      call MPI_SEND(fidx, 1, MPI_integer, 0, MERGE_FIXED, MPI_comm_world, ierr)
      if (ierr /= MPI_SUCCESS) &
          call fatal_error("streamlines", "MPI_SEND could not send")
      if (fidx > 0) then
        buffer_tmp = transpose(fixed_points)
        call MPI_SEND(buffer_tmp, fidx*3, MPI_real, 0, MERGE_FIXED, MPI_comm_world, ierr)              
        if (ierr /= MPI_SUCCESS) &
            call fatal_error("streamlines", "MPI_SEND could not send")
      endif
    endif
!
!   redistribute the fixed points to the appropriate cores
!
    if (iproc == 0) then 
      call MPI_Bcast(fidx_all, 1, MPI_integer, 0, MPI_comm_world, ierr)
      buffer_tmp = transpose(fixed_points_all)
      call MPI_Bcast(buffer_tmp, fidx_all*3, MPI_real, 0, MPI_comm_world, ierr)   
    else
      call MPI_Bcast(fidx_all, 1, MPI_integer, 0, MPI_comm_world, ierr)
      call MPI_Bcast(buffer_tmp, fidx_all*3, MPI_real, 0, MPI_comm_world, ierr)
      fixed_points_all(1:fidx_all,:) = transpose(buffer_tmp(:,1:fidx_all))
    endif
!
    fidx = 0
    do j = 1,fidx_all
!     find the corresponding core
!     NB: the whole thing only works for nprocz = 1
      x_proc = floor((fixed_points_all(j,1)-x0)*real(nprocx)/Lx)
      y_proc = floor((fixed_points_all(j,2)-y0)*real(nprocy)/Ly)
      if ((x_proc == ipx) .and. (y_proc == ipy)) then
        fidx = fidx + 1
        fixed_points(fidx,:) = fixed_points_all(j,:)
      endif      
    enddo
!
!   this is kept for the moment, but will be removed soon
    write(str_tmp, "(I10.1,A)") iproc, '/fixed_points.dat'
    write(filename, *) 'data/proc', adjustl(trim(str_tmp))
    open(unit = 1, file = adjustl(trim(filename)), form = "unformatted", position = "append")
    write(1) tfixed_points_write
    write(1) float(fidx)
    do j=1,fidx
      write(1) fixed_points(j,:)
    enddo
    close(1)
!
    deallocate(tracer_tmp)
    deallocate(vv)
  end subroutine wfixed_points
!***********************************************************************
endmodule Fixed_point
