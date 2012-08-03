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
  use Messages
  use Streamlines
!
  implicit none
!
! a few constants
  integer :: FINISHED_FIXED = 98
! the arrays with the values for x, y, and z for all cores (tracer xyz)
  real, pointer, dimension (:) :: xt, yt, zt
! fixed points array
  real, dimension (1000,3) :: fixed_points
! total number of fixed points of this core
  integer :: fidx
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
!
!  Output fixed points-data in 'tfixed_points' time intervals
!
    file = trim(datadir)//'/tfixed_points.dat'
    if (ifirst==0) then
      call read_snaptime(file,tfixed_points,ntracers,dtracers,t)
!     Read the previous fixed points from the file.
      write(str_tmp, "(I10.1,A)") iproc, '/fixed_points.dat'
      write(filename, *) 'data/proc', adjustl(trim(str_tmp))
      open(unit = 1, file = adjustl(trim(filename)), form = "unformatted")
!     loop until we find the last entry
      IOstatus = 0
!
      read(1,iostat = IOstatus) tfixed_points_write
      if (IOstatus == 0) then
        read(1) fidx_read
        fidx = int(fidx_read)
        allocate(fixed_tmp(fidx,3))
        do j=1,fidx
          read(1) fixed_tmp(j,:)
        enddo
        fixed_points(1:fidx,:) = fixed_tmp(:,:)
        deallocate(fixed_tmp)
        close(1)
      endif
      ifirst=1
    endif
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
!
    real, dimension (mx,my,mz,mfarray) :: f
    real :: point(2), point_old(2), fixed_point(2)
!   the integrated quantity along the field line
    real :: q
!     real, pointer, dimension(2) :: fixed_point
    real, pointer, dimension (:,:,:,:) :: vv
!   the points for the extrapolation
    real :: der(2)
    real, pointer, dimension (:,:) :: tracers
    real :: diff(13)
    integer :: iter
    real :: lambda, dl
!
    intent(out) :: fixed_point, q
!
    allocate(tracers(13,7))
!
    iter = 0
    dl = min(dx,dy)/30.
    do
!     trace the necessary field lines for the gradient
      tracers(1,:) = (/point(1),point(2),point(1),point(2),z(1+nghost)-ipz*nz*dz+dz,0.,1./)
      tracers(2,:) = (/point(1)-3*dl,point(2),point(1)-3*dl,point(2),z(1+nghost)-ipz*nz*dz+dz,0.,2./)
      tracers(3,:) = (/point(1)-2*dl,point(2),point(1)-2*dl,point(2),z(1+nghost)-ipz*nz*dz+dz,0.,2./)
      tracers(4,:) = (/point(1)-1*dl,point(2),point(1)-1*dl,point(2),z(1+nghost)-ipz*nz*dz+dz,0.,2./)
      tracers(5,:) = (/point(1)+1*dl,point(2),point(1)+1*dl,point(2),z(1+nghost)-ipz*nz*dz+dz,0.,2./)
      tracers(6,:) = (/point(1)+2*dl,point(2),point(1)+2*dl,point(2),z(1+nghost)-ipz*nz*dz+dz,0.,2./)
      tracers(7,:) = (/point(1)+3*dl,point(2),point(1)+3*dl,point(2),z(1+nghost)-ipz*nz*dz+dz,0.,2./)
      tracers(8,:) = (/point(1),point(2)-3*dl,point(1),point(2)-3*dl,z(1+nghost)-ipz*nz*dz+dz,0.,4./)
      tracers(9,:) = (/point(1),point(2)-2*dl,point(1),point(2)-2*dl,z(1+nghost)-ipz*nz*dz+dz,0.,4./)
      tracers(10,:) = (/point(1),point(2)-1*dl,point(1),point(2)-1*dl,z(1+nghost)-ipz*nz*dz+dz,0.,4./)
      tracers(11,:) = (/point(1),point(2)+1*dl,point(1),point(2)+1*dl,z(1+nghost)-ipz*nz*dz+dz,0.,4./)
      tracers(12,:) = (/point(1),point(2)+2*dl,point(1),point(2)+2*dl,z(1+nghost)-ipz*nz*dz+dz,0.,4./)
      tracers(13,:) = (/point(1),point(2)+3*dl,point(1),point(2)+3*dl,z(1+nghost)-ipz*nz*dz+dz,0.,4./)
      call trace_streamlines(f,tracers,13,vv)
      diff(1:13) = sqrt((tracers(1:13,3) - tracers(1:13,1))**2 + (tracers(1:13,4) - tracers(1:13,2))**2)
      q = tracers(1,7)
!     determine the 6th order gradient at this point
      der(1) = (-diff(2)/60. + 3./20.*diff(3) - 3./4.*diff(4) + 3./4.*diff(5) - 3./20.*diff(6) +diff(7)/60.)/dl
      der(2) = (-diff(8)/60. + 3./20.*diff(9) - 3./4.*diff(10) + 3./4.*diff(11) - 3./20.*diff(12) +diff(13)/60.)/dl
!     if the gradient is 0 we have already reached the fixed points
      if ((der(1) == 0 .and. der(2) == 0) .or. (der(1)**2+der(2)**2 == 0)) then
        fixed_point = point
        exit
      else
!       gradient descent method
        point_old = point
        lambda = diff(1)/(der(1)**2+der(2)**2)
!       avoid that the step length becomes too large
        if (lambda > min(dx,dy)) lambda = min(dx,dy)
        point = point - der*lambda
      endif
      if (iter > 50) then
        fixed_point = point
        exit
      endif
      iter = iter + 1
    enddo
!
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
!  trace stream lines of the vetor field stored in f(:,:,:,iaa)
!
!   13-feb-12/simon: coded
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
        if (abs(poincare) >= 6.2) then
          call get_fixed_point(f,(/(tracers2(j+(l-1)*(nx*trace_sub+addx),1)+tracers2(j+1+(l-1)*(nx*trace_sub+addx),1))/2., &
              (tracers2(j+(l-1)*(nx*trace_sub+addx),2)+tracers2(j+l*(nx*trace_sub+addx),2))/2./), &
              fixed_points(fidx,1:2), fixed_points(fidx,3), vv)
! !  Fixed point lies outside the physical domain.
!           if ((fixed_points(fidx,1) < xt(1)) .or. (fixed_points(fidx,1) > xt(nxgrid*trace_sub)) .or. &
!               (fixed_points(fidx,2) < yt(1)) .or. (fixed_points(fidx,2) > yt(nygrid*trace_sub))) then
!             write(*,*) iproc, "warning: fixed point lies outside the domain"
!  Fixed point lies outside the cell.
          if ((fixed_points(fidx,1) < tracers2(j+(l-1)*(nx*trace_sub+addx),1)) .or. &
              (fixed_points(fidx,1) > tracers2(j+1+(l-1)*(nx*trace_sub+addx),1)) .or. &
              (fixed_points(fidx,2) < tracers2(j+(l-1)*(nx*trace_sub+addx),2)) .or. &
              (fixed_points(fidx,2) > tracers2(j+l*(nx*trace_sub+addx),2))) then
            write(*,*) iproc, "warning: fixed point lies outside the domain"
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
    integer :: j, flag
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
!
!   allocate memory for the traced field
    allocate(vv(nx,ny,nz,3))
!
    allocate(tracer_tmp(1,7))
!
    call keep_compiler_quiet(path)
!
! !   mpi debug output for finding dead locks
!     write(str_tmp, "(I10.1,A)") iproc, '/mpi_debug.out'
!     write(filename, *) 'data/proc', adjustl(trim(str_tmp))
!     open(unit = 11, file = adjustl(trim(filename)), access="sequential", form="formatted")
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
! !   mpi debug output for finding dead locks
!     rewind(11)
!     write(11,*) iproc, "fixed_points: starting fixedp oint finding at", fixed_points(:,1:2)
    do j=1,fidx
      point = fixed_points(j,1:2)
      call get_fixed_point(f, point, fixed_points(j,1:2), fixed_points(j,3), vv)
    enddo
!
!   Tell every other core that we have finished.
! !   mpi debug output for finding dead locks
!     rewind(11)
!     write(11,*) iproc, "fixed_points: we have finished"
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
! !   mpi debug output for finding dead locks
!     rewind(11)
!     write(11,*) iproc, "fixed_points: wait for other cores"
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
! !   mpi debug output for finding dead locks
!     rewind(11)
!     write(11,*) iproc, "fixed_points: other cores have finished"
!
!   Wait for other cores. This ensures that uncomplete fixed points wont get written out.
    call MPI_BARRIER(MPI_comm_world, ierr)
! !   mpi debug output for finding dead locks
!     rewind(11)
!     write(11,*) iproc, "fixed_points: MPI_BARRIER passed"
!
    write(str_tmp, "(I10.1,A)") iproc, '/fixed_points.dat'
    write(filename, *) 'data/proc', adjustl(trim(str_tmp))
    open(unit = 1, file = adjustl(trim(filename)), form = "unformatted", position = "append")
    write(1) tfixed_points_write
    write(1) float(fidx)
    do j=1,fidx
      write(1) fixed_points(j,:)
    enddo
    close(1)
!     close(11)
!
    deallocate(tracer_tmp)
    deallocate(vv)
  end subroutine wfixed_points
!***********************************************************************
endmodule Fixed_point
