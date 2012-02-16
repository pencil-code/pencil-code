! $Id: streamlines.f90 17621 2011-09-05 07:05:20Z iomsn $
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
!***************************************************************
module Streamlines
!
  use Cdata
  use Cparam
  use Messages
!
  implicit none
!
  contains
!  
!*********************************************************************** 
  subroutine get_grid_pos(phys_pos, grid_pos, outside)
!
! Determines the grid cell for the physical location.
!
! 13-feb-12/simon: coded
!
    real, dimension(3) :: phys_pos
    integer, dimension(3) :: grid_pos
    real :: delta
    integer :: j, outside
!
    intent(in) :: phys_pos
    intent(out) :: grid_pos
!
    outside = 0
!
    delta = Lx
    do j=1,nxgrid+2*nghost
      if (abs(phys_pos(1) - x(j)) < delta) then
        grid_pos(1) = j
        delta = abs(phys_pos(1) - x(j))
      endif
    enddo
!   check if the point lies outside the domain
    if (delta > dx) outside = 1
!      
    delta = Ly
    do j=1,nygrid+2*nghost
      if (abs(phys_pos(2) - y(j)) < delta) then
        grid_pos(2) = j
        delta = abs(phys_pos(2) - y(j))
      endif
    enddo
!   check if the point lies outside the domain
    if (delta > dy) outside = 1
!
    delta = Lz
    do j=1,nzgrid+2*nghost
      if (abs(phys_pos(3) - z(j)) < delta) then
        grid_pos(3) = j
        delta = abs(phys_pos(3) - z(j))
      endif
    enddo
!   check if the point lies outside the domain
    if (delta > dz) outside = 1
!
!   consider the processor indices
    grid_pos(1) = grid_pos(1) - nx*ipx - nghost
    grid_pos(2) = grid_pos(2) - ny*ipy - nghost
    grid_pos(3) = grid_pos(3) - nz*ipz - nghost
!
  endsubroutine get_grid_pos
!***********************************************************************
  subroutine trace_streamlines(f,tracers,n_tracers,h_max,h_min,l_max,tol)
!
!  trace stream lines of the vetor field stored in f(:,:,:,iaa)
!
!   13-feb-12/simon: coded
!
    use Sub
!
    real, dimension (mx,my,mz,mfarray) :: f
!     real, dimension(:,:), allocatable :: tracers
    real, dimension (nx*ny,7) :: tracers
!     real, dimension (1024*16,7) :: tracers
!     real, dimension (16384,7) :: tracers
    integer :: n_tracers, tracer_idx
!   the auxiliary magnetic field B
    real, dimension (nx,ny,nz,3) :: bb_aux
    real :: h_max, h_min, l_max, tol, dh, dist2
!   auxilliary vectors for the tracing
    real, dimension(3) :: x_mid, x_single, x_half, x_double
!   current position on the grid
    integer, dimension(3) :: grid_pos
    integer :: loop_count, outside = 0
!
    intent(in) :: f,n_tracers,h_max,h_min,l_max,tol
!
!     write(*,*) "allocating memory for the tracers"
!     allocate (tracers(n_tracers,7))
!
!   tracing stream lines
    write(*,*) "tracing ", n_tracers, " stream lines"
!
!   convert the magnetic vector potential into the magnetic field
    do m=1,ny
      do n=1,nz
        call curl(f,iaa,bb_aux(:,m,n,:))
      enddo
    enddo
!
!     do n=1,nz
!       do m=1,ny
!         do l=1,nx
!           write(*,*) x(l+nghost), y(m+nghost), z(n+nghost), bb_aux(l,m,n,:)
!         enddo
!       enddo
!     enddo
!
!
!   open the destination file
    open(unit = 1, file = "tracers.dat", form = "unformatted")
!
    do tracer_idx=1,n_tracers
      tracers(tracer_idx, 6) = 0.
!     initial step length dh
      dh = sqrt(h_max*h_min/2.0)
      loop_count = 0
!
!       call get_grid_pos(tracers(tracer_idx,3:5), grid_pos)
!       write(*,*) "tracer = ", tracers(tracer_idx,3:5), " grid_pos = ", grid_pos
      do
!       (a) Single step (midpoint method):
        call get_grid_pos(tracers(tracer_idx,3:5),grid_pos,outside)
        if (outside == 1) exit
        if (any(grid_pos <= 0) .or. (grid_pos(1) > nx) .or. &
            (grid_pos(2) > ny) .or. (grid_pos(3) > nz)) exit
        x_mid = tracers(tracer_idx,3:5) + 0.5*dh*bb_aux(grid_pos(1),grid_pos(2),grid_pos(3),:)
        call get_grid_pos(x_mid,grid_pos,outside)
        if (outside == 1) exit
        if (any(grid_pos <= 0) .or. (grid_pos(1) > nx) .or. &
            (grid_pos(2) > ny) .or. (grid_pos(3) > nz)) exit
        x_single = tracers(tracer_idx,3:5) + dh*bb_aux(grid_pos(1),grid_pos(2),grid_pos(3),:)
!
!       (b) Two steps with half stepsize:
        call get_grid_pos(tracers(tracer_idx,3:5),grid_pos,outside)
        if (outside == 1) exit
        if (any(grid_pos <= 0) .or. (grid_pos(1) > nx) .or. &
            (grid_pos(2) > ny) .or. (grid_pos(3) > nz)) exit
        x_mid = tracers(tracer_idx,3:5) + 0.25*dh*bb_aux(grid_pos(1),grid_pos(2),grid_pos(3),:)
        call get_grid_pos(x_mid,grid_pos,outside)
        if (outside == 1) exit
        if (any(grid_pos <= 0) .or. (grid_pos(1) > nx) .or. &
            (grid_pos(2) > ny) .or. (grid_pos(3) > nz)) exit
        x_half = tracers(tracer_idx,3:5) + 0.5*dh*bb_aux(grid_pos(1),grid_pos(2),grid_pos(3),:)
        call get_grid_pos(x_half,grid_pos,outside)
        if (outside == 1) exit
        if (any(grid_pos <= 0) .or. (grid_pos(1) > nx) .or. &
            (grid_pos(2) > ny) .or. (grid_pos(3) > nz)) exit
        x_mid = x_half + 0.25*dh*bb_aux(grid_pos(1),grid_pos(2),grid_pos(3),:)
        call get_grid_pos(x_mid,grid_pos,outside)
        if (outside == 1) exit
        if (any(grid_pos <= 0) .or. (grid_pos(1) > nx) .or. &
            (grid_pos(2) > ny) .or. (grid_pos(3) > nz)) exit
        x_double = x_half + 0.5*dh*bb_aux(grid_pos(1),grid_pos(2),grid_pos(3),:)
!
!       (c) Check error (difference between methods):
        dist2 = dot_product((x_single-x_double),(x_single-x_double))
        if (dist2 > tol**2) then
          dh = 0.5*dh
          if (abs(dh) <= h_min) then
!             write(*,*) 'Error: stepsize underflow'
            exit
          endif
        else
!           dist2 = sqrt(dot_product((x_double - tracers(tracer_idx,3:5)),(x_double - tracers(tracer_idx,3:5))))
          tracers(tracer_idx,3:5) = x_double
!           write(*,*) "tracer = ", tracers(tracer_idx,3:5), dh, loop_count
!           write(*,*) "bb_aux = ", bb_aux(grid_pos(1),grid_pos(2),grid_pos(3),:)
          tracers(tracer_idx, 6) = tracers(tracer_idx, 6) + dh
          if (abs(dh) < h_min) dh = 2*dh
!           if (dist2 < (tol/1000.)**2) dh = 2*dh
        endif
!
        if (tracers(tracer_idx, 6) >= l_max) then
!           write(*,*) "dh = ", dh, " line length = ", tracers(tracer_idx, 6), &
!               " loop_count = ", loop_count, " at grid_pos ", grid_pos
          exit
        endif
!
        loop_count = loop_count + 1
      enddo
!
!       if (any(grid_pos <= 0) .or. (grid_pos(1) > nx) .or. &
!           (grid_pos(2) > ny) .or. (grid_pos(3) > nz)) then
!         write(*,*) "tracer reached a boundary for tracer number ", &
!         tracer_idx, " loop_count = ", loop_count, " at grid_pos ", grid_pos
!       endif
      if (tracers(tracer_idx, 6) == 0) &
!       write(*,*) "length 0 at grid_pos ", grid_pos, " outside = ", outside
      write(1) tracers(tracer_idx,:)
      write(*,*) tracers(tracer_idx,:)
!
!     pass the current position to the next core
    enddo
    close(1)
!
  endsubroutine trace_streamlines
!***********************************************************************
endmodule Streamlines
