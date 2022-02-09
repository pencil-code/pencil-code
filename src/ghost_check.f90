! $Id$
!
!  This module tests the ghost zones consistency.
!  For more documentation about the ghost zones, please refer to Mpicomm.f90.
!
module Ghost_check
!
  use Cdata
  use Cparam
!
  implicit none
!
  contains
!
!***********************************************************************
    function blocks_equal(msg,a,b)
!
!  Helper routine to check the equality of two arrays.
!
!  07-mar-2011/Bourdin.KIS: coded
!
      use Mpicomm, only: stop_it
!
      character (len=*), intent(in) :: msg
      real, dimension (:,:,:,:), intent(in) :: a, b
      logical :: blocks_equal
!
      integer :: num_x, num_y, num_z, num_a, px, py, pz, pa
!
      num_x = size (a, 1)
      num_y = size (a, 2)
      num_z = size (a, 3)
      num_a = size (a, 4)
!
      if (num_x /= size (b, 1)) call stop_it ('blocks_equal: size mismatch in X')
      if (num_y /= size (b, 2)) call stop_it ('blocks_equal: size mismatch in Y')
      if (num_z /= size (b, 3)) call stop_it ('blocks_equal: size mismatch in Z')
      if (num_a /= size (b, 4)) call stop_it ('blocks_equal: size mismatch in A')
!
      blocks_equal = .true.
      do pa = 1, min (num_a, mcom)
        do pz = 1, num_z
          do py = 1, num_y
            do px = 1, num_x
              if (a(px,py,pz,pa) /= b(px,py,pz,pa)) then
                write (100+iproc,*) msg, ' => ', px, py, pz, pa, ' : ', a(px,py,pz,pa), b(px,py,pz,pa)
                blocks_equal = .false.
              endif
            enddo
          enddo
        enddo
      enddo
!
    endfunction blocks_equal
!***********************************************************************
    subroutine check_ghosts_consistency(f, msg)
!
!  Helper routine to check the consistency of the ghost cell values.
!
!  07-mar-2011/Bourdin.KIS: coded
!
      use Mpicomm, only: stop_it, collect_xy, collect_z
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      character (len=*), intent(in) :: msg
!
      real, dimension (:,:,:,:), allocatable :: buffer, global, lower, middle, upper
      integer :: px, py, pz, lx, ux, ly, uy, lz, uz, alloc_err, mwx
      logical :: ok
!
      if (lpencil_check_at_work) return
!
      if (lroot) mwx=mx*nprocx    ! to avoid compile-time length calculation in allocation of "global" below

      if (lfirst_proc_xy) then
        allocate (buffer(mx*nprocx,my*nprocy,mz,mfarray), stat=alloc_err)
        if (alloc_err > 0) call stop_it ('check_ghosts_consistency: could not allocate buffer memory.')
      endif
      if (lroot) then
        allocate (global(mwx,my*nprocy,mz*nprocz,mfarray), stat=alloc_err)
        if (alloc_err > 0) call stop_it ('check_ghosts_consistency: could not allocate global memory.')
        allocate (lower(mx,my,mz,mfarray), middle(mx,my,mz,mfarray), upper(mx,my,mz,mfarray), stat=alloc_err)
        if (alloc_err > 0) call stop_it ('check_ghosts_consistency: could not allocate box memory.')
      endif
!
      call collect_xy (f, buffer)
      if (lfirst_proc_xy) call collect_z (buffer, global)
!
      if (lfirst_proc_xy) deallocate (buffer)
      if (.not. lroot) return
!
      ok = .true.
!
      ! check in X:
      write (100+iproc,*) 'Checking X-direction:'
      do py = 1, nprocy
        do pz = 1, nprocz
          if (nprocx == 1) then
            if (lperi(1)) then
              write (100+iproc,*) 'Subdomain ', 1, py, pz
              middle = global(:,(py-1)*my+1:py*my,(pz-1)*mz+1:pz*mz,:)
              ok = ok .and. blocks_equal ("X  l2i:l2 <>  1:l1-1", middle(l2i:l2,m1:m2,n1:n2,:), middle(1:l1-1,m1:m2,n1:n2,:))
              ok = ok .and. blocks_equal ("X l2+1:mx <> l1:l1i ", middle(l2+1:mx,m1:m2,n1:n2,:), middle(l1:l1i,m1:m2,n1:n2,:))
            endif
          else
            do px = 1, nprocx
              write (100+iproc,*) 'Subdomain ', px, py, pz
              lx = px - 1
              if (lperi(1) .and. (lx < 1)) lx = lx + nprocx
              ux = px + 1
              if (lperi(1) .and. (ux > nprocx)) ux = ux - nprocx
              middle = global((px-1)*mx+1:px*mx,(py-1)*my+1:py*my,(pz-1)*mz+1:pz*mz,:)
              ! check lower neighbor
              if (lx > 0) then
                lower = global((lx-1)*mx+1:lx*mx,(py-1)*my+1:py*my,(pz-1)*mz+1:pz*mz,:)
                ok = ok .and. blocks_equal ("X  l2i:l2 <>  1:l1-1", lower(l2i:l2,m1:m2,n1:n2,:), middle(1:l1-1,m1:m2,n1:n2,:))
              endif
              ! check upper neighbor
              if (ux <= nprocx) then
                upper = global((ux-1)*mx+1:ux*mx,(py-1)*my+1:py*my,(pz-1)*mz+1:pz*mz,:)
                ok = ok .and. blocks_equal ("X l2+1:mx <> l1:l1i ", middle(l2+1:mx,m1:m2,n1:n2,:), upper(l1:l1i,m1:m2,n1:n2,:))
              endif
            enddo
          endif
        enddo
      enddo
!
      ! check in Y:
      write (100+iproc,*) 'Checking Y-direction:'
      do px = 1, nprocx
        do pz = 1, nprocz
          if (nprocy == 1) then
            if (lperi(2)) then
              write (100+iproc,*) 'Subdomain ', px, 1, pz
              middle = global((px-1)*mx+1:px*mx,:,(pz-1)*mz+1:pz*mz,:)
              ok = ok .and. blocks_equal ("Y  m2i:m2 <>  1:m1-1", middle(l1:l2,m2i:m2,n1:n2,:), middle(l1:l2,1:m1-1,n1:n2,:))
              ok = ok .and. blocks_equal ("Y m2+1:my <> m1:m1i ", middle(l1:l2,m2+1:my,n1:n2,:), middle(l1:l2,m1:m1i,n1:n2,:))
            endif
          else
            do py = 1, nprocy
              write (100+iproc,*) 'Subdomain ', px, py, pz
              ly = ipy - 1
              if (lperi(2) .and. (ly < 1)) ly = ly + nprocy
              uy = ipy + 1
              if (lperi(2) .and. (uy > nprocy)) uy = uy - nprocy
              middle = global((px-1)*mx+1:px*mx,(py-1)*my+1:py*my,(pz-1)*mz+1:pz*mz,:)
              ! check lower neighbor
              if (ly > 0) then
                lower = global((px-1)*mx+1:px*mx,(ly-1)*my+1:ly*my,(pz-1)*mz+1:pz*mz,:)
                ok = ok .and. blocks_equal ("Y  m2i:m2 <>  1:m1-1", lower(l1:l2,m2i:m2,n1:n2,:), middle(l1:l2,1:m1-1,n1:n2,:))
              endif
              ! check upper neighbor
              if (uy <= nprocy) then
                upper = global((px-1)*mx+1:px*mx,(uy-1)*my+1:uy*my,(pz-1)*mz+1:pz*mz,:)
                ok = ok .and. blocks_equal ("Y m2+1:my <> m1:m1i ", middle(l1:l2,m2+1:my,n1:n2,:), upper(l1:l2,m1:m1i,n1:n2,:))
              endif
            enddo
          endif
        enddo
      enddo
!
      ! check in z:
      write (100+iproc,*) 'Checking Z-direction:'
      do px = 1, nprocx
        do py = 1, nprocy
          if (nprocz == 1) then
            if (lperi(3)) then
              write (100+iproc,*) 'Subdomain ', px, py, 1
              middle = global((px-1)*mx+1:px*mx,(py-1)*my+1:py*my,:,:)
              ok = ok .and. blocks_equal ("Z  n2i:n2 <>  1:n1-1", middle(l1:l2,m1:m2,n2i:n2,:), middle(l1:l2,m1:m2,1:n1-1,:))
              ok = ok .and. blocks_equal ("Z n2+1:mz <> n1:n1i ", middle(l1:l2,m1:m2,n2+1:mz,:), middle(l1:l2,m1:m2,n1:n1i,:))
            endif
          else
            do pz = 1, nprocz
              write (100+iproc,*) 'Subdomain ', px, py, pz
              lz = pz - 1
              if (lperi(3) .and. (lz < 1)) lz = lz + nprocz
              uz = pz + 1
              if (lperi(3) .and. (uz > nprocz)) uz = uz - nprocz
              middle = global((px-1)*mx+1:px*mx,(py-1)*my+1:py*my,(pz-1)*mz+1:pz*mz,:)
              ! check lower neighbor
              if (lz > 0) then
                lower = global((px-1)*mx+1:px*mx,(py-1)*my+1:py*my,(lz-1)*mz+1:lz*mz,:)
                ok = ok .and. blocks_equal ("Z  n2i:n2 <>  1:n1-1", lower(l1:l2,m1:m2,n2i:n2,:), middle(l1:l2,m1:m2,1:n1-1,:))
              endif
              ! check upper neighbor
              if (uz <= nprocz) then
                upper = global((px-1)*mx+1:px*mx,(py-1)*my+1:py*my,(uz-1)*mz+1:uz*mz,:)
                ok = ok .and. blocks_equal ("Z n2+1:mz <> n1:n1i ", middle(l1:l2,m1:m2,n2+1:mz,:), upper(l1:l2,m1:m2,n1:n1i,:))
              endif
            enddo
          endif
        enddo
      enddo
!
      if (.not. ok) then
        write (*,*) '=> ERROR: found inconsistency in ghost cells!'
        write (*,*) '=> SUBROUTINE: ', msg
        !call sleep (1)
        stop
      endif
!
      if (lroot) deallocate (global, lower, middle, upper)
!
    endsubroutine check_ghosts_consistency
!***********************************************************************
    subroutine check_ghosts_consistency_nompi(f,routine)
!
!  Non-MPI routine to check the consistency of all ghost cell values.
!  This routine will be removed, and the regular routine will be extended.
!
!  07-mar-2011/Bourdin.KIS: coded
!
      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      character (len=*), intent(in) :: routine
!
      logical :: ok
!
      if (lpencil_check_at_work) return
!
      ok = .true.
!
      ! check in X:
      if (lperi(1)) then
        ok = ok .and. blocks_equal ("X  l2i:l2 <>  1:l1-1", f(l2i:l2,m1:m2,n1:n2,:), f(1:l1-1,m1:m2,n1:n2,:))
        ok = ok .and. blocks_equal ("X l2+1:mx <> l1:l1i ", f(l2+1:mx,m1:m2,n1:n2,:), f(l1:l1i,m1:m2,n1:n2,:))
      endif
!
      ! check in Y:
      if (lperi(2)) then
        ok = ok .and. blocks_equal ("Y  m2i:m2 <>  1:m1-1", f(l1:l2,m2i:m2,n1:n2,:), f(l1:l2,1:m1-1,n1:n2,:))
        ok = ok .and. blocks_equal ("Y m2+1:my <> m1:m1i ", f(l1:l2,m2+1:my,n1:n2,:), f(l1:l2,m1:m1i,n1:n2,:))
      endif
!
      ! check in z:
      if (lperi(2)) then
        ok = ok .and. blocks_equal ("Z  n2i:n2 <>  1:n1-1", f(l1:l2,m1:m2,n2i:n2,:), f(l1:l2,m1:m2,1:n1-1,:))
        ok = ok .and. blocks_equal ("Z n2+1:mz <> n1:n1i ", f(l1:l2,m1:m2,n2+1:mz,:), f(l1:l2,m1:m2,n1:n1i,:))
      endif
!
      if (.not. ok) then
        write (*,*) '=> ERROR: found inconsistency in ghost cells!'
        write (*,*) '=> SUBROUTINE: ', routine
        !call sleep (1)
        stop
      endif
!
    endsubroutine check_ghosts_consistency_nompi
!***********************************************************************
endmodule Ghost_check
