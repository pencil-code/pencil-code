! $Id: boundcond.f90,v 1.20 2002-07-06 20:29:17 brandenb Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   boundcond.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  Module for boundary conditions. Extracted from (no)mpicomm, since
!!!  all non-periodic (external) boundary conditions require the same
!!!  code for serial and parallel runs.

module Boundcond

  use Mpicomm
 
  implicit none
  
  contains

!***********************************************************************
    subroutine boundconds(f)
!
!  Physical boundary conditions in x except for periodic stuff.
!
      use Cdata
      use Entropy
      use Magnetic
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: fder,cs2_2d
      integer :: i,j
!
      if(ldebug) print*,'ENTER: boundconds'
!
!  Boundary conditions in x
!
      !
      !  shearing sheet boundary condition (default)
      !  can still use other boundary conditions (even with shear)
      !
      if (bcx1(1)=='she') then
         if (headtt) print*,'use shearing sheet boundary condition'
         call initiate_shearing(f)
         if (nprocy>1 .OR. (.NOT. lmpicomm)) call finalise_shearing(f)
      else
        do j=1,mvar
          !
          ! `lower' bdry
          !
          if (ldebug) write(*,'(A,I2,A,A)') ' bcx1(',j,')=',bcx1(j)
          if (ipx == 0) then
            select case(bcx1(j))
            case ('p')          ! periodic
              if (nprocx==1) then
                f(1:l1-1,:,:,j) = f(l2i:l2,:,:,j)
                ! so far edges are copied twice, corners three times...
              endif
            case ('s')          ! symmetry
              do i=1,nghost; f(l1-i,:,:,j) = f(l1+i,:,:,j); enddo
            case ('a')          ! antisymmetry
              f(l1,:,:,j) = 0.  ! ensure bdry value=0 (indep.of initial cond.)
              do i=1,nghost; f(l1-i,:,:,j) = -f(l1+i,:,:,j); enddo
            case ('a2')         ! antisymmetry relative to boundary value
              do i=1,nghost; f(l1-i,:,:,j) = 2*f(l1,:,:,j)-f(l1+i,:,:,j); enddo
            case default
              if (lroot) &
                   print*, "No such boundary condition bcx1 = ", &
                           bcx1(j), " for j=", j
              call stop_it("")
            endselect
          endif
          !
          ! `upper' bdry
          !
          if (ldebug) write(*,'(A,I2,A,A)') ' bcx2(',j,')=',bcx2(j)
          if (ipx == nprocx-1) then
            select case(bcx2(j))
            case ('p')          ! periodic
              if (nprocx==1) then
                f(l2+1:mx,:,:,j) = f(l1:l1i,:,:,j)
                ! so far edges are copied twice, corners three times...
              endif
            case ('s')          ! symmetry
              do i=1,nghost; f(l2+i,:,:,j) = f(l2-i,:,:,j); enddo
            case ('a')          ! antisymmetry
              f(l2,:,:,j) = 0.  ! ensure bdry value=0 (indep.of initial cond.)
              do i=1,nghost; f(l2+i,:,:,j) = -f(l2-i,:,:,j); enddo
            case ('a2')             ! antisymmetry relative to boundary value
              do i=1,nghost; f(l2+i,:,:,j) = 2*f(l2,:,:,j)-f(l2-i,:,:,j); enddo
            case default
              if (lroot) &
                   print*, "No such boundary condition bcx2 = ", &
                           bcx2(j), " for j=", j
              call stop_it("")              
            endselect
          endif
        enddo
      endif
!
!  Boundary conditions in y
!
      do j=1,mvar 
        !
        ! `lower' bdry
        !
        if (ldebug) write(*,'(A,I2,A,A)') ' bcy1(',j,')=',bcy1(j)
        if (ipy == 0) then
          select case(bcy1(j))
          case ('p')              ! periodic
            if (nprocy==1) then
              f(:,1:m1-1,:,j) = f(:,m2i:m2,:,j)
              ! so far edges are copied twice, corners three times...
            endif
          case ('s')              ! symmetry
            do i=1,nghost; f(:,m1-i,:,j) = f(:,m1+i,:,j); enddo
          case ('a')              ! antisymmetry
            f(:,m1,:,j) = 0.      ! ensure bdry value=0 (indep.of initial cond.)
            do i=1,nghost; f(:,m1-i,:,j) = -f(:,m1+i,:,j); enddo
          case ('a2')             ! antisymmetry relative to boundary value
            do i=1,nghost; f(:,m1-i,:,j) = 2*f(:,m1,:,j)-f(:,m1+i,:,j); enddo
          case default
            if (lroot) &
                 print*, "No such boundary condition bcy1 = ", &
                         bcy1(j), " for j=", j
            call stop_it("")
          endselect
        endif
        !
        ! `upper' bdry
        !
        if (ldebug) write(*,'(A,I2,A,A)') ' bcy2(',j,')=',bcy2(j)
        if (ipy == nprocy-1) then
          select case(bcy2(j))
          case ('p')              ! periodic
            if (nprocy==1) then
              f(:,m2+1:my,:,j) = f(:,m1:m1i,:,j)
            endif
          case ('s')              ! symmetry
            do i=1,nghost; f(:,m2+i,:,j) = f(:,m2-i,:,j); enddo
          case ('a')              ! antisymmetry
            f(:,m2,:,j) = 0.      ! ensure bdry value=0 (indep.of initial cond.)
            do i=1,nghost; f(:,m2+i,:,j) = -f(:,m2-i,:,j); enddo
          case ('a2')             ! antisymmetry relative to boundary value
            do i=1,nghost; f(:,m2+i,:,j) = 2*f(:,m2,:,j)-f(:,m2-i,:,j); enddo
          case default
            if (lroot) &
                 print*, "No such boundary condition bcy2 = ", &
                         bcy2(j), " for j=", j
            call stop_it("")
          endselect
        endif
      enddo
!
!  Boundary conditions in z
!
      do j=1,mvar
        !
        ! `lower' bdry
        !
        if (ldebug) write(*,'(A,I2,A,A)') ' bcz1(',j,')=',bcz1(j)
        if (ipz == 0) then
          select case(bcz1(j))
          case ('p')              ! periodic
            if (nprocz==1) then
              f(:,:,1:n1-1,j) = f(:,:,n2i:n2,j)
            endif
          case ('s')              ! symmetry
            do i=1,nghost; f(:,:,n1-i,j) = f(:,:,n1+i,j); enddo
          case ('a')              ! antisymmetry
            f(:,:,n1,j) = 0.      ! ensure bdry value=0 (indep.of initial cond.)
            do i=1,nghost; f(:,:,n1-i,j) = -f(:,:,n1+i,j); enddo
          case ('a2')             ! antisymmetry relative to boundary value
            do i=1,nghost; f(:,:,n1-i,j) = 2*f(:,:,n1,j)-f(:,:,n1+i,j); enddo
          case ('c1')             ! complex (processed in its own routine)
          case ('c2')             ! complex (processed in its own routine)
            ! handled below
          case ('db')
            !
            !  Set ghost zone to reproduce one-sided boundary consition 
            !  (2nd order):
            !  Finding the derivatives on the boundary using a one 
            !  sided final difference method. This derivative is being 
            !  used to calculate the boundary points. This will probably
            !  only be used for ln(rho)
            !
            do i=1,nghost
               fder=(-3*f(:,:,n1-i+1,j)+4*f(:,:,n1-i+2,j)&
                    -f(:,:,n1-i+3,j))/(2*dz)
               f(:,:,n1-i,j)=f(:,:,n1-i+2,j)-2*dz*fder
            end do
         case ('ce')
            !
            !  This makes cs2 (temperature) constant on the boundaries
            !
            if (lentropy) then 
               cs2_2d=cs20*exp(gamma1*f(:,:,n2,ilnrho)+gamma*f(:,:,n2,ient))
            else
               cs2_2d=cs20;
            end if
            do i=1,nghost
               f(:,:,n1-i,j)=1./gamma*(-gamma1*f(:,:,n1-i,ilnrho)-log(cs20)&
                    +log(cs2_2d))
            end do
          case default
            if (lroot) &
                 print*, "No such boundary condition bcz1 = ", &
                         bcz1(j), " for j=", j
            call stop_it("")
          endselect
        endif
        !
        ! `upper' bdry
        !
        if (ldebug) write(*,'(A,I2,A,A)') ' bcz2(',j,')=',bcz2(j)
        if (ipz == nprocz-1) then
          select case(bcz2(j))
          case ('p')              ! periodic
            if (nprocz==1) then
              f(:,:,n2+1:mz,j) = f(:,:,n1:n1i,j)
            endif
          case ('s')              ! symmetry
            do i=1,nghost; f(:,:,n2+i,j) = f(:,:,n2-i,j); enddo
          case ('a')              ! antisymmetry
            f(:,:,n2,j) = 0.      ! ensure bdry value=0 (indep.of initial cond.)
            do i=1,nghost; f(:,:,n2+i,j) = -f(:,:,n2-i,j); enddo
          case ('a2')             ! antisymmetry relative to boundary value
            do i=1,nghost; f(:,:,n2+i,j) = 2*f(:,:,n2,j)-f(:,:,n2-i,j); enddo
          case ('c1')             ! complex (processed in its own routine)
          case ('c2')             ! complex (processed in its own routine)
            ! handled below
          case ('db')
            !
            !  Finding the derivatives on the boundary using a one 
            !  sided final difference method. This derivative is being 
            !  used to calculate the boundary points. This will probably
            !  only be used for ln(rho).
            !
            do i=1,nghost
               fder=(3*f(:,:,n2+i-1,j)-4*f(:,:,n2+i-2,j)&
                    +f(:,:,n2+i-3,j))/(2*dz)
               f(:,:,n2+i,j)=f(:,:,n2+i-2,j)+2*dz*fder
            end do
         case ('ce')
            !
            !  This makes cs2 (temperature) constant on the boundaries
            !
            if (lentropy) then 
               cs2_2d=cs20*exp(gamma1*f(:,:,n2,ilnrho)+gamma*f(:,:,n2,ient))
            else
               cs2_2d=cs20;
            end if
            do i=1,nghost
               f(:,:,n2+i,j)=1./gamma*(-gamma1*f(:,:,n2+i,ilnrho)-log(cs20)&
                    +log(cs2_2d))
            end do
          case default
            if (lroot) &
                 print*, "No such boundary condition bcz2 = ", &
                         bcz2(j), " for j=", j
              call stop_it("")
          endselect
        endif
      enddo
!
      if (lentropy)  call bc_ss(f)
      if (lmagnetic) call bc_aa(f)
!
    endsubroutine boundconds
!***********************************************************************
 
endmodule Boundcond
