! $Id: boundcond.f90,v 1.30 2002-09-21 16:35:50 dobler Exp $

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
    subroutine boundconds_x(f)
!
!  Physical boundary conditions in x except for periodic stuff.
!  For the x-direction, the routine needs to be called immediately;
!  for the y- and z-directions is needs to be called after all
!  communication is done, because it needs to overwrite the otherwise
!  periodic boundary conditions that have already been solved under mpi.
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!
      use Cdata
      use Entropy
      use Magnetic
      use Radiation
!
      real, dimension (mx,my,mz,mvar) :: f
      integer :: i,j
!
      if(ldebug) print*,'ENTER: boundconds'
!
!  Boundary conditions in x
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
            case ('cT')         ! constant temp. (processed in own routine)
              if (j==ient) call bc_ss_temp_x(f,'bot')
            case ('sT')         ! symmetric temp. (processed in own routine)
              if (j==ient) call bc_ss_stemp_x(f,'bot')
            case ('in')
              if (j==ie) call bc_ee_inflow_x(f,'bot')
            case ('out')
              if (j==ie) call bc_ee_outflow_x(f,'bot')
            case ('db')
               call bc_db_x(f,'bot',j)
            case ('osc')
!               if (j==ilnrho) call bc_lnrho_osc_x(f,'bot')
               call bc_osc_x(f,'bot',j)
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
            case ('a2')            ! antisymmetry relative to boundary value
              do i=1,nghost; f(l2+i,:,:,j) = 2*f(l2,:,:,j)-f(l2-i,:,:,j); enddo
            case ('cT')            ! constant temp. (processed in own routine)
              if (j==ient) call bc_ss_temp_x(f,'top')
            case ('sT')            ! symmetric temp. (processed in own routine)
              if (j==ient) call bc_ss_stemp_x(f,'top')
            case ('in')
               if (j==ie) call bc_ee_inflow_x(f,'top')
            case ('out')
               if (j==ie) call bc_ee_outflow_x(f,'top')
            case ('db')
               call bc_db_x(f,'top',j)
            case ('osc')
!               if (j==ilnrho) call bc_lnrho_osc_x(f,'top')
               call bc_osc_x(f,'top',j)
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
    endsubroutine boundconds_x
!***********************************************************************
    subroutine boundconds_y(f)
!
!  Physical boundary conditions in y except for periodic stuff.
!  For the x-direction, the routine needs to be called immediately;
!  for the y- and z-directions is needs to be called after all
!  communication is done, because it needs to overwrite the otherwise
!  periodic boundary conditions that have already been solved under mpi.
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!
      use Cdata
      use Entropy
      use Magnetic
!
      real, dimension (mx,my,mz,mvar) :: f
      integer :: i,j
!
      if(ldebug) print*,'ENTER: boundconds'
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
          case ('cT')             ! constant temp. (processed in own routine)
            if (j==ient) call bc_ss_temp_y(f,'bot')
          case ('sT')             ! symmetric temp. (processed in own routine)
            if (j==ient) call bc_ss_stemp_y(f,'bot')
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
          case ('cT')             ! constant temp. (processed in own routine)
            if (j==ient) call bc_ss_temp_y(f,'top')
          case ('sT')             ! symmetric temp. (processed in own routine)
            if (j==ient) call bc_ss_stemp_y(f,'top')
          case default
            if (lroot) &
                 print*, "No such boundary condition bcy2 = ", &
                         bcy2(j), " for j=", j
            call stop_it("")
          endselect
        endif
      enddo
!
    endsubroutine boundconds_y
!***********************************************************************
    subroutine boundconds_z(f)
!
!  Physical boundary conditions in z except for periodic stuff.
!  For the x-direction, the routines needs to be called immediately;
!  for the y- and z-directions is needs to be called after all
!  communication is done, because it needs to overwrite the otherwise
!  periodic boundary conditions that have already been solved under mpi.
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!
      use Cdata
      use Entropy
      use Magnetic
      use Density
!
      real, dimension (mx,my,mz,mvar) :: f
      integer :: i,j
!
      if(ldebug) print*,'ENTER: boundconds'
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
            if (j==ient) call bc_ss_flux(f,'bot')
            if (j==iaa)  call bc_aa_pot(f,'bot')
          case ('cT')             ! constant temp. (processed in own routine)
            if (j==ient) call bc_ss_temp_z(f,'bot')
          case ('sT')             ! symmetric temp. (processed in own routine)
            if (j==ient) call bc_ss_stemp_z(f,'bot')
          case ('c2')             ! complex (processed in its own routine)
            if (j==ient) call bc_ss_temp_old(f,'bot')
          case ('db')             ! complex (processed in its own routine)
            call bc_db_z(f,'bot',j) 
          case ('ce')             ! complex (processed in its own routine) 
             if (j==ient) call bc_ss_energy(f,'bot')
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
            if (j==ient) call bc_ss_flux(f,'top')
            if (j==iaa)  call bc_aa_pot(f,'top')
          case ('cT')             ! constant temp. (processed in own routine)
            if (j==ient) call bc_ss_temp_z(f,'top')
          case ('sT')             ! symmetric temp. (processed in own routine)
            if (j==ient) call bc_ss_stemp_z(f,'top')
          case ('c2')             ! complex (processed in its own routine)
            if (j==ient) call bc_ss_temp_old(f,'top')
          case ('db')             ! complex (processed in its own routine)
             call bc_db_z(f,'top',j) 
          case ('ce')             ! complex (processed in its own routine)
             if (j==ient) call bc_ss_energy(f,'top')
          case default
            if (lroot) &
                 print*, "No such boundary condition bcz2 = ", &
                         bcz2(j), " for j=", j
            call stop_it("")
          endselect
        endif
      enddo
!
    endsubroutine boundconds_z
!***********************************************************************
    subroutine bc_db_z(f,topbot,j)
!
!  boundary condition for density
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the density module
!  13-aug-2002/nils: moved into boundcond
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: fder
      integer :: i,j
      !
      !  Set ghost zone to reproduce one-sided boundary condition 
      !  (2nd order):
      !  Finding the derivatives on the boundary using a one 
      !  sided final difference method. This derivative is being 
      !  used to calculate the boundary points. This will probably
      !  only be used for ln(rho)
      !
    select case(topbot)
!
! Bottom boundary
!
    case('bot')
       do i=1,nghost
          fder=(-3*f(:,:,n1-i+1,j)+4*f(:,:,n1-i+2,j)&
               -f(:,:,n1-i+3,j))/(2*dz)
          f(:,:,n1-i,j)=f(:,:,n1-i+2,j)-2*dz*fder
       end do
    case('top')
       do i=1,nghost
          fder=(3*f(:,:,n2+i-1,j)-4*f(:,:,n2+i-2,j)&
               +f(:,:,n2+i-3,j))/(2*dz)
          f(:,:,n2+i,j)=f(:,:,n2+i-2,j)+2*dz*fder
       end do
    case default
       if(lroot) print*,"invalid argument for 'bc_db_z'"
    endselect
!
  end subroutine bc_db_z
!*********************************************************************** 
    subroutine bc_db_x(f,topbot,j)
!
!  boundary condition for density
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the density module
!  13-aug-2002/nils: moved into boundcond
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (my,mz) :: fder
      integer :: i,j
      !
      !  Set ghost zone to reproduce one-sided boundary condition 
      !  (2nd order):
      !  Finding the derivatives on the boundary using a one 
      !  sided final difference method. This derivative is being 
      !  used to calculate the boundary points. This will probably
      !  only be used for ln(rho)
      !
      select case(topbot)
!
! Bottom boundary
!
      case('bot')
        do i=1,nghost
          fder=(-3*f(l1-i+1,:,:,j)+4*f(l1-i+2,:,:,j)&
               -f(l1-i+3,:,:,j))/(2*dx)
          f(l1-i,:,:,j)=f(l1-i+2,:,:,j)-2*dx*fder
        enddo
      case('top')
        do i=1,nghost
          fder=(3*f(l2+i-1,:,:,j)-4*f(l2+i-2,:,:,j)&
               +f(l2+i-3,:,:,j))/(2*dx)
          f(l2+i,:,:,j)=f(l2+i-2,:,:,j)+2*dx*fder
        enddo
      case default
        if(lroot) print*,"invalid argument for 'bc_db_x'"
      endselect
!
    endsubroutine bc_db_x
!***********************************************************************
    subroutine bc_osc_x(f,topbot,j)
!
!  12-aug-02/nils: coded
!  14-aug-02/nils: moved to boundcond
!
      use Cdata
      use Density
      use Hydro
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar) :: f
      real :: ampl_osc,frec
      integer :: i,pnts=10,j
!
      if (j==ilnrho) then
        ampl_osc=ampl_osc_lnrho
        frec=frec_lnrho
      elseif (j==iux) then
        ampl_osc=ampl_osc_ux
        frec=frec_ux
      else
        if(lroot) print*,"invalid argument for 'bc_osc_x'"
      endif
!         
      if (topbot=='bot') then
        do i=1,pnts
          f(i,:,:,j)=ampl_osc*sin(t*frec)*cos(2*pi*x(i)*mx/(Lx*pnts))
        enddo
      else
        do i=1,pnts
          f(mx+1-i,:,:,j)=ampl_osc*sin(t*frec)*cos(2*pi*x(mx+1-i)*mx/(pnts*Lx))
        enddo
      endif
    endsubroutine bc_osc_x
!*********************************************************************** 
    subroutine update_ghosts(a)
!
!  update all ghost zones of a
!  21-sep-02/wolf: extracted from wsnaps
!
      use Cparam
      use Mpicomm
!
      real, dimension (mx,my,mz,mvar) :: a
!
      call boundconds_x(a)
      call boundconds_y(a)
      call boundconds_z(a)
      call initiate_isendrcv_bdry(a)
      call finalise_isendrcv_bdry(a)
!
    endsubroutine update_ghosts
!***********************************************************************

endmodule Boundcond
