! $Id: boundcond.f90,v 1.2 2002-05-21 07:25:50 brandenb Exp $

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   boundcond.f90   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!  Module for boundary conditions. Extracted from (no)mpicomm, since
!!!  all non-periodic (external) boundary conditions require the same
!!!  code for serial and parallel runs.

module Boundcond

  implicit none

  contains

!***********************************************************************
    subroutine boundconds(f,errmesg)
!
!  Physical boundary conditions except for periodic stuff.
!  If errmesg is set to a non-empty string, the calling routine will call
!  stop_it(errmesg). This deferred abort mechanism is necessary because
!  Boundcond is used by Mpicomm, so we cannot call stop_it (from Mpicomm)
!  here.
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: tmp_xy
      real, dimension (7) :: lnrho
      real :: dlnrho
      integer :: i,j,k
      character (len=*) :: errmesg
!
      errmesg=""
!
!  Boundary conditions in x
!
      do j=1,mvar
        !
        ! `lower' bdry
        !
        select case(bcx1(j))
        case ('p')              ! periodic
          if (nprocx==1) then
            f(1:l1-1,:,:,j) = f(l2i:l2,:,:,j)
          endif
        case ('s')              ! symmetry
          do i=1,nghost; f(l1-i,:,:,j) = f(l1+i,:,:,j); enddo
        case ('a')              ! antisymmetry
          f(l1,:,:,j) = 0.      ! ensure bdry value=0 (indep.of initial cond.)
          do i=1,nghost; f(l1-i,:,:,j) = -f(l1+i,:,:,j); enddo
        case ('a2')             ! antisymmetry relative to boundary value
          do i=1,nghost; f(l1-i,:,:,j) = 2*f(l1,:,:,j)-f(l1+i,:,:,j); enddo
        case default
          if (lroot) &
               print*, "No such boundary condition bcx1 = ", &
                       bcx1(j), " for j=", j
          STOP
        endselect
        !
        ! `upper' bdry
        !
        select case(bcx2(j))
        case ('p')              ! periodic
          if (nprocx==1) then
            f(l2+1:mx,:,:,j) = f(l1:l1i,:,:,j)
          endif
        case ('s')              ! symmetry
          do i=1,nghost; f(l2+i,:,:,j) = f(l2-i,:,:,j); enddo
        case ('a')              ! antisymmetry
          f(l2,:,:,j) = 0.      ! ensure bdry value=0 (indep.of initial cond.)
          do i=1,nghost; f(l2+i,:,:,j) = -f(l2-i,:,:,j); enddo
        case ('a2')             ! antisymmetry relative to boundary value
          do i=1,nghost; f(l2+i,:,:,j) = 2*f(l2,:,:,j)-f(l2-i,:,:,j); enddo
        case default
          if (lroot) &
               print*, "No such boundary condition bcx2 = ", &
                       bcx2(j), " for j=", j
          STOP
        endselect
      enddo
!
!  Boundary conditions in y
!
      do j=1,mvar 
        !
        ! `lower' bdry
        !
        select case(bcy1(j))
        case ('p')              ! periodic
          if (nprocy==1) then
            f(:,1:m1-1,:,j) = f(:,m2i:m2,:,j)
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
          STOP
        endselect
        !
        ! `upper' bdry
        !
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
          STOP
        endselect
      enddo
!
!  Boundary conditions in z
!
      do j=1,mvar
        !
        ! `lower' bdry
        !
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
          ! handled below
        case default
          if (lroot) &
               print*, "No such boundary condition bcz1 = ", &
                       bcz1(j), " for j=", j
          STOP
        endselect
        !
        ! `upper' bdry
        !
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
        case default
          if (lroot) &
               print*, "No such boundary condition bcz2 = ", &
                       bcz2(j), " for j=", j
          STOP
        endselect
      enddo
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!
      if (lentropy) then
        if (bcz1(ient) == "c1") then
          if (bcz1(ilnrho) /= "a2") &
               errmesg = "BOUNDCONDS: Inconsistent boundary conditions 1."
          tmp_xy = gamma1/gamma & ! 1/T_0 (i.e. 1/T at boundary)
                   * exp(-gamma*f(:,:,n1,ient) - gamma1*f(:,:,n1,ilnrho))
          tmp_xy = Fheat/(hcond0*hcond1) * tmp_xy ! F_heat/(lambda T_0)
          do i=1,nghost
            f(:,:,n1-i,ient) = &
                 (2*i*dx*tmp_xy &
                  + 2*gamma1*(f(:,:,n1+i,ilnrho)-f(:,:,n1,ilnrho)) &
                 )/gamma &
                 + f(:,:,n1+i,ient)
          enddo
        endif
!
        if (bcz2(ient) == "c1") then
          if (bcz2(ilnrho) /= "a2") &
               errmesg = "BOUNDCONDS: Inconsistent boundary conditions 2."
          tmp_xy = gamma1/gamma & ! 1/T_0 (i.e. 1/T at boundary)
                   * exp(-gamma*f(:,:,n2,ient) - gamma1*f(:,:,n2,ilnrho))
          tmp_xy = Fheat/(hcond0*hcond2) * tmp_xy ! F_heat/(lambda T_0)
          do i=1,nghost
            f(:,:,n2+i,ient) = &
                 (-2*i*dx*tmp_xy &
                  + 2*gamma1*(f(:,:,n2-i,ilnrho)-f(:,:,n2,ilnrho)) &
                 )/gamma &
                 + f(:,:,n2-i,ient)
          enddo
        endif
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for
!  entropy and density.
!  NB: Sound speed is set to cs0, so this is mostly useful for top boundary.  
!
        if (bcz2(ient) == "c2") then
          if (bcz1(ilnrho) /= "a2") &
               errmesg = "BOUNDCONDS: Inconsistent boundary conditions 4."
          tmp_xy = (-gamma1*f(:,:,n2,ilnrho) + alog(cs20/gamma)) / gamma
          f(:,:,n2,ient) = tmp_xy
          do i=1,nghost
            f(:,:,n2+i,ient) = 2*tmp_xy - f(:,:,n2-i,ient)
          enddo
        endif
!
      endif                     ! (lentropy)
!
    endsubroutine boundconds
!***********************************************************************

endmodule Boundcond
