! $Id: nompicomm.f90,v 1.21 2002-05-08 17:47:28 dobler Exp $

!!!!!!!!!!!!!!!!!!!!!!!
!!!  nompicomm.f90  !!!
!!!!!!!!!!!!!!!!!!!!!!!

!!!  Module with dummy MPI stuff.
!!!  This allows code to be run on single cpu machine

module Mpicomm

  use Cparam
  use Cdata, only: lroot

  implicit none

  interface mpibcast_real       ! Overload the `mpibcast_real' function
    module procedure mpibcast_real_scl
    module procedure mpibcast_real_arr
  endinterface

  integer, dimension (ny*nz) :: mm,nn
  integer :: ierr,imn
  integer :: iproc,root=0
  integer :: ipx=0,ipy=0,ipz=0
  logical, dimension (nx*ny) :: necessary=.false.

  contains

!***********************************************************************
    subroutine mpicomm_init()
!
      use General
      use Cdata, only: lmpicomm,directory
!
!  sets iproc in order that we write in the correct directory
!
      character (len=4) :: chproc
      integer :: m,n
!
!  for single cpu machine, set processor to root value
!
      lmpicomm = .false.
      iproc=root
      lroot = .true.
!
!  produce index-array for the sequence of points to be worked through.
!  Trivial here, since no communication.
!
      imn=1
      do n=n1,n2
        do m=m1,m2
          mm(imn)=m
          nn(imn)=n
          imn=imn+1
        enddo
      enddo
!
!  produce a directory name for the output
!
      call chn(iproc,chproc)
      directory='tmp/proc'//chproc
!
    endsubroutine mpicomm_init
!***********************************************************************
    subroutine initiate_isendrcv_bdry(f)
!
      use Cdata
!
!  for one processor, use periodic boundary conditions
!  but in this dummy routine this is done in finalise_isendrcv_bdry
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my) :: tmp_xy
      real, dimension (7) :: lnrho
      real :: dlnrho
      integer :: i,j,k
!
!  Boundary conditions in x
!
      do j=1,mvar
        !
        ! `lower' bdry
        !
        select case(bcx1(j))
        case ('p')              ! periodic
          f( 1:l1-1,:,:,j) = f(l2i:l2,:,:,j)
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
          f(l2+1:mx,:,:,j) = f(l1:l1i,:,:,j)
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
          f(:, 1:m1-1,:,:) = f(:,m2i:m2,:,:)
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
          f(:,m2+1:my,:,:) = f(:,m1:m1i,:,:)
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
          f(:,:, 1:n1-1,j) = f(:,:,n2i:n2,j)
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
          f(:,:,n2+1:mz,j) = f(:,:,n1:n1i,j)
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
               call stop_it("Inconsistent boundary conditions 1.")
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
               call stop_it("Inconsistent boundary conditions 2.")
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
               call stop_it("Inconsistent boundary conditions 4.")
          tmp_xy = (-gamma1*f(:,:,n2,ilnrho) + alog(cs20/gamma)) / gamma
          f(:,:,n2,ient) = tmp_xy
          do i=1,nghost
            f(:,:,n2+i,ient) = 2*tmp_xy - f(:,:,n2-i,ient)
          enddo
        endif
!
      endif                     ! (lentropy)
!
    endsubroutine initiate_isendrcv_bdry
!***********************************************************************
    subroutine finalise_isendrcv_bdry(f)
!
      use General
!
!  do periodic boundary conditions
!
      real, dimension (mx,my,mz,mvar) :: f
      real :: dummy
!
!HPF$ ALIGN (:,:,:,*) WITH tmpl :: f
!HPF$ SHADOW(0,3,3,0) :: f
!
      dummy=f(1,1,1,1)  !(prevent compiler warning "unused variable ...")
!
    endsubroutine finalise_isendrcv_bdry
!***********************************************************************
    subroutine mpibcast_int(ibcast_array,nbcast_array)
!
      integer :: nbcast_array
      integer, dimension(nbcast_array) :: ibcast_array,dummy
!    
      dummy=ibcast_array
    endsubroutine mpibcast_int
!***********************************************************************
    subroutine mpibcast_real_scl(bcast_array,nbcast_array)
!
      integer :: nbcast_array
      real :: bcast_array,dummy
!
      if (nbcast_array/=1) stop "problem in mpibcast_real_scl"
      dummy=bcast_array
    endsubroutine mpibcast_real_scl
!***********************************************************************
    subroutine mpibcast_real_arr(bcast_array,nbcast_array)
!
      integer :: nbcast_array
      real, dimension(nbcast_array) :: bcast_array,dummy
!
      dummy=bcast_array
    endsubroutine mpibcast_real_arr
!***********************************************************************
    subroutine mpireduce_max(fmax_tmp,fmax,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fmax_tmp, fmax
!
      fmax=fmax_tmp
    endsubroutine mpireduce_max
!***********************************************************************
    subroutine mpireduce_sum(fsum_tmp,fsum,nreduce)
!
      integer :: nreduce
      real, dimension(nreduce) :: fsum_tmp,fsum
!
      fsum=fsum_tmp
    endsubroutine mpireduce_sum
!***********************************************************************
    subroutine mpifinalize
    endsubroutine mpifinalize
!***********************************************************************
    subroutine stop_it(msg)
!
!  Print message and stop
!  6-nov-01/wolf: coded
!
      character (len=*) :: msg
!      
      if (lroot) write(0,'(A,A)') 'STOPPED: ', msg
      call mpifinalize
      STOP
    endsubroutine stop_it
!***********************************************************************

endmodule Mpicomm
