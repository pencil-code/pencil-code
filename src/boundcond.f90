! $Id: boundcond.f90,v 1.36 2002-11-12 16:07:58 dobler Exp $

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
!  Apply boundary conditions in all three directions.
!  Note that we _must_ call boundconds_{x,y,z} in this order, or edges and
!  corners will not be OK.
!
!  10-oct-02/wolf: coded
!
      use Cparam
!
      real, dimension (mx,my,mz,mvar) :: f
!
      call boundconds_x(f)      ! Do not change this order.
      call boundconds_y(f)
      call boundconds_z(f)      
!
    endsubroutine boundconds
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
      integer :: i,j,k,ip_ok
      character (len=bclen), dimension(mvar) :: bc12
      character (len=3) :: topbot
!
      if(ldebug) print*,'ENTER: boundconds_x'
!
      select case(nxgrid)
!
      case(1)
        if(headtt) print*,'no x-boundary'
!
!  Boundary conditions in x
!  shearing sheet boundary condition (default)
!  can still use other boundary conditions (even with shear)
!
      case default
        if (bcx1(1)=='she') then
          if (headtt) print*,'use shearing sheet boundary condition'
          call initiate_shearing(f)
          if (nprocy>1 .OR. (.NOT. lmpicomm)) call finalise_shearing(f)
        else
          do k=1,2                ! loop over 'bot','top'
            if (k==1) then
              topbot='bot'; bc12=bcx1; ip_ok=0
            else
              topbot='top'; bc12=bcx2; ip_ok=nprocx-1
            endif
            !
            do j=1,mvar
              if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcx',k,'(',j,')=',bc12(j)
              if (ipx == ip_ok) then
                select case(bc12(j))
                case ('p')        ! periodic
                  call bc_per_x(f,topbot,j)
                case ('s')        ! symmetry
                  call bc_sym_x(f,+1,topbot,j)
                case ('a')        ! antisymmetry
                  call bc_sym_x(f,-1,topbot,j)
                case ('a2')       ! antisymmetry relative to boundary value
                  call bc_sym_x(f,-1,topbot,j,REL=.true.)
                case ('cT')       ! constant temp.
                  if (j==ient) call bc_ss_temp_x(f,topbot)
                case ('sT')       ! symmetric temp.
                  if (j==ient) call bc_ss_stemp_x(f,topbot)
                case ('in')
                  if (j==ie) call bc_ee_inflow_x(f,topbot)
                case ('out')
                  if (j==ie) call bc_ee_outflow_x(f,topbot)
                case ('db')
                  call bc_db_x(f,topbot,j)
                case ('osc')
!                 if (j==ilnrho) call bc_lnrho_osc_x(f,topbot)
                  call bc_osc_x(f,topbot,j)
                case default
                  if (lroot) &
                       print*, "No such boundary condition bcx1/2 = ", &
                               bc12(j), " for j=", j
                  call stop_it("")
                endselect
              endif
            enddo
          enddo
        endif
      endselect
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
      integer :: i,j,k,ip_ok
      character (len=bclen), dimension(mvar) :: bc12
      character (len=3) :: topbot
!
      if(ldebug) print*,'ENTER: boundconds_y'
!
      select case(nygrid)
!
      case(1)
        if(headtt) print*,'no y-boundary'
!
!  Boundary conditions in y
!
      case default
        do k=1,2                ! loop over 'bot','top'
          if (k==1) then
            topbot='bot'; bc12=bcy1; ip_ok=0
          else
            topbot='top'; bc12=bcy2; ip_ok=nprocy-1
          endif
          !
          do j=1,mvar 
            if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcy',k,'(',j,')=',bc12(j)
            if (ipy == ip_ok) then
              select case(bc12(j))
              case ('p')        ! periodic
                call bc_per_y(f,topbot,j)
              case ('s')        ! symmetry
                call bc_sym_y(f,+1,topbot,j)
              case ('a')        ! antisymmetry
                call bc_sym_y(f,-1,topbot,j)
              case ('a2')       ! antisymmetry relative to boundary value
                call bc_sym_y(f,-1,topbot,j,REL=.true.)
              case ('cT')       ! constant temp.
                if (j==ient) call bc_ss_temp_y(f,topbot)
              case ('sT')       ! symmetric temp.
                if (j==ient) call bc_ss_stemp_y(f,topbot)
              case default
                if (lroot) &
                     print*, "No such boundary condition bcy1/2 = ", &
                             bc12(j), " for j=", j
                call stop_it("")
              endselect
            endif
          enddo
        enddo
      endselect
!
    endsubroutine boundconds_y
!***********************************************************************
    subroutine boundconds_z(f)
!
!  Physical boundary conditions in z except for periodic stuff.
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
      use Density
!
      real, dimension (mx,my,mz,mvar) :: f
      integer :: i,j,k,ip_ok
      character (len=bclen), dimension(mvar) :: bc12
      character (len=3) :: topbot
!
      if(ldebug) print*,'ENTER: boundconds_z'
!
      select case(nzgrid)
!
      case(1)
        if(headtt) print*,'no z-boundary'
!
!  Boundary conditions in z
!
      case default
        do k=1,2                ! loop over 'bot','top'
          if (k==1) then
            topbot='bot'; bc12=bcz1; ip_ok=0
          else
            topbot='top'; bc12=bcz2; ip_ok=nprocz-1
          endif
          !
          do j=1,mvar
            if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcz',k,'(',j,')=',bc12(j)
            if (ipz == ip_ok) then
              select case(bc12(j))
              case ('p')        ! periodic
                call bc_per_z(f,topbot,j)
              case ('s')        ! symmetry
                call bc_sym_z(f,+1,topbot,j)
              case ('a')        ! antisymmetry
                call bc_sym_z(f,-1,topbot,j)
              case ('a2')       ! antisymmetry relative to boundary value
                call bc_sym_z(f,-1,topbot,j,REL=.true.)
              case ('c1')       ! complex
                if (j==ient) call bc_ss_flux(f,topbot)
                if (j==iaa)  call bc_aa_pot(f,topbot)
              case ('cT')       ! constant temp.
                if (j==ilnrho) call bc_lnrho_temp_z(f,topbot)
                if (j==ient)   call bc_ss_temp_z(f,topbot)
              case ('sT')       ! symmetric temp.
                if (j==ient) call bc_ss_stemp_z(f,topbot)
              case ('c2')       ! complex
                if (j==ient) call bc_ss_temp_old(f,topbot)
              case ('db')       ! complex
                call bc_db_z(f,topbot,j) 
              case ('ce')       ! complex
                if (j==ient) call bc_ss_energy(f,topbot)
              case ('')         ! do nothing; assume that everything is set
              case default
                if (lroot) &
                     print*, "No such boundary condition bcz1/2 = ", &
                             bc12(j), " for j=", j
                call stop_it("")
              endselect
            endif
          enddo
        enddo
      endselect
!
    endsubroutine boundconds_z
!***********************************************************************
    subroutine boundconds_yz_corner(f)
!
!  Boundary condition for the yz corners. Needs to be called after the
!  communication has finished.
!
!   12-nov-02/wolf: coded
!
      use Cdata
      use Entropy
      use Magnetic
!
      real, dimension (mx,my,mz,mvar) :: f
      integer, dimension(2*nghost) :: idy,idz
      integer :: i,j,k,ip_ok
      character (len=bclen), dimension(mvar) :: bc12
      character (len=3) :: topbot
!
      if(ldebug) print*,'ENTER: boundconds_y'
!
!  Anything to be done?
!      
      if ((nygrid == 1) .or. (nzgrid ==1)) then
        if(headtt) print*,'no y- or z-boundary'
        return
        print*, '================== NEVER GOT HERE ====================='
      endif
!
!  Set up array indices for the yz corners
!
      do i=1,nghost
        idy(i)=i; idy(nghost+i)=nghost+ny+i
        idz(i)=i; idz(nghost+i)=nghost+nz+i
      enddo
!
!  Boundary conditions in y
!
      do k=1,2                ! loop over 'bot','top'
        if (k==1) then
          topbot='bot'; bc12=bcy1; ip_ok=0
        else
          topbot='top'; bc12=bcy2; ip_ok=nprocy-1
        endif
        !
        do j=1,mvar 
          if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcy',k,'(',j,')=',bc12(j)
          if (ipy == ip_ok) then
            select case(bc12(j))
            case ('p')        ! periodic
              call bc_per_y(f(:,:,idz,:),topbot,j)
            case ('s')        ! symmetry
              call bc_sym_y(f(:,:,idz,:),+1,topbot,j)
            case ('a')        ! antisymmetry
              call bc_sym_y(f(:,:,idz,:),-1,topbot,j)
            case ('a2')       ! antisymmetry relative to boundary value
              call bc_sym_y(f(:,:,idz,:),-1,topbot,j,REL=.true.)
            case ('cT')       ! constant temp.
              if (j==ient) call bc_ss_temp_y(f(:,:,idz,:),topbot)
            case ('sT')       ! symmetric temp.
              if (j==ient) call bc_ss_stemp_y(f(:,:,idz,:),topbot)
            case default
              if (lroot) &
                   print*, "No such boundary condition bcy1/2 = ", &
                   bc12(j), " for j=", j
              call stop_it("")
            endselect
          endif
        enddo
      enddo
!
!  Boundary conditions in z
!
      do k=1,2                ! loop over 'bot','top'
        if (k==1) then
          topbot='bot'; bc12=bcz1; ip_ok=0
        else
          topbot='top'; bc12=bcz2; ip_ok=nprocz-1
        endif
        !
        do j=1,mvar
          if (ldebug) write(*,'(A,I1,A,I2,A,A)') ' bcz',k,'(',j,')=',bc12(j)
          if (ipz == ip_ok) then
            select case(bc12(j))
            case ('p')        ! periodic
              call bc_per_z(f(:,idy,:,:),topbot,j)
            case ('s')        ! symmetry
              call bc_sym_z(f(:,idy,:,:),+1,topbot,j)
            case ('a')        ! antisymmetry
              call bc_sym_z(f(:,idy,:,:),-1,topbot,j)
            case ('a2')       ! antisymmetry relative to boundary value
              call bc_sym_z(f(:,idy,:,:),-1,topbot,j,REL=.true.)
            case ('c1')       ! complex
              if (j==ient) call bc_ss_flux(f(:,idy,:,:),topbot)
              if (j==iaa)  call bc_aa_pot(f(:,idy,:,:),topbot)
            case ('cT')       ! constant temp.
              if (j==ilnrho) call bc_lnrho_temp_z(f(:,idy,:,:),topbot)
              if (j==ient)   call bc_ss_temp_z(f(:,idy,:,:),topbot)
            case ('sT')       ! symmetric temp.
              if (j==ient) call bc_ss_stemp_z(f(:,idy,:,:),topbot)
            case ('c2')       ! complex
              if (j==ient) call bc_ss_temp_old(f(:,idy,:,:),topbot)
            case ('db')       ! complex
              call bc_db_z(f(:,idy,:,:),topbot,j) 
            case ('ce')       ! complex
              if (j==ient) call bc_ss_energy(f(:,idy,:,:),topbot)
            case ('')         ! do nothing; assume that everything is set
            case default
              if (lroot) &
                   print*, "No such boundary condition bcz1/2 = ", &
                   bc12(j), " for j=", j
              call stop_it("")
            endselect
          endif
        enddo
      enddo
!
    endsubroutine boundconds_yz_corner
!***********************************************************************
    subroutine bc_per_x(f,topbot,j)
!
!  periodic boundary condition
!  11-nov-02/wolf: coded
!
      use Cdata
!
      real, dimension(:,:,:,:) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case(topbot)

      case('bot')               ! bottom boundary
        if (nprocx==1) f(1:l1-1,:,:,j) = f(l2i:l2,:,:,j)

      case('top')               ! top boundary
        if (nprocx==1) f(l2+1:mx,:,:,j) = f(l1:l1i,:,:,j)

      case default
        if(lroot) print*, topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_per_x
!***********************************************************************
    subroutine bc_per_y(f,topbot,j)
!
!  periodic boundary condition
!  11-nov-02/wolf: coded
!
      use Cdata
!
      real, dimension (:,:,:,:) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case(topbot)

      case('bot')               ! bottom boundary
        if (nprocy==1) f(:,1:m1-1,:,j) = f(:,m2i:m2,:,j)

      case('top')               ! top boundary
        if (nprocy==1) f(:,m2+1:my,:,j) = f(:,m1:m1i,:,j)

      case default
        if(lroot) print*, topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_per_y
!***********************************************************************
    subroutine bc_per_z(f,topbot,j)
!
!  periodic boundary condition
!  11-nov-02/wolf: coded
!
      use Cdata
!
      real, dimension (:,:,:,:) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case(topbot)

      case('bot')               ! bottom boundary
        if (nprocz==1) f(:,:,1:n1-1,j) = f(:,:,n2i:n2,j)

      case('top')               ! top boundary
        if (nprocz==1) f(:,:,n2+1:mz,j) = f(:,:,n1:n1i,j)

      case default
        if(lroot) print*, topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_per_z
!***********************************************************************
    subroutine bc_sym_x(f,sgn,topbot,j,rel)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension(:,:,:,:) :: f
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif

      select case(topbot)

      case('bot')               ! bottom boundary
        if (relative) then
          do i=1,nghost; f(l1-i,:,:,j)=2*f(l1,:,:,j)+sgn*f(l1+i,:,:,j); enddo
        else
          do i=1,nghost; f(l1-i,:,:,j)=              sgn*f(l1+i,:,:,j); enddo
          if (sgn<0) f(l1,:,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case('top')               ! top boundary
        if (relative) then
          do i=1,nghost; f(l2+i,:,:,j)=2*f(l2,:,:,j)+sgn*f(l2-i,:,:,j); enddo
        else
          do i=1,nghost; f(l2+i,:,:,j)=              sgn*f(l2-i,:,:,j); enddo
          if (sgn<0) f(l2,:,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case default
        if(lroot) print*, topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_sym_x
!***********************************************************************
    subroutine bc_sym_y(f,sgn,topbot,j,rel)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif

      select case(topbot)

      case('bot')               ! bottom boundary
        if (relative) then
          do i=1,nghost; f(:,m1-i,:,j)=2*f(:,m1,:,j)+sgn*f(:,m1+i,:,j); enddo
        else
          do i=1,nghost; f(:,m1-i,:,j)=              sgn*f(:,m1+i,:,j); enddo
          if (sgn<0) f(:,m1,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case('top')               ! top boundary
        if (relative) then
          do i=1,nghost; f(:,m2+i,:,j)=2*f(:,m2,:,j)+sgn*f(:,m2-i,:,j); enddo
        else
          do i=1,nghost; f(:,m2+i,:,j)=              sgn*f(:,m2-i,:,j); enddo
          if (sgn<0) f(:,m2,:,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case default
        if(lroot) print*, topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_sym_y
!***********************************************************************
    subroutine bc_sym_z(f,sgn,topbot,j,rel)
!
!  Symmetry boundary conditions.
!  (f,-1,topbot,j)            --> antisymmetry             (f  =0)
!  (f,+1,topbot,j)            --> symmetry                 (f' =0)
!  (f,-1,topbot,j,REL=.true.) --> generalized antisymmetry (f''=0)
!  Don't combine rel=T and sgn=1, that wouldn't make much sense.
!
!  11-nov-02/wolf: coded
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: sgn,i,j
      logical, optional :: rel
      logical :: relative
!
      if (present(rel)) then; relative=rel; else; relative=.false.; endif

      select case(topbot)

      case('bot')               ! bottom boundary
        if (relative) then
          do i=1,nghost; f(:,:,n1-i,j)=2*f(:,:,n1,j)+sgn*f(:,:,n1+i,j); enddo
        else
          do i=1,nghost; f(:,:,n1-i,j)=              sgn*f(:,:,n1+i,j); enddo
          if (sgn<0) f(:,:,n1,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case('top')               ! top boundary
        if (relative) then
          do i=1,nghost; f(:,:,n2+i,j)=2*f(:,:,n2,j)+sgn*f(:,:,n2-i,j); enddo
        else
          do i=1,nghost; f(:,:,n2+i,j)=              sgn*f(:,:,n2-i,j); enddo
          if (sgn<0) f(:,:,n2,j) = 0. ! set bdry value=0 (indep of initcond)
        endif

      case default
        if(lroot) print*, topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_sym_z
!***********************************************************************
    subroutine bc_db_z(f,topbot,j)
!
!  ``One-sided'' boundary condition for density.
!  Set ghost zone to reproduce one-sided boundary condition 
!  (2nd order):
!  Finding the derivatives on the boundary using a one 
!  sided final difference method. This derivative is being 
!  used to calculate the boundary points. This will probably
!  only be used for ln(rho)
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the density module
!  13-aug-2002/nils: moved into boundcond
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension(:,:,:,:) :: f
      real, dimension (mx,my) :: fder
      integer :: i,j
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
        enddo
      case('top')
        do i=1,nghost
          fder=(3*f(:,:,n2+i-1,j)-4*f(:,:,n2+i-2,j)&
               +f(:,:,n2+i-3,j))/(2*dz)
          f(:,:,n2+i,j)=f(:,:,n2+i-2,j)+2*dz*fder
        enddo
      case default
        if(lroot) print*,"invalid argument for 'bc_db_z'"
      endselect
!
    endsubroutine bc_db_z
!*********************************************************************** 
    subroutine bc_db_x(f,topbot,j)
!
!  ``One-sided'' boundary condition for density.
!  Set ghost zone to reproduce one-sided boundary condition 
!  (2nd order):
!  Finding the derivatives on the boundary using a one 
!  sided final difference method. This derivative is being 
!  used to calculate the boundary points. This will probably
!  only be used for ln(rho)
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the density module
!  13-aug-2002/nils: moved into boundcond
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension(:,:,:,:) :: f
      real, dimension (my,mz) :: fder
      integer :: i,j
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
      real, dimension(:,:,:,:) :: f
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
      call boundconds(a)
      call initiate_isendrcv_bdry(a)
      call finalise_isendrcv_bdry(a)
!!!call boundconds_yz_corner(a)
!
    endsubroutine update_ghosts
!***********************************************************************

endmodule Boundcond
