! $Id: boundcond.f90,v 1.48 2003-07-01 13:37:03 brandenb Exp $

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
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      call boundconds_x(f)      ! Do not change this order.
      call boundconds_y(f)
      call boundconds_z(f)      
!
    endsubroutine boundconds
!***********************************************************************
    subroutine boundconds_x(f)
!
!  Boundary conditions in x except for periodic part handled by communication.
!  boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-`corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!
      use Cdata
      use Entropy
      use Magnetic
      use Radiation
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j,k,ip_ok
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
!  Boundary conditions in x except for periodic part handled by communication.
!  boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-`corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!
      use Cdata
      use Entropy
      use Magnetic
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j,k,ip_ok
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
!  Boundary conditions in x except for periodic part handled by communication.
!  boundconds_x() needs to be called before communicating (because we
!  communicate the x-ghost points), boundconds_[yz] after communication
!  has finished (they need some of the data communicated for the edges
!  (yz-`corners').
!
!   8-jul-02/axel: split up into different routines for x,y and z directions
!  11-nov-02/wolf: unified bot/top, now handled by loop
!
      use Cdata
      use Entropy
      use Magnetic
      use Density
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j,k,ip_ok
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
              case ('a3')       ! a2 - wiggles
                call bc_asym3(f,topbot,j)
              case ('1s')        ! one-sided
                call bc_onesided_z(f,topbot,j)
              case ('c1')       ! complex
                if (j==ient) call bc_ss_flux(f,topbot)
                if (j==iaa)  call bc_aa_pot(f,topbot)
              case ('cT')       ! constant temp.
                if (j==ilnrho) call bc_lnrho_temp_z(f,topbot)
                if (j==ient)   call bc_ss_temp_z(f,topbot)
              case ('cp')       ! constant pressure
                if (j==ilnrho) call bc_lnrho_pressure_z(f,topbot)
              case ('sT')       ! symmetric temp.
                if (j==ient) call bc_ss_stemp_z(f,topbot)
              case ('c2')       ! complex
                if (j==ient) call bc_ss_temp_old(f,topbot)
              case ('db')       ! complex
                call bc_db_z(f,topbot,j) 
              case ('ce')       ! complex
                if (j==ient) call bc_ss_energy(f,topbot)
              case ('e1')       ! extrapolation
                call bc_extrap_2_1(f,topbot,j)
              case ('e2')       ! extrapolation
                call bc_extrap_2_2(f,topbot,j)
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
    subroutine bc_per_x(f,topbot,j)
!
!  periodic boundary condition
!  11-nov-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
    subroutine bc_asym3(f,topbot,j)
!
!  Generalized antisymmetric bc (a al `a2') with removal of Nyquist wiggles
!  Does not seem to help against wiggles -- use upwinding instead
!
!  TEMPORARY HACK: Commented put calculation of Nyquist, as this creates
!  problems for some 2D runs and this boundary condition was not really
!  helpful so far. Will either have to find a better solution or remove
!  this altogether. wd, 21-jun-2003
!
!  17-jun-03/wolf: coded
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my) :: Nyquist=impossible
      integer :: j
!
      select case(topbot)

      case('bot')               ! bottom boundary
        ! Nyquist = 0.25*(f(:,:,n1,j)-2*f(:,:,n1+1,j)+f(:,:,n1+2,j))
        ! Nyquist = 0.0625*(     f(:,:,n1  ,j)+f(:,:,n1+4,j) &
        !                   - 4*(f(:,:,n1+1,j)+f(:,:,n1+3,j)) &
        !                   + 6* f(:,:,n1+2,j) )
        f(:,:,n1-1,j) = 2*f(:,:,n1,j) - f(:,:,n1+1,j) -4*Nyquist
        f(:,:,n1-2,j) = 2*f(:,:,n1,j) - f(:,:,n1+2,j)
        f(:,:,n1-3,j) = 2*f(:,:,n1,j) - f(:,:,n1+3,j) -4*Nyquist

      case('top')               ! top boundary
        ! Nyquist = 0.25*(f(:,:,n2,j)-2*f(:,:,n2-1,j)+f(:,:,n2-2,j))
        ! Nyquist = 0.0625*(     f(:,:,n2  ,j)+f(:,:,n2-4,j) &
        !                   - 4*(f(:,:,n2-1,j)+f(:,:,n2-3,j)) &
        !                   + 6* f(:,:,n2-2,j) )
        f(:,:,n2+1,j) = 2*f(:,:,n2,j) - f(:,:,n2-1,j) -4*Nyquist
        f(:,:,n2+2,j) = 2*f(:,:,n2,j) - f(:,:,n2-2,j)
        f(:,:,n2+3,j) = 2*f(:,:,n2,j) - f(:,:,n2-3,j) -4*Nyquist

      case default
        if(lroot) print*, topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_asym3
!***********************************************************************
    subroutine bc_onesided_z(f,topbot,j)
!
!  One-sided conditions.
!  These expressions result from combining Eqs(207)-(210), astro-ph/0109497,
!  corresponding to (9.207)-(9.210) in Ferriz-Mas proceedings.
!
!   5-apr-03/axel: coded
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i,j,k
!
      select case(topbot)

      case('bot')               ! bottom boundary
        do i=1,nghost
          k=n1-i
          f(:,:,k,j)=7*f(:,:,k+1,j) &
                   -21*f(:,:,k+2,j) &
                   +35*f(:,:,k+3,j) &
                   -35*f(:,:,k+4,j) &
                   +21*f(:,:,k+5,j) &
                    -7*f(:,:,k+6,j) &
                      +f(:,:,k+7,j)
        enddo

      case('top')               ! top boundary
        do i=1,nghost
          k=n2+i
          f(:,:,k,j)=7*f(:,:,k-1,j) &
                   -21*f(:,:,k-2,j) &
                   +35*f(:,:,k-3,j) &
                   -35*f(:,:,k-4,j) &
                   +21*f(:,:,k-5,j) &
                    -7*f(:,:,k-6,j) &
                      +f(:,:,k-7,j)
        enddo

      case default
        if(lroot) print*, topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_onesided_z
!***********************************************************************
    subroutine bc_extrap_2_1(f,topbot,j)
!
!  Extrapolation boundary condition.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!   19-jun-03/wolf: coded
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j
!
      select case(topbot)

      case('bot')               ! bottom boundary
        f(:,:,n1-1,j)=0.25*(  9*f(:,:,n1,j)- 3*f(:,:,n1+1,j)- 5*f(:,:,n1+2,j)+ 3*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=0.05*( 81*f(:,:,n1,j)-43*f(:,:,n1+1,j)-57*f(:,:,n1+2,j)+39*f(:,:,n1+3,j))
        f(:,:,n1-3,j)=0.05*(127*f(:,:,n1,j)-81*f(:,:,n1+1,j)-99*f(:,:,n1+2,j)+73*f(:,:,n1+3,j))

      case('top')               ! top boundary
        f(:,:,n2+1,j)=0.25*(  9*f(:,:,n2,j)- 3*f(:,:,n2-1,j)- 5*f(:,:,n2-2,j)+ 3*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=0.05*( 81*f(:,:,n2,j)-43*f(:,:,n2-1,j)-57*f(:,:,n2-2,j)+39*f(:,:,n2-3,j))
        f(:,:,n2+3,j)=0.05*(127*f(:,:,n2,j)-81*f(:,:,n2-1,j)-99*f(:,:,n2-2,j)+73*f(:,:,n2-3,j))

      case default
        if(lroot) print*, topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_extrap_2_1
!***********************************************************************
    subroutine bc_extrap_2_2(f,topbot,j)
!
!  Extrapolation boundary condition.
!  Correct for polynomials up to 2nd order, determined 2 further degrees
!  of freedom by minimizing L2 norm of coefficient vector.
!
!   19-jun-03/wolf: coded
!    1-jul-03/axel: introduced abbreviations n1p4,n2m4
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j,n1p4,n2m4
!
!  abbreviations, because otherwise the ifc compiler complains
!  for 1-D runs without vertical extent
!
      n1p4=n1+4
      n2m4=n2-4
!
      select case(topbot)

      case('bot')               ! bottom boundary
        f(:,:,n1-1,j)=0.2   *(  9*f(:,:,n1,j)                 -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
        f(:,:,n1-2,j)=0.2   *( 15*f(:,:,n1,j)- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
        f(:,:,n1-3,j)=1./35.*(157*f(:,:,n1,j)-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))

      case('top')               ! top boundary
        f(:,:,n2+1,j)=0.2   *(  9*f(:,:,n2,j)                 -  4*f(:,:,n2-2,j)- 3*f(:,:,n2-3,j)+ 3*f(:,:,n2m4,j))
        f(:,:,n2+2,j)=0.2   *( 15*f(:,:,n2,j)- 2*f(:,:,n2-1,j)-  9*f(:,:,n2-2,j)- 6*f(:,:,n2-3,j)+ 7*f(:,:,n2m4,j))
        f(:,:,n2+3,j)=1./35.*(157*f(:,:,n2,j)-33*f(:,:,n2-1,j)-108*f(:,:,n2-2,j)-68*f(:,:,n2-3,j)+87*f(:,:,n2m4,j))

      case default
        if(lroot) print*, topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_extrap_2_2
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      real, dimension (mx,my,mz,mvar+maux) :: a
!
      call boundconds_x(a)
      call initiate_isendrcv_bdry(a)
      call finalise_isendrcv_bdry(a)
      call boundconds_y(a)
      call boundconds_z(a)
!
    endsubroutine update_ghosts
!***********************************************************************

endmodule Boundcond
