! $Id: boundcond.f90,v 1.62 2004-04-26 16:05:16 dobler Exp $

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
      use Ionization
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j,k,ip_ok
      character (len=bclen), dimension(mvar) :: bc12
      character (len=3) :: topbot
!
      if(ldebug) print*,'boundconds_x: ENTER: boundconds_x'
!
      select case(nxgrid)
!
      case(1)
        if(headtt) print*,'boundconds_x: no x-boundary'
!
!  Boundary conditions in x
!  shearing sheet boundary condition (default)
!  can still use other boundary conditions (even with shear)
!
      case default
        if (bcx1(1)=='she') then
          if (ip<12.and.headtt) print*, &
               'boundconds_x: use shearing sheet boundary condition'
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
                  if (j==iss) call bc_ss_temp_x(f,topbot)
                case ('sT')       ! symmetric temp.
                  if (j==iss) call bc_ss_stemp_x(f,topbot)
                case ('in')
                  if (j==ie) call bc_ee_inflow_x(f,topbot)
                case ('out')
                  if (j==ie) call bc_ee_outflow_x(f,topbot)
                case ('db')
                  call bc_db_x(f,topbot,j)
                case ('1')        ! f=1 (for debugging)
                  call bc_one_x(f,topbot,j)
                case default
                  print*, &
                       "boundconds_x: No such boundary condition bcx1/2 = ", &
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
      use Ionization
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j,k,ip_ok
      character (len=bclen), dimension(mvar) :: bc12
      character (len=3) :: topbot
!
      if(ldebug) print*,'boundconds_y: ENTER: boundconds_y'
!
      select case(nygrid)
!
      case(1)
        if(headtt) print*,'boundconds_y: no y-boundary'
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
                if (j==iss) call bc_ss_temp_y(f,topbot)
              case ('sT')       ! symmetric temp.
                if (j==iss) call bc_ss_stemp_y(f,topbot)
              case ('1')        ! f=1 (for debugging)
                call bc_one_y(f,topbot,j)
              case default
                print*, "boundconds_y: No such boundary condition bcy1/2 = ", &
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
      use Entropy, only: hcond0,hcond1,Fbot,FbotKbot,Ftop,FtopKtop,chi, &
                         lmultilayer,lcalc_heatcond_constchi
      use Magnetic
      use Density
      use Ionization
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real :: Ftopbot,FtopbotK
      integer :: j,k,ip_ok
      character (len=bclen), dimension(mvar) :: bc12
      character (len=3) :: topbot
!
      if(ldebug) print*,'boundconds_z: ENTER: boundconds_z'
!
      select case(nzgrid)
!
      case(1)
        if(headtt) print*,'boundconds_z: no z-boundary'
!
!  Boundary conditions in z
!
      case default
        do k=1,2                ! loop over 'bot','top'
          if (k==1) then
            topbot='bot'
            bc12=bcz1
            ip_ok=0
            Ftopbot=Fbot
            FtopbotK=FbotKbot
          else
            topbot='top'
            bc12=bcz2
            ip_ok=nprocz-1
            Ftopbot=Ftop
            FtopbotK=FtopKtop
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
              case ('v3')       ! vanishing third derivative
                call bc_van3rd_z(f,topbot,j)
              case ('1s')       ! one-sided
                call bc_onesided_z(f,topbot,j)
              case ('c1')       ! complex
                if (j==iss) call bc_ss_flux(f,topbot,hcond0,hcond1,Ftopbot,FtopbotK,chi, &
                                  lmultilayer,lcalc_heatcond_constchi)
                if (j==iaa) call bc_aa_pot(f,topbot)
              case ('cT')       ! constant temp.
                if (j==ilnrho) call bc_lnrho_temp_z(f,topbot)
                if (j==iss)   call bc_ss_temp_z(f,topbot)
              case ('cp')       ! constant pressure
                if (j==ilnrho) call bc_lnrho_pressure_z(f,topbot)
              case ('sT')       ! symmetric temp.
                if (j==iss) call bc_ss_stemp_z(f,topbot)
              case ('c2')       ! complex
                if (j==iss) call bc_ss_temp_old(f,topbot)
              case ('db')       ! complex
                call bc_db_z(f,topbot,j) 
              case ('ce')       ! complex
                if (j==iss) call bc_ss_energy(f,topbot)
              case ('e1')       ! extrapolation
                call bc_extrap_2_1(f,topbot,j)
              case ('e2')       ! extrapolation
                call bc_extrap_2_2(f,topbot,j)
              case ('b1')       ! extrapolation with zero value (improved 'a')
                call bc_extrap0_2_0(f,topbot,j)
              case ('b2')       ! extrapolation with zero value (improved 'a')
                call bc_extrap0_2_1(f,topbot,j)
              case ('b3')       ! extrapolation with zero value (improved 'a')
                call bc_extrap0_2_2(f,topbot,j)
              case ('f')        ! freeze value
                ! tell other modules not to change boundary value
                call bc_freeze_var_z(topbot,j)
                call bc_sym_z(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('fB')       ! frozen-in B-field
                ! tell other modules not to change boundary value
                call bc_frozen_in_bb_z(topbot)
                call bc_sym_z(f,-1,topbot,j,REL=.true.) ! antisymm wrt boundary
              case ('g')        ! set to given value(s) or function
                call bc_force_z(f,topbot,j)
              case ('1')        ! f=1 (for debugging)
                call bc_one_z(f,topbot,j)
              case ('')         ! do nothing; assume that everything is set
              case default
                print*, "boundconds_z: No such boundary condition bcz1/2 = ", &
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
        print*, "bc_per_x: ", topbot, " should be `top' or `bot'"

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
        print*, "bc_per_y: ", topbot, " should be `top' or `bot'"

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
        print*, "bc_per_z: ", topbot, " should be `top' or `bot'"

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
        print*, "bc_sym_x: ", topbot, " should be `top' or `bot'"

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
        print*, "bc_sym_y: ", topbot, " should be `top' or `bot'"

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
        print*, "bc_sym_z: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_sym_z
!***********************************************************************
    subroutine bc_van3rd_z(f,topbot,j)
!
!  Boundary condition with vanishing 3rd derivative
!  (useful for vertical hydrostatic equilibrium in discs)
!
!  19-aug-03/anders: coded
!
    use Cdata
!
    character (len=3) :: topbot
    real, dimension (mx,my,mz,mvar+maux) :: f
    real, dimension (mx,my) :: cpoly0,cpoly1,cpoly2
    integer :: i,j

    select case(topbot)

    case('bot')
      cpoly0(:,:)=f(:,:,n1,j)
      cpoly1(:,:)=-(3*f(:,:,n1,j)-4*f(:,:,n1+1,j)+f(:,:,n1+2,j))/(2*dz)
      cpoly2(:,:)=-(-f(:,:,n1,j)+2*f(:,:,n1+1,j)-f(:,:,n1+2,j)) /(2*dz**2)
      do i=1,nghost
        f(:,:,n1-i,j) = cpoly0(:,:) - cpoly1(:,:)*i*dz + cpoly2(:,:)*(i*dz)**2
      enddo

    case('top')
      cpoly0(:,:)=f(:,:,n2,j)
      cpoly1(:,:)=-(-3*f(:,:,n2,j)+4*f(:,:,n2-1,j)-f(:,:,n2-2,j))/(2*dz)
      cpoly2(:,:)=-(-f(:,:,n2,j)+2*f(:,:,n2-1,j)-f(:,:,n2-2,j))/(2*dz**2)
      do i=1,nghost
        f(:,:,n2+i,j) = cpoly0(:,:) + cpoly1(:,:)*i*dz + cpoly2(:,:)*(i*dz)**2
      enddo

    endselect

    endsubroutine bc_van3rd_z
!***********************************************************************
    subroutine bc_asym3(f,topbot,j)
!
!  Generalized antisymmetric bc (a al `a2') with removal of Nyquist wiggles
!  Does not seem to help against wiggles -- use upwinding instead
!
!  TEMPORARY HACK: Commented out calculation of Nyquist, as this creates
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
        print*, "bc_asym3: ", topbot, " should be `top' or `bot'"

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
        print*, "bc_onesided_z ", topbot, " should be `top' or `bot'"

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
        print*, "bc_extrap_2_1: ", topbot, " should be `top' or `bot'"

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
        print*, "bc_extrap_2_2: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_extrap_2_2
!***********************************************************************
    subroutine bc_extrap0_2_0(f,topbot,j)
!
!  Extrapolation boundary condition for f(bdry)=0.
!  Correct for polynomials up to 2nd order, determined no further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!    9-oct-03/wolf: coded
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j
!
      select case(topbot)

!       case('bot')               ! bottom boundary
!         f(:,:,n1  ,j)= 0.       ! set bdry value=0 (indep of initcond)
!         f(:,:,n1-1,j)=- 3*f(:,:,n1+1,j)+  f(:,:,n1+2,j)
!         f(:,:,n1-2,j)=- 8*f(:,:,n1+1,j)+3*f(:,:,n1+2,j)
!         f(:,:,n1-3,j)=-15*f(:,:,n1+1,j)+6*f(:,:,n1+2,j)

!       case('top')               ! top boundary
!         f(:,:,n2  ,j)= 0.       ! set bdry value=0 (indep of initcond)
!         f(:,:,n2+1,j)=- 3*f(:,:,n2-1,j)+  f(:,:,n2-2,j)
!         f(:,:,n2+2,j)=- 8*f(:,:,n2-1,j)+3*f(:,:,n2-2,j)
!         f(:,:,n2+3,j)=-15*f(:,:,n2-1,j)+6*f(:,:,n2-2,j)


!! Nyquist-filtering
      case('bot')               ! bottom boundary
        f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=(1/11.)*(-17*f(:,:,n1+1,j)- 9*f(:,:,n1+2,j)+ 8*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=      2*(- 2*f(:,:,n1+1,j)-   f(:,:,n1+2,j)+   f(:,:,n1+3,j))
        f(:,:,n1-3,j)=(3/11.)*(-27*f(:,:,n1+1,j)-13*f(:,:,n1+2,j)+14*f(:,:,n1+3,j))

      case('top')               ! top boundary
        f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=(1/11.)*(-17*f(:,:,n2-1,j)- 9*f(:,:,n2-2,j)+ 8*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=      2*(- 2*f(:,:,n2-1,j)-   f(:,:,n2-2,j)+   f(:,:,n2-3,j))
        f(:,:,n2+3,j)=(3/11.)*(-27*f(:,:,n2-1,j)-13*f(:,:,n2-2,j)+14*f(:,:,n2-3,j))


! !! Nyquist-transparent
!       case('bot')               ! bottom boundary
!         f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
!         f(:,:,n1-1,j)=(1/11.)*(-13*f(:,:,n1+1,j)-14*f(:,:,n1+2,j)+10*f(:,:,n1+3,j))
!         f(:,:,n1-2,j)=(1/11.)*(-48*f(:,:,n1+1,j)-17*f(:,:,n1+2,j)+20*f(:,:,n1+3,j))
!         f(:,:,n1-3,j)=         - 7*f(:,:,n1+1,j)- 4*f(:,:,n1+2,j)+ 4*f(:,:,n1+3,j)

!       case('top')               ! top boundary
!         f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
!         f(:,:,n2+1,j)=(1/11.)*(-13*f(:,:,n2-1,j)-14*f(:,:,n2-2,j)+10*f(:,:,n2-3,j))
!         f(:,:,n2+2,j)=(1/11.)*(-48*f(:,:,n2-1,j)-17*f(:,:,n2-2,j)+20*f(:,:,n2-3,j))
!         f(:,:,n2+3,j)=         - 7*f(:,:,n2-1,j)- 4*f(:,:,n2-2,j)+ 4*f(:,:,n2-3,j)

      case default
        print*, "bc_extrap0_2_0: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_extrap0_2_0
!***********************************************************************
    subroutine bc_extrap0_2_1(f,topbot,j)
!
!  Extrapolation boundary condition for f(bdry)=0.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  NOTE: This is not the final formula, but just bc_extrap_2_1() with f(bdry)=0
!
!    9-oct-03/wolf: coded
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j
!
      select case(topbot)

      case('bot')               ! bottom boundary
        f(:,:,n1  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=0.25*(- 3*f(:,:,n1+1,j)- 5*f(:,:,n1+2,j)+ 3*f(:,:,n1+3,j))
        f(:,:,n1-2,j)=0.05*(-43*f(:,:,n1+1,j)-57*f(:,:,n1+2,j)+39*f(:,:,n1+3,j))
        f(:,:,n1-3,j)=0.05*(-81*f(:,:,n1+1,j)-99*f(:,:,n1+2,j)+73*f(:,:,n1+3,j))

      case('top')               ! top boundary
        f(:,:,n2  ,j)=0.        ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=0.25*(- 3*f(:,:,n2-1,j)- 5*f(:,:,n2-2,j)+ 3*f(:,:,n2-3,j))
        f(:,:,n2+2,j)=0.05*(-43*f(:,:,n2-1,j)-57*f(:,:,n2-2,j)+39*f(:,:,n2-3,j))
        f(:,:,n2+3,j)=0.05*(-81*f(:,:,n2-1,j)-99*f(:,:,n2-2,j)+73*f(:,:,n2-3,j))

      case default
        print*, "bc_extrap0_2_1: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_extrap0_2_1
!***********************************************************************
    subroutine bc_extrap0_2_2(f,topbot,j)
!
!  Extrapolation boundary condition for f(bdry)=0.
!  Correct for polynomials up to 2nd order, determined 1 further degree
!  of freedom by minimizing L2 norm of coefficient vector.
!
!  NOTE: This is not the final formula, but just bc_extrap_2_2() with f(bdry)=0
!
!    9-oct-03/wolf: coded
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
        f(:,:,n1  ,j)= 0.       ! set bdry value=0 (indep of initcond)
        f(:,:,n1-1,j)=0.2   *(                 -  4*f(:,:,n1+2,j)- 3*f(:,:,n1+3,j)+ 3*f(:,:,n1p4,j))
        f(:,:,n1-2,j)=0.2   *(- 2*f(:,:,n1+1,j)-  9*f(:,:,n1+2,j)- 6*f(:,:,n1+3,j)+ 7*f(:,:,n1p4,j))
        f(:,:,n1-3,j)=1./35.*(-33*f(:,:,n1+1,j)-108*f(:,:,n1+2,j)-68*f(:,:,n1+3,j)+87*f(:,:,n1p4,j))

      case('top')               ! top boundary
        f(:,:,n2  ,j)= 0.       ! set bdry value=0 (indep of initcond)
        f(:,:,n2+1,j)=0.2   *(                 -  4*f(:,:,n2-2,j)- 3*f(:,:,n2-3,j)+ 3*f(:,:,n2m4,j))
        f(:,:,n2+2,j)=0.2   *(- 2*f(:,:,n2-1,j)-  9*f(:,:,n2-2,j)- 6*f(:,:,n2-3,j)+ 7*f(:,:,n2m4,j))
        f(:,:,n2+3,j)=1./35.*(-33*f(:,:,n2-1,j)-108*f(:,:,n2-2,j)-68*f(:,:,n2-3,j)+87*f(:,:,n2m4,j))

      case default
        print*, "bc_extrap0_2_2: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_extrap0_2_2
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
        print*,"bc_db_z: invalid argument for 'bc_db_z'"
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
        print*,"bc_db_x: invalid argument for 'bc_db_x'"
      endselect
!
    endsubroutine bc_db_x
!***********************************************************************
    subroutine bc_force_z(f,topbot,j)
!
!  Force values of j-th variable on vertical boundary topbot.
!  This can either be used for freezing variables at the boundary, or for
!  enforcing a certain time-dependent function of (x,y).
!
!  Currently this is hard-coded for velocity components (ux,uy) and quite
!  useless. Plan is to read time-dependent velocity field from disc and
!  apply it as boundary condition here.
!
!  26-apr-2004/wolf: coded
!
      use Cdata
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: i,j
!
      if (topbot /= 'bot') &
           call stop_it("BC_FORCE_Z: only implemented for lower boundary yet")
!
      select case (force_lower_bound)
      case ('uxy_sin-cos')
        call bc_force_uxy_sin_cos(f,n1,j)
      case ('uxy_convection')
        call stop_it("BC_FORCE_Z: uxy_convection not implemented yet")
      case default
        if (lroot) print*, "No such value for force_lower_bound: <", &
             trim(force_lower_bound),">"
        call stop_it("")
      endselect
!
!  Now fill ghost zones imposing antisymmetry w.r.t. the values just set:
!
      do i=1,nghost; f(:,:,n1-i,j)=2*f(:,:,n1,j)-f(:,:,n1+i,j); enddo
!
    endsubroutine bc_force_z
!***********************************************************************
    subroutine bc_force_uxy_sin_cos(f,iiz,j)
!
!  Set (ux, uy) = (cos y, sin x) in vertical layer
!
!  26-apr-2004/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: iiz,j
      real :: kx,ky
!
      if (iuz == 0) call stop_it("BC_FORCE_UXY_SIN_COS: Bad idea...")
!
      if (j==iux) then
        if (Ly>0) then; ky=2*pi/Ly; else; ky=0.; endif 
        f(:,:,iiz,j) = spread(cos(ky*y),1,mx)
      elseif (j==iuy) then
        if (Lx>0) then; kx=2*pi/Lx; else; kx=0.; endif 
        f(:,:,iiz,j) = spread(sin(kx*x),2,my)
      elseif (j==iuz) then
        f(:,:,n1,j) = 0.
      endif
!
    endsubroutine bc_force_uxy_sin_cos
!***********************************************************************
    subroutine bc_one_x(f,topbot,j)
!
!  Set bdry values to 1 for debugging purposes
!
!  11-jul-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case(topbot)

      case('bot')               ! bottom boundary
          f(1:l1-1,:,:,j)=1.

      case('top')               ! top boundary
          f(l2+1:mx,:,:,j)=1.

      case default
        print*, "bc_one_x: ",topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_one_x
!***********************************************************************
    subroutine bc_one_y(f,topbot,j)
!
!  Set bdry values to 1 for debugging purposes
!
!  11-jul-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case(topbot)

      case('bot')               ! bottom boundary
          f(:,1:m1-1,:,j)=1.

      case('top')               ! top boundary
          f(:,m2+1:my,:,j)=1.

      case default
        print*, "bc_one_y: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_one_y
!***********************************************************************
    subroutine bc_one_z(f,topbot,j)
!
!  Set bdry values to 1 for debugging purposes
!
!  11-jul-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      integer :: j
      character (len=3) :: topbot
!
      select case(topbot)

      case('bot')               ! bottom boundary
          f(:,:,1:n1-1,j)=1.

      case('top')               ! top boundary
          f(:,:,n2+1:mz,j)=1.

      case default
        print*, "bc_one_z: ", topbot, " should be `top' or `bot'"

      endselect
!
    endsubroutine bc_one_z
!***********************************************************************
    subroutine bc_freeze_var_z(topbot,j)
!
!  Tell other modules that variable with slot j is to be frozen in on
!  given boundary
!
      use Cdata
!
      integer :: j
      character (len=3) :: topbot
!
      lfrozen_bcs_z = .true.    ! set flag

      select case(topbot)
      case('bot')               ! bottom boundary
        lfrozen_bot_var_z(j) = .true.
      case('top')               ! top boundary
        lfrozen_top_var_z(j) = .true.
      case default
        print*, "bc_freeze_var_z: ", topbot, " should be `top' or `bot'"
      endselect
!
    endsubroutine bc_freeze_var_z
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
