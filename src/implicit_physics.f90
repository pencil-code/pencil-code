! $Id$
!
! This module solves the radiative diffusion implicitly thanks
! to an Alternate Direction Implicit Scheme (ADI) in a D'Yakonov
! form
!     lambda_x T(n+1/2) = lambda_x + lambda_z
!     lambda_z T(n+1) = T(n+1/2)
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
! COMMUNICATED AUXILIARIES 1
!
!***************************************************************
module ImplicitPhysics
!
  use Cdata
  use Cparam
  use Messages, only: svn_id, fatal_error
  use Sub, only: keep_compiler_quiet
  use General, only: tridag, cyclic
!
  implicit none
!
  include 'implicit_physics.h'
!
  interface heatcond_TT ! Overload subroutine `hcond_TT' function
    module procedure heatcond_TT_0d  ! get one value (hcond, dhcond)
    module procedure heatcond_TT_1d  ! get 1d-arrays (hcond, dhcond)
    module procedure heatcond_TT_2d  ! get 2d-arrays (hcond, dhcond)
  end interface
!
  real, pointer :: hcond0, Fbot
  logical, pointer :: lADI_mixed
  real :: Tbump, Kmax, Kmin, hole_slope, hole_width, hole_alpha
  real :: dx_2, dz_2, cp1
  logical :: lyakonov=.true.
!
  contains
!***********************************************************************
    subroutine register_implicit_physics()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrho; increase nvar accordingly.
!
!  03-mar-2010/dintrans: coded
!
      use FArrayManager, only: farray_register_auxiliary
!
      call farray_register_auxiliary('TTold',iTTold,communicated=.true.)
      print*, 'iTTold=', iTTold
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
       "$Id$")
!
    endsubroutine register_implicit_physics
!***********************************************************************
    subroutine initialize_implicit_physics(f)
!
      use SharedVariables, only: get_shared_variable
      use MpiComm, only: stop_it
      use EquationOfState, only: get_cp1
!
      implicit none
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: ierr
      real, dimension(:), pointer :: hole_params
!
      call get_shared_variable('hcond0', hcond0, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                 "there was a problem when getting hcond0")
      call get_shared_variable('Fbot', Fbot, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                "there was a problem when getting Fbot")
      call get_shared_variable('lADI_mixed', lADI_mixed, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                "there was a problem when getting lADI_mixed")
      call get_shared_variable('hole_params', hole_params, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                "there was a problem when getting the hole_params array")
      Tbump=hole_params(1)
      Kmax=hole_params(2)
      Kmin=hole_params(3)
      hole_slope=hole_params(4)
      hole_width=hole_params(5)
      hole_alpha=(Kmax-Kmin)/(pi/2.+atan(hole_slope*hole_width**2))
      if (lroot .and. ldebug) then
        print*, '************ hole parameters ************'
        print*,'Tbump, Kmax, Kmin, hole_slope, hole_width, hole_alpha=', &
               Tbump, Kmax, Kmin, hole_slope, hole_width, hole_alpha
        print*, '*****************************************'
      endif
!
      if (lrun) then
! hcondADI is dynamically shared with boundcond() for the 'c3' BC
        call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
      else
        hcondADI=spread(Kmax, 1, mx)
      endif
!
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
!
    endsubroutine initialize_implicit_physics
!***********************************************************************
    subroutine calc_heatcond_ADI(f)
!
!  10-sep-07/gastine+dintrans: wrapper to the two possible ADI subroutines
!  ADI_Kconst: constant radiative conductivity
!  ADI_Kprof: radiative conductivity depends on T, i.e. hcond(T)
!
      implicit none
!
      real, dimension(mx,my,mz,mfarray) :: f
      integer :: ierr
!
      if (hcond0 /= impossible) then
        if (nx == 1) then
          call ADI_Kconst_1d(f)
        else
          if (nprocz>1) then
            call ADI_Kconst_MPI(f)
          else
            if (lyakonov) then
              call ADI_Kconst_yakonov(f)
            else
              call ADI_Kconst(f)
            endif
          endif
        endif
      else
        if (nx == 1) then
          if (lADI_mixed) then
            call ADI_Kprof_1d_mixed(f)
          else
            call ADI_Kprof_1d(f)
          endif
        else
          if (nprocz>1) then
            if (lADI_mixed) then
              call ADI_Kprof_MPI_mixed(f)
            else
              call ADI_Kprof_MPI(f)
            endif
          else
            if (lADI_mixed) then
              call ADI_Kprof_mixed(f)
            else
              call ADI_Kprof(f)
            endif
          endif
        endif
      endif
!
    endsubroutine calc_heatcond_ADI
!***********************************************************************
    subroutine ADI_Kconst(f)
!
!  08-Sep-07/gastine+dintrans: coded
!  2-D ADI scheme for the radiative diffusion term (see
!  Peaceman & Rachford 1955). Each direction is solved implicitly:
!
!    (1-dt/2*Lambda_x)*T^(n+1/2) = (1+dt/2*Lambda_y)*T^n + source/2
!    (1-dt/2*Lambda_y)*T^(n+1)   = (1+dt/2*Lambda_x)*T^(n+1/2) + source/2
!
!  where Lambda_x and Lambda_y denote diffusion operators and the source
!  term comes from the explicit advance.
!
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, mpoly0
      use Boundcond, only: update_ghosts
!
      implicit none
!
      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz) :: finter, source, rho, TT
      real, dimension(nx)    :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nz)    :: az, bz, cz, wz, rhsz, workz
      real    :: aalpha, bbeta
!
!  first update all the ghost zones in the f-array
!
      call update_ghosts(f)
!
      TT=f(:,4,:,iTTold)
      source=(f(:,4,:,ilnTT)-TT)/dt
      if (ldensity) then
        rho=exp(f(:,4,:,ilnrho))
      else
        rho=1.
      endif
!
!  row dealt implicitly
!
      do j=n1,n2
        wx=dt*gamma*hcond0*cp1/rho(l1:l2,j)
        ax=-wx*dx_2/2.
        bx=1.+wx*dx_2
        cx=ax
        rhsx=TT(l1:l2,j)+wx*dz_2/2.*                         &
             (TT(l1:l2,j+1)-2.*TT(l1:l2,j)+TT(l1:l2,j-1))    &
             +dt/2.*source(l1:l2,j)
!
! x boundary conditions: periodic
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax, bx, cx, aalpha, bbeta, rhsx, workx, nx)
        finter(l1:l2,j)=workx
      enddo
!
! finter must be periodic in the x-direction
!
      finter(1:l1-1,:)=finter(l2i:l2,:)
      finter(l2+1:mx,:)=finter(l1:l1i,:)
!
!  columns dealt implicitly
!
      do i=l1,l2
        wz=dt*gamma*hcond0*cp1/rho(i,n1:n2)
        az=-wz*dz_2/2.
        bz=1.+wz*dz_2
        cz=az
        rhsz=finter(i,n1:n2)+wz*dx_2/2.*                               &
             (finter(i+1,n1:n2)-2.*finter(i,n1:n2)+finter(i-1,n1:n2))  &
             +dt/2.*source(i,n1:n2)
        !
        ! z boundary conditions
        ! Always constant temperature at the top
        !
        bz(nz)=1. ; az(nz)=0.
        rhsz(nz)=cs2top/gamma_m1
        select case (bcz1(ilnTT))
          ! Constant temperature at the bottom
          case ('cT')
            bz(1)=1.  ; cz(1)=0. 
            rhsz(1)=cs2bot/gamma_m1
          ! Constant flux at the bottom
          case ('c1')
            bz(1)=1.   ; cz(1)=-1
            rhsz(1)=dz*Fbot/hcond0
! we can use here the second-order relation for the first derivative: 
! (T_{j+1}-T_{j_1})/2dz = dT/dz --> T_{j-1} = T_{j+1} - 2*dz*dT/dz 
! and insert this expression in the difference relation to eliminate T_{j-1}:
! a_{j-1}*T_{j-1} + b_j T_j + c_{j+1}*T_{j+1} = RHS
!           cz(1)=cz(1)+az(1)
!           rhsz(1)=rhsz(1)-2.*az(1)*dz*Fbot/hcond0
          case default 
           call fatal_error('ADI_Kconst','bcz on TT must be cT or c1')
        endselect
!
        call tridag(az, bz, cz, rhsz, workz)
        f(i,4,n1:n2,ilnTT)=workz
      enddo
!
    endsubroutine ADI_Kconst
!***********************************************************************
    subroutine ADI_Kprof(f)
!
!  10-Sep-07/gastine+dintrans: coded
!  2-D ADI scheme for the radiative diffusion term where the radiative
!  conductivity depends on T (uses heatcond_TT to compute hcond _and_
!  dhcond). The ADI scheme is of Yakonov's form:
!
!    (1-dt/2*J_x)*lambda = f_x(T^n) + f_y(T^n) + source
!    (1-dt/2*J_y)*beta   = lambda
!    T^(n+1) = T^n + dt*beta
!
!    where J_x and J_y denote Jacobian matrices df/dT.
!
      use EquationOfState, only: gamma
!
      implicit none
!
      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz) :: source, hcond, dhcond, finter, val, TT, rho
      real, dimension(nx)    :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nz)    :: az, bz, cz, wz, rhsz, workz
      real    :: aalpha, bbeta
!
      source=(f(:,4,:,ilnTT)-f(:,4,:,iTTold))/dt
! BC important not for the x-direction (always periodic) but for 
! the z-direction as we must impose the 'c3' BC at the 2nd-order
! before going in the implicit stuff
      call heatcond_TT(f(:,4,:,iTTold), hcond, dhcond)
      call boundary_ADI(f(:,4,:,iTTold), hcond(:,n1))
      TT=f(:,4,:,iTTold)
      if (ldensity) then
        rho=exp(f(:,4,:,ilnrho))
      else
        rho=1.
      endif
!
!  rows dealt implicitly
!
      do j=n1,n2
       wx=cp1*gamma/rho(l1:l2,j)
! ax=-dt/2*J_x for i=i-1 (lower diagonal)
       ax=-dt*wx*dx_2/4.*(dhcond(l1-1:l2-1,j)    &
         *(TT(l1-1:l2-1,j)-TT(l1:l2,j))          &
         +hcond(l1-1:l2-1,j)+hcond(l1:l2,j))
! bx=1-dt/2*J_x for i=i (main diagonal)
       bx=1.+dt*wx*dx_2/4.*(dhcond(l1:l2,j)      &
         *(2.*TT(l1:l2,j)-TT(l1-1:l2-1,j)        &
         -TT(l1+1:l2+1,j))+2.*hcond(l1:l2,j)     &
         +hcond(l1+1:l2+1,j)+hcond(l1-1:l2-1,j))
! cx=-dt/2*J_x for i=i+1 (upper diagonal)
       cx=-dt*wx*dx_2/4.*(dhcond(l1+1:l2+1,j)    &
          *(TT(l1+1:l2+1,j)-TT(l1:l2,j))         &
          +hcond(l1:l2,j)+hcond(l1+1:l2+1,j))
! rhsx=f_y(T^n) + f_x(T^n) (Eq. 3.6)
! do first f_y(T^n)
       rhsx=wx*dz_2/2.*((hcond(l1:l2,j+1)        &
           +hcond(l1:l2,j))*(TT(l1:l2,j+1)       &
           -TT(l1:l2,j))-(hcond(l1:l2,j)         &
           +hcond(l1:l2,j-1))                    &
           *(TT(l1:l2,j)-TT(l1:l2,j-1)))
! then add f_x(T^n)
       rhsx=rhsx+wx*dx_2/2.*((hcond(l1+1:l2+1,j)         &
         +hcond(l1:l2,j))*(TT(l1+1:l2+1,j)-TT(l1:l2,j))  &
           -(hcond(l1:l2,j)+hcond(l1-1:l2-1,j))          &
           *(TT(l1:l2,j)-TT(l1-1:l2-1,j)))+source(l1:l2,j)
!
! x boundary conditions: periodic
       aalpha=cx(nx) ; bbeta=ax(1)
       call cyclic(ax,bx,cx,aalpha,bbeta,rhsx,workx,nx)
       finter(l1:l2,j)=workx(1:nx)
      enddo
!
!  columns dealt implicitly
!
      do i=l1,l2
       wz=dt*cp1*gamma*dz_2/rho(i,n1:n2)
       az=-wz/4.*(dhcond(i,n1-1:n2-1)   &
         *(TT(i,n1-1:n2-1)-TT(i,n1:n2)) &
         +hcond(i,n1-1:n2-1)+hcond(i,n1:n2))
!
       bz=1.+wz/4.*(dhcond(i,n1:n2)*             &
         (2.*TT(i,n1:n2)-TT(i,n1-1:n2-1)         &
         -TT(i,n1+1:n2+1))+2.*hcond(i,n1:n2)     &
         +hcond(i,n1+1:n2+1)+hcond(i,n1-1:n2-1))
!
       cz=-wz/4.*(dhcond(i,n1+1:n2+1)            &
         *(TT(i,n1+1:n2+1)-TT(i,n1:n2))          &
         +hcond(i,n1:n2)+hcond(i,n1+1:n2+1))
!
       rhsz=finter(i,n1:n2)
!
! z boundary conditions
! Constant temperature at the top: T^(n+1)-T^n=0
       bz(nz)=1. ; az(nz)=0.
       rhsz(nz)=0.
! bottom
       select case (bcz1(ilnTT))
! Constant temperature at the bottom: T^(n+1)-T^n=0
         case ('cT')
          bz(1)=1. ; cz(1)=0.
          rhsz(1)=0.
! Constant flux at the bottom
         case ('c3')
          bz(1)=1. ; cz(1)=-1.
          rhsz(1)=0.
         case default 
          call fatal_error('ADI_Kprof','bcz on TT must be cT or c3')
       endselect
!
       call tridag(az,bz,cz,rhsz,workz)
       val(i,n1:n2)=workz(1:nz)
      enddo
!
      f(:,4,:,ilnTT)=f(:,4,:,iTTold)+dt*val
!
! update hcond used for the 'c3' condition in boundcond.f90
!
      call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof
!***********************************************************************
    subroutine ADI_Kprof_MPI(f)
!
!  15-jan-10/gastine: coded
!  2-D ADI scheme for the radiative diffusion term where the radiative
!  conductivity depends on T (uses heatcond_TT to compute hcond _and_
!  dhcond). The ADI scheme is of Yakonov's form:
!
!    (1-dt/2*J_x)*lambda = f_x(T^n) + f_z(T^n) + source
!    (1-dt/2*J_z)*beta   = lambda
!    T^(n+1) = T^n + dt*beta
!
!    where J_x and J_z denote Jacobian matrices df/dT.
!  08-mar-2010/dintrans: added the case of a non-square domain (ibox-loop)
!  21-aug-2010/dintrans: simplified version that uses Anders' original
!    transp_xz and transp_zx subroutines
!
      use EquationOfState, only: gamma
      use Mpicomm, only: transp_xz, transp_zx
      use Boundcond, only: update_ghosts
!
      implicit none
!
      integer, parameter :: mzt=nzgrid+2*nghost
      integer, parameter :: n1t=nghost+1, n2t=n1t+nzgrid-1
      integer, parameter :: nxt=nx/nprocz
      integer :: i ,j
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz)   :: source, hcond, dhcond, finter, TT, rho
      real, dimension(mzt,nxt) :: hcondt, dhcondt, fintert, TTt, rhot, valt
      real, dimension(nx,nz)   :: val
      real, dimension(nx)      :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nzgrid)  :: az, bz, cz, wz, rhsz, workz
      real :: aalpha, bbeta
!
!  It is necessary to communicate ghost-zones points between
!  processors to ensure a correct transposition of these ghost
!  zones. It is needed by rho,rhot and source,sourcet.
!
      call update_ghosts(f)
      source=(f(:,4,:,ilnTT)-f(:,4,:,iTTold))/dt
!
! BC important not for the x-direction (always periodic) but for 
! the z-direction as we must impose the 'c3' BC at the 2nd-order
! before going in the implicit stuff
!
      TT=f(:,4,:,iTTold)
      call heatcond_TT(TT, hcond, dhcond)
      call boundary_ADI(TT, hcond(:,n1))
      if (ldensity) then
        rho=exp(f(:,4,:,ilnrho))
      else
        rho=1.
      endif
!
! rows dealt implicitly
!
      do j=n1,n2
        wx=cp1*gamma/rho(l1:l2,j)
! ax=-dt/2*J_x for i=i-1 (lower diagonal)
        ax=-dt*wx*dx_2/4.*(dhcond(l1-1:l2-1,j)     &
           *(TT(l1-1:l2-1,j)-TT(l1:l2,j))          &
           +hcond(l1-1:l2-1,j)+hcond(l1:l2,j))
! bx=1-dt/2*J_x for i=i (main diagonal)
        bx=1.+dt*wx*dx_2/4.*(dhcond(l1:l2,j)       &
           *(2.*TT(l1:l2,j)-TT(l1-1:l2-1,j)        &
           -TT(l1+1:l2+1,j))+2.*hcond(l1:l2,j)     &
           +hcond(l1+1:l2+1,j)+hcond(l1-1:l2-1,j))
! cx=-dt/2*J_x for i=i+1 (upper diagonal)
        cx=-dt*wx*dx_2/4.*(dhcond(l1+1:l2+1,j)     &
           *(TT(l1+1:l2+1,j)-TT(l1:l2,j))          &
           +hcond(l1:l2,j)+hcond(l1+1:l2+1,j))
! rhsx=f_z(T^n) + f_x(T^n) (Eq. 3.6)
! do first f_z(T^n)
        rhsx=wx*dz_2/2.*((hcond(l1:l2,j+1)         &
             +hcond(l1:l2,j))*(TT(l1:l2,j+1)       &
             -TT(l1:l2,j))-(hcond(l1:l2,j)         &
             +hcond(l1:l2,j-1))                    &
             *(TT(l1:l2,j)-TT(l1:l2,j-1)))
! then add f_x(T^n)
        rhsx=rhsx+wx*dx_2/2.*((hcond(l1+1:l2+1,j)            &
             +hcond(l1:l2,j))*(TT(l1+1:l2+1,j)-TT(l1:l2,j))  &
             -(hcond(l1:l2,j)+hcond(l1-1:l2-1,j))            &
             *(TT(l1:l2,j)-TT(l1-1:l2-1,j)))+source(l1:l2,j)
!
! periodic boundary conditions in the x-direction
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax, bx, cx, aalpha, bbeta, rhsx, workx, nx)
        finter(l1:l2,j)=workx(1:nx)
      enddo
!
! do the transpositions x <--> z
!
      call transp_xz(finter(l1:l2,n1:n2), fintert(n1t:n2t,:))
      call transp_xz(rho(l1:l2,n1:n2), rhot(n1t:n2t,:))
      call transp_xz(TT(l1:l2,n1:n2), TTt(n1t:n2t,:))
      call heatcond_TT(TTt, hcondt, dhcondt)
!
      do i=1,nxt
        wz=dt*cp1*gamma*dz_2/rhot(n1t:n2t,i)
        az=-wz/4.*(dhcondt(n1t-1:n2t-1,i)            &
           *(TTt(n1t-1:n2t-1,i)-TTt(n1t:n2t,i))      &
           +hcondt(n1t-1:n2t-1,i)+hcondt(n1t:n2t,i))
!
        bz=1.+wz/4.*(dhcondt(n1t:n2t,i)*                 &
           (2.*TTt(n1t:n2t,i)-TTt(n1t-1:n2t-1,i)         &
           -TTt(n1t+1:n2t+1,i))+2.*hcondt(n1t:n2t,i)     &
           +hcondt(n1t+1:n2t+1,i)+hcondt(n1t-1:n2t-1,i))
!
        cz=-wz/4.*(dhcondt(n1t+1:n2t+1,i)            &
           *(TTt(n1t+1:n2t+1,i)-TTt(n1t:n2t,i))      &
           +hcondt(n1t:n2t,i)+hcondt(n1t+1:n2t+1,i))
!
        rhsz=fintert(n1t:n2t,i)
!
! z boundary conditions
! Constant temperature at the top: T^(n+1)-T^n=0
!
        bz(nzgrid)=1. ; az(nzgrid)=0.
        rhsz(nzgrid)=0.
! bottom
        select case (bcz1(ilnTT))
! Constant temperature at the bottom: T^(n+1)-T^n=0
          case ('cT')
            bz(1)=1. ; cz(1)=0.
            rhsz(1)=0.
! Constant flux at the bottom
          case ('c3')
            bz(1)=1. ; cz(1)=-1.
            rhsz(1)=0.
          case default 
            call fatal_error('ADI_Kprof','bcz on TT must be cT or c3')
        endselect
        call tridag(az, bz, cz, rhsz, workz)
        valt(n1t:n2t,i)=workz(1:nzgrid)
      enddo ! i
!
! come back on the grid (x,z)
!
      call transp_zx(valt(n1t:n2t,:), val)
      f(l1:l2,4,n1:n2,ilnTT)=f(l1:l2,4,n1:n2,iTTold)+dt*val
!
! update hcond used for the 'c3' condition in boundcond.f90
!
      if (iproc==0) call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof_MPI
!***********************************************************************
    subroutine boundary_ADI(f_2d, hcond)
!
! 13-Sep-07/gastine: computed two different types of boundary 
! conditions for the implicit solver:
!     - Always periodic in x-direction
!     - Possibility to choose between 'cT' and 'c3' in z direction
! Note: 'c3' means that the flux is constant at the *bottom only*
!
      implicit none
!
      real, dimension(mx,mz) :: f_2d
      real, dimension(mx), optional :: hcond
!
! x-direction: periodic
!
      f_2d(1:l1-1,:)=f_2d(l2i:l2,:)
      f_2d(l2+1:mx,:)=f_2d(l1:l1i,:)
!
! top boundary condition z=z(n2): always constant temperature
!
      if (llast_proc_z) then
        f_2d(:,n2+1)=2.*f_2d(:,n2)-f_2d(:,n2-1)
      endif
!
! bottom bondary condition z=z(n1): constant T or imposed flux dT/dz
!
      if (iproc==0) then
      select case (bcz1(ilnTT))
        case ('cT') ! constant temperature
          f_2d(:,n1-1)=2.*f_2d(:,n1)-f_2d(:,n1+1)
        case ('c3') ! constant flux
          if (.not. present(hcond)) then
            f_2d(:,n1-1)=f_2d(:,n1+1)+2.*dz*Fbot/hcond0
          else 
            f_2d(:,n1-1)=f_2d(:,n1+1)+2.*dz*Fbot/hcond(:)
          endif
      endselect
      endif
!
    endsubroutine boundary_ADI
!***********************************************************************
    subroutine ADI_Kconst_1d(f)
!
! 18-sep-07/dintrans: coded
! Implicit Crank Nicolson scheme in 1-D for a constant K (not 
! really an ADI but keep the generic name for commodity).
!
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top
!
      implicit none
!
      integer :: j, jj
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mz) :: rho, TT
      real, dimension(nz) :: a, b, c, rhs, work
      real :: wz
!
      TT=f(4,4,:,ilnTT)
      rho=exp(f(4,4,:,ilnrho))
!
      do j=n1,n2
        wz=dt*gamma*hcond0*cp1/rho(j)
        jj=j-nghost
        a(jj)=-wz*dz_2/2.
        b(jj)=1.+wz*dz_2
        c(jj)=a(jj)
!
!       rhs(jj)=TT(j)+wz*dz_2/2.*(TT(j+1)-2.*TT(j)+TT(j-1))+dt*source(j)
        rhs(jj)=TT(j)+wz*dz_2/2.*(TT(j+1)-2.*TT(j)+TT(j-1))
      enddo
! apply the boundary conditions *outside* the j-loop
! Always constant temperature at the top
      b(nz)=1. ; a(nz)=0.
      rhs(nz)=cs2top/gamma_m1
      if (bcz1(ilnTT)=='cT') then
! Constant temperature at the bottom
        b(1)=1. ; c(1)=0. 
        rhs(1)=cs2bot/gamma_m1
      else
! Constant flux at the bottom
        b(1)=1.  ; c(1)=-1.
        rhs(1)=dz*Fbot/hcond0
      endif
      call tridag(a, b, c, rhs, work)
      f(4,4,n1:n2,ilnTT)=work
!
    endsubroutine ADI_Kconst_1d
!***********************************************************************
    subroutine ADI_Kprof_1d(f)
!
! 18-sep-07/dintrans: coded
! Implicit 1-D case for a temperature-dependent conductivity K(T).
! Not really an ADI but keep the generic name for commodity.
!
      use EquationOfState, only: gamma
!
      implicit none
!
      integer :: j, jj
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mz) :: source, rho, TT, hcond, dhcond
      real, dimension(nz) :: a, b, c, rhs, work
      real :: wz, hcondp, hcondm
!
      source=(f(4,4,:,ilnTT)-f(4,4,:,iTTold))/dt
      rho=exp(f(4,4,:,ilnrho))
!
! need to set up the 'c3' BC at the 2nd-order before the implicit stuff
!
      call heatcond_TT(f(4,4,:,iTTold), hcond, dhcond)
      hcondADI=spread(hcond(1), 1, mx)
      call boundary_ADI(f(:,4,:,iTTold), hcondADI)
      TT=f(4,4,:,iTTold)
!
      do j=n1,n2
        jj=j-nghost
        wz=dt*dz_2*gamma*cp1/rho(j)
        hcondp=hcond(j+1)+hcond(j)
        hcondm=hcond(j)+hcond(j-1)
!
        a(jj)=-wz/4.*(hcondm-dhcond(j-1)*(TT(j)-TT(j-1)))
        b(jj)=1.-wz/4.*(-hcondp-hcondm+dhcond(j)*(TT(j+1)-2.*TT(j)+TT(j-1)))
        c(jj)=-wz/4.*(hcondp+dhcond(j+1)*(TT(j+1)-TT(j)))
        rhs(jj)=wz/2.*(hcondp*(TT(j+1)-TT(j))-hcondm*(TT(j)-TT(j-1))) &
                +dt*source(j)
!
! Always constant temperature at the top: T^(n+1)-T^n = 0
!
        b(nz)=1. ; a(nz)=0.
        rhs(nz)=0.
        if (bcz1(ilnTT)=='cT') then
! Constant temperature at the bottom
          b(1)=1. ; c(1)=0. 
          rhs(1)=0.
        else
! Constant flux at the bottom: d/dz [T^(n+1)-T^n] = 0
          b(1)=1.  ; c(1)=-1.
          rhs(1)=0.
        endif
      enddo
      call tridag(a, b, c, rhs, work)
      f(4,4,n1:n2,ilnTT)=f(4,4,n1:n2,iTTold)+work
!
! Update the bottom value of hcond used for the 'c3' BC in boundcond
!
      call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof_1d
!***********************************************************************
    subroutine ADI_Kconst_MPI(f)
!
!  04-sep-2009/dintrans: coded
!  parallel version of the ADI scheme for the K=cte case
!
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top
      use Mpicomm, only: transp_xz, transp_zx
      use Boundcond, only: update_ghosts
!
      implicit none
!
      integer, parameter :: nxt=nx/nprocz
      integer, parameter :: n1t=nghost+1, n2t=n1t+nzgrid-1
      integer, parameter :: mzt=nzgrid+2*nghost
      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz)   :: finter, source, rho, TT
      real, dimension(mzt,nxt) :: fintert, rhot, sourcet, wtmp
      real, dimension(nx)      :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nzgrid)  :: az, bz, cz, wz, rhsz, workz
      real  :: aalpha, bbeta
!
      call update_ghosts(f)
!
      TT=f(:,4,:,iTTold)
      source=(f(:,4,:,ilnTT)-TT)/dt
      if (ldensity) then
        rho=exp(f(:,4,:,ilnrho))
      else
        rho=1.
      endif
!
! Rows dealt implicitly
!
      do j=n1,n2
        wx=dt*gamma*hcond0*cp1/rho(l1:l2,j)
        ax=-wx*dx_2/2
        bx=1.+wx*dx_2
        cx=ax
        rhsx=TT(l1:l2,j)+wx*dz_2/2*                      &
             (TT(l1:l2,j+1)-2*TT(l1:l2,j)+TT(l1:l2,j-1)) &
             +dt/2*source(l1:l2,j)
!
! x boundary conditions: periodic
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax, bx, cx, aalpha, bbeta, rhsx, workx, nx)
        finter(l1:l2,j)=workx(1:nx)
      enddo
!
! Do the transpositions x <--> z
!
       call transp_xz(finter(l1:l2,n1:n2), fintert(n1t:n2t,:))
       call transp_xz(rho(l1:l2,n1:n2), rhot(n1t:n2t,:))
       call transp_xz(source(l1:l2,n1:n2), sourcet(n1t:n2t,:))
!
! Columns dealt implicitly
!
      do i=1,nxt
        wz=dt*gamma*hcond0*cp1/rhot(n1t:n2t,i)
        az=-wz*dz_2/2
        bz=1.+wz*dz_2
        cz=az
        rhsz=fintert(n1t:n2t,i)+wz*dx_2/2*                        &
          (fintert(n1t:n2t,i+1)-2*fintert(n1t:n2t,i)+fintert(n1t:n2t,i-1))  &
          +dt/2*sourcet(n1t:n2t,i)
        !
        ! z boundary conditions
        ! Always constant temperature at the top
        !
        bz(nzgrid)=1. ; az(nzgrid)=0.
        rhsz(nzgrid)=cs2top/gamma_m1
        select case (bcz1(ilnTT))
          ! Constant temperature at the bottom
          case ('cT')
            bz(1)=1.  ; cz(1)=0. 
            rhsz(1)=cs2bot/gamma_m1
          ! Constant flux at the bottom
          case ('c1')
            bz(1)=1.   ; cz(1)=-1
            rhsz(1)=dz*Fbot/hcond0
! we can use here the second-order relation for the first derivative: 
! (T_{j+1}-T_{j_1})/2dz = dT/dz --> T_{j-1} = T_{j+1} - 2*dz*dT/dz 
! and insert this expression in the difference relation to eliminate T_{j-1}:
! a_{j-1}*T_{j-1} + b_j T_j + c_{j+1}*T_{j+1} = RHS
!
!           cz(1)=cz(1)+az(1)
!           rhsz(1)=rhsz(1)-2.*az(1)*dz*Fbot/hcond0
          case default 
           call fatal_error('ADI_Kconst','bcz on TT must be cT or c1')
        endselect
!
        call tridag(az, bz, cz, rhsz, workz)
        wtmp(n1t:n2t,i)=workz(1:nzgrid)
      enddo
      call transp_zx(wtmp(n1t:n2t,:), f(l1:l2,4,n1:n2,ilnTT))
!
    endsubroutine ADI_Kconst_MPI
!***********************************************************************
    subroutine heatcond_TT_2d(TT, hcond, dhcond)
!
! 07-Sep-07/gastine: computed 2-D radiative conductivity hcond(T) with
! its derivative dhcond=dhcond(T)/dT.
!
      implicit none
!
      real, dimension(:,:), intent(in) :: TT
      real, dimension(:,:), intent(out) :: hcond
      real, dimension(:,:), optional :: dhcond
!
      hcond=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      if (present(dhcond)) &
        dhcond=2.*hole_alpha/(1.+hcond**2)*hole_slope*(TT-Tbump)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(hcond))
!
    endsubroutine heatcond_TT_2d
!***********************************************************************
    subroutine heatcond_TT_1d(TT, hcond, dhcond)
!
! 18-Sep-07/dintrans: computed 1-D radiative conductivity 
! hcond(T) with its derivative dhcond=dhcond(T)/dT.
!
      implicit none
!
      real, dimension(:), intent(in) :: TT
      real, dimension(:), intent(out) :: hcond
      real, dimension(:), optional :: dhcond
!
      hcond=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      if (present(dhcond)) &
        dhcond=2.*hole_alpha/(1.+hcond**2)*hole_slope*(TT-Tbump)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(hcond))
!
    endsubroutine heatcond_TT_1d
!***********************************************************************
    subroutine heatcond_TT_0d(TT, hcond, dhcond)
!
! 07-Sep-07/gastine: computed the radiative conductivity hcond(T)
! with its derivative dhcond=dhcond(T)/dT at a given temperature.
!
      implicit none
!
      real, intent(in) :: TT
      real, intent(out) :: hcond
      real, optional :: dhcond
!
      hcond=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      if (present(dhcond)) &
        dhcond=2.*hole_alpha/(1.+hcond**2)*hole_slope*(TT-Tbump)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(hcond))
!
    endsubroutine heatcond_TT_0d
!***********************************************************************
    subroutine ADI_Kprof_1d_mixed(f)
!
! 28-feb-10/dintrans: coded
! Simpler version where a part of the radiative diffusion term is
! computed during the explicit advance.
!
      use EquationOfState, only: gamma
!
      implicit none
!
      integer :: j, jj
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mz) :: source, TT, hcond, dhcond, dLnhcond, chi
      real, dimension(nz) :: a, b, c, rhs, work
      real :: wz
!
      source=(f(4,4,:,ilnTT)-f(4,4,:,iTTold))/dt
      call heatcond_TT(f(4,4,:,iTTold), hcond, dhcond)
!
! need to set up the 'c3' BC at the 2nd-order before the implicit stuff
!
      hcondADI=spread(hcond(1), 1, mx)
      call boundary_ADI(f(:,4,:,iTTold), hcondADI)
      TT=f(4,4,:,iTTold)
      if (ldensity) then
        chi=cp1*hcond/exp(f(4,4,:,ilnrho))
      else
        chi=cp1*hcond
      endif
      dLnhcond=dhcond/hcond
!
      do j=n1,n2
        jj=j-nghost
        wz=dt*dz_2*gamma*chi(j)
!
        a(jj)=-wz/2.
        b(jj)=1.-wz/2.*(-2.+dLnhcond(j)*(TT(j+1)-2.*TT(j)+TT(j-1)))
        c(jj)=-wz/2.
        rhs(jj)=wz*(TT(j+1)-2.*TT(j)+TT(j-1))+dt*source(j)
!
! Always constant temperature at the top: T^(n+1)-T^n = 0
!
        b(nz)=1. ; a(nz)=0.
        rhs(nz)=0.
        if (bcz1(ilnTT)=='cT') then
! Constant temperature at the bottom
          b(1)=1. ; c(1)=0. 
          rhs(1)=0.
        else
! Constant flux at the bottom: d/dz [T^(n+1)-T^n] = 0
          b(1)=1.  ; c(1)=-1.
          rhs(1)=0.
        endif
      enddo
!
      call tridag(a, b, c, rhs, work)
      f(4,4,n1:n2,ilnTT)=f(4,4,n1:n2,iTTold)+work
!
! Update the bottom value of hcond used for the 'c3' BC in boundcond
!
      call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof_1d_mixed
!***********************************************************************
    subroutine ADI_Kprof_mixed(f)
!
!  28-fev-2010/dintrans: coded
!  simpler version where one part of the radiative diffusion term is
!  computed during the explicit step. The implicit part remains 
!  of Yakonov's form:
!
!    (1-dt/2*J_x)*lambda = f_x(T^n) + f_z(T^n) + source
!    (1-dt/2*J_y)*beta   = lambda
!    T^(n+1) = T^n + dt*beta
!
!    where J_x and J_y denote Jacobian matrices df/dT.
!
      use EquationOfState, only: gamma
      use Boundcond, only: update_ghosts
!
      implicit none
!
      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz) :: source, hcond, dhcond, finter, val, TT, &
                                rho, chi, dLnhcond
      real, dimension(nx)    :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nz)    :: az, bz, cz, wz, rhsz, workz
      real :: aalpha, bbeta
!
      call update_ghosts(f)
!
      source=(f(:,4,:,ilnTT)-f(:,4,:,iTTold))/dt
! BC important not for the x-direction (always periodic) but for 
! the z-direction as we must impose the 'c3' BC at the 2nd-order
! before going in the implicit stuff
      call heatcond_TT(f(:,4,:,iTTold), hcond, dhcond)
      call boundary_ADI(f(:,4,:,iTTold), hcond(:,n1))
      TT=f(:,4,:,iTTold)
      if (ldensity) then
        chi=cp1*hcond/exp(f(:,4,:,ilnrho))
!        chi=cp1*hcond0/exp(f(:,4,:,ilnrho))
      else
        chi=cp1*hcond
      endif
      dLnhcond=dhcond/hcond
!      dLnhcond=0.
!
! rows in the x-direction dealt implicitly
!
      do j=n1,n2
        wx=gamma*chi(l1:l2,j)
        ax=-dt/2.*wx*dx_2
        bx=1.-dt/2.*wx*dx_2*(-2.+dLnhcond(l1:l2,j)* &
           (TT(l1+1:l2+1,j)-2.*TT(l1:l2,j)+TT(l1-1:l2-1,j)))
        cx=-dt/2.*wx*dx_2
! rhsx=f_x(T^n) + f_z(T^n) + source
! do first f_z(T^n)
        rhsx=wx*dz_2*(TT(l1:l2,j+1)-2.*TT(l1:l2,j)+TT(l1:l2,j-1))
! then add f_x(T^n) + source
        rhsx=rhsx+wx*dx_2*(TT(l1+1:l2+1,j)-2.*TT(l1:l2,j)+TT(l1-1:l2-1,j)) &
             +source(l1:l2,j)
!
! periodic boundary conditions in x --> cyclic matrix
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax,bx,cx,aalpha,bbeta,rhsx,workx,nx)
        finter(l1:l2,j)=workx
      enddo
!
! columns in the z-direction dealt implicitly
!
      do i=l1,l2
        wz=dt*gamma*dz_2*chi(i,n1:n2)
        az=-wz/2.
        bz=1.-wz/2.*(-2.+dLnhcond(i,n1:n2)*    &
          (TT(i,n1+1:n2+1)-2.*TT(i,n1:n2)+TT(i,n1-1:n2-1)))
        cz=-wz/2.
        rhsz=finter(i,n1:n2)
!
! z boundary conditions
! Constant temperature at the top: T^(n+1)-T^n=0
       bz(nz)=1. ; az(nz)=0.
       rhsz(nz)=0.
! bottom
       select case (bcz1(ilnTT))
! Constant temperature at the bottom: T^(n+1)-T^n=0
         case ('cT')
          bz(1)=1. ; cz(1)=0.
          rhsz(1)=0.
! Constant flux at the bottom
         case ('c3')
          bz(1)=1. ; cz(1)=-1.
          rhsz(1)=0.
         case default 
          call fatal_error('ADI_Kprof_mixed','bcz on TT must be cT or c3')
       endselect
!
       call tridag(az,bz,cz,rhsz,workz)
       val(i,n1:n2)=workz(1:nz)
      enddo
!
      f(:,4,:,ilnTT)=f(:,4,:,iTTold)+dt*val
!
! update hcond used for the 'c3' condition in boundcond.f90
!
      call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof_mixed
!***********************************************************************
    subroutine ADI_Kprof_MPI_mixed(f)
!
!  01-mar-2010/dintrans: coded
!  parallel version of the ADI_Kprof_mixed subroutine
!
      use EquationOfState, only: gamma
      use Mpicomm, only: transp_xz, transp_zx
      use Boundcond, only: update_ghosts
!
      implicit none
!
      integer, parameter :: mxt=nx/nprocz+2*nghost, mzt=nzgrid+2*nghost
      integer, parameter :: n1t=nghost+1, n2t=n1t+nzgrid-1
      integer :: i, j, ibox, ll1, ll2
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz)  :: source, hcond, dhcond, finter, val, TT, &
                                 chi, dLnhcond
      real, dimension(mzt,mz) :: fintert, TTt, chit, dLnhcondt, valt
      real, dimension(nx)     :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nzgrid) :: az, bz, cz, wz, rhsz, workz
      real :: aalpha, bbeta
!
! needed for having the correct ghost zones for ilnTT
!
      call update_ghosts(f)
!
      source=(f(:,4,:,ilnTT)-f(:,4,:,iTTold))/dt
! BC important not for the x-direction (always periodic) but for 
! the z-direction as we must impose the 'c3' BC at the 2nd-order
! before going in the implicit stuff
      TT=f(:,4,:,iTTold)
      call heatcond_TT(TT, hcond, dhcond)
      call boundary_ADI(TT, hcond(:,n1))
      if (ldensity) then
        chi=cp1*hcond/exp(f(:,4,:,ilnrho))
!        chi=cp1*hcond0/exp(f(:,4,:,ilnrho))
      else
        chi=cp1*hcond
      endif
      dLnhcond=dhcond/hcond
!      dLnhcond=0.
!
! rows in the x-direction dealt implicitly
!
      do j=n1,n2
        wx=gamma*chi(l1:l2,j)
        ax=-dt/2.*wx*dx_2
        bx=1.-dt/2.*wx*dx_2*(-2.+dLnhcond(l1:l2,j)* &
           (TT(l1+1:l2+1,j)-2.*TT(l1:l2,j)+TT(l1-1:l2-1,j)))
        cx=-dt/2.*wx*dx_2
! rhsx=f_x(T^n) + f_z(T^n) + source
! do first f_z(T^n)
        rhsx=wx*dz_2*(TT(l1:l2,j+1)-2.*TT(l1:l2,j)+TT(l1:l2,j-1))
! then add f_x(T^n) + source
        rhsx=rhsx+wx*dx_2*(TT(l1+1:l2+1,j)-2.*TT(l1:l2,j)+TT(l1-1:l2-1,j)) &
             +source(l1:l2,j)
!
! periodic boundary conditions in x --> cyclic matrix
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax, bx, cx, aalpha, bbeta, rhsx, workx, nx)
        finter(l1:l2,j)=workx
      enddo
!
! do the transpositions x <--> z
!
      do ibox=1,nxgrid/nzgrid
        ll1=l1+(ibox-1)*nzgrid
        ll2=l1+ibox*nzgrid-1

!        call transp_xz(finter(ll1:ll2,n1:n2), fintert(n1t:n2t,n1:n2))
!        call transp_xz(chi(ll1:ll2,n1:n2), chit(n1t:n2t,n1:n2))
!        call transp_xz(dLnhcond(ll1:ll2,n1:n2), dLnhcondt(n1t:n2t,n1:n2))
!        call transp_xz(TT(ll1:ll2,n1:n2), TTt(n1t:n2t,n1:n2))
!
! columns in the z-direction dealt implicitly
! be careful! we still play with the l1,l2 and n1,n2 indices but that applies
! on *transposed* arrays
!
        do i=n1,n2
          wz=dt*gamma*dz_2*chit(n1t:n2t,i)
          az=-wz/2.
          bz=1.-wz/2.*(-2.+dLnhcondt(n1t:n2t,i)*    &
             (TTt(n1t+1:n2t+1,i)-2.*TTt(n1t:n2t,i)+TTt(n1t-1:n2t-1,i)))
          cz=-wz/2.
          rhsz=fintert(n1t:n2t,i)
!
! z boundary conditions
! Constant temperature at the top: T^(n+1)-T^n=0
!
          bz(nzgrid)=1. ; az(nzgrid)=0.
          rhsz(nzgrid)=0.
! bottom
          select case (bcz1(ilnTT))
! Constant temperature at the bottom: T^(n+1)-T^n=0
            case ('cT')
              bz(1)=1. ; cz(1)=0.
              rhsz(1)=0.
! Constant flux at the bottom
            case ('c3')
              bz(1)=1. ; cz(1)=-1.
              rhsz(1)=0.
            case default 
              call fatal_error('ADI_Kprof_MPI_mixed','bcz on TT must be cT or c3')
          endselect
          call tridag(az, bz, cz, rhsz, workz)
          valt(n1t:n2t,i)=workz
        enddo ! i
!
! come back on the grid (x,z)
!
!        call transp_xz(valt(n1t:n2t,n1:n2), val(ll1:ll2,n1:n2))
        f(ll1:ll2,4,n1:n2,ilnTT)=f(ll1:ll2,4,n1:n2,iTTold)+dt*val(ll1:ll2,n1:n2)
      enddo ! ibox
!
! update hcond used for the 'c3' condition in boundcond.f90
!
      if (iproc==0) call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof_MPI_mixed
!***********************************************************************
    subroutine ADI_Kconst_yakonov(f)
!
!  26-Jan-2011/dintrans: coded
!  2-D ADI scheme for the radiative diffusion term for a constant
!  radiative conductivity K. The ADI scheme is of Yakonov's form:
!
!    (1-dt/2*Lamba_x)*T^(n+1/2) = Lambda_x(T^n) + Lambda_z(T^n) + source
!    (1-dt/2*Lamba_z)*T^(n+1)   = T^(n+1/2)
!
!  where Lambda_x and Lambda_y denote diffusion operators and the source
!  term comes from the explicit advance.
!  Note: this form is more adapted for a parallelisation compared the 
!  Peaceman & Rachford one.
!
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top
      use Boundcond, only: update_ghosts
!
      implicit none
!
      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,mz) :: source, finter, val, TT, rho
      real, dimension(nx)    :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nz)    :: az, bz, cz, wz, rhsz, workz
      real :: aalpha, bbeta
!
      call update_ghosts(f)
!
      source=(f(:,4,:,ilnTT)-f(:,4,:,iTTold))/dt
      TT=f(:,4,:,iTTold)
      if (ldensity) then
        rho=exp(f(:,4,:,ilnrho))
      else
        rho=1.
      endif
!
!  rows dealt implicitly
!
      do j=n1,n2
        wx=dt*cp1*gamma*hcond0/rho(l1:l2,j)
        ax=-wx*dx_2/2.
        bx=1.+wx*dx_2
        cx=ax
        rhsx=TT(l1:l2,j)+ &
             wx*dz_2/2.*(TT(l1:l2,j+1)-2.*TT(l1:l2,j)+TT(l1:l2,j-1))
        rhsx=rhsx+wx*dx_2/2.*                                 &
             (TT(l1+1:l2+1,j)-2.*TT(l1:l2,j)+TT(l1-1:l2-1,j)) &
             +dt*source(l1:l2,j)
!
! x boundary conditions: periodic
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax,bx,cx,aalpha,bbeta,rhsx,workx,nx)
        finter(l1:l2,j)=workx
      enddo
!
!  columns dealt implicitly
!
      do i=l1,l2
        wz=dt*cp1*gamma*hcond0*dz_2/rho(i,n1:n2)
        az=-wz/2.
        bz=1.+wz
        cz=az
        rhsz=finter(i,n1:n2)
!
! z boundary conditions
!
! Constant temperature at the top
        bz(nz)=1. ; az(nz)=0.
        rhsz(nz)=cs2top/gamma_m1
! bottom
        select case (bcz1(ilnTT))
          ! Constant temperature at the bottom
          case ('cT')
            bz(1)=1. ; cz(1)=0.
            rhsz(1)=cs2bot/gamma_m1
          ! Constant flux at the bottom: c1 condition
          case ('c1')
            bz(1)=1.   ; cz(1)=-1
            rhsz(1)=dz*Fbot/hcond0
          case default 
            call fatal_error('ADI_Kprof_yakonov','bcz on TT must be cT or c1')
        endselect
!
        call tridag(az,bz,cz,rhsz,workz)
        f(i,4,n1:n2,ilnTT)=workz
      enddo
!
    endsubroutine ADI_Kconst_yakonov
!***********************************************************************
endmodule ImplicitPhysics
