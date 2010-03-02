! $Id$
!
! This module solves the radiative diffusion implicitly thanks
! to an Alternate Direction Implicit Scheme (ADI) in a D'Yakonov
! form
!     lambda_x T(n+1/2) = lambda_x + lambda_z
!     lambda_z T(n+1) = T(n+1/2)
!
!***************************************************************
module ImplicitPhysics
!
  use Cdata
  use Cparam
  use EquationOfState, only: mpoly0
  use Messages
  use Sub, only: keep_compiler_quiet
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
  real, pointer :: hcond0, Fbot, Tbump, Kmax, hole_slope, hole_width, &
                   hole_alpha
  logical, pointer :: lADI_mixed
!
  contains
!***********************************************************************
    subroutine init_param_ADI()
!
      use SharedVariables, only: get_shared_variable
      use MpiComm, only: stop_it
!
      implicit none
!
      integer :: ierr
!
      call get_shared_variable('hcond0', hcond0, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                 "there was a problem when getting hcond0")
      call get_shared_variable('Fbot', Fbot, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                "there was a problem when getting Fbot")
!
      call get_shared_variable('Tbump', Tbump, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                "there was a problem when getting Tbump")
      call get_shared_variable('Kmax', Kmax, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                "there was a problem when getting Kmax")
      call get_shared_variable('hole_slope', hole_slope, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                "there was a problem when getting hole_slope")
      call get_shared_variable('hole_alpha', hole_alpha, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                "there was a problem when getting hole_alpha")
      call get_shared_variable('hole_width', hole_width, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                "there was a problem when getting hole_width")
      call get_shared_variable('lADI_mixed', lADI_mixed, ierr)
      if (ierr/=0) call stop_it("implicit_physics: "//&
                "there was a problem when getting lADI_mixed")
      print*, 'lADI_mixed in implicit_physics=', lADI_mixed
!
    endsubroutine init_param_ADI
!***********************************************************************
    subroutine calc_heatcond_ADI(finit,f)
!
!  10-sep-07/gastine+dintrans: wrapper to the two possible ADI subroutines
!  ADI_Kconst: constant radiative conductivity
!  ADI_Kprof: radiative conductivity depends on T, i.e. hcond(T)
!
!
      implicit none
!
      real, dimension(mx,my,mz,mfarray) :: finit, f
      integer :: ierr
!
      if (hcond0 /= impossible) then
        if (nx == 1) then
          call ADI_Kconst_1d(finit,f)
        else
          if (nprocz>1) then
            call ADI_Kconst_MPI(finit,f)
          else
            call ADI_Kconst(finit,f)
          endif
        endif
      else
        if (nx == 1) then
          if (lADI_mixed) then
            call ADI_Kprof_1d_mixed(finit,f)
          else
            call ADI_Kprof_1d(finit,f)
          endif
        else
          if (nprocz>1) then
            if (lADI_mixed) then
              call ADI_Kprof_MPI_mixed(finit,f)
            else
              call ADI_Kprof_MPI(finit,f)
            endif
          else
            if (lADI_mixed) then
              call ADI_Kprof_mixed(finit,f)
            else
              call ADI_Kprof(finit,f)
            endif
          endif
        endif
      endif
!
    endsubroutine calc_heatcond_ADI
!***********************************************************************
    subroutine ADI_Kconst(finit,f)
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
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, get_cp1
      use General, only: tridag
      use Gravity, only: gravz
!
      implicit none
!
      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(nx,nz) :: finter, source, rho, TT
      real, dimension(nx)    :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nz)    :: az, bz, cz, wz, rhsz, workz
      real    :: aalpha, bbeta, cp1, dx_2, dz_2, tmp_flux
!
      TT=finit(l1:l2,4,n1:n2,ilnTT)
      source=(f(l1:l2,4,n1:n2,ilnTT)-TT)/dt
      if (ldensity) then
        rho=exp(f(l1:l2,4,n1:n2,ilnrho))
      else
        rho=1.
      endif
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
!
!  row dealt implicitly
!
      do j=1,nz
        wx=dt*gamma*hcond0*cp1/rho(:,j)
        ax=-wx*dx_2/2.
        bx=1.+wx*dx_2
        cx=ax
!
        if (j==1) then
          if (bcz1(ilnTT)=='cT') then
            !constant temperature: T(j-1)=2*T(j)-T(j+1)
            rhsx=TT(:,j)+dt/2.*source(:,j)
          else
            !imposed flux: T(j-1)=T(j+1)-2*dz*tmp_flux with tmp_flux=-Fbot/hcond0
            Fbot=-gamma/(gamma-1.)*hcond0*gravz/(mpoly0+1.)
            tmp_flux=-Fbot/hcond0
            rhsx=TT(:,j)+wx*dz_2/2.*                            &
                (TT(:,j+1)-2.*TT(:,j)+TT(:,j+1)-2.*dz*tmp_flux) &
                +dt/2.*source(:,j)
          endif
        elseif (j==nz) then
          !constant temperature: T(j+1)=2*T(j)-T(j-1)
          rhsx=TT(:,j)+dt/2.*source(:,j)
        else
          rhsx=TT(:,j)+wx*dz_2/2.*                         &
              (TT(:,j+1)-2.*TT(:,j)+TT(:,j-1))             &
               +dt/2.*source(:,j)
        endif
!
! x boundary conditions: periodic
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax, bx, cx, aalpha, bbeta, rhsx, workx, nx)
        finter(:,j)=workx
      enddo
!
!  columns dealt implicitly
!
      do i=1,nx
        wz=dt*gamma*hcond0*cp1/rho(i,:)
        az=-wz*dz_2/2.
        bz=1.+wz*dz_2
        cz=az
!
        if (i==1) then
          rhsz=finter(i,:)+wz*dx_2/2.*                       &
              (finter(i+1,:)-2.*finter(i,:)+finter(nx,:))    &
              +dt/2.*source(i,:)
        elseif (i==nx) then
          rhsz=finter(i,:)+wz*dx_2/2.*                       &
              (finter(1,:)-2.*finter(i,:)+finter(i-1,:))     &
              +dt/2.*source(i,:)
        else
          rhsz=finter(i,:)+wz*dx_2/2.*                       &
              (finter(i+1,:)-2.*finter(i,:)+finter(i-1,:))   &
              +dt/2.*source(i,:)
        endif
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
!           bz(1)=1.   ; cz(1)=-1
!           rhsz(1)=dz*Fbot/hcond0
! we can use here the second-order relation for the first derivative: 
! (T_{j+1}-T_{j_1})/2dz = dT/dz --> T_{j-1} = T_{j+1} - 2*dz*dT/dz 
! and insert this expression in the difference relation to eliminate T_{j-1}:
! a_{j-1}*T_{j-1} + b_j T_j + c_{j+1}*T_{j+1} = RHS
            cz(1)=cz(1)+az(1)
            rhsz(1)=rhsz(1)-2.*az(1)*dz*Fbot/hcond0
          case default 
           call fatal_error('ADI_Kconst','bcz on TT must be cT or c1')
        endselect
!
        call tridag(az, bz, cz, rhsz, workz)
        f(i+nghost,4,n1:n2,ilnTT)=workz
      enddo
!
    endsubroutine ADI_Kconst
!***********************************************************************
    subroutine ADI_Kprof(finit,f)
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
      use EquationOfState, only: gamma,get_cp1
      use General, only: tridag
!
      implicit none
!
      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: finit, f
      real, dimension(mx,mz) :: source, hcond, dhcond, finter, val, TT, rho
      real, dimension(nx)    :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nz)    :: az, bz, cz, wz, rhsz, workz
      real    :: aalpha, bbeta
      real    :: dx_2, dz_2, cp1
!
      source=(f(:,4,:,ilnTT)-finit(:,4,:,ilnTT))/dt
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
! BC important not for the x-direction (always periodic) but for 
! the z-direction as we must impose the 'c3' BC at the 2nd-order
! before going in the implicit stuff
      call heatcond_TT(finit(:,4,:,ilnTT), hcond, dhcond)
      call boundary_ADI(finit(:,4,:,ilnTT), hcond(:,n1))
      TT=finit(:,4,:,ilnTT)
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
      f(:,4,:,ilnTT)=finit(:,4,:,ilnTT)+dt*val
!
! update hcond used for the 'c3' condition in boundcond.f90
!
      call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof
!***********************************************************************
    subroutine ADI_Kprof_MPI(finit,f)
!
!  15-jan-10/gastine: coded
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
      use EquationOfState, only: gamma,get_cp1
      use General, only: tridag
      use Mpicomm, only: transp_mxmz, transp_xz, &
          transp_zx, initiate_isendrcv_bdry, finalize_isendrcv_bdry
!
      implicit none
!
      integer, parameter :: mxt=nx/nprocz+2*nghost
      integer, parameter :: mzt=nzgrid+2*nghost
      integer, parameter :: l1t=nghost+1, n1t=nghost+1
      integer, parameter :: l2t=l1t+nx/nprocz-1, n2t=n1t+nzgrid-1
      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: finit, f, ftmp
      real, dimension(mzt,my,mxt,mfarray) :: ftmpt
      real, dimension(mx,mz) :: source, hcond, dhcond, finter, TT, rho, val
      real, dimension(mzt,mxt) :: hcondt, dhcondt, fintert, TTt, rhot, valt
      real, dimension(nx)    :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nzgrid)    :: az, bz, cz, wz, rhsz, workz
      real :: dx_2, dz_2, cp1, aalpha, bbeta
!
!  It is necessary to communicate ghost-zones points between
!  processors to ensure a correct transposition of these ghost
!  zones. It is needed by rho,rhot and source,sourcet.
!
      call initiate_isendrcv_bdry(finit, ilnTT, ilnTT)
      call finalize_isendrcv_bdry(finit, ilnTT, ilnTT)
      source=(f(:,4,:,ilnTT)-finit(:,4,:,ilnTT))/dt
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
!
! BC important not for the x-direction (always periodic) but for 
! the z-direction as we must impose the 'c3' BC at the 2nd-order
! before going in the implicit stuff
!
      TT=finit(:,4,:,ilnTT)
      call heatcond_TT(TT, hcond, dhcond)
      call boundary_ADI(TT, hcond(:,n1))
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
! rhsx=f_y(T^n) + f_x(T^n) (Eq. 3.6)
! do first f_y(T^n)
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
! x boundary conditions: periodic
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax, bx, cx, aalpha, bbeta, rhsx, workx, nx)
        finter(l1:l2,j)=workx(1:nx)
      enddo
!
! do the transpositions x <--> z
!
      call transp_xz(finter(l1:l2,n1:n2), fintert(l1:l2,n1:n2))
      call transp_xz(rho(l1:l2,n1:n2), rhot(l1:l2,n1:n2))
      call transp_xz(hcond(l1:l2,n1:n2), hcondt(l1:l2,n1:n2))
      call transp_xz(dhcond(l1:l2,n1:n2), dhcondt(l1:l2,n1:n2))
      call transp_mxmz(TT, TTt)
!
      do i=l1t,l2t
        wz=dt*cp1*gamma*dz_2/rhot(n1t:n2t,i)
        az=-wz/4.*(dhcondt(n1t-1:n2t-1,i)   &
           *(TTt(n1t-1:n2t-1,i)-TTt(n1t:n2t,i)) &
           +hcondt(n1t-1:n2t-1,i)+hcondt(n1t:n2t,i))
!
        bz=1.+wz/4.*(dhcondt(n1t:n2t,i)*             &
           (2.*TTt(n1t:n2t,i)-TTt(n1t-1:n2t-1,i)         &
           -TTt(n1t+1:n2t+1,i))+2.*hcondt(n1t:n2t,i)     &
           +hcondt(n1t+1:n2t+1,i)+hcondt(n1t-1:n2t-1,i))
!
        cz=-wz/4.*(dhcondt(n1t+1:n2t+1,i)            &
           *(TTt(n1t+1:n2t+1,i)-TTt(n1t:n2t,i))          &
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
!
        call tridag(az, bz, cz, rhsz, workz)
        valt(n1t:n2t,i)=workz(1:nzgrid)
      enddo
      call transp_zx(valt(l1:l2,n1:n2), val(l1:l2,n1:n2))
      f(:,4,:,ilnTT)=finit(:,4,:,ilnTT)+dt*val
!
! update hcond used for the 'c3' condition in boundcond.f90
!
      if (iproc==0) then
        call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
      endif
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
      if (iproc == nprocz-1) then
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
    subroutine cyclic(a,b,c,alpha,beta,r,x,n)
!
! 08-Sep-07/gastine+dintrans: coded
! inversion of a tridiagonal matrix with periodic BC (alpha and beta
! coefficients); used in the ADI scheme
!
      use General, only: tridag
!
      implicit none
!
      integer :: i,n
      integer, parameter    :: NMAX=1200
      real    :: alpha, beta,gamma,fact      
      real, dimension(n)    :: a,b,c,r,x,bb,u,z
!     real, dimension(NMAX) :: bb,u,z
!
      if (n <= 2) call fatal_error('cyclic', 'n too small in cyclic')
      if (n > NMAX) call fatal_error('cyclic', 'NMAX too small in cyclic')
      gamma=-b(1)
      bb(1)=b(1)-gamma
      bb(n)=b(n)-alpha*beta/gamma
      do i=2,n-1
        bb(i)=b(i)
      enddo
      call tridag(a,bb,c,r,x)
      u(1)=gamma
      u(n)=alpha
      do i=2,n-1
        u(i)=0.
      enddo
      call tridag(a,bb,c,u,z)
      fact=(x(1)+beta*x(n)/gamma)/(1.+z(1)+beta*z(n)/gamma)
      do i=1,n
        x(i)=x(i)-fact*z(i)
      enddo
!
      return
    endsubroutine cyclic
!***********************************************************************
    subroutine ADI_Kconst_1d(finit,f)
!
! 18-sep-07/dintrans: coded
! Implicit Crank Nicolson scheme in 1-D for a constant K (not 
! really an ADI but keep the generic name for commodity).
!
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, get_cp1
      use General, only: tridag
!
      implicit none
!
      integer :: j, jj
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(mz) :: rho, TT
      real, dimension(nz) :: a, b, c, rhs, work
      real  :: cp1, dz_2, wz
!
!     source=(f(4,4,:,ilnTT)-finit(4,4,:,ilnTT))/dt
!     TT=finit(4,4,:,ilnTT)
      TT=f(4,4,:,ilnTT)
      rho=exp(f(4,4,:,ilnrho))
      call get_cp1(cp1)
      dz_2=1./dz**2
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
      call keep_compiler_quiet(finit)
!
    endsubroutine ADI_Kconst_1d
!***********************************************************************
    subroutine ADI_Kprof_1d(finit,f)
!
! 18-sep-07/dintrans: coded
! Implicit 1-D case for a temperature-dependent conductivity K(T).
! Not really an ADI but keep the generic name for commodity.
!
      use EquationOfState, only: gamma,get_cp1
      use General, only: tridag
!
      implicit none
!
      integer :: j, jj
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(mz) :: source, rho, TT, hcond, dhcond
      real, dimension(nz) :: a, b, c, rhs, work
      real  :: cp1, dz_2, wz, hcondp, hcondm
!
      source=(f(4,4,:,ilnTT)-finit(4,4,:,ilnTT))/dt
      call get_cp1(cp1)
      dz_2=1./dz**2
      rho=exp(f(4,4,:,ilnrho))
! need to set up the 'c3' BC at the 2nd-order before the implicit stuff
      call heatcond_TT(finit(4,4,:,ilnTT), hcond, dhcond)
      hcondADI=spread(hcond(1), 1, mx)
      call boundary_ADI(finit(:,4,:,ilnTT), hcondADI)
      TT=finit(4,4,:,ilnTT)
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
!
        rhs(jj)=wz/2.*(hcondp*(TT(j+1)-TT(j))-hcondm*(TT(j)-TT(j-1))) &
          +dt*source(j)
!
! Always constant temperature at the top: T^(n+1)-T^n = 0
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
      call tridag(a,b,c,rhs,work)
      f(4,4,n1:n2,ilnTT)=work+TT(n1:n2)
!
! Update the bottom value of hcond used for the 'c3' BC in boundcond
!
      call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof_1d
!***********************************************************************
    subroutine ADI_Kconst_MPI(finit,f)
!
!  04-sep-2009/dintrans: coded
!  parallel version of the ADI scheme for the K=cte case
!
      use EquationOfState, only: gamma, gamma_m1, cs2bot, cs2top, get_cp1
      use General, only: tridag
      use Mpicomm, only: transp_xz, transp_zx, MPI_adi_x, MPI_adi_z
!
      implicit none
!
      integer :: i,j
      integer, parameter :: nxt=nx/nprocz
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(nx,nz)  :: finter, source, rho, TT
      real, dimension(nx)     :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nzgrid) :: az, bz, cz, wz, rhsz, workz
      real, dimension(nzgrid,nxt) :: fintert, rhot, sourcet, wtmp
      real, dimension(nx)     :: tmpx1, tmpx2, send_bufx1, send_bufx2
      real, dimension(nzgrid) :: tmpz1, tmpz2, send_bufz1, send_bufz2
      real  :: aalpha, bbeta, cp1, dx_2, dz_2
!
      TT=finit(l1:l2,4,n1:n2,ilnTT)
      source=(f(l1:l2,4,n1:n2,ilnTT)-TT)/dt
      if (ldensity) then
        rho=exp(f(l1:l2,4,n1:n2,ilnrho))
      else
        rho=1.
      endif
      call get_cp1(cp1)
      dx_2=1/dx**2
      dz_2=1/dz**2
!
! Communicate the first and last pencils of size nx
!
      send_bufx1=TT(:,1)
      send_bufx2=TT(:,nz)
      call MPI_adi_x(tmpx1, tmpx2, send_bufx1, send_bufx2)
!
! Rows dealt implicitly
!
      do j=1,nz
        wx=dt*gamma*hcond0*cp1/rho(:,j)
        ax=-wx*dx_2/2
        bx=1.+wx*dx_2
        cx=ax
!
        if (j==1) then
          rhsx=TT(:,j)+wx*dz_2/2*          &
               (TT(:,j+1)-2*TT(:,j)+tmpx2) &
               +dt/2*source(:,j)
        elseif (j==nz) then
          rhsx=TT(:,j)+wx*dz_2/2*          &
               (tmpx1-2*TT(:,j)+TT(:,j-1)) &
               +dt/2*source(:,j)
        else
          rhsx=TT(:,j)+wx*dz_2/2*             &
              (TT(:,j+1)-2*TT(:,j)+TT(:,j-1)) &
               +dt/2*source(:,j)
        endif
!
! x boundary conditions: periodic
!
        aalpha=cx(nx) ; bbeta=ax(1)
        call cyclic(ax, bx, cx, aalpha, bbeta, rhsx, workx, nx)
        finter(:,j)=workx
      enddo
!
! Do the transpositions x <--> z
!
      call transp_xz(finter, fintert)
      call transp_xz(rho, rhot)
      call transp_xz(source, sourcet)
!
! Communicate the first and last pencils of size nzgrid
!
      send_bufz1=fintert(:,1)
      send_bufz2=fintert(:,nxt)
      call MPI_adi_z(tmpz1, tmpz2, send_bufz1, send_bufz2)
!
! Columns dealt implicitly
!
      do i=1,nxt
        wz=dt*gamma*hcond0*cp1/rhot(:,i)
        az=-wz*dz_2/2
        bz=1.+wz*dz_2
        cz=az
!
        if (i==1) then
          rhsz=fintert(:,i)+wz*dx_2/2*               &
              (fintert(:,i+1)-2*fintert(:,i)+tmpz2)  &
              +dt/2*sourcet(:,i)
        elseif (i==nxt) then
          rhsz=fintert(:,i)+wz*dx_2/2*               &
              (tmpz1-2*fintert(:,i)+fintert(:,i-1))  &
              +dt/2*sourcet(:,i)
        else
          rhsz=fintert(:,i)+wz*dx_2/2*                        &
              (fintert(:,i+1)-2*fintert(:,i)+fintert(:,i-1))  &
              +dt/2*sourcet(:,i)
        endif
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
!           bz(1)=1.   ; cz(1)=-1
!           rhsz(1)=dz*Fbot/hcond0
! we can use here the second-order relation for the first derivative: 
! (T_{j+1}-T_{j_1})/2dz = dT/dz --> T_{j-1} = T_{j+1} - 2*dz*dT/dz 
! and insert this expression in the difference relation to eliminate T_{j-1}:
! a_{j-1}*T_{j-1} + b_j T_j + c_{j+1}*T_{j+1} = RHS
!
            cz(1)=cz(1)+az(1)
            rhsz(1)=rhsz(1)-2.*az(1)*dz*Fbot/hcond0
          case default 
           call fatal_error('ADI_Kconst','bcz on TT must be cT or c1')
        endselect
!
        call tridag(az, bz, cz, rhsz, workz)
        wtmp(:,i)=workz
      enddo
      call transp_zx(wtmp, f(l1:l2,4,n1:n2,ilnTT))
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
      real, dimension(mx,mz) :: TT, arg, hcond
      real, dimension(mx,mz), optional :: dhcond
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
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
      real, dimension(:)           :: TT, hcond
      real, dimension(:), optional :: dhcond
      real, dimension(size(TT,1))  :: arg
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
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
      real :: TT, arg, hcond
      real, optional :: dhcond
!
      arg=hole_slope*(TT-Tbump-hole_width)*(TT-Tbump+hole_width)
      hcond=Kmax+hole_alpha*(-pi/2.+atan(arg))
      if (present(dhcond)) dhcond=2.*hole_alpha/(1.+arg**2)*hole_slope*(TT-Tbump)
!
    endsubroutine heatcond_TT_0d
!***********************************************************************
    subroutine ADI_Kprof_1d_mixed(finit,f)
!
! 28-feb-10/dintrans: coded
! Simpler version where a part of the radiative diffusion term is
! computed during the explicit advance
!
      use EquationOfState, only: gamma,get_cp1
      use General, only: tridag
!
      implicit none
!
      integer :: j, jj
      real, dimension(mx,my,mz,mfarray) :: finit,f
      real, dimension(mz) :: source, rho, TT, hcond, dhcond
      real, dimension(nz) :: a, b, c, rhs, work
      real  :: cp1, dz_2, wz, hcondp, hcondm
!
      source=(f(4,4,:,ilnTT)-finit(4,4,:,ilnTT))/dt
      call get_cp1(cp1)
      dz_2=1./dz**2
      rho=exp(f(4,4,:,ilnrho))
! need to set up the 'c3' BC at the 2nd-order before the implicit stuff
      call heatcond_TT(finit(4,4,:,ilnTT), hcond, dhcond)
      hcondADI=spread(hcond(1), 1, mx)
      call boundary_ADI(finit(:,4,:,ilnTT), hcondADI)
      TT=finit(4,4,:,ilnTT)
!
      do j=n1,n2
        jj=j-nghost
        wz=dt*dz_2*gamma*hcond(j)*cp1/rho(j)
!
        a(jj)=-wz/2.
        b(jj)=1.-wz/2.*(-2.+dhcond(j)/hcond(j)*(TT(j+1)-2.*TT(j)+TT(j-1)))
        c(jj)=-wz/2.
        rhs(jj)=wz*(TT(j+1)-2.*TT(j)+TT(j-1))+dt*source(j)
! Always constant temperature at the top: T^(n+1)-T^n = 0
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
      call tridag(a,b,c,rhs,work)
      f(4,4,n1:n2,ilnTT)=work+TT(n1:n2)
!
! Update the bottom value of hcond used for the 'c3' BC in boundcond
!
      call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof_1d_mixed
!***********************************************************************
    subroutine ADI_Kprof_mixed(finit,f)
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
      use EquationOfState, only: gamma,get_cp1
      use General, only: tridag
      use Boundcond, only: update_ghosts
!
      implicit none
!
      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: finit, f
      real, dimension(mx,mz) :: source, hcond, dhcond, finter, val, TT, &
                                rho, chi, dLnhcond
      real, dimension(nx)    :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nz)    :: az, bz, cz, wz, rhsz, workz
      real :: dx_2, dz_2, cp1, aalpha, bbeta
!
      call update_ghosts(finit)
!
      source=(f(:,4,:,ilnTT)-finit(:,4,:,ilnTT))/dt
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
! BC important not for the x-direction (always periodic) but for 
! the z-direction as we must impose the 'c3' BC at the 2nd-order
! before going in the implicit stuff
      call heatcond_TT(finit(:,4,:,ilnTT), hcond, dhcond)
      call boundary_ADI(finit(:,4,:,ilnTT), hcond(:,n1))
      TT=finit(:,4,:,ilnTT)
      if (ldensity) then
        rho=exp(f(:,4,:,ilnrho))
      else
        rho=1.
      endif
      chi=cp1*hcond/rho
      dLnhcond=dhcond/hcond
!      chi=cp1*hcond0/rho
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
      f(:,4,:,ilnTT)=finit(:,4,:,ilnTT)+dt*val
!
! update hcond used for the 'c3' condition in boundcond.f90
!
      call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
!
    endsubroutine ADI_Kprof_mixed
!***********************************************************************
    subroutine ADI_Kprof_MPI_mixed(finit,f)
!
!  01-mar-2010/dintrans: coded
!  parallel version of the ADI_Kprof_mixed subroutine
!
      use EquationOfState, only: gamma, get_cp1
      use General, only: tridag
      use Mpicomm, only: transp_mxmz, transp_xz, transp_zx
      use Boundcond, only: update_ghosts
!
      implicit none
!
      integer :: i,j
      real, dimension(mx,my,mz,mfarray) :: finit, f
      real, dimension(mx,mz) :: source, hcond, dhcond, finter, val, TT, &
            rho, chi, dLnhcond, fintert, TTt, chit, dLnhcondt, valt
      real, dimension(nx)     :: ax, bx, cx, wx, rhsx, workx
      real, dimension(nzgrid) :: az, bz, cz, wz, rhsz, workz
      real :: dx_2, dz_2, cp1, aalpha, bbeta
!
! needed for having the correct ghost zones for ilnTT
!
      call update_ghosts(finit)
!
      source=(f(:,4,:,ilnTT)-finit(:,4,:,ilnTT))/dt
      call get_cp1(cp1)
      dx_2=1./dx**2
      dz_2=1./dz**2
! BC important not for the x-direction (always periodic) but for 
! the z-direction as we must impose the 'c3' BC at the 2nd-order
! before going in the implicit stuff
      TT=finit(:,4,:,ilnTT)
      call heatcond_TT(TT, hcond, dhcond)
      call boundary_ADI(TT, hcond(:,n1))
      if (ldensity) then
        rho=exp(f(:,4,:,ilnrho))
      else
        rho=1.
      endif
      chi=cp1*hcond/rho
      dLnhcond=dhcond/hcond
!      chi=cp1*hcond0/rho
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
      call transp_xz(finter(l1:l2,n1:n2), fintert(l1:l2,n1:n2))
      call transp_xz(chi(l1:l2,n1:n2), chit(l1:l2,n1:n2))
      call transp_xz(dLnhcond(l1:l2,n1:n2), dLnhcondt(l1:l2,n1:n2))
      call transp_mxmz(TT, TTt)
!
! columns in the z-direction dealt implicitly
! be careful! we still play with the l1,l2 and n1,n2 indices but that applies
! on *transposed* arrays
!
      do i=n1,n2
        wz=dt*gamma*dz_2*chit(l1:l2,i)
        az=-wz/2.
        bz=1.-wz/2.*(-2.+dLnhcondt(l1:l2,i)*    &
           (TTt(l1+1:l2+1,i)-2.*TTt(l1:l2,i)+TTt(l1-1:l2-1,i)))
        cz=-wz/2.
        rhsz=fintert(l1:l2,i)
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
!
        call tridag(az, bz, cz, rhsz, workz)
        valt(l1:l2,i)=workz
      enddo
!
! come back on the grid (x,z)
!
      call transp_zx(valt(l1:l2,n1:n2), val(l1:l2,n1:n2))
      f(l1:l2,4,n1:n2,ilnTT)=TT(l1:l2,n1:n2)+dt*val(l1:l2,n1:n2)
!
! update hcond used for the 'c3' condition in boundcond.f90
!
      if (iproc==0) then
        call heatcond_TT(f(:,4,n1,ilnTT), hcondADI)
      endif
!
    endsubroutine ADI_Kprof_MPI_mixed
!***********************************************************************
endmodule ImplicitPhysics
