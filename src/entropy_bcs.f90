  module EnergyBcs

    use Cdata
    use DensityMethods
    use EquationOfState, only: get_gamma_etc, cs20, lnrho0, cs2bot, cs2top
    use Messages
!
    private

    real, dimension(:,:), pointer :: reference_state

    real :: gamma, gamma1, gamma_m1, cp, cv, cp1, cv1
    logical, pointer :: lheatc_kramers, lheatc_chiconst
    real, pointer :: chi,chi_t,hcondzbot,hcondztop,sigmaSBt
    real, pointer :: hcondxbot, hcondxtop    !, cs2bot, cs2top
    real, pointer :: chit_prof1,chit_prof2,hcond0_kramers, nkramers
    real, pointer :: Fbot, Ftop, FtopKtop,FbotKbot, mpoly
!
    character (len=labellen) :: meanfield_Beq_profile
    real, pointer :: meanfield_Beq, chit_quenching, uturb
    real, dimension(:), pointer :: B_ext
!
    integer, parameter :: XBOT=1, XTOP=nx

    include 'energy_bcs.h'
    include 'eos_params.h'

    contains
!**************************************************************************************************
    subroutine initialize_energy_bcs

      use EquationOfState, only: get_gamma_etc
      use SharedVariables, only: get_shared_variable

      call get_gamma_etc(gamma,cp,cv)
      gamma1=1./gamma; gamma_m1=gamma-1.
      cp1=1./cp; cv1=1./cv
!
!  Get the shared variables
!
      if (lreference_state) call get_shared_variable('reference_state',reference_state, &
                                                     caller='initialize_energy_bcs')
      call get_shared_variable('sigmaSBt',sigmaSBt,caller='initialize_energy_bcs')
      if (.not.lthermal_energy) then
        call get_shared_variable('Fbot',Fbot)
        call get_shared_variable('Ftop',Ftop)
        call get_shared_variable('FbotKbot',FbotKbot)
        call get_shared_variable('FtopKtop',FtopKtop)
        call get_shared_variable('lheatc_chiconst',lheatc_chiconst)
        call get_shared_variable('lheatc_kramers',lheatc_kramers)
        call get_shared_variable('chi_t',chi_t)
        call get_shared_variable('chit_prof1',chit_prof1)
        call get_shared_variable('chit_prof2',chit_prof2)
        if (lheatc_kramers) then
          call get_shared_variable('hcond0_kramers',hcond0_kramers)
          call get_shared_variable('nkramers',nkramers)
        endif
        call get_shared_variable('hcondzbot',hcondzbot)
        call get_shared_variable('hcondztop',hcondztop)
        call get_shared_variable('hcondxbot',hcondxbot)
        call get_shared_variable('hcondxtop',hcondxtop)
      else
        allocate(Fbot,Ftop,FbotKbot,FtopKtop)
        Fbot=0.; Ftop=0.; FbotKbot=0.; FbotKbot=0.
        allocate(lheatc_chiconst,lheatc_kramers); lheatc_chiconst=.false.; lheatc_kramers=.false.
        allocate(chi_t,chit_prof1,chit_prof2); chi_t=0.
        allocate(hcondzbot,hcondztop,hcondxbot,hcondxtop)
        hcondzbot=0.; hcondztop=0.; hcondxbot=0.; hcondxtop=0.
      endif
      !call get_shared_variable('cs2bot',cs2bot)
      !call get_shared_variable('cs2top',cs2top)
      call get_shared_variable('chi',chi)
    
!
      if (ldensity.and..not.lstratz) then
        call get_shared_variable('mpoly',mpoly)
      else
        call warning('initialize_eos','mpoly not obtained from density, set impossible')
        allocate(mpoly); mpoly=impossible
      endif
!
      if (lrun .and. lmagn_mf) then
        !call get_shared_variable('meanfield_Beq_profile',dummy)
        !meanfield_Beq_profile=dummy
        call get_shared_variable('meanfield_Beq',meanfield_Beq)
        call get_shared_variable('chit_quenching',chit_quenching)
        call get_shared_variable('uturb',uturb)
        call get_shared_variable('B_ext',B_ext)
      endif

    endsubroutine initialize_energy_bcs
!**************************************************************************************************
    subroutine bc_ss_flux(f,topbot,lone_sided)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!  13-mar-2011/pete: c1 condition for z-boundaries with Kramers' opacity
!   4-jun-2015/MR: factor cp added in front of tmp_xy
!  30-sep-2016/MR: changes for use of one-sided BC formulation (chosen by setting new optional switch lone_sided)
!
      use DensityMethods, only: getdlnrho_z, getderlnrho_z
      use Deriv, only: bval_from_neumann, set_ghosts_for_onesided_ders
      use General, only: loptest
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      logical, optional :: lone_sided
!
      real, dimension (size(f,1),size(f,2)) :: tmp_xy,cs2_xy,rho_xy
      integer :: i
!
      if (ldebug) print*,'bc_ss_flux: ENTER - cs20=',cs20
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!  ===============
!
      case(BOT)
!
!  calculate Fbot/(K*cs2)
!
        if (pretend_lnTT) then
          tmp_xy=-FbotKbot/exp(f(:,:,n1,iss))
          do i=1,nghost
            f(:,:,n1-i,iss)=f(:,:,n1+i,iss)-dz2_bound(-i)*tmp_xy
          enddo
        else
!
          call getrho(f(:,:,n1,ilnrho),rho_xy)
          cs2_xy = f(:,:,n1,iss)         ! here cs2_xy = entropy
          if (lreference_state) &
            cs2_xy(l1:l2,:) = cs2_xy(l1:l2,:) + spread(reference_state(:,iref_s),2,my)
!
          if (ldensity_nolog) then
            cs2_xy=cs20*exp(gamma_m1*(log(rho_xy)-lnrho0)+cv1*cs2_xy)
          else
            cs2_xy=cs20*exp(gamma_m1*(f(:,:,n1,ilnrho)-lnrho0)+cv1*cs2_xy)
          endif
!
!  Check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!  Check also whether Kramers opacity is used, then hcond itself depends
!  on density and temperature.
!
          if (lheatc_chiconst) then
            tmp_xy=Fbot/(rho_xy*chi*cs2_xy)
          else if (lheatc_kramers) then
            tmp_xy=Fbot*rho_xy**(2*nkramers)*(cp*gamma_m1)**(6.5*nkramers) &
                   /(hcond0_kramers*cs2_xy**(6.5*nkramers+1.))
          else
            tmp_xy=cp*FbotKbot/cs2_xy      ! when FbotKbot= Fbot/Kbot then factor cp tb added
          endif
!
!  enforce ds/dz + (cp-cv)*dlnrho/dz = - cp*(cp-cv)*Fbot/(Kbot*cs2)
!
          if (loptest(lone_sided)) then
            call not_implemented('bc_ss_flux', 'one-sided BC')
            call getderlnrho_z(f,n1,rho_xy)                           ! rho_xy=d_z ln(rho)
            call bval_from_neumann(f,topbot,iss,3,rho_xy)
            call set_ghosts_for_onesided_ders(f,topbot,iss,3,.true.)
          else
            do i=1,nghost
              call getdlnrho_z(f(:,:,:,ilnrho),n1,i,rho_xy)           ! rho_xy=del_z ln(rho)
              f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+(cp-cv)*(rho_xy+dz2_bound(-i)*tmp_xy)   !! factor cp removed
            enddo
          endif
        endif
!
!  top boundary
!  ============
!
      case(TOP)
!
!  calculate Ftop/(K*cs2)
!
        if (pretend_lnTT) then
          tmp_xy=-FtopKtop/exp(f(:,:,n2,iss))
          do i=1,nghost
             f(:,:,n2-i,iss)=f(:,:,n2+i,iss)-dz2_bound(i)*tmp_xy
          enddo
        else
!
          call getrho(f(:,:,n2,ilnrho),rho_xy)
          cs2_xy = f(:,:,n2,iss)             ! here cs2_xy = entropy
          if (lreference_state) &
            cs2_xy(l1:l2,:) = cs2_xy(l1:l2,:) + spread(reference_state(:,iref_s),2,my)
!
          if (ldensity_nolog) then
            cs2_xy=cs20*exp(gamma_m1*(log(rho_xy)-lnrho0)+cv1*cs2_xy)
          else
            cs2_xy=cs20*exp(gamma_m1*(f(:,:,n2,ilnrho)-lnrho0)+cv1*cs2_xy)
          endif
!
!  Check whether we have chi=constant at top, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!  Check also whether Kramers opacity is used, then hcond itself depends
!  on density and temperature.
!
          if (lheatc_chiconst) then
            tmp_xy=Ftop/(rho_xy*chi*cs2_xy)
          else if (lheatc_kramers) then
            tmp_xy=Ftop*rho_xy**(2*nkramers)*(cp*gamma_m1)**(6.5*nkramers) &
                   /(hcond0_kramers*cs2_xy**(6.5*nkramers+1.))
          else
            tmp_xy=cp*FtopKtop/cs2_xy        !! factor cp added
          endif
!
!  enforce ds/dz + (cp-cv)*dlnrho/dz = - cp*(cp-cv)*Ftop/(K*cs2)
!
          if (loptest(lone_sided)) then
            call not_implemented('bc_ss_flux', 'one-sided BC')
            call getderlnrho_z(f,n2,rho_xy)                           ! rho_xy=d_z ln(rho)
            call bval_from_neumann(f,topbot,iss,3,rho_xy)
            call set_ghosts_for_onesided_ders(f,topbot,iss,3,.true.)
          else
            do i=1,nghost
              call getdlnrho_z(f(:,:,:,ilnrho),n2,i,rho_xy)        ! rho_xy=del_z ln(rho)
              f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+(cp-cv)*(-rho_xy-dz2_bound(i)*tmp_xy)   !! factor cp removed!
            enddo
          endif
        endif
!
      case default
        call fatal_error('bc_ss_flux','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_flux_turb(f,topbot)
!
!  Black body boundary condition for entropy (called when bcz='Fgs')
!  setting F = sigmaSBt*T^4 where sigmaSBt is related to the
!  Stefan-Boltzmann constant.
!
!   04-may-2009/axel: adapted from bc_ss_flux
!   31-may-2010/pete: replaced sigmaSB by a `turbulent' sigmaSBt
!    4-jun-2015/MR: corrected sign of dsdz_xy for bottom boundary;
!                   added branches for Kramers heat conductivity (using sigmaSBt!)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      real, dimension (size(f,1),size(f,2)) :: dsdz_xy,cs2_xy,lnrho_xy,rho_xy, &
                                               TT_xy,dlnrhodz_xy,chi_xy
      integer :: i
!
      if (ldebug) print*,'bc_ss_flux_turb: ENTER - cs20=',cs20
!
!  Get the shared variables for magnetic quenching effect in a
!  mean-field description of a radiative boundary condition.
!  Ideally, one would like this to reside in magnetic/meanfield,
!  but this leads currently to circular dependencies.
!
!
!  lmeanfield_chitB and chi_t0
!
!     if (lmagnetic) then
!       call get_shared_variable('lmeanfield_chitB',lmeanfield_chitB)
!       if (lmeanfield_chitB) then
!         call get_shared_variable('chi_t0',chi_t0)
!       else
!         lmeanfield_chitB=.false.
!       endif
!     endif
!
      select case (topbot)
!
!  Check whether we want to do top or bottom (this is precessor dependent)
!
!  bottom boundary
!  ===============
!
      case(BOT)
!
!  set ghost zones such that dsdz_xy obeys
!  - chi_t rho T dsdz_xy - hcond gTT = sigmaSBt*TT^4
!
        call getlnrho(f(:,:,n1,ilnrho),lnrho_xy)
        cs2_xy=gamma_m1*(lnrho_xy-lnrho0)+cv1*f(:,:,n1,iss)
        if (lreference_state) cs2_xy(l1:l2,:)=cs2_xy(l1:l2,:)+cv1*spread(reference_state(:,iref_s),2,my)
        cs2_xy=cs20*exp(cs2_xy)
!
        call getrho(f(:,:,n1,ilnrho),rho_xy)
        TT_xy=cs2_xy/(gamma_m1*cp)
!
        dlnrhodz_xy= coeffs_1_z(1,1)*(f(:,:,n1+1,ilnrho)-f(:,:,n1-1,ilnrho)) &
                    +coeffs_1_z(2,1)*(f(:,:,n1+2,ilnrho)-f(:,:,n1-2,ilnrho)) &
                    +coeffs_1_z(3,1)*(f(:,:,n1+3,ilnrho)-f(:,:,n1-3,ilnrho))
!
        if (lheatc_kramers) then
          dsdz_xy= cv*( (sigmaSBt/hcond0_kramers)*TT_xy**(3-6.5*nkramers)*rho_xy**(2.*nkramers) &
                       +gamma_m1*dlnrhodz_xy)      ! no turbulent diffusivity considered here!
        else
          dsdz_xy=(sigmaSBt*TT_xy**3+hcondzbot*gamma_m1*dlnrhodz_xy)/ &
                  (chit_prof1*chi_t*rho_xy+hcondzbot/cv)
        endif
!
!  enforce ds/dz=-(sigmaSBt*T^3 + hcond*(gamma-1)*glnrho)/(chi_t*rho+hcond/cv)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+dz2_bound(-i)*dsdz_xy
        enddo
!
!  top boundary
!  ============
!
      case(TOP)
!
!  set ghost zones such that dsdz_xy obeys
!  - chi_t rho T dsdz_xy - hcond gTT = sigmaSBt*TT^4
!
        if (ilnrho/=0) then
          call getlnrho(f(:,:,n2,ilnrho),lnrho_xy)            ! here rho_xy=log(rho)
          cs2_xy=gamma_m1*(lnrho_xy-lnrho0)+cv1*f(:,:,n2,iss) ! this is lncs2
          if (lreference_state) &
            cs2_xy(l1:l2,:)=cs2_xy(l1:l2,:) + cv1*spread(reference_state(:,iref_s),2,my)
          cs2_xy=cs20*exp(cs2_xy)
          call getrho(f(:,:,n2,ilnrho),rho_xy)                ! here rho_xy=rho
        endif
!
!  This TT_xy is used for b.c. used in BB14
!
        TT_xy=cs2_xy/(gamma_m1*cp)
        dlnrhodz_xy = coeffs_1_z(1,2)*(f(:,:,n2+1,ilnrho)-f(:,:,n2-1,ilnrho)) &
                     +coeffs_1_z(2,2)*(f(:,:,n2+2,ilnrho)-f(:,:,n2-2,ilnrho)) &
                     +coeffs_1_z(3,2)*(f(:,:,n2+3,ilnrho)-f(:,:,n2-3,ilnrho))
!
        if (ldensity_nolog) dlnrhodz_xy=dlnrhodz_xy/rho_xy   ! (not used)
!
!  Set chi_xy=chi, which sets also the ghost zones.
!  chi_xy consists of molecular and possibly turbulent values.
!  The turbulent value can be quenched (but not in ghost zones).
!
      chi_xy=chi
      if (lmagnetic) then
        !if (lmeanfield_chitB) then
        !  n=n2
        !  do m=m1,m2
        !    call bdry_magnetic(f,quench,'meanfield_chitB')
        !  enddo
        !  chi_xy(l1:l2,m)=chi+chi_t0*quench
        !endif
      endif
!
!  Select to use either: sigmaSBt*TT^4 = - K dT/dz - chi_t*rho*T*ds/dz,
!                    or: sigmaSBt*TT^4 = - chi_xy*rho*cp dT/dz - chi_t*rho*T*ds/dz.
!
      if (lheatc_kramers) then
        dsdz_xy=-cv*(sigmaSBt*TT_xy**3+hcond0_kramers*TT_xy**(6.5*nkramers)*rho_xy**(-2.*nkramers)* &
                gamma_m1*dlnrhodz_xy)/(hcond0_kramers*TT_xy**(6.5*nkramers)*rho_xy**(-2.*nkramers) &
                       + chit_prof2*chi_t*rho_xy/gamma)
      elseif (hcondztop==impossible) then
        dsdz_xy=-(sigmaSBt*TT_xy**3+chi_xy*rho_xy*cp*gamma_m1*dlnrhodz_xy)/ &
                 (chit_prof2*chi_t*rho_xy+chi_xy*rho_xy*cp/cv)
      else
!
!  This is what was used in the BB14 paper.
!
        dsdz_xy=-(sigmaSBt*TT_xy**3+hcondztop*gamma_m1*dlnrhodz_xy)/ &
                 (chit_prof2*chi_t*rho_xy+hcondztop/cv)
      endif
!
!  Apply condition;
!  enforce ds/dz=-(sigmaSBt*T^3 + hcond*(gamma-1)*glnrho)/(chi_t*rho+hcond/cv)
!
      do i=1,nghost
        f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+dz2_bound(i)*dsdz_xy
      enddo
!
!  capture undefined entries
!
      case default
        call fatal_error('bc_ss_flux_turb','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_flux_turb
!***********************************************************************
    subroutine bdry_magnetic(f,quench,task)
!
!  Calculate magnetic properties needed for z boundary conditions.
!  This routine contains calls to more specialized routines.
!
!   8-jun-13/axel: coded, originally in magnetic, but cyclic dependence
!  21-jan-15/MR  : changes for use of reference state.
!
      use Sub, only: curl, dot2
      !use Boundcond, only: boundconds_x, boundconds_y, boundconds_z
      !use Mpicomm, only: initiate_isendrcv_bdry, finalize_isendrcv_bdry
      !use Magnetic_meanfield, only: meanfield_chitB
!
      real, dimension (:,:,:,:), intent(in) :: f
      real, dimension (:),       intent(out):: quench
      character (len=*), intent(in) :: task
!
      real, dimension (size(quench),3) :: bb
      real, dimension (size(quench)) :: rho,b2
      integer :: j
!
      !character (len=linelen), pointer :: dummy
!
      select case (task)
!
      case ('meanfield_chitB')
!
        !call boundconds_x(f,iax,iaz)
        !call initiate_isendrcv_bdry(f,iax,iaz)
        !call finalize_isendrcv_bdry(f,iax,iaz)
        !call boundconds_y(f,iax,iaz)
        !call boundconds_z(f,iax,iaz)
!
!  Add the external field.
!
        call curl(f,iaa,bb)
        do j=1,3
          bb(:,j)=bb(:,j)!+B_ext(j)
        enddo
        call dot2(bb,b2)
        call getrho(f(:,m,n,ilnrho),rho)
!
!  Call mean-field routine.
!
        call meanfield_chitB(rho,b2,quench)
!
!  capture undefined entries
!
      case default
        call fatal_error('bdry_magnetic','invalid value of topbot')
      endselect
!
    endsubroutine bdry_magnetic
!***********************************************************************
    subroutine meanfield_chitB(rho,b2,quench)
!
!  Calculate magnetic properties needed for z boundary conditions.
!  This routine contails calls to more specialized routines.
!
!   8-jun-13/axel: coded
!
      real, dimension(:), intent(IN) :: rho,b2
      real, dimension(:), intent(OUT):: quench
!
      real, dimension(size(rho)) :: Beq21
!
!  compute Beq21 = 1/Beq^2
!XX
!     select case (meanfield_Beq_profile)
!     case ('uturbconst');
        Beq21=mu01/(rho*uturb**2)
!     case default;
!       Beq21=1./meanfield_Beq**2
!     endselect
!
!  compute chit_quenching
!
      quench=1./(1.+chit_quenching*b2*Beq21)
!
    endsubroutine meanfield_chitB
!***********************************************************************
    subroutine bc_ss_flux_turb_x(f,topbot)
!
!  Black body boundary condition for entropy (called when bcx='Fgs')
!  setting F = sigmaSBt*T^4 where sigmaSBt is related to the
!  Stefan-Boltzmann constant.
!
!   31-may-2010/pete: adapted from bc_ss_flux_turb
!   20-jul-2010/pete: expanded to take into account hcond/=0
!   21-jan-2015/MR: changes for reference state.
!   22-jan-2015/MR: corrected bug in branches for pretend_lnTT=T
!    6-jun-2015/MR: added branches for Kramers heat conductivity (using sigmaSBt!)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,2),size(f,3)) :: dsdx_yz,cs2_yz,rho_yz,dlnrhodx_yz,TT_yz
      real, dimension (size(f,2),size(f,3)) :: hcond_total
      integer :: i
!
      if (ldebug) print*,'bc_ss_flux_turb: ENTER - cs20=',cs20
!
      select case (topbot)
!
!  Check whether we want to do top or bottom (this is processor dependent)
!
!  bottom boundary
!  ===============
!
      case(BOT)
!
! For the case of pretend_lnTT=T, set glnTT=-sigma*T^3/hcond
!
        if (pretend_lnTT) then
          do i=1,nghost
            f(l1-i,:,:,iss)=f(l1+i,:,:,iss) + &     ! MR: corrected, plus sign correct?
                dx2_bound(-i)*sigmaSBt*exp(f(l1,:,:,iss))**3/hcondxbot
          enddo
        else
!
!  set ghost zones such that dsdx_yz obeys
!  - chi_t rho T dsdx_yz - hcond gTT = sigmaSBt*TT^4
!
          call getrho(f(l1,:,:,ilnrho),XBOT,rho_yz)
!
          if (ldensity_nolog) then
            cs2_yz=f(l1,:,:,iss)                    ! here cs2_yz = entropy
            if (lreference_state) cs2_yz = cs2_yz+reference_state(XBOT,iref_s)
            cs2_yz=cs20*exp(gamma_m1*(log(rho_yz)-lnrho0)+cv1*cs2_yz)
          else
            cs2_yz=cs20*exp(gamma_m1*(f(l1,:,:,ilnrho)-lnrho0)+cv1*f(l1,:,:,iss))
          endif
!
          TT_yz=cs2_yz/(gamma_m1*cp)
!
!  Calculate   d rho/d x    or   d ln(rho) / dx
!
          dlnrhodx_yz= coeffs_1_x(1,1)*(f(l1+1,:,:,ilnrho)-f(l1-1,:,:,ilnrho)) &
                      +coeffs_1_x(2,1)*(f(l1+2,:,:,ilnrho)-f(l1-2,:,:,ilnrho)) &
                      +coeffs_1_x(3,1)*(f(l1+3,:,:,ilnrho)-f(l1-3,:,:,ilnrho))
!
          if (ldensity_nolog) then
!
!  Add gradient of reference density and divide by total density
!
            if (lreference_state) then
              dlnrhodx_yz=dlnrhodx_yz + reference_state(XBOT,iref_grho)
              dlnrhodx_yz=dlnrhodx_yz/(rho_yz + reference_state(XBOT,iref_rho))
            else
              dlnrhodx_yz=dlnrhodx_yz/rho_yz
            endif
!
          endif
!
          if (lheatc_kramers) then
            dsdx_yz=-cv*( (sigmaSBt/hcond0_kramers)*TT_yz**(3-6.5*nkramers)*rho_yz**(2.*nkramers) &
                         +gamma_m1*dlnrhodx_yz)      ! no turbulent diffusivity considered here!
          else
            dsdx_yz=-(sigmaSBt*TT_yz**3+hcondxbot*gamma_m1*dlnrhodx_yz)/ &
                     (chit_prof1*chi_t*rho_yz+hcondxbot/cv)
          endif
!
!  Substract gradient of reference entropy.
!
          if (lreference_state) dsdx_yz = dsdx_yz - reference_state(XBOT,iref_gs)
!
!  enforce ds/dx = - (sigmaSBt*T^3 + hcond*(gamma-1)*glnrho)/(chi_t*rho+hcond/cv)
!
          do i=1,nghost
            f(l1-i,:,:,iss)=f(l1+i,:,:,iss)-dx2_bound(-i)*dsdx_yz
          enddo
        endif
!
!  top boundary
!  ============
!
      case(TOP)
!
!  For the case of pretend_lnTT=T, set glnTT=-sigma*T^3/hcond
!
        if (pretend_lnTT) then
          do i=1,nghost
            f(l2+i,:,:,iss)=f(l2-i,:,:,iss) + &     ! MR: corrected, plus sign correct?
                            dx2_bound(i)*sigmaSBt*exp(f(l2,:,:,iss))**3/hcondxtop
          enddo
        else
!
!  set ghost zones such that dsdx_yz obeys
!  - chi_t rho T dsdx_yz - hcond gTT = sigmaSBt*TT^4
!
          call getrho(f(l2,:,:,ilnrho),XTOP,rho_yz)
!
          if (ldensity_nolog) then
            cs2_yz=f(l2,:,:,iss)                         ! here cs2_yz = entropy
            if (lreference_state) cs2_yz = cs2_yz+reference_state(XTOP,iref_s)
            cs2_yz=cs20*exp(gamma_m1*(log(rho_yz)-lnrho0)+cv1*cs2_yz)
          else
            cs2_yz=cs20*exp(gamma_m1*(f(l2,:,:,ilnrho)-lnrho0)+cv1*f(l2,:,:,iss))
          endif
!
          TT_yz=cs2_yz/(gamma_m1*cp)
!
!  Calculate d rho/d x  or  d ln(rho) / d x
!
          dlnrhodx_yz= coeffs_1_x(1,2)*(f(l2+1,:,:,ilnrho)-f(l2-1,:,:,ilnrho)) &
                      +coeffs_1_x(2,2)*(f(l2+2,:,:,ilnrho)-f(l2-2,:,:,ilnrho)) &
                      +coeffs_1_x(3,2)*(f(l2+3,:,:,ilnrho)-f(l2-3,:,:,ilnrho))
!
          if (ldensity_nolog) then
!
!  Add gradient of reference density to d rho/d x and divide by total density
!
            if (lreference_state) then
              dlnrhodx_yz=dlnrhodx_yz + reference_state(XTOP,iref_grho)
              dlnrhodx_yz=dlnrhodx_yz/(rho_yz + reference_state(XTOP,iref_rho))
            else
              dlnrhodx_yz=dlnrhodx_yz/rho_yz
            endif
!
          endif
!
!  Compute total heat conductivity
!
          if (lheatc_kramers .or. (hcondxtop /= 0.0)) then
            hcond_total = hcondxtop
            if (lheatc_kramers) hcond_total = hcond_total + &
                hcond0_kramers*TT_yz**(6.5*nkramers)*rho_yz**(-2.*nkramers)
!
            dsdx_yz = -(sigmaSBt*TT_yz**3+hcond_total*gamma_m1*dlnrhodx_yz) / &
                       (chit_prof2*chi_t*rho_yz+hcond_total/cv)
!
!  Substract gradient of reference entropy.
!
            if (lreference_state) dsdx_yz = dsdx_yz - reference_state(XTOP,iref_gs)
!
!  enforce ds/dx = - (sigmaSBt*T^3 + hcond*(gamma-1)*glnrho)/(chi_t*rho+hcond/cv)
!
            do i=1,nghost
              f(l2+i,:,:,iss)=f(l2-i,:,:,iss)+dx2_bound(i)*dsdx_yz
            enddo
          endif
        endif
!
!  capture undefined entries
!
      case default
        call fatal_error('bc_ss_flux_turb_x','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_flux_turb_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_x(f,topbot)
!
!   Constant conductive + turbulent flux through the surface
!
!   08-apr-2014/pete: coded
!   21-jan-2015/MR: changes for reference state.
!    4-jun-2015/MR: added missing factor cp in RB;
!                   added branches for Kramers heat conductivity
!
      use DensityMethods, only: getdlnrho_x
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: dsdx_yz, TT_yz, rho_yz, dlnrhodx_yz, Kxbot
      integer :: i
!
      if (ldebug) print*,'bc_ss_flux_condturb: ENTER - cs20=',cs20
!
      select case (topbot)
!
!  Check whether we want to do top or bottom (this is precessor dependent)
!
!  bottom boundary
!  ===============
!
      case(BOT)
!
! Do the pretend_lnTT=T case first
!
        if (pretend_lnTT) then
          call not_implemented("bc_ss_flux_condturb_x","for pretend_lnTT=T")
        else
!
!  Set ghost zones such that -chi_t*rho*T*grads -hcond*gTT = Fbot
!
          call getrho(f(l1,:,:,ilnrho),XBOT,rho_yz)
!
          if (ldensity_nolog) then
            TT_yz=f(l1,:,:,iss)                                      ! TT_yz here entropy
            if (lreference_state) TT_yz = TT_yz+reference_state(XBOT,iref_s)
            TT_yz=cs20*exp(gamma_m1*(log(rho_yz)-lnrho0)+cv1*TT_yz)  ! TT_yz here cs2
          else
            TT_yz=cs20*exp(gamma_m1*(f(l1,:,:,ilnrho)-lnrho0)+cv1*f(l1,:,:,iss))
          endif
          TT_yz=TT_yz/(cp*gamma_m1)                                  ! TT_yz here temperature
!
!  Calculate d rho/d x   or   d ln(rho)/ d x
!
          dlnrhodx_yz= coeffs_1_x(1,1)*(f(l1+1,:,:,ilnrho)-f(l1-1,:,:,ilnrho)) &
                      +coeffs_1_x(2,1)*(f(l1+2,:,:,ilnrho)-f(l1-2,:,:,ilnrho)) &
                      +coeffs_1_x(3,1)*(f(l1+3,:,:,ilnrho)-f(l1-3,:,:,ilnrho))
!
          if (ldensity_nolog) then
!
!  Add gradient of reference density and divide by total density (but not used later!).
!
            if (lreference_state) then
              dlnrhodx_yz=dlnrhodx_yz + reference_state(XBOT,iref_grho)
              dlnrhodx_yz=dlnrhodx_yz/(rho_yz + reference_state(XBOT,iref_rho))
            else
              dlnrhodx_yz=dlnrhodx_yz/rho_yz
            endif
!
          endif
!
          if (lheatc_kramers) then
            Kxbot=hcond0_kramers*TT_yz**(6.5*nkramers)/rho_yz**(2.*nkramers)
          else
            Kxbot=hcondxbot
          endif
!
          dsdx_yz=(Fbot/TT_yz)/(chit_prof1*chi_t*rho_yz + Kxbot*cv1)
!
!  Add gradient of reference entropy.
!
          if (lreference_state) dsdx_yz = dsdx_yz + reference_state(XBOT,iref_gs)
!
!  Enforce ds/dx = -(cp*gamma_m1*Fbot/cs2 + K*gamma_m1*glnrho)/(gamma*K+chi_t*rho)
!
          do i=1,nghost
            call getdlnrho_x(f(:,:,:,ilnrho),l1,i,rho_yz,dlnrhodx_yz)
            f(l1-i,:,:,iss)=f(l1+i,:,:,iss) + Kxbot*gamma_m1/(Kxbot*cv1+chit_prof1*chi_t*rho_yz)* &
                            dlnrhodx_yz + dx2_bound(-i)*dsdx_yz
          enddo
        endif
!
!  top boundary
!  ============
!
      case(TOP)
!
         call not_implemented("bc_ss_flux_condturb_x","for the top boundary")
!
!  capture undefined entries
!
      case default
        call fatal_error('bc_ss_flux_condturb_x','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_flux_condturb_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_mean_x(f,topbot)
!
!   Constant conductive + turbulent flux through the surface applied on
!   the spherically symmetric part, zero gradient for the fluctuation
!   at the boundary.
!
!   18-dec-2014/pete: coded
!
      use Mpicomm, only: mpiallreduce_sum
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (my,mz) :: dsdx_yz
      real, dimension (-nghost:nghost) :: lnrmx, lnrmx_tmp
      real :: cs2mx, cs2mx_tmp
      real :: fact, dlnrmxdx, tmp1
      real, dimension(ny) :: tmp2
      integer :: i,l,lf
!
      if (ldebug) print*,'bc_ss_flux_condturb_mean_x: ENTER - cs20=',cs20
!
      select case (topbot)
!
!  Check whether we want to do top or bottom (this is precessor dependent)
!
!  bottom boundary
!  ===============
!
      case(BOT)
!
! Do the pretend_lnTT=T case first
!
        if (pretend_lnTT) then
           call not_implemented("bc_ss_flux_condturb_mean_x","for pretend_lnTT=T")
        else
!
!  Compute yz-averaged log density and sound speed
!
          fact=1./((cos(y0)-cos(y0+Lxyz(2)))*Lxyz(3))
          cs2mx=0.
          do n=n1,n2
            call getlnrho(f(l1,:,n,ilnrho),XBOT,tmp2)
            tmp2 = gamma_m1*(tmp2-lnrho0) + cv1*f(l1,m1:m2,n,iss)
            if (lreference_state) tmp2 = tmp2 + cv1*reference_state(XBOT,iref_s)
            cs2mx = cs2mx+sum(cs20*exp(tmp2)*dVol_y(m1:m2))*dVol_z(n)
          enddo
          cs2mx=fact*cs2mx
!
          lnrmx=0.
          fact=1./((cos(y0)-cos(y0+Lxyz(2)))*Lxyz(3))
          do l=-nghost,nghost
            tmp1=0.
            lf=l+nghost+1
            do n=n1,n2
              call getlnrho(f(lf,:,n,ilnrho),XBOT,tmp2)         ! doubled call to getlnrho not yet optimal
              tmp1=tmp1+sum(tmp2*dVol_y(m1:m2))*dVol_z(n)
            enddo
            lnrmx(l)=lnrmx(l)+tmp1
          enddo
          lnrmx=fact*lnrmx
!
!  Communication over all processors in the yz plane.
!
          if (nprocy>1.or.nprocz>1) then
            call mpiallreduce_sum(lnrmx,lnrmx_tmp,2*nghost+1,idir=23)
            call mpiallreduce_sum(cs2mx,cs2mx_tmp,idir=23)
            lnrmx=lnrmx_tmp
            cs2mx=cs2mx_tmp
          endif
!
          do i=1,nghost; lnrmx(-i)=2.*lnrmx(0)-lnrmx(i); enddo    !???
!
!  Compute x-derivative of mean lnrho
!
          dlnrmxdx= coeffs_1_x(1,1)*(lnrmx(1)-lnrmx(-1)) &
                   +coeffs_1_x(2,1)*(lnrmx(2)-lnrmx(-2)) &
                   +coeffs_1_x(3,1)*(lnrmx(3)-lnrmx(-3))
!
!  Set ghost zones such that -chi_t*rho*T*grads -hcond*gTT = Fbot, i.e.
!  enforce:
!    ds/dx = -(cp*gamma_m1*Fbot/cs2 + K*gamma_m1*glnrho)/(gamma*K+chi_t*rho)
!
          dsdx_yz=(-cp*gamma_m1*Fbot/cs2mx)/ &
                  (chit_prof1*chi_t*exp(lnrmx(0)) + hcondxbot*gamma) - gamma_m1/gamma*dlnrmxdx
!
          if (lreference_state) dsdx_yz = dsdx_yz - reference_state(XBOT,iref_gs)
!
          do i=1,nghost
            f(l1-i,:,:,iss)=f(l1+i,:,:,iss)-dx2_bound(-i)*dsdx_yz
          enddo
        endif
!
!  top boundary
!  ============
!
      case(TOP)
!
         call not_implemented("bc_ss_flux_condturb_mean_x","for the top boundary")
!
!  capture undefined entries
!
      case default
        call fatal_error('bc_ss_flux_condturb_mean_x','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_flux_condturb_mean_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_z(f,topbot)
!
!   Constant conductive + turbulent flux through the surface
!
!   15-jul-2014/pete: adapted from bc_ss_flux_condturb_x
!    4-jun-2015/MR: added missing factor cp in RB
!                   added branches for Kramers heat conductivity
!
      use DensityMethods, only: getdlnrho_z
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: dsdz_xy, TT_xy, rho_xy
      integer :: i
!
      if (ldebug) print*,'bc_ss_flux_turb: ENTER - cs20=',cs20
!
      select case (topbot)
!
!  Check whether we want to do top or bottom (this is precessor dependent)
!
!  bottom boundary
!  ===============
!
      case(BOT)
!
! Do the pretend_lnTT=T case first
!
        if (pretend_lnTT) then
           call not_implemented("bc_ss_flux_condturb_z","for pretend_lnTT=T")
        else
!
!  Set ghost zones such that -chi_t*rho*T*grads -hcond*gTT = Fbot
!
          call getlnrho(f(:,:,n1,ilnrho),rho_xy)             ! here rho_xy = ln(rho)
          TT_xy=f(:,:,n1,iss)                                ! here TT_xy = entropy
          if (lreference_state) TT_xy(l1:l2,:) = TT_xy(l1:l2,:) + spread(reference_state(:,iref_s),2,my)
          TT_xy=cs20*exp(gamma_m1*(rho_xy-lnrho0)+cv1*TT_xy) ! here TT_xy = cs2
          TT_xy=TT_xy/(cp*gamma_m1)                          ! here TT_xy temperature
!
          call getrho(f(:,:,n1,ilnrho),rho_xy)               ! here rho_xy = rho
!
          !fac=(1./60)*dz_1(n1)
          !dlnrhodz_xy=fac*(+ 45.0*(f(:,:,n1+1,ilnrho)-f(:,:,n1-1,ilnrho)) &
          !                 -  9.0*(f(:,:,n1+2,ilnrho)-f(:,:,n1-2,ilnrho)) &
          !                 +      (f(:,:,n1+3,ilnrho)-f(:,:,n1-3,ilnrho)))
!
          if (lheatc_kramers) then
            dsdz_xy=gamma1*(Fbot/hcond0_kramers)*rho_xy**(2.*nkramers)/TT_xy**(6.5*nkramers+1.)
                                     ! no turbulent diffusivity considered here!
            rho_xy=1.-gamma1         ! not efficient
          elseif (lheatc_chiconst) then
            dsdz_xy= (Fbot/TT_xy)/(rho_xy*(chit_prof1*chi_t + cp*gamma*chi))
            rho_xy = chi*gamma_m1/(chit_prof1*chi_t/cp+gamma*chi)
          else
            dsdz_xy= (Fbot/TT_xy)/(chit_prof1*chi_t*rho_xy + hcondzbot*gamma)
            rho_xy = hcondzbot*gamma_m1/(chit_prof1*chi_t*rho_xy+gamma*hcondzbot)
          endif
!
!  Enforce ds/dz = -(cp*gamma_m1*Fbot/cs2 + K*gamma_m1*glnrho)/(gamma*K+chi_t*rho)
!
          do i=1,nghost
            call getdlnrho_z(f(:,:,:,ilnrho),n1,i,TT_xy)             ! here TT_xy = d_z ln(rho)
            f(:,:,n1-i,iss)=f(:,:,n1+i,iss) + cp*(rho_xy*TT_xy+dz2_bound(-i)*dsdz_xy)
          enddo
!
        endif
!
!  top boundary
!  ============
!
      case(TOP)
!
         call not_implemented("bc_ss_flux_condturb_z","for top boundary")
!
!  capture undefined entries
!
      case default
        call fatal_error('bc_ss_flux_condturb_z','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_flux_condturb_z
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  23-jun-2003/tony: implemented for leos_fixed_ionization
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,1),size(f,2)) :: tmp_xy
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp_old: ENTER - cs20=',cs20
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for entropy.
!  This assumes that the density is already set (ie density must register
!  first!)
!  tmp_xy = s(x,y) on the boundary.
!  gamma*s/cp = [ln(cs2/cs20)-(gamma-1)ln(rho/rho0)]
!
      if ((bcz12(ilnrho,topbot) /= 'a2') .and. (bcz12(ilnrho,topbot) /= 'a3')) &
        call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3.')
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (ldebug) print*, 'bc_ss_temp_old: set bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*,'bc_ss_temp_old: cannot have cs2bot = ', cs2bot, ' <= 0'
!
        call getlnrho(f(:,:,n1,ilnrho),tmp_xy)             ! here tmp_xy = log(rho)
        tmp_xy = (-gamma_m1*(tmp_xy-lnrho0) + log(cs2bot/cs20)) / gamma
!
        if (lreference_state) &
          tmp_xy(l1:l2,:) = tmp_xy(l1:l2,:) - spread(reference_state(:,iref_s),2,my)
        f(:,:,n1,iss) = tmp_xy
!
        do i=1,nghost
          f(:,:,n1-i,iss) = 2*tmp_xy - f(:,:,n1+i,iss)     ! reference_state?
        enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, 'bc_ss_temp_old: set top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, 'bc_ss_temp_old: cannot have cs2top = ',cs2top, ' <= 0'
!
  !     if (bcz12(ilnrho,1) /= 'a2') &
  !          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 4.')
        call getlnrho(f(:,:,n2,ilnrho),tmp_xy)
        tmp_xy = (-gamma_m1*(tmp_xy-lnrho0) + log(cs2top/cs20)) / gamma
!
        if (lreference_state) &
          tmp_xy(l1:l2,:) = tmp_xy(l1:l2,:) - spread(reference_state(:,iref_s),2,my)
        f(:,:,n2,iss) = tmp_xy
!
        do i=1,nghost
          f(:,:,n2+i,iss) = 2*tmp_xy - f(:,:,n2-i,iss)     ! reference_state?
        enddo
!
      case default
        call fatal_error('bc_ss_temp_old','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_temp_old
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso_energ(f,topbot)
!
!  Boundary condition for density *and* entropy.
!
!  This sets
!    \partial_{z} \ln\rho
!  such that
!    \partial_{z} p = \rho g_{z},
!  i.e. it enforces hydrostatic equlibrium at the boundary.
!
!  Currently this is only correct if
!    \partial_{z} lnT = 0
!  at the boundary.
!
!  12-Juil-2006/dintrans: coded
!
      use EquationOfState, only: eoscalc
      use Gravity, only: potential, gravz
      use Sub, only: div
!
      real, dimension (:,:,:,:), intent (inout) :: f
      integer, intent(IN) :: topbot
!
      real, dimension (size(f,1),size(f,2)) :: cs2
      real, dimension (l2-l1+1) :: divu
      real :: rho,ss,dlnrhodz, dssdz, cs2_point
      real :: potp,potm
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary
!
      case(BOT)
!
!  The following might work for anelastic
!
        if (ldensity) then
          if (bcz12(iss,1)/='hs') &
            call fatal_error("bc_lnrho_hds_z_iso", &
              "This boundary condition for density is "// &
              "currently only correct for bcz1(iss)='hs'")

          if (bcz12(ilnrho,1)/='nil') &
            call fatal_error("bc_lnrho_hds_z_iso", &
              "To avoid a double computation, this boundary condition "// &
              "should be used only with bcz1(ilnrho)='nil' for density")

!
          rho=getrho_s(f(l1,m1,n1,ilnrho),l1)
          ss=f(l1,m1,n1,iss)
          if (lreference_state) ss=ss+reference_state(XBOT,iref_s)
          call eoscalc(irho_ss,rho,ss, cs2=cs2_point)
!
          dlnrhodz = gamma *gravz/cs2_point
          if (ldensity_nolog) dlnrhodz=dlnrhodz*rho    ! now dlnrhodz = d rho/d z
!
          dssdz = -gamma_m1*gravz/cs2_point
!
          do i=1,nghost
            f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) - dz2_bound(-i)*dlnrhodz
            f(:,:,n1-i,iss   ) = f(:,:,n1+i,iss   ) - dz2_bound(-i)*dssdz
          enddo
        elseif (lanelastic) then
          if (bcz12(iss_b,1)/='hs') &
            call fatal_error("bc_lnrho_hds_z_iso", &
              "This boundary condition for density is "// &
              "currently only correct for bcz1(iss)='hs'")
          
          if (bcz12(irho_b,1)/='nil') &
            call fatal_error("bc_lnrho_hds_z_iso", &
              "To avoid a double computation, this boundary condition "// &
              "should be used only with bcz1(irho_b)='nil' for density.")

          call eoscalc(ipp_ss,log(f(l1,m1,n1,irho_b)),f(l1,m1,n1,iss_b),cs2=cs2_point)
!
          dlnrhodz = gamma *gravz/cs2_point
          dssdz    = gamma_m1*gravz/cs2_point
!
          do i=1,nghost
            f(:,:,n1-i,irho_b) = f(:,:,n1+i,irho_b) - dz2_bound(-i)*dlnrhodz*f(:,:,n1+1,irho_b)
            f(:,:,n1-i,iss_b ) = f(:,:,n1+i,iss_b ) - dz2_bound(-i)*dssdz
          enddo
        else
          call not_implemented("bc_lnrho_hds_z_iso","at bottom for ldensity=F and lanelastic=F")
        endif
!
!  Top boundary
!
      case(TOP)
!
        if (ldensity) then
!
          if (bcz12(iss,2)/='hs') &
            call fatal_error("bc_lnrho_hds_z_iso", &
                "This boundary condition for density is "//&
                "currently only correct for bcz2(iss)='hs'")

          if (bcz12(ilnrho,2)/='nil') &
            call fatal_error("bc_lnrho_hds_z_iso", &
              "To avoid a double computation, this boundary condition "// &
              "should be used only with bcz2(ilnrho)='nil' for density.")

!
          rho=getrho_s(f(l2,m2,n2,ilnrho),l2)
          ss=f(l2,m2,n2,iss)
          if (lreference_state) ss=ss+reference_state(XTOP,iref_s)
          call eoscalc(irho_ss,rho,ss,cs2=cs2_point)
!
          dlnrhodz = gamma *gravz/cs2_point
          if (ldensity_nolog) dlnrhodz=dlnrhodz*rho    ! now dlnrhodz = d rho/d z
!
          dssdz    = -gamma_m1*gravz/cs2_point
!
          do i=1,nghost
            f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) + dz2_bound(i)*dlnrhodz
            f(:,:,n2+i,iss   ) = f(:,:,n2-i,iss   ) + dz2_bound(i)*dssdz
          enddo
        else
          call not_implemented("bc_lnrho_hds_z_iso","at top for ldensity=F")
        endif
!
      case default
!
      endselect
!
    endsubroutine bc_lnrho_hds_z_iso_energ
!***********************************************************************
    subroutine bc_ss_temp_x(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      real, dimension(my,mz) :: lnrho_yz
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp_x: cs20=',cs20
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (ldebug) print*, 'bc_ss_temp_x: set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) call fatal_error('bc_ss_temp_x','cannot have cs2bot<=0')
!
        if (.not. pretend_lnTT) then
          tmp = 2*cv*log(cs2bot/cs20)
          call getlnrho(f(l1,:,:,ilnrho),XBOT,lnrho_yz)
          f(l1,:,:,iss) =   0.5*tmp - (cp-cv)*(lnrho_yz - lnrho0)
!
          if (lreference_state) f(l1,:,:,iss) = f(l1,:,:,iss) - reference_state(XBOT,iref_s)
!
          if (ldensity_nolog) then
            do i=1,nghost
              if (lreference_state) then
!
! Reference state assumed symmetric about boundary point.
!
                f(l1-i,:,:,iss) = - f(l1+i,:,:,iss) + tmp - 2*reference_state(i,iref_s) &
                                  - (cp-cv)*(log( (f(l1+i,:,:,irho)+reference_state(i,iref_rho)) &
                                   *(f(l1-i,:,:,irho)+reference_state(i,iref_rho)) ) - 2*lnrho0)
              else
                f(l1-i,:,:,iss) = - f(l1+i,:,:,iss) + tmp &
                                  - (cp-cv)*(log(f(l1+i,:,:,irho)*f(l1-i,:,:,irho)) - 2*lnrho0)
              endif
            enddo
          else
            do i=1,nghost
              f(l1-i,:,:,iss) = - f(l1+i,:,:,iss) + tmp &
                                - (cp-cv)*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
            enddo
          endif
!
        elseif (pretend_lnTT) then
           f(l1,:,:,iss) = log(cs2bot/gamma_m1)
           do i=1,nghost; f(l1-i,:,:,iss)=2*f(l1,:,:,iss)-f(l1+i,:,:,iss); enddo
        endif
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, 'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) call fatal_error('bc_ss_temp_x','cannot have cs2top<=0')
!
        if (.not. pretend_lnTT) then
!
          tmp = 2*cv*log(cs2top/cs20)
          call getlnrho(f(l2,:,:,ilnrho),XTOP,lnrho_yz)
          f(l2,:,:,iss) = 0.5*tmp - (cp-cv)*(lnrho_yz - lnrho0)
!
          if (lreference_state) f(l2,:,:,iss) = f(l2,:,:,iss) - reference_state(XTOP,iref_s)
!
!  Distinguish cases for linear and logarithmic density
!
          if (ldensity_nolog) then
            do i=1,nghost
              if (lreference_state) then
!
! Reference state assumed symmetric about boundary point.
!
                f(l2+i,:,:,iss) =-f(l2-i,:,:,iss) + tmp - 2.*reference_state(nx-i,iref_s) &
                                 - (cp-cv)*(log((f(l2-i,:,:,irho)+reference_state(XTOP-i,iref_rho)) &
                                               *(f(l2+i,:,:,irho)+reference_state(XTOP-i,iref_rho)))-2*lnrho0)
              else
                f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
                                  -(cp-cv)*(log(f(l2-i,:,:,irho)*f(l2+i,:,:,irho))-2*lnrho0)
              endif
            enddo
          else
            do i=1,nghost
              f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
                                - (cp-cv)*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
            enddo
          endif
        elseif (pretend_lnTT) then
          f(l2,:,:,iss) = log(cs2top/gamma_m1)
          do i=1,nghost; f(l2+i,:,:,iss)=2*f(l2,:,:,iss)-f(l2-i,:,:,iss); enddo
        endif
!
      case default
        call fatal_error('bc_ss_temp_x','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_y(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      integer :: i
      real, dimension(mx,mz) :: lnrho_xz
!
      if (ldebug) print*,'bc_ss_temp_y: cs20=',cs20
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (ldebug) print*, &
                   'bc_ss_temp_y: set y bottom temperature - cs2bot=',cs2bot
        if (cs2bot<=0.) call fatal_error('bc_ss_temp_y','cannot have cs2bot<=0')
!
        tmp = 2*cv*log(cs2bot/cs20)
        call getlnrho(f(:,m1,:,ilnrho),lnrho_xz)
!
        f(:,m1,:,iss) = 0.5*tmp - (cp-cv)*(lnrho_xz-lnrho0)
        if (lreference_state) &
          f(l1:l2,m1,:,iss) = f(l1:l2,m1,:,iss) - spread(reference_state(:,iref_s),2,mz)
!
        do i=1,nghost
          if (ldensity_nolog) then
            if (lreference_state) then
              ! not yet fully implemented
              f(l1:l2,m1-i,:,iss) =- f(l1:l2,m1+i,:,iss) + tmp &
              - (cp-cv)*(log((f(l1:l2,m1+i,:,irho)+spread(reference_state(:,iref_rho),2,mz))* &
                             (f(l1:l2,m1-i,:,irho)+spread(reference_state(:,iref_rho),2,mz)))-2*lnrho0)
            else
              f(:,m1-i,:,iss) =- f(:,m1+i,:,iss) + tmp &
                               - (cp-cv)*(log(f(:,m1+i,:,irho)*f(:,m1-i,:,irho))-2*lnrho0)
            endif
          else
            f(:,m1-i,:,iss) =- f(:,m1+i,:,iss) + tmp &
                             - (cp-cv)*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
          endif
        enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, &
                     'bc_ss_temp_y: set y top temperature - cs2top=',cs2top
        if (cs2top<=0.) call fatal_error('bc_ss_temp_y','cannot have cs2top<=0')
!
        tmp = 2*cv*log(cs2top/cs20)
        call getlnrho(f(:,m2,:,ilnrho),lnrho_xz)
!
        f(:,m2,:,iss) = 0.5*tmp - (cp-cv)*(lnrho_xz-lnrho0)
        if (lreference_state) &
          f(l1:l2,m2,:,iss) = f(l1:l2,m2,:,iss) - spread(reference_state(:,iref_s),2,mz)
!
        do i=1,nghost
          if (ldensity_nolog) then
            if (lreference_state) then
              ! not yet fully implemented
              f(l1:l2,m2+i,:,iss) =- f(l1:l2,m2-i,:,iss) + tmp &
              - (cp-cv)*(log((f(l1:l2,m2-i,:,irho)+spread(reference_state(:,iref_rho),2,mz))* &
                             (f(l1:l2,m2+i,:,irho)+spread(reference_state(:,iref_rho),2,mz)))-2*lnrho0)
            else
              f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
                   - (cp-cv)*(log(f(:,m2-i,:,irho)*f(:,m2+i,:,irho))-2*lnrho0)
            endif
          else
            f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
                 - (cp-cv)*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
          endif
        enddo
!
      case default
        call fatal_error('bc_ss_temp_y','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot,lone_sided)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!  11-oct-2016/MR: changes for use of one-sided BC formulation (chosen by setting new optional switch lone_sided)
!
      use General, only: loptest
      use Deriv, only: set_ghosts_for_onesided_ders
!
      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
      logical, optional :: lone_sided
!
      real :: tmp
      integer :: i
      real, dimension(mx,my) :: lnrho_xy
!
      if (ldebug) print*,'bc_ss_temp_z: cs20=',cs20
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (ldebug) print*, 'bc_ss_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) call fatal_error('bc_ss_temp_z','cannot have cs2bot<=0')

        if (pretend_lnTT) then
!
          tmp = 2*cv*log(cs2bot/cs20)
          call getlnrho(f(:,:,n1,ilnrho),lnrho_xy)
          f(:,:,n1,iss) = 0.5*tmp - (cp-cv)*(lnrho_xy-lnrho0)
!
          if (lreference_state) &
            f(l1:l2,:,n1,iss) = f(l1:l2,:,n1,iss) - spread(reference_state(:,iref_s),2,my)
!
!  Distinguish cases for linear and logarithmic density
!
          if (loptest(lone_sided)) then
            call set_ghosts_for_onesided_ders(f,topbot,iss,3,.true.)
          else
            if (ldensity_nolog) then
!
              do i=1,nghost
                f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
!                   - (cp-cv)*(log(f(:,:,n1+i,irho)*f(:,:,n1-i,irho))-2*lnrho0)
!AB: this could be better
!MR: but is not equivalent  !! tb corrected
!    Why anyway different from cases below, which set *s* antisymmtric w.r.t. boundary value?
                   - 2*(cp-cv)*(lnrho_xy-lnrho0)
              enddo
            else
              do i=1,nghost
                f(:,:,n1-i,iss) = -f(:,:,n1+i,iss) + tmp &
                    - (cp-cv)*(f(:,:,n1+i,ilnrho)+f(:,:,n1-i,ilnrho)-2*lnrho0)
              enddo
            endif
          endif
!
        elseif (pretend_lnTT) then
          f(:,:,n1,iss) = log(cs2bot/gamma_m1)
          if (loptest(lone_sided)) then
            call set_ghosts_for_onesided_ders(f,topbot,iss,3,.true.)
          else
            do i=1,nghost; f(:,:,n1-i,iss)=2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
          endif
        endif
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*,'bc_ss_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) call fatal_error('bc_ss_temp_z','cannot have cs2top<=0')

!DM+PC next two lines need to be looked into.
!AB: This was implemented in revision: 17029 dhruba.mitra, but it works!
        if (lread_oldsnap) &
          cs2top=cs20*exp(gamma*f(l2,m2,n2,iss)/cp+gamma_m1*(f(l2,m2,n2,ilnrho)-lnrho0))
!
        if (.not. pretend_lnTT) then
          tmp = 2*cv*log(cs2top/cs20)
          call getlnrho(f(:,:,n2,ilnrho),lnrho_xy)
          f(:,:,n2,iss) = 0.5*tmp - (cp-cv)*(lnrho_xy-lnrho0)
          if (lreference_state) &
            f(l1:l2,:,n2,iss) = f(l1:l2,:,n2,iss) - spread(reference_state(:,iref_s),2,my)
!
          if (loptest(lone_sided)) then
            call set_ghosts_for_onesided_ders(f,topbot,iss,3,.true.)
          else
!
!  Distinguish cases for linear and logarithmic density
!
            if (ldensity_nolog) then
              do i=1,nghost
                f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
                     !- (cp-cv)*(log(f(:,:,n2-i,irho)*f(:,:,n2+i,irho))-2*lnrho0)
!AB: this could be better
                     - 2*(cp-cv)*(lnrho_xy-lnrho0)
              enddo
            else
              do i=1,nghost
                f(:,:,n2+i,iss) = -f(:,:,n2-i,iss) + tmp &
                     - (cp-cv)*(f(:,:,n2-i,ilnrho)+f(:,:,n2+i,ilnrho)-2*lnrho0)
              enddo
            endif
          endif
        elseif (pretend_lnTT) then
          f(:,:,n2,iss) = log(cs2top/gamma_m1)
          if (loptest(lone_sided)) then
            call set_ghosts_for_onesided_ders(f,topbot,iss,3,.true.)
          else
            do i=1,nghost; f(:,:,n2+i,iss)=2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
          endif
        endif
!
      case default
        call fatal_error('bc_ss_temp_z','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_lnrho_temp_z(f,topbot)
!
!  boundary condition for lnrho *and* ss: constant temperature
!
!  27-sep-2002/axel: coded
!  19-aug-2005/tobi: distributed across ionization modules
!
      use Gravity, only: gravz
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      integer :: i
      real, dimension(mx,my) :: lnrho_xy
!
      if (ldebug) print*,'bc_lnrho_temp_z: cs20=',cs20
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (ldebug) print*, 'bc_lnrho_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) call fatal_error('bc_lnrho_temp_z','cannot have cs2bot<=0')

        tmp = cv*log(cs2bot/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        call getlnrho(f(:,:,n1,ilnrho),lnrho_xy)
        f(:,:,n1,iss) = tmp - (cp-cv)*(lnrho_xy-lnrho0)
        if (lreference_state) &
          f(l1:l2,:,n1,iss) = f(l1:l2,:,n1,iss) - spread(reference_state(:,iref_s),2,my)
!
        do i=1,nghost; f(:,:,n1-i,iss) = 2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2bot
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=-gravz/cs2bot
        do i=1,nghost
          f(:,:,n1-i,ilnrho)=f(:,:,n1+i,ilnrho) + cp1*(f(:,:,n1+i,iss)-f(:,:,n1-i,iss))+dz2_bound(-i)*tmp !check
        enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, 'bc_lnrho_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) call fatal_error('bc_lnrho_temp_z','cannot have cs2top<=0')
        tmp = cv*log(cs2top/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        call getlnrho(f(:,:,n2,ilnrho),lnrho_xy)
        f(:,:,n2,iss) = tmp - (cp-cv)*(lnrho_xy-lnrho0)
        if (lreference_state) &
          f(l1:l2,:,n2,iss) = f(l1:l2,:,n2,iss) - spread(reference_state(:,iref_s),2,my)
!
        do i=1,nghost; f(:,:,n2+i,iss) = 2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2top
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=gravz/cs2top
        do i=1,nghost
          f(:,:,n2+i,ilnrho)=f(:,:,n2-i,ilnrho) + cp1*(f(:,:,n2-i,iss)-f(:,:,n2+i,iss))+dz2_bound(i)*tmp   !check
        enddo
!
      case default
        call fatal_error('bc_lnrho_temp_z','invalid value of topbot')
      endselect
!
    endsubroutine bc_lnrho_temp_z
!***********************************************************************
    subroutine bc_lnrho_pressure_z(f,topbot)
!
!  boundary condition for lnrho: constant pressure
!
!   4-apr-2003/axel: coded
!   1-may-2003/axel: added the same for top boundary
!  19-aug-2005/tobi: distributed across ionization modules
!
      use Gravity, only: lnrho_bot,lnrho_top,ss_bot,ss_top
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'bc_lnrho_pressure_z: cs20=',cs20
!
!  Constant pressure, i.e. antisymmetric
!  This assumes that the entropy is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*,'bc_lnrho_pressure_z: lnrho_top,ss_top=',lnrho_top,ss_top
!
!  fix entropy if inflow (uz>0); otherwise leave s unchanged
!  afterwards set s antisymmetrically about boundary value
!
        if (lentropy) then
!         do m=m1,m2
!         do l=l1,l2
!           if (f(l,m,n1,iuz)>=0) then
!             f(l,m,n1,iss)=ss_bot
!           else
!             f(l,m,n1,iss)=f(l,m,n1+1,iss)
!           endif
!         enddo
!         enddo
          f(:,:,n2,iss)=ss_top
          if (lreference_state) &
            f(l1:l2,:,n2,iss) = f(l1:l2,:,n2,iss) - spread(reference_state(:,iref_s),2,my)
!
          do i=1,nghost; f(:,:,n2+i,iss)=2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n2,ilnrho)=lnrho_top+cp1*(ss_top-f(:,:,n2,iss))
        else
          f(:,:,n2,ilnrho)=lnrho_top
        endif
!
!  make density antisymmetric about boundary
!  another possibility might be to enforce hydrostatics
!  ie to set dlnrho/dz=-g/cs^2, assuming zero entropy gradient
!
        do i=1,nghost
          f(:,:,n2+i,ilnrho)=2*f(:,:,n2,ilnrho)-f(:,:,n2-i,ilnrho)
        enddo
!
!  bottom boundary
!
      case(BOT)
        if (ldebug) print*,'bc_lnrho_pressure_z: lnrho_bot,ss_bot=',lnrho_bot,ss_bot
!
!  fix entropy if inflow (uz>0); otherwise leave s unchanged
!  afterwards set s antisymmetrically about boundary value
!
        if (lentropy) then
!         do m=m1,m2
!         do l=l1,l2
!           if (f(l,m,n1,iuz)>=0) then
!             f(l,m,n1,iss)=ss_bot
!           else
!             f(l,m,n1,iss)=f(l,m,n1+1,iss)
!           endif
!         enddo
!         enddo
          f(:,:,n1,iss)=ss_bot
          if (lreference_state) &
            f(l1:l2,:,n1,iss) = f(l1:l2,:,n1,iss) - spread(reference_state(:,iref_s),2,my)
!
          do i=1,nghost; f(:,:,n1-i,iss)=2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n1,ilnrho)=lnrho_bot+cp1*(ss_bot-f(:,:,n1,iss))   !! factor cp1 added, check
        else
          f(:,:,n1,ilnrho)=lnrho_bot
        endif
!
!  make density antisymmetric about boundary
!  another possibility might be to enforce hydrostatics
!  ie to set dlnrho/dz=-g/cs^2, assuming zero entropy gradient
!
        do i=1,nghost
          f(:,:,n1-i,ilnrho)=2*f(:,:,n1,ilnrho)-f(:,:,n1+i,ilnrho)
        enddo
!
      case default
        call fatal_error('bc_lnrho_pressure_z','invalid value of topbot')
      endselect
!
    endsubroutine bc_lnrho_pressure_z
!***********************************************************************
    subroutine bc_ss_temp2_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      real :: tmp
      real, dimension(mx,my) :: lnrho_xy
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp2_z: cs20=',cs20
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (ldebug) print*, 'bc_ss_temp2_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) call fatal_error('bc_ss_temp2_z','cannot have cs2bot<=0')
!
        tmp = cv*log(cs2bot/cs20)
        do i=0,nghost
          call getlnrho(f(:,:,n1-i,ilnrho),lnrho_xy)
          f(:,:,n1-i,iss) = tmp - (cp-cv)*(lnrho_xy-lnrho0)
        enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, 'bc_ss_temp2_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) call fatal_error('bc_ss_temp2_z','cannot have cs2top<=0')
!
        tmp = cv*log(cs2top/cs20)
        do i=0,nghost
          call getlnrho(f(:,:,n2+i,ilnrho),lnrho_xy)
          f(:,:,n2+i,iss) = tmp - (cp-cv)*(lnrho_xy-lnrho0)
        enddo
      case default
        call fatal_error('bc_ss_temp2_z','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_temp3_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  22-jan-2013/axel: coded to impose cs2bot and dcs2bot at bottom
!
      use Gravity, only: gravz
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      real :: tmp,dcs2bot
      integer :: i
      real, dimension(mx,my) :: lnrho_xy
!
      if (ldebug) print*,'bc_ss_temp3_z: cs20=',cs20
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
!  Not yet adapted to reference_state
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        dcs2bot=gamma*gravz/(mpoly+1.)
        if (ldebug) print*, 'bc_ss_temp3_z: set cs2bot,dcs2bot=',cs2bot,dcs2bot
        if (cs2bot<=0.) call fatal_error('bc_ss_temp3_z','cannot have cs2bot<=0')
!
        do i=0,nghost
          call getlnrho(f(:,:,n1-i,ilnrho),lnrho_xy)
          f(:,:,n1-i,iss) = cv*log((cs2bot-0.5*dz2_bound(-i)*dcs2bot)/cs20)-(cp-cv)*(lnrho_xy-lnrho0)
        enddo
!
!  top boundary
!
      case(TOP)
        if (ldebug) print*, 'bc_ss_temp3_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) call fatal_error('bc_ss_temp3_z','cannot have cs2top<=0')
!
        tmp = cv*log(cs2top/cs20)
        do i=0,nghost
          call getlnrho(f(:,:,n2+i,ilnrho),lnrho_xy)
          f(:,:,n2+i,iss) = tmp - (cp-cv)*(lnrho_xy-lnrho0)
        enddo
!
      case default
        call fatal_error('bc_ss_temp3_z','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_temp3_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      use DensityMethods, only: getdlnrho_x
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
      real, dimension(:,:), allocatable :: rho_yz,dlnrho
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
      if (lreference_state) allocate(dlnrho(my,mz),rho_yz(my,mz))
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (cs2bot<=0.) print*, 'bc_ss_stemp_x: cannot have cs2bot<=0'
!
        if (lreference_state) call getrho(f(l1,:,:,ilnrho),XBOT,rho_yz)
!
        do i=1,nghost
          if (ldensity_nolog) then
!
            if (lreference_state) then
              call getdlnrho_x(f(:,:,:,ilnrho),l1,i,rho_yz,dlnrho)     ! dlnrho = d_x ln(rho)
              f(l1-i,:,:,iss) =  f(l1+i,:,:,iss) + dx2_bound(-i)*reference_state(XBOT,iref_gs) &
                               + (cp-cv)*dlnrho
            else
              f(l1-i,:,:,iss) =  f(l1+i,:,:,iss) + (cp-cv)*(log(f(l1+i,:,:,irho)/f(l1-i,:,:,irho)))
            endif
          else
            f(l1-i,:,:,iss) =  f(l1+i,:,:,iss) + (cp-cv)*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
          endif
        enddo
!
!  top boundary
!
      case(TOP)
        if (cs2top<=0.) print*, 'bc_ss_stemp_x: cannot have cs2top<=0'
!
        if (lreference_state) call getrho(f(l2,:,:,ilnrho),XTOP,rho_yz)
!
        do i=1,nghost
          if (ldensity_nolog) then
            if (lreference_state) then
              call getdlnrho_x(f(:,:,:,ilnrho),l2,i,rho_yz,dlnrho)    ! dlnrho = d_x ln(rho)
              f(l2+i,:,:,iss) =  f(l2-i,:,:,iss) - dx2_bound(i)*reference_state(XTOP,iref_gs) &
                               - (cp-cv)*dlnrho
            else
              f(l2+i,:,:,iss) = f(l2-i,:,:,iss) + (cp-cv)*log(f(l2-i,:,:,irho)/f(l2+i,:,:,irho))
            endif
          else
            f(l2+i,:,:,iss) = f(l2-i,:,:,iss) + (cp-cv)*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
          endif
        enddo
!
      case default
        call fatal_error('bc_ss_stemp_x','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_stemp_x
!***********************************************************************
    subroutine bc_ss_stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      use DensityMethods, only: getdlnrho_y
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      integer :: i
      real, dimension(mx,mz) :: dlnrho
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (cs2bot<=0.) call fatal_error('bc_ss_stemp_y','cannot have cs2bot<=0')

        do i=1,nghost
          call getdlnrho_y(f(:,:,:,ilnrho),m1,i,dlnrho)    ! dlnrho = d_y ln(rho)
          f(:,m1-i,:,iss) = f(:,m1+i,:,iss) + (cp-cv)*dlnrho
        enddo
!
!  top boundary
!
      case(TOP)
        if (cs2top<=0.) call fatal_error('bc_ss_stemp_y','cannot have cs2top<=0')

        do i=1,nghost
          call getdlnrho_y(f(:,:,:,ilnrho),m2,i,dlnrho)    ! dlnrho = d_y ln(rho)
          f(:,m2+i,:,iss) = f(:,m2-i,:,iss) - (cp-cv)*dlnrho
        enddo
!
      case default
        call fatal_error('bc_ss_stemp_y','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      use DensityMethods, only: getdlnrho_z
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
      real, dimension(mx,my) :: dlnrho
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        if (cs2bot<=0.) call fatal_error('bc_ss_stemp_z','cannot have cs2bot<=0')

        do i=1,nghost
          call getdlnrho_z(f(:,:,:,ilnrho),n1,i,dlnrho)     ! dlnrho = d_z ln(rho)
          f(:,:,n1-i,iss) = f(:,:,n1+i,iss) + (cp-cv)*dlnrho
        enddo
!
!  top boundary
!
      case(TOP)
        if (cs2top<=0.) call fatal_error('bc_ss_stemp_z','cannot have cs2top<=0')

        do i=1,nghost
          call getdlnrho_z(f(:,:,:,ilnrho),n2,i,dlnrho)     ! dlnrho = d_z ln(rho)
          f(:,:,n2+i,iss) = f(:,:,n2-i,iss) - (cp-cv)*dlnrho
        enddo
      case default
        call fatal_error('bc_ss_stemp_z','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_a2stemp_x(f,topbot)
!
!  Boundary condition for entropy: adopt boundary value for temperature in
!  the ghost zone to handle shock profiles in interstellar with steep +ve
!  1st derivative in cooled remnant shells, followed by steep -ve 1st
!  derivative inside remnant.
!  s or a2 for temperature both unstable and unphysical as the unshocked
!  exterior ISM will be comparatively homogeneous, hence allowing the ghost
!  zone to fluctuate matching the boundary values is a reasonable approx
!  of the physical flow, whilst avoiding unphysical spikes to wreck the
!  calculation.
!
!  25-2010/fred: adapted from bc_ss_stemp_z
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
!  Uniform temperature/sound speed condition for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
        case(BOT)
          if (cs2bot<=0.) call fatal_error('bc_ss_a2stemp_x','cannot have cs2bot<=0')

          do i=1,nghost
            f(l1-i,:,:,iss) = f(l1+1-i,:,:,iss)+(cp-cv)*(f(l1+1-i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
          enddo
!
!  top boundary
!
        case(TOP)
          if (cs2top<=0.) call fatal_error('bc_ss_a2stemp_x','cannot have cs2top<=0')

          do i=1,nghost
            f(l2+i,:,:,iss) = f(l2-1+i,:,:,iss)+(cp-cv)*(f(l2-1+i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
          enddo
!
        case default
          call fatal_error('bc_ss_a2stemp_x','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_a2stemp_x
!***********************************************************************
    subroutine bc_ss_a2stemp_y(f,topbot)
!
!  Boundary condition for entropy: adopt boundary value for temperature in
!  the ghost zone to handle shock profiles in interstellar with steep +ve
!  1st derivative in cooled remnant shells, followed by steep -ve 1st
!  derivative inside remnant.
!  s or a2 for temperature both unstable and unphysical as the unshocked
!  exterior ISM will be comparatively homogeneous, hence allowing the ghost
!  zone to fluctuate matching the boundary values is a reasonable approx
!  of the physical flow, whilst avoiding unphysical spikes to wreck the
!  calculation.
!
!  25-2010/fred: adapted from bc_ss_stemp_z
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
!  Uniform temperature/sound speed condition for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
        case(BOT)
          if (cs2bot<=0.) call fatal_error('bc_ss_a2stemp_y','cannot have cs2bot<=0')

          do i=1,nghost
            f(:,m1-i,:,iss) = f(:,m1+1-i,:,iss)+(cp-cv)*(f(:,m1+1-i,:,ilnrho)-f(:,m1-i,:,ilnrho))
          enddo
!
!  top boundary
!
        case(TOP)
          if (cs2top<=0.) call fatal_error('bc_ss_a2stemp_y','cannot have cs2top<=0')
              
          do i=1,nghost
            f(:,m2+i,:,iss) = f(:,m2-1+i,:,iss)+(cp-cv)*(f(:,m2-1+i,:,ilnrho)-f(:,m2+i,:,ilnrho))
          enddo
!
        case default
          call fatal_error('bc_ss_a2stemp_y','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_a2stemp_y
!***********************************************************************
    subroutine bc_ss_a2stemp_z(f,topbot)
!
!  Boundary condition for entropy: adopt boundary value for temperature in
!  the ghost zone to handle shock profiles in interstellar with steep +ve
!  1st derivative in cooled remnant shells, followed by steep -ve 1st
!  derivative inside remnant.
!  s or a2 for temperature both unstable and unphysical as the unshocked
!  exterior ISM will be comparatively homogeneous, hence allowing the ghost
!  zone to fluctuate matching the boundary values is a reasonable approx
!  of the physical flow, whilst avoiding unphysical spikes to wreck the
!  calculation.
!
!  25-2010/fred: adapted from bc_ss_stemp_z
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
!  Uniform temperature/sound speed condition for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
        case(BOT)
          if (cs2bot<=0.) call fatal_error('bc_ss_a2stemp_z','cannot have cs2bot<=0')

          do i=1,nghost
            f(:,:,n1-i,iss) = f(:,:,n1+1-i,iss) + (cp-cv)*(f(:,:,n1+1-i,ilnrho)-f(:,:,n1-i,ilnrho))
          enddo
!
!  top boundary
!
        case(TOP)
          if (cs2top<=0.) call fatal_error('bc_ss_a2stemp_z','cannot have cs2top<=0')

          do i=1,nghost
            f(:,:,n2+i,iss) = f(:,:,n2-1+i,iss) + (cp-cv)*(f(:,:,n2-1+i,ilnrho)-f(:,:,n2+i,ilnrho))
          enddo
        case default
          call fatal_error('bc_ss_a2stemp_z','invalid value of topbot')
      endselect
!
    endsubroutine bc_ss_a2stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
!  boundary condition for entropy
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the entropy module
!  26-aug-2003/tony: distributed across ionization modules
!
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f

      real, dimension (size(f,1),size(f,2)) :: TT_2d
      integer :: i
!
!  The 'ce' boundary condition for entropy makes the energy constant at
!  the boundaries.
!  This assumes that the density is already set (ie density must register
!  first!)
!
    select case (topbot)
!
! Bottom boundary
!
    case(BOT)
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      if (ldensity_nolog) then
        TT_2d=exp(gamma_m1*log(f(:,:,n1,irho)) +cv1*f(:,:,n1,iss))
      else
        TT_2d=exp(gamma_m1*    f(:,:,n1,ilnrho)+cv1*f(:,:,n1,iss))
      endif
      do i=1,nghost
        if (ldensity_nolog) then
          f(:,:,n1-i,iss)=cv*(-gamma_m1*log(f(:,:,n1-i,irho)) +log(TT_2d))
        else
          f(:,:,n1-i,iss)=cv*(-gamma_m1*    f(:,:,n1-i,ilnrho)+log(TT_2d))
        endif 
     enddo
!
! Top boundary
!
    case(TOP)
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      if (ldensity_nolog) then
        TT_2d=exp(gamma_m1*log(f(:,:,n2,irho)) +cv1*f(:,:,n2,iss))
      else
        TT_2d=exp(gamma_m1*    f(:,:,n2,ilnrho)+cv1*f(:,:,n2,iss))
      endif
      do i=1,nghost
        if (ldensity_nolog) then
          f(:,:,n2+i,iss)=cv*(-gamma_m1*log(f(:,:,n2+i,irho)) +log(TT_2d))
        else
          f(:,:,n2+i,iss)=cv*(-gamma_m1*    f(:,:,n2+i,ilnrho)+log(TT_2d))
        endif
      enddo
    case default
      call fatal_error('bc_ss_energy','invalid of topbot')
    endselect
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_ism_energ(f,topbot,j)
!
!  30-nov-15/fred: Replaced bc_ctz and bc_cdz.
!  Apply observed scale height locally from Reynolds 1991, Manchester & Taylor
!  1981 for warm ionized gas - dominant scale height above 500 parsecs.
!  Apply constant local temperature across boundary for entropy.
!  Motivation to prevent numerical spikes in shock fronts, which cannot be
!  absorbed in only three ghost cells, but boundary thermodynamics still
!  responsive to interior dynamics.
!
      use EquationOfState, only: get_cp_cv

      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,k,stat

      real :: density_scale1, density_scale
      real, dimension (:,:), allocatable :: cv_,cp_
!
      if (j/=iss) call fatal_error('bc_ism_energ','only for entropy')

      density_scale=density_scale_factor
      density_scale1=1./density_scale
!
      if (leos_ionization.and.iyH/=0) then
        allocate(cp_(size(f,1),size(f,2)),stat=stat)
        if (stat>0) call fatal_error('bc_ism_energ','could not allocate cp_')
        allocate(cv_(size(f,1),size(f,2)),stat=stat)
        if (stat>0) call fatal_error('bc_ism_energ','could not allocate cv_')
      else
        allocate(cp_(1,1),cv_(1,1)); cp_=cp; cv_=cv
      endif
!
      select case (topbot)
!
      case(BOT)               ! bottom boundary

        if (leos_ionization.and.iyH/=0) call get_cp_cv(f,0,0,n1,cp_,cv_)

        do k=1,nghost

          if (leos_chemistry.or.iyH==0) then
            f(:,:,k,iss)=f(:,:,n1,iss) + (z(n1)-z(k))*density_scale1
          else
            if (ldensity_nolog) then
              f(:,:,n1-k,iss)=f(:,:,n1,iss)+(cp_-cv_)*(log(f(:,:,n1,irho))-log(f(:,:,n1-k,irho))) + &
                              cv_*log((z(n1)-z(n1-k))*density_scale+1.)
            else
              f(:,:,n1-k,iss)=f(:,:,n1,iss)+(cp_-cv_)*(f(:,:,n1,ilnrho)-f(:,:,n1-k,ilnrho))+ &
                              cv_*log((z(n1)-z(n1-k))*density_scale+1.)
            endif
          endif

        enddo
!
      case(TOP)               ! top boundary
!
        if (leos_ionization.and.iyH/=0) call get_cp_cv(f,0,0,n2,cp_,cv_)

        do k=1,nghost

          if (leos_chemistry.or.iyH==0) then
            f(:,:,n2+k,iss)=f(:,:,n2,iss) + (z(n2+k)-z(n2))*density_scale1
          else
            if (ldensity_nolog) then
              f(:,:,n2+k,iss)=f(:,:,n2,iss)+(cp_-cv_)*(log(f(:,:,n2,irho))-log(f(:,:,n2+k,irho)))+ &
                              cv_*log((z(n2+k)-z(n2))*density_scale+1.)
            else
              f(:,:,n2+k,iss)=f(:,:,n2,iss)+(cp_-cv_)*(f(:,:,n2,ilnrho)-f(:,:,n2+k,ilnrho))+ &
                              cv_*log((z(n2+k)-z(n2))*density_scale+1.)
            endif
          endif

        enddo
!
      case default
        print*, "bc_ism_energ topbot should be BOT or TOP"
      endselect
!
    endsubroutine bc_ism_energ
!***********************************************************************
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'energy_bcs_dummies.inc'
!***********************************************************************
  endmodule EnergyBcs
