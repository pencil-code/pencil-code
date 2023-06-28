  module EnergyBcs

    use Cdata
    use DensityMethods
    use EquationOfState, only: get_gamma_etc, cs20, lnrho0, ipp_ss, irho_ss, cs2bot, cs2top
    use Messages
!
    private
    real, dimension(:,:), pointer :: reference_state

    real :: gamma, gamma1, gamma_m1, cp, cv, cp1, cv1
    logical, pointer :: lheatc_chiconst
    real, pointer :: chi,chi_t,hcondzbot,hcondztop
    real, pointer :: hcondxbot, hcondxtop    !, cs2bot, cs2top
    real, pointer :: Fbot, Ftop, FtopKtop,FbotKbot, mpoly
!
    character (len=labellen) :: meanfield_Beq_profile
    real, pointer :: meanfield_Beq, chit_quenching, uturb
    real, dimension(:), pointer :: B_ext
!
    integer, parameter :: XBOT=1, XTOP=nx

    include 'energy_bcs.h'

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
      call get_shared_variable('Fbot',Fbot)
      call get_shared_variable('Ftop',Ftop)
      call get_shared_variable('FbotKbot',FbotKbot)
      call get_shared_variable('FtopKtop',FtopKtop)
      !call get_shared_variable('cs2bot',cs2bot)
      !call get_shared_variable('cs2top',cs2top)
      call get_shared_variable('chi',chi)
      call get_shared_variable('lheatc_chiconst',lheatc_chiconst)

      call get_shared_variable('chi_t',chi_t)
      call get_shared_variable('hcondzbot',hcondzbot)
      call get_shared_variable('hcondztop',hcondztop)
!
      call get_shared_variable('hcondxbot',hcondxbot)
      call get_shared_variable('hcondxtop',hcondxtop)

      call get_shared_variable('lheatc_chiconst',lheatc_chiconst)
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
!   3-oct-16/MR: added new optional switch lone_sided
!
      use SharedVariables,only: get_shared_variable
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_ss_flux','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_flux_turb(f,topbot)
!
!   4-may-2009/axel: dummy
!
      use General, only: keep_compiler_quiet

      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_ss_flux_turb','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb
!***********************************************************************
    subroutine bc_ss_flux_turb_x(f,topbot)
!
!   31-may-2010/pete: dummy
!
      use General, only: keep_compiler_quiet

      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_ss_flux_turb_x','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_x(f,topbot)
!
!   23-apr-2014/pete: dummy
!
      use General, only: keep_compiler_quiet

      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call not_implemented('bc_ss_flux_condturb_x','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_mean_x(f,topbot)
!
!   07-jan-2015/pete: dummy
!
      use General, only: keep_compiler_quiet

      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call not_implemented('bc_ss_flux_condturb_mean_x','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_mean_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_z(f,topbot)
!
!   15-jul-2014/pete: dummy
!
      use General, only: keep_compiler_quiet

      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call not_implemented('bc_ss_flux_condturb_z','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_z
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
      use General, only: keep_compiler_quiet

      real, dimension (:,:,:,:), intent(inout) :: f
      integer, intent(IN) :: topbot
!
      call not_implemented('bc_ss_temp_old','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_old
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

      call not_implemented('bc_ss_temp_x','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)

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

      call not_implemented('bc_ss_temp_y','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)

    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot,lone_sided)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!   3-oct-16/MR: added new optional switch lone_sided
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      logical, optional :: lone_sided

      call not_implemented('bc_ss_temp_z','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f

      call not_implemented('bc_lnrho_temp_z','in noeos')

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_temp_z
!***********************************************************************
    subroutine bc_lnrho_pressure_z(f,topbot)    ! essentially = eos_idealgas
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
!  Constant pressure, i.e. antisymmetric
!  This assumes that the entropy is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
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
          do i=1,nghost; f(:,:,n2+i,iss)=2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n2,ilnrho)=lnrho_top+ss_top-f(:,:,n2,iss)   !!! diff to eos_idealgas
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
!  top boundary
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
          do i=1,nghost; f(:,:,n1-i,iss)=2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n1,ilnrho)=lnrho_bot+ss_bot-f(:,:,n1,iss)
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
        call fatal_error('bc_lnrho_pressure_z','invalid argument')
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
      call not_implemented('bc_ss_temp2_z','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_temp3_z(f,topbot)
!
!  31-jan-2013/axel: coded to impose cs2bot and dcs2bot at bottom
!
      use General, only: keep_compiler_quiet

      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_ss_temp3_z','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_ss_stemp_x','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)

    endsubroutine bc_ss_stemp_x
!***********************************************************************
    subroutine bc_ss_stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_ss_stemp_y','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  Boundary condition for entropy: symmetric temperature.
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_ss_stemp_z','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_a2stemp_x(f,topbot)
!
!  boundary condition for entropy: asymmetric temperature vanishing 2nd deriv
!
!  22-sep-2010/fred: adapted from bc_ss_stemp_z
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
!  Symmetric temperature/sound speed for entropy.
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
            if (ldensity_nolog) then
              f(l1-i,:,:,iss) = min( &
                  2*f(l1+1-i,:,:,iss)-f(l1+2-i,:,:,iss)+gamma_m1/gamma* &
                  log(f(l1+1-i,:,:,irho)**2/f(l1+2-i,:,:,irho)/f(l1-i,:,:,irho)), &
                  f(l1+i,:,:,iss)+gamma_m1/gamma* &
                  log(f(l1+i,:,:,irho)/f(l1-i,:,:,irho)))
            else
              f(l1-i,:,:,iss) = min( &
                  2*f(l1+1-i,:,:,iss)-f(l1+2-i,:,:,iss)+gamma_m1/gamma* &
                  (2*f(l1+1-i,:,:,ilnrho)-f(l1+2-i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)), &
                  f(l1+i,:,:,iss)+gamma_m1/gamma* &
                  (f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)))
            endif
          enddo
!
!  top boundary
!
        case(TOP)
          if (cs2top<=0.) call fatal_error('bc_ss_a2stemp_x','cannot have cs2top<=0')
          do i=1,nghost
            if (ldensity_nolog) then
              f(l2+i,:,:,iss) = min( &
                  2*f(l2-1+i,:,:,iss)-f(l2+2-i,:,:,iss)+gamma_m1/gamma* &
                  log(f(l2-1+i,:,:,irho)**2/f(l2+2-i,:,:,irho)/f(l2+i,:,:,irho)), &
                  f(l2+i,:,:,iss)+gamma_m1/gamma* &
                  log(f(l2-i,:,:,irho)/f(l2+i,:,:,irho)))
            else
              f(l2+i,:,:,iss) = min( &
                  2*f(l2-1+i,:,:,iss)-f(l2+2-i,:,:,iss)+gamma_m1/gamma* &
                  (2*f(l2-1+i,:,:,ilnrho)-f(l2+2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho)), &
                  f(l2+i,:,:,iss)+gamma_m1/gamma* &
                  (f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho)))
            endif
          enddo
!
        case default
          call fatal_error('bc_ss_a2stemp_x','invalid argument')
      endselect
!
    endsubroutine bc_ss_a2stemp_x
!***********************************************************************
    subroutine bc_ss_a2stemp_y(f,topbot)
!
!  boundary condition for entropy: asymmetric temperature vanishing 2nd deriv
!
!  22-sep-2010/fred: adapted from bc_ss_stemp_y
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
!  Symmetric temperature/sound speed for entropy.
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
            if (ldensity_nolog) then
              f(:,m1-i,:,iss) = min( &
                  2*f(:,m1+1-i,:,iss)-f(:,m1+2-i,:,iss)+gamma_m1/gamma* &
                  log(f(:,m1+1-i,:,irho)**2/f(:,m1+2-i,:,irho)/f(:,m1-i,:,irho)), &
                  f(:,m1-i,:,iss)+gamma_m1/gamma* &
                  log(f(:,m1+i,:,irho)/f(:,m1-i,:,irho)))
            else
              f(:,m1-i,:,iss) = min( &
                  2*f(:,m1+1-i,:,iss)-f(:,m1+2-i,:,iss)+gamma_m1/gamma* &
                  (2*f(:,m1+1-i,:,ilnrho)-f(:,m1+2-i,:,ilnrho)-f(:,m1-i,:,ilnrho)), &
                  f(:,m1-i,:,iss)+gamma_m1/gamma* &
                  (f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho)))
            endif
          enddo
!
!  top boundary
!
        case(TOP)
          if (cs2top<=0.) call fatal_error('bc_ss_a2stemp_y','cannot have cs2top<=0')
          do i=1,nghost
            if (ldensity_nolog) then
              f(:,m2+i,:,iss) = min( &
                  2*f(:,m2-1+i,:,iss)-f(:,m2-2+i,:,iss)+gamma_m1/gamma* &
                  log(f(:,m2-1+i,:,irho)**2/f(:,m2-2+i,:,irho)/f(:,m2+i,:,irho)), &
                  f(:,m2-i,:,iss)+gamma_m1/gamma* &
                  log(f(:,m2-i,:,irho)/f(:,m2+i,:,irho)))
            else
              f(:,m2+i,:,iss) = min( &
                  2*f(:,m2-1+i,:,iss)-f(:,m2-2+i,:,iss)+gamma_m1/gamma* &
                  (2*f(:,m2-1+i,:,ilnrho)-f(:,m2-2+i,:,ilnrho)-f(:,m2+i,:,ilnrho)), &
                  f(:,m2-i,:,iss)+gamma_m1/gamma* &
                  (f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho)))
            endif
          enddo
!
        case default
          call fatal_error('bc_ss_a2stemp_y','invalid argument')
      endselect
!
    endsubroutine bc_ss_a2stemp_y
!***********************************************************************
    subroutine bc_ss_a2stemp_z(f,topbot)
!
!  boundary condition for entropy: asymmetric temperature vanishing 2nd deriv
!
!  22-sep-2010/fred: adapted from bc_ss_stemp_z
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
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
          if (cs2bot<=0.) call fatal_error('bc_ss_a2stemp_z','cannot have cs2bot<=0')
          do i=1,nghost
            if (ldensity_nolog) then
              f(:,:,n1-i,iss) = min( &
                  2*f(:,:,n1+1-i,iss)-f(:,:,n1+2-i,iss) + gamma_m1/gamma* &
                  log(f(:,:,n1+1-i,irho)**2/f(:,:,n1+2-i,irho)/f(:,:,n1-i,irho)), &
                  f(:,:,n1+i,iss)+gamma_m1/gamma* &
                  log(f(:,:,n1+i,irho)/f(:,:,n1-i,irho)))
            else
              f(:,:,n1-i,iss) = min( &
                  2*f(:,:,n1+1-i,iss)-f(:,:,n1+2-i,iss) + gamma_m1/gamma* &
                  (2*f(:,:,n1+1-i,ilnrho)-f(:,:,n1+2-i,ilnrho)-f(:,:,n1-i,ilnrho)), &
                  f(:,:,n1+i,iss)+gamma_m1/gamma* &
                  (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)))
            endif
          enddo
!
!  top boundary
!
        case(TOP)
          if (cs2top<=0.) call fatal_error('bc_ss_a2stemp_z','cannot have cs2top<=0')
          do i=1,nghost
            if (ldensity_nolog) then
              f(:,:,n2+i,iss) = min( &
                  2*f(:,:,n2-1+i,iss)-f(:,:,n2-2+i,iss) + gamma_m1/gamma* &
                  log(f(:,:,n2-1+i,irho)**2/f(:,:,n2-2+i,irho)/f(:,:,n2+i,irho)), &
                  f(:,:,n2-i,iss)+gamma_m1/gamma* &
                  log(f(:,:,n2-i,irho)/f(:,:,n2+i,irho)))
            else
              f(:,:,n2+i,iss) = min( &
                  2*f(:,:,n2-1+i,iss)-f(:,:,n2-2+i,iss) + gamma_m1/gamma* &
                  (2*f(:,:,n2-1+i,ilnrho)-f(:,:,n2-2+i,ilnrho)-f(:,:,n2+i,ilnrho)), &
                  f(:,:,n2-i,iss)+gamma_m1/gamma* &
                  (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)))
            endif
          enddo
        case default
          call fatal_error('bc_ss_a2stemp_z','invalid argument')
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
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_ss_energy','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      use General, only: keep_compiler_quiet

      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_stellar_surface','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)
!
      use General, only: keep_compiler_quiet

      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call not_implemented('bc_lnrho_cfb_r_iso','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso(f,topbot)
!
      use General, only: keep_compiler_quiet

      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call not_implemented('bc_lnrho_hds_z_iso','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hds_z_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)
!
      use General, only: keep_compiler_quiet

      real, dimension (:,:,:,:) :: f
      integer, intent(IN) :: topbot
!
      call not_implemented('bc_lnrho_hdss_z_iso','in noeos')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_stellar_surface","in noeos")

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_stellar_surface
!***********************************************************************
  endmodule EnergyBcs

