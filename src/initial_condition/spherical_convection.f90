! $Id$
!
!  Isentropic initial condition (density, entropy) for convection in
!  spherical coordinates. Produces an isentropic stratification with
!  given surface gravity and surface temperature, and a heat
!  conduction profile proportional to r^-15. The setup was originally
!  introduced in Astron. Nachr., 332, 883.
!
!  02-sep-13/pete: coded
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
  use EquationOfState
  use Sub, only: step, der_step
!
  implicit none
!
  include '../initial_condition.h'
!
  real :: star_luminosity=1.0, Rstar=1.0, Rtran=1.2, chit0=1e8
  real :: xi0=1.0, npoly1=1.5, npoly_jump=1.0, nad=1.5, rbot=1.0, wbot=0.1
  real :: npoly_fac=1.0, npoly_exp=1.0, r_ss=1.0,chiSGS_top=1.0
  real :: Fbottom, wtran=0.02,Tcor_jump=1.0, kramers_hcond0=0.0
  logical :: lcorona=.false., lwrite_cooling_profile=.false.
  logical :: lwrite_hcond_profile=.true., lkappa_constchi=.false.
  character (len=labellen) :: strat_type='polytropic'
!
  namelist /initial_condition_pars/ &
      star_luminosity, Rstar, nad, npoly1, npoly_jump, xi0, & 
      lcorona, Rtran, wtran, Tcor_jump, strat_type, r_ss, npoly_fac, &
      npoly_exp, chiSGS_top, chit0, lwrite_cooling_profile, &
      lwrite_hcond_profile, lkappa_constchi, rbot, wbot
!
  contains
!***********************************************************************
    subroutine register_initial_condition()
!
!  Register variables associated with this module; likely none.
!
!  07-may-09/wlad: coded
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_initial_condition
!***********************************************************************
    subroutine initialize_initial_condition(f)
!
!  Initialize any module variables which are parameter dependent.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_initial_condition
!***********************************************************************
    subroutine initial_condition_all(f,profiles)
!
!  Initializes all the f arrays in one call. This subroutine is called last.
!
!  02-sep-13/pete: coded
!  06-oct-13/joern: add coronal layer
!  25-jan-15/MR: extended case 'polytropic' to linear density
!  10-feb-15/MR: added optional parameter 'profiles' (intended to replace f):
!                profiles are returned instead of set into f.
!                Profiles hcond and glhcond not written then!
!
      use SharedVariables, only: get_shared_variable
      use EquationOfState, only: gamma, rho0, cs20
      use General, only: safe_character_assign
      use Mpicomm, only: stop_it, mpiallreduce_sum
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (nx,*),             optional, intent(out)  :: profiles
!
      real, dimension (mx) :: TT, TTc, rho_prof
      real, dimension (mx) :: lnrho, ss_prof, cs2_prof
      real, dimension (nxgrid) :: kappa, gkappa, npoly2, gnpoly2
      real, dimension (nxgrid) :: rho_global, TT_global, TTc_global, dTdr_global, dTdrc_global
      real, dimension (nxgrid) :: dlnTdr_global, dlnrhodr_global, lnrho_global
      real, dimension (nxgrid) :: ss_global, cs2_global, npoly_global
      real, dimension (nxgrid) :: drhodr_global, del2rho_global
      real :: T00, rho00, Rsurf, Tsurf, coef1, L00, sigma, cs2_surf, cs2_top
      real :: cs2_bot
      real :: Tcor, Rmin, wmin, cs2_cor, rho_surf
      real :: Lsun=3.84e26, Rsun=7e8, Omsun=2.7e-6, Msun=2e30, cvsun=20786.1
      real :: GG=6.67348e-11, rhosun=200., fluxratio, Omsim, gratio, rratio
      real :: T00sun=2.23e6
      real :: volume, total_mass, tmp, tau_KH
      real, pointer :: gravx, cp, cv
      integer :: i, j, n, m, ix, ierr, nsurf, nsurf_global
      integer, parameter :: unit=1

      character (len=120) :: wfile
!
!     Retrieve cp, cv, and gravx
!
      call get_shared_variable('cp', cp, caller='initial_condition_all')
      call get_shared_variable('cv', cv)
      call get_shared_variable('gravx', gravx)
!
!  Select type of stratification
!
      select case (strat_type)
!
!  Single polytrope for the convection zone
!
      case ('polytropic')
!
!  Surface of the star
!
      if (lcorona) then
         Rsurf=Rstar
         if (x0+Lxyz(1)<=Rstar) then
           write(unit=errormsg,fmt=*) &
           'initial_condition: your wedge has to have a radius bigger than Rstar'//& 
           ' for using a corona'
           call fatal_error('initial_condition',errormsg)
         endif
         Rtran=Rtran*Rstar
         wtran=wtran*Rstar
         Rmin=Rsurf+(Rtran-Rsurf)/6.
         wmin=wtran/2.0
         if (lroot) &
           print*,'initial_condition: you are using a coronal envelope'
         do i=1,nx
           if (x(l1+i)>=Rsurf) then
             nsurf=i-1
             exit
           endif
         enddo
         do j=1, nxgrid
           if (xglobal(nghost+j)>=Rsurf) then
             nsurf_global=j
             exit
           endif
         enddo
      else
         Rsurf=x0+Lxyz(1)
         nsurf=l2-l1
         nsurf_global=nxgrid
      endif
!
!  Temperature using a constant polytropic index npoly1
!
      TT=1.
      TT(l1:l2)=gravx/(cv*(gamma-1.))*(xi0/Rstar + 1./(npoly1+1.)*(1./x(l1:l2) - 1./Rsurf))
      TT_global=gravx/(cv*(gamma-1.))*(xi0/Rstar + 1./(npoly1+1.)*(1./xglobal(nghost+1:nxgrid+nghost) - 1./Rsurf))
      dTdr_global=-gravx/xglobal(nghost+1:nxgrid+nghost)**2./(cv*(gamma-1)*(npoly1+1.))

      T00=gravx/(cv*(gamma-1.))*(xi0/Rstar + 1./(npoly1+1.)*(1./x0 - 1./Rsurf))
      Tsurf=gravx/(cv*(gamma-1.))*xi0/Rstar
!
      if (lcorona) then
!
!  Using a step function for the temperature profile in the corona,
!  and another step function to make a smooth transition.
!
        Tcor=Tcor_jump*T00
        TTc=1.
        TTc(l1+nsurf:mx)=Tsurf+(Tcor-Tsurf)*step(x(l1+nsurf:mx),Rtran,wtran)
        TT(l1+nsurf:mx)=TT(l1+nsurf:mx)+(TTc(l1+nsurf:mx)-TT(l1+nsurf:mx))*step(x(l1+nsurf:mx), Rmin, wmin)
!
!  global temperature
!
        TTc_global(nsurf_global:nxgrid)=Tsurf+(Tcor-Tsurf)*step(xglobal(nghost+nsurf_global:nxgrid+nghost),Rtran,wtran)
        TT_global(nsurf_global:nxgrid)=TT_global(nsurf_global:nxgrid) + & 
                                      (TTc_global(nsurf_global:nxgrid)-TT_global(nsurf_global:nxgrid))* & 
                                      step(xglobal(nghost+nsurf_global:nxgrid+nghost), Rmin, wmin)
!
!  global temperature derivative
!
        dTdrc_global=0
        dTdrc_global(nsurf_global:nxgrid)=(Tcor-Tsurf)*der_step(xglobal(nghost+nsurf_global:nxgrid+nghost),Rtran,wtran)
        dTdr_global(nsurf_global:nxgrid)=dTdr_global(nsurf_global:nxgrid)+(dTdrc_global(nsurf_global:nxgrid) - & 
             dTdr_global(nsurf_global:nxgrid))*step(xglobal(nghost+nsurf_global:nxgrid+nghost), Rmin, wmin) + &
             (TTc_global(nsurf_global:nxgrid)-TT_global(nsurf_global:nxgrid)) * &
             der_step(xglobal(nghost+nsurf_global:nxgrid+nghost),Rmin, wmin)
        dlnTdr_global=dTdr_global/TT_global
      endif
!
!  Density stratification assuming an isentropic atmosphere with ss=0. 
!
      rho_prof=rho0*(TT/T00)**(1./(gamma-1.))
      rho00=rho0*(T00/T00)**(1./(gamma-1.))
      rho_surf=rho0*(Tsurf/T00)**(1./(gamma-1.))
      rho_global=rho0*(TT_global/T00)**(1./(gamma-1.))
!
      drhodr_global=(1./(gamma-1.))*rho_global*dTdr_global/TT_global
      del2rho_global=(1./(gamma-1.))*(drhodr_global*dTdr_global/TT_global - &
           rho_global*(dTdr_global/TT_global)**2. + &
           rho_global*2.*gravx/(xglobal(nghost+1:nxgrid+nghost)**3.)/ &
           (cv*(gamma-1.)*(npoly1+1.)*TT_global)) + &
           (2./xglobal(nghost+1:nxgrid+nghost))*drhodr_global
!
      lnrho=log(rho_prof/rho0)
      lnrho_global=log(rho_global/rho0)
!
      if (lcorona .or. lkappa_constchi) then
        dlnrhodr_global=-dlnTdr_global-gravx/xglobal(nghost+1:nxgrid+nghost)**2/(cv*(gamma-1)*TT_global)
        if (lcorona) then
          do j=nsurf_global-nsurf_global/10, nxgrid
            lnrho_global(j)=lnrho_global(j-1)+dlnrhodr_global(j-1)*(xglobal(nghost+j)-xglobal(nghost+j-1))
          enddo
          lnrho(l1:l2)=lnrho_global(ipx*nx+1:(ipx+1)*nx)
        endif
      else
        dlnrhodr_global=0.
      endif
!
!  Renormalize entropy with rho0 and cs20
!
      cs2_prof=cs20*TT*cv*gamma*(gamma-1.)
      cs2_global=cs20*TT_global*cv*gamma*(gamma-1.)
      ss_prof=log(cs2_prof/cs20)/gamma - & 
              (gamma-1.)/(gamma)*(lnrho-log(rho0))
      ss_global=log(cs2_global/cs20)/gamma - & 
              (gamma-1.)/(gamma)*(lnrho_global-log(rho0))
!
!  Put lnrho and ss into the f-array
!
      if (present(profiles)) then
!
         profiles(:,iref_rho  ) = rho_prof(l1:l2)
         profiles(:,iref_grho ) = drhodr_global (ipx*nx+1:(ipx+1)*nx)
         profiles(:,iref_d2rho) = del2rho_global(ipx*nx+1:(ipx+1)*nx)
         profiles(:,iref_s    ) = ss_global(1)
!
         if (lroot) then 
           call safe_character_assign(wfile,'reference_state.dat')
           open(unit,file=wfile,status='unknown')
           do ix=1,nxgrid
              write(unit,'(4(2x,1pe15.8))') rho_global(ix),drhodr_global(ix),del2rho_global(ix),ss_global(1)
           enddo
           close(unit)
         endif
!
      elseif (lreference_state) then
         f(l1:l2,:,:,irho) = 0.
         f(l1:l2,:,:,iss) = 0.
      else
!
        do m=m1,m2
        do n=n1,n2
          if (ldensity_nolog) then
            f(l1:l2,m,n,irho) = exp(lnrho(l1:l2))
          else
            f(l1:l2,m,n,ilnrho) = lnrho(l1:l2)
          endif
          f(l1:l2,m,n,iss) = ss_prof(l1:l2)
        enddo
        enddo
!
      endif
!
!  Compute hcond and glhcond using global x-array
!
      coef1=star_luminosity*rho0*sqrt(gravx*Rstar)*cv*(gamma-1.)/(4.*pi)
!
      npoly2=0.
      gnpoly2=0.
!
      do n=1,nxgrid
!
!     set kappa that chi is constant with initial stratification
!
        if (lkappa_constchi) then
          npoly2(n)=nad*exp(lnrho_global(n)-lnrho_global(1))
          gnpoly2(n)=npoly2(n)*dlnrhodr_global(n)
        else
          npoly2(n)=npoly_jump*(xglobal(nghost+n)/x0)**(-15.)+nad-npoly_jump
          gnpoly2(n)=15./xglobal(nghost+n)*(nad-npoly_jump-npoly2(n))
          if ((xglobal(nghost+n)>=Rstar) .and. &
              ((npoly2(n)+1.)/exp(lnrho_global(n))>2.*(npoly2(1)+1.)/exp(lnrho_global(1)))) then
            npoly2(n)=2.*(npoly2(1)+1)*exp(lnrho_global(n)-lnrho_global(1))-1.
            gnpoly2(n)=2.*(npoly2(1)+1)*exp(lnrho_global(n)-lnrho_global(1))*dlnrhodr_global(n)
          endif
        endif
      enddo
!
      kappa=coef1*(npoly2+1.)
      gkappa=coef1*gnpoly2
!
!  Thermal equilibrium condition
!
      case ('thermal_equilibrium')
!
!  Compute hcond and glhcond using global x-array
!
        coef1=star_luminosity*rho0*sqrt(gravx*Rstar)*cv*(gamma-1.)/(4.*pi)
!
        do n=1,nxgrid
          npoly2(n)=npoly_jump*(xglobal(nghost+n)/x0)**(-15.)+nad-npoly_jump
          gnpoly2(n)=15./xglobal(nghost+n)*(nad-npoly_jump-npoly2(n))
        enddo
!
        kappa=coef1*(npoly2+1.)
        gkappa=coef1*gnpoly2
!
      case ('piecewise-poly')
!
!  Piecewise polytropic thermodynamic stratification (hcond not computed!)
!
      Rsurf=x0+Lxyz(1)
      Tsurf=gravx/(cv*(gamma-1.))*xi0/Rstar
!
!  Compute depth dependent 'polytropic index' npoly
!
      npoly_global  = npoly1 - (npoly_jump)*(step(xglobal(nghost+1:nxgrid+nghost), rbot, wbot))
!
!  Temperature gradient using a constant 'polytropic index' npoly
!
      dTdr_global=-gravx/xglobal(nghost+1:nxgrid+nghost)**2./(cv*(gamma-1)*(npoly_global+1.))
!
!  Integrate temperature from the known surface value to the interior
!
      TT_global(nxgrid)=Tsurf
      do j=1,nxgrid-1
        TT_global(nxgrid-j)=TT_global(nxgrid-j+1)-dTdr_global(nxgrid-j+1)*(xglobal(nxgrid+nghost-j+1)-xglobal(nxgrid+nghost-j))
      enddo
      TT(l1:l2)=TT_global(ipx*nx+1:(ipx+1)*nx)
      T00=TT_global(1)
      dlnTdr_global=dTdr_global/TT_global
!
!  Density gradient assuming hydrostatic equilibrium
!
      dlnrhodr_global=-dlnTdr_global-gravx/xglobal(nghost+1:nxgrid+nghost)**2/(cv*(gamma-1)*TT_global)

      lnrho_global(1)=log(10*rho0) ! Dangerous hard-coded value for time being
      do j=2,nxgrid
        lnrho_global(j)=lnrho_global(j-1)+dlnrhodr_global(j-1)*(xglobal(nghost+j)-xglobal(nghost+j-1))
      enddo
!
!  Renormalize density such that rho=rho0 at r=rbot
!
      do j=1,nxgrid
         if (xglobal(nghost+j) < rbot) then
           del2rho_global(j)=lnrho_global(j) ! del2rho_global is not lnrho_global(r < rbot)
         else
           del2rho_global(j)=maxval(lnrho_global)
         endif
      enddo
      lnrho_global=lnrho_global-minval(del2rho_global)
      rho_global=exp(lnrho_global)
      lnrho(l1:l2)=lnrho_global(ipx*nx+1:(ipx+1)*nx)
      rho00=rho0
      rho_surf=exp(lnrho_global(nxgrid))
!
!  Renormalize entropy with rho0 and cs20
!
      cs2_prof=cs20*TT*cv*gamma*(gamma-1.)
      cs2_global=cs20*TT_global*cv*gamma*(gamma-1.)
      ss_prof=log(cs2_prof/cs20)/gamma - &
              (gamma-1.)/(gamma)*(lnrho-log(rho0))
      ss_global=log(cs2_global/cs20)/gamma - &
              (gamma-1.)/(gamma)*(lnrho_global-log(rho0))
!
!  Put lnrho and ss into the f-array
!
      do m=m1,m2
      do n=n1,n2
        if (ldensity_nolog) then
          f(l1:l2,m,n,irho) = exp(lnrho(l1:l2))
        else
          f(l1:l2,m,n,ilnrho) = lnrho(l1:l2)
        endif
        f(l1:l2,m,n,iss) = ss_prof(l1:l2)
      enddo
      enddo
!
!  Fail if reference state is used
!
      if (lreference_state) then
         write(unit=errormsg,fmt=*) &
           'initial_condition: piecewise-poly initial condition not'//&
           ' implemented with reference state'
         call fatal_error('initial_condition',errormsg)
      endif
!
      endselect
!
!  Write kappa and gkappa to file to be read by run.in
!
      if (lwrite_hcond_profile) then
        if (lroot) then 
          call safe_character_assign(wfile,'hcond_glhc.dat')
          open(unit,file=wfile,status='unknown')
          do ix=1,nxgrid
            write(unit,'(2(2x,1pe12.5))') kappa(ix),gkappa(ix)
          enddo
          close(unit)
        endif
      endif
!
!  Write cs2 to a file to be used in cooling
!
      if (lwrite_cooling_profile) then
        if (lroot) then 
          call safe_character_assign(wfile,'cooling_profile.dat')
          open(unit,file=wfile,status='unknown')
          do ix=1,nxgrid
            write(unit,'(3(2x,1pe17.10))') cs2_global(ix)
          enddo
          close(unit)
        endif
      endif
!
!  Compute flux at the bottom and a modified Stefan-Boltzmann constant
!  assuming the the total flux at the outer boundary radiates through
!  the surface for the Fgs boundary condition.
!
      L00=star_luminosity*rho0*gravx**1.5*sqrt(Rstar)
      Fbottom=L00/(4.*pi*x0**2)
      sigma=(L00/(4.*pi*Rsurf**2))/Tsurf**4
      cs2_bot=cs20*T00*cv*gamma*(gamma-1.)
      cs2_top=cs20*Tsurf*cv*gamma*(gamma-1.)
      if (lcorona) then
        cs2_top=cs20*TT_global(nxgrid)*cv*gamma*(gamma-1.)
        cs2_surf=cs20*Tsurf*cv*gamma*(gamma-1.)
        cs2_cor=cs20*Tcor*cv*gamma*(gamma-1.)
      endif

!
!  Compute the ratio of the dimensionless luminosity in the simulation
!  versus the Sun in order to determine a rotation rate for the
!  simulation which reproduces the rotational effect in the Sun.
!
      fluxratio=star_luminosity*rhosun*(GG*Msun)**1.5*sqrt(Rsun)/Lsun
      gratio=gravx/(GG*Msun)
      rratio=Rsun/Rstar
      Omsim=fluxratio**(1./3.)*sqrt(gratio)*rratio**(1.5)*Omsun
      chiSGS_top=sqrt(gratio)*fluxratio**(1./3.)*chit0/sqrt(rratio)
!
!  Compute hcond0_kramers for runs with Kramers opacity
!
      kramers_hcond0=-Fbottom/(rho_global(1)**(-2)*TT_global(1)**6.5*dTdr_global(1))
!
!  Compute total mass.
!
      total_mass=0.
!
      do n=n1,n2
        tmp=0.
        do m=m1,m2
          tmp=tmp+sum(exp(lnrho(l1:l2))*dVol_x(l1:l2))*dVol_y(m)
        enddo
        total_mass=total_mass+tmp*dVol_z(n)
      enddo
!
      if (ncpus>1) then
        call mpiallreduce_sum(total_mass,tmp)
        total_mass=tmp
      endif
!
      volume=((x0+Lxyz(1))**3-x0**3)*(cos(y0)-cos(y0+Lxyz(2)))*((z0+Lxyz(3))-z0)/3.
!
!  Compute Kelvin-Helmholtz time scale: GM*M_domain/(2*Rsurf*Luminosity)
!
!      tau_KH=gravx*total_mass/L00/(Rsurf-x0)/2.
! PJK: Assuming the convection zone as a point mass at r=R now.
! PJK: The previous formulation is effectively the same but putting
! PJK: the CZ as a point mass at r=0.3R.
      tau_KH=gravx*total_mass/(2.*Rsurf*L00)
!
      if (lroot) then
         print*,''
         print*,'initial_condition: Fbottom    =',Fbottom
         print*,'initial_condition: SigmaSBt   =',sigma
         print*,'initial_condition: cs2bot     =',cs2_bot
         print*,'initial_condition: cs2top     =',cs2_top
         print*,'initial_condition: fluxratio  =',fluxratio
         print*,'initial_condition: Omsim      =',Omsim
         print*,'initial_condition: chiSGS_top =',chiSGS_top
         print*,'initial_condition: gratio     =',gratio
         print*,'initial_condition: rratio     =',rratio
         print*,'initial_condition: volume     =',volume
         print*,'initial_condition: total_mass =',total_mass
         print*,'initial_condition: number of density scale heights =',lnrho_global(1)-lnrho_global(nxgrid)
         print*,'initial_condition: hcond0_kramers =',kramers_hcond0
         print*,'initial_condition: Kelvin-Helmholtz time in PC units =', tau_KH
         print*,'initial_condition: Kelvin-Helmholtz time in years =', tau_KH*Omsim/Omsun/3.15360e+07
         if (lcorona) then
           print*, ''
           print*,'initial_condition: rcool      =',Rsurf+(Rtran-Rsurf)/6.
           print*,'initial_condition: wcool      =',wmin/1.5
           print*,'initial_condition: cs2cool    =',cs2_surf*0.85
           print*,'initial_condition: rcool2     =',Rtran
           print*,'initial_condition: wcool2     =',wtran
           print*,'initial_condition: cs2cool2   =',cs2_cor
         endif
         print*,''
!
!  Compute temperatures at the surface and in the corona assuming solar
!  temperature at the base of the convection zone.
!
         print*,'initial_condition: Temperature at the surface (of domain)     =',Tsurf*T00sun/T00, 'K'
         print*,'initial_condition: Temperature at the bottom                  =',T00sun, 'K'
         if (lcorona) then
           print*,'initial_condition: Temperature in the corona                  =',Tcor*T00sun/T00, 'K'
         endif
         print*,'initial_condition: Density stratification in convection zone  =',rho00/rho_surf
         if (lcorona) then
           print*,'initial_condition: Density stratification in the corona       =',&
                exp(lnrho_global(nsurf_global)-lnrho_global(nxgrid))
           print*,'initial_condition: Density stratification with corona         =',exp(log(rho00)-lnrho_global(nxgrid))
         endif
         print*,'initial_condition: Turbulent heat conductivity at the surface =',chit0, 'm^2/s'
         print*,''
         if (.not. lequidist(1)) then
           print*, 'initial_condition: using a nonequidistant grid in x with dx_bot/dx_top =',&
                    dx1grid(nxgrid)/dx1grid(1)
           print*, ''
         endif
      endif
!
    endsubroutine initial_condition_all
!***********************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uu
!***********************************************************************
    subroutine initial_condition_lnrho(f)
!
!  Initialize logarithmic density. init_lnrho will take care of
!  converting it to linear density if you use ldensity_nolog.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=initial_condition_pars, IOSTAT=iostat)
!
    endsubroutine read_initial_condition_pars
!***********************************************************************
    subroutine write_initial_condition_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=initial_condition_pars)
!
    endsubroutine write_initial_condition_pars
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from noinitial_condition.f90 for any    **
!**  InitialCondition routines not implemented in this file        **
!**                                                                **
    include '../initial_condition_dummies.inc'
!********************************************************************
endmodule InitialCondition
