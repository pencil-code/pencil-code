! $Id$
!
!  This module includes physics related to planet atmospheres.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! variables global to this module
!
  real, dimension(my,mz) :: mu_ss=0.
  real, dimension(my) :: lat  ! latitude in [rad]
  real, dimension(mz) :: lon  ! longitude in [rad]
  real, dimension(mx,my,mz) :: rr,rr1,siny,cosy  !  r,1/r,sin(th) and cos(th)
  real, dimension(mx,my,mz,3) :: Bext=0.,Jext=0.  !  time-dependent external fields
  real, dimension(mx) :: fact_near_topbot=0.
  real, dimension(my) :: fact_near_polar=0.
!  constants for unit conversion
  real :: gamma=1.
  real :: r2m=1., rho2kg_m3=1., u2m_s=1., cp2si=1.
  real :: pp2Pa=1., TT2K=1., tt2s=1., g2m3_s2=1., eta2si=1.
!
! variables in the temperature reference profile
!
  real :: dlog10p_T_ref, log10p_T_ref_min, log10p_T_ref_max
  real, dimension(:), allocatable :: tau_rad,logp_ref,dTeq,Teq_night
  integer :: n_T_ref
!
! variables in the magnetic diffusivity reference profile
!
  real :: ref_eta_dlog10p,ref_eta_log10p_min,ref_eta_log10p_max
  real :: ref_eta_dT
  real, dimension(:), allocatable :: ref_eta_log10p,ref_eta_T,ref_etaP,ref_etaPT,ref_log10etaPT
  integer :: ref_eta_nP, ref_eta_nT
!
! Init parameters
!
  !  reference values for unit conversion
  real :: R_planet=9.5e7   !  unit: [m]
  real :: rho_ref=1.e-2    !  unit: [kg/m3]
  real :: cs_ref=2.e3      !  unit: [m/s]
  real :: cp_ref=1.44e4    !  unit: [J/(kg*K)]
  real :: T_ref=1533       !  unit: [K]
  !
  real :: lon_ss=0., lat_ss=0.            ! unit: [degree]
  real :: dTeqbot=0., dTeqtop=100.        ! unit: [K]
  real :: peqtop=1.d2, peqbot=1.d6        ! unit: [Pa]
  real :: tauradtop=1.d4, tauradbot=1.d7  ! unit: [s]
  real :: pradtop=1.d3, pradbot=1.d6      ! unit: [Pa]
  real :: pbot0=1.e7                      ! unit: [Pa]
  real :: q_drag=0.0
  real :: q_sponge_r=0.0
  real :: q_sponge_polar=0.0
  real :: frac_sponge_r=0
  real :: frac_sponge_polar=0
  !
  logical :: linit_equilibrium=.false.
  logical :: leta_planet_as_aux=.false.
!
! Run parameters
!
  real :: tau_slow_heating=-1.,t0_slow_heating=0.,dTeq_max=1000.
  real :: Bext_ampl=0.
  real :: eta_floor=-1., eta_ceiling=-1.
  character (len=labellen) :: iBext='nothing',ieta_PT='nothing'
!
!
!
  namelist /special_init_pars/ &
      R_planet,rho_ref,cs_ref,cp_ref,T_ref,pbot0,&
      lon_ss,lat_ss,peqtop,peqbot,tauradtop,tauradbot,&
      pradtop,pradbot,dTeqbot,dTeqtop,linit_equilibrium,&
      Bext_ampl,iBext,frac_sponge_r,frac_sponge_polar,q_sponge_polar,&
      q_drag,q_sponge_r,leta_planet_as_aux
!
  namelist /special_run_pars/ &
      tau_slow_heating,t0_slow_heating,Bext_ampl,iBext,frac_sponge_r,&
      q_drag,q_sponge_r,&
      dTeq_max,dTeqtop,dTeqbot,ieta_PT,eta_floor,eta_ceiling,leta_planet_as_aux
!
!
! Declare index of new variables in f array (if any).
!
!!   integer :: ispecial=0
!!   integer :: ispecaux=0
!
!! Diagnostic variables (needs to be consistent with reset list below).
!
!!   integer :: idiag_POSSIBLEDIAGNOSTIC=0
!
  contains
!****************************************************************************
    subroutine initialize_special(f)
!
      use Sub, only: register_report_aux
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (.not.ltemperature_nolog) call fatal_error('initialize_special', &
          'special/planet_atmosphere is formulated in TT only')
!
      if (leta_planet_as_aux) call register_report_aux('eta_planet', ieta_planet)
!
!  convert y and z to latitude and longitude in rad
!
      lat=0.5*pi-y
      lon=z-pi
!
!  3d coordinates for convenience
!
      rr = spread(spread(x,2,my),3,mz)
      rr1 = 1./rr
      siny = spread(spread(sin(y),1,mx),3,mz)
      cosy = spread(spread(cos(y),1,mx),3,mz)
!
!  unit conversion
!
      call prepare_unit_conversion
!
!  read in the reference P-T profile, and calculate Teq and tau_rad etc.
!
      call prepare_Tref_and_tau
!
!  read in the reference eta(P,T) profile
!
      if (lmagnetic) then
        if (ieta_PT=='eta_P') call read_eta_P_file
        if (ieta_PT=='eta_PT') call read_eta_PT_file
      endif
!
!  define sponge layer factors
!
      if (frac_sponge_r>0.) then
!  bottom:
        where ( (x-xyz0(1) <= frac_sponge_r*Lxyz(1)) ) &
            fact_near_topbot = 0.5 + 0.5*cos(pi*(x-xyz0(1))/(frac_sponge_r*Lxyz(1)))
!  top:
        where ( (xyz1(1)-x <= frac_sponge_r*Lxyz(1)) ) &
            fact_near_topbot = 0.5 + 0.5*cos(pi*(xyz1(1)-x)/(frac_sponge_r*Lxyz(1)))
      endif
      !
      if (frac_sponge_polar>0.) then
!  north:
        where ( (y-xyz0(2) <= frac_sponge_polar*Lxyz(2)) ) &
            fact_near_polar = 0.5 + 0.5*cos(pi*(y-xyz0(2))/(frac_sponge_polar*Lxyz(2)))
!  south:
        where ( (xyz1(2)-y <= frac_sponge_polar*Lxyz(2)) ) &
            fact_near_polar = 0.5 + 0.5*cos(pi*(xyz1(2)-y)/(frac_sponge_polar*Lxyz(2)))
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialize with self-consistent equilibrium state
!  23-oct-23/hongzhe: coded
!
      use Gravity, only: g0
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      real :: p_tmp,dp,rhoeq_tmp,Teq_tmp,tau_tmp
      integer, parameter :: nsub=24
      integer :: i,j,k,isub
!
      if (linit_equilibrium) then
!
!  Compute the initial equilibrium state for pressure and density.
!  dP = -rho_eq(P) * g(r) * dr. Integrate from bottom.
!
        call get_mu_ss(mu_ss,lon_ss,lat_ss)
        do j=1,my
        do k=1,mz
          p_tmp = pbot0  !  in [Pa]
          do i=1,(ipx+1)*nx
            do isub=1,nsub  !  use small length step, dx/nsub
              call calc_Teq_tau_pmn(Teq_tmp,tau_tmp,p_tmp,j,k,1.) ! Teq in [K]
              rhoeq_tmp = p_tmp/Teq_tmp/(cp_ref*(gamma-1.)/gamma) / rho2kg_m3  !  in code unit
              dp = -rhoeq_tmp * g0/xglobal(i+nghost)**2 * dx/nsub * pp2Pa  ! in [Pa]
              p_tmp = p_tmp+dp
              if (p_tmp<0.) call fatal_error('init_special', &
                  'failed to compute initial state, probably because of too low Nx')
            enddo
!
            if (i>=ipx*nx+1) then
              f(i-ipx*nx+nghost,j,k,ilnrho) = log(rhoeq_tmp)
              f(i-ipx*nx+nghost,j,k,iTT)    =  Teq_tmp/TT2K
            endif
          enddo
        enddo
        enddo
      endif
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
    lpenc_requested(i_rho)=.true.
    lpenc_requested(i_pp)=.true.
    lpenc_requested(i_uu)=.true.
!
    if (lmagnetic .and. iBext/='nothing') then
      lpenc_requested(i_rho1)=.true.
      lpenc_requested(i_bb)=.true.
      lpenc_requested(i_jj)=.true.
    endif
!
    if (lmagnetic .and. ieta_PT/='nothing') then
      lpenc_requested(i_rho1)=.true.
      lpenc_requested(i_cv1)=.true.
      lpenc_requested(i_TT)=.true.
      lpenc_requested(i_j2)=.true.
      lpenc_requested(i_del2a)=.true.
    endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!  Add damping layers near the top and bottom boundaries. Ref: Dowling (1998).
!
!  24-nov-23/kuan,hongzhe: coded
!  21-feb-24/kuan: revised sponge layer
!
      use Sub, only: cross_mn,multsv_mn
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx,3) :: jtot,btot,jxbtot,jxbtotr
      integer :: i, j
!
!  Add top or bottom sponge layer
!
      do i=iux,iuz
        df(l1:l2,m,n,i) = df(l1:l2,m,n,i) - fact_near_topbot(l1:l2) * q_sponge_r * f(l1:l2,m,n,i)
      enddo
!
!  Add polar sponge layer of ur
!
      df(l1:l2,m,n,iux) = df (l1:l2,m,n,iux) - fact_near_polar(m) * q_sponge_polar * f(l1:l2,m,n,iux)
!
!  Add velocity drag, dudt = ... - q_drag * u
!
      df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - q_drag * f(l1:l2,m,n,iux:iuz)
!
!  Add Lorentz force from the background field
!
      if (lmagnetic) then
        jtot = p%jj + Jext(l1:l2,m,n,:)
        btot = p%bb + Bext(l1:l2,m,n,:)
        call cross_mn(jtot,btot,jxbtot)
        call multsv_mn(p%rho1,jxbtot,jxbtotr)
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) + jxbtotr
      endif
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
!  08-sep-23/hongzhe: specific things for planet atmospheres.
!                     Reference: Rogers & Komacek (2014)
!  27-sep-23/hongzhe: outsourced from temperature_idealgas.f90
!  26-feb-24/kuan: Possibility of slowly turning on the heating term
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(nx) :: Teq_x,tau_rad_x
      real :: f_slow_heating
!
!  Calculate the local equilibrium temperature Teq_x,
!  and the local radiative cooling time tau_rad_x,
!  for all l at m,n, given the local pressure.
!  Teq_x is in [K] and tau_rad_x is in [s].
!
      if (tau_slow_heating>0) then
        f_slow_heating = min(1.0d0*dTeq_max/dTeqtop,1.0+(t-t0_slow_heating)/tau_slow_heating*(dTeq_max-dTeqtop)/dTeqtop)
      else
        f_slow_heating = 1.
      endif
      if (t<t0_slow_heating) f_slow_heating=1.
!
      call calc_Teq_tau_mn(Teq_x,tau_rad_x,p%pp*pp2Pa,m,n,f_slow_heating)
!
      df(l1:l2,m,n,iTT) = df(l1:l2,m,n,iTT) - (p%TT-Teq_x/TT2K)/(tau_rad_x/tt2s)
!
    endsubroutine special_calc_energy
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!  09-oct-23/hongzhe: add a background dipole field
!
      use Sub, only: cross_mn
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx,3) :: uxb_ext
      real, dimension (nx) :: eta_x
!
      call cross_mn(p%uu,Bext(l1:l2,m,n,:),uxb_ext)
      df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + uxb_ext
!
!  add customized eta profile
!
      select case (ieta_PT)
      case ('nothing')
        ! do nothing
      case ('eta_P')
        call interpolation1d(ref_eta_log10p,log10(ref_etaP), &
                             ref_eta_nP,ref_eta_dlog10p, &
                             log10(p%pp*pp2Pa),eta_x )
        eta_x = 10.**eta_x
      case ('eta_PT')
        call interpolation2d(ref_eta_log10p,ref_eta_T,ref_log10etaPT, &
                             ref_eta_nP,ref_eta_nT,ref_eta_dlog10p,ref_eta_dT, &
                             log10(p%pp),p%TT,eta_x )
        eta_x = 10.**eta_x
        if (eta_floor>0.) then
          where (eta_x<eta_floor) eta_x = eta_floor
        endif
      !
        if (eta_ceiling>0.) then
          where (eta_x>eta_ceiling) eta_x = eta_ceiling
        endif
      case ('Perna+2010')
        !  reference: Perna+2010, Eq. (1) and Menou+2012 Sec. 3.3
        eta_x = 214.545 * exp(25188./p%TT/TT2K) * &
            sqrt(cp_ref*p%rho*rho2kg_m3) * (p%TT*TT2K)**(-0.25) / eta2si
        where (eta_x > eta_ceiling)
          eta_x = eta_ceiling
        elsewhere (eta_x < eta_floor)
          eta_x = eta_floor
        endwhere
      case default
        call fatal_error('special_calc_magnetic','no such ieta_PT')
      endselect
!
      if (leta_planet_as_aux) f(l1:l2,m,n,ieta_planet) = eta_x
!
! Apply customized eta profile to the induction and heat equations
!
      if (ieta_PT/='nothing') then
        df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + spread(eta_x,2,3) * p%del2a
        df(l1:l2,m,n,iTT)   = df(l1:l2,m,n,iTT) + p%cv1*p%rho1*mu0 * eta_x * p%j2
        if (lfirst.and.ldt) maxdiffus=max(maxdiffus,eta_x*dxyz_2)
      endif
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
!  compute cos(angle between the substellar point)
!  could be time-dependent
!
      call get_mu_ss(mu_ss,lon_ss,lat_ss)
!
      call calc_Bext
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine prepare_unit_conversion
!
!  Calculate constants that convert code units to SI.
!
!  28-sep-23/hongzhe: coded
!
      use EquationOfState, only: rho0,cs0,get_gamma_etc
!
      real :: cp
!
      if (unit_system/='SI') call fatal_error('prepare_unit_conversion','please use SI system')
!
      call get_gamma_etc(gamma,cp=cp)
!
      r2m = R_planet/xyz0(1)         !  length to [m]
      rho2kg_m3 = rho_ref/rho0       !  density to [kg/m3]
      u2m_s = cs_ref/cs0             !  velocisty to [m/s]
      cp2si = cp_ref/cp              !  cp to [J/(kg*K)]
!
      pp2Pa = rho2kg_m3 * u2m_s**2.  !  pressure to [Pa]
      TT2K = u2m_s**2. / cp2si       !  temperature to [K]
      tt2s = r2m / u2m_s             !  time to [s]
      g2m3_s2 = r2m * u2m_s**2.      !  GM to [m3/s2]
      eta2si = r2m**2./tt2s          !  eta to m^2/s
!
      if (lroot) then
        print*,'Constants for uit conversion: gamma,r2m,rho2kg_m3,u2m_s,cp2si= ', &
                gamma,r2m,rho2kg_m3,u2m_s,cp2si
        print*,'Constants for uit conversion: pp2Pa,TT2K,tt2s, g2m3_s2,eta2si= ', &
                pp2Pa,TT2K,tt2s,g2m3_s2,eta2si
      endif
!
    endsubroutine  prepare_unit_conversion
!***********************************************************************
    subroutine prepare_Tref_and_tau
!
!   Read the reference profile and calculated tau_rad.
!   All quantities in this subroutine are in SI units.
!
!   28-sep-23/xianyu,hongzhe: coded
!
      real, dimension(:), allocatable :: p_temp_ref,temp_ref
      real :: alpha
      integer :: i
      logical :: lTref_file_exists
!
!  read in Tref, in physical unit.
!
!!HZ: need use a less specific file name
      inquire(FILE='iro-teq-tint100K-regrid-Pa.txt', EXIST=lTref_file_exists)
      if (.not.lTref_file_exists) call fatal_error('initialize_special', &
          'Must provide a Tref file')
!
      open(1,file='iro-teq-tint100K-regrid-Pa.txt')
      read(1,*)
      read(1,*) dlog10p_T_ref, log10p_T_ref_min, log10p_T_ref_max
      read(1,*) n_T_ref
      if(allocated(logp_ref)) deallocate(logp_ref)
      allocate(logp_ref(n_T_ref),temp_ref(n_T_ref))
      read(1,*) logp_ref,temp_ref
      if (lroot) then
        print*, 'n_T_ref=',n_T_ref
        print *,'Here is the baseline radiative equil T profile'
        print *,'used in the Newtonian relaxation:'
        print *,'p [Pa] and T_eq[K]:'
        do i=1,n_T_ref
          print*, 'logp_ref,temp_ref=',logp_ref(i),temp_ref(i)
        enddo
      endif
      close(1)
!
!  convert units, and calculate tau_rad, dTeq, and Teq_night
!
      if(allocated(tau_rad))    deallocate(tau_rad)
      if(allocated(dTeq))       deallocate(dTeq)
      if(allocated(Teq_night))  deallocate(Teq_night)
      allocate( p_temp_ref(n_T_ref), tau_rad(n_T_ref), dTeq(n_T_ref),Teq_night(n_T_ref) )
      p_temp_ref = 10.**logp_ref
  !
      where (p_temp_ref>=peqtop .and. p_temp_ref<=peqbot)
        dTeq = dTeqbot + (dTeqtop - dTeqbot) * log(p_temp_ref/peqbot)/log(peqtop/peqbot)
      elsewhere (p_temp_ref < peqtop)
        dTeq = dTeqtop
      elsewhere
        dTeq = dTeqbot
      endwhere
!
      alpha=log(tauradtop/tauradbot)/log(pradtop/pradbot)
      where (p_temp_ref>=pradtop .and. p_temp_ref<=pradbot)
        tau_rad = tauradbot * (p_temp_ref/pradbot)**alpha
      elsewhere (p_temp_ref < pradtop)
        tau_rad = tauradtop
      elsewhere
        tau_rad = tauradbot
      endwhere
!
      Teq_night = temp_ref - 0.5*dTeq
!
      deallocate(temp_ref,p_temp_ref)
!
!  for debug purpose, output Teq at day- and night-points
!
      if (lroot) then
        open(1,file=trim(datadir)//'/Teq_night.dat',position='append')
        write(1,*) Teq_night
        write(1,*) logp_ref
        close(1)
!
        open(1,file=trim(datadir)//'/Teq_day.dat',position='append')
        write(1,*) Teq_night+dTeq
        write(1,*) logp_ref
        close(1)
      endif
!
    endsubroutine  prepare_Tref_and_tau
!***********************************************************************
    subroutine get_mu_ss(mu_ss,lonss,latss)
!
!  Given the substellar longitude lonss and latitude latss (both in
!  degree), compute the cosine of the angle of the sun from overhead.
!  This will be used in future RT calculation.
!  It also tells you whether you're on the dayside, since
!  mu_ss>0 for dayside, mu_ss=0 on the terminator, and mu_ss<0 on the nightside.
!
!  28-sep-23/xianyu,hongzhe: coded
!
      use General, only: itoa
!
      real, dimension(my,mz), intent(out) :: mu_ss
      real, intent (in) :: lonss, latss
      integer :: j,k
      character (len=1) :: chproc
!
      real, PARAMETER :: deg2rad=pi/180.
!
      mu_ss = spread(cos(lon),1,my) * spread(cos(lat),2,mz) &
             * cos(lonss*deg2rad) * cos(latss*deg2rad) + &
              spread(sin(lon),1,my)*spread(cos(lat),2,mz) &
             * sin(lonss*deg2rad) * cos(latss*deg2rad) + &
              spread(sin(lat),2,mz)*sin(latss*deg2rad)
!
!  debug
!
!      chproc=itoa(iproc)
!      open(1,file=trim(datadir)//'/mu_proc'//chproc//'.dat',position='append')
!      do j=m1,m2
!      do k=n1,n2
!        write(1,*) y(j),z(k),mu_ss(j,k)
!      enddo
!      enddo
!      close(1)      
!
    endsubroutine  get_mu_ss
!***********************************************************************
    subroutine calc_Bext
!
!  Calculate the external B field originated from the planet interior.
!  This contributes an extra uxb term in dA/dt.
!
!  13-nov-23/hongzhe: coded
!
      real :: r0,r1,th0,th1
      select case (iBext)
      case ('nothing')
        Bext = 0.
        Jext = 0.
      case ('dipole')
        ! dipole = mu0/(4pi) * ( 3*rhat*(rhat dot m)-m ) / |r|^3
        !        = mu0/(4pi) * m * { 2*costh/r^3, sinth/r^3, 0 }
        Bext(:,:,:,1) = Bext_ampl * 2.*cosy * (xyz0(1)*rr1)**3.
        Bext(:,:,:,2) = Bext_ampl * siny    * (xyz0(1)*rr1)**3.
        Bext(:,:,:,3) = 0.
        Jext = 0.
      case ('dipole2')
        ! mimic initaa='dipole' to account for b.c.
        ! A_phi ~ (r-r0)(r-r1)(sin(theta)-sin(theta0))(sin(theta)-sin(theta1))
        r0 = xyz0(1)
        r1 = xyz1(1)
        th0 = xyz0(2)
        th1 = xyz1(2)
        Bext(:,:,:,1) = Bext_ampl * cosy*(r0-rr)*(r1-rr)*(siny*(3.*siny-2.*sin(th1))+sin(th0)*(-2.*siny+sin(th1)))/rr/siny
        Bext(:,:,:,2) = Bext_ampl * (r0*r1-2.*(r0+r1)*rr+3.*rr**2)*(-siny+sin(th0))*(siny-sin(th1))/rr
        Bext(:,:,:,3) = 0.
        call fatal_error('calc_Bext','Jext not implemented for iBext=dipole2')
      case default
        call fatal_error('calc_Bext','no such iBext')
      endselect
!
    endsubroutine  calc_Bext
!***********************************************************************
    subroutine calc_Teq_tau_pmn(T_local,tau_local,press,m,n,f_slow)
!
!  Given the local pressure p and position m,n, calculate the equilibrium
!  temperature T_local and the radiative cooling time tau_local, by
!  interpolating Tref. Reference: Komacek+Showman2016.
!  The output T_local is in [K] and tau_local is in [s], and
!  the input press is in [Pa].
!
!  23-oct-23/hongzhe: outsourced from special_calc_energy
!  26-feb-24/kuan: added slow heating term
!
      real, intent(out) :: T_local,tau_local
      real, intent(in) :: press,f_slow
      integer, intent(in) :: m,n
!
      real :: log10p_local,Teq_local1,Teq_local2
      integer :: ip
!
!  Index of the logp_ref that is just smaller than log10(pressure)
!
      log10p_local = log10(press)
      ip = 1+floor((log10p_local-log10p_T_ref_min)/dlog10p_T_ref)
!
!  Interpolation for T_local and tau_local
!
      if (ip>=n_T_ref) then
        T_local = (Teq_night(n_T_ref)-0.5*(f_slow-1.)*dTeq(n_T_ref)) + dTeq(n_T_ref)*max(0.,mu_ss(m,n))*f_slow
        tau_local = tau_rad(n_T_ref)
      elseif (ip<=1) then
        T_local = (Teq_night(1)-0.5*(f_slow-1.)*dTeq(1)) + dTeq(1)*max(0.,mu_ss(m,n))*f_slow
        tau_local = tau_rad(1)
      else
!
!  The two closest values of equilibrium T given press,m,n
!
        Teq_local1 = (Teq_night(ip)  -0.5*(f_slow-1.)*dTeq(ip))   + dTeq(ip)  *max(0.,mu_ss(m,n))*f_slow
        Teq_local2 = (Teq_night(ip+1)-0.5*(f_slow-1.)*dTeq(ip+1)) + dTeq(ip+1)*max(0.,mu_ss(m,n))*f_slow
        T_local = Teq_local1+(Teq_local2-Teq_local1)*   &
                  (log10p_local-logp_ref(ip))/   &
                  (logp_ref(ip+1)-logp_ref(ip))   ! unit of K
!
        tau_local = log10(tau_rad(ip))+  &
                    (log10(tau_rad(ip+1))-log10(tau_rad(ip)))*   &
                    (log10p_local-logp_ref(ip))/   &
                    (logp_ref(ip+1)-logp_ref(ip))
        tau_local = 10.**tau_local  ! unit of s
      endif
!
    endsubroutine calc_Teq_tau_pmn
!***********************************************************************
    subroutine calc_Teq_tau_mn(T_local,tau_local,press,m,n,f_slow)
!
!  Same as calc_Teq_tau_pmn, but for an array of pressure values
!
!  23-oct-23/hongzhe: coded
!
    real, dimension(nx), intent(out) :: T_local,tau_local
    real, dimension(nx), intent(in) :: press
    real, intent(in) :: f_slow
    integer, intent(in) :: m,n
!
    integer :: i
!
    do i=1,nx
      call calc_Teq_tau_pmn( T_local(i),tau_local(i), press(i),m,n,f_slow)
    enddo
!
    endsubroutine calc_Teq_tau_mn
!***********************************************************************
    subroutine read_eta_P_file
!
!   Read the reference eta profile.
!   All quantities in this subroutine are in SI units.
!
!   13-mar-24/hongzhe: coded
!
      integer :: i
      logical :: leta_file_exists
!
!  read in eta, in SI unit
!
      inquire(FILE='eta_P.txt', EXIST=leta_file_exists)
      if (.not.leta_file_exists) call fatal_error('read_eta_P_file', &
          'Must provide an eta(P) file')
!
      open(1,file='eta_P.txt')
      read(1,*)
      read(1,*) ref_eta_dlog10p, ref_eta_log10p_min, ref_eta_log10p_max
      read(1,*) ref_eta_nP
      if(allocated(ref_eta_log10p)) deallocate(ref_eta_log10p)
      if(allocated(ref_etaP)) deallocate(ref_etaP)
      allocate(ref_eta_log10p(ref_eta_nP),ref_etaP(ref_eta_nP))
      read(1,*) ref_eta_log10p,ref_etaP
      if (lroot) then
        print*, 'ref_eta_nP=',ref_eta_nP
        print *,'Here is the eta(P) profile:'
        print *,'p [Pa] and eta[m^2/s]:'
        do i=1,ref_eta_nP
          print*, 'ref_eta_log10p,ref_etaP=',ref_eta_log10p(i),ref_etaP(i)
        enddo
      endif
      close(1)
!
!  convert to code unit
!
      ref_etaP=ref_etaP/eta2si
!
!  for debug purpose, output the result
!
      if (lroot) then
        open(1,file=trim(datadir)//'/eta_P_normalized.dat',status='replace')
        write(1,*) ref_eta_log10p
        write(1,*) ref_etaP
        close(1)
      endif
!
    endsubroutine  read_eta_P_file
!***********************************************************************
    subroutine read_eta_PT_file
!
!   Read the reference eta profile.
!   All quantities in this subroutine are in SI units.
!
!   2-jul-24/hongzhe: coded
!
      integer :: i,j
      logical :: leta_file_exists
!
!  read in eta(P,T) in SI unit and convert to code unit
!
      inquire(FILE='eta_PT.txt', EXIST=leta_file_exists)
      if (.not.leta_file_exists) call fatal_error('read_eta_PT_file', &
          'Must provide an eta(P,T) file')
!
      open(1,file='eta_PT.txt')
      read(1,*) ! skip the comment line
      ! pressure
      read(1,*) ref_eta_nP
      if (allocated(ref_eta_log10p)) deallocate(ref_eta_log10p)
      allocate(ref_eta_log10p(ref_eta_nP))
      read(1,*) ref_eta_log10p
      ref_eta_log10p = log10(ref_eta_log10p/pp2Pa)
      ref_eta_dlog10p = (ref_eta_log10p(ref_eta_nP) - ref_eta_log10p(1)) / (ref_eta_nP-1.)
      ! temperature
      read(1,*) ref_eta_nT
      if (allocated(ref_eta_T)) deallocate(ref_eta_T)
      allocate(ref_eta_T(ref_eta_nT))
      read(1,*) ref_eta_T
      ref_eta_T = ref_eta_T/TT2K
      ref_eta_dT = (ref_eta_T(ref_eta_nT)-ref_eta_T(1)) / (ref_eta_nT-1.)
      ! eta, reshaping into 1d array
      if (allocated(ref_etaPT)) deallocate(ref_etaPT)
      allocate(ref_etaPT(ref_eta_nP*ref_eta_nT))
      read(1,*) ref_etaPT
      ref_etaPT=ref_etaPT/eta2si
!
      close(1)
!
!  store in log10
!
      if (allocated(ref_log10etaPT)) deallocate(ref_log10etaPT)
      allocate(ref_log10etaPT(ref_eta_nP*ref_eta_nT))
      ref_log10etaPT = log10(ref_etaPT)
!
!  for debug purpose, output the result
!
      if (lroot) then
        open(1,file=trim(datadir)//'/eta_PT_normalized.dat',status='replace')
        write(1,*) ref_eta_log10p
        write(1,*) ref_eta_T
        write(1,*) reshape(ref_etaPT,(/ref_eta_nP,ref_eta_nT/))
        close(1)
      endif
!
    endsubroutine  read_eta_PT_file
!***********************************************************************
    subroutine interpolation1d(xref,yref,nref,dxref,xlocal,ylocal)
!
!  Linear interpolation from the reference table (xref,yref).
!
      real, dimension(:), intent(in) :: xref,yref
      integer, intent(in) :: nref
      real, intent(in) :: dxref
      real, dimension(nx), intent(in) :: xlocal
      real, dimension(nx), intent(out) :: ylocal
!
      real, dimension(nx) :: xlocal2,x1,x2,y1,y2
      integer, dimension(nx) :: ix
!
      xlocal2 = min(max(xlocal,xref(1)),xref(nref))
!
!  Find the index so that xref(ix) is just smaller than xlocal
!
      ix = 1+floor((xlocal2-xref(1))/dxref)
!
!  Linear interpolation
!
      x1 = xref(ix)
      x2 = xref(ix+1)
      y1 = yref(ix)
      y2 = yref(ix+1)
      ylocal = (-x2*y1+x1*y2)/(x1-x2) + (y1-y2)/(x1-x2)*xlocal2
!
    endsubroutine interpolation1d
!***********************************************************************
    subroutine interpolation2d(xref,yref,zref,nxref,nyref,dxref,dyref,xlocal,ylocal,zlocal)
!
!  Linear interpolation from the reference table (xref,yref,zref).
!
      real, dimension(:), intent(in) :: xref,yref,zref
      integer, intent(in) :: nxref,nyref
      real, intent(in) :: dxref,dyref
      real, dimension(nx), intent(in) :: xlocal,ylocal
      real, dimension(nx), intent(out) :: zlocal
!
      real, dimension(nx) :: x1,x2,y1,y2,z11,z12,z21,z22
      real, dimension(nx) :: xlocal2,ylocal2,fact,a0,a1,a2,a3
      integer, dimension(nx) :: ix,iy
!
      xlocal2 = min(max(xlocal,xref(1)),xref(nxref))
      ylocal2 = min(max(ylocal,yref(1)),yref(nyref))
!
!  Find the index so that xref(ix) is just smaller than xlocal,
!  and similarly to y
!
      ix = 1+floor((xlocal2-xref(1))/dxref)
      iy = 1+floor((ylocal2-yref(1))/dyref)
!
!  Polynomial interpolation: z=a0+a1*x+a2*y+a3*x*y
!
      x1 = xref(ix)
      x2 = xref(ix+1)
      y1 = yref(iy)
      y2 = yref(iy+1)
      z11 = zref((iy-1)*nxref + ix)
      z12 = zref((iy  )*nxref + ix)
      z21 = zref((iy-1)*nxref + ix+1)
      z22 = zref((iy  )*nxref + ix+1)
      !
      fact = 1./(x1-x2)/(y1-y2)
      a0 = (x2*y2*z11-x2*y1*z12-x1*y2*z21+x1*y1*z22)*fact
      a1 = (y2*(z21-z11)+y1*(z12-z22))*fact
      a2 = (x2*(z12-z11)+x1*(z21-z22))*fact
      a3 = (z11-z12-z21+z22)*fact
      !
      zlocal = a0 + a1*xlocal2 + a2*ylocal2 + a3*xlocal2*ylocal2
!
    endsubroutine interpolation2d
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copied dummy routines from nospecial.f90 for any              **
!**  routine not implemented in this file.                         **
!**                                                                **
    include '../special_dummies.inc'
!*********************************************************************** 
endmodule Special
