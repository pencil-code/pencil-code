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
  use Cparam
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
  real, dimension(mx,my,mz,3) :: Bext=0.  !  time-dependent external field
!  constants for unit conversion
  real :: gamma=1.
  real :: r2m=1., rho2kg_m3=1., u2m_s=1., cp2si=1.
  real :: pp2Pa=1., TT2K=1., tt2s=1., g2m3_s2=1.
!
! variables in the temperature reference profile
!
  real :: dlog10p_T_ref, log10p_T_ref_min, log10p_T_ref_max
  real, dimension(:), allocatable :: tau_rad,logp_ref,dTeq,Teq_night
  integer :: n_T_ref
!
! variables in the magnetic diffusivity reference profile
!
  real :: dlog10p_eta_ref, log10p_eta_ref_min, log10p_eta_ref_max
  real, dimension(:), allocatable :: log10p_eta_ref, eta_ref
  integer :: n_eta_ref
!
! Init parameters
!
  !  reference values for unit conversion
  real :: R_planet=9.5e7   !  unit: [m]
  real :: rho_ref=1.e-2    !  unit: [kg/m3]
  real :: cs_ref=2.e3      !  unit: [m/s]
  real :: cp_ref=1.44e4    !  unit: [J/(kg*K)]
  real :: T_ref=1533       !  unit: [K]
  real :: eta0_ref=1.12314e6  !  [m^2/s], eta at T_ref
  !
  real :: lon_ss=0., lat_ss=0.            ! unit: [degree]
  real :: dTeqbot=0., dTeqtop=100.        ! unit: [K]
  real :: peqtop=1.d2, peqbot=1.d6        ! unit: [Pa]
  real :: tauradtop=1.d4, tauradbot=1.d7  ! unit: [s]
  real :: pradtop=1.d3, pradbot=1.d6      ! unit: [Pa]
  real :: pbot0=1.e7                      ! unit: [Pa]
  real :: q_drag=0.0
  real :: q_sponge=0.0
  !
  integer :: n_sponge=0
  !
  logical :: linit_equilibrium=.false.
  logical :: lsponge_top=.false.,lsponge_bottom=.false.,lvelocity_drag=.false.,lsponge_dt=.false.
!
! Run parameters
!
  real :: tau_slow_heating=-1.,t0_slow_heating=0.,dTeq_max=1000.
  real :: Bext_ampl=0.,f_eta=0.
  character (len=labellen) :: iBext='nothing',ietaPT='nothing',ieta_order='1'
!
!
!
  namelist /special_init_pars/ &
      R_planet,rho_ref,cs_ref,cp_ref,T_ref,pbot0,&
      lon_ss,lat_ss,peqtop,peqbot,tauradtop,tauradbot,&
      pradtop,pradbot,dTeqbot,dTeqtop,linit_equilibrium,&
      Bext_ampl,iBext,n_sponge,lsponge_top,lsponge_bottom,&
      lvelocity_drag,q_drag,q_sponge,lsponge_dt,eta0_ref
!
  namelist /special_run_pars/ &
      tau_slow_heating,t0_slow_heating,Bext_ampl,iBext,n_sponge,&
      lsponge_top,lsponge_bottom,lvelocity_drag,q_drag,q_sponge,lsponge_dt,&
      dTeq_max,dTeqtop,dTeqbot,ietaPT,f_eta,ieta_order
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
      real, dimension (mx,my,mz,mfarray) :: f
!
      if (.not.ltemperature_nolog) call fatal_error('initialize_special', &
          'special/planet_atmosphere is formulated in TT only')
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
        select case (ietaPT)
        case ('nothing')
          !  do nothing
        case ('eta_P')
          call prepare_eta_P
        case default
          call fatal_error('initialize_special','no such ietaPT')
        endselect
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
    if (lmagnetic .and. ietaPT/='nothing') then
      select case (ieta_order)
      case ('1')
        lpenc_requested(i_del2a)=.true.
      case ('3')
        lpenc_requested(i_del6a)=.true.
      case default
        call fatal_error('pencil_criteria_special','no such ieta_order')
      endselect
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
!  ksp in his work is equal to n_sponge.
!  q(j) is correspond to mu0(n_sponge+1-j) in eq(57).
!
!  24-nov-23/kuan,hongzhe: coded
!  21-feb-24/kuan: revised sponge layer
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension(:), allocatable :: q
      integer :: i, j
!
      if (n_sponge>0 .and. (lsponge_top.or.lsponge_bottom)) then
        allocate(q(n_sponge))
        q=0.
      endif
!
!  Add top or bottom sponge layer
!
      if (lsponge_top .and. llast_proc_x .and. it>1) then
        do j=1,n_sponge
          q=q_sponge*(1.-cos(pi*j/n_sponge))
        enddo
        if (lsponge_dt) q=q/dt
        !
        do i=iux,iuz
          df(l2-n_sponge+1:l2,m,n,i) = df(l2-n_sponge+1:l2,m,n,i) - q * f(l2-n_sponge+1:l2,m,n,i)
        enddo
      endif
!
      if (lsponge_bottom .and. lfirst_proc_x .and. it>1) then
        do j=1,n_sponge
          q=q_sponge*(1.-cos(pi*(n_sponge-j+1)/n_sponge))
        enddo
        if (lsponge_dt) q=q/dt
        !
        do i=iux,iuz
          df(l1:l1+n_sponge-1,m,n,i) = df (l1:l1+n_sponge-1,m,n,i) - q * f(l1:l1+n_sponge-1,m,n,i)
        enddo
      endif
!
! Add velocity drag, dudt = ... - q_drag * u
!
      if (lvelocity_drag) &
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz) - q_drag * f(l1:l2,m,n,iux:iuz)
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
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!
      real, dimension (nx,3) :: uxb_ext,fres
      real, dimension (nx) :: eta_x
!
      call cross_mn(p%uu,Bext(l1:l2,m,n,:),uxb_ext)
      df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + uxb_ext
!
!  add customized eta profile
!
      select case (ietaPT)
      case ('nothing')
        ! do nothing
      case ('eta_P')
        call calc_eta_p(eta_x,p%pp*pp2Pa)
      case default
        call fatal_error('special_calc_magnetic','no such ietaPT')
      endselect
!
      if (ietaPT/='nothing') then
        select case (ieta_order)
        case ('1')
          df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + &
              f_eta * p%del2a * spread(eta_x,2,3)
        case ('3')
          df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + &
              f_eta * p%del6a * spread(eta_x,2,3)
        case default
          call fatal_error('special_calc_magnetic','no such ieta_order')
        endselect
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
!
      if (lroot) then
        print*,'Constants for uit conversion: gamma,r2m,rho2kg_m3,u2m_s,cp2si = ', &
                gamma,r2m,rho2kg_m3,u2m_s,cp2si
        print*,'Constants for uit conversion: pp2Pa,TT2K,tt2s, g2m3_s2 = ', &
                pp2Pa,TT2K,tt2s,g2m3_s2
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
!
      select case (iBext)
      case ('nothing')
        Bext = 0.
      case ('dipole')
        ! dipole = mu0/(4pi) * ( 3*rhat*(rhat dot m)-m ) / |r|^3
        !        = mu0/(4pi) * m * { 2*costh/r^3, sinth/r^3, 0 }
        Bext(:,:,:,1) = Bext_ampl * 2.*cosy * (xyz0(1)*rr1)**3.
        Bext(:,:,:,2) = Bext_ampl * siny    * (xyz0(1)*rr1)**3.
        Bext(:,:,:,3) = 0.
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
      real :: log10pp,Teq_local1,Teq_local2
      integer :: ip
!
!  Index of the logp_ref that is just smaller than log10(pressure)
!
      log10pp = log10(press)
      ip = 1+floor((log10pp-log10p_T_ref_min)/dlog10p_T_ref)
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
                  (log10pp-logp_ref(ip))/   &
                  (logp_ref(ip+1)-logp_ref(ip))   ! unit of K
!
        tau_local = log10(tau_rad(ip))+  &
                    (log10(tau_rad(ip+1))-log10(tau_rad(ip)))*   &
                    (log10pp-logp_ref(ip))/   &
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
    subroutine prepare_eta_P
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
      if (.not.leta_file_exists) call fatal_error('prepare_eta_P', &
          'Must provide an eta(P) file')
!
      open(1,file='eta_P.txt')
      read(1,*)
      read(1,*) dlog10p_eta_ref, log10p_eta_ref_min, log10p_eta_ref_max
      read(1,*) n_eta_ref
      if(allocated(log10p_eta_ref)) deallocate(log10p_eta_ref)
      if(allocated(eta_ref)) deallocate(eta_ref)
      allocate(log10p_eta_ref(n_eta_ref),eta_ref(n_eta_ref))
      read(1,*) log10p_eta_ref,eta_ref
      if (lroot) then
        print*, 'n_eta_ref=',n_eta_ref
        print *,'Here is the eta(P) profile:'
        print *,'p [Pa] and eta[m^2/s]:'
        do i=1,n_eta_ref
          print*, 'log10p_eta_ref,eta_ref=',log10p_eta_ref(i),eta_ref(i)
        enddo
      endif
      close(1)
!
!  normalize by eta at T=T_ref
!
      eta_ref=eta_ref/eta0_ref
!
!  for debug purpose, output the result
!
      if (lroot) then
        open(1,file=trim(datadir)//'/eta_P_normalized.dat',status='replace')
        write(1,*) log10p_eta_ref
        write(1,*) eta_ref
        close(1)
      endif
!
    endsubroutine  prepare_eta_P
!***********************************************************************
    subroutine calc_eta_p(eta_local,press)
!
!  Given the local pressure p, calculate the magnetic diffusivity
!  eta_local by interpolation. The output eta_local is dimensionless.
!
!  15-mar-24/hongzhe: coded
!
      real, dimension(nx), intent(out) :: eta_local
      real, dimension(nx), intent(in) :: press
!
      real :: log10pp
      integer :: ix,ip
!
!  Index of the logp_ref that is just smaller than log10(pressure)
!
      do ix=1,nx
        log10pp = log10(press(ix))
        ip = 1+floor((log10pp-log10p_eta_ref_min)/dlog10p_eta_ref)
!
!  Interpolation for eta_local
!
        if (ip>=n_eta_ref) then
          eta_local(ix) = eta_ref(n_eta_ref)
        elseif (ip<=1) then
          eta_local(ix) = eta_ref(1)
        else
!
!  The two closest values of eta(P) given the pressure
!
          eta_local(ix) = log10(eta_ref(ip))+  &
                      (log10(eta_ref(ip+1))-log10(eta_ref(ip))) * &
                      (log10pp-log10p_eta_ref(ip))/   &
                      (log10p_eta_ref(ip+1)-log10p_eta_ref(ip))
          eta_local(ix) = 10.**eta_local(ix)
        endif
      enddo
!
    endsubroutine calc_eta_p
!***********************************************************************
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
