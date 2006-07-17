! $Id: interstellar.f90,v 1.124 2006-07-17 11:26:21 mee Exp $
!
!  This modules contains the routines for SNe-driven ISM simulations.
!  Still in development. 
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linterstellar = .true.
!
! MAUX CONTRIBUTION 1
!
!***************************************************************

module Interstellar
!
  use Cparam !, only: nx,labellen,impossible
!  use Cdata, only: ilnrho,iss,ip,x,y,z
  use Messages
!
  implicit none
!
  include 'interstellar.h'
  include 'record_types.h'
!
  type SNRemnant
    real :: x, y, z, t    ! Time and location
    double precision :: EE, MM ! Mass and energy injected
    double precision :: rhom   ! Local mean density at explosion time
    real :: radius        ! Injection radius
    real :: t_sedov
    real :: t_damping
    real :: heat_energy
    real :: damping_factor
    real :: energy_loss
    integer :: l,m,n      ! Grid position
    integer :: iproc,ipy,ipz
    integer :: SN_type
    integer :: state
  endtype
!
  type ExplosionSite
    real :: rho, lnrho, yH, lnTT, TT, ss, ee
  endtype
!
! Enumeration of Supernovae types.
!
  integer, parameter :: SNtype_I = 1, SNtype_II = 2
!
! Enumeration of Supernovae states.
!
  integer, parameter :: SNstate_invalid = 0
  integer, parameter :: SNstate_waiting = 1
  integer, parameter :: SNstate_exploding = 2
  integer, parameter :: SNstate_damping = 3
!
!
  real :: xsi_sedov=1.17   ! Estimated value for the similarity variable at shock
!
! 'Current' SN Explosion site parameters
!
  type (SNRemnant) :: SNR
!
  type (ExplosionSite) :: SNsite
!
  integer :: icooling=0
!
! Squared distance to the SNe site along the current pencil
!
  double precision, dimension(nx) :: dr2_SN     ! Pencil storing radius to SN
  double precision, dimension(nx,3) :: outward_normal_SN     ! Pencil storing radius to SN
!
!  Save space for last SNI time
!
  real :: t_next_SNI=0.0
  real :: t_interval_SNI=impossible
!
! normalisation factors for 1-d, 2-d, and 3-d profiles like exp(-r^6)
! ( 1d: 2    int_0^infty exp(-(r/a)^6)     dr) / a
!   2d: 2 pi int_0^infty exp(-(r/a)^6) r   dr) / a^2
!   3d: 4 pi int_0^infty exp(-(r/a)^6) r^2 dr) / a^3 )
! ( cf. 3.128289613 -- from where ?!? )
! NB: 1d and 2d results just from numerical integration -- calculate
!      exact integrals at some point...
!
  double precision, parameter, dimension(3) :: &
             cnorm_gaussian_SN = (/ 0.8862269255, 3.141592654, 5.568327998 /) 
  double precision, parameter, dimension(3) :: &
             cnorm_gaussian2_SN = (/ 0.9064024771, 2.784163999, 3.849760109 /)
  double precision, parameter, dimension(3) :: &
             cnorm_SN = (/ 1.855438667 , 2.805377875 , 3.712218666 /) 
  double precision, parameter, dimension(3) :: &
             cnorm_para_SN = (/  1.33333333,  1.5707963, 1.6755161 /) 
  double precision, parameter, dimension(3) :: &
             cnorm_quar_SN = (/  0.,  2.0943951, 0. /) 
!            3-D  was 3.71213666 but replaced with Maple result.... 
!
!  cp1=1/cp used to convert TT (and ss) into interstellar code units
!  (useful, as many conditions conveniently expressed in terms of TT)
!  code units based on:
!    [length]  = 1kpc  = 3.09 10^21 cm
!    [time]    = 1Gyr  = 3.15 10^16 s             !no, on [u]=1km/s...
!    [rho]     =       = 1.00 10^-24 g/cm^3
!  Lambdaunits converts coolH into interstellar code units.
!
  double precision :: unit_Lambda
!
! Remnants of the old Galaxy code formulation
!  real, parameter :: cp1=27.8   !=R * gamma / (mu * (gamma-1))  27.8 
!  real, parameter :: TTunits=46.6
!
! Minimum resulting central temperature of a SN explosion. 
! If this is not reached then consider moving mass to acheive this.
!
  real, parameter :: TT_SN_min_cgs=1.e7
  real :: TT_SN_min=impossible
  real :: TT_cutoff_cgs=100.
  real :: TT_cutoff=impossible
!
! SNe placement limitations (for code stability)
!
  double precision, parameter :: rho_SN_min_cgs=1E-28 
  real, parameter :: TT_SN_max_cgs=5E8
  real :: rho_SN_min=impossible, TT_SN_max=impossible
!
! SNI per (x,y)-area explosion rate 
!
  double precision, parameter :: SNI_area_rate_cgs=1.330982784D-56
  real :: SNI_area_rate=impossible
!
! Some useful constants 
!
  double precision, parameter :: kpc_cgs=3.086d+21         ! [cm]
  real, parameter :: yr_cgs=3.155692E7 
  double precision, parameter :: solar_mass_cgs=1.989e33
  real :: solar_mass=impossible
!
! Scale heights for SNI/II with Gaussian z distributions
!
  real, parameter :: h_SNI_cgs=1.00295e21,h_SNII_cgs=2.7774e20
  real :: h_SNI=impossible,h_SNII=impossible
!
! Self regulating SNII explosion coefficients
!
  real, parameter :: cloud_rho_cgs=1.e-24,cloud_TT_cgs=4000.
  real, parameter :: cloud_tau_cgs=2.E7 * yr_cgs
  double precision, parameter :: mass_SN_progenitor_cgs=10.*solar_mass_cgs
  real, parameter :: frac_converted=0.02,frac_heavy=0.10
!  real, parameter :: tosolarMkpc3=1.483e7
  real :: cloud_rho=impossible,cloud_TT=impossible
  real :: cloud_tau=impossible   !was =2e-2 
  real :: mass_SN_progenitor=impossible 
!
! Total SNe energy
!
  double precision, parameter :: ampl_SN_cgs=1D51
  real :: frac_ecr=0.1, frac_eth=0.9
  real :: ampl_SN=impossible
!
! SNe composition
!
  logical :: lSN_eth=.true., lSN_ecr=.true.,lSN_mass=.true., lSN_velocity=.false.
!
! Total mass added by a SNe
!
  double precision, parameter :: mass_SN_cgs=10.*solar_mass_cgs
  real :: mass_SN=impossible
  real :: velocity_SN=impossible
!
! Size of SN insertion site (energy and mass) and shell in mass movement
!
  real :: sigma_SN, sigma_SN1
  real, parameter :: width_SN_cgs=3.086E19
  real :: energy_width_ratio=1. !308533
  real :: mass_width_ratio=2.
  real :: velocity_width_ratio=1.
  real :: outer_shell_proportion = 2.
  real :: inner_shell_proportion = 1.5
  real :: width_SN=impossible
!
! Parameters for 'averaged'-SN heating
!
  real :: r_SNI_yrkpc2=4.e-6, r_SNII_yrkpc2=3.e-5
  real :: r_SNI=3.e+4, r_SNII=4.e+3
  real :: average_SNI_heating=0., average_SNII_heating=0.
!
! Limit placed of minimum density resulting from cavity creation
! 
  real, parameter :: rho_min=1.e-6
!
! Cooling timestep limiter coefficient
!
  real :: cdt_tauc=0.05
!
! COMMENT ME
!
  real :: last_SN_t=0.
!
! Time to wait before allowing SN to begin firing
!
  real :: t_settle=0.
!
! Number of randomly placed SNe to generate to mix the initial medium
!
  integer :: initial_SNI=0
  integer :: initial_SNI_done=0
!
! Parameters for UV heating of Wolfire et al.
! 
  real, parameter :: rhoUV_cgs=0.1
  real, parameter :: GammaUV_cgs=0.0147
  real, parameter :: TUV_cgs=7000.,T0UV_cgs=12000.,cUV_cgs=5.e-4
!
  real :: GammaUV=impossible,T0UV=impossible,cUV=impossible
  double precision, dimension(7) :: coolT_cgs, coolH_cgs
  real, dimension(7) :: coolB, lncoolH, lncoolT
!
  real :: coolingfunction_scalefactor=1.
  real :: heatingfunction_scalefactor=1.
!
  real :: heating_rate = 0.015
  real :: heating_rate_code = impossible
!
!
!
  logical :: ltestSN = .false.  ! If set .true. SN are only exploded at the
                                ! origin and ONLY the type I scheme is used
                                ! Used in kompaneets test etc.
!
! SN, heating, cooling type flags  
!
  logical :: lSNI=.true., lSNII=.false.
!
! Damp central regions of SNRs at early times
!
  logical :: lSNR_damping = .false.
  real :: SNR_damping = 0.
  real :: SNR_damping_time_cgs = 5E4*yr_cgs
  real :: SNR_damping_rate_cgs = 1E3*yr_cgs
  real :: SNR_damping_time = impossible
  real :: SNR_damping_rate = impossible
!
  logical :: lsmooth_coolingfunc = .false.
  logical :: laverage_SN_heating = .false.
  logical :: lheating_UV         = .true.

  logical :: lforce_locate_SNI=.false.
  logical :: uniform_zdist_SNI = .false.
!
! Requested SNe location (used for test SN)
!
  real :: center_SN_x = impossible
  real :: center_SN_y = impossible
  real :: center_SN_z = impossible 
!
! Cooling time diagnostic
!
  integer :: idiag_taucmin=0
  integer :: idiag_Hmax=0
  integer :: idiag_Lamm=0
  integer :: idiag_nrhom=0
  integer :: idiag_rhoLm=0
  integer :: idiag_Gamm=0
!
! Heating function, cooling function and mass movement
! method selection.
!
  character (len=labellen) :: cooling_select  = 'RB'
  character (len=labellen) :: heating_select  = 'wolfire'
  character (len=labellen) :: thermal_profile = 'gaussian3'
  character (len=labellen) :: velocity_profile= 'lineartanh'
!  character (len=labellen) :: velocity_profile= 'quadratictanh'
  character (len=labellen) :: mass_profile    = 'gaussian3'
  character (len=labellen) :: mass_movement   = 'off'
  character (len=labellen) :: cavity_profile  = 'gaussian3'
!
! start parameters
!
  namelist /interstellar_init_pars/ &
      cooling_select, heating_select, heating_rate
!
! run parameters
!
  namelist /interstellar_run_pars/ &
      ampl_SN, mass_SN, mass_SN_progenitor,cloud_tau, mass_width_ratio, &
      lSN_velocity, velocity_SN, velocity_width_ratio, &
      uniform_zdist_SNI, ltestSN, mass_movement, &
      lSNI, lSNII, laverage_SN_heating, coolingfunction_scalefactor, &
      lsmooth_coolingfunc, heatingfunction_scalefactor, &
      width_SN, inner_shell_proportion, outer_shell_proportion, &
      center_SN_x, center_SN_y, center_SN_z, &
      frac_ecr, frac_eth, lSN_eth, lSN_ecr, lSN_mass, &
      h_SNI, h_SNII, SNI_area_rate, rho_SN_min, TT_SN_max, &
      lheating_UV, cdt_tauc, initial_SNI, initial_SNI_done, &
      cooling_select, heating_select, heating_rate, t_settle, &
      lforce_locate_SNI, thermal_profile,velocity_profile, mass_profile, &
      t_next_SNI, SNR_damping

  contains

!***********************************************************************
    subroutine register_interstellar()
!
!  19-nov-02/tony: coded
!
      use Cdata, only: linterstellar,ip,nvar,lroot,naux,naux_com,maux,aux_count,varname,aux_var !,mvar,nvar
      use Mpicomm, only: stop_it
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_interstellar called twice')
      first = .false.
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_interstellar: ENTER'
      endif
!
      icooling = mvar+naux+1             ! indices to access uu
      naux = naux+1             ! added 3 variables
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_interstellar: naux = ', naux
        print*, 'register_interstellar: icooling = ', icooling
      endif
!
!  Put variable names in array
!
      varname(icooling) = 'cooling'
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: interstellar.f90,v 1.124 2006-07-17 11:26:21 mee Exp $")
!
! Check we aren't registering too many auxiliary variables
!
      if (naux > maux) then
        if (lroot) write(0,*) 'naux = ', naux, ', maux= ', maux
        call stop_it('register_interstellar: naux > maux')
      endif
!
!  Writing files for use with IDL
!
      if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=',cooling $'
      if (naux+naux_com  == maux+maux_com) aux_var(aux_count)=',cooling'
      aux_count=aux_count+1
      if (lroot) write(15,*) 'cooling = fltarr(mx,my,mz)*one'
!
    endsubroutine register_interstellar
!***********************************************************************
    subroutine initialize_interstellar(lstarting)
!
!  Perform any post-parameter-read initialization eg. set derived 
!  parameters
!
!  24-nov-02/tony: coded
!
!  read parameters from seed.dat and interstellar.dat
!
      use Cdata, only: datadir, directory, seed, nseed, lroot, pi, Lxyz, &
                       unit_energy, unit_density, unit_length, unit_mass, &
                       unit_temperature, unit_time, unit_velocity, unit_system, &
                       dxmax
      use General, only: random_seed_wrapper
      use Sub, only: inpui,inpup
      use Mpicomm, only: stop_it
      use EquationOfState, only: getmu
!
      logical, save :: first=.true.
      logical :: lstarting
      logical :: exist
      real :: mu
!
      if (lroot) print*,'initialize_interstellar: t_next_SNI',t_next_SNI
!
      if (lroot.and.uniform_zdist_SNI) then
         print*,'initialize_interstellar: using UNIFORM z-distribution of SNI'
      endif
!
!
      call getmu(mu) 
      if (unit_system=='cgs') then
!
!  this Lambda as such enters as n^2*Lambda(T) on the rhs of the
!  energy equation per unit volume
!
      unit_Lambda = unit_velocity**2 / unit_density / unit_time 
      elseif (unit_system=='SI') then
        call stop_it('initialize_interstellar: SI unit conversions not implemented')
      endif
      if (lroot) print*,'initialize_interstellar: unit_Lambda',unit_Lambda

!
! Mara: Initialize cooling parameters according to selection
! Default selection 'RB' Rosen & Bregman (1993)
! Alternative selection 'SS' Sanchez-Salcedo et al. (2002)
! Turn off cooling: 'no'
! cooling_select in interstellar_init_pars added
!
      if (cooling_select == 'RB') then
         if (lroot) print*,'initialize_interstellar: default RB cooling fct'
         coolT_cgs=(/ 100.D0, 2000.D0, 8000.D0, 1.D5, 1.D6, 1.D9, 0D0 /)
         coolH_cgs=(/ 2.2380D-32, 1.0012D-30, 4.6240D-36, 1.7800D-18, 3.2217D-27, tiny(0.D0), tiny(0.D0) /) / ( m_p_cgs )**2 
         coolB=(/ 2.,       1.5,      2.867,    -.65,    0.5,      tiny(0.), tiny(0.) /)
      else if (cooling_select == 'RBr') then
         if (lroot) print*,'initialize_interstellar: RB cooling fct (revised)'
         coolT_cgs=(/ 100.D0, 2000.D0, 8000.D0, 1.D5, 1.D6, 1.D9, 0D0 /)
         coolH_cgs=(/ 2.2380D-32, 1.0012D-30, 4.6240D-36, 1.7800D-18, 3.2217D-27, tiny(0.D0), tiny(0.D0) /) / ( m_p_cgs )**2 
         coolB=(/ 2.,       1.5,      2.867,    -.65,    0.5,      tiny(0.), tiny(0.) /)
      else if (cooling_select == 'SS') then
         ! These are the SS et al (2002) coefficients multiplied by m_proton**2
         ! to obtain same units as RB above
         if (lroot) print*,'initialize_interstellar: SS cooling function'
         coolT_cgs=(/   10D0,   141D0,   313D0,  6102D0,      1D5,       1D9, 0D0 /)
         coolH_cgs=(/3.42D16, 9.10D18, 1.11D20,  2.00D8, tiny(0D0), tiny(0D0), tiny(0D0) /)
         coolB    =(/   2.12,     1.0,    0.56,    3.67,     -.65 , tiny(0.), tiny(0.) /)
      else if (cooling_select == 'SSr') then
         ! revised to make continuous
         if (lroot) print*,'initialize_interstellar: revised SS cooling fct'
         coolT_cgs=(/   10D0,   141D0,    313D0, 6102D0,     1D5,       1D9, 0D0 /)
         coolH_cgs=(/3.70D16, 9.46D18, 1.185D20, 2.00D8, 7.96D29, tiny(0D0), tiny(0D0) /)
         coolB    =(/   2.12,     1.0,     0.56,   3.67,   -0.65, tiny(0.) , tiny(0.) /)
      else if (cooling_select == 'SSrr') then
         ! revised to make continuous
         if (lroot) print*,'initialize_interstellar: revised SS cooling fct'
         coolT_cgs=(/   10D0,   141D0,    313D0, 6102D0,     1D5,       4D7, 1D10 /)
         coolH_cgs=(/3.703109927416290D16, &
                     9.455658188464892D18, &
                     1.185035244783337D20, &
                     1.9994576479D8, &
                     7.96D29, &
                     1.440602814622207D21, &
                     0D0 /)
         coolB    =(/   2.12,     1.0,     0.56,   3.67,   -0.65, 0.5, tiny(0.)  /)
      else if (cooling_select == 'off') then
         if (lroot) print*,'initialize_interstellar: no cooling applied'
         coolT_cgs=tiny(0.D0)
         coolH_cgs=tiny(0.D0)
         coolB=tiny(0.)
      end if
!

      lncoolH = real(dlog(coolH_cgs / unit_Lambda * (unit_temperature**coolB) * coolingfunction_scalefactor))
      lncoolT = real(dlog(coolT_cgs / unit_temperature))
!
      heating_rate_code=heating_rate*real(unit_length/unit_velocity**3)
!
      if (unit_system=='cgs') then
        if (TT_SN_max==impossible) TT_SN_max=TT_SN_max_cgs / unit_temperature
        if (rho_SN_min==impossible) rho_SN_min=rho_SN_min_cgs / unit_density
        TT_SN_min=TT_SN_min_cgs / unit_temperature
        TT_cutoff=TT_cutoff / unit_temperature
        if (SNI_area_rate==impossible) SNI_area_rate=SNI_area_rate_cgs * unit_length**2 * unit_time
        if (h_SNI==impossible)         h_SNI=h_SNI_cgs / unit_length
        h_SNII=h_SNII_cgs / unit_length
        solar_mass=solar_mass_cgs / unit_mass
        if (lroot) print*,'initialize_interstellar: solar_mass (code) =',solar_mass
        cloud_rho=cloud_rho_cgs / unit_density
        cloud_TT=cloud_TT_cgs / unit_temperature
        r_SNI =r_SNI_yrkpc2  * (unit_time/yr_cgs) * (unit_length/kpc_cgs)**2
        r_SNII=r_SNII_yrkpc2 * (unit_time/yr_cgs) * (unit_length/kpc_cgs)**2
        T0UV=T0UV_cgs / unit_temperature
        cUV=cUV_cgs / unit_temperature
                !GammaUV=GammaUV_cgs !/ (unit_velocity**2 * unit_time)
        if (GammaUV == impossible) &
                GammaUV=GammaUV_cgs * real(unit_length/unit_velocity**3)
        if (ampl_SN == impossible) ampl_SN=ampl_SN_cgs / unit_energy 
        if (lroot) print*,'initialize_interstellar: ampl_SN = ',ampl_SN
        if (cloud_tau == impossible) cloud_tau=cloud_tau_cgs / unit_time 
        if (mass_SN == impossible) mass_SN=mass_SN_cgs / unit_mass 
        if (mass_SN_progenitor == impossible) &
                mass_SN_progenitor=mass_SN_progenitor_cgs / unit_mass 
        if (width_SN == impossible) width_SN=width_SN_cgs / unit_length 
        if (SNR_damping_time == impossible) SNR_damping_time=SNR_damping_time_cgs / unit_time 
        if (SNR_damping_rate == impossible) SNR_damping_rate=SNR_damping_rate_cgs / unit_time 
      else
        call stop_it('initialize_interstellar: SI unit conversions not implemented')
      endif
!
!  SNRdamping factor
!
      if (SNR_damping/=0.) lSNR_damping=.true.
!
!  Slopeyness used for tanh rounding profiles etc.
!
      sigma_SN=dxmax*3
      sigma_SN1=1./sigma_SN
!
      t_interval_SNI = 1./(SNI_area_rate * Lxyz(1) * Lxyz(2))
      average_SNI_heating =r_SNI *ampl_SN/(sqrt(pi)*h_SNI )*heatingfunction_scalefactor
      average_SNII_heating=r_SNII*ampl_SN/(sqrt(pi)*h_SNII)*heatingfunction_scalefactor
      if (lroot) print*,'initialize_interstellar: t_interval_SNI =',t_interval_SNI,Lxyz(1),Lxyz(2),SNI_area_rate
!
      if (lroot.and.ip<14) then
        print*,'initialize_interstellar: nseed,seed',nseed,seed(1:nseed)
        print*,'initialize_interstellar: finished'
      endif
!
      if (lroot.and. (.not. lstarting)) then
         open(1,file=trim(datadir)//'/sn_series.dat',position='append')
         write(1,'("#",3A)')  &
          '---it-----t----------itype---iproc----l--m--n----', &
          '----x-----------y--------z-----', &
          '--rho--TT----EE----'
         close(1)
      endif
!
      if (ltestSN) then
        t_interval_SNI=1.E10
        t_next_SNI=0.
      endif
!
      if (lroot) print*,"initialize_interstellar: t_next_SNI=",t_next_SNI
!
    endsubroutine initialize_interstellar
!***********************************************************************
    subroutine input_persistent_interstellar(id,lun,done)
!
!  Read in the stored time of the next SNI
!
      use Cdata, only: lroot
!
      integer :: id,lun
      logical :: done
!
      if (id==id_record_T_NEXT_SN) then
        read (lun) t_next_SNI
        done=.true.
      endif
      if (lroot) print*,'input_persistent_interstellar: ', t_next_SNI
!
    endsubroutine input_persistent_interstellar
!***********************************************************************
    subroutine output_persistent_interstellar(lun)
!
!  Writes out the time of the next SNI
!
      use Cdata, only: lroot
!
      integer :: lun
!
      if (lroot.and.lSNI.and.lSNII) &
        print*,'output_persistent_interstellar: ', t_next_SNI
      write (lun) id_record_T_NEXT_SN
      write (lun) t_next_SNI
!
    endsubroutine output_persistent_interstellar
!***********************************************************************
    subroutine rprint_interstellar(lreset,lwrite)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Cdata, only: nname, cname, cform
      use Sub, only: parse_name
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        SNR_damping=0.
        idiag_taucmin=0
        idiag_Hmax=0
        idiag_Lamm=0
        idiag_nrhom=0
        idiag_rhoLm=0
        idiag_Gamm=0

!        TT_SN_max=impossible 
!        rho_SN_min=impossible 
!        SNI_area_rate=impossible 
!        h_SNI=impossible 
!        GammaUV=impossible 
!        width_SN=impossible 
!        ampl_SN=impossible 
!        mass_SN=impossible 
!        mass_SN_progenitor=impossible 
!        cloud_tau=impossible
     endif

     lpenc_requested(i_ee)=.true.
     lpenc_requested(i_lnTT)=.true.
     lpenc_requested(i_TT1)=.true.
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'taucmin',idiag_taucmin)
        call parse_name(iname,cname(iname),cform(iname),'Hmax',idiag_Hmax)
        call parse_name(iname,cname(iname),cform(iname),'Lamm',idiag_Lamm)
        call parse_name(iname,cname(iname),cform(iname),'nrhom',idiag_nrhom)
        call parse_name(iname,cname(iname),cform(iname),'rhoLm',idiag_rhoLm)
        call parse_name(iname,cname(iname),cform(iname),'Gamm',idiag_Gamm)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_taucmin=',idiag_taucmin
        write(3,*) 'i_Hmax=',idiag_Hmax
        write(3,*) 'i_Lamm=',idiag_Lamm
        write(3,*) 'i_nrhom=',idiag_nrhom
        write(3,*) 'i_rhoLm=',idiag_rhoLm
        write(3,*) 'i_Gamm=',idiag_Gamm
        write(3,*) 'icooling=',icooling
      endif
!
    endsubroutine rprint_interstellar
!***********************************************************************
    subroutine read_interstellar_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
                                                                                                   
      if (present(iostat)) then
        read(unit,NML=interstellar_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=interstellar_init_pars,ERR=99)
      endif
                                                                                                   
                                                                                                   
99    return
    endsubroutine read_interstellar_init_pars
!***********************************************************************
    subroutine write_interstellar_init_pars(unit)
      integer, intent(in) :: unit
                                                                                                   
      write(unit,NML=interstellar_init_pars)
                                                                                                   
    endsubroutine write_interstellar_init_pars
!***********************************************************************
    subroutine read_interstellar_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
                                                                                                   
      if (present(iostat)) then
        read(unit,NML=interstellar_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=interstellar_run_pars,ERR=99)
      endif
                                                                                                   
                                                                                                   
99    return
    endsubroutine read_interstellar_run_pars
!***********************************************************************
    subroutine write_interstellar_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=interstellar_run_pars)
!
    endsubroutine write_interstellar_run_pars
!***********************************************************************
    subroutine pencil_criteria_interstellar()
! 
!  All pencils that the Interstellar module depends on are specified here.
! 
!  26-mar-05/tony: coded
!  11-mar-06/axel: added idiag_nrhom
!
      use Cdata
!
      lpenc_requested(i_ee)=.true.
      lpenc_requested(i_lnrho)=.true.
      lpenc_requested(i_lnTT)=.true.
      lpenc_requested(i_TT1)=.true.
!
!  diagnostic pencils
!
!AB:  if (idiag_nrhom/=0) lpenc_diagnos(i_rho)=.true.
!
    endsubroutine pencil_criteria_interstellar
!***********************************************************************
    subroutine calc_heat_cool_interstellar(f,df,p,Hmax)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  We may want to move it to the entropy module for good, because its use
!  is not restricted to interstellar runs (could be used for solar corona).
!  Also, it doesn't pose an extra load on memory usage or compile time.
!  (We should allow that UV heating can be turned off; so rhoUV should
!  be made an input parameter.)
!
!  19-nov-02/graeme: adapted from calc_heat_cool
!  10-aug-03/axel: TT is used as input
!   3-apr-06/axel: add ltemperature switch
!
      use Cdata, only: lroot, headtt, lfirst, ldiagnos, ldt, m, n,  &
                       iss, ilnTT, unit_length, unit_velocity, dt1_max, z, &
                       datadir, pretend_lnTT, t
!
      use Sub, only: max_mn_name, sum_mn_name
      use EquationOfState, only: gamma, gamma11
!
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case) :: p
      
      real, dimension (nx), intent(inout) :: Hmax
      real, dimension (nx) :: heat,cool, cool_beta
      real, dimension (nx) :: damp_profile
      real :: minqty
!      real ::  mu
      integer :: i
!
!  identifier
!
      if(headtt) print*,'calc_heat_cool_interstellar: ENTER'
!
!  define T in K, for calculation of both UV heating and radiative cooling
!
!  add T-dept radiative cooling, from Rosen et al., ApJ, 413, 137, 1993 ('RB')
!  OR
!  Sanchez-Salcedo et al. ApJ, 577, 768, 2002 ('SS').
!  cooling is Lambda*rho^2, with (eq 7)
!     Lambda=coolH(i)*TT*coolB(i),   for coolT(i) <= T < coolT(i+1)
!  nb: our coefficients coolH(i) differ from those in Rosen et al. by
!   factor (mu mp)^2, with mu=1.2, since Rosen works in number density, n.
!   (their cooling = Lambda*n^2,  rho=mu mp n.)
!  The factor Lambdaunits converts from cgs units to code units.
!
!  [Currently, coolT(1) is not modified, but this may be necessary
!  to avoid creating gas too cold to resolve.]
!
      cool=0.0
      cool_beta=1.0
cool_loop: do i=1,6
        if (lncoolT(i) .ge. lncoolT(i+1)) exit cool_loop
        where (lncoolT(i) <= p%lnTT .and. p%lnTT < lncoolT(i+1))
               cool=cool+exp(lncoolH(i)+p%lnrho+p%lnTT*coolB(i))
               cool_beta=coolB(i)
        endwhere
      enddo cool_loop

!
!  possibility of temporal smoothing of cooling function
!
      if (lsmooth_coolingfunc) cool=(cool+f(l1:l2,m,n,icooling))*0.5
!
!      open(1,file=trim(datadir)//'/cooling.dat',position='append')
!      do i=1,nx
!        write(1,'(4e15.8)') exp(p%lnTT(i)), cool(i)/exp(p%lnrho(i)), cool(i), exp(p%lnrho)
!      enddo
!      close(1)
!print*,'min,max interstellar heating: ',minval(cool),maxval(cool)
!
!  add UV heating, cf. Wolfire et al., ApJ, 443, 152, 1995
!  with the values above, this gives about 0.012 erg/g/s (T < ~1.e4 K)
!  nb: need rho0 from density_[init/run]_pars, if i want to implement
!      the the arm/interarm scaling.
!
!      heat=0.0
!tony: DISABLE UV HEATING -- Requires reformulation
!      if (lheating_UV) then
!        heat(1:nx)=GammaUV*0.5*(1.0+tanh(cUV*(T0UV-exp(lnTT(1:nx)))))
!      endif
!
!tony: need to do unit_system stuff with scale heights etc.
!
! Mara:UV heating off, heating_select='off' in interstellar_init_pars
!      Constant heating with a rate heating_rate[erg/g/s]: 'cst'
!      Heating balancing the initial cooling function: 'eql'
!      Default: 'off' with heating_rate = 0.015
      If (heating_select == 'cst') Then
         heat = heating_rate_code
      Else If (heating_select == 'wolfire') Then
        heat(1:nx)=GammaUV*0.5*(1.0+tanh(cUV*(T0UV-exp(p%lnTT))))
      Else If (heating_select == 'wolfire_min') Then
        heat(1:nx)=GammaUV*0.5*(1.0+tanh(cUV*(T0UV-exp(p%lnTT))))
        heat = max(heat,heating_rate_code)
      Else If (heating_select == 'eql') Then
         If (headtt) Then 
            heating_rate = cool(1)
            If (lroot) Then
               Print*,'Heating balancing the initial cooling profile (eql)'
               Print*,'Note: works only for unstratified cases when cool is z-independent!'
               Print*,'      heating_rate is overwritten'
            End If
         End If
         heat = heating_rate
      Else If (heating_select == 'off') Then
         heat = 0.
      End If

      if (SNR%state==SNstate_damping) then
        call proximity_SN()        
        minqty=0.5*(1.+tanh((t-SNR%t_damping)/SNR_damping_rate))
        damp_profile=((1.-minqty)*0.5*(1.+tanh((sqrt(dr2_SN)-(SNR%radius*2.))*sigma_SN1-2.))+minqty)
        cool=cool*damp_profile
        heat=heat*damp_profile
      endif
!      where (p%lnTT < TT_cutoff) cool=min(cool,heat)
!      Print*,'Debug: cool, heat',cool,heat
!
!  Average SN heating (due to SNI and SNII)
!  The amplitudes of both types is assumed the same (=ampl_SN)
!
 
      if (laverage_SN_heating) then
        heat=heat+average_SNI_heating *exp(-(z(n)/h_SNI )**2)
        heat=heat+average_SNII_heating*exp(-(z(n)/h_SNII)**2)
      endif
!
!  prepare diagnostic output
!
      if(ldiagnos) then
        if(idiag_Hmax/=0) &
          call max_mn_name(heat,idiag_Hmax)
        if(idiag_taucmin/=0) &
          call max_mn_name(cool/p%ee,idiag_taucmin,lreciprocal=.true.)  
        if(idiag_Lamm/=0) &
          call sum_mn_name(cool,idiag_Lamm)  
        if(idiag_nrhom/=0) &
          call sum_mn_name(cool/p%ee,idiag_nrhom)  
!--       call sum_mn_name(cool*p%rho/p%ee,idiag_nrhom)  
!AB: the factor rho is already included in cool, so cool=rho*Lambda
        if(idiag_rhoLm/=0) &
          call sum_mn_name(p%rho*cool,idiag_rhoLm)  
        if(idiag_Gamm/=0) &
          call sum_mn_name(p%rho*heat,idiag_Gamm)  
      endif

! Limit timestep by the cooling time (having subtracted any heating) 
!    dt1_max=max(dt1_max,cdt_tauc*(cool)/ee,cdt_tauc*(heat)/ee)

    if (lfirst.and.ldt) then
        !dt1_max=max(dt1_max,(cool-heat)/(p%ee*cdt_tauc))
        dt1_max=max(dt1_max,cool/(p%ee*cdt_tauc))
        !dt1_max=max(dt1_max,(cool-heat)/(p%ee*cdt_tauc)*max(abs(cool_beta-1.),1.))
      !where ((heat-cool) > 0.) &
      !  dt1_max=max(dt1_max,heat/(p%ee*cdt_tauc))
!      dt1_max=max(dt1_max,(cool-0.99*heat)/(p%ee*cdt_tauc)*max(abs(cool_beta-1.)**4,1.))
!print*,'dt1_max(1),p%ee(1)=',dt1_max(1),p%ee(1)
            !cdt_tauc*(heat)/ee)
      Hmax=Hmax+heat
    endif
!
!  For clarity we have constructed the rhs in erg/s/g [=T*Ds/Dt]
!  so therefore we now need to multiply by TT1.
!
       f(l1:l2,m,n,icooling)=cool
       if (ltemperature) then
         df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+p%TT1*(heat-cool)*gamma
       elseif (pretend_lnTT) then
         df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+p%TT1*(heat-cool)*gamma
       else
         df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+p%TT1*(heat-cool)
       endif
!      
    endsubroutine calc_heat_cool_interstellar
!***********************************************************************
    subroutine check_SN(f,df)
!
!  Checks for SNe, and implements appropriately:
!   relevant subroutines in entropy.f90
!
    use Cdata, only: headtt, t
!
    real, dimension(mx,my,mz,mvar+maux) :: f
    real, dimension(mx,my,mz,mvar) :: df
    logical :: l_SNI=.false.   !only allow SNII if no SNI this step
                               !(may not be worth keeping)
!
    intent(inout) :: f
!
!  identifier  
!
      if(headtt) print*,'check_SN: ENTER'
!
!  Do separately for SNI (simple scheme) and SNII (Boris' scheme)
!
    if (t < t_settle) return
    call calc_snr_damping_factor(f)
    if (lSNI) call check_SNI (f,df,l_SNI)
    if (lSNII) call check_SNII(f,df,l_SNI)
    call calc_snr_damping_add_heat(f)
!
    endsubroutine check_SN
!***********************************************************************
    subroutine check_SNI(f,df,l_SNI)
!
!  If time for next SNI, then implement, and calculate time of subsequent SNI
!
    use Cdata, only: headtt, lroot, t
    use General, only: random_number_wrapper
!
    real, dimension(mx,my,mz,mvar+maux) :: f
    real, dimension(mx,my,mz,mvar) :: df
    real, dimension(1) :: franSN
    logical :: l_SNI
    integer :: try_count
!
    intent(inout) :: f,l_SNI
!
!  identifier
!
    if(headtt) print*,'check_SNI: ENTER'
!
    l_SNI=.false.
    if (initial_SNI_done < initial_SNI) then
      SNR%t=t
      SNR%state=SNstate_exploding
      SNR%radius=width_SN
      if (uniform_zdist_SNI) then
        call position_SN_uniformz(f)
      else
        call position_SN_gaussianz(f,h_SNII)
      endif
      call explode_SN(f,df,1)
      initial_SNI_done=initial_SNI_done+1
    elseif (t >= t_next_SNI) then
      SNsite%TT=1E20
      SNsite%rho=0.
      SNR%t=t
      SNR%radius=width_SN
      SNR%state=SNstate_waiting
      try_count=100
      do while ((SNsite%rho .lt. rho_SN_min) .or. (SNsite%TT .gt. TT_SN_max))
        if (ltestSN) then  
          call position_SN_testposition(f)
        elseif (uniform_zdist_SNI) then
          call position_SN_uniformz(f)
        else
          call position_SN_gaussianz(f,h_SNI)
        endif

        if (lforce_locate_SNI.and. &
             ((SNsite%rho .lt. rho_SN_min) .or. &
              (SNsite%TT .gt. TT_SN_max))) then
          call find_nearest_SNI(f)
        endif

        try_count=try_count-1
        if (try_count.eq.0) then
          if (lroot) print*,"check_SNI: 100 RETRIES OCCURED - skipping SNI insertion"
          exit
        endif
      enddo
       if (try_count.ne.0) then
          SNR%state=SNstate_exploding
          call explode_SN(f,df,1)
       !  pre-determine time for next SNI
          if (lroot.and.ip<14) print*,"check_SNI: Old t_next_SNI=",t_next_SNI
          call random_number_wrapper(franSN)   
          t_next_SNI=t + (1.0 + 0.4*(franSN(1)-0.5)) * t_interval_SNI
          if (lroot.and.ip<20) print*,'check_SNI: Next SNI at time = ',t_next_SNI
       l_SNI=.true.
      endif
    endif
!
    endsubroutine check_SNI
!***********************************************************************
    subroutine check_SNII(f,df,l_SNI)
!
!  Check for SNII, via self-regulating scheme.
!
    use Cdata, only: dx, dy, dz, t, ilnrho, headtt, lroot
    use General, only: random_number_wrapper
    use Mpicomm, only: mpireduce_sum, mpibcast_real
    use EquationOfState, only: eoscalc
! 
    real, dimension(mx,my,mz,mvar+maux) :: f
    real, dimension(mx,my,mz,mvar) :: df
    real, dimension(nx) :: rho,rho_cloud,ss,lnTT,TT,yH
    real :: cloud_mass,cloud_mass_dim,freq_SNII,prob_SNII,dv
    real, dimension(1) :: franSN,fsum1,fsum1_tmp,fmpi1
    real, dimension(ncpus) :: cloud_mass_byproc
    integer :: icpu, m, n
    logical :: l_SNI
    real :: dtsn
!
    intent(in) :: l_SNI
    intent(inout) :: f
!
!  identifier
!
    if(lroot.and.headtt.and.ip<14) print*,'check_SNII: ENTER'
!
    if (l_SNI) return         ! only do if no SNI this step
!
    cloud_mass=0.0
    do n=n1,n2
      do m=m1,m2
        rho(1:nx)=exp(f(l1:l2,m,n,ilnrho))
        call eoscalc(f,nx,yH=yH,lnTT=lnTT)
        TT=exp(lnTT)

        rho_cloud(1:nx)=0.0
        !print*,'min,max TT:  ', minval(TT), maxval(TT)
        !print*,'min,max rho: ', minval(rho), maxval(rho)

        where (rho(1:nx) >= cloud_rho .and. TT(1:nx) <= cloud_TT)   &
          rho_cloud(1:nx)=rho(1:nx)
        cloud_mass=cloud_mass+sum(rho_cloud(1:nx))
      enddo
    enddo
    fsum1_tmp=(/ cloud_mass /)
!       print*,'check_SNII, iproc,fsum1_tmp:',iproc,fsum1_tmp(1)
    call mpireduce_sum(fsum1_tmp,fsum1,1) 
    call mpibcast_real(fsum1,1)
    dv=1.
    if (nxgrid/=1) dv=dv*dx
    if (nygrid/=1) dv=dv*dy
    if (nzgrid/=1) dv=dv*dz
    cloud_mass_dim=fsum1(1)*dv/solar_mass
    if (lroot .and. ip < 14) &
          print*,'check_SNII: cloud_mass_dim,fsum(1),dv,solar_mass:',cloud_mass_dim,fsum1(1),dv,solar_mass
    !if (franSN(1) <= prob_SNII) then
    !  print*,'check_SNII: iproc,fsum1:',iproc,fsum1(1)
    ! need convert to dimensional units, for rate/probability calculation only. 
    ! don't overwrite cloud_mass (on individual processors), as it's re-used.
    !if (lroot .and. ip < 14) &
    !     print*,'check_SNII: cloud_mass_dim:',cloud_mass_dim
    !
    dtsn=t-last_SN_t
    freq_SNII= &
      frac_heavy*frac_converted*cloud_mass_dim/mass_SN_progenitor/cloud_tau
    prob_SNII=freq_SNII*dtsn
    call random_number_wrapper(franSN)   
    !if (lroot .and. ip < 20) then
    if (cloud_mass .gt. 0.) then
      print*,'check_SNII: freq,prob,rnd,dtsn:',freq_SNII,prob_SNII,franSN(1),dtsn
      print*,'check_SNII: frac_heavy,frac_converted,cloud_mass_dim,mass_SN,cloud_tau',&
              frac_heavy,frac_converted,cloud_mass_dim,mass_SN,cloud_tau
    endif
    !endif
    if (franSN(1) <= prob_SNII) then
      SNR%state=SNstate_waiting
      !  position_SN_bycloudmass needs the cloud_masses for each processor;  
      !   communicate and store them here, to avoid recalculation.
      cloud_mass_byproc(:)=0.0
      ! use non-root broadcasts for the communication...
      do icpu=1,ncpus
        fmpi1=cloud_mass
        call mpibcast_real(fmpi1,1,icpu-1)
        cloud_mass_byproc(icpu)=fmpi1(1)
      enddo
      ! if (lroot.and.ip<14) print*,'check_SNII: cloud_mass_byproc:',cloud_mass_byproc
      call position_SN_bycloudmass(f,cloud_mass_byproc)
      SNR%t=t
      SNR%radius=width_SN
      SNR%state=SNstate_exploding
      call explode_SN(f,df,2)
      last_SN_t=0.
    endif
    !
    endsubroutine check_SNII
!***********************************************************************
    subroutine position_SN_testposition(f)
!
!   determine position for next SN (w/ fixed scale-height)
!
    use Cdata, only: headtt, dx, dy, dz, lroot, lperi, xyz0
!    use General
!
    real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
!
    real :: z00, x00, y00
    real, dimension(3) :: fran3
    integer :: i  
!
    if(headtt) print*,'position_SN_testposition: ENTER'
!    
!  Calculate the global (nzgrid) lower z-coordinate
!
    if (lperi(1)) then; x00=xyz0(1)-.5*dx; else; x00=xyz0(1); endif
    if (lperi(2)) then; y00=xyz0(2)-.5*dy; else; y00=xyz0(2); endif
    if (lperi(3)) then; z00=xyz0(3)-.5*dz; else; z00=xyz0(3); endif
!
!  Pick SN position (SNR%l,SNR%m,SNR%n)
!
    if (lroot) then
      if (center_SN_x.eq.impossible) then
        i=max(int(nxgrid/2)+1,1)
      else
        i=int((center_SN_x-x00)/dx)
      endif
      SNR%l=i+nghost
!
      if (center_SN_y.eq.impossible) then
        i=max(int(nygrid/2)+1,1)
      else
        i=int((center_SN_y-y00)/dy)
      endif
      SNR%ipy=(i-1)/ny  ! uses integer division
      SNR%m=i-(SNR%ipy*ny)+nghost
!
      if (center_SN_z.eq.impossible) then
        i=max(int(nzgrid/2)+1,1)
      else
        i=int((center_SN_z-z00)/dz)
      endif
      SNR%ipz=(i-1)/nz   ! uses integer division
      SNR%n=i-(SNR%ipz*nz)+nghost
      SNR%iproc=SNR%ipz*nprocy + SNR%ipy
    endif
    call share_SN_parameters(f)
!
    endsubroutine position_SN_testposition
!***********************************************************************
    subroutine position_SN_gaussianz(f,h_SN)
!
!   determine position for next SN (w/ fixed scale-height)
!
    use Cdata, only: headtt, dx, dy, dz, lroot, lperi, xyz0
    use General, only: random_number_wrapper
!    use Mpicomm
!    use General
!
    real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
    real, intent(in) :: h_SN
!
    real, dimension(nzgrid) :: cum_prob_SN
    real :: zn, z00, x00, y00
    real, dimension(3) :: fran3
    integer :: i, nzskip=10   !prevent SN from being too close to boundaries
!
!
!
    if(headtt) print*,'position_SN_gaussianz: ENTER'
!    
!  Calculate the global (nzgrid) lower z-coordinate
!
    if (lperi(1)) then; x00=xyz0(1)+.5*dx; else; x00=xyz0(1); endif
    if (lperi(2)) then; y00=xyz0(2)+.5*dy; else; y00=xyz0(2); endif
    if (lperi(3)) then; z00=xyz0(3)+.5*dz; else; z00=xyz0(3); endif
!
!  Pick SN position (SNR%l,SNR%m,SNR%n)
!
    call random_number_wrapper(fran3)    ! get 3 random numbers
!                                        ! on all processors to keep
!                                        ! rnd. generators in sync
    if (lroot) then
      i=int(fran3(1)*nxgrid)+1
      SNR%l=i+nghost
!
      i=int(fran3(2)*nygrid)+1
      SNR%ipy=(i-1)/ny  ! uses integer division
      SNR%m=i-(SNR%ipy*ny)+nghost
!
!  Cumulative probability function in z currently calculated each time.
!  It's constant, and could be stored (and calculated in init)
!
      cum_prob_SN(1:nzskip)=0.0
      do i=nzskip+1,nzgrid-nzskip
        zn=z00+(i-1)*dz
        cum_prob_SN(i)=cum_prob_SN(i-1)+exp(-(zn/h_SN)**2)
      enddo
      cum_prob_SN=cum_prob_SN/cum_prob_SN(nzgrid-nzskip)
!
!  The following should never be needed, but just in case floating point 
!  errors ever lead to cum_prob_SNI(nzgrid-nzskip) < rnd < 1.
!
      cum_prob_SN(nzgrid-nzskip+1:nzgrid)=1.0   
          
      do i=nzskip+1,nzgrid-nzskip
        if (cum_prob_SN(i-1) <= fran3(3) .and. fran3(3) < cum_prob_SN(i)) &
           then
          SNR%ipz=(i-1)/nz  ! uses integer division
          SNR%n=i-(SNR%ipz*nz)+nghost
          exit
        endif
      enddo
      SNR%iproc=SNR%ipz*nprocy + SNR%ipy
    endif
!
    call share_SN_parameters(f)
!
    endsubroutine position_SN_gaussianz
!***********************************************************************
    subroutine position_SN_uniformz(f)
!
!   determine position for next SN (w/ fixed scale-height)
!
    use Cdata, only: headtt, dx, dy, dz, lroot, lperi, xyz0
    use General, only: random_number_wrapper
!    use Mpicomm
!    use General
!
    real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
!
    real :: z00, x00, y00
    real, dimension(3) :: fran3
    integer :: i   !prevent SN from being too close to boundaries
!
    if(headtt) print*,'position_SN_uniformz: ENTER'
!    
!  Calculate the global (nzgrid) lower z-coordinate
!
    if (lperi(1)) then; x00=xyz0(1)+.5*dx; else; x00=xyz0(1); endif
    if (lperi(2)) then; y00=xyz0(2)+.5*dy; else; y00=xyz0(2); endif
    if (lperi(3)) then; z00=xyz0(3)+.5*dz; else; z00=xyz0(3); endif
!
!  Pick SN position (SNR%l,SNR%m,SNR%n)
!
    call random_number_wrapper(fran3)    ! get 3 random numbers
!                                         ! on all processors to keep
!                                         ! rnd. generators in sync
    if (lroot) then
      i=int(fran3(1)*nxgrid)+1
      if (nxgrid==1) i=1
      SNR%l=i+nghost

      i=int(fran3(2)*nygrid)+1
      if (nygrid==1) i=1
      SNR%ipy=(i-1)/ny  ! uses integer division
      SNR%m=i-(SNR%ipy*ny)+nghost

      i=int(fran3(3)*nzgrid)+1
      if (nzgrid==1) i=1
      SNR%ipz=(i-1)/nz   ! uses integer division
      SNR%n=i-(SNR%ipz*nz)+nghost
      SNR%iproc=SNR%ipz*nprocy + SNR%ipy
    endif
!
    call share_SN_parameters(f)
!
    endsubroutine position_SN_uniformz
!***********************************************************************
    subroutine position_SN_bycloudmass(f,cloud_mass_byproc)
!
!  Determine position for next SNII (using Boris' scheme)
!  It seems impractical to sort all high density points across all processors;
!  instead, we just construct cumulative pdfs that allow us to pick a processor,
!  and then a point on that processor, with probability proportional to rho.
!  As a result, the SN position is *not* independent of ncpus (or of nprocy 
!  and nprocz).  (It is repeatable given fixed nprocy/z though.)
!
    use Cdata, only: lroot, iproc, ilnrho
    use General, only: random_number_wrapper
!    use Mpicomm
    use EquationOfState, only: eoscalc
!
    real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
    real, intent(in) , dimension(ncpus) :: cloud_mass_byproc
!
    real, dimension(0:ncpus) :: cum_prob_byproc
    real, dimension(1) :: franSN
    real :: cloud_mass,cum_mass,cum_prob_onproc
    real, dimension(nx) :: lnrho,rho,ss,lnTT,TT,yH
    integer :: icpu,l,m,n
!
!
!  identifier
!
      if(lroot.and.ip<14) print*,'position_SN_bycloudmass: ENTER'
!
!  Construct cumulative distribution function, using cloud_mass_byproc.
!  NB: icpu=iproc+1 (iproc in [0,ncpus-1], icpu in [1,ncpus] )
!
    cum_prob_byproc=0.0
    do icpu=1,ncpus
      cloud_mass=cloud_mass_byproc(icpu)
      cum_prob_byproc(icpu)=cum_prob_byproc(icpu-1)+cloud_mass_byproc(icpu)
    enddo
    cum_prob_byproc(:)=cum_prob_byproc(:)/cum_prob_byproc(ncpus)
    if (lroot.and.ip<14) then
      print*,'position_SN_bycloudmass: cloud_mass_byproc=',cloud_mass_byproc
      print*,'position_SN_bycloudmass: cum_prob_byproc=',cum_prob_byproc
      print*,'position_SN_bycloudmass: cloud_mass=',cloud_mass
    endif
!     
!  Use random number to detemine which processor SN is on.
!  (Use root processor for rand, to ensure repeatability.)
!
    call random_number_wrapper(franSN)   
    do icpu=1,ncpus
      if (cum_prob_byproc(icpu-1) <= franSN(1) .and.                      &
           franSN(1) < cum_prob_byproc(icpu)) then
        SNR%iproc=icpu-1 
        exit
      endif
    enddo
    if (lroot.and.ip<14) &
          print*, 'position_SN_bycloudmass: franSN(1),SNR%iproc=',franSN(1),SNR%iproc
!
!  Use random number to pick SNII location on the right processor.
!  (No obvious reason to re-use the original random number for this.)
!    franSN(1)=(franSN(1)-cum_prob_byproc(SNR%iproc)) /                      &
!              (cum_prob_byproc(SNR%iproc+1)-cum_prob_byproc(SNR%iproc))
!
    call random_number_wrapper(franSN)   
    if (iproc == SNR%iproc) then
      cum_mass=0.0
find_SN: do n=n1,n2
      do m=m1,m2
        lnrho=f(l,m,n,ilnrho)
        rho=exp(lnrho)
        call eoscalc(f,nx,yH=yH,lnTT=lnTT)
        TT=exp(lnTT)
        do l=1,nx
          if (rho(l) >= cloud_rho .and. TT(l) <= cloud_TT) then
            cum_mass=cum_mass+rho(l)
            cum_prob_onproc=cum_mass/cloud_mass
            if (franSN(1) <= cum_prob_onproc) then
              SNR%l=l+l1-1; SNR%m=m; SNR%n=n
              if (lroot.and.ip<14) &
                print*,'position_SN_bycloudmass: cum_mass,cum_prob_onproc,franSN(1)=', &
                                 cum_mass,cum_prob_onproc,franSN(1)
              exit find_SN
            endif
          endif
        enddo
      enddo
      enddo find_SN
    endif
!
    call share_SN_parameters(f)
!
    endsubroutine position_SN_bycloudmass
!***********************************************************************
    subroutine find_nearest_SNI(f)
!
!   Given a presently unsuitable SNI explosion site... Find the nearest 
!   suitable location
!
    use Cdata, only: headtt,iproc,ilnrho
    use EquationOfState, only: eoscalc
    use Mpicomm, only: mpibcast_int
!    use Mpicomm
    use General, only: random_number_wrapper
!    use EquationOfState
!
    real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
!
    real, dimension(nx) :: rho_test, lnTT_test, TT_test
    real, dimension(1) :: fran_location
    integer, dimension(4) :: new_lmn
    integer :: ii
    integer :: deltar2, deltar2_test
    integer :: nfound=0, chosen_site
    integer :: m,n
!
!
    if(headtt) print*,'find_nearest_SNI: ENTER'
!    
    call random_number_wrapper(fran_location)   
    if (iproc==SNR%iproc) then
      deltar2=nx**2+ny**2+nx**2
      do n=n1,n2
      do m=m1,m2
        rho_test=exp(f(l1:l2,m,n,ilnrho))
        call eoscalc(f,nx,lnTT=lnTT_test)
        TT_test=exp(lnTT_test)
!   
        do ii=l1,l2
          if ((SNsite%rho .gt. rho_SN_min) .and. &
              (SNsite%TT .lt. TT_SN_max)) then
            deltar2_test=((ii-SNR%l)**2+(m-SNR%m)**2+(n-SNR%n)**2)
            if (deltar2_test .lt. deltar2) then
              nfound=1
              deltar2=deltar2_test
              new_lmn=(/ nfound, ii, m, n /)
            elseif (deltar2==deltar2_test) then
              nfound=nfound+1
            endif
          endif
        enddo
      enddo
      enddo
!
      if (nfound==0) then
        new_lmn=(/ nfound, SNR%l, SNR%m, SNR%n /)
      elseif (nfound.gt.1) then
        chosen_site=int(nfound*fran_location(1)+0.5) 
        nfound=0
  search_two: do n=n1,n2
        do m=m1,m2
          rho_test=exp(f(l1:l2,m,n,ilnrho))
          call eoscalc(f,nx,lnTT=lnTT_test)
          TT_test=exp(lnTT_test)
! 
          do ii=l1,l2
            if ((SNsite%rho .gt. rho_SN_min) .and. &
                (SNsite%TT .lt. TT_SN_max)) then
              deltar2_test=((ii-SNR%l)**2+(m-SNR%m)**2+(n-SNR%n)**2)
              if (deltar2==deltar2_test) then
                nfound=nfound+1
                if (nfound==chosen_site) then
                  new_lmn=(/ 1, SNR%l, SNR%m, SNR%n /)
                  exit search_two
                endif
              endif
            endif
          enddo
        enddo
        enddo search_two
      endif
    endif    
!
!
    call mpibcast_int(new_lmn,4,SNR%iproc)
    nfound=new_lmn(1)   
! 
    if (nfound.gt.0) then 
      SNR%l=new_lmn(2)
      SNR%m=new_lmn(3)
      SNR%n=new_lmn(4)
      call share_SN_parameters(f)
    endif
!
    endsubroutine find_nearest_SNI
!***********************************************************************
    subroutine share_SN_parameters(f)
!   
!   Handle common SN positioning processor communications
!
!
!   27-aug-2003/tony: coded
!    
    use Cdata
    use EquationOfState, only: eoscalc,ilnrho_lnTT
    use Mpicomm, only: mpibcast_int, mpibcast_real
!      
    real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
!    
    real, dimension(nx) :: lnTT
    real, dimension(5) :: fmpi5
    integer, dimension(4) :: impi4
!
!  Broadcast position to all processors from root;
!  also broadcast SNR%iproc, needed for later broadcast of SNsite%rho.
!
!
    impi4=(/ SNR%iproc, SNR%l, SNR%m, SNR%n /)
    call mpibcast_int(impi4,4)
    SNR%iproc=impi4(1)
    SNR%l=impi4(2)
    SNR%m=impi4(3)
    SNR%n=impi4(4)
!
!  With current SN scheme, we need rho at the SN location.
!
    if (iproc==SNR%iproc) then
      !SNsite%lnrho=alog(sum(exp(f(SNR%l:SNR%l+1,SNR%m:SNR%m+1,SNR%n:SNR%n+1,ilnrho)))/4.)
      !SNsite%ss=sum(f(SNR%l:SNR%l+1,SNR%m:SNR%m+1,SNR%n:SNR%n+1,iss))/4.
      SNsite%lnrho=f(SNR%l,SNR%m,SNR%n,ilnrho)
      m=SNR%m
      n=SNR%n
      call eoscalc(f,nx,lnTT=lnTT)
      SNsite%lnTT=lnTT(SNR%l-l1+1)
      SNR%x=0.; SNR%y=0.; SNR%z=0.
      if (nxgrid/=1) SNR%x=x(SNR%l)+dx/2.
      if (nygrid/=1) SNR%y=y(SNR%m)+dy/2.
      if (nzgrid/=1) SNR%z=z(SNR%n)+dz/2.
    print*, &
 'share_SN_parameters: (MY SNe) SNR%iproc,x_SN,y_SN,z_SN,SNR%l,SNR%m,SNR%n,SNsite%rho,SNsite%lnTT = ' &
          ,SNR%iproc,SNR%x,SNR%y,SNR%z,SNR%l,SNR%m,SNR%n,SNsite%rho,SNsite%lnTT
    else
      ! Better initialise these to something on the other processors
      SNsite%lnrho=0.
      SNsite%lnTT=0.
      SNR%x=0.
      SNR%y=0.
      SNR%z=0.
    endif
!
!  Broadcast to all processors.
!
    fmpi5=(/ SNR%x, SNR%y, SNR%z, SNsite%lnrho, SNsite%lnTT /)
    call mpibcast_real(fmpi5,5,SNR%iproc)
!
    SNR%x=fmpi5(1); SNR%y=fmpi5(2); SNR%z=fmpi5(3); 
    SNsite%lnrho=fmpi5(4); SNsite%lnTT=fmpi5(5)
!
    SNsite%rho=exp(SNsite%lnrho);
!
    call eoscalc(ilnrho_lnTT,SNsite%lnrho,SNsite%lnTT, &
                    yH=SNsite%yH,ss=SNsite%ss,ee=SNsite%ee)
    SNsite%TT=exp(SNsite%lnTT)
!
    if (lroot.and.ip<54) print*, &
 'share_SN_parameters: SNR%iproc,x_SN,y_SN,z_SN,SNR%l,SNR%m,SNR%n,SNsite%rho,SNsite%ss,SNsite%TT = ' &
          ,SNR%iproc,SNR%x,SNR%y,SNR%z,SNR%l,SNR%m,SNR%n,SNsite%rho,SNsite%ss,SNsite%TT
!
    endsubroutine share_SN_parameters
!***********************************************************************
    subroutine explode_SN(f,df,itype_SN)
      !
      !  Implement SN (of either type), at pre-calculated position
      !  (This can all be made more efficient, after debugging.)
      !
      !  ??-nov-02/grs : coded from GalaxyCode                        
      !  20-may-03/tony: pencil formulation and broken into subroutines
      !
      use Cdata, only: lroot,dx,dy,dz,dimensionality, &
                       ilnrho,iecr,ilntt,iyh, &
                       iuu, iux, iuy, iuz, &
                       it,t,datadir,m,n, pi
      use EquationOfState, only: ilnrho_ee, eoscalc, getdensity, eosperturb
      use Mpicomm, only: mpireduce_max, mpibcast_real, mpibcast_double, mpireduce_sum_double
!      use Sub, only: update_snaptime
!      use Slices, only: tvid, nvid
!
      real, intent(inout), dimension(mx,my,mz,mvar+maux) :: f
      real, intent(inout), dimension(mx,my,mz,mvar) :: df
      integer, intent(in) :: itype_SN
!
      double precision :: c_SN,cmass_SN,cvelocity_SN
      double precision :: mass_shell
      real :: rho_SN_lowest
      double precision :: width_energy, width_mass, width_velocity
      double precision :: cavity_depth, r_cavity, rhom, ekintot
      double precision ::  rhom_new, ekintot_new
      real :: rho_SN_new,lnrho_SN_new,yH_SN_new,lnTT_SN_new,ee_SN_new
      real :: TT_SN_new,dv, uu_sedov
!      
      double precision, dimension(nx) :: deltarho, deltaEE
      double precision, dimension(nx,3) :: deltauu
      real, dimension(1) :: fmpi1, fmpi1_tmp
      double precision, dimension(2) :: dmpi2, dmpi2_tmp
      real, dimension(nx) ::  lnrho, yH, lnTT, TT, rho_old, ee_old
      real, dimension(nx,3) :: uu
      character (len=4) :: ch
!
      logical :: lmove_mass=.false.
!!
!! For taking video slices centred on the most recent SN explosion.
!! DOES NOT WORK FOR PARALLEL JOBS - so outcommented
!!
!!      ! If slices are following the SN explosions then set the slice position
!!      ! and take a snapshot now.
!!      if (slice_position=='S') theN
!!        ix=SNR%l
!!        iy=SNR%m
!!        iz=SNR%n
!!        iz2=n2
!!        lwrite_slice_xy2=(SNR%ipz==ipz)
!!        lwrite_slice_xy=(SNR%ipz==ipz)
!!        lwrite_slice_xz=(SNR%ipy==ipy)
!!        lwrite_slice_yz=.true.
!!        ! And do a snapshot immediately! 
!!        call update_snaptime(trim(datadir)//'/tvid.dat',tvid,nvid,dvid,t,lvid,ch,ENUM=.false.)
!!      endif
!!      !
      !  identifier
      !
      if(lroot.and.ip<12) print*,'explode_SN: itype_SN=',itype_SN
      !
      width_energy   = SNR%radius*energy_width_ratio
      width_mass     = SNR%radius*mass_width_ratio
      width_velocity = SNR%radius*velocity_width_ratio
!
! Energy insertion normalization 
!
      if (thermal_profile=="gaussian3") then
        c_SN=ampl_SN/(cnorm_SN(dimensionality)*width_energy**dimensionality)
      elseif (thermal_profile=="gaussian2") then
        c_SN=ampl_SN/(cnorm_gaussian2_SN(dimensionality)*width_energy**dimensionality)
      elseif (thermal_profile=="gaussian") then
        c_SN=ampl_SN/(cnorm_gaussian_SN(dimensionality)*width_energy**dimensionality)
      elseif (thermal_profile=="quadratictanh") then
        c_SN=ampl_SN/(cnorm_para_SN(dimensionality)*width_energy**dimensionality)
      elseif (thermal_profile=="quartictanh") then
        c_SN=ampl_SN/(cnorm_quar_SN(dimensionality)*width_energy**dimensionality)
      elseif (thermal_profile=="tanh") then
        if (dimensionality==1) then
          c_SN=ampl_SN/( 2.*width_energy )
        elseif (dimensionality==2) then
          c_SN=ampl_SN/( pi*(width_energy)**2 )
        elseif (dimensionality==3) then
          c_SN=ampl_SN/( 4./3.*pi*(width_energy)**3 )
        endif
      endif

      if (lroot) print*,'explode_SN: c_SN         =',c_SN
!
! Mass insertion normalization 
!
      if (lSN_mass) then
        if (mass_profile=="gaussian3") then
          cmass_SN=mass_SN/(cnorm_SN(dimensionality)*width_mass**dimensionality)
        elseif (mass_profile=="quadratic") then
          cmass_SN=mass_SN/(cnorm_para_SN(dimensionality)*width_mass**dimensionality)
        elseif (mass_profile=="tanh") then
          if (dimensionality==1) then
            cmass_SN=mass_SN/( 2.*width_mass )
          elseif (dimensionality==2) then
            cmass_SN=mass_SN/( pi*(width_mass)**2 )
          elseif (dimensionality==3) then
            cmass_SN=mass_SN/( 4./3.*pi*(width_mass)**3 )
          endif
        endif
        if (lroot) print*,'explode_SN: cmass_SN     =',cmass_SN
      else 
        cmass_SN=0.
      endif

      !
      ! Calculate cross over point between mass addition and removal
      ! if mass movement is used
      !
      r_cavity=width_energy*(dimensionality*log(outer_shell_proportion/inner_shell_proportion)   &
                          /((1./inner_shell_proportion**6)      &
                            - (1./outer_shell_proportion**6))                       &
                        )**(1./6.)
      if (lroot) print*,'explode_SN: dimensionality,r_cavity',dimensionality,r_cavity
      if (lroot) print*,'explode_SN: shell_(inner, outer)_prop.=',inner_shell_proportion,outer_shell_proportion
!
!
      dv=1.
      if (nxgrid/=1) dv=dv*dx
      if (nygrid/=1) dv=dv*dy
      if (nzgrid/=1) dv=dv*dz
      if (lroot) print*,'explode_SN: dv=',dv

      if (lroot.and.ip<14) print*,'explode_SN: width_energy,c_SN,SNsite%rho=', width_energy,c_SN,SNsite%rho
        
      !
      !  Now deal with (if nec.) mass relocation
      !

      if (lroot) print*,'explode_SN: rho_new,SNsite%ee=',SNsite%rho,SNsite%ee
      ee_SN_new = (SNsite%ee+frac_eth*c_SN/(SNsite%rho+cmass_SN))
      if (lroot) print*,'explode_SN: rho_SN_new,ee_SN_new=',SNsite%rho+cmass_SN,ee_SN_new
      call eoscalc(ilnrho_ee,real(log(SNsite%rho+cmass_SN)),ee_SN_new, &
                              lnTT=lnTT_SN_new,yH=yH_SN_new)
      TT_SN_new=exp(lnTT_SN_new)
!
!  Calculate explosion site mean density
!
      call get_properties(f,SNR,rhom,ekintot)
      SNR%rhom=rhom
!
!  Calculate effective Sedov evolution time
!
      SNR%t_sedov = sqrt((SNR%radius/xsi_sedov)**5*SNR%rhom/ampl_SN)
      uu_sedov = 0.4*SNR%radius/SNR%t_sedov
!
! Velocity insertion normalization 
!
      if (lSN_velocity) then
        cvelocity_SN=uu_sedov  !/(cnorm_SN(dimensionality)*(width_SN*velocity_width_ratio)**dimensionality)
        !cvelocity_SN=velocity_SN !/(cnorm_SN(dimensionality)*(width_SN*velocity_width_ratio)**dimensionality)
        if (lroot) print*,'explode_SN: cvelocity_SN =',cvelocity_SN
      else 
        cvelocity_SN=0.
      endif

      if(lroot.and.ip<32) print*, &
         'explode_SN: SNsite%TT, TT_SN_new, TT_SN_min, SNsite%ee =', &
                                SNsite%TT,TT_SN_new,TT_SN_min, SNsite%ee
      if (lroot) print*,'explode_SN: yH_SN_new =',yH_SN_new

      if ((TT_SN_new < TT_SN_min).or.(mass_movement=='constant')) then
         if (lroot) print*,'explode_SN: SN will be too cold!'
         lmove_mass=.not.(mass_movement == 'off')
         ! lmove_mass=.false.  ! use to switch off for debug...

         ! The bit that BREAKS the pencil formulation...
         ! must know the total moved mass BEFORE attempting mass relocation 

         ! ASSUME: SN will fully ionize the gas at its centre
         if (lmove_mass) then
           if (lroot) print*,'explode_SN: moving mass to compensate.'
           call getdensity(real((SNsite%ee*SNsite%rho)+c_SN),TT_SN_min,1.,rho_SN_new)
           if (mass_movement=='rho-cavity') then 
             call get_lowest_rho(f,r_cavity,rho_SN_lowest)
             cavity_depth=SNsite%rho-rho_SN_new
             if (cavity_depth .gt. rho_SN_lowest-rho_min) then
               cavity_depth=rho_SN_lowest-rho_min
               if (cavity_depth .le. 0.) then
                 cavity_depth=0.
                 lmove_mass=.false.
               endif
               if (lroot) print*,"Reduced cavity from:,", &
                                  SNsite%rho-rho_SN_new," to: ", &
                                  cavity_depth
               rho_SN_new=SNsite%rho-cavity_depth
               lnrho_SN_new=log(rho_SN_new)
             endif
           elseif (mass_movement=='Galaxycode') then
             lnrho_SN_new=log(rho_SN_new-cmass_SN)
             cavity_depth=max(SNsite%lnrho-lnrho_SN_new,0.)
             cavity_profile="gaussian3log"
           elseif (mass_movement=='constant') then
             lnrho_SN_new=log(cmass_SN)
             cavity_depth=cmass_SN
             cavity_profile="tanh"
           endif
         endif
           
         if (lmove_mass) then
           ee_SN_new=(SNsite%ee*SNsite%rho+frac_eth*c_SN)/rho_SN_new

           call eoscalc(ilnrho_ee,lnrho_SN_new,ee_SN_new, &
                                 lnTT=lnTT_SN_new,yH=yH_SN_new)
           TT_SN_new=exp(lnTT_SN_new)

           if(lroot.and.ip<32) print*, &
              'explode_SN: Relocate mass... TT_SN_new, rho_SN_new=', &
                                                     TT_SN_new,rho_SN_new

           if (mass_movement=='rho_cavity') then 
             ! Do nowt.
           elseif (mass_movement=='Galaxycode') then 
             call calc_cavity_mass_lnrho(f,width_energy,cavity_depth,mass_shell)
             if (lroot.and.ip<32) &
               print*, 'explode_SN: mass_shell=',mass_shell
           elseif (mass_movement=='constant') then 
             call calc_cavity_mass_lnrho(f,width_energy,cavity_depth,mass_shell)
             if (lroot.and.ip<32) &
               print*, 'explode_SN: mass_shell=',mass_shell
           endif
         endif
      endif
      

      SNR%EE=0. 
      SNR%MM=0.
      !EE_SN2=0.
      do n=n1,n2
         do m=m1,m2
            
            ! Calculate the distances to the SN origin for all points
            ! in the current pencil and store in the dr2_SN global array
            call proximity_SN
            ! Get the old energy
            lnrho=f(l1:l2,m,n,ilnrho)
            rho_old=exp(lnrho) 
            deltarho=0.

            call eoscalc(f,nx,yH=yH,lnTT=lnTT,ee=ee_old)
            TT=exp(lnTT)

            ! Apply perturbations
            call injectenergy_SN(deltaEE,width_energy,c_SN,SNR%EE)
            if (lmove_mass) then
              if (mass_movement=='rho_cavity') then 
                if (lSN_mass) then
                  call make_cavity_rho(deltarho,width_energy,cavity_depth, &
                            cnorm_SN(dimensionality),SNR%MM)
                else
                  call make_cavity_rho(deltarho,width_energy,cavity_depth, &
                            cnorm_SN(dimensionality),SNR%MM)
                endif
                lnrho=log(rho_old(1:nx)+deltarho(1:nx))
              elseif (mass_movement=='Galaxycode') then 
                if (lSN_mass) then
                  call make_cavity_lnrho(lnrho,width_energy,cavity_depth, &
                      (mass_shell+mass_SN),cnorm_SN(dimensionality),SNR%MM)
                else
                  call make_cavity_lnrho(lnrho,width_energy,cavity_depth, &
                        mass_shell,cnorm_SN(dimensionality),SNR%MM)
                endif
              elseif (mass_movement=='constant') then 
                call make_cavity_lnrho(lnrho,width_mass,cmass_SN, &
                      mass_shell,cnorm_SN(dimensionality),SNR%MM)
              endif
            else
              if (lSN_mass) then
                call injectmass_SN(deltarho,width_mass,cmass_SN,SNR%MM)
                lnrho=log(rho_old(1:nx)+deltarho(1:nx))
              endif
            endif

            if (lSN_velocity) then
              uu=f(l1:l2,m,n,iux:iuz)
              call injectvelocity_SN(deltauu,width_velocity,cvelocity_SN)
              f(l1:l2,m,n,iux:iuz)=uu+deltauu
            endif
!  lnrho=log(rho_min)
!    where (rho_old(1:nx)+deltarho(1:nx) .gt. rho_min) lnrho=log(rho_old(1:nx)+deltarho(1:nx))
  
            TT=exp(lnTT)

            if (lcosmicray.and.lSN_ecr) then 
              f(l1:l2,m,n,iecr) = f(l1:l2,m,n,iecr) + (deltaEE * frac_ecr) 
            endif
!
!  Save changes to f-array
!
            f(l1:l2,m,n,ilnrho)=lnrho
            if (lSN_eth) then
              call eosperturb(f,nx,ee=real((ee_old*rho_old+deltaEE*frac_eth) &
                                                    /exp(lnrho)))
            endif
            call eoscalc(f,nx,lnTT=lnTT,yH=yH)
            lnTT=log(TT)
            if (lentropy.and.ilnTT.ne.0) f(l1:l2,m,n,ilnTT)=lnTT
            if (iyH.ne.0) f(l1:l2,m,n,iyH)=yH
!
       enddo
      enddo

      call get_properties(f,SNR,rhom_new,ekintot_new)
      print*,"TOTAL KINETIC ENERGY CHANGE:",ekintot_new-ekintot
!
!  Sum and share diagnostics etc. amongst processors
!
      dmpi2_tmp=(/ SNR%MM, SNR%EE /)
      call mpireduce_sum_double(dmpi2_tmp,dmpi2,2) 
      call mpibcast_double(dmpi2,2)
      SNR%MM=dmpi2(1)*dv
      SNR%EE=dmpi2(2)*dv
! Extra debug - no longer calculated 
!      EE2_SN=fmpi2(3)*dv; 
!print*,'EE2_SN = ',EE2_SN
!
      if (lroot.and.ip<20) print*, &
           'explode_SN: SNR%MM=',SNR%MM

      if (lroot) then
         open(1,file=trim(datadir)//'/sn_series.dat',position='append')
         !write(1,*)  &
        print*, 'explode_SN:    step, time = ', it,t
        print*, 'explode_SN:            dv = ', dv
        print*, 'explode_SN:       SN type = ', itype_SN 
        print*, 'explode_SN: proc, l, m, n = ', SNR%iproc, SNR%l,SNR%m,SNR%n
        print*, 'explode_SN:       x, y, z = ', SNR%x,SNR%y,SNR%z
        print*, 'explode_SN:       rho, TT = ', SNsite%rho,SNsite%TT
        print*, 'explode_SN:  Mean density = ', SNR%rhom
        print*, 'explode_SN:  Total energy = ', SNR%EE
        print*, 'explode_SN:    Total mass = ', SNR%MM
        print*, 'explode_SN:    Sedov time = ', SNR%t_sedov
        print*, 'explode_SN:    Shell velocity  = ', uu_sedov
         write(1,'(i10,e13.5,2i4,3i10,6e13.5)')  &
          it,t, &
          itype_SN, &
          SNR%iproc, SNR%l,SNR%m,SNR%n, &
          SNR%x,SNR%y,SNR%z, &
          SNsite%rho,SNsite%TT,SNR%EE
         close(1)
      endif

      if (lSNR_damping) then
        SNR%state=SNstate_damping
        SNR%t_damping=SNR_damping_time-SNR%t_sedov+t
      else
        SNR%state=SNstate_invalid
      endif

    endsubroutine explode_SN
!!***********************************************************************
!    subroutine calc_interstellar_SNR_smooth(f)
!!
!      use Cdata
!!
!      real, intent(inout), dimension(mx,my,mz,mvar+maux) :: f
!      real, dimension(nx) :: ee, rho, profile_smooth
!      real, dimension(nx) :: ee_new, rho_new
!      real, dimension(nx,3) :: uu, uu_new
!      real :: EE_tot, MM_tot, EKIN_tot
!      real :: EE_local, MM_local, EKIN_local
!      real, dimension(3) :: MOM_local, MOM_local_new
!      real :: EE_local_new, MM_local_new, EKIN_local_new
!!
!      if (SNR%state/=SNstate_damping) return
!!
!      EE_tot=0.
!      MM_tot=0.
!      MOM_tot=0.
!      do n=n1,n2; do m=m1,m2
!        call proximity_SN(dr2_SN)
!        profile_smooth=(exp(-(dr2_SN(1:nx)/width**2)))
!        call eoscalc(f,nx,ee=ee)
!        rho = exp(f(l1:l2,m,n,ilnrho))
!        rho_new=rho*
!        EE_tot=EE_tot+sum(rho*ee)
!        EKIN_tot=EKIN_tot+sum(sqrt(dot2(uu))*rho)
!        MOM_tot=MOM_tot+sum(sqrt(dot2(uu))*rho)
!        MM_tot=MM_tot+sum(rho)
!
!        where (dr2_SN .le. radius2)
!          EE_local=EE_local+sum(rho*ee)
!          MOM_local=MOM_local+sum(sqrt(dot2(uu))*rho)
!          MM_local=MM_local+sum(rho)
!        endwhere 
!      enddo; enddo
!print*,"SNR smooth: Total thermal (energy, momentum, mass)", EE_tot, MOM_tot, MM_tot
!
!
!!
!!  Have to divide by dxmin**2 to compensate for the * dxmin**2 in
!!  the shock code!
!!
!      penc=max(penc,SNR_damping*exp(-(dr2_SN_mx/SNR%radius**2))/dxmin**2)
!!
!    endsubroutine calc_interstellar_SNRdamping
!***********************************************************************
    subroutine calc_snr_damping_factor(f)
!
      use Cdata
      use Mpicomm
      use Sub
!
      real, intent(inout), dimension(mx,my,mz,mvar+maux) :: f
      integer, dimension(3) :: worst
      real, dimension(nx,3) :: r_vec, r_hat, uu
      real, dimension(nx) :: uur, ttc
      real, dimension(nx) :: r2, lnrho
      real :: radius2, fac, fac2
      integer :: l=0
!
      if (SNR%state/=SNstate_damping) return

      if ((t-SNR%t_damping)/SNR_damping_rate >= 2.5) then
        print*,"No more damping!!"
        SNR%state=SNstate_invalid
        return
      endif


!      SNR%damping_factor = SNR_damping * 0.5 * (1. - tanh((t-SNR%t_damping)/SNR_damping_rate))
      

      fac=0.
      radius2=SNR%radius**2
      do n=n1,n2; do m=m1,m2
        call proximity_SN()
        uu = f(l1:l2,m,n,iux:iuz)
        lnrho = f(l1:l2,m,n,ilnrho)
        r_vec=0.
        r_vec(:,1) = x(l1:l2) - SNR%x
        r_vec(:,2) = y(m) - SNR%y
        r_vec(:,3) = z(n) - SNR%z
        call dot2(r_vec,r2)
!        r_hat=r_vec/spread(sqrt(r2),2,3)
        call dot(r_vec,uu,uur)
        uur=-uur/r2
        where (r2 .gt. radius2)
          uur=0.
        endwhere 
        fac=max(fac,maxval(uur))
      enddo; enddo
!      fac=fac*2. 
     call mpiallreduce_max(fac,fac2)

!     print*,"Choose damping factor of:",fac2*2E-4
      SNR%damping_factor = fac2 * 2E-4 * (1. - tanh((t-SNR%t_damping)/SNR_damping_rate))
!

    endsubroutine calc_snr_damping_factor
!***********************************************************************
    subroutine calc_snr_unshock(penc)
!
      use Cdata
!
!      real, intent(inout), dimension(mx,my,mz,mvar+maux) :: f
!      integer, intent(in) :: ivar
      real, dimension(mx), intent(inout) :: penc
      real, dimension(mx) :: dr2_SN_mx
!
      if (SNR%state/=SNstate_damping) return
      call proximity_SN_mx(dr2_SN_mx)
      penc = penc*(1.-exp(-(dr2_SN_mx/SNR%radius**2)))
!
!
    endsubroutine calc_snr_unshock
!***********************************************************************
    subroutine calc_snr_damping(p)
!
      use Cdata
      use Sub, only: multsv, multsv_add, dot
!
!      real, intent(inout), dimension(mx,my,mz,mvar+maux) :: f
!      integer, intent(in) :: ivar
      type (pencil_case) :: p
      real, dimension(nx) :: profile
      real, dimension(nx,3) :: tmp, tmp2, gprofile
!
      if (SNR%state/=SNstate_damping) return

      if (lfirst) SNR%energy_loss=0.
      call proximity_SN
      profile = SNR%damping_factor*exp(-(dr2_SN/SNR%radius**2))
      gprofile = -SNR%damping_factor*spread(x(l1:l2)*(dr2_SN)**2/(SNR%radius**(2.5))*exp(-(dr2_SN/SNR%radius**2)),2,3)
      if (ldensity) then
        call multsv(p%divu,p%glnrho,tmp2)
        tmp=tmp2 + p%graddivu
        call multsv(profile,tmp,tmp2)
        call multsv_add(tmp2,p%divu,gprofile,tmp)
        if (lfirst.and.ldt) p%diffus_total=p%diffus_total+sum(profile)
        SNR%energy_loss=SNR%energy_loss-sum(profile*p%rho*p%divu**2)
        
        p%fvisc=p%fvisc+tmp
      endif
!
!
!  Have to divide by dxmin**2 to compensate for the * dxmin**2 in
!  the shock code!
!
!      penc=max(penc,SNR_damping*exp(-(dr2_SN_mx/SNR%radius**2))/dxmin**2)
!
    endsubroutine calc_snr_damping
!***********************************************************************
    subroutine calc_snr_damp_int(int_dt)
!
      use Cdata
      use Sub, only: multsv, multsv_add
      use EquationOfState, only: eoscalc, eosperturb
!
      real :: int_dt, dv
!
      if (SNR%state/=SNstate_damping) return
!
      dv=1.
      if (nxgrid/=1) dv=dv*dx
      if (nygrid/=1) dv=dv*dy
      if (nzgrid/=1) dv=dv*dz
!
      SNR%heat_energy=SNR%heat_energy+int_dt*SNR%energy_loss*dv*int_dt
!
    endsubroutine calc_snr_damp_int
!***********************************************************************
    subroutine calc_snr_damping_add_heat(f)
!
      use Cdata
      use Mpicomm
      use Sub, only: multsv, multsv_add
      use EquationOfState, only: eoscalc, eosperturb
!
      real, intent(inout), dimension(mx,my,mz,mvar+maux) :: f
      real, dimension(nx) :: profile, ee_old, rho, lnrho
      real :: factor
      real, dimension(1) :: fmpi, fmpi_tmp
!
      if (SNR%state/=SNstate_damping) return

      fmpi_tmp(1)=SNR%heat_energy 
      call mpireduce_sum(fmpi_tmp,fmpi,1) 
      call mpibcast_real(fmpi,1) 
      SNR%heat_energy=fmpi(1)

      factor=SNR%heat_energy/(cnorm_gaussian_SN(dimensionality)*SNR%radius**dimensionality)
      do n=n1,n2; do m=m1,m2
        call proximity_SN
        call eoscalc(f,nx,ee=ee_old)
        lnrho=f(l1:l2,m,n,ilnrho)
        rho=exp(lnrho)
        profile = factor*exp(-(dr2_SN/SNR%radius**2))
        call eosperturb(f,nx,ee=real((ee_old*rho+profile) &
                                                   /exp(lnrho)))
      enddo; enddo

      SNR%heat_energy=0.

      if (t >= SNR%t_damping) then
        SNR%state=SNstate_invalid
      endif
!
    endsubroutine calc_snr_damping_add_heat
!***********************************************************************
    subroutine get_properties(f,remnant,rhom,ekintot)
!
!  Calculate integral of mass cavity profile  
!
!  22-may-03/tony: coded
!
      use Cdata
      use Sub
      use Mpicomm
!
      real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
      type (SNRemnant) :: remnant
      double precision :: radius2
      double precision :: rhom, ekintot
      integer :: point_count
      real, dimension(nx) :: rho, u2
      integer, dimension(nx) :: mask
      double precision, dimension(3) :: tmp,tmp2
      real :: dv 
!
!  Obtain distance to SN
!
      dv=1.
      if (nxgrid/=1) dv=dv*dx
      if (nygrid/=1) dv=dv*dy
      if (nzgrid/=1) dv=dv*dz
!
!  Sum all points inside 1.5 * SNR radius and divide by number of points 
!
      radius2 = (remnant%radius)**2
      tmp=0.
      do n=n1,n2
         do m=m1,m2
            call proximity_SN
            mask=1
            rho=exp(f(l1:l2,m,n,ilnrho))
            call dot2(f(l1:l2,m,n,iuu:iuu+2),u2)
            tmp(3)=tmp(3)+sum(rho*u2)
            where (dr2_SN(1:nx) .gt. radius2) 
              rho(1:nx)=0.
              mask(1:nx)=0
            endwhere
            tmp(1)=tmp(1)+sum(rho)
            tmp(2)=tmp(2)+sum(mask)
         enddo
      enddo
!
      call mpireduce_sum_double(tmp,tmp2,3)
      call mpibcast_double(tmp2,3)
      ekintot=tmp2(3)*dv
      rhom=tmp2(1)/tmp2(2) 
!
    endsubroutine get_properties
!***********************************************************************
    subroutine get_lowest_rho(f,radius,rho_lowest)
!
!  Calculate integral of mass cavity profile  
!
!  22-may-03/tony: coded
!
      use Cdata, only: ilnrho,m,n
      use Mpicomm
!
      real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
      real, intent(out) :: rho_lowest
      real :: tmp
      double precision, intent(in) :: radius
      double precision :: radius2
      real, dimension(nx) :: rho
!
!  Find lowest rho value in the surronding cavity
!
      rho_lowest=1E10 
      radius2 = radius**2
      do n=n1,n2
         do m=m1,m2
            call proximity_SN
            rho=f(l1:l2,m,n,ilnrho)
            where (dr2_SN(1:nx) .gt. radius2) rho=1E10
            rho_lowest=min(rho_lowest,minval(rho(1:nx)))
         enddo
      enddo
!
      tmp=-exp(rho_lowest)
      call mpireduce_max(tmp,rho_lowest) 
      call mpibcast_real(rho_lowest,2)
!
    endsubroutine get_lowest_rho
!***********************************************************************
    subroutine proximity_SN()
!
!  Calculate pencil of distance to SN explosion site
!
!  20-may-03/tony: extracted from explode_SN code written by grs
!  22-may-03/tony: pencil formulation
!
      use Cparam
      use Cdata, only: Lx,Ly,Lz,x,y,z,lperi,m,n
!
      double precision,dimension(nx) :: dx_SN, dr_SN
      double precision :: dy_SN 
      double precision :: dz_SN 
      integer :: j
!
!  Obtain distance to SN
!
         dx_SN=x(l1:l2)-SNR%x
         if (lperi(1)) then
           where (dx_SN .gt. Lx/2.) dx_SN=dx_SN-Lx
           where (dx_SN .lt. -Lx/2.) dx_SN=dx_SN+Lx
         endif
!
         dy_SN=y(m)-SNR%y
         if (lperi(2)) then
           if (dy_SN .gt. Ly/2.) dy_SN=dy_SN-Ly
           if (dy_SN .lt. -Ly/2.) dy_SN=dy_SN+Ly
         endif
!
         dz_SN=z(n)-SNR%z
         if (lperi(3)) then
           if (dz_SN .gt. Lz/2.) dz_SN=dz_SN-Lz
           if (dz_SN .lt. -Lz/2.) dz_SN=dz_SN+Lz
         endif
!
         dr2_SN=dx_SN**2 + dy_SN**2 + dz_SN**2           
!
         if (lSN_velocity) then
           dr_SN=dsqrt(dr2_SN)
           outward_normal_SN(:,1)=dx_SN/dr_SN
           outward_normal_SN(:,2)=dy_SN/dr_SN
           outward_normal_SN(:,3)=dz_SN/dr_SN
         endif
!       
    endsubroutine proximity_SN
!***********************************************************************
    subroutine proximity_SN_mx(dr2_SN_mx)
!
!  Calculate pencil of distance to SN explosion site
!
!  20-may-03/tony: extracted from explode_SN code written by grs
!  22-may-03/tony: pencil formulation
!
      use Cparam
      use Cdata, only: Lx,Ly,Lz,x,y,z,lperi,m,n
!
      real,dimension(mx), intent(out) :: dr2_SN_mx
      real,dimension(mx) :: dx_SN, dr_SN
      real :: dy_SN 
      real :: dz_SN 
      integer :: j
!
!  Obtain distance to SN
!
         dx_SN=x-SNR%x
         if (lperi(1)) then
           where (dx_SN .gt. Lx/2.) dx_SN=dx_SN-Lx
           where (dx_SN .lt. -Lx/2.) dx_SN=dx_SN+Lx
         endif
!
         dy_SN=y(m)-SNR%y
         if (lperi(2)) then
           if (dy_SN .gt. Ly/2.) dy_SN=dy_SN-Ly
           if (dy_SN .lt. -Ly/2.) dy_SN=dy_SN+Ly
         endif
!
         dz_SN=z(n)-SNR%z
         if (lperi(3)) then
           if (dz_SN .gt. Lz/2.) dz_SN=dz_SN-Lz
           if (dz_SN .lt. -Lz/2.) dz_SN=dz_SN+Lz
         endif
!
         dr2_SN_mx=dx_SN**2 + dy_SN**2 + dz_SN**2           
!
    endsubroutine proximity_SN_mx
!***********************************************************************
    subroutine calc_cavity_mass_lnrho(f,width,depth,mass_removed)
!
!  Calculate integral of mass cavity profile  
!
!  22-may-03/tony: coded
!
      use Cdata, only: dx,dy,dz,ilnrho,m,n
      use Mpicomm, only: mpibcast_double, mpireduce_sum_double
!
      real, intent(in), dimension(mx,my,mz,mvar+maux) :: f
      double precision, intent(in) :: width, depth
      double precision, intent(out) :: mass_removed
      real, dimension(nx) :: lnrho, lnrho_old
      real, dimension(nx) :: rho
      double precision, dimension(1) :: dmpi1, dmpi1_tmp
      double precision, dimension(nx) :: profile_cavity
      real :: dv
!
!  Obtain distance to SN
!
      dv=1.
      if (nxgrid/=1) dv=dv*dx
      if (nygrid/=1) dv=dv*dy
      if (nzgrid/=1) dv=dv*dz

      !mass_start=0.
      !mass_end=0.
      mass_removed=0.
      do n=n1,n2
        do m=m1,m2
          call proximity_SN()

          lnrho_old=f(l1:l2,m,n,ilnrho)
          if (cavity_profile=="gaussian3log") then
            profile_cavity=(depth*exp(-(dr2_SN(1:nx)/width**2)**3))
            lnrho=lnrho_old - profile_cavity
            mass_removed=mass_removed+sum(exp(lnrho_old)-exp(lnrho))
          elseif (cavity_profile=="gaussian3") then
            profile_cavity=(depth*exp(-(dr2_SN(1:nx)/width**2)**3))
            lnrho=lnrho_old - profile_cavity
            mass_removed=mass_removed+sum(exp(lnrho_old)-exp(lnrho))
          elseif (cavity_profile=="gaussian2") then
            profile_cavity=(depth*exp(-(dr2_SN(1:nx)/width**2)**2))
            lnrho=lnrho_old - profile_cavity
            mass_removed=mass_removed+sum(exp(lnrho_old)-exp(lnrho))
          elseif (cavity_profile=="gaussian") then
            profile_cavity=(depth*exp(-(dr2_SN(1:nx)/width**2)))
            lnrho=lnrho_old - profile_cavity
            mass_removed=mass_removed+sum(exp(lnrho_old)-exp(lnrho))
          elseif (cavity_profile=="tanh") then
            profile_cavity=(1.-tanh( (width-sqrt(dr2_SN(1:nx)) ) *sigma_SN1 ))*0.5
            rho=exp(lnrho_old)*profile_cavity
            mass_removed=mass_removed+sum(exp(lnrho_old)-rho)
          endif

        enddo
      enddo
      dmpi1_tmp=(/ mass_removed /)
      call mpireduce_sum_double(dmpi1_tmp,dmpi1,1) 
      call mpibcast_double(dmpi1,1)
      mass_removed=dmpi1(1)*dv

    endsubroutine calc_cavity_mass_lnrho
!***********************************************************************
    subroutine make_cavity_rho(deltarho,width,depth, &
                             cnorm_dim,MMtot_SN)   
!
      use Cdata, only: dimensionality
!
      double precision, intent(in) :: width, depth, cnorm_dim
      double precision, intent(inout) :: MMtot_SN
      double precision, intent(out), dimension(nx) :: deltarho
!
      double precision, dimension(nx) :: profile_shell_outer,profile_shell_inner
      double precision :: width_shell_outer, width_shell_inner, c_shell
!
      width_shell_outer=outer_shell_proportion*width
      width_shell_inner=inner_shell_proportion*width
!
!      deltarho(1:nx) =  -depth*exp(-(dr2_SN(1:nx)/width**2)**3)
!      
      c_shell=-depth*cnorm_dim / ( (1. / width_shell_outer**dimensionality) -   &
                                   (1. / width_shell_inner**dimensionality))
!
!  add missing mass back into shell
!
      profile_shell_outer(1:nx)=                              &
           exp(-(dr2_SN(1:nx)/width_shell_outer**2)**3)       &
              /cnorm_dim/width_shell_outer**dimensionality
      profile_shell_inner(1:nx)=                              &
           exp(-(dr2_SN(1:nx)/width_shell_inner**2)**3)       &
              /cnorm_dim/width_shell_inner**dimensionality     
      deltarho(1:nx)= c_shell *      &
           (profile_shell_outer(1:nx) - profile_shell_inner(1:nx))
      MMtot_SN=MMtot_SN + sum(deltarho(1:nx))  
!
    endsubroutine make_cavity_rho
!***********************************************************************
    subroutine make_cavity_lnrho(lnrho,width,depth,mass_shell, &
                             cnorm_dim,MMtot_SN)   
!
      use Cdata, only: dimensionality
!
      double precision, intent(in) :: width, depth, mass_shell, cnorm_dim
      double precision, intent(inout) :: MMtot_SN
      real, intent(inout), dimension(nx) :: lnrho
!
      double precision, dimension(nx) :: profile_shell_outer,profile_shell_inner, profile_cavity
      double precision :: width_shell_outer, width_shell_inner, c_shell
      double precision :: mass_before, mass_after
!
      width_shell_outer=outer_shell_proportion*width
      width_shell_inner=inner_shell_proportion*width
!      
      c_shell=mass_shell / (cnorm_dim * (width_shell_outer**dimensionality - width_shell_inner**dimensionality))
!
      profile_shell_outer(1:nx)=exp(-(dr2_SN(1:nx)/width_shell_outer**2)**3)
      profile_shell_inner(1:nx)=exp(-(dr2_SN(1:nx)/width_shell_inner**2)**3)
!      
      mass_before=sum(exp(lnrho(1:nx)))  
      if (cavity_profile=="gaussian3log") then
        profile_cavity=(depth*exp(-(dr2_SN(1:nx)/width**2)**3))
        lnrho = lnrho(1:nx) - profile_cavity
        lnrho = log(exp(lnrho(1:nx)) + c_shell *    &
           (profile_shell_outer(1:nx) - profile_shell_inner(1:nx)))
      elseif (cavity_profile=="gaussian3") then
        profile_cavity=(depth*exp(-(dr2_SN(1:nx)/width**2)**3))
        lnrho = lnrho(1:nx) - profile_cavity
        lnrho = log(exp(lnrho(1:nx)) + c_shell *    &
           (profile_shell_outer(1:nx) - profile_shell_inner(1:nx)))
      elseif (cavity_profile=="gaussian2") then
        profile_cavity=(depth*exp(-(dr2_SN(1:nx)/width**2)**2))
        lnrho = lnrho(1:nx) - profile_cavity
        lnrho = log(exp(lnrho(1:nx)) + c_shell *    &
           (profile_shell_outer(1:nx) - profile_shell_inner(1:nx)))
      elseif (cavity_profile=="gaussian") then
        profile_cavity=(depth*exp(-(dr2_SN(1:nx)/width**2)))
        lnrho = lnrho(1:nx) - profile_cavity
        lnrho = log(exp(lnrho(1:nx)) + c_shell *    &
           (profile_shell_outer(1:nx) - profile_shell_inner(1:nx)))
      elseif (cavity_profile=="tanh") then
        profile_cavity=(1.-tanh( (width- sqrt(dr2_SN(1:nx)) )*sigma_SN1 ))*0.5
        lnrho = log(exp(lnrho(1:nx))*profile_cavity + c_shell *    &
           (  profile_shell_outer(1:nx) - profile_shell_inner(1:nx)) &
            + depth*(1.-tanh( (sqrt(dr2_SN(1:nx)) - width )*sigma_SN1 ))*0.5 &
           )
      endif
      mass_after=sum(exp(lnrho(1:nx)))  
      MMtot_SN=MMtot_SN + (mass_after-mass_before)
!      
    endsubroutine make_cavity_lnrho
!***********************************************************************
    subroutine injectenergy_SN(deltaEE,width,c_SN,EEtot_SN)
!      use Cdata
      !
      double precision, intent(in) :: width,c_SN
      double precision, intent(inout) :: EEtot_SN
      double precision, intent(out), dimension(nx) :: deltaEE
      !
      double precision, dimension(nx) :: profile_SN
      
      ! Whether mass moved or not, inject energy.
      !

      if (thermal_profile=="gaussian3") then
        profile_SN=exp(-(dr2_SN(1:nx)/width**2)**3)
      elseif (thermal_profile=="gaussian2") then
        profile_SN=exp(-(dr2_SN(1:nx)/width**2)**2)
      elseif (thermal_profile=="gaussian") then
        profile_SN=exp(-(dr2_SN(1:nx)/width**2))
      elseif (thermal_profile=="quadratictanh") then
        profile_SN=max(1d0-(dr2_SN(1:nx)/width**2),0d0) &
            *0.5*(1.-tanh((sqrt(dr2_SN)-width)*sigma_SN1))
      elseif (thermal_profile=="quartictanh") then
        profile_SN=max(1d0-(dr2_SN(1:nx)/width**2)**2,0d0) &
            *0.5*(1.-tanh((sqrt(dr2_SN)-width)*sigma_SN1))
      elseif (thermal_profile=="tanh") then
        profile_SN=(1.-tanh((sqrt(dr2_SN(1:nx))-width)*sigma_SN1))*0.5
      endif

      deltaEE(1:nx)=c_SN*profile_SN(1:nx) ! spatial energy density 
      EEtot_SN=EEtot_SN+sum(deltaEE(1:nx))   

    endsubroutine injectenergy_SN
!***********************************************************************
    subroutine injectmass_SN(deltarho,width,cmass_SN,MMtot_SN)
!      use Cdata
      !
      double precision, intent(in) :: width,cmass_SN
      double precision, intent(inout) :: MMtot_SN
      double precision, intent(out), dimension(nx) :: deltarho
      !
      double precision, dimension(nx) :: profile_SN
      
      ! Inject mass.
      !

      if (mass_profile=="gaussian3") then
        profile_SN=exp(-(dr2_SN(1:nx)/width**2)**3)
      elseif (mass_profile=="gaussian2") then
        profile_SN=exp(-(dr2_SN(1:nx)/width**2)**2)
      elseif (mass_profile=="gaussian") then
        profile_SN=exp(-(dr2_SN(1:nx)/width**2))
      elseif (mass_profile=="quadratic") then
        profile_SN=max(1d0-(dr2_SN(1:nx)/width**2),0D0)
      elseif (mass_profile=="tanh") then
      !
      !  This is normally handled in the mass movement section
      !
        profile_SN=(1.-tanh((sqrt(dr2_SN(1:nx))-width)*sigma_SN1))*0.5
      endif

      deltarho(1:nx)=cmass_SN*profile_SN(1:nx) ! spatial mass density 
      MMtot_SN=MMtot_SN+sum(deltarho(1:nx))   

    endsubroutine injectmass_SN
!***********************************************************************
    subroutine injectvelocity_SN(deltauu,width,cvelocity_SN)
!
      use Cdata
!
      double precision, intent(in) :: width,cvelocity_SN
      double precision, intent(out), dimension(nx,3) :: deltauu
!
      double precision, dimension(nx) :: profile_SN
      
      integer :: j
!
! Calculate deltauu
!
      if (velocity_profile=="quintictanh") then
        profile_SN=((sqrt(dr2_SN)/width)**5) &
                   *0.5*(1.-tanh((sqrt(dr2_SN)-(1.1*width))/(0.08*width)))
!
      elseif (velocity_profile=="lineartanh") then
        profile_SN=max(sqrt(dr2_SN)/width,1.)*0.5*(1.-tanh((sqrt(dr2_SN)-width)*sigma_SN1-2.))
!
      elseif (velocity_profile=="quadratictanh") then
        profile_SN=min((dr2_SN/(width**2)),0.5*(1.-tanh((sqrt(dr2_SN)-width)*sigma_SN1-2.)))
!
      elseif (velocity_profile=="gaussian3") then
        profile_SN=exp(-(dr2_SN(1:nx)/width**2)**3)
!
      elseif (velocity_profile=="cubictanh") then
        profile_SN=(sqrt(dr2_SN/width)**3) &
               * (1.-tanh((sqrt(dr2_SN)-(1.1*width))*sigma_SN1))
!
      elseif (velocity_profile=="gaussian3der") then
        profile_SN=(((sqrt(dr2_SN)**5)/width**6*(1./35.)) &
                  * exp(-(dr2_SN(1:nx)/width**2)**3)) ! &
                  !* (1.-tanh((sqrt(dr2_SN)-(1.1*width))/sigma_SN))
!
      elseif (velocity_profile=="quadratic") then
        profile_SN=dr2_SN(1:nx)/width**2
        where (dr2_SN.gt.(width**2)) profile_SN=0.
      endif
!
      do j=1,3
        deltauu(1:nx,j)=cvelocity_SN*profile_SN(1:nx)*outward_normal_SN(1:nx,j) ! spatial mass density 
      enddo
    endsubroutine injectvelocity_SN
!
endmodule interstellar
!***********************************************************************
