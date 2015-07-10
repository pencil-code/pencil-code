! $Id$
!
!  This modules contains the routines for SNe-driven ISM simulations.
!  Still in development.
!
!***************** AUTOMATIC CPARAM.INC GENERATION ***************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lheat_cool = .true.
!
! MAUX CONTRIBUTION 2
! COMMUNICATED AUXILIARIES 1
!
!*****************************************************************************
module Heat_Cool
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'heat_cool.h'
  include '../record_types.h'
!
!
!  integer :: icooling=0
!  integer :: inetcool=0
!
!  cp1=1/cp used to convert TT (and ss) into heat_cool code units
!  (useful, as many conditions conveniently expressed in terms of TT)
!  code units based on:
!    [length]  = 1kpc  = 3.09 10^21 cm
!    [time]    = 1Gyr  = 3.15 10^16 s             !no, on [u]=1km/s...
!    [rho]     =       = 1.00 10^-24 g/cm^3
!  Lambdaunits converts coolH into heat_cool code units.
!
  real :: unit_Lambda
!
!  Parameters for 'averaged'-SN heating
!
  real :: r_SNI_yrkpc2=4.E-6, r_SNII_yrkpc2=3.E-5
  real :: r_SNI=3.E+4, r_SNII=4.E+3
  real :: average_SNI_heating=0., average_SNII_heating=0.
!
!  Cooling timestep limiter coefficient
!  (This value 0.08 is overly restrictive. cdt_tauc=0.5 is a better value.)
!
  real :: cdt_tauc=0.5
!
!  Parameters for UV heating of Wolfire et al.
!
  real, parameter :: rhoUV_cgs=0.1
  real, parameter :: GammaUV_cgs=0.0147
  real, parameter :: TUV_cgs=7000., T0UV_cgs=20000., cUV_cgs=5.E-4
  real :: GammaUV=impossible, T0UV=impossible, cUV=impossible
!
!  04-jan-10/fred:
!  Amended cool dim from 7 to 11 to accomodate WSW dimension.
!  Appended null last term to all arrays for RBN and SS cooling
!
  double precision, dimension(11) :: coolH_cgs
  real, dimension(11) :: coolT_cgs
  real, dimension(11) :: coolB, lncoolH, lncoolT
  integer :: ncool
!
!  TT & z-dependent uv-heating profile
!
  real, dimension(mz) :: heat_z, zrho
  logical :: lthermal_hse=.false., lheatz_min=.true.
!
  real :: coolingfunction_scalefactor=1.
  real :: heatingfunction_scalefactor=1.
!
  real :: heating_rate = 0.015
  real :: heating_rate_code = impossible
!
!  Cooling & heating flags
!
  logical :: lsmooth_coolingfunc = .false.
  logical :: laverage_SN_heating = .false.
  logical :: lheating_UV         = .true.
!
!  Cooling time diagnostic
!
  integer :: idiag_taucmin=0
  integer :: idiag_Hmax=0
  integer :: idiag_Lamm=0
  integer :: idiag_nrhom=0
  integer :: idiag_rhoLm=0
  integer :: idiag_Gamm=0
!
!  Heating function, cooling function
!  method selection.
!
  character (len=labellen) :: cooling_select  = 'RBN'
  character (len=labellen) :: heating_select  = 'wolfire'
!
!  start parameters
!
  namelist /heat_cool_init_pars/ &
      init_heat_cool, &
      cooling_select, heating_select, heating_rate, &
      lthermal_hse, lheatz_min
!
! run parameters
!
  namelist /heat_cool_run_pars/ &
      cdt_tauc, &
      laverage_SN_heating, coolingfunction_scalefactor,  lforce_locate_SNI,&
      heatingfunction_scalefactor, &
      lheating_UV, cooling_select, heating_select, heating_rate, &
      lthermal_hse, lheatz_min
!
  contains
!
!***********************************************************************
    subroutine register_heat_cool()
!
!  19-nov-02/tony: coded
!
      use FArrayManager
!
      call farray_register_auxiliary('netcool',inetcool,communicated=.true.)
      call farray_register_auxiliary('cooling',icooling)
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Writing files for use with IDL
!
      if (naux+naux_com <  maux+maux_com) aux_var(aux_count)=',netcool $'
      if (naux+naux_com == maux+maux_com) aux_var(aux_count)=',netcool'
      aux_count=aux_count+1
      if (lroot) write(15,*) 'netcool = fltarr(mx,my,mz)*one'
!
    endsubroutine register_heat_cool
!***********************************************************************
    subroutine initialize_heat_cool(f)
!
!  Perform any post-parameter-read initialization eg. set derived
!  parameters
!
!  24-nov-02/tony: coded
!
!  read parameters from seed.dat and heat_cool.dat
!
      use General, only: random_seed_wrapper
      use Mpicomm, only: stop_it
      use EquationOfState, only: getmu
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real :: mu
!
      f(:,:,:,icooling)=0.0
      f(:,:,:,inetcool)=0.0
!
      call getmu(f,mu)
      if (unit_system=='cgs') then
!
!  this Lambda as such enters as n^2*Lambda(T) on the rhs of the
!  energy equation per unit volume
!
        unit_Lambda = unit_velocity**2 / unit_density / unit_time
      elseif (unit_system=='SI') then
        call stop_it('initialize_heat_cool: SI unit conversions not implemented')
      endif
      if (lroot) print*,'initialize_heat_cool: unit_Lambda',unit_Lambda
!
!  Mara: Initialize cooling parameters according to selection
!  Default selection 'RBN' Rosen & Bregman (1993)
!  Alternative selection 'SS' Sanchez-Salcedo et al. (2002)
!  Turn off cooling: cooling_select='off'
!  cooling_select in heat_cool_init_pars added
!
      if (cooling_select == 'RBN') then
        if (lroot) print*,'initialize_heat_cool: default RBN cooling function'
        coolT_cgs = (/  10.0,         &
                        2000.0,       &
                        8000.0,       &
                        1E5,          &
                        1E6,          &
                        1E17,         &
                        tiny(0E0),    &
                        tiny(0E0),    &
                        tiny(0E0),    &
                        tiny(0E0),    &
                        tiny(0E0) /)
        coolH_cgs = (/  2.2380D-32,   &
                        1.0012D-30,   &
                        4.6240D-36,   &
                        1.7783524D-18,&
                        2.238814D-25, &
                        tiny(0D0),    &
                        tiny(0D0),    &
                        tiny(0D0),    &
                        tiny(0D0),    &
                        tiny(0D0),    &
                        tiny(0D0) /)  / ( m_p_cgs )**2
        coolB =     (/  2.0,          &
                        1.5,          &
                        2.867,        &
                       -0.65,         &
                        0.5,          &
                        tiny(0.),     &
                        tiny(0.),     &
                        tiny(0.),     &
                        tiny(0.),     &
                        tiny(0.),     &
                        tiny(0.)  /)
        ncool=5
      else if (cooling_select == 'SS') then
!
!  These are the SS et al (2002) coefficients multiplied by m_proton**2
!
        if (lroot) print*,'initialize_heat_cool: revised SS cooling fct'
        coolT_cgs = (/  10.0,                 &
                        141.0,                &
                        313.0,                &
                        6102.0,               &
                        1E5,                  &
                        4E7,                  &
                        1E17,                 &
                        tiny(0E0),            &
                        tiny(0E0),            &
                        tiny(0E0),            &
                        tiny(0.0)  /)
        coolH_cgs = (/  3.703109927416290D16, &
                        9.455658188464892D18, &
                        1.185035244783337D20, &
                        1.9994576479D8,       &
                        7.96D29,              &
                        1.440602814622207D21, &
                        tiny(0D0),            &
                        tiny(0D0),            &
                        tiny(0D0),            &
                        tiny(0D0),            &
                        tiny(0D0)  /)
        coolB =     (/  2.12,                 &
                        1.0,                  &
                        0.56,                 &
                        3.67,                 &
                        -0.65,                &
                        0.5,                  &
                        tiny(0.),             &
                        tiny(0.),             &
                        tiny(0.),             &
                        tiny(0.),             &
                        tiny(0.)   /)
        ncool=6
!
!  26-Jan-10/fred
!  Combines Sanchez-Salcedo (2002) with Slyz et al (2005) above 1E5K
!  as Gressel simulation (2008) with constants revised for continuity
!
      else if (cooling_select == 'WSW') then
        if (lroot) print*,'initialize_heat_cool: WSW cooling fct'
        coolT_cgs = (/  10.0,                 &
                        141.0,                &
                        313.0,                &
                        6102.0,               &
                        1E5,                  &
                        2.88E5,               &
                        4.73E5,               &
                        2.11E6,               &
                        3.98E6,               &
                        2.0E7,                &
                        1.0E17      /)
        coolH_cgs = (/  3.703109927416290D16, &
                        9.455658188464892D18, &
                        1.185035244783337D20, &
                        1.102120336D10,       &
                        1.236602671D27,       &
                        2.390722374D42,       &
                        4.003272698D26,       &
                        1.527286104D44,       &
                        1.608087849D22,       &
                        9.228575532D20,       &
                        tiny(0D0) /)
        coolB =     (/  2.12,                 &
                        1.0,                  &
                        0.56,                 &
                        3.21,                 &
                       -0.20,                 &
                       -3.0,                  &
                       -0.22,                 &
                       -3.00,                 &
                        0.33,                 &
                        0.50,                 &
                        tiny(0.) /)
        ncool=10
      else if (cooling_select == 'off') then
        if (lroot) print*,'initialize_heat_cool: no cooling applied'
        coolT_cgs=tiny(0.0)
        coolH_cgs=tiny(0.0)
        coolB=tiny(0.)
      endif
!
! BEGIN TEMPORARY
      if (any(coolH_cgs(1:ncool) == 0) &
        .or. any(coolT_cgs(1:ncool+1) == 0)) then
        call fatal_error('initialize_heat_cool', &
        'Calculating lncoolH and lncoolT: One of the cooling coefficient is zero')
      endif
! END TEMPORARY
      lncoolH(1:ncool) = real(log(coolH_cgs(1:ncool)) - log(unit_Lambda) &
                              + log(unit_temperature**coolB(1:ncool)) &
                              + log(coolingfunction_scalefactor))
      lncoolT(1:ncool+1) = real(log(coolT_cgs(1:ncool+1) / unit_temperature))
!
      heating_rate_code=heating_rate*real(unit_length/unit_velocity**3)
!
      if (unit_system=='cgs') then
        T0UV=T0UV_cgs / unit_temperature
        cUV=cUV_cgs * unit_temperature
        if (GammaUV==impossible) &
            GammaUV=GammaUV_cgs * real(unit_length/unit_velocity**3)
      else
        call stop_it('initialize_heat_cool: SI unit conversions not implemented')
      endif
!
      if (heating_select == 'thermal-hs') then
        call thermal_hs(f,zrho)
        call heat_heat_cool(f,heat_z,zrho)
      endif
!
      average_SNI_heating = &
          r_SNI *ampl_SN/(sqrt(pi)*h_SNI )*heatingfunction_scalefactor
      average_SNII_heating= &
          r_SNII*ampl_SN/(sqrt(pi)*h_SNII)*heatingfunction_scalefactor
      if (lroot) print*,'initialize_heat_cool: t_interval_SNI =', &
          t_interval_SNI,Lxyz(1),Lxyz(2),SNI_area_rate
!
      if (lroot .and. (ip<14)) then
        print*,'initialize_heat_cool: nseed,seed',nseed,seed(1:nseed)
        print*,'initialize_heat_cool: finished'
      endif
!
!  Write unit_Lambda to pc_constants file
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,'(a,1pd26.16)') 'unit_Lambda=',unit_Lambda
        close (1)
      endif
!
    endsubroutine initialize_heat_cool
!*****************************************************************************
    subroutine rprint_heat_cool(lreset,lwrite)
!
!  Reads and registers print parameters relevant to heat_cool
!
!  01-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset
!  (This needs to be consistent with what is defined above!)
!
      if (lreset) then
        SNR_damping=0.
        idiag_taucmin=0
        idiag_Hmax=0
        idiag_Lamm=0
        idiag_nrhom=0
        idiag_rhoLm=0
        idiag_Gamm=0
!
     endif
!
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
!  Write column in which each heat_cool variable is stored
!
      if (lwr) then
        write(3,*) 'icooling=',icooling
        write(3,*) 'inetcool=',inetcool
      endif
!
    endsubroutine rprint_heat_cool
!*****************************************************************************
    subroutine get_slices_heat_cool(f,slices)
!
!  Write slices for animation of Interstellar variables.
!
!  26-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Cooling profile
!
        case ('cooling')
          slices%yz=f(ix_loc,m1:m2 ,n1:n2  ,icooling)
          slices%xz=f(l1:l2 ,iy_loc,n1:n2  ,icooling)
          slices%xy=f(l1:l2 ,m1:m2 ,iz_loc ,icooling)
          slices%xy2=f(l1:l2,m1:m2 ,iz2_loc,icooling)
          slices%ready = .true.
        case ('netcool')
          slices%yz=f(ix_loc,m1:m2 ,n1:n2  ,inetcool)
          slices%xz=f(l1:l2 ,iy_loc,n1:n2  ,inetcool)
          slices%xy=f(l1:l2 ,m1:m2 ,iz_loc ,inetcool)
          slices%xy2=f(l1:l2,m1:m2 ,iz2_loc,inetcool)
          slices%ready = .true.
!
      endselect
!
    endsubroutine get_slices_heat_cool
!***********************************************************************
    subroutine read_heat_cool_init_pars(iostat)
!
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include '../parallel_unit.h'
!
      read(parallel_unit, NML=heat_cool_init_pars, IOSTAT=iostat)
!
    endsubroutine read_heat_cool_init_pars
!***********************************************************************
    subroutine write_heat_cool_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=heat_cool_init_pars)
!
    endsubroutine write_heat_cool_init_pars
!***********************************************************************
    subroutine read_heat_cool_run_pars(iostat)
!
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include '../parallel_unit.h'
!
      read(parallel_unit, NML=heat_cool_run_pars, IOSTAT=iostat)
!
    endsubroutine read_heat_cool_run_pars
!***********************************************************************
    subroutine write_heat_cool_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=heat_cool_run_pars)
!
    endsubroutine write_heat_cool_run_pars
!***********************************************************************
    subroutine init_heat_cool(f)
!
!  Initialise some explosions etc.
!  24-nov-2002/tony: coded
!
      use General, only: itoa
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      logical :: lnothing=.true.
      character (len=intlen) :: iinit_str
      integer :: i,j,iSNR
!
      intent(inout) :: f
!
      do j=1,ninit
!
      if (initheat_cool(j)/='nothing') then
!
      lnothing=.false.
      iinit_str=itoa(j)
!
!  Select different initial conditions
!
      select case (initheat_cool(j))
!
!        case ('single')
!
        case default
!
!  Catch unknown values
!
          write(unit=errormsg,fmt=*) 'No such value for initheat_cool(' &
                           //trim(iinit_str)//'): ',trim(initheat_cool(j))
          call fatal_error('init_heat_cool',errormsg)
!
      endselect
!
      if (lroot) print*,'init_heat_cool: initheat_cool(' &
                        //trim(iinit_str)//') = ',trim(initheat_cool(j))
      endif
!
      enddo
!
      if (lnothing.and.lroot) print*,'init_heat_cool: nothing'
!
      call tidy_SNRs
!
    endsubroutine init_heat_cool
!*****************************************************************************
    subroutine pencil_criteria_heat_cool()
!
!  All pencils that the Interstellar module depends on are specified here.
!
!  26-mar-05/tony: coded
!  11-mar-06/axel: added idiag_nrhom
!
      lpenc_requested(i_ee)=.true.
      lpenc_requested(i_lnrho)=.true.
      lpenc_requested(i_lnTT)=.true.
      lpenc_requested(i_TT1)=.true.
      lpenc_requested(i_rho1)=.true.
!
!  Diagnostic pencils
!
!  AB:
!      if (idiag_nrhom/=0) lpenc_diagnos(i_rho)=.true.
!
    endsubroutine pencil_criteria_heat_cool
!*****************************************************************************
    subroutine heat_cool_before_boundary(f)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  01-aug-06/tony: coded
!
      use Diagnostics, only: max_mn_name, sum_mn_name
      use EquationOfState, only: gamma, gamma1, eoscalc
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real, dimension (nx) :: heat,cool,lnTT, lnrho
      real, dimension (nx) :: damp_profile
      real :: minqty
      integer :: i,iSNR
!
      if (.not.(lcooltime_smooth.or.lcooltime_despike)) return
!
!  Identifier
!
      if (headtt) print*,'heat_cool_before_boundary: ENTER'
!
!  Precalculate radiative cooling function
!
    do n=n1,n2
    do m=m1,m2
      if (ldensity_nolog) then
        lnrho = log(f(l1:l2,m,n,irho))
      else
        lnrho = f(l1:l2,m,n,ilnrho)
      endif
      call eoscalc(f,nx,lnTT=lnTT)
!
      call calc_cool_func(cool,lnTT,lnrho)
      call calc_heat(heat,lnTT)
!
      if (nSNR==0) then
!
!  05-sep-10/fred
!  NB The applied net heating/cooling: (heat-cool)/temp is stored in
!  inetcool. The radiative cooling rho*Lambda is stored in icooling.
!  Both are diagnostic.
!
        if (ltemperature) then
          f(l1:l2,m,n,inetcool)=exp(-lnTT)*(heat-cool)*gamma
        elseif (pretend_lnTT) then
          f(l1:l2,m,n,inetcool)=exp(-lnTT)*(heat-cool)*gamma
        else
          f(l1:l2,m,n,inetcool)=exp(-lnTT)*(heat-cool)
        endif
      else
        damp_profile=1.
        do i=1,nSNR
          iSNR=SNR_index(i)
          if (SNRs(iSNR)%state==SNstate_damping) then
            call proximity_SN(SNRs(iSNR))
            minqty=0.5*(1.+tanh((t-SNRs(iSNR)%t_damping)/SNR_damping_rate))
            damp_profile = damp_profile * ( (1.-minqty) * 0.5 * (1.+ &
                tanh((sqrt(dr2_SN)-(SNRs(iSNR)%radius*2.))*sigma_SN1-2.)) &
                + minqty )
          endif
        enddo
        if (ltemperature) then
          f(l1:l2,m,n,inetcool)=exp(-lnTT)*(heat-cool)*gamma*damp_profile
        elseif (pretend_lnTT) then
          f(l1:l2,m,n,inetcool)=exp(-lnTT)*(heat-cool)*gamma*damp_profile
        else
          f(l1:l2,m,n,inetcool)=exp(-lnTT)*(heat-cool)*damp_profile
        endif
      endif
    enddo
    enddo
!
!
    endsubroutine heat_cool_before_boundary
!*****************************************************************************
    subroutine calc_heat_cool(f,df,p,Hmax)
!
!  This routine calculates and applies the optically thin cooling function
!  together with UV heating.
!
!  We may want to move it to the entropy module for good, because its use
!  is not restricted to heat_cool runs (could be used for solar corona).
!  Also, it doesn't pose an extra load on memory usage or compile time.
!  (We should allow that UV heating can be turned off; so rhoUV should
!  be made an input parameter.)
!
!  19-nov-02/graeme: adapted from calc_heat_cool
!  10-aug-03/axel: TT is used as input
!   3-apr-06/axel: add ltemperature switch
!  05-sep-10/fred: added TT_cutoff option, comments, revised diagnostics.
!
      use Diagnostics, only: max_mn_name, sum_mn_name
      use EquationOfState, only: gamma, gamma1
      use Sub, only: smooth_kernel, despike, dot2
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case) :: p
!
      real, dimension (nx), intent(inout) :: Hmax
      real, dimension (nx) :: heat,cool,heatcool,netheat,netcool
      real :: minqty
      integer :: i, iSNR
!
!  Identifier
!
      if (headtt) print*,'calc_heat_cool_heat_cool: ENTER'
!
!  05-sep-10/fred
!  NB redistributing the applied cooling/heating using smooth_kernel or
!  despike was found to add to the thermal instability at low temperatures.
!  Since heatcool is divided by TT for the entropy equation, heatcool is
!  shared with neighbours with low temperatures ~0.001 they are rapidly
!  amplified producing both superfluids and hyper-heating to crash the code.
!  I therefore recommend not using them at all.
!
      call calc_cool_func(cool,p%lnTT,p%lnrho)
!
      call calc_heat(heat,p%lnTT)
!
!  For clarity we have constructed the rhs in erg/s/g [=T*Ds/Dt] so therefore
!  we now need to multiply by TT1. At very low temperatures this becomes
!  unstable, so limit denominator to TT_cutoff to prevent unresolvable
!  supercooling or heating spikes.
!  This may have been caused by the use of smoothing and despiking discussed
!  above, however am leaving this switch as an option until a long enough run
!  has demonstrated whether or not it is required. Fred
!
      if (ltemperature) then
        heatcool=p%TT1*(heat-cool)*gamma
        if (lTT_cutoff) then
          where (p%TT1>TT_cutoff1) heatcool=TT_cutoff1*(heat-cool)*gamma
        endif
      elseif (pretend_lnTT) then
        heatcool=p%TT1*(heat-cool)*gamma
        if (lTT_cutoff) then
          where (p%TT1>TT_cutoff1) heatcool=TT_cutoff1*(heat-cool)*gamma
        endif
      else
        heatcool=p%TT1*(heat-cool)
        if (lTT_cutoff) then
          where (p%TT1>TT_cutoff1) heatcool=TT_cutoff1*(heat-cool)
        endif
      endif
!
!  Save result in diagnostic aux variable
!  cool=rho*Lambda, heatcool=(Gamma-rho*Lambda)/TT
!
      f(l1:l2,m,n,icooling)=cool
      f(l1:l2,m,n,inetcool)=heatcool
!
!  Prepare diagnostic output
!  Since these variables are divided by Temp when applied it is useful to
!  monitor the actual applied values for diagnostics so TT1 included.
!
      if (ldiagnos) then
        if (idiag_Hmax/=0) then
          netheat=heatcool
          where (heatcool<0.0) netheat=0.0
          call max_mn_name(netheat/p%ee,idiag_Hmax)
        endif
        if (idiag_taucmin/=0) then
          netcool=-heatcool
          where (heatcool>=0.0) netcool=1.0
          call max_mn_name(netcool/p%ee,idiag_taucmin,lreciprocal=.true.)
        endif
        if (idiag_Lamm/=0) &
          call sum_mn_name(p%rho1*cool*p%TT1,idiag_Lamm)
        if (idiag_nrhom/=0) &
          call sum_mn_name(cool/p%ee,idiag_nrhom)
!  --       call sum_mn_name(cool*p%rho/p%ee,idiag_nrhom)
!  AB: the factor rho is already included in cool, so cool=rho*Lambda
        if (idiag_rhoLm/=0) &
          call sum_mn_name(p%TT1*cool,idiag_rhoLm)
        if (idiag_Gamm/=0) &
          call sum_mn_name(p%TT1*heat,idiag_Gamm)
      endif
!
!  Limit timestep by the cooling time (having subtracted any heating)
!  dt1_max=max(dt1_max,cdt_tauc*(cool)/ee,cdt_tauc*(heat)/ee)
!
      if (lfirst.and.ldt) then
        dt1_max=max(dt1_max,(-heatcool)/(p%ee*cdt_tauc))
        where (heatcool>0.0) Hmax=Hmax+heatcool
        dt1_max=max(dt1_max,Hmax/(p%ee*cdt_tauc))
      endif
!
!  Apply heating/cooling to temperature/entropy variable
!
      if (ltemperature) then
        df(l1:l2,m,n,ilnTT)=df(l1:l2,m,n,ilnTT)+heatcool
      else
        df(l1:l2,m,n,iss)=df(l1:l2,m,n,iss)+heatcool
      endif
!
    endsubroutine calc_heat_cool_heat_cool
!*****************************************************************************
    subroutine calc_cool_func(cool,lnTT,lnrho)
!
!  This routine calculates the temperature dependent radiative cooling.
!  Applies Rosen et al., ApJ, 413, 137, 1993 ('RBN') OR
!  Sanchez-Salcedo et al. ApJ, 577, 768, 2002 ('SS') OR
!  Slyz et al MNRAS, 356 2005 ('WSW') fit of Wolfire with Sarazin & White
!  ApJ, 443:152-168, 1985 and ApJ, 320:32-48, 1987
!
!
!  Cooling is Lambda*rho^2, with (eq 7) & Lambda=coolH(i)*TT**coolB(i),
!  for coolT(i) <= TT < coolT(i+1). Nb: our coefficients coolH(i) differ from
!  those in Rosen et al. by a factor (mu mp)^2, with mu=1.2, since Rosen works
!  in number density, n. (their cooling = Lambda*n^2,  rho=mu mp n)
!  The factor Lambdaunits converts from cgs units to code units.
!
!  [Currently, coolT(1) is not modified, but this may be necessary
!  to avoid creating gas too cold to resolve.]
!
      real, dimension (nx), intent(out) :: cool
      real, dimension (nx), intent(in) :: lnTT, lnrho
      integer :: i
!
      cool=0.0
      cool_loop: do i=1,ncool
        if (lncoolT(i) >= lncoolT(i+1)) exit cool_loop
        where (lncoolT(i) <= lnTT .and. lnTT < lncoolT(i+1))
          cool=cool+exp(lncoolH(i)+lnrho+lnTT*coolB(i))
        endwhere
      enddo cool_loop
    endsubroutine calc_cool_func
!*****************************************************************************
    subroutine calc_heat(heat,lnTT)
!
!  This routine adds UV heating, cf. Wolfire et al., ApJ, 443, 152, 1995
!  with the values above, this gives about 0.012 erg/g/s (T < ~1.E4 K)
!  Nb: need rho0 from density_[init/run]_pars, to implement the arm/interarm
!  scaling.
!
!  Control with heating_select in heat_cool_init_pars/run_pars.
!  Default heating_rate GammaUV = 0.015.
!
      real, dimension (nx), intent(out) :: heat
      real, dimension (nx), intent(in) :: lnTT
!
!  Constant heating with a rate heating_rate[erg/g/s].
!
      if (heating_select == 'cst') then
         heat = heating_rate_code
      else if (heating_select == 'wolfire') then
        heat(1:nx)=GammaUV*0.5*(1.0+tanh(cUV*(T0UV-exp(lnTT))))
      else if (heating_select == 'wolfire_min') then
        heat(1:nx)=GammaUV*0.5*(1.0+tanh(cUV*(T0UV-exp(lnTT))))
        heat = max(heat,heating_rate_code)
!
!  If using thermal-hs in initial entropy this must also be specified for
!  thermal equilibrium and applies for vertically stratified density supported
!  by vertical gravity profile 'Ferriere'.
!
      else if (heating_select == 'thermal-hs') then
        heat(1:nx) = heat_z(n)*0.5*(1.0+tanh(cUV*(T0UV-exp(lnTT))))
      else if (heating_select == 'off') then
        heat = 0.
      endif
!
    endsubroutine calc_heat
!*****************************************************************************
 endmodule Heat_Cool
