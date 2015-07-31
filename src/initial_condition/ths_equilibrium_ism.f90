! $Id$
!
!  This module provide a way for users to specify custom initial
!  conditions.
!
!  The module provides a set of standard hooks into the Pencil Code
!  and currently allows the following customizations:
!
!   Description                               | Relevant function call
!  ------------------------------------------------------------------------
!   Initial condition registration            | register_initial_condition
!     (pre parameter read)                    |
!   Initial condition initialization          | initialize_initial_condition
!     (post parameter read)                   |
!                                             |
!   Initial condition for momentum            | initial_condition_uu
!   Initial condition for density             | initial_condition_lnrho
!   Initial condition for entropy             | initial_condition_ss
!   Initial condition for magnetic potential  | initial_condition_aa
!                                             |
!   Initial condition for all in one call     | initial_condition_all
!     (called last)                           |
!
!   And a similar subroutine for each module with an "init_XXX" call.
!   The subroutines are organized IN THE SAME ORDER THAT THEY ARE CALLED.
!   First uu, then lnrho, then ss, then aa, and so on.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: linitial_condition = .true.
!
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!    linitial_condition = .true.
! to enable use of custom initial conditions.
!
! The rest of this file may be used as a template for your own initial
! conditions. Simply fill out the prototypes for the features you want
! to use.
!
! Save the file with a meaningful name, e.g. mhs_equilibrium.f90, and
! place it in the $PENCIL_HOME/src/initial_condition directory. This
! path has been created to allow users to optionally check their
! contributions in to the Pencil Code SVN repository. This may be
! useful if you are working on/using an initial condition with
! somebody else or may require some assistance from one from the main
! Pencil Code team. HOWEVER, less general initial conditions should
! not go here (see below).
!
! You can also place initial condition files directly in the run
! directory. Simply create the folder 'initial_condition' at the same
! level as the *.in files and place an initial condition file there.
! With pc_setupsrc this file is linked automatically into the local
! src directory. This is the preferred method for initial conditions
! that are not very general.
!
! To use your additional initial condition code, edit the
! Makefile.local in the src directory under the run directory in which
! you wish to use your initial condition. Add a line that says e.g.
!
!    INITIAL_CONDITION =   initial_condition/mhs_equilibrium
!
! Here mhs_equilibrium is replaced by the filename of your new file,
! not including the .f90 extension.
!
! This module is based on Tony's special module.
!
module InitialCondition
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../initial_condition.h'
!
!  Gravitational parameters from Kuijken & Gilmore
!
  real :: a_S, a_D
  real, parameter :: a_S_cgs=4.4e-9, a_D_cgs=1.7e-9
  real :: z_S, z_D
  real, parameter :: z_S_cgs=6.172e20, z_D_cgs=3.086e21
!
!  Observed number density per cubic cm parameters from Dickey & Lockman
!  Includes neutral hydrogen and warm ionized hydrogen plus helium 
!  proportionately. Multiply by m_u_cgs for gas density
!
  real, parameter, dimension(5) :: nfraction_cgs = (/0.6541, 0.1775, 0.1028, 0.0245, 0.0411/) ! particles per cm cubed 
  real, parameter, dimension(5) :: hscale_cgs = (/3.919e20, 9.813e20, 1.244e21, 2.160e20, 2.777e21/) ! scale height in cm
  real, dimension(5) :: rho_fraction, hscale 
!
!  Parameters for UV heating of Wolfire et al.
!
  real, parameter :: rhoUV_cgs=0.1
  real, parameter :: GammaUV_cgs=0.0147
  real, parameter :: TUV_cgs=7000., T0UV_cgs=20000., cUV_cgs=5.E-4
  real :: GammaUV=impossible, T0UV=impossible, cUV=impossible
  real :: unit_Lambda, unit_Gamma, coolingfunction_scalefactor=1.
!
!  04-jan-10/fred:
!  Amended cool dim from 7 to 11 to accomodate WSW dimension.
!  Appended null last term to all arrays for RBN and SS cooling
!
  double precision, dimension(11) :: lncoolH, coolH_cgs
  real, dimension(11) :: coolT_cgs
  real, dimension(11) :: coolB, lncoolT
  integer :: ncool
!
!  Heating function, cooling function and mass movement
!  method selection.
!
  character (len=labellen) :: cooling_select  = 'WSW'
  character (len=labellen) :: heating_select  = 'wolfire'
!
!
  real, parameter :: T0hs_cgs=1E3, zT_cgs = 2.785E21
  real :: T0hs=impossible, zT=impossible, rho_max, const_T_frac=1.25
!
!
!  TT & z-dependent uv-heating profile
!
  real, dimension(mz) :: rho, TT, rho_obs, drhodz, d2rhodz2, dTdz, d2Tdz2
  real, dimension(mz) :: eta_mz, chi_mz, lnTT
!
!  Magnetic profile - the magnetic field is propto sqrt density
!  Diffusivities propto sqrt(T)
!
  real, parameter :: amplaa_cgs=1e-21  ! 1 nano Gauss (G g^{1/2} cm^-{3/2})
  real, parameter :: chi_cgs = 8.425e23  ! 1 cm^2/s/K^0.5
  real, parameter :: eta_cgs = 8.425e18  ! 1 cm^2/s(/K^0.5 if cspeed)
  real :: amplaa=impossible, chi_th=impossible, eta_cspeed=impossible
  real :: eta=impossible, ybias_aa = 1.0 ! adjust angle of field defalt 45deg
!
!  multiple options for resistivity are possible in magnetic.f90, but here
!  only constant or sound speed dependent eta exclusively
!
  character (len=labellen) :: iresistivity=''
  character (len=labellen), dimension(ninit) :: initaa='nothing'
!
!  start parameters
!
  namelist /initial_condition_pars/ &
      T0hs, const_T_frac, amplaa, zT, cooling_select, &
      heating_select, chi_th, eta_cspeed, eta, iresistivity, initaa, &
      ybias_aa
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
!  Initializes all the f arrays in one call.  This subroutine is called last.
!
!  21-dec-10/ccyang: coded
!  15-feb-15/MR: optional parameter 'profiles' added
!
      use Messages, only: fatal_error
      use EquationOfState, only: getmu, eoscalc
!
      real, dimension (mx,my,mz,mfarray), optional, intent(inout):: f
      real, dimension (:,:),              optional, intent(out)  :: profiles
! 
!  SAMPLE IMPLEMENTATION
!
      

      call keep_compiler_quiet(f)
      if (present(profiles)) then
        call fatal_error('initial_condition_all', &
          'If profiles are asked for, a real initial condition must be specified')
        call keep_compiler_quiet(profiles)
      endif
!
    endsubroutine initial_condition_all
!*****************************************************************************
    subroutine initial_condition_uu(f)
!
!  Initialize the velocity field.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      f(:,:,:,iuu:iuu+2) = 0.
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
!  07-jul-15/fred: coded
!
      use EquationOfState, only: getmu, get_cp1
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: cp1, muhs, Rgas
      real, dimension(mz) :: Gamma_z, Lambda_z, gravz
      integer :: i 
!
!  Set up physical units.
!
      call select_cooling(cooling_select,lncoolT,lncoolH,coolB)
      unit_Gamma  = unit_velocity**3 / unit_length
      if (unit_system=='cgs') then
        a_S = a_S_cgs/unit_velocity*unit_time
        z_S = z_S_cgs/unit_length
        a_D = a_D_cgs/unit_velocity*unit_time
        z_D = z_D_cgs/unit_length
!        rho_fraction = nfraction * m_u_cgs/unit_density
!        hscale = hscale_cgs/unit_length       
!        if (T0hs == impossible) T0hs=T0hs_cgs/unit_temperature
!        if (zT == impossible) zT=zT_cgs/unit_length
        if (GammaUV == impossible) GammaUV=GammaUV_cgs/unit_Gamma
        if (amplaa == impossible) amplaa = &
            amplaa_cgs/unit_magnetic/sqrt(unit_density)
        if (eta == impossible) eta=eta_cgs/unit_velocity/unit_length
        if (eta_cspeed == impossible) eta_cspeed=&
            eta_cgs/unit_velocity/unit_length*sqrt(unit_temperature)
        if (chi_th == impossible) chi_th= &
            chi_cgs/unit_velocity/unit_length*sqrt(unit_temperature)
      else if (unit_system=='SI') then
        call fatal_error('initial_condition_lnrho','SI unit conversions not inplemented')
      endif
!
!  Stellar and dark halo gravity profile 
!
      gravz = -a_S*z/sqrt(z_S**2 + z**2) - a_D * z/z_D
!
 !     call temperature_ansatz(TT,dTdz,d2Tdz2,rho_obs,drhodz,d2rhodz2)
      call temperature_ansatz()
!
!  Cooling function Lambda_i T^beta_i for T_i<T<T_(i+1) 
!
      Lambda_z=0.0
      cool_loop: do i=1,ncool
        if (lncoolT(i) >= lncoolT(i+1)) exit cool_loop
        where (lncoolT(i) <= lnTT .and. lnTT < lncoolT(i+1))
          Lambda_z=Lambda_z+exp(lncoolH(i)+lnTT*coolB(i))
        endwhere
      enddo cool_loop
!
!  UV heating requires heating_select = 'wolfire' in interstellar_run_pars
!  and lthermal_hs = F in interstellar_run_pars
!

      if (heating_select /= 'wolfire') then
        call fatal_error('initial_condition_lnrho','heating_select "'&
           // trim(heating_select) // '" not inplemented')
      else
        Gamma_z = GammaUV*0.5*(1.0+tanh(cUV*(T0UV-exp(lnTT))))
      endif
      call getmu(f,muhs)
      call get_cp1(cp1)
      Rgas = k_B/m_u
      if (.not. lmagnetic) then
        amplaa = 0.
        eta_mz = 0.
      else
        if (iresistivity=='eta-const') then
          if (lroot) print*, 'resistivity: constant eta'
          eta_mz = eta
        else if (iresistivity=='eta-cspeed') then
          if (lroot) print*, 'resistivity: sound speed dependent SN driven ISM'
          eta_mz = eta_cspeed*sqrt(TT)
        else
          if (lroot) print*, 'No such value for iresistivity: ', &
              trim(iresistivity)
          call fatal_error('initialize_magnetic','')
        endif
      endif
      chi_mz=chi_th*sqrt(TT)
!
!  rho is an even function of z, so multiply odd components by sign(1.0,z)
!
      rho = Gamma_z/Lambda_z + chi_mz/cp1/Lambda_z * (sign(1.0,z)* &
            d2Tdz2 + (0.5/TT - mu0*Rgas/(mu0*Rgas*TT + amplaa**2*muhs) &
            ) * dTdz**2 + dTdz*gravz*mu0*muhs/(mu0*Rgas*TT + amplaa**2*muhs) &
            ) + eta_mz*mu0*amplaa**2/(mu0*Rgas*TT + amplaa**2*muhs)**2 * &
            (muhs**2*gravz**2 - 2.*muhs*gravz*Rgas*dTdz + Rgas**2*dTdz**2)/ &
            (4.*Lambda_z)
      do n=n1,n2
        f(:,:,n,ilnrho)=log(rho(n))
        f(:,:,n,inetcool) = rho(n)*Lambda_z(n) - Gamma_z(n)
      enddo
!
    endsubroutine initial_condition_lnrho
!***********************************************************************
    subroutine temperature_ansatz()
!
!  Initialize vertical profile for temperature using empirical density
!  models from which to derive initial density and entropy distributuins
!
!  07-jul-15/fred: coded
!
!  Set up physical units.
!
      if (unit_system=='cgs') then
        rho_fraction = nfraction_cgs * m_u_cgs/unit_density
        hscale = hscale_cgs/unit_length       
        if (T0hs == impossible) T0hs=T0hs_cgs/unit_temperature
        if (zT == impossible) zT=zT_cgs/unit_length
      else if (unit_system=='SI') then
        call fatal_error('initial_condition_lnrho','SI unit conversions not inplemented')
      endif
!
!  Observed density profile and 1st/2nd derivatives used to specify initial
!  temperature profile and calculate initial density for thermo-hydrostatic
!  balance, subject to radiative transfer
!
      rho_obs = rho_fraction(1) * exp(-z**2/hscale(1)**2) &
              + rho_fraction(2) * exp(-z**2/hscale(2)**2) &
              + rho_fraction(3) * exp(-abs(z)/hscale(3)) &
              + rho_fraction(4) * exp(-abs(z)/hscale(4)) &
              + rho_fraction(5) * exp(-abs(z)/hscale(5)) 
      rho_max = sum(rho_fraction)
!
!  rho_obs is an even function, but its derivatives must be odd functions of z
!
      drhodz = -2*z/hscale(1)**2 * rho_fraction(1) * exp(-z**2/hscale(1)**2) &
               -2*z/hscale(2)**2 * rho_fraction(2) * exp(-z**2/hscale(2)**2) &
              -sign(1.0,z)*rho_fraction(3) * exp(-abs(z)/hscale(3))/hscale(3) &
              -sign(1.0,z)*rho_fraction(4) * exp(-abs(z)/hscale(4))/hscale(4) &
              -sign(1.0,z)*rho_fraction(5) * exp(-abs(z)/hscale(5))/hscale(5) 
!
      d2rhodz2 = &
        -sign(1.0,z)*2/hscale(1)**2 * rho_fraction(1) * exp(-z**2/hscale(1)**2) &
        +4*sign(1.0,z)*z**2/hscale(1)**4*rho_fraction(1)*exp(-z**2/hscale(1)**2) &
        -sign(1.0,z)*2/hscale(2)**2 * rho_fraction(2) * exp(-z**2/hscale(2)**2) &
        +4*sign(1.0,z)*z**2/hscale(2)**4*rho_fraction(2)*exp(-z**2/hscale(2)**2) &
        +sign(1.0,z)* rho_fraction(3) * exp(-abs(z)/hscale(3))/hscale(3)**2 &
        +sign(1.0,z)* rho_fraction(4) * exp(-abs(z)/hscale(4))/hscale(4)**2 &
        +sign(1.0,z)* rho_fraction(5) * exp(-abs(z)/hscale(5))/hscale(5)**2 
!
!  Temperature Ansatz, constructed to preserve positive density for all z
!  and column density same order as observation 
!  1st/2nd derivatives required for calulation of initial density for
!  thermo-hydrostatic balance, subject to radiative transfer
!
      TT = T0hs * rho_max/rho_obs * exp(-abs(z)/zT) + const_T_frac * T0hs
      lnTT = log(TT)
!
!  TT is an even function, but its derivatives must be odd functions of z
!
      dTdz =  -TT * (sign(1.0,z)/zT + drhodz/rho_obs)
      d2Tdz2 = TT * (sign(1.0,z)/zT**2 + 2*drhodz/(zT*rho_obs) &
                   + 2*sign(1.0,z)*drhodz**2/rho_obs**2 - d2rhodz2/rho_obs)
!
    endsubroutine temperature_ansatz
!***********************************************************************
    subroutine initial_condition_ss(f)
!
!  Initialize entropy.
!
!  07-may-09/wlad: coded
!
      use EquationOfState, only: eoscalc, ilnrho_lnTT 
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: lnrho, ss
!
      do n=n1,n2
        lnrho=log(rho(n))
!
        call eoscalc(ilnrho_lnTT,lnrho,lnTT(n),ss=ss)
!
        f(:,:,n,iss)=ss
!
      enddo
!
    endsubroutine initial_condition_ss
!***********************************************************************
    subroutine initial_condition_aa(f)
!
!  Initialize the magnetic vector potential.
!
!  This routine sets up an initial magnetic field y-parallel(azimuthal) 
!  with a magnitude directly proportional to the density. 
!
      use Mpicomm, only: mpireduce_sum, mpibcast_real
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      integer :: i ,j, l, m, n
!
      if (.not. lmagnetic) then
!
        call keep_compiler_quiet(f)
      else     
        do j=1,ninit
!  
          select case (initaa(j))
          case ('nothing'); if (lroot .and. j==1) print*,'init_aa: nothing'
          case ('zero', '0'); f(:,:,:,iax:iaz) = 0.0
          case ('uniform-x')
            f(:,:,:,iax:iay) = 0.0
            do l=l1,l2
            do n=n1,n2
              f(l,m1:m2,n,iaz) = amplaa * y(m1:m2) * sqrt(rho(n))
            enddo
            enddo
          case ('uniform-y')
            f(:,:,:,iax:iay) = 0.0
            do m=m1,m2
            do n=n1,n2
              f(l1:l2,m,n,iaz) = amplaa * x(l1:l2) * sqrt(rho(n))
            enddo
            enddo
          case ('uniform-x+y')
            f(:,:,:,iax:iay) = 0.0
            do l=l1,l2
            do m=m1,m2
            do n=n1,n2
              f(l,m,n,iaz) = amplaa * (x(l) + ybias_aa*y(m)) * sqrt(rho(n)) / &
                             sqrt(1 + ybias_aa**2)
            enddo
            enddo
            enddo
          case default
!  
!  Catch unknown values.
!  
            call fatal_error('inititial_condition_aa', &
                'init_aa value "' // trim(initaa(j)) // '" not recognised')
!  
          endselect
!  
!  End loop over initial conditions.
!  
        enddo
!
      endif
!
    endsubroutine initial_condition_aa
!***********************************************************************
    subroutine initial_condition_aatest(f)
!
!  Initialize testfield.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_aatest
!***********************************************************************
    subroutine initial_condition_uutest(f)
!
!  Initialize testflow.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uutest
!***********************************************************************
    subroutine initial_condition_lncc(f)
!
!  Initialize passive scalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lncc
!***********************************************************************
    subroutine initial_condition_chiral(f)
!
!  Initialize chiral.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chiral
!***********************************************************************
    subroutine initial_condition_chemistry(f)
!
!  Initialize chemistry.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_chemistry
!***********************************************************************
    subroutine initial_condition_uud(f)
!
!  Initialize dust fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uud
!***********************************************************************
    subroutine initial_condition_nd(f)
!
!  Initialize dust fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_nd
!***********************************************************************
    subroutine initial_condition_uun(f)
!
!  Initialize neutral fluid velocity.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_uun
!***********************************************************************
    subroutine initial_condition_lnrhon(f)
!
!  Initialize neutral fluid density.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_lnrhon
!***********************************************************************
    subroutine initial_condition_ecr(f)
!
!  Initialize cosmic rays.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_ecr
!***********************************************************************
    subroutine initial_condition_fcr(f)
!
!  Initialize cosmic ray flux.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_fcr
!***********************************************************************
    subroutine initial_condition_solid_cells(f)
!
!  Initialize solid cells.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_solid_cells
!***********************************************************************
    subroutine initial_condition_cctest(f)
!
!  Initialize testscalar.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine initial_condition_cctest
!***********************************************************************
    subroutine initial_condition_xxp(f,fp)
!
!  Initialize particles' positions.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_xxp
!*****************************************************************************
    subroutine select_cooling(cooling_select,lncoolT,lncoolH,coolB)
!
!  Routine for selecting parameters for temperature dependent cooling
!  Lambda. 
!
      character (len=labellen), intent(IN) :: cooling_select  
      real, dimension (:), intent(OUT)  :: lncoolT, coolB
      double precision, dimension (:), intent(OUT)  :: lncoolH
!
!
      if (cooling_select == 'RBNr') then
        if (lroot) print*,'initialize_interstellar: RBN cooling fct (revised)'
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
      else if (cooling_select == 'WSW') then
        if (lroot) print*,'initialize_interstellar: WSW cooling fct'
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
        if (lroot) print*,'initialize_interstellar: no cooling applied'
        coolT_cgs=tiny(0.0)
        coolH_cgs=tiny(0.0)
        coolB=tiny(0.)
      endif
      unit_Lambda = unit_velocity**2 / unit_density / unit_time
      lncoolH(1:ncool) = real(log(coolH_cgs(1:ncool)) - log(unit_Lambda) &
                              + log(unit_temperature**coolB(1:ncool)) &
                              + log(coolingfunction_scalefactor))
      lncoolT(1:ncool+1) = real(log(coolT_cgs(1:ncool+1) / unit_temperature))
!
    endsubroutine select_cooling
!***********************************************************************
    subroutine initial_condition_vvp(f,fp)
!
!  Initialize particles' velocities.
!
!  07-may-09/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (:,:), intent(inout) :: fp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
!
    endsubroutine initial_condition_vvp
!***********************************************************************
    subroutine read_initial_condition_pars(iostat)
!
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include '../parallel_unit.h'
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
    subroutine initial_condition_clean_up
!
!  04-may-11/dhruba: coded
! dummy
!
    endsubroutine initial_condition_clean_up
!***********************************************************************
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
