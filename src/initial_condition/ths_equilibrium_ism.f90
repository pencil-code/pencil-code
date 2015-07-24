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
  real :: a_S, a_D
  real, parameter :: a_S_cgs=4.4e-9, a_D_cgs=1.7e-9
  real :: z_S, z_D
  real, parameter :: z_S_cgs=6.172e20, z_D_cgs=3.086e21
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
  double precision, dimension(11) :: lncoolH, coolH_cgs
  real, dimension(11) :: coolT_cgs
  real, dimension(11) :: coolB, lncoolT
  integer :: ncool
!
!  Heating function, cooling function and mass movement
!  method selection.
!
  character (len=labellen) :: cooling_select  = 'RBN'
  character (len=labellen) :: heating_select  = 'wolfire'
!
!
  real, parameter :: rho0ts_cgs=3.5E-24, T0hs_cgs=7.088E2
  real :: rho0ts=impossible, T0hs=impossible
!
!
!  TT & z-dependent uv-heating profile
!
  real, dimension(mz) :: TT, zrho
  logical :: lthermal_hse=.true.
!
!  Magnetic profile
!
  real, parameter :: amplaa_cgs=1e-9 ! 1 nano Gauss
  real :: amplaa=impossible
!
!  start parameters
!
  namelist /initial_condition_pars/ &
      rho0ts, T0hs, lthermal_hse, amplaa
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
      use EquationOfState, only: getmu
      use Interstellar, only: select_cooling 
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: logrho, muhs
!
!  Set up physical units.
!
      call select_cooling(cooling_select,lncoolT,lncoolH,coolB)
      if (unit_system=='cgs') then
        a_S = a_S_cgs/unit_velocity*unit_time
        z_S = z_S_cgs/unit_length
        a_D = a_D_cgs/unit_velocity*unit_time
        z_D = z_D_cgs/unit_length
        if (rho0ts == impossible) rho0ts=rho0ts_cgs/unit_density
        if (T0hs == impossible) T0hs=T0hs_cgs/unit_temperature
      else if (unit_system=='SI') then
        call fatal_error('initial_condition_lnrho','SI unit conversions not inplemented')
      endif
!
      call getmu(f,muhs)
!
      cool=0.0
      cool_loop: do i=1,ncool
        if (lncoolT(i) >= lncoolT(i+1)) exit cool_loop
        where (lncoolT(i) <= lnTT .and. lnTT < lncoolT(i+1))
          cool=cool+exp(lncoolH(i)+lnrho+lnTT*coolB(i))
        endwhere
      enddo cool_loop
      if (lroot) print*, 'thermal-hs: '// &
          'hydrostatic thermal equilibrium density and entropy profiles'
!
!  requires heating_select = 'wolfire' in interstellar_run_pars
!
      Gamma_z = GammaUV*0.5*(1.0+tanh(cUV*(T0UV-exp(lnTT))))
      do n=1,mz
        if (lthermal_hse) then
          logrho = log(rho0ts)+(a_S*z_S*m_u*muhs/k_B/T0hs)*(log(T0hs)- &
              log(T0hs/(a_S*z_S)* &
              (a_S*sqrt(z_S**2+(z(n))**2)+0.5*a_D*(z(n))**2/z_D)))
        else
          logrho = log(rho0ts)-0.015*(- &
              a_S*z_S+ &
              a_S*sqrt(z_S**2+(z(n))**2)+0.5*a_D*(z(n))**2/z_D)
        endif
        if (logrho < -40.0) logrho=-40.0
        zrho(n)=exp(logrho)
      enddo
      do n=n1,n2
        f(:,:,n,ilnrho)=log(zrho(n))
      enddo
      if (lroot) print*, 'zhro =',minval(zrho),maxval(zrho),'T0hs =',T0hs,'rho0ts =',rho0ts
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
      use EquationOfState, only: eoscalc, ilnrho_lnTT 
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real :: TT, lnTT, lnrho, ss
!
!  Set up physical units.
!
      if (unit_system=='cgs') then
        a_S = a_S_cgs/unit_velocity*unit_time
        z_S = z_S_cgs/unit_length
        a_D = a_D_cgs/unit_velocity*unit_time
        z_D = z_D_cgs/unit_length
        if (T0hs == impossible) T0hs=T0hs_cgs/unit_temperature
      else if (unit_system=='SI') then
        call fatal_error('initial_condition_lnrho','SI unit conversions not inplemented')
      endif
!
      if (lroot) print*, 'thermal-hs: '// &
          'hydrostatic thermal equilibrium density and entropy profiles'
      do n=1,mz
        if (ldensity_nolog) then
          lnrho=log(f(1,1,n,irho))
        else
          lnrho=f(1,1,n,ilnrho)
        endif
!
        TT=T0hs/(a_S*z_S)* &
            (a_S*sqrt(z_S**2+(z(n))**2)+0.5*a_D*(z(n))**2/z_D)
        lnTT=log(TT)
!
        call eoscalc(ilnrho_lnTT,lnrho,lnTT,ss=ss)
!
        f(:,:,n,iss)=ss
!
      enddo
!
!
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
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
      integer :: i ,icpu
      real :: tmp1,tmp3
      real, dimension(ncpus)::sumtmp,tmp2
!
      if (unit_system=='cgs') then
        if (amplaa == impossible) amplaa = amplaa_cgs/unit_magnetic
      else if (unit_system=='SI') then
        call fatal_error('initial_condition_magnetic','SI unit conversions not inplemented')
      endif
!
      tmp2(:)=0.0
      sumtmp(:)=0.0
      if (ldensity_nolog) then
        tmp1=sum(f(l1:l2,m1,n1:n2,irho))
      else
        tmp1=sum(exp(f(l1:l2,m1,n1:n2,ilnrho)))
      endif
!
!  Calculate the total mass on each processor tmp1 and identify it with the
!  appropriate processor in the array tmp2
!
      do icpu=1,ncpus
        tmp3=tmp1
        call mpibcast_real(tmp3,icpu-1)
        tmp2(icpu)=tmp3
      enddo
!
!  If nprocz is 1 then start summing mass below from zero (sumtmp above).
!  Otherwise sum the masses on the processors below from which to start
!  summing the mass on this processor.
!
      if (ncpus>nprocy) then
        do icpu=nprocy+1,ncpus
          sumtmp(icpu)=sumtmp(icpu-nprocy)+tmp2(icpu-nprocy)
        enddo
      endif
      if (lroot) print*,'sumtmp =',sumtmp
      print*,'sumtmp on iproc =',sumtmp(iproc+1),iproc
      if (amplaa==0) then
        f(:,:,:,iaa:iaa+2)=0
        if (lroot) print*,'ferriere_uniform_y: set B field to zero; i=',i
      else
        print*,'ferriere_uniform_y: uniform y-field approx rho; i=',i
        if ((ip<=16).and.lroot) print*,'uniform_y: amplaa=',amplaa
        do n=n1,n2
        do m=m1,m2
          if (ldensity_nolog) then
            f(l1:l2,m,n,iaa)=amplaa*(sumtmp(iproc+1)+&
                sum(f(l1:l2,m,n1:n,irho)))*dx*dz
          else
            f(l1:l2,m,n,iaa)=amplaa*(sumtmp(iproc+1)+&
                sum(exp(f(l1:l2,m,n1:n,ilnrho))))*dx*dz
          endif
          f(l1:l2,m,n,iaa+1)=0.0
          f(l1:l2,m,n,iaa+2)=0.0
        enddo
        enddo
      endif
!  SAMPLE IMPLEMENTATION
!
      call keep_compiler_quiet(f)
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
