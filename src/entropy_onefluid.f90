! $Id$
!
!  This module takes care of entropy (initial condition
!  and time advance) for a fluid consisting of gas and perfectly
!  coupled pressureless dust.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .true.
! CPARAM logical, parameter :: ltemperature = .false.
! CPARAM logical, parameter :: lthermal_energy = .false.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cp1tilde; glnTT(3)
! PENCILS PROVIDED TT; TT1; Ma2; ugss; hss(3,3); hlnTT(3,3)
! PENCILS PROVIDED del2ss; del6ss; fpres(3); sglnTT(3); advec_cs2
!
!***************************************************************
module Energy
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Interstellar
  use Messages
!
  implicit none
!
  include 'energy.h'
!
  real :: ss_const=0.0
  real :: T0=1.0
  logical :: lupw_ss=.false.
  logical, pointer :: lpressuregradient_gas
  logical :: lviscosity_heat=.true.
  logical :: ladvection_energy=.true.
  logical :: lenergy_slope_limited=.false.
  character (len=labellen), dimension(ninit) :: initss='nothing'
  character (len=intlen) :: iinit_str
  real, dimension(3) :: beta_glnrho_global=0.
!
  namelist /entropy_init_pars/ &
      initss, ss_const, T0, beta_glnrho_global, ladvection_energy
!
  namelist /entropy_run_pars/ &
      lupw_ss, beta_glnrho_global, ladvection_energy
!
  integer :: idiag_dtc=0,idiag_ethm=0,idiag_ethdivum=0,idiag_ssm=0
  integer :: idiag_ugradpm=0,idiag_ethtot=0,idiag_ssmphi=0
  integer :: idiag_yHm=0,idiag_yHmax=0,idiag_TTm=0,idiag_TTmax=0,idiag_TTmin=0
  integer :: idiag_fconvz=0,idiag_dcoolz=0,idiag_fradz=0,idiag_fturbz=0
  integer :: idiag_ssmz=0,idiag_ssmy=0,idiag_ssmx=0,idiag_TTmz=0
!
  real, dimension(:), pointer :: beta_glnrho_scaled
!
  contains
!***********************************************************************
    subroutine register_energy
!
!  Initialise variables which should know that we solve an entropy
!  equation: iss, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
      call farray_register_pde('ss',iss)
!
!  Identify version number.
!
      if (lroot) call svn_id( &
          "$Id$")
!
      call put_shared_variable('lviscosity_heat',lviscosity_heat,caller='register_energy')
!
      if (.not.ldensity.or.lboussinesq) &
        call put_shared_variable('beta_glnrho_global',beta_glnrho_global)
!
    endsubroutine register_energy
!***********************************************************************
    subroutine initialize_energy(f)
!
!  Called by run.f90 after reading parameters, but before the time loop.
!
!  21-jul-2002/wolf: coded
!
      use EquationOfState
      use Gravity, only: gravz,g0
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: i
!
!  check any module dependencies.
!
      if (.not. leos) call fatal_error('initialize_energy', &
          'EOS=noeos but entropy requires an EQUATION OF STATE for the fluid')
      call select_eos_variable('ss',iss)
!
!  For global density gradient beta=H/r*dlnrho/dlnr, calculate actual
!  gradient dlnrho/dr = beta/H.
!
      call get_shared_variable('beta_glnrho_scaled',beta_glnrho_scaled)
!
!  Turn off pressure gradient term and advection for 0-D runs.
!
      if (nxgrid*nygrid*nzgrid==1) then
        lpressuregradient_gas=.false.
        ladvection_energy=.false.
        if (lroot) then
          print*, 'initialize_energy: 0-D run, turned off pressure gradient term'
          print*, 'initialize_energy: 0-D run, turned off advection of entropy'
        endif
      endif
!
!  logical variable lpressuregradient_gas shared with hydro modules
!
      call get_shared_variable('lpressuregradient_gas',lpressuregradient_gas,caller='initialize_energy')
!
      if (lenergy_slope_limited) call not_implemented('initialize_energy','slope-limited diffusion')
!
      call keep_compiler_quiet(f)
!
      endsubroutine initialize_energy
!***********************************************************************
    subroutine read_energy_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_init_pars, IOSTAT=iostat)
!
    endsubroutine read_energy_init_pars
!***********************************************************************
    subroutine write_energy_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_init_pars)
!
    endsubroutine write_energy_init_pars
!***********************************************************************
    subroutine read_energy_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=entropy_run_pars, IOSTAT=iostat)
!
    endsubroutine read_energy_run_pars
!***********************************************************************
    subroutine write_energy_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=entropy_run_pars)
!
    endsubroutine write_energy_run_pars
!***********************************************************************
    subroutine init_energy(f)
!
!  Initialise entropy; called from start.f90.
!
      use Sub
      use Gravity
      use General, only: itoa
      use Initcond
      use InitialCondition, only: initial_condition_ss
      use EquationOfState,  only: rho0, lnrho0, isothermal_entropy, &
                                  isothermal_lnrho_ss, eoscalc, ilnrho_pp
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: j
!
      intent(inout) :: f
!
      do j=1,ninit
!
        if (initss(j)/='nothing') then
!
          iinit_str=itoa(j)
!
!  select different initial conditions
!
          select case (initss(j))
!
            case ('zero', '0'); f(:,:,:,iss) = 0.
            case ('const_ss'); f(:,:,:,iss)=f(:,:,:,iss)+ss_const
            case ('isothermal'); call isothermal_entropy(f(:,:,:,ilnrho),T0,f(:,:,:,ss))
            case ('isothermal_lnrho_ss')
              if (lroot) print*, 'init_energy: Isothermal density and entropy stratification'
              call isothermal_lnrho_ss(f(:,:,:,ilnrho),T0,rho0,f(:,:,:,iss))
            case default
!
!  Catch unknown values
!
              call fatal_error('init_energy','no such initss('//trim(iinit_str)//'): '//trim(initss(j)))
          endselect
!
        endif
!
      enddo
!
!  Interface for user's own initial condition
!
      if (linitial_condition) call initial_condition_ss(f)
!
    endsubroutine init_energy
!***********************************************************************
    subroutine pencil_criteria_energy
!
!  All pencils that the Entropy module depends on are specified here.
!
!  20-11-04/anders: coded
!
      if (ldt) lpenc_requested(i_cs2)=.true.
      if (lpressuregradient_gas) then
        lpenc_requested(i_cs2)=.true.
        lpenc_requested(i_cp1tilde)=.true.
        lpenc_requested(i_glnrho)=.true.
        lpenc_requested(i_gss)=.true.
      endif
      if (ladvection_energy) lpenc_requested(i_ugss)=.true.
      if (pretend_lnTT) lpenc_requested(i_divu)=.true.
      if (lpressuregradient_gas) lpenc_requested(i_cp1tilde)=.true.
!
      if (any(beta_glnrho_scaled/=0.0)) lpenc_requested(i_cs2)=.true.
!
      lpenc_diagnos2d(i_ss)=.true.
!
      if (idiag_ethdivum/=0) lpenc_diagnos(i_divu)=.true.
      if (idiag_ssm/=0 .or. idiag_ssmz/=0 .or. idiag_ssmy/=0.or.idiag_ssmx/=0) &
          lpenc_diagnos(i_ss)=.true.
      if (idiag_ethm/=0 .or. idiag_ethtot/=0 .or. idiag_ethdivum/=0) then
          lpenc_diagnos(i_rho)=.true.
          lpenc_diagnos(i_ee)=.true.
      endif
      if (idiag_fconvz/=0 .or. idiag_fturbz/=0 ) then
          lpenc_diagnos(i_rho)=.true.
          lpenc_diagnos(i_TT)=.true.  !(to be replaced by enthalpy)
      endif
      if (idiag_TTm/=0 .or. idiag_TTmz/=0 .or. idiag_TTmax/=0 .or. idiag_TTmin/=0) &
          lpenc_diagnos(i_TT)=.true.
      if (idiag_yHm/=0 .or. idiag_yHmax/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_dtc/=0) lpenc_diagnos(i_cs2)=.true.
!
    endsubroutine pencil_criteria_energy
!***********************************************************************
    subroutine pencil_interdep_energy(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_del2lnTT)) then
        lpencil_in(i_del2lnrho)=.true.
        lpencil_in(i_del2ss)=.true.
      endif
      if (lpencil_in(i_glnTT)) then
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_gss)=.true.
      endif
      if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
      if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
      if (lpencil_in(i_ugss)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_gss)=.true.
      endif
      if (pretend_lnTT .and. lpencil_in(i_glnTT)) lpencil_in(i_gss)=.true.
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif
      if (lpencil_in(i_hlnTT)) then
        if (pretend_lnTT) then
          lpencil_in(i_hss)=.true.
        else
          lpencil_in(i_hlnrho)=.true.
          lpencil_in(i_hss)=.true.
        endif
      endif
!  The pencils cs2 and cp1tilde come in a bundle, so enough to request one.
      if (lpencil_in(i_cs2) .and. lpencil_in(i_cp1tilde)) &
          lpencil_in(i_cp1tilde)=.false.
!
    endsubroutine pencil_interdep_energy
!***********************************************************************
    subroutine calc_pencils_energy(f,p)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
!  16-05-12/MR: dead branch for pressure force calculation added
!
      use EquationOfState
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension (nx,3) :: glnrhod, glnrhotot
      real, dimension (nx) :: eps
      integer :: i
!
      intent(in) :: f
      intent(inout) :: p
! ss
      if (lpencil(i_ss)) p%ss=f(l1:l2,m,n,iss)
! gss
      if (lpencil(i_gss)) call grad(f,iss,p%gss)
! pp
      if (lpencil(i_pp)) call eoscalc(f,nx,pp=p%pp)
! ee
      if (lpencil(i_ee)) call eoscalc(f,nx,ee=p%ee)
! lnTT
      if (lpencil(i_lnTT)) call eoscalc(f,nx,lnTT=p%lnTT)
! TT
      if (lpencil(i_TT)) p%TT=exp(p%lnTT)
! TT1
      if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
! cs2 and cp1tilde
      if (lpencil(i_cs2) .or. lpencil(i_cp1tilde)) call pressure_gradient(f,p%cs2,p%cp1tilde)
! Ma2
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
! glnTT
      if (lpencil(i_glnTT)) then
        if (pretend_lnTT) then
           p%glnTT=p%gss
        else
          call temperature_gradient(f,p%glnrho,p%gss,p%glnTT)
        endif
      endif
! ugss
      if (lpencil(i_ugss)) call u_dot_grad(f,iss,p%gss,p%uu,p%ugss,UPWIND=lupw_ss)
!ajwm Should probably combine the following two somehow.
! hss
      if (lpencil(i_hss)) call g2ij(f,iss,p%hss)
! del2ss
      if (lpencil(i_del2ss)) call del2(f,iss,p%del2ss)
! del2lnTT
      if (lpencil(i_del2lnTT)) call temperature_laplacian(f,p)
! del6ss
      if (lpencil(i_del6ss)) call del6(f,iss,p%del6ss)
! hlnTT
      if (lpencil(i_hlnTT)) then
         if (pretend_lnTT) then
           p%hlnTT=p%hss
         else
           call temperature_hessian(f,p%hlnrho,p%hss,p%hlnTT)
         endif
       endif
! uud (for dust continuity equation in one fluid approximation)
      if (lpencil(i_uud)) p%uud(:,:,1)=p%uu
! divud (for dust continuity equation in one fluid approximation)
      if (lpencil(i_divud)) p%divud(:,1)=p%divu
!
      if (lpencil(i_fpres)) &
        call not_implemented('calc_pencils_energy', &
                 'calculation of pressure force for entropy_onefluid')
! sglnTT 
      if (lpencil(i_sglnTT)) &
        call not_implemented('calc_pencils_energy','pencil sglnTT for entropy_onefluid')
!
!  ``cs2/dx^2'' for timestep
!
      if (lhydro.and.ldensity.and.lfirst.and.ldt) then
        p%advec_cs2=p%cs2*dxyz_2
        if (headtt.or.ldebug) print*,'calc_pencils_energy: max(advec_cs2) =',maxval(p%advec_cs2)
      endif
!
    endsubroutine calc_pencils_energy
!***********************************************************************
    subroutine denergy_dt(f,df,p)
!
!  Calculate right hand side of entropy equation.
!
!      use Conductivity, only: heat_conductivity
      use Special, only: special_calc_energy
      use Sub
      use Viscosity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: Hmax
      integer :: j,ju
!
      intent(inout)  :: f,p
      intent(out) :: df
!
      Hmax = 0.0
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'denergy_dt: SOLVE denergy_dt'
      if (headtt) call identify_bcs('ss',iss)
!
      if (lhydro) then
!
!  Pressure term in momentum equation (setting lpressuregradient_gas to
!  .false. allows suppressing pressure term for test purposes).
!
        if (lpressuregradient_gas) then
          do j=1,3
            ju=j+iuu-1
            df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) - p%cs2*(p%glnrho(:,j) + p%cp1tilde*p%gss(:,j))
          enddo
!
!  Add pressure force from global density gradient.
!  WARNING (AJ): This may be implemented inconsistently, since we have
!  here linearised rho and P independently.
!
          if (any(beta_glnrho_scaled/=0.0)) then
            if (headtt) print*, 'denergy_dt: adding global pressure gradient force'
            do j=1,3
              df(l1:l2,m,n,(iux-1)+j) = df(l1:l2,m,n,(iux-1)+j) - p%cs2*beta_glnrho_scaled(j)
            enddo
          endif
        endif
!
!  Advection term.
!
        if (ladvection_energy) df(l1:l2,m,n,iss) = df(l1:l2,m,n,iss) - p%ugss
      endif
!
!  Calculate viscous contribution to entropy.
!
      if (lviscosity) call calc_viscous_heat(df,p,Hmax)
!
!  Thermal conduction.
!
!      if (lconductivity) call heat_conductivity(f,df,p)
!
!  Entry possibility for "personal" entries.
!
      if (lspecial) call special_calc_energy(f,df,p)
!
      if (lhydro.and.ldensity.and.lfirst.and.ldt) advec_cs2=p%advec_cs2
!
      call calc_diagnostics_energy(f,p)

    endsubroutine denergy_dt
!***********************************************************************
    subroutine calc_diagnostics_energy(f,p)

      use Diagnostics

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p
!
!  Calculate entropy related diagnostics.
!
      call keep_compiler_quiet(f)

      if (ldiagnos) then
        call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_dtc/=0) call max_mn_name(sqrt(p%advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        if (idiag_ethm/=0) call sum_mn_name(p%rho*p%ee,idiag_ethm)
        if (idiag_ethtot/=0) call integrate_mn_name(p%rho*p%ee,idiag_ethtot)
        if (idiag_ethdivum/=0) call sum_mn_name(p%rho*p%ee*p%divu,idiag_ethdivum)
        call sum_mn_name(p%ss,idiag_ssm)
      endif
!
!  1D averages. Happens at every it1d timesteps, NOT at every it1
!  idiag_fradz is done in the calc_heatcond routine
!
      if (l1davgfirst) then
        if (idiag_fconvz/=0) call xysum_mn_name_z(p%rho*p%uu(:,3)*p%TT,idiag_fconvz)
        call xysum_mn_name_z(p%ss,idiag_ssmz)
        call xzsum_mn_name_y(p%ss,idiag_ssmy)
        call yzsum_mn_name_x(p%ss,idiag_ssmx)
        call xysum_mn_name_z(p%TT,idiag_TTmz)
      endif
!
    endsubroutine calc_diagnostics_energy
!***********************************************************************
    subroutine rprint_energy(lreset,lwrite)
!
!  Reads and registers print parameters relevant to entropy.
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
!
      integer :: iname,inamez,inamey,inamex,irz
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset.
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_dtc=0; idiag_ethm=0; idiag_ethdivum=0; idiag_ssm=0
        idiag_ugradpm=0; idiag_ethtot=0; idiag_ssmphi=0
        idiag_yHmax=0; idiag_yHm=0; idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_fconvz=0; idiag_dcoolz=0; idiag_fradz=0; idiag_fturbz=0
        idiag_ssmz=0; idiag_ssmy=0; idiag_ssmx=0; idiag_TTmz=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'ethtot',idiag_ethtot)
        call parse_name(iname,cname(iname),cform(iname),'ethdivum',idiag_ethdivum)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',idiag_ugradpm)
        call parse_name(iname,cname(iname),cform(iname),'yHm',idiag_yHm)
        call parse_name(iname,cname(iname),cform(iname),'yHmax',idiag_yHmax)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
      enddo
!
!  Check for those quantities for which we want xy-averages.
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fturbz',idiag_fturbz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fconvz',idiag_fconvz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'dcoolz',idiag_dcoolz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'fradz',idiag_fradz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'ssmz',idiag_ssmz)
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'TTmz',idiag_TTmz)
      enddo
!
!  Check for those quantities for which we want xz-averages.
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'ssmy',idiag_ssmy)
      enddo
!
!  Check for those quantities for which we want yz-averages.
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'ssmx',idiag_ssmx)
      enddo
!
!  Check for those quantities for which we want phi-averages.
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ssmphi',idiag_ssmphi)
      enddo
!
!  Write column where which entropy variable is stored.
!
      if (lwr) then
        call farray_index_append('iyH',iyH)
        call farray_index_append('ilnTT',ilnTT)
        call farray_index_append('iTT',iTT)
      endif
!
    endsubroutine rprint_energy
!***********************************************************************    
    subroutine update_char_vel_energy(f)
!
! TB implemented.
!
!   25-sep-15/MR+joern: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f

      call keep_compiler_quiet(f)

      call warning('update_char_vel_energy', &
           'characteristic velocity not yet implemented for entropy_onefluid')

    endsubroutine update_char_vel_energy
!***********************************************************************
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'energy_dummies.inc'
!***********************************************************************
endmodule Energy
