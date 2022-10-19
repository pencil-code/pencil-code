! $Id$
!
!  This module is for systems with spatially fixed entropy
!  distribution. This implies Ds/Dt=u.grads only, which is used
!  in Ds/Dt=(1/gamma)*Dlnp/Dt-Dlnrho/Dt, for example.
!  This procedure has been used in the context of accretion discs
!  (see von Rekowski et al., 2003, A&A 398, 825). The shock jump
!  relations are modified (see Sect 9.3.6 of Brandenburg 2003, in
!  "Computational aspects...", ed. Ferriz-Mas & Nunez, Taylor & Francis,
!  or astro-ph/0109497.
!
!  This current implementation has temporarily been used to imitate
!  a corona in solar context.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ltemperature = .false.
! CPARAM logical, parameter :: lthermal_energy = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Energy
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  real :: dummy,xss_corona,wss_corona,ss_corona
  real, dimension(nx) :: profxss,dprofxss
  character (len=labellen) :: iss_profile='nothing'
!
  ! parameters necessary for consistency with other modules,
  ! but not directly used.
  real :: hcond0=0.,hcond1=impossible,chi=impossible
  real :: Fbot=impossible,FbotKbot=impossible,Kbot=impossible
  real :: Ftop=impossible,FtopKtop=impossible
  real :: ss_const=1.
  logical :: lmultilayer=.true.
  logical :: lcalc_heatcond_constchi=.false.
  logical :: lenergy_slope_limited=.false.
  character (len=labellen), dimension(ninit) :: initss='nothing'
  integer :: pushpars2c  ! should be procedure pointer (F2003)
!
  ! input parameters
  namelist /entropy_init_pars/ &
      initss, ss_const
!
  ! run parameters
  namelist /entropy_run_pars/ &
      dummy
!
  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_dtc=0,idiag_ethm=0,idiag_ethdivum=0,idiag_ssm=0
  integer :: idiag_ugradpm=0,idiag_ethtot=0,idiag_dtchi=0,idiag_ssmphi=0
!
  contains
!***********************************************************************
    subroutine register_energy
!
!  initialise variables which should know that we solve an energy
!  equation: iss, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      lentropy = .true.
!
!  identify version number
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_energy
!***********************************************************************
    subroutine initialize_energy
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!
      use Cdata
!
      integer :: ierr
!
  ! run parameters
  namelist /entropy_initialize_pars/ &
       iss_profile,xss_corona,wss_corona,ss_corona
!
!  read initialization parameters
!
      if (lroot) print*,'read energy_initialize_pars'
      open(1,FILE='run.in',FORM='formatted',STATUS='old')
      read(1,NML=entropy_initialize_pars,ERR=99,IOSTAT=ierr)
      close(1)
!
!  set profiles
!
      if (iss_profile=='nothing') then
        profxss=1.
      elseif (iss_profile=='corona') then
        profxss=ss_corona*.5*(1.+tanh((x(l1:l2)-xss_corona)/wss_corona))
        dprofxss=ss_corona*.5/(wss_corona*cosh((x(l1:l2)-xss_corona)/wss_corona)**2)
      endif
!
      if (lroot) then
        print*,'initialize_energy: iss_profile=',iss_profile
        print*,'initialize_energy: wss_corona,ss_corona=',xss_corona,wss_corona,ss_corona
      endif
!
      return
99    if (lroot) print*,'i/o error'
!
    endsubroutine initialize_energy
!***********************************************************************
    subroutine init_energy(f)
!
!  initialise energy; called from start.f90
!  07-nov-2001/wolf: coded
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      do j=1,ninit
        select case (initss(j))
!
          case ('zero', '0'); f(:,:,:,iss) = 0.
          case ('const_ss'); f(:,:,:,iss)=f(:,:,:,iss)+ss_const
!
        endselect
      enddo
!
    endsubroutine init_energy
!***********************************************************************
    subroutine energy_after_boundary(f)

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f

      call keep_compiler_quiet(f)

      if (lenergy_slope_limited) &
        call fatal_error('energy_after_boundary', &
                         'Slope-limited diffusion not implemented')

    endsubroutine energy_after_boundary
!***********************************************************************
    subroutine denergy_dt(f,df,uu,glnrho,divu,rho1,lnrho,cs2,TT1,shock,gshock,bb,bij)
!
!  28-mar-02/axel: dummy routine, adapted from entropy.f90 of 6-nov-01.
!  19-may-02/axel: added isothermal pressure gradient
!   9-jun-02/axel: pressure gradient added to du/dt already here
!
      use Density
      use Diagnostics
      use Sub
      use EquationOfState, only: pressure_gradient
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: bij
      real, dimension (nx,3) :: uu,glnrho,gshock,bb
      real, dimension (nx) :: lnrho,rho1,divu,cs2,TT1,uglnrho,shock,rho
      real, dimension (nx) :: ss,cp1tilde
      integer :: j,ju
!
      intent(in) :: f,uu,glnrho,rho1,shock,gshock
      intent(out) :: cs2,TT1  !(df is dummy)
!
!  sound speed squared and inverse temperature
!
      TT1=0.
      ss=profxss
      call pressure_gradient(f,lnrho,cs2,cp1tilde)
      print*, it, m, n, maxval(cs2), maxval(ss)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=cs2*dxyz_2
      if (headtt.or.ldebug) print*,'denergy_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  subtract isothermal/polytropic pressure gradient term in momentum equation
!
      if (lhydro) then
        df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-cs2*(glnrho(:,1)+dprofxss)
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-cs2*(glnrho(:,2))
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-cs2*(glnrho(:,3))
      endif
!
!  Calculate energy related diagnostics
!
      if (ldiagnos) then
        if (idiag_dtc/=0) &
            call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        if (idiag_ugradpm/=0) then
          call dot_mn(uu,glnrho,uglnrho)
          call sum_mn_name(rho*cs2*uglnrho,idiag_ugradpm)
        endif
      endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(uu)
      call keep_compiler_quiet(divu)
      call keep_compiler_quiet(rho1)
      call keep_compiler_quiet(shock)
      call keep_compiler_quiet(gshock)
      call keep_compiler_quiet(bb)
      call keep_compiler_quiet(bij)
!
    endsubroutine denergy_dt
!***********************************************************************
    subroutine calc_diagnostics_energy(f,p)

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)

    endsubroutine calc_diagnostics_energy
!***********************************************************************
    subroutine energy_after_boundary(f)
!
!  dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine energy_after_boundary
!***********************************************************************
    subroutine rprint_energy(lreset,lwrite)
!
!  reads and registers print parameters relevant to energy
!
!   1-jun-02/axel: adapted from magnetic fields
!
      use Cdata
      use Diagnostics
      use FArrayManager, only: farray_index_append
!
      integer :: iname,irz
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
        idiag_dtc=0; idiag_ethm=0; idiag_ethdivum=0; idiag_ssm=0
        idiag_ugradpm=0; idiag_ethtot=0; idiag_dtchi=0; idiag_ssmphi=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'ethtot',idiag_ethtot)
        call parse_name(iname,cname(iname),cform(iname),&
            'ethdivum',idiag_ethdivum)
        call parse_name(iname,cname(iname),cform(iname),'ethm',idiag_ethm)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'ugradpm',idiag_ugradpm)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),'ssmphi',idiag_ssmphi)
      enddo
!
!  write column where which energy variable is stored
!
      if (lwr) then
        call farray_index_append('iss',iss)
        call farray_index_append('iyH',iyH)
        call farray_index_append('ilnTT',ilnTT)
        call farray_index_append('iTT',iTT)
      endif
!
    endsubroutine rprint_energy
!***********************************************************************
    subroutine fill_farray_pressure(f)
!
!  18-feb-10/anders: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine fill_farray_pressure
!***********************************************************************
    subroutine impose_energy_floor(f)
!
!  Dummy subroutine
!
      real, dimension(mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_energy_floor
!***********************************************************************
    subroutine dynamical_thermal_diffusion(uc)
!
!  Dummy subroutine
!
      real, intent(in) :: uc
!
      call keep_compiler_quiet(uc)
!
    endsubroutine dynamical_thermal_diffusion
!***********************************************************************
    subroutine split_update_energy(f)
!
!  Dummy subroutine
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine
!***********************************************************************
    subroutine expand_shands_energy
!
!  Presently dummy, for possible use
!
    endsubroutine expand_shands_energy
!***********************************************************************
    subroutine energy_after_timestep(f,df,dtsub)
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real :: dtsub
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dtsub)
!
    endsubroutine energy_after_timestep
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
           'characteristic velocity not yet implemented for entropy_const')

    endsubroutine update_char_vel_energy
!***********************************************************************
endmodule Energy

