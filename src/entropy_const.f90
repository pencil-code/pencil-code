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
! CPARAM logical, parameter :: lentropy = .true.
! CPARAM logical, parameter :: ltemperature = .false.
! CPARAM logical, parameter :: lthermal_energy = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED Ma2; fpres(3); sglnTT(3); advec_cs2
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
!
  ! input parameters
  namelist /entropy_init_pars/ &
      initss, ss_const
!
  ! run parameters
  namelist /entropy_run_pars/ &
      dummy
!
! diagnostics (needs to be consistent with reset list below)
!
  integer :: idiag_dtc=0,idiag_ethm=0,idiag_ethdivum=0,idiag_ssm=0
  integer :: idiag_ugradpm=0,idiag_ethtot=0,idiag_dtchi=0,idiag_ssmphi=0
!
  real, dimension(3) :: beta_glnrho_global=0.

  contains
!***********************************************************************
    subroutine register_energy
!
!  initialise variables which should know that we solve an energy
!  equation: iss, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
!  identify version number
!
      use SharedVariables, only: put_shared_variable

      if (lroot) call svn_id( &
           "$Id$")
!
      if (.not.ldensity.or.lboussinesq) call put_shared_variable('beta_glnrho_global',beta_glnrho_global)

    endsubroutine register_energy
!***********************************************************************
    subroutine initialize_energy(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: ierr
!
! run parameters
!
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
      call keep_compiler_quiet(f)
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
      real, dimension (mx,my,mz,mfarray) :: f
!
      integer :: j

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
    subroutine energy_before_boundary(f)

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f

      call keep_compiler_quiet(f)

    endsubroutine energy_before_boundary
!***********************************************************************
    subroutine energy_after_boundary(f)

      real, dimension (mx,my,mz,mfarray), intent(INOUT) :: f

      call keep_compiler_quiet(f)

      if (lenergy_slope_limited) &
        call not_implemented('energy_after_boundary','slope-limited diffusion')

    endsubroutine energy_after_boundary
!***********************************************************************
    subroutine pencil_criteria_energy
!
    endsubroutine pencil_criteria_energy
!***********************************************************************
    subroutine pencil_interdep_energy(lpencil_in)
!
!  Interdependency among pencils from the Energy module is specified here.
!
      logical, dimension (npencils) :: lpencil_in

      call keep_compiler_quiet(lpencil_in)

    endsubroutine pencil_interdep_energy
!***********************************************************************
    subroutine calc_pencils_energy(f,p)
!
!  Calculate Energy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) then
        p%advec_cs2=cs2*dxyz_2
        if (headtt.or.ldebug) print*,'calc_pencils_energy: max(advec_cs2) =',maxval(p%advec_cs2)
      endif

      call keep_compiler_quiet(f)

    endsubroutine calc_pencils_energy
!***********************************************************************
    subroutine denergy_dt(f,df,p)
!
!  28-mar-02/axel: dummy routine, adapted from entropy.f90 of 6-nov-01.
!  19-may-02/axel: added isothermal pressure gradient
!   9-jun-02/axel: pressure gradient added to du/dt already here
!
      use DensityMethods, only: getlnrho
      use EquationOfState, only: pressure_gradient
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type(pencil_case) :: p
 
      real, dimension (nx,3) :: glnrho
      real, dimension (nx) :: lnrho,cs2,cp1tilde
      integer :: j,ju
!
      intent(in) :: f
!
!  sound speed squared and inverse temperature
!
      call pressure_gradient(f,cs2,cp1tilde)
      call getlnrho(f(:,m,n,ilnrho),lnrho)

      if (headtt) print*, it, m, n, maxval(cs2), maxval(profxss)
!
!  subtract isothermal/polytropic pressure gradient term in momentum equation
!
      if (lhydro) then
        df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)-cs2*(glnrho(:,1)+dprofxss)
        df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)-cs2*(glnrho(:,2))
        df(l1:l2,m,n,iuz)=df(l1:l2,m,n,iuz)-cs2*(glnrho(:,3))
      endif

      if (lfirst.and.ldt) advec_cs2 = p%advec_cs2

      call calc_diagnostics_energy(f,p)

    endsubroutine denergy_dt
!***********************************************************************
    subroutine calc_diagnostics_energy(f,p)

      use Diagnostics
      use Sub, only: dot_mn

      real, dimension (mx,my,mz,mfarray) :: f
      type(pencil_case) :: p
!
!  Calculate energy related diagnostics
!
      if (ldiagnos) then
        if (lfirst.and.ldt.and.idiag_dtc/=0) call max_mn_name(sqrt(p%advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        if (idiag_ugradpm/=0) then
          call dot_mn(p%uu,p%glnrho,p%uglnrho)
          call sum_mn_name(p%rho*p%cs2*p%uglnrho,idiag_ugradpm)
        endif
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_diagnostics_energy
!***********************************************************************
    subroutine update_char_vel_energy(f)
!
!  Updates characteristic velocity for slope-limited diffusion.

      real, dimension (mx,my,mz,mfarray) :: f

      call keep_compiler_quiet(f)

    endsubroutine update_char_vel_energy
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
        call parse_name(iname,cname(iname),cform(iname),'ethdivum',idiag_ethdivum)
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

