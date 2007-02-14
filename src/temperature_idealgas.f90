! $Id: temperature_idealgas.f90,v 1.1 2007-02-14 13:04:17 wlyra Exp $

!  This module replaces the entropy module by using lnT as dependent
!  variable. For a perfect gas with constant coefficients (no ionization)
!  we have (1-1/gamma) * cp*T = cs02 * exp( (gamma-1)*ln(rho/rho0)-gamma*s/cp )
!  At a later point we may want to rename the module Entropy into Energy

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lentropy = .false.
! CPARAM logical, parameter :: ltemperature = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED Ma2,uglnTT,fpres
!
!***************************************************************

module Entropy

  use Cparam
  use Cdata
  use Messages
  use Interstellar

  implicit none

  include 'entropy.h'

  real :: radius_lnTT=0.1,ampl_lnTT=0.,widthlnTT=2*epsi
  real :: lnTT_left=1.0,lnTT_right=1.0,lnTT_const=0.0,TT_const=1.0
  real :: kx_lnTT=1.0,ky_lnTT=1.0,kz_lnTT=1.0
  real :: chi=0.0,heat_uniform=0.0
  real :: zbot=0.0,ztop=0.0
  real :: tau_heat_cor=-1.0,tau_damp_cor=-1.0,zcor=0.0,TT_cor=0.0
  logical :: lpressuregradient_gas=.true.,ladvection_temperature=.true.
  logical :: lupw_lnTT=.false.,lcalc_heat_cool=.false.
  logical :: lheatc_chiconst=.false.,lheatc_chiconst_accurate=.false.
  character (len=labellen), dimension(ninit) :: initlnTT='nothing'
  character (len=4) :: iinit_str

! Delete (or use) me asap!                                                       
  real :: hcond0,hcond1,Fbot,FbotKbot,Ftop,Kbot,FtopKtop
  logical :: lmultilayer

  ! input parameters
  namelist /entropy_init_pars/ &
      initlnTT,radius_lnTT,ampl_lnTT,widthlnTT, &
      lnTT_left,lnTT_right,lnTT_const,TT_const, &
      kx_lnTT,ky_lnTT,kz_lnTT

  ! run parameters
  namelist /entropy_run_pars/ &
       lupw_lnTT,lpressuregradient_gas,ladvection_temperature, &
      heat_uniform,chi,tau_heat_cor,tau_damp_cor,zcor,TT_cor, &
      lheatc_chiconst_accurate,lcalc_heat_cool
  ! other variables (needs to be consistent with reset list below)
  integer :: idiag_TTmax=0,idiag_TTmin=0,idiag_TTm=0
  integer :: idiag_yHmax=0,idiag_yHmin=0,idiag_yHm=0
  integer :: idiag_eth=0,idiag_ssm=0
  integer :: idiag_dtchi=0,idiag_dtc=0
  integer :: idiag_eem=0,idiag_ppm=0,idiag_csm=0

  contains

!***********************************************************************
    subroutine register_entropy()
!
!  initialise variables which should know that we solve an entropy
!  equation: ilnTT, etc; increase nvar accordingly
!
!  6-nov-01/wolf: coded
!
      logical, save :: first=.true.
!
      if (.not. first) call fatal_error('register_entropy','module registration called twice')
      first = .false.
!
      ilnTT = nvar+1             ! index to access temperature
      nvar = nvar+1
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_entropy: nvar = ', nvar
        print*, 'register_entropy: ilnTT = ', ilnTT
      endif
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: temperature_idealgas.f90,v 1.1 2007-02-14 13:04:17 wlyra Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call fatal_error('register_entropy','nvar > mvar')
      endif
!
!  Put variable name in array
!
      varname(ilnTT) = 'lnTT'
!
!  Writing files for use with IDL
!
      if (lroot) then
         if (maux == 0) then
            if (nvar < mvar) write(4,*) ',lnTT $'
            if (nvar == mvar) write(4,*) ',lnTT'
         else
            write(4,*) ',lnTT $'
         endif
         write(15,*) 'lnTT = fltarr(mx,my,mz)*one'
      endif
!                                       
    endsubroutine register_entropy
!***********************************************************************
    subroutine initialize_entropy(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-jul-2002/wolf: coded
!
      use Cdata
      use FArrayManager
      use Gravity, only: g0
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
      type (pencil_case) :: p
!
      if (.not. leos) then
         call fatal_error('initialize_entropy','EOS=noeos but entropy requires an EQUATION OF STATE for the fluid')
      endif
!
      call select_eos_variable('lnTT',ilnTT)
!
    endsubroutine initialize_entropy
!***********************************************************************
    subroutine init_ss(f,xx,yy,zz)
!
!  initialise lnTT; called from start.f90
!
!  13-dec-2002/axel+tobi: adapted from init_ss
!
!  initialise entropy; called from start.f90
!  07-nov-2001/wolf: coded 
!  24-nov-2002/tony: renamed for consistancy (i.e. init_[variable name])
!
      use General, only: chn
      use Sub, only: blob
      use Initcond, only: jump
!
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mx,my,mz), intent (in) :: xx,yy,zz
      logical :: lnothing=.true.
!
      do iinit=1,ninit
!
      if (initlnTT(iinit)/='nothing') then
!
      lnothing=.false.

      call chn(iinit,iinit_str)
!
!  select different initial conditions
!
      select case(initlnTT(iinit))
      case('zero', '0'); f(:,:,:,ilnTT) = 0.
      case('const_lnTT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+lnTT_const
      case('const_TT'); f(:,:,:,ilnTT)=f(:,:,:,ilnTT)+log(TT_const)
         
         !
        case default
          !
          !  Catch unknown values
          !
          write(unit=errormsg,fmt=*) 'No such value for initss(' &
                           //trim(iinit_str)//'): ',trim(initlnTT(iinit))
          call fatal_error('init_ss',errormsg)

      endselect

      if (lroot) print*,'init_ss: initss(' &
                        //trim(iinit_str)//') = ',trim(initlnTT(iinit))
      endif
      enddo
      if (lnothing.and.lroot) print*,'init_ss: nothing'
!
      if (NO_WARN) print*,xx,yy  !(to keep compiler quiet)        

    endsubroutine init_ss
!***********************************************************************
    subroutine pencil_criteria_entropy()
!
      use Cdata
!
      if (ldt) lpenc_requested(i_cs2)=.true.
      if (lpressuregradient_gas) lpenc_requested(i_fpres)=.true.
!
      if (lviscosity) then
        lpenc_requested(i_cv1)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_visc_heat)=.true.
      endif
!
      !if (ldensity) then
      !   lpenc_requested(i_delta)=.true.
      !   lpenc_requested(i_divu)=.true.
      !endif
!
      if (lcalc_heat_cool) then
        lpenc_requested(i_rho1)=.true.
        lpenc_requested(i_TT)=.true.
        lpenc_requested(i_TT1)=.true.
        lpenc_requested(i_cv1)=.true.
      endif
!
      if (ladvection_temperature) lpenc_requested(i_uglnTT)=.true.
!
      !if (lheatc_shock) then
      !   lpenc_requested(i_glnrho)=.true.
      !   lpenc_requested(i_gss)=.true.
      !   lpenc_requested(i_del2lnTT)=.true.
      !   lpenc_requested(i_gshock)=.true.
      !   lpenc_requested(i_shock)=.true.
      !   lpenc_requested(i_glnTT)=.true.
      !endif
!
!  Diagnostics
!
      if (idiag_TTmax/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_TTmin/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_TTm/=0) lpenc_diagnos(i_TT)=.true.
      if (idiag_yHmax/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_yHmin/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_yHm/=0) lpenc_diagnos(i_yH)=.true.
      if (idiag_eth/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_ee)=.true.
      endif
      if (idiag_ssm/=0) lpenc_diagnos(i_ss)=.true.
      if (idiag_dtchi/=0) then
        lpenc_diagnos(i_rho1)=.true.
        lpenc_diagnos(i_cv1)=.true.
      endif
      if (idiag_dtchi/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_csm/=0) lpenc_diagnos(i_cs2)=.true.
      if (idiag_eem/=0) lpenc_diagnos(i_ee)=.true.
      if (idiag_ppm/=0) lpenc_diagnos(i_pp)=.true.
!
    endsubroutine pencil_criteria_entropy
!***********************************************************************
    subroutine pencil_interdep_entropy(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_Ma2)) then
        lpencil_in(i_u2)=.true.
        lpencil_in(i_cs2)=.true.
      endif

      if (lpencil_in(i_uglnTT)) then
        lpencil_in(i_glnTT)=.true.
      endif
!
      if (lpencil_in(i_fpres)) then
        lpencil_in(i_cs2)=.true.
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_gss)=.true.
        lpencil_in(i_glnTT)=.true.
      endif
!
    endsubroutine pencil_interdep_entropy
!*********************************************************************** 
    subroutine calc_pencils_entropy(f,p)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  20-11-04/anders: coded
! 
      use EquationOfState
      use Sub

      real, dimension (mx,my,mz,mfarray), intent (in) :: f
      type (pencil_case), intent (inout) :: p
      integer :: j
!
!  Mach Speed
!
      if (lpencil(i_Ma2)) p%Ma2=p%u2/p%cs2
!
!  Temperature advection
!  (Needs to be here because of lupw_lnTT)
!
      if (lpencil(i_uglnTT)) then
        call u_dot_grad(f,ilnTT,p%glnTT,p%uu,p%uglnTT,UPWIND=lupw_lnTT)
      endif
!
      ! fpres
      if (lpencil(i_fpres)) then
        if (leos_idealgas) then
          do j=1,3
            p%fpres(:,j)=-p%cs2*(p%glnrho(:,j) + p%glnTT(:,j))*gamma11
          enddo
!  TH: The following would work if one uncomments the intrinsic operator
!  extensions in sub.f90. Please Test.
!  p%fpres      =-p%cs2*(p%glnrho + p%glnTT)*gamma11
        else
          do j=1,3
            p%fpres(:,j)=-p%cs2*(p%glnrho(:,j) + p%cp1tilde*p%gss(:,j))
          enddo
        endif
      endif
!
    endsubroutine calc_pencils_entropy
!**********************************************************************
    subroutine dss_dt(f,df,p)
!
!  calculate right hand side of entropy equation
!  heat condution is currently disabled until old stuff,
!  which in now in calc_heatcond, has been reinstalled.
!  DlnTT/Dt = -gamma1*divu + gamma*TT1*RHS
!
!  13-dec-02/axel+tobi: adapted from entropy
!
      use Cdata
      use Mpicomm
      use Sub
      use Global
      use Viscosity, only: calc_viscous_heat
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: Hmax=0.
      integer :: j,ju
!
      intent(inout) :: f,p
      intent(out) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'SOLVE dlnTT_dt'
      if (headtt) call identify_bcs('lnTT',ilnTT)
      if (headtt) print*,'dss_dt: lnTT,cs2=', p%lnTT(1), p%cs2(1)
!
!  entropy gradient: needed for advection and pressure gradient
!
      !call grad(f,ilnTT,glnTT)
!
!  sound speed squared
!
      if (headtt) print*,'dss_dt: cs20=',p%cs2(1)
!
!  ``cs2/dx^2'' for timestep
!
      if (lfirst.and.ldt) advec_cs2=p%cs2*dxyz_2
      if (headtt.or.ldebug) print*,'dss_dt: max(advec_cs2) =',maxval(advec_cs2)
!
!  subtract pressure gradient term in momentum equation
!
      if (lhydro) then
         if (lpressuregradient_gas) then
            do j=1,3
               ju=j+iuu-1
               df(l1:l2,m,n,ju) = df(l1:l2,m,n,ju) + p%fpres(:,j)
            enddo
         endif
      endif
!
!  advection term and PdV-work
!
        if (ladvection_temperature) then
        df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%uglnTT
      endif
!
!  Calculate viscous contribution to temperature
!
      if (lviscosity) call calc_viscous_heat(f,df,p,Hmax)
!
!  Various heating conduction contributions
!
      if (lcalc_heat_cool)  call calc_heat_cool(f,df,p)
!
!  Need to add left-hand-side of the continuity equation (see manual)
!  Check this

      !if (ldensity) then
      !  df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%gamma1*p%divu/p%delta
      !endif
!
!  Calculate entropy related diagnostics
!
      if (ldiagnos) then
        if (idiag_TTmax/=0) call max_mn_name(p%TT,idiag_TTmax)
        if (idiag_TTmin/=0) call max_mn_name(-p%TT,idiag_TTmin,lneg=.true.)
        if (idiag_TTm/=0)   call sum_mn_name(p%TT,idiag_TTm)
        if (idiag_yHmax/=0) call max_mn_name(p%yH,idiag_yHmax)
        if (idiag_yHmin/=0) call max_mn_name(-p%yH,idiag_yHmin,lneg=.true.)
        if (idiag_yHm/=0)   call sum_mn_name(p%yH,idiag_yHm)
        if (idiag_eth/=0)   call sum_mn_name(p%ee/p%rho1,idiag_eth)
        if (idiag_ssm/=0)   call sum_mn_name(p%ss,idiag_ssm)
        if (idiag_dtc/=0) then
          call max_mn_name(sqrt(advec_cs2)/cdt,idiag_dtc,l_dt=.true.)
        endif
        if (idiag_eem/=0) call sum_mn_name(p%ee,idiag_eem)
        if (idiag_ppm/=0) call sum_mn_name(p%pp,idiag_ppm)
        if (idiag_csm/=0) call sum_mn_name(p%cs2,idiag_csm,lsqrt=.true.)
      endif
!
    endsubroutine dss_dt
!***********************************************************************
    subroutine calc_heat_cool(f,df,p)

      use EquationOfState, only: gamma,gamma1
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx) :: tau,cooling,kappa,a1,a3,pTT
      real :: a2,kappa0,kappa0_cgs
!
!  Initialize
!
      intent(in) :: p
      intent(out) :: df
!
      if (headtt) print*,'enter calc_heat_cool'
!
      kappa0_cgs=2e-6  !cm2/g
      kappa0=kappa0_cgs*unit_density*unit_length
      kappa=kappa0*p%TT**2
!
!  Optical Depth tau=kappa*rho*H
!  If we are using 2D, the pencil value p%rho is actually
!   sigma, the column density, sigma=rho*2*H
!
      if (nzgrid==1) then
         tau = .5*kappa*p%rho
      else
         call fatal_error("calc_heat_cool","opacity not yet implemented for 3D")
      endif
!
! Analytical gray description of Hubeny (1990)
! a1 is the optically thick contribution,
! a3 the optically thin one.
!
      a1=0.375*tau ; a2=0.433013 ; a3=0.25/tau
!
! cooling for Energy: 2*sigmaSB*p%TT**4/(a1+a2+a3)
!
      cooling = 2*sigmaSB*p%rho1*p%TT**4/(a1) !+a2+a3)
!
!  this cooling has dimension of energy over time
!
      df(l1:l2,m,n,ilnTT) = df(l1:l2,m,n,ilnTT) - p%cv1*p%TT1*cooling
!
    endsubroutine calc_heat_cool
!***********************************************************************
    subroutine read_entropy_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=entropy_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=entropy_init_pars,ERR=99)
      endif


99    return
    endsubroutine read_entropy_init_pars
!***********************************************************************
    subroutine write_entropy_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=entropy_init_pars)

    endsubroutine write_entropy_init_pars
!***********************************************************************
    subroutine read_entropy_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=entropy_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=entropy_run_pars,ERR=99)
      endif


99    return
    endsubroutine read_entropy_run_pars
!***********************************************************************
    subroutine write_entropy_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=entropy_run_pars)

    endsubroutine write_entropy_run_pars
!***********************************************************************
    subroutine rprint_entropy(lreset,lwrite)
!
!  reads and registers print parameters relevant to entropy
!
!   1-jun-02/axel: adapted from magnetic fields
!
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
        idiag_TTmax=0; idiag_TTmin=0; idiag_TTm=0
        idiag_yHmax=0; idiag_yHmin=0; idiag_yHm=0
        idiag_eth=0; idiag_ssm=0
        idiag_dtchi=0; idiag_dtc=0
        idiag_eem=0; idiag_ppm=0; idiag_csm=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'TTmax',idiag_TTmax)
        call parse_name(iname,cname(iname),cform(iname),'TTmin',idiag_TTmin)
        call parse_name(iname,cname(iname),cform(iname),'TTm',idiag_TTm)
        call parse_name(iname,cname(iname),cform(iname),'yHmax',idiag_yHmax)
        call parse_name(iname,cname(iname),cform(iname),'yHmin',idiag_yHmin)
        call parse_name(iname,cname(iname),cform(iname),'yHm',idiag_yHm)
        call parse_name(iname,cname(iname),cform(iname),'eth',idiag_eth)
        call parse_name(iname,cname(iname),cform(iname),'ssm',idiag_ssm)
        call parse_name(iname,cname(iname),cform(iname),'dtchi',idiag_dtchi)
        call parse_name(iname,cname(iname),cform(iname),'dtc',idiag_dtc)
        call parse_name(iname,cname(iname),cform(iname),'eem',idiag_eem)
        call parse_name(iname,cname(iname),cform(iname),'ppm',idiag_ppm)
        call parse_name(iname,cname(iname),cform(iname),'csm',idiag_csm)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'ilnTT=',ilnTT
        write(3,*) 'iyH=',iyH
        write(3,*) 'iss=',iss
        write(3,*) 'i_TTmax=',idiag_TTmax
        write(3,*) 'i_TTmin=',idiag_TTmin
        write(3,*) 'i_TTm=',idiag_TTm
        write(3,*) 'i_yHmax=',idiag_yHmax
        write(3,*) 'i_yHmin=',idiag_yHmin
        write(3,*) 'i_yHm=',idiag_yHm
        write(3,*) 'i_eth=',idiag_eth
        write(3,*) 'i_ssm=',idiag_ssm
        write(3,*) 'i_dtchi=',idiag_dtchi
        write(3,*) 'i_dtc=',idiag_dtc
        write(3,*) 'i_eem=',idiag_eem
        write(3,*) 'i_ppm=',idiag_ppm
        write(3,*) 'i_csm=',idiag_csm
      endif
!
    endsubroutine rprint_entropy
!***********************************************************************
endmodule Entropy
