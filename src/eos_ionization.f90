! $Id$
!
!  This modules contains the routines for simulation with
!  simple hydrogen ionization.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .true., leos_ionization = .true., leos_temperature_ionization=.false.
! CPARAM logical, parameter :: leos_idealgas = .false., leos_chemistry = .false.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 2
! COMMUNICATED AUXILIARIES 1
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cp; cp1; cp1tilde
! PENCILS PROVIDED glnTT(3); TT; TT1; gTT(3); yH; hss(3,3); hlnTT(3,3); del2TT; del6TT; del6lnTT
! PENCILS PROVIDED del2ss; del6ss; del2lnTT; cv; cv1; glnmumol(3); ppvap; csvap2
! PENCILS PROVIDED rho_anel; rho1gpp(3)
!
!***************************************************************
module EquationOfState
!
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'eos.h'
  include 'eos_params.h'

!  secondary parameters calculated in initialize
  real :: TT_ion,lnTT_ion,TT_ion_,lnTT_ion_
  real :: ss_ion,ee_ion,kappa0,xHe_term,ss_ion1,Srad0
  real :: lnrho_e,lnrho_e_,lnrho_H,lnrho_He,Rgas,mu1yHxHe
! namelist parameters
  real, parameter :: yHmin=tiny(TT_ion), yHmax=1-epsilon(TT_ion)
  real :: xHe=0.1
  real :: yMetals=0
  real :: yHacc=1e-5
!ajwm  Moved here from Density.f90
!ajwm  Completely irrelevant to eos_ionization but density and entropy need
!ajwm  reworking to be independent of these things first
!ajwm  can't use impossible else it breaks reading param.nml
!ajwm  SHOULDN'T BE HERE... But can wait till fully unwrapped
  real :: cs0=impossible, rho0=impossible
  real :: cs20=impossible, lnrho0=impossible
  logical :: lpp_as_aux=.false., lcp_as_aux=.false.
  real :: gamma=impossible
  real :: lnTT0=impossible, TT0=impossible
  real :: cs2bot=impossible, cs2top=impossible
! input parameters
  namelist /eos_init_pars/ xHe,yMetals,yHacc,lpp_as_aux,lcp_as_aux
! run parameters
  namelist /eos_run_pars/ xHe,yMetals,yHacc,lpp_as_aux,lcp_as_aux
!
  integer :: imass=0, ivars_mod
!
  real :: Cp_const=impossible
  real :: Pr_number=0.7
  logical :: lpres_grad=.false.
!
  contains
!***********************************************************************
    subroutine register_eos
!
!   2-feb-03/axel: adapted from Interstellar module
!   13-jun-03/tobi: re-adapted from visc_shock module
!
      use FArrayManager
      use Sub
      use SharedVariables,only: put_shared_variable
!
!  Set indices for auxiliary variables.
!
      !call farray_register_auxiliary('yH',iyH,communicated=.true.)
      call farray_register_auxiliary('yH',iyH)
      if (.not.ltemperature.or.ltemperature_nolog) &
        call farray_register_auxiliary('lnTT',ilnTT,communicated=.true.)
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Writing files for use with IDL.
!
      aux_var(aux_count)=',yh $'
      aux_count=aux_count+1
      if (naux < maux)  aux_var(aux_count)=',lnTT $'
      if (naux == maux) aux_var(aux_count)=',lnTT'
      aux_count=aux_count+1
      if (lroot) then
        write(15,*) 'yH = fltarr(mx,my,mz)*one'
        write(15,*) 'lnTT = fltarr(mx,my,mz)*one'
      endif

      call put_shared_variable('gamma',gamma,caller='register_eos')

      if (.not.ldensity) then
        call put_shared_variable('rho0',rho0)
        call put_shared_variable('lnrho0',lnrho0)
      endif
!
    endsubroutine register_eos
!***********************************************************************
    subroutine units_eos
!
!  If unit_temperature hasn't been specified explictly in start.in,
!  set it to 1 (Kelvin).
!
!  24-jun-06/tobi: coded
!
!  11-feb-23/fred lfix_unit_std seeks to adopt a numerically stable value
!                 as in ideal gas, yet to explore here
!
      if (unit_temperature==impossible) then
        if (lfix_unit_std) then
          unit_temperature=unit_density*unit_velocity**2/k_B_cgs*13.6
        else
          unit_temperature=1.
        endif
      endif
      if (lroot) print*,'unit temperature',unit_temperature
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos(f)
!
!  Perform any post-parameter-read initialization, e.g. set derived
!  parameters.
!
!   2-feb-03/axel: adapted from Interstellar module
!
      use General
      use Sub, only: register_report_aux
!
      real, dimension (mx,my,mz,mfarray),intent(INOUT):: f

      if (lroot) print*,'initialize_eos: ENTER'
!
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
      if (pretend_lnTT) then
        call warning('initialize_eos','pretend_lnTT is not used with ionization')
        pretend_lnTT=.false.
      endif
      Rgas=k_B/m_p
      mu1yHxHe=1+3.97153*xHe
      TT_ion=chiH/k_B
      lnTT_ion=log(TT_ion)
      TT_ion_=chiH_/k_B
      lnTT_ion_=log(chiH_/k_B)
      lnrho_e=1.5*log((m_e/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_H=1.5*log((m_H/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_He=1.5*log((m_He/hbar)*(chiH/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      lnrho_e_=1.5*log((m_e/hbar)*(chiH_/hbar)/2./pi)+log(m_H)+log(mu1yHxHe)
      ss_ion=k_B/m_H/mu1yHxHe
      ss_ion1=1/ss_ion
      ee_ion=ss_ion*TT_ion
      kappa0=sigmaH_/m_H/mu1yHxHe/4.0
      Srad0=sigmaSB*TT_ion**4.0D0/pi
!
      if (xHe>0) then
        xHe_term=xHe*(log(xHe)-lnrho_He)
      elseif (xHe<0) then
        call fatal_error('initialize_eos','xHe < 0 makes no sense')
      else
        xHe_term=0
      endif
!
!  pressure and cp as optional auxiliary variable
!
      if (lpp_as_aux) call register_report_aux('pp',ipp)
      if (lcp_as_aux) call register_report_aux('cp',icp)
!
      if (lrun) call init_eos(f)

      if (lroot.and.ip<14) then
        print*,'initialize_eos: reference values for ionization'
        print*,'initialize_eos: yHmin,yHmax,yMetals=',yHmin,yHmax,yMetals
        print*,'initialize_eos: TT_ion,ss_ion,kappa0=', &
                TT_ion,ss_ion,kappa0
        print*,'initialize_eos: lnrho_e,lnrho_H,lnrho_He,lnrho_e_=', &
                lnrho_e,lnrho_H,lnrho_He,lnrho_e_
      endif
!
!  write scale non-free constants to file; to be read by idl
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,*) 'TT_ion=',TT_ion
        write (1,*) 'lnTT_ion=',lnTT_ion
        write (1,*) 'TT_ion_=',TT_ion_
        write (1,*) 'lnTT_ion_=',lnTT_ion_
        write (1,*) 'lnrho_e=',lnrho_e
        write (1,*) 'lnrho_H=',lnrho_H
        write (1,*) 'lnrho_p=',lnrho_H
        write (1,*) 'lnrho_He=',lnrho_He
        write (1,*) 'lnrho_e_=',lnrho_e_
        write (1,*) 'ss_ion=',ss_ion
        write (1,*) 'ee_ion=',ee_ion
        write (1,*) 'kappa0=',kappa0
        write (1,*) 'Srad0=',Srad0
        write (1,*) 'k_B=',k_B
        write (1,*) 'm_H=',m_H
        close (1)
      endif
!
    endsubroutine initialize_eos
!***********************************************************************
    subroutine select_eos_variable(variable,findex)
!
!  Calculate average particle mass in the gas relative to
!
!   02-apr-06/tony: implemented
!
      character (len=*), intent(in) :: variable
      integer, intent(in) :: findex
!
      call keep_compiler_quiet(variable)
      call keep_compiler_quiet(findex)
!  DUMMY ideagas version below
!!      integer :: this_var=0
!!      integer, save :: ieosvar=0
!!      integer, save :: ieosvar_selected=0
!!      integer, parameter :: ieosvar_lnrho = 1
!!      integer, parameter :: ieosvar_rho   = 2
!!      integer, parameter :: ieosvar_ss    = 4
!!      integer, parameter :: ieosvar_lnTT  = 8
!!      integer, parameter :: ieosvar_TT    = 16
!!      integer, parameter :: ieosvar_cs2   = 32
!!      integer, parameter :: ieosvar_pp    = 64
!!!
!!      if (ieosvar>=2) &
!!        call fatal_error("select_eos_variable", &
!!             "2 thermodynamic quantities have already been defined while attempting to add a 3rd: ") !//variable)
!!
!!      ieosvar=ieosvar+1
!!
!!!      select case (variable)
!!      if (variable=='ss') then
!!          this_var=ieosvar_ss
!!          if (findex<0) then
!!            leos_isentropic=.true.
!!          endif
!!      elseif (variable=='cs2') then
!!          this_var=ieosvar_cs2
!!          if (findex==-2) then
!!            leos_localisothermal=.true.
!!          elseif (findex<0) then
!!            leos_isothermal=.true.
!!          endif
!!      elseif (variable=='lnTT') then
!!          this_var=ieosvar_lnTT
!!          if (findex<0) then
!!            leos_isothermal=.true.
!!          endif
!!      elseif (variable=='lnrho') then
!!          this_var=ieosvar_lnrho
!!          if (findex<0) then
!!            leos_isochoric=.true.
!!          endif
!!      elseif (variable=='rho') then
!!          this_var=ieosvar_rho
!!          if (findex<0) then
!!            leos_isochoric=.true.
!!          endif
!!      elseif (variable=='pp') then
!!          this_var=ieosvar_pp
!!          if (findex<0) then
!!            leos_isobaric=.true.
!!          endif
!!      else
!!        call fatal_error("select_eos_variable", &
!!             "unknown thermodynamic variable")
!!      endif
!!      if (ieosvar==1) then
!!        ieosvar1=findex
!!        ieosvar_selected=ieosvar_selected+this_var
!!        return
!!      endif
!!!
!!! Ensure the indexes are in the correct order.
!!!
!!      if (this_var<ieosvar_selected) then
!!        ieosvar2=ieosvar1
!!        ieosvar1=findex
!!      else
!!        ieosvar2=findex
!!      endif
!!      ieosvar_selected=ieosvar_selected+this_var
!!      select case (ieosvar_selected)
!!        case (ieosvar_lnrho+ieosvar_ss)
!!          ieosvars=ilnrho_ss
!!        case (ieosvar_lnrho+ieosvar_lnTT)
!!          ieosvars=ilnrho_lnTT
!!        case (ieosvar_lnrho+ieosvar_cs2)
!!          ieosvars=ilnrho_cs2
!!        case default
!!          print*,"select_eos_variable: Thermodynamic variable combination, ieosvar_selected= ",ieosvar_selected
!!          call fatal_error("select_eos_variable", &
!!             "This thermodynamic variable combination is not implemented: ")
!!      endselect
!
    endsubroutine select_eos_variable
!***********************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='lnTT'.or.cnamev=='yH') cformv='DEFINED'
      endif
!
    endsubroutine rprint_eos
!***********************************************************************
    subroutine get_slices_eos(f,slices)
!
!  Write slices for animation of Eos variables.
!
!  26-jul-06/tony: coded
!
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices.
!
      select case (trim(slices%name))
!
!  Temperature.
!
        case ('lnTT'); call assign_slices_scal(slices,f,ilnTT)
!
!  Degree of ionization.
!
        case ('yH'); call assign_slices_scal(slices,f,iyH)
!
      endselect
!
    endsubroutine get_slices_eos
!***********************************************************************
    subroutine pencil_criteria_eos
!
!  All pencils that the EquationOfState module depends on are specified here.
!
!  02-apr-06/tony: coded
!
!  EOS is a pencil provider but evolves nothing so it is unlokely that
!  it will require any pencils for it's own use.
!
!  pp pencil if lpp_as_aux
!
      if (lpp_as_aux) lpenc_requested(i_pp)=.true.
      if (lcp_as_aux) lpenc_requested(i_cp1tilde)=.true.
!
    endsubroutine pencil_criteria_eos
!***********************************************************************
    subroutine pencil_interdep_eos(lpencil_in)
!
!  Interdependency among pencils from the Entropy module is specified here.
!
!  20-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_cv).or.lpencil_in(i_cp).or.&
          lpencil_in(i_cv1).or.lpencil_in(i_cp1)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_TT1)=.true.
      endif
!
      if (lpencil_in(i_gTT)) then
        lpencil_in(i_glnTT)=.true.
        lpencil_in(i_TT)=.true.
      endif
      !if (lpencil_in(i_del2lnTT)) then
      !  lpencil_in(i_del2lnrho)=.true.
      !  lpencil_in(i_del2ss)=.true.
      !endif
      !if (lpencil_in(i_glnTT)) then
      !  lpencil_in(i_glnrho)=.true.
      !  lpencil_in(i_gss)=.true.
      !endif
      if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
      if (lpencil_in(i_TT1)) lpencil_in(i_lnTT)=.true.
!
      if (lpencil_in(i_hlnTT)) then
        lpencil_in(i_hss)=.true.
        if (.not.pretend_lnTT) lpencil_in(i_hlnrho)=.true.
      endif
!
    endsubroutine pencil_interdep_eos
!***********************************************************************
    subroutine calc_pencils_eos_std(f,p)
!
! Envelope adjusting calc_pencils_eos_pencpar to the standard use with
! lpenc_loc=lpencil
!
!  9-oct-15/MR: coded
!
      real, dimension (mx,my,mz,mfarray),intent(INOUT):: f
      type (pencil_case),                intent(INOUT):: p
!
      call calc_pencils_eos_pencpar(f,p,lpencil)
!
    endsubroutine calc_pencils_eos_std
!***********************************************************************
    subroutine calc_pencils_eos_pencpar(f,p,lpenc_loc)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  02-apr-06/tony: coded
!  09-oct-15/MR: added mask parameter lpenc_loc
!  13-feb-23/FG: added pencils for cp to cv1
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      intent(inout):: f,p
      intent(in) :: lpenc_loc
!
      integer :: i
!
! THE FOLLOWING 2 ARE CONCEPTUALLY WRONG
! FOR pretend_lnTT since iss actually contain lnTT NOT entropy!
! The code is not wrong however since this is correctly
! handled by the eos module.
! ss
      if (lpenc_loc(i_ss)) p%ss=f(l1:l2,m,n,iss)
! gss
      if (lpenc_loc(i_gss)) call grad(f,iss,p%gss)
! pp
      if (lpenc_loc(i_pp)) call eoscalc(f,nx,pp=p%pp)
! ee
      if (lpenc_loc(i_ee)) call eoscalc(f,nx,ee=p%ee)
! lnTT
      if (lpenc_loc(i_lnTT)) p%lnTT=f(l1:l2,m,n,ilnTT)
!      if (lpenc_loc(i_lnTT)) call eoscalc(f,nx,lnTT=p%lnTT)
! yH
      if (lpenc_loc(i_yH)) p%yH=f(l1:l2,m,n,iyH)
!      if (lpenc_loc(i_yH)) call eoscalc(f,nx,yH=p%yH)
! TT
      if (lpenc_loc(i_TT)) p%TT=exp(p%lnTT)
! TT1
      if (lpenc_loc(i_TT1)) p%TT1=exp(-p%lnTT)
! cs2 and cp1tilde
      if (lpenc_loc(i_cs2) .or. lpenc_loc(i_cp1tilde)) call pressure_gradient(f,p%cs2,p%cp1tilde)
! glnTT
      if (lpenc_loc(i_glnTT)) call grad(f,ilnTT,p%glnTT)
!        call temperature_gradient(f,p%glnrho,p%gss,p%glnTT)
! gTT
      if (lpenc_loc(i_gTT)) then
        do i=1,3; p%gTT(:,i)=p%glnTT(:,i)*p%TT; enddo
      endif
! hss
      if (lpenc_loc(i_hss)) call g2ij(f,iss,p%hss)
! del2ss
      if (lpenc_loc(i_del2ss)) call del2(f,iss,p%del2ss)
! del2lnTT
      if (lpenc_loc(i_del2lnTT)) call del2(f,ilnTT,p%del2lnTT)
!        call temperature_laplacian(f,p)
! del2TT
      if (lpenc_loc(i_del2TT)) call not_implemented('calc_pencils_eos_pencpar','del2TT')
! del6lnTT
      if (lpenc_loc(i_del6lnTT)) call del6(f,ilnTT,p%del6lnTT)
! del6TT
      if (lpenc_loc(i_del6TT)) call not_implemented('calc_pencils_eos_pencpar','del6TT')
! del6ss
      if (lpenc_loc(i_del6ss)) call del6(f,iss,p%del6ss)
! hlnTT
      if (lpenc_loc(i_hlnTT)) call temperature_hessian(f,p%hlnrho,p%hss,p%hlnTT)
!
      if (lpenc_loc(i_glnmumol)) p%glnmumol(:,:)=0.
!
!  cv/cp pencils are computed following A&A 587, A90 (2016)
!  DOI: 10.1051/0004-6361/201425396
!  TBD Helium single and double ionization states to be included
!
! cp
      if (lpenc_loc(i_cp)) &
          p%cp=(2.5+p%yH*(1-p%yH)/((2-p%yH)*xHe+2)*(2.5+p%TT1*TT_ion)**2)* &
          Rgas*mu1yHxHe/(1+xHe+p%yH)
! cv
      if (lpenc_loc(i_cv)) &
          p%cv=(1.5+p%yH*(1-p%yH)/((2-p%yH)*(1+p%yH+xHe))*(1.5+p%TT1*TT_ion)**2)* &
          Rgas*mu1yHxHe/(1+xHe+p%yH)
! cp1
      if (lpenc_loc(i_cp1)) &
          p%cp1=(1+xHe+p%yH)/((2.5+p%yH*(1-p%yH)/((2-p%yH)*xHe+2)* &
          (2.5+p%TT1*TT_ion)**2)*Rgas*mu1yHxHe)
! cv1
      if (lpenc_loc(i_cv1)) &
          p%cv1=(1+xHe+p%yH)/((1.5+p%yH*(1-p%yH)/((2-p%yH)*(1+p%yH+xHe))* &
          (1.5+p%TT1*TT_ion)**2)*Rgas*mu1yHxHe)
!
!  pressure and cp as optional auxiliary pencils
!
      if (lpp_as_aux) f(l1:l2,m,n,ipp)=p%pp
      if (lcp_as_aux) f(l1:l2,m,n,icp)=p%cp1tilde
!
    endsubroutine calc_pencils_eos_pencpar
!***********************************************************************
    subroutine getmu(f,mu_tmp)
!
!  Calculate average particle mass.
!  Note that the particles density is N = nHI + nHII + ne + nHe
!  = (1-y)*nH + y*nH + y*nH + xHe*nH = (1 + yH + xHe) * nH, where
!  nH is the number of protons per cubic centimeter.
!  The number of particles per mole is therefore 1 + yH + xHe.
!  The mass per mole is M=1.+3.97153*xHe, so the mean molecular weight
!  per particle is M/N = (1.+3.97153*xHe)/(1 + yH + xHe).
!
!   12-aug-03/tony: implemented
!   13-feb-23/fred: mu_tmp now included in pencils for cp, etc.
!                   call fatal error for getmu
!
      real, dimension (mx,my,mz,mfarray), optional :: f
      real, optional, intent(out) :: mu_tmp
!
!  for variable ionization, it doesn't make sense to calculate
!  just a single value of mu, because it must depend on position.
!  Therefore, call fatal error.
!
      call not_implemented('getmu','for eos_ionization, use cp/cv pencils')
!
! tobi: the real mean molecular weight would be:
!
! mu_tmp=(1.+3.97153*xHe)/(1+yH+xHe)
!
    endsubroutine getmu
!***********************************************************************
    subroutine init_eos(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      f(:,:,:,iyH) = 0.5*(yHmax-yHmin)
      call ioncalc(f)
!
    endsubroutine init_eos
!***********************************************************************
    subroutine ioncalc(f)
!
!   calculate degree of ionization and temperature
!   This routine is called from equ.f90 and operates on the full 3-D array.
!
!   13-jun-03/tobi: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx) :: lnrho,ss,yH,lnTT
!
      do n=1,mz
      do m=1,my
        if (ldensity_nolog) then
          lnrho=log(f(:,m,n,ilnrho))
        else
          lnrho=f(:,m,n,ilnrho)
        endif
        ss=f(:,m,n,iss)
        yH=f(:,m,n,iyH)
        call rtsafe_pencil(lnrho,ss,yH)
        f(:,m,n,iyH)=yH
        lnTT=(ss/ss_ion+(1-yH)*(log(1-yH+epsi)-lnrho_H) &
              +yH*(2*log(yH)-lnrho_e-lnrho_H)+xHe_term)/(1+yH+xHe)
        lnTT=(2.0/3.0)*(lnTT+lnrho-2.5)+lnTT_ion
        f(:,m,n,ilnTT)=lnTT
      enddo
      enddo
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)
!
!  calculate density. Is currently only being used by the interstellar
!  module. I guess we can/should replace this now by a call to eoscalc.
!
      real, intent(in) :: EE,TT,yH
      real, intent(out) :: rho
!
      rho=EE/(1.5*(1.+yH+xHe)*ss_ion*TT+yH*ee_ion)
!
    endsubroutine getdensity
!***********************************************************************
    subroutine get_gamma_etc(gamma,cp,cv)
!
      real, intent(OUT) :: gamma
      real, optional, intent(OUT) :: cp,cv
!
      if (headt) call warning('get_gamma_etc','gamma, cp, and cv are not constant in eos_ionization.'// &
                              achar(10)//'The values provided are for one-atomic ideal gas. Use at own risk')
      gamma=5./3.
      if (present(cp)) cp=1.
      if (present(cv)) cv=3./5.

    endsubroutine get_gamma_etc
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      use General, only: ioptest
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: cs2
      real, dimension(nx), intent(out), optional :: cp1tilde

      real, dimension(nx) :: lnrho,yH,lnTT
      real, dimension(nx) :: R,dlnTTdy,dRdy,temp
      real, dimension(nx) :: dlnPPdlnrho,fractions,fractions1
      real, dimension(nx) :: dlnPPdss,TT1

      lnrho=f(l1:l2,m,n,ilnrho)
      yH=f(l1:l2,m,n,iyH)
      lnTT=f(l1:l2,m,n,ilnTT)
!
      TT1=exp(-lnTT)
      fractions=1+yH+xHe
      fractions1=1/fractions
!
      R=lnrho_e-lnrho+1.5*(lnTT-lnTT_ion)-TT_ion*TT1+log(1-yH+epsi)-2*log(yH)
      dlnTTdy=(2*(-R-TT_ion*TT1)-3)/3*fractions1
      dRdy=dlnTTdy*(1.5+TT_ion*TT1)-1/(1-yH+epsi)-2/yH
      temp=(dlnTTdy+fractions1)/dRdy
      dlnPPdlnrho=(5-2*TT_ion*TT1*temp)/3
      dlnPPdss=ss_ion1*fractions1*(dlnPPdlnrho-temp-1)
      cs2=fractions*ss_ion*dlnPPdlnrho/TT1
      if (present(cp1tilde)) cp1tilde=dlnPPdss/dlnPPdlnrho
!
    endsubroutine pressure_gradient_farray
!***********************************************************************
    subroutine pressure_gradient_point(lnrho,ss,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: cs2,cp1tilde
      real :: yH,lnTT
      real :: R,dlnTTdy,dRdy,temp
      real :: dlnPPdlnrho,fractions,fractions1
      real :: dlnPPdss,TT1
!
      yH=0.5
      call rtsafe(ilnrho_ss,lnrho,ss,yHmin,yHmax,yH)
      fractions=(1+yH+xHe)
      fractions1=1/fractions
      lnTT=(2.0/3.0)*((ss/ss_ion+(1-yH)*(log(1-yH+epsi)-lnrho_H) &
                       +yH*(2*log(yH)-lnrho_e-lnrho_H) &
                       +xHe_term)/fractions+lnrho-2.5)+lnTT_ion
      TT1=exp(-lnTT)
!
      R=lnrho_e-lnrho+1.5*(lnTT-lnTT_ion)-TT_ion*TT1+log(1-yH+epsi)-2*log(yH)
      dlnTTdy=(2*(-R-TT_ion*TT1)-3)/3*fractions1
      dRdy=dlnTTdy*(1.5+TT_ion*TT1)-1/(1-yH+epsi)-2/yH
      temp=(dlnTTdy+fractions1)/dRdy
      dlnPPdlnrho=(5-2*TT_ion*TT1*temp)/3
      dlnPPdss=ss_ion1*fractions1*(dlnPPdlnrho-temp-1)
      cs2=fractions*ss_ion*dlnPPdlnrho/TT1
      cp1tilde=dlnPPdss/dlnPPdlnrho
!
    endsubroutine pressure_gradient_point
!***********************************************************************
    subroutine temperature_gradient(f,glnrho,gss,glnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(in) :: glnrho,gss
      real, dimension(nx,3), intent(out) :: glnTT
      real, dimension(nx) :: lnrho,yH,lnTT,TT1,fractions1
      real, dimension(nx) :: R,dlnTTdy,dRdy
      real, dimension(nx) :: dlnTTdydRdy,dlnTTdlnrho,dlnTTdss
      integer :: j
!
      lnrho=f(l1:l2,m,n,ilnrho)
      yH=f(l1:l2,m,n,iyH)
      lnTT=f(l1:l2,m,n,ilnTT)
      TT1=exp(-lnTT)
      fractions1=1/(1+yH+xHe)
!
      R=lnrho_e-lnrho+1.5*(lnTT-lnTT_ion)-TT_ion*TT1+log(1-yH+epsi)-2*log(yH)
      dlnTTdy=((2.0/3.0)*(-R-TT_ion*TT1)-1)*fractions1
      dRdy=dlnTTdy*(1.5+TT_ion*TT1)-1/(1-yH+epsi)-2/yH
      dlnTTdydRdy=dlnTTdy/dRdy
      dlnTTdlnrho=(2.0/3.0)*(1-TT_ion*TT1*dlnTTdydRdy)
      dlnTTdss=(dlnTTdlnrho-dlnTTdydRdy)*fractions1*ss_ion1
      do j=1,3
        glnTT(:,j)=dlnTTdlnrho*glnrho(:,j)+dlnTTdss*gss(:,j)
      enddo
!
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally hlnPP and hlnTT
!   hP/rho=cs2*(hlnrho+cp1tilde*hss)
!
!   10-apr-04/axel: adapted from temperature_gradient
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3,3), intent(in) :: hlnrho,hss
      real, dimension(nx,3,3), intent(out) :: hlnTT
      real, dimension(nx) :: lnrho,yH,lnTT,TT1,fractions1
      real, dimension(nx) :: R,dlnTTdy,dRdy
      real, dimension(nx) :: dlnTTdydRdy,dlnTTdlnrho,dlnTTdss
      integer :: i,j
!
      lnrho=f(l1:l2,m,n,ilnrho)
      yH=f(l1:l2,m,n,iyH)
      lnTT=f(l1:l2,m,n,ilnTT)
      TT1=exp(-lnTT)
      fractions1=1/(1+yH+xHe)
!
      R=lnrho_e-lnrho+1.5*(lnTT-lnTT_ion)-TT_ion*TT1+log(1-yH+epsi)-2*log(yH)
      dlnTTdy=((2.0/3.0)*(-R-TT_ion*TT1)-1)*fractions1
      dRdy=dlnTTdy*(1.5+TT_ion*TT1)-1/(1-yH+epsi)-2/yH
      dlnTTdydRdy=dlnTTdy/dRdy
      dlnTTdlnrho=(2.0/3.0)*(1-TT_ion*TT1*dlnTTdydRdy)
      dlnTTdss=(dlnTTdlnrho-dlnTTdydRdy)*fractions1*ss_ion1
      do j=1,3
      do i=1,3
        hlnTT(:,i,j)=dlnTTdlnrho*hlnrho(:,i,j)+dlnTTdss*hss(:,i,j)
      enddo
      enddo
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eosperturb(f,psize,ee,pp,ss)
!
!  Set f(l1:l2,m,n,iss), depending on the values of ee and pp
!
!  20-jan-15/MR: changes for use of reference state
!
  use SharedVariables, only: get_shared_variable
  use DensityMethods, only: getlnrho,getrho
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(in), optional :: ee, pp, ss
!
      real, dimension(psize) :: lnrho_,ss_,rho_,TT_,lnTT_,yH_,fractions
      real, dimension(:,:), pointer :: reference_state
      integer :: i
!
      if (psize==nx) then
        if (lreference_state) &
          call get_shared_variable('reference_state',reference_state,caller='eosperturb')
        if (ldensity_nolog) then
          if (lreference_state) then
            call getrho(f(:,m,n,irho),rho_)
            lnrho_=log(rho_)+log(f(l1:l2,m,n,irho))
          else
            lnrho_=log(f(l1:l2,m,n,irho))
          endif
        else
          if (lreference_state) then
            call getlnrho(f(:,m,n,ilnrho),lnrho_)
            lnrho_=lnrho_+log(f(l1:l2,m,n,ilnrho))
          else
            lnrho_=f(l1:l2,m,n,ilnrho)
          endif
        endif
        if (present(ee)) then
          call eoscalc(ilnrho_ee,lnrho_,ee,ss=ss_)
        elseif (present(pp)) then
          call eoscalc(ilnrho_pp,lnrho_,pp,ss=ss_)
        elseif (present(ss)) then
          ss_=ss
        endif
!
        f(l1:l2,m,n,iss) = ss_
        if (lreference_state) f(l1:l2,m,n,iss) = f(l1:l2,m,n,iss) - reference_state(:,iref_s)
!
      elseif (psize==mx) then
!
!  Reference state not yet considered in this branch as undefined in ghost zones.
!
        if (ldensity_nolog) then
          lnrho_=log(f(:,m,n,irho))
        else
          lnrho_=f(:,m,n,ilnrho)
        endif
        if (present(ee)) then
          yH_=yHmax
          yH_=0.5*min(ee/ee_ion,yH_)
          ivars_mod=ilnrho_ee
          do i=1,nx
            call rtsafe(ilnrho_ee,lnrho_(i),ee(i),yHmin,yHmax*min(ee(i)/ee_ion,1.0),yH_(i))
          enddo
          !call rtsafe_elem(lnrho_,ee,0.*yH_+yHmin,yHmax*min(ee/ee_ion,1.0),yH_)
          fractions=(1+yH_+xHe)
          TT_=(ee-yH_*ee_ion)/(1.5*fractions*ss_ion)
          lnTT_=log(TT_)
          ss_=ss_ion*(fractions*(1.5*(lnTT_-lnTT_ion)-lnrho_+2.5) &
                      -yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                      -(1-yH_)*(log(1-yH_+epsi)-lnrho_H)-xHe_term)
        elseif (present(pp)) then
          yH_=0.5*yHmax
          ivars_mod=ilnrho_pp
          do i=1,nx
            call rtsafe(ilnrho_pp,lnrho_(i),pp(i),yHmin,yHmax,yH_(i))
          enddo
          !call rtsafe_elem(lnrho_,pp,0.*yH_+yHmin,0.*yH_+yHmax,yH_)
          fractions=(1+yH_+xHe)
          rho_=exp(lnrho_)
          TT_=pp/(fractions*ss_ion*rho_)
          lnTT_=log(TT_)
          ss_=ss_ion*(fractions*(1.5*(lnTT_-lnTT_ion)-lnrho_+2.5) &
                     -yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                     -(1-yH_)*(log(1-yH_+epsi)-lnrho_H)-xHe_term)
        elseif (present(ss)) then
          ss_=ss
        endif
        f(:,m,n,iss) = ss_
      endif
!
    endsubroutine eosperturb
!!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,yH,lnTT,ee,pp,cs2,kapparho)
!
!   Calculate thermodynamical quantities
!
!    2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!   10-feb-23/fred: call to pressure_gradient to yield cs2
!
      use Sub
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(out), optional :: lnrho
      real, dimension(psize), intent(out), optional :: yH,lnTT
      real, dimension(psize), intent(out), optional :: ee,pp,kapparho
      real, dimension(psize), optional :: cs2
      real, dimension(psize) :: lnrho_,yH_,lnTT_,TT,fractions,exponent
      integer :: i1,i2
!
      select case (psize)
        case (nx); i1=l1; i2=l2
        case (mx); i1=1; i2=mx
        case default; call fatal_error("eoscalc_farray","no such pencil size")
      end select

      lnrho_=f(i1:i2,m,n,ilnrho)
      yH_   =f(i1:i2,m,n,iyH)
      lnTT_ =f(i1:i2,m,n,ilnTT)

      if (present(lnrho)) lnrho=lnrho_
      if (present(yH)) yH=yH_
      if (present(lnTT)) lnTT=lnTT_
      if (present(cs2)) call pressure_gradient(f,cs2(i1:i2))
!
      if (present(ee).or.present(pp).or.present(kapparho)) TT=exp(lnTT_)
      if (present(ee).or.present(pp)) fractions=1+yH_+xHe
!
      if (present(ee)) ee=1.5*fractions*ss_ion*TT+yH_*ee_ion
      if (present(pp)) pp=fractions*exp(lnrho_)*TT*ss_ion
!
!  Hminus opacity
!
      if (present(kapparho)) then
!        lnchi=2*lnrho_-lnrho_e_+1.5*(lnTT_ion_-lnTT_) &
!             +TT_ion_/TT_+log(yH_+yMetals)+log(1-yH_+epsi)+lnchi0
        exponent = (2*lnrho_-lnrho_e_+1.5*(lnTT_ion_-lnTT_)+TT_ion_/TT)
        !
        ! Ensure exponentiation and successive multiplication with
        ! numbers up to ~ 1e4 avoids overflow:
        !
        exponent = min(exponent,log(huge1)-5.)
        kapparho = exp(exponent + alog(yH_+yMetals))*(1-yH_)*kappa0
      endif
!
    endsubroutine eoscalc_farray
!***********************************************************************
    subroutine eoscalc_point_(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
    
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp,cs2

      select case (ivars)
      case (ipp_ss, irho_eth, ilnrho_eth)
        call not_implemented("eoscalc_point_","thermodynamic variable combinations ipp_ss, irho_eth, ilnrho_eth")
      case default
        call fatal_error("eoscalc_pencil","unknown independent variables combination")
      end select
        
      ivars_mod=ivars

      !call eoscalc_elem(var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)

    endsubroutine eoscalc_point_
!***********************************************************************
    subroutine eoscalc_pencil_(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)

      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), intent(out), optional :: lnrho,ss
      real, dimension(nx), intent(out), optional :: yH,lnTT
      real, dimension(nx), intent(out), optional :: ee,pp,cs2

      select case (ivars)
      case (ipp_ss, irho_eth, ilnrho_eth)
        call not_implemented("eoscalc_pencil_","thermodynamic variable combinations ipp_ss, irho_eth, ilnrho_eth")
      case default
        !!!call fatal_error("eoscalc_pencil","unknown independent variables combination")
      end select

      ivars_mod=ivars

      !call eoscalc_elem(var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)

    endsubroutine eoscalc_pencil_
!***********************************************************************
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!   13-mar-04/tony: modified
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp,cs2
      real :: lnrho_,ss_,yH_,lnTT_,TT_,TT1_,rho_,ee_,pp_
      real :: fractions,rhs,sqrtrhs
!
      ivars_mod=ivars
      select case (ivars)
!
      case (ilnrho_ss,irho_ss)
        if (ivars==ilnrho_ss) then
          lnrho_=var1
          if (present(pp)) rho_=exp(var1)
        else
          lnrho_=alog(var1)
          if (present(pp)) rho_=var1
        endif
        ss_=var2
        yH_=0.5*yHmax
        call rtsafe(ivars,lnrho_,ss_,yHmin,yH_+yHmax,yH_)

        if (present(ee).or.present(lnTT).or.present(pp)) then
          fractions=1+yH_+xHe
          lnTT_=(2.0/3.0)*((ss_/ss_ion+(1-yH_)*(log(1-yH_+epsi)-lnrho_H) &
                            +yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                            +xHe_term)/fractions+lnrho_-2.5)+lnTT_ion
          if (present(ee).or.present(pp)) then
            TT_=exp(lnTT_)
            if (present(ee)) ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion
            if (present(pp)) pp_=fractions*rho_*TT_*ss_ion
          endif
        endif
!
      case (ilnrho_lnTT,ilnrho_TT,irho_TT)
        if (ivars==irho_TT) then
          if (present(pp)) rho_=var1
          lnrho_=log(var1)
        else
          lnrho_=var1
          if (present(pp)) rho_=exp(var1)
        endif
        if (ivars==ilnrho_lnTT) then
          lnTT_=var2
          TT_=exp(lnTT_)
        else
          TT_=var2
          lnTT_=log(TT_)
        endif
        rhs=exp(lnrho_e-lnrho_+1.5*(lnTT_-lnTT_ion)-TT_ion/TT_)
        rhs=max(rhs,tini)       ! avoid log(0.) below
        sqrtrhs=sqrt(rhs)
        yH_=2*sqrtrhs/(sqrtrhs+sqrt(4+rhs))

        if (present(ee).or.present(ss).or.present(pp)) then
          fractions=1+yH_+xHe
          if (present(ss)) ss_= get_ss(lnTT_,lnrho_,yH_,fractions)
          if (present(ee)) ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion
          if (present(pp)) pp_=fractions*rho_*TT_*ss_ion
        endif
!
      case (ilnrho_ee,irho_ee)
        if (ivars==ilnrho_ee) then
          lnrho_=var1
          rho_=exp(lnrho_)
        else
          rho_=var1
          lnrho_=log(rho_)
        endif
        ee_=var2
        yH_=yHmax
        yH_=0.5*min(ee_/ee_ion,yH_)
        call rtsafe(ivars,lnrho_,ee_,0.*yH_+yHmin,yHmax*min(ee_/ee_ion,1.0),yH_)

        if (present(lnTT).or.present(ss).or.present(pp)) then
          fractions=1+yH_+xHe
          TT_=(ee_-yH_*ee_ion)/(1.5*fractions*ss_ion)
          if (present(lnTT).or.present(ss)) then
            lnTT_=log(TT_)
            if (present(ss)) ss_= get_ss(lnTT_,lnrho_,yH_,fractions)
          endif
          if (present(pp)) pp_=fractions*rho_*TT_*ss_ion
        endif
!
      case (ilnrho_pp,irho_pp)
        if (ivars==ilnrho_pp) then
          lnrho_=var1
          rho_=exp(lnrho_)
        else
          rho_=var1
          lnrho_=log(rho_)
        endif
        pp_=var2
        yH_=0.5*yHmax
        call rtsafe(ivars,lnrho_,pp_,0.*yH_+yHmin,0.*yH_+yHmax,yH_)

        if (present(lnTT).or.present(ss).or.present(ee)) then
          fractions=1+yH_+xHe
          TT_=pp_/(fractions*ss_ion*rho_)
          if (present(lnTT).or.present(ss)) then
            lnTT_=log(TT_)
            if (present(ss)) ss_= get_ss(lnTT_,lnrho_,yH_,fractions)
          endif
          if (present(ee)) ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion   !simplify
        endif
!
      case (ipp_ss, irho_eth, ilnrho_eth)
        call not_implemented("eoscalc_point","thermodynamic variable combinations ipp_ss, irho_eth, ilnrho_eth")
      case default
        call fatal_error("eoscalc_pencil","unknown independent variables combination")
      end select
!
      if (present(lnrho)) lnrho=lnrho_
      if (present(ss)) ss=ss_
      if (present(yH)) yH=yH_
      if (present(lnTT)) lnTT=lnTT_
      if (present(ee)) ee=ee_
      if (present(pp)) pp=pp_
      if (present(cs2)) cs2=impossible
!
    endsubroutine eoscalc_point
!***********************************************************************
    function get_ss_pencil(lnTT,lnrho,yH,fractions) result(ss)
!
!  Calculates entropy pencil from lnTT,lnrho,yH,fractions pencils.
!
! 11-nov-2023/MR: coded
!
      real, dimension(nx), intent(IN) :: lnTT, lnrho, yH, fractions
      real, dimension(nx) :: ss

      ss = ss_ion*(fractions*(1.5*(lnTT-lnTT_ion)-lnrho+2.5) &
          -yH*(2*log(yH)-lnrho_e-lnrho_H) &
          -(1-yH)*(log(1-yH+epsi)-lnrho_H)-xHe_term)

    endfunction get_ss_pencil
!***********************************************************************
    elemental function get_ss(lnTT,lnrho,yH,fractions) result(ss)
!
!  Calculates entropy from lnTT,lnrho,yH,fractions at a point
!  Not efficient for pencils.
!
! 11-nov-2023/MR: coded
!
      real, intent(IN) :: lnTT, lnrho, yH, fractions
      real :: ss

      ss = ss_ion*(fractions*(1.5*(lnTT-lnTT_ion)-lnrho+2.5) &
          -yH*(2*log(yH)-lnrho_e-lnrho_H) &
          -(1-yH)*(log(1-yH+epsi)-lnrho_H)-xHe_term)

    endfunction get_ss
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp,cs2)
!
!   Calculate thermodynamical quantities
!
!   i13-mar-04/tony: modified
!
      use Mpicomm, only: stop_it
!
      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), intent(out), optional :: lnrho,ss
      real, dimension(nx), intent(out), optional :: yH,lnTT
      real, dimension(nx), intent(out), optional :: ee,pp,cs2
      real, dimension(nx) :: lnrho_,ss_,yH_,lnTT_,TT_,TT1_,rho_,ee_,pp_
      real, dimension(nx) :: fractions,rhs,sqrtrhs
      integer :: i
!
      ivars_mod=ivars
      select case (ivars)
!
      case (ilnrho_ss,irho_ss)
        if (ivars==ilnrho_ss) then
          lnrho_=var1
        else
          lnrho_=alog(var1)
        endif
        ss_=var2
        yH_=0.5*yHmax
        !call rtsafe_elem(lnrho_,ss_,0.*yH_+yHmin,0.*yH_+yHmax,yH_)
        do i=1,nx
          call rtsafe(ilnrho_ss,lnrho_(i),ss_(i),yHmin,yHmax,yH_(i))
        enddo

        if (present(ee).or.present(lnTT).or.present(pp)) then
          fractions=1+yH_+xHe
          lnTT_=(2.0/3.0)*((ss_/ss_ion+(1-yH_)*(log(1-yH_+epsi)-lnrho_H) &
                            +yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                            +xHe_term)/fractions+lnrho_-2.5)+lnTT_ion
          if (present(ee).or.present(pp)) then
            TT_=exp(lnTT_)
            ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion
            pp_=fractions*exp(lnrho_)*TT_*ss_ion
          endif
        endif
!
      case (ilnrho_lnTT,ilnrho_TT,irho_TT)
        if (ivars==irho_TT) then
          if (present(pp)) rho_=var1
          lnrho_=log(var1)
        else
          lnrho_=var1
          if (present(pp)) rho_=exp(var1)
        endif
        if (ivars==ilnrho_lnTT) then
          lnTT_=var2
          TT_=exp(lnTT_)
          TT1_=exp(-lnTT_)
        else
          TT_=var2
          lnTT_=log(TT_)
          TT1_=1./TT_
        endif
        rhs=exp(lnrho_e-lnrho_+1.5*(lnTT_-lnTT_ion)-TT_ion*TT1_)
        rhs=max(rhs,tini)       ! avoid log(0.) below
        sqrtrhs=sqrt(rhs)
        yH_=2*sqrtrhs/(sqrtrhs+sqrt(4+rhs))

        if (present(ee).or.present(ss).or.present(pp)) then
          fractions=1+yH_+xHe
          if (present(ss)) ss_= get_ss_pencil(lnTT_,lnrho_,yH_,fractions)
          if (present(ee)) ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion
          if (present(pp)) pp_=fractions*rho_*TT_*ss_ion
        endif
!
      case (ilnrho_ee,irho_ee)
        if (ivars==ilnrho_ee) then
          lnrho_=var1
          rho_=exp(lnrho_)
        else
          rho_=var1
          lnrho_=log(rho_)
        endif
        ee_=var2
        yH_=yHmax
        yH_=0.5*min(ee_/ee_ion,yH_)
        !call rtsafe_elem(lnrho_,ee_,0.*yH_+yHmin,yHmax*min(ee_/ee_ion,1.0),yH_)
        do i=1,nx
          call rtsafe(ilnrho_ee,lnrho_(i),ee_(i),yHmin,yHmax,yH_(i))
        enddo

        if (present(lnTT).or.present(ss).or.present(pp)) then
          fractions=1+yH_+xHe
          TT_=(ee_-yH_*ee_ion)/(1.5*fractions*ss_ion)
          if (present(lnTT).or.present(ss)) then
            lnTT_=log(TT_)
            if (present(ss)) ss_= get_ss_pencil(lnTT_,lnrho_,yH_,fractions)
          endif
          if (present(pp)) pp_=fractions*rho_*TT_*ss_ion
        endif
!
      case (ilnrho_pp,irho_pp)
        if (ivars==ilnrho_pp) then
          lnrho_=var1
          rho_=exp(lnrho_)
        else
          rho_=var1
          lnrho_=log(rho_)
        endif
        pp_=var2
        yH_=0.5*yHmax
        !call rtsafe_elem(lnrho_,pp_,0.*yH_+yHmin,0.*yH_+yHmax,yH_)
        do i=1,nx
          call rtsafe(ilnrho_pp,lnrho_(i),pp_(i),yHmin,yHmax,yH_(i))
        enddo

        if (present(lnTT).or.present(ss).or.present(ee)) then
          fractions=1+yH_+xHe
          TT_=pp_/(fractions*ss_ion*rho_)
          if (present(lnTT).or.present(ss)) then
            lnTT_=log(TT_)
            if (present(ss)) ss_= get_ss_pencil(lnTT_,lnrho_,yH_,fractions)
          endif
          if (present(ee)) ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion   !simplify
        endif
!
      case (ipp_ss, irho_eth, ilnrho_eth)
        call not_implemented("eoscalc_point","thermodynamic variable combinations ipp_ss, irho_eth, ilnrho_eth")
      case default
        call fatal_error("eoscalc_pencil","unknown combination of thermodynamic variables")
      end select
!
      if (present(lnrho)) lnrho=lnrho_
      if (present(ss)) ss=ss_
      if (present(yH)) yH=yH_
      if (present(lnTT)) lnTT=lnTT_
      if (present(ee)) ee=ee_
      if (present(pp)) pp=pp_
      if (present(cs2)) cs2=impossible
!
    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine read_eos_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=eos_init_pars, IOSTAT=iostat)
!
    endsubroutine read_eos_init_pars
!***********************************************************************
    subroutine write_eos_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=eos_init_pars)
!
    endsubroutine write_eos_init_pars
!***********************************************************************
    subroutine read_eos_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=eos_run_pars, IOSTAT=iostat)
!
    endsubroutine read_eos_run_pars
!***********************************************************************
    subroutine write_eos_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=eos_run_pars)
!
    endsubroutine write_eos_run_pars
!***********************************************************************
    subroutine rtsafe_pencil(lnrho,ss,yH)
!
!   safe newton raphson algorithm (adapted from NR) !
!   09-apr-03/tobi: changed to subroutine
!
      real, dimension(mx), intent(in) :: lnrho,ss
      real, dimension(mx), intent(inout) :: yH
!
      real, dimension(mx) :: dyHold,dyH,yHlow,yHhigh,ff,dff
      real, dimension(mx) :: lnTT_,dlnTT_,TT1_,fractions1
      logical, dimension(mx) :: found
      integer             :: i
      integer, parameter  :: maxit=1000
!
      yHlow=yHmin
      yHhigh=yHmax
      dyH=yHhigh-yHlow
      dyHold=dyH
!
      found=.false.
!
      fractions1=1/(1+yH+xHe)
      lnTT_=(2.0/3.0)*((ss/ss_ion+(1-yH)*(log(1-yH+epsi)-lnrho_H) &
                         +yH*(2*log(yH)-lnrho_e-lnrho_H) &
                         +xHe_term)*fractions1+lnrho-2.5)
      TT1_=exp(-lnTT_)
      ff=lnrho_e-lnrho+1.5*lnTT_-TT1_+log(1-yH+epsi)-2*log(yH)
      dlnTT_=((2.0/3.0)*(-ff-TT1_)-1)*fractions1
! wd: Need to add epsi at the end, as otherwise (yH-yHlow)*df  will
! wd: eventually yield an overflow.
! wd: Something like  sqrt(tini)  would probably also do instead of epsi,
! wd: but even epsi does not affect the auto-tests so far.
!      df=dlnTT_*(1.5+TT1_)-1/(1-yH+epsi)-2/yH
      dff=dlnTT_*(1.5+TT1_)-1/(1-yH+epsi)-2/(yH+epsi)
!
      do i=1,maxit
        where (.not.found)
          where (      sign(1.,((yH-yHlow)*dff-ff)) &
                    == sign(1.,((yH-yHhigh)*dff-ff)) &
                  .or. abs(2*ff) > abs(dyHold*dff) )
            !
            !  Bisection
            !
            dyHold=dyH
            dyH=0.5*(yHhigh-yHlow)
            yH=yHhigh-dyH
          elsewhere
            !
            !  Newton-Raphson
            !
            dyHold=dyH
            dyH=ff/dff
            ! Apply floor to dyH (necessary to avoid negative yH in samples
            ! /0d-tests/heating_ionize)
            dyH=min(dyH,yH-yHmin)
            dyH=max(dyH,yH-yHmax)    ! plausibly needed as well
            yH=yH-dyH
          endwhere
        endwhere

        where (abs(dyH)>max(yHacc,1e-31)*max(yH,1e-31))     ! use max to avoid underflow
          fractions1=1/(1+yH+xHe)
          lnTT_=(2.0/3.0)*((ss/ss_ion+(1-yH)*(log(1-yH+epsi)-lnrho_H) &
                             +yH*(2*log(yH)-lnrho_e-lnrho_H) &
                             +xHe_term)*fractions1+lnrho-2.5)
          TT1_=exp(-lnTT_)
          ff=lnrho_e-lnrho+1.5*lnTT_-TT1_+log(1-yH+epsi)-2*log(yH)
          dlnTT_=((2.0/3.0)*(-ff-TT1_)-1)*fractions1
          dff=dlnTT_*(1.5+TT1_)-1/(1-yH+epsi)-2/yH
          where (ff<0)
            yHhigh=yH
          elsewhere
            yHlow=yH
          endwhere
        elsewhere
          found=.true.
        endwhere
        if (all(found)) return
      enddo
!
    endsubroutine rtsafe_pencil
!***********************************************************************
    subroutine rtsafe(ivars,var1,var2,yHlb,yHub,yH,rterror,rtdebug)
!
!   safe newton raphson algorithm (adapted from NR) !
!   09-apr-03/tobi: changed to subroutine
!
      integer, intent(in)            :: ivars
      real, intent(in)               :: var1,var2
      real, intent(in)               :: yHlb,yHub
      real, intent(inout)            :: yH
      logical, intent(out), optional :: rterror
      logical, intent(in), optional  :: rtdebug
!
      real               :: dyHold,dyH,yHl,yHh,f,df,temp
      integer            :: i
      integer, parameter :: maxit=1000
!
      if (present(rterror)) rterror=.false.
      if (present(rtdebug)) then
        if (rtdebug) print*,'rtsafe: i,yH=',0,yH
      endif
!
      yHl=yHlb
      yHh=yHub
      dyH=1
      dyHold=dyH
!
      call saha(ivars,var1,var2,yH,f,df)
!
      do i=1,maxit
        if (present(rtdebug)) then
          if (rtdebug) print*,'rtsafe: i,yH=',i,yH
        endif
        if (        sign(1.,((yH-yHl)*df-f)) &
                 == sign(1.,((yH-yHh)*df-f)) &
              .or. abs(2*f) > abs(dyHold*df) ) then
          dyHold=dyH
          dyH=0.5*(yHl-yHh)
          yH=yHh+dyH
          if (yHh==yH) return
        else
          dyHold=dyH
          dyH=f/df
          temp=yH
          yH=yH-dyH
          if (temp==yH) return
        endif
        if (abs(dyH)<yHacc*yH) return
        call saha(ivars,var1,var2,yH,f,df)
        if (f<0) then
          yHh=yH
        else
          yHl=yH
        endif
      enddo
!
      if (present(rterror)) rterror=.true.
!
    endsubroutine rtsafe
!***********************************************************************
    subroutine saha(ivars,var1,var2,yH,f,df)
!
!   we want to find the root of f
!
!   23-feb-03/tobi: errors fixed
!
      integer, intent(in)          :: ivars
      real, intent(in)             :: var1,var2,yH
      real, intent(out)            :: f,df
!
      real :: lnrho,ss,ee,pp
      real :: lnTT_,dlnTT_,TT1_,fractions1
!
      fractions1=1/(1+yH+xHe)
!
      select case (ivars)
      case (ilnrho_ss)
        lnrho=var1
        ss=var2
        lnTT_=(2.0/3.0)*((ss/ss_ion+(1-yH)*(log(1-yH+epsi)-lnrho_H) &
                         +yH*(2*log(yH)-lnrho_e-lnrho_H) &
                         +xHe_term)*fractions1+lnrho-2.5)
      case (ilnrho_ee)
        lnrho=var1
        ee=var2
        !print*,'saha: TT_',2.0/3.0*(ee-yH*ee_ion)*fractions1
        !lnTT_=log(2.0/3.0*(ee/ee_ion-yH)*fractions1)
      !  if (ee<yH*ee_ion) then
      !    lnTT_=-25.
      !  else
        lnTT_=log(2.0/3.0*(ee-yH*ee_ion)*fractions1)
      !  endif
      case (ilnrho_pp)
        lnrho=var1
        pp=var2
        lnTT_=log(pp/ss_ion*fractions1)-lnrho
      case default
        call fatal_error("saha","unknown thermodynamic variable combination")
        lnTT_=0.
      end select
!
      TT1_=exp(-lnTT_)
      f=lnrho_e-lnrho+1.5*lnTT_-TT1_+log(1-yH+epsi)-2*log(yH)
      dlnTT_=((2.0/3.0)*(-f-TT1_)-1)*fractions1
      df=dlnTT_*(1.5+TT1_)-1/(1-yH+epsi)-2/yH
!
    endsubroutine saha
!***********************************************************************
    elemental subroutine rtsafe_elem(var1,var2,yHlb,yHub,yH)
!
!   safe newton raphson algorithm (adapted from NR) !
!   09-apr-03/tobi: changed to subroutine
!
      real, intent(in)    :: var1,var2
      real, intent(in)    :: yHlb,yHub
      real, intent(inout) :: yH
!
      real               :: dyHold,dyH,yHl,yHh,f,df,temp
      integer            :: i
      integer, parameter :: maxit=1000
!
      yHl=yHlb
      yHh=yHub
      dyH=1
      dyHold=dyH
!
      call saha_elem(var1,var2,yH,f,df)
!
      do i=1,maxit
        if (        sign(1.,((yH-yHl)*df-f)) &
                 == sign(1.,((yH-yHh)*df-f)) &
              .or. abs(2*f) > abs(dyHold*df) ) then
          dyHold=dyH
          dyH=0.5*(yHl-yHh)
          yH=yHh+dyH
          if (yHh==yH) return
        else
          dyHold=dyH
          dyH=f/df
          temp=yH
          yH=yH-dyH
          if (temp==yH) return
        endif
        if (abs(dyH)<yHacc*yH) return
        call saha_elem(var1,var2,yH,f,df)
        if (f<0) then
          yHh=yH
        else
          yHl=yH
        endif
      enddo
!
    endsubroutine rtsafe_elem
!***********************************************************************
    elemental subroutine saha_elem(var1,var2,yH,f,df)
!
!   We want to find the root of f.
!
!   23-feb-03/tobi: errors fixed
!
      real, intent(in)  :: var1,var2,yH
      real, intent(out) :: f,df
!
      real :: lnrho,ss,ee,pp
      real :: lnTT,dlnTT,TT1,fractions1
!
      fractions1=1/(1+yH+xHe)
!
      select case (ivars_mod)
      case (ilnrho_ss)
        lnrho=var1
        ss=var2
        lnTT=(2.0/3.0)*((ss/ss_ion+(1-yH)*(log(1-yH+epsi)-lnrho_H) &
                        +yH*(2*log(yH)-lnrho_e-lnrho_H) &
                        +xHe_term)*fractions1+lnrho-2.5)
      case (ilnrho_ee)
        lnrho=var1
        ee=var2
        !print*,'saha: TT',2.0/3.0*(ee-yH*ee_ion)*fractions1
        !lnTT=log(2.0/3.0*(ee/ee_ion-yH)*fractions1)
      !  if (ee<yH*ee_ion) then
      !    lnTT=-25.
      !  else
        lnTT=log(2.0/3.0*(ee-yH*ee_ion)*fractions1)
      !  endif
      case (ilnrho_pp)
        lnrho=var1
        pp=var2
        lnTT=log(pp/ss_ion*fractions1)-lnrho
      end select
!
      TT1=exp(-lnTT)
      f=lnrho_e-lnrho+1.5*lnTT-TT1+log(1-yH+epsi)-2*log(yH)
      dlnTT=((2.0/3.0)*(-f-TT1)-1)*fractions1
      df=dlnTT*(1.5+TT1)-1/(1-yH+epsi)-2/yH
!
    endsubroutine saha_elem
!***********************************************************************
    subroutine isothermal_entropy(lnrho_arr,T0,ss_arr)
!
!  Isothermal stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!  Sound speed (and hence Temperature), is
!  initialised to the reference value:
!           sound speed: cs^2_0            from start.in
!           density: rho0 = exp(lnrho0)
!
!  11-jun-03/tony: extracted from isothermal routine in Density module
!                  to allow isothermal condition for arbitrary density
!  17-oct-03/nils: works also with leos_ionization=T
!  18-oct-03/tobi: distributed across ionization modules
!
      real, intent(in) :: T0
      real, dimension(mx,my,mz), intent(in) :: lnrho_arr
      real, dimension(mx,my,mz), intent(out):: ss_arr

      real, dimension(nx) :: lnrho,yH,K,sqrtK,yH_term,one_yH_term
!
      do n=n1,n2
      do m=m1,m2

        lnrho=lnrho_arr(l1:l2,m,n)
        K=exp(lnrho_e-lnrho-TT_ion/T0)*(T0/TT_ion)**1.5
        sqrtK=sqrt(K)
        yH=2*sqrtK/(sqrtK+sqrt(4+K))
!
        where (yH>0)
          yH_term=yH*(2*log(yH)-lnrho_e-lnrho_H)
        elsewhere
          yH_term=0
        endwhere
!
        where (yH<1)
          one_yH_term=(1-yH)*(log(1-yH+epsi)-lnrho_H)
        elsewhere
          one_yH_term=0
        endwhere
!
        ss_arr(l1:l2,m,n)=ss_ion*((1+yH+xHe)*(1.5*log(T0/TT_ion)-lnrho+2.5)-yH_term-one_yH_term-xHe_term)
!
      enddo
      enddo
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine bc_ss_flux(f,topbot,lone_sided)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!   3-oct-16/MR: added new optional switch lone_sided
!
      use Gravity
      use SharedVariables,only: get_shared_variable
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      logical, pointer :: lmultilayer, lheatc_chiconst
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      logical, optional :: lone_sided
      real, dimension (size(f,1),size(f,2)) :: tmp_xy,TT_xy,rho_xy,yH_xy
      integer :: i
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
!  Get the shared variables
!
      call get_shared_variable('hcond0',hcond0,caller='bc_ss_flux')
      call get_shared_variable('hcond1',hcond1)
      call get_shared_variable('Fbot',Fbot)
      call get_shared_variable('Ftop',Ftop)
      call get_shared_variable('FbotKbot',FbotKbot)
      call get_shared_variable('FtopKtop',FtopKtop)
      call get_shared_variable('chi',chi)
      call get_shared_variable('lmultilayer',lmultilayer)
      call get_shared_variable('lheatc_chiconst',lheatc_chiconst)
!
      select case (topbot)
!
!  bottom boundary
!  ---------------
!
      case(BOT)
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Fbot,hcond=',Fbot,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n1,ilnrho))
        TT_xy=exp(f(:,:,n1,ilnTT))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!
        if (lheatc_chiconst) then
          tmp_xy=Fbot/(rho_xy*chi*TT_xy)
        else
          tmp_xy=FbotKbot/TT_xy
        endif
!
!  get ionization fraction at bottom boundary
!
        yH_xy=f(:,:,n1,iyH)
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+ss_ion*(1+yH_xy+xHe)* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+3*i*dz*tmp_xy)
        enddo
!
!  top boundary
!  ------------
!
      case(TOP)
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0
        endif
!
!  calculate Fbot/(K*cs2)
!
        rho_xy=exp(f(:,:,n2,ilnrho))
        TT_xy=exp(f(:,:,n2,ilnTT))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!
        if (lheatc_chiconst) then
          tmp_xy=Ftop/(rho_xy*chi*TT_xy)
        else
          tmp_xy=FtopKtop/TT_xy
        endif
!
!  get ionization fraction at top boundary
!
        yH_xy=f(:,:,n2,iyH)
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n2+i,iss)= f(:,:,n2-i,iss)+ss_ion*(1+yH_xy+xHe)* &
                          (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)-3*i*dz*tmp_xy)
        enddo
!
      case default
        call fatal_error("bc_ss_flux","topbot should be BOT or TOP")
      endselect
      call keep_compiler_quiet(present(lone_sided))
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_flux_turb(f,topbot)
!
!  dummy routine
!
!   4-may-2009/axel: dummy routine
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_flux_turb","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb
!***********************************************************************
    subroutine bc_ss_flux_turb_x(f,topbot)
!
!  dummy routine
!
!   31-may-2010/pete: dummy routine
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_flux_turb_x","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_x(f,topbot)
!
!   23-apr-2014/pete: dummy
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_flux_condturb_x","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_mean_x(f,topbot)
!
!   07-jan-2015/pete: dummy
!
      integer, intent(IN) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call not_implemented("bc_ss_flux_condturb_z_mean_x","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_mean_x
!***********************************************************************
    subroutine bc_ss_flux_condturb_z(f,topbot)
!
!   15-jul-2014/pete: dummy
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_flux_condturb_z","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_condturb_z
!***********************************************************************
    subroutine bc_ss_temp_old(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  23-jun-2003/tony: implemented for leos_fixed_ionization
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_temp_old","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_old
!***********************************************************************
    subroutine bc_ss_temp_x(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_temp_x","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_x
!***********************************************************************
    subroutine bc_ss_temp_y(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_temp_y","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot,lone_sided)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      logical, optional :: lone_sided
!
      call not_implemented("bc_ss_temp_z","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
      call keep_compiler_quiet(present(lone_sided))
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_lnrho_temp_z(f,topbot)
!
!  boundary condition for density: constant temperature
!
!  19-aug-2005/tobi: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_lnrho_temp_z","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_temp_z
!***********************************************************************
    subroutine bc_lnrho_pressure_z(f,topbot)
!
!  boundary condition for density: constant pressure
!
!  19-aug-2005/tobi: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_lnrho_pressure_z","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_pressure_z
!***********************************************************************
    subroutine bc_ss_temp2_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_temp2_z","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_temp3_z(f,topbot)
!
!  31-jan-2013/axel: coded to impose cs2bot and dcs2bot at bottom
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented('bc_ss_temp3_z','in eos_ionization')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp3_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_stemp_x","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_x
!***********************************************************************
    subroutine bc_ss_stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_stemp_y","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!  26-sep-2003/tony: coded
!
      use Gravity
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,1),size(f,2),nghost) :: lnrho,ss,yH,lnTT,TT,K,sqrtK,yH_term,one_yH_term
      integer :: i
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case(BOT)
        do i=1,nghost
          f(:,:,n1-i,ilnTT) = f(:,:,n1+i,ilnTT)
        enddo
!
        lnrho=f(:,:,1:n1-1,ilnrho)
        lnTT=f(:,:,1:n1-1,ilnTT)
        TT=exp(lnTT)
!
        K=exp(lnrho_e-lnrho-TT_ion/TT)*(TT/TT_ion)**1.5
        sqrtK=sqrt(K)
        yH=2*sqrtK/(sqrtK+sqrt(4+K))
!
        where (yH>0)
          yH_term=yH*(2*log(yH)-lnrho_e-lnrho_H)
        elsewhere
          yH_term=0
        endwhere
!
        where (yH<1)
          one_yH_term=(1-yH)*(log(1-yH+epsi)-lnrho_H)
        elsewhere
          one_yH_term=0
        endwhere
!
        ss=ss_ion*((1+yH+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5) &
                    -yH_term-one_yH_term-xHe_term)
!
        f(:,:,1:n1-1,iyH)=yH
        f(:,:,1:n1-1,iss)=ss
!
!  top boundary
!
      case(TOP)
        do i=1,nghost
          f(:,:,n2+i,ilnTT) = f(:,:,n2-i,ilnTT)
        enddo
!
        lnrho=f(:,:,n2+1:,ilnrho)
        lnTT =f(:,:,n2+1:,ilnTT)
        TT=exp(lnTT)
!
        K=exp(lnrho_e-lnrho-TT_ion/TT)*(TT/TT_ion)**1.5
        sqrtK=sqrt(K)
        yH=2*sqrtK/(sqrtK+sqrt(4+K))
!
        where (yH>0)
          yH_term=yH*(2*log(yH)-lnrho_e-lnrho_H)
        elsewhere
          yH_term=0
        endwhere
!
        where (yH<1)
          one_yH_term=(1-yH)*(log(1-yH+epsi)-lnrho_H)
        elsewhere
          one_yH_term=0
        endwhere
!
        ss=ss_ion*((1+yH+xHe)*(1.5*log(TT/TT_ion)-lnrho+2.5) &
                    -yH_term-one_yH_term-xHe_term)
!
        f(:,:,n2+1:,iyH)=yH
        f(:,:,n2+1:,iss)=ss
!
      case default
        call fatal_error('bc_ss_stemp_z','topbot should be BOT or TOP')
      endselect
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_a2stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_a2stemp_x","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_x
!***********************************************************************
    subroutine bc_ss_a2stemp_y(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_a2stemp_y","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_y
!***********************************************************************
    subroutine bc_ss_a2stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_ss_a2stemp_z","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_a2stemp_z
!***********************************************************************
    subroutine bc_ss_energy(f,topbot)
!
!  boundary condition for entropy
!
!  may-2002/nils: coded
!  11-jul-2002/nils: moved into the entropy module
!  26-aug-2003/tony: distributed across ionization modules
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
!  The 'ce' boundary condition for entropy makes the energy constant at
!  the boundaries.
!  This assumes that the density is already set (ie density must register
!  first!)
!
      call not_implemented("bc_ss_stemp_y","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_stellar_surface","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_lnrho_cfb_r_iso","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_lnrho_hds_z_iso","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hds_z_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call not_implemented("bc_lnrho_hdss_z_iso","in eos_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine bc_ism(f,topbot,j)
!
!  30-nov-15/fred: Replaced bc_ctz and bc_cdz.
!  Apply observed scale height locally from Reynolds 1991, Manchester & Taylor
!  1981 for warm ionized gas - dominant scale height above 500 parsecs.
!  Apply constant local temperature across boundary for entropy.
!  Motivation to prevent numerical spikes in shock fronts, which cannot be
!  absorbed in only three ghost cells, but boundary thermodynamics still
!  responsive to interior dynamics.
!  06-jun-22/fred update to allow setting scale height in start.in or run.in
!  default is density_scale_factor=impossible so that scale_factor is 0.9, assuming
!  unit_length = 1 kpc and scale is 400 pc. To change scale height add to
!  start_pars or run_pars density_scale_factor=... in dimensionless units
!
      integer, intent(IN) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: j,k,stat
      real :: density_scale1, density_scale
      real, dimension (:,:), allocatable :: cv,cp
!
      if (density_scale_factor==impossible) then
        density_scale=density_scale_cgs/unit_length
      else
        density_scale=density_scale_factor
      endif
      density_scale1=1./density_scale
!
      if (j==iss.and..not.ltemperature) then
        allocate(cp(size(f,1),size(f,2)),stat=stat)
        if (stat>0) call fatal_error('bc_ism','could not allocate cp')
        allocate(cv(size(f,1),size(f,2)),stat=stat)
        if (stat>0) call fatal_error('bc_ism','could not allocate cv')
      endif
!
      select case (topbot)
!
      case(BOT)               ! bottom boundary
        do k=1,nghost
          if (j==irho .or. j==ilnrho) then
            if (ldensity_nolog) then
              f(:,:,k,j)=f(:,:,n1,j)*exp(-(z(n1)-z(k))*density_scale1)
            else
              f(:,:,k,j)=f(:,:,n1,j) - (z(n1)-z(k))*density_scale1
            endif
          else if (j==iss) then
            if (.not.ltemperature) then !case for entropy
              cp=(2.5+f(:,:,n1,iyH)*(1-f(:,:,n1,iyH))/((2-f(:,:,n1,iyH))*xHe+2)* &
                  (2.5+TT_ion/exp(f(:,:,n1,ilnTT)))**2)* &
                  Rgas*mu1yHxHe/(1+xHe+f(:,:,n1,iyH))
              cv=(1.5+f(:,:,n1,iyH)*(1-f(:,:,n1,iyH))/((2-f(:,:,n1,iyH))* &
                  (1+f(:,:,n1,iyH)+xHe))*(1.5+TT_ion/exp(f(:,:,n1,ilnTT)))**2)* &
                  Rgas*mu1yHxHe/(1+xHe+f(:,:,n1,iyH))
              if (ldensity_nolog) then
                f(:,:,n1-k,j)=f(:,:,n1,j)+(cp-cv) * &
                    (log(f(:,:,n1,j-1))-log(f(:,:,n1-k,j-1))) + &
                    cv*log((z(n1)-z(n1-k))*density_scale+1.)
              else
                f(:,:,n1-k,j)=f(:,:,n1,j)+(cp-cv)*&
                    (f(:,:,n1,j-1)-f(:,:,n1-k,j-1))+&
                    cv*log((z(n1)-z(n1-k))*density_scale+1.)
              endif
            else !case for lnTT
              f(:,:,n1-k,j)=f(:,:,n1,j)+log((z(n1)-z(n1-k))*density_scale+1.)
            endif
          else
            call fatal_error('bc_ism','only for irho, ilnrho, iuz or iss')
          endif
        enddo
!
      case(TOP)               ! top boundary
        do k=1,nghost
          if (j==irho .or. j==ilnrho) then
            if (ldensity_nolog) then
              f(:,:,n2+k,j)=f(:,:,n2,j)*exp(-(z(n2+k)-z(n2))*density_scale1)
            else
              f(:,:,n2+k,j)=f(:,:,n2,j) - (z(n2+k)-z(n2))*density_scale1
            endif
          else if (j==iss) then
            if (.not.ltemperature) then !case for entropy
              cp=(2.5+f(:,:,n2,iyH)*(1-f(:,:,n2,iyH))/((2-f(:,:,n2,iyH))*xHe+2)* &
                  (2.5+TT_ion/exp(f(:,:,n2,ilnTT)))**2)* &
                  Rgas*mu1yHxHe/(1+xHe+f(:,:,n2,iyH))
              cv=(1.5+f(:,:,n2,iyH)*(1-f(:,:,n2,iyH))/((2-f(:,:,n2,iyH))* &
                  (1+f(:,:,n2,iyH)+xHe))*(1.5+TT_ion/exp(f(:,:,n2,ilnTT)))**2)* &
                  Rgas*mu1yHxHe/(1+xHe+f(:,:,n2,iyH))
              if (ldensity_nolog) then
                f(:,:,n2+k,j)=f(:,:,n2,j)+(cp-cv)*&
                    (log(f(:,:,n2,j-1))-log(f(:,:,n2+k,j-1)))+&
                    cv*log((z(n2+k)-z(n2))*density_scale+1.)
              else
                f(:,:,n2+k,j)=f(:,:,n2,j)+(cp-cv)*&
                    (f(:,:,n2,j-1)-f(:,:,n2+k,j-1))+&
                    cv*log((z(n2+k)-z(n2))*density_scale+1.)
              endif
            else !case for lnTT
              f(:,:,n1-k,j)=f(:,:,n1,j)+log((z(n1)-z(n1-k))*density_scale+1.)
            endif
          else
            call fatal_error('bc_ism','only for irho, ilnrho, iuz or iss')
          endif
        enddo
!
      case default
        call fatal_error("bc_ism","topbot should be BOT or TOP")
!
      endselect
!
    endsubroutine bc_ism
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(cs20,p_par(1))

    endsubroutine pushpars2c
!***********************************************************************
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include 'eos_dummies.inc'
!***********************************************************************
endmodule EquationOfState
