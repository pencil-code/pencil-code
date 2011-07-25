! $Id$
!
!  This modules contains the routines for simulation with
!  simple hydrogen ionization.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 2
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; cp; cp1; cp1tilde
! PENCILS PROVIDED glnTT(3); TT; TT1; gTT(3); yH; hss(3,3); hlnTT(3,3)
! PENCILS PROVIDED del2ss; del6ss; del2lnTT; cv1; glnmumol(3); ppvap; csvap2
! PENCILS PROVIDED rho_anel
!
!***************************************************************
module EquationOfState
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'eos.h'
! integers specifying which independent variables to use in eoscalc
  integer, parameter :: ilnrho_ss=1,ilnrho_ee=2,ilnrho_pp=3,ilnrho_lnTT=4
  integer, parameter :: irho_ss=7, ilnrho_TT=9, irho_TT=10, ipp_ss=11
  integer, parameter :: ipp_cs2=12
!
  interface eoscalc              ! Overload subroutine eoscalc
    module procedure eoscalc_farray
    module procedure eoscalc_pencil
    module procedure eoscalc_point
  endinterface
!
  interface pressure_gradient    ! Overload subroutine pressure_gradient
    module procedure pressure_gradient_farray
    module procedure pressure_gradient_point
  endinterface
!  secondary parameters calculated in initialize
  real :: TT_ion,lnTT_ion,TT_ion_,lnTT_ion_
  real :: ss_ion,ee_ion,kappa0,xHe_term,ss_ion1,Srad0
  real :: lnrho_e,lnrho_e_,lnrho_H,lnrho_He
  integer :: l
! namelist parameters
  real, parameter :: yHmin=tiny(TT_ion), yHmax=1-epsilon(TT_ion)
  real :: xHe=0.1
  real :: yMetals=0
  real :: yHacc=1e-5
! input parameters
  namelist /eos_init_pars/ xHe,yMetals,yHacc
! run parameters
  namelist /eos_run_pars/ xHe,yMetals,yHacc
!ajwm  Moved here from Density.f90
!ajwm  Completely irrelevant to eos_ionization but density and entropy need
!ajwm  reworking to be independent of these things first
!ajwm  can't use impossible else it breaks reading param.nml
!ajwm  SHOULDN'T BE HERE... But can wait till fully unwrapped
  real :: cs0=impossible, rho0=impossible, cp=impossible
  real :: cs20=impossible, lnrho0=impossible
  logical :: lcalc_cp = .false.
  real :: gamma=impossible, gamma_m1=impossible,gamma_inv=impossible
  real :: cs2top_ini=impossible, dcs2top_ini=impossible
  real :: cs2bot=impossible, cs2top=impossible
!ajwm  Not sure this should exist either...
  real :: cs2cool=0.
  real :: mpoly=1.5, mpoly0=1.5, mpoly1=1.5, mpoly2=1.5
  integer :: isothtop=0
  real, dimension (3) :: beta_glnrho_global=0.0,beta_glnrho_scaled=0.0
!
  character (len=labellen) :: ieos_profile='nothing'
  real, dimension(mz) :: profz_eos=1.,dprofz_eos=0.
!
  real, dimension(nchemspec,18) :: species_constants
  real, dimension(nchemspec,7)     :: tran_data
  real, dimension(nchemspec)  :: Lewis_coef, Lewis_coef1
!
  contains
!***********************************************************************
    subroutine register_eos()
!
!   2-feb-03/axel: adapted from Interstellar module
!   13-jun-03/tobi: re-adapted from visc_shock module
!
      use FArrayManager
      use Sub
!
      leos_ionization=.true.
!
!  Set indices for auxiliary variables.
!
      call farray_register_auxiliary('yH',iyH)
      call farray_register_auxiliary('lnTT',ilnTT)
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
!
    endsubroutine register_eos
!***********************************************************************
    subroutine units_eos()
!
!  If unit_temperature hasn't been specified explictly in start.in,
!  set it to 1 (Kelvin).
!
!  24-jun-06/tobi: coded
!
      if (unit_temperature==impossible) unit_temperature=1.
!
    endsubroutine units_eos
!***********************************************************************
    subroutine initialize_eos()
!
!  Perform any post-parameter-read initialization, e.g. set derived
!  parameters.
!
!   2-feb-03/axel: adapted from Interstellar module
!
      use General
      use Mpicomm, only: stop_it
!
      real :: mu1yHxHe
!
      if (lroot) print*,'initialize_eos: ENTER'
!
!  ionization parameters
!  since m_e and chiH, as well as hbar are all very small
!  it is better to divide m_e and chiH separately by hbar.
!
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
        call stop_it('initialize_eos: xHe lower than zero makes no sense')
      else
        xHe_term=0
      endif
!
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
!  Writes iyH and ilnTT to index.pro file
!
!  14-jun-03/axel: adapted from rprint_radiation
!  21-11-04/anders: moved diagnostics to entropy
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_eos
!***********************************************************************
    subroutine get_slices_eos(f,slices)
!
!  Write slices for animation of Eos variables.
!
!  26-jul-06/tony: coded
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
        case ('lnTT')
          slices%yz =f(ix_loc,m1:m2,n1:n2,ilnTT)
          slices%xz =f(l1:l2,iy_loc,n1:n2,ilnTT)
          slices%xy =f(l1:l2,m1:m2,iz_loc,ilnTT)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,ilnTT)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,ilnTT)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,ilnTT)
          slices%ready=.true.
!
!  Degree of ionization.
!
        case ('yH')
          slices%yz =f(ix_loc,m1:m2,n1:n2,iyH)
          slices%xz =f(l1:l2,iy_loc,n1:n2,iyH)
          slices%xy =f(l1:l2,m1:m2,iz_loc,iyH)
          slices%xy2=f(l1:l2,m1:m2,iz2_loc,iyH)
          if (lwrite_slice_xy3) slices%xy3=f(l1:l2,m1:m2,iz3_loc,iyH)
          if (lwrite_slice_xy4) slices%xy4=f(l1:l2,m1:m2,iz4_loc,iyH)
          slices%ready=.true.
!
      endselect
!
    endsubroutine get_slices_eos
!***********************************************************************
    subroutine pencil_criteria_eos()
!
!  All pencils that the EquationOfState module depends on are specified here.
!
!  02-04-06/tony: coded
!
!  EOS is a pencil provider but evolves nothing so it is unlokely that
!  it will require any pencils for it's own use.
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
      if (lpencil_in(i_gTT)) then
        lpencil_in(i_glnTT)=.true.
        lpencil_in(i_TT)=.true.
      endif
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
!
      if (lpencil_in(i_hlnTT)) then
        lpencil_in(i_hss)=.true.
        if (.not.pretend_lnTT) lpencil_in(i_hlnrho)=.true.
      endif
!
    endsubroutine pencil_interdep_eos
!***********************************************************************
    subroutine calc_pencils_eos(f,p)
!
!  Calculate Entropy pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  02-04-06/tony: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      integer :: i
!
      intent(in) :: f
      intent(inout) :: p
!
! THE FOLLOWING 2 ARE CONCEPTUALLY WRONG
! FOR pretend_lnTT since iss actually contain lnTT NOT entropy!
! The code is not wrong however since this is correctly
! handled by the eos module.
! ss
      if (lpencil(i_ss)) p%ss=f(l1:l2,m,n,iss)
!
! gss
      if (lpencil(i_gss)) call grad(f,iss,p%gss)
! pp
      if (lpencil(i_pp)) call eoscalc(f,nx,pp=p%pp)
! ee
      if (lpencil(i_ee)) call eoscalc(f,nx,ee=p%ee)
! lnTT
      if (lpencil(i_lnTT)) call eoscalc(f,nx,lnTT=p%lnTT)
! yH
      if (lpencil(i_yH)) call eoscalc(f,nx,yH=p%yH)
! TT
      if (lpencil(i_TT)) p%TT=exp(p%lnTT)
! TT1
      if (lpencil(i_TT1)) p%TT1=exp(-p%lnTT)
! cs2 and cp1tilde
      if (lpencil(i_cs2) .or. lpencil(i_cp1tilde)) &
          call pressure_gradient(f,p%cs2,p%cp1tilde)
! glnTT
      if (lpencil(i_glnTT)) then
        call temperature_gradient(f,p%glnrho,p%gss,p%glnTT)
      endif
! gTT
      if (lpencil(i_gTT)) then
        do i=1,3; p%gTT(:,i)=p%glnTT(:,i)*p%TT; enddo
      endif
! hss
      if (lpencil(i_hss)) then
        call g2ij(f,iss,p%hss)
      endif
! del2ss
      if (lpencil(i_del2ss)) then
        call del2(f,iss,p%del2ss)
      endif
! del2lnTT
      if (lpencil(i_del2lnTT)) then
          call temperature_laplacian(f,p)
      endif
! del6ss
      if (lpencil(i_del6ss)) then
        call del6(f,iss,p%del6ss)
      endif
! hlnTT
      if (lpencil(i_hlnTT)) then
        call temperature_hessian(f,p%hlnrho,p%hss,p%hlnTT)
      endif
!
      if (lpencil(i_glnmumol)) p%glnmumol(:,:)=0.
!
    endsubroutine calc_pencils_eos
!***********************************************************************
    subroutine getmu(f,mu)
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
!
      real, dimension (mx,my,mz,mfarray), optional :: f
      real, intent(out) :: mu
!
      mu=1.+3.97153*xHe
!
! tobi: the real mean molecular weight would be:
!
! mu=(1.+3.97153*xHe)/(1+yH+xHe)
!
      call keep_compiler_quiet(present(f))
!
    endsubroutine getmu
!***********************************************************************
    subroutine getmu_array(f,mu1_full_tmp)
!
!  dummy routine to calculate mean molecular weight
!
!   16-mar-10/natalia
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: mu1_full_tmp
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(mu1_full_tmp)
!
    endsubroutine getmu_array
!***********************************************************************
    subroutine ioninit(f)
!
!  the ionization fraction has to be set to a value yH0 < yH < yHmax before
!  rtsafe is called for the first time
!
!  12-jul-03/tobi: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      f(:,:,:,iyH) = 0.5*(yHmax-yHmin)
      call ioncalc(f)
!
    endsubroutine ioninit
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
  subroutine gettemperature(f,TT_tmp)
!
     real, dimension (mx,my,mz,mfarray) :: f
     real, dimension (mx,my,mz), intent(out) :: TT_tmp
!
     call keep_compiler_quiet(f)
     call keep_compiler_quiet(TT_tmp)
!
   endsubroutine gettemperature
!***********************************************************************
 subroutine getpressure(pp_tmp)
!
     real, dimension (mx,my,mz), intent(out) :: pp_tmp
!
     call keep_compiler_quiet(pp_tmp)
!
   endsubroutine getpressure
!***********************************************************************
    subroutine get_cp1(cp1_)
!
!  04-nov-06/axel: added to alleviate spurious use of pressure_gradient
!
!  return the value of cp1 to outside modules
!
      real, intent(out) :: cp1_
!
!  for variable ionization, it doesn't make sense to calculate
!  just a single value of cp1, because it must depend on position.
!  Therefore, return impossible, so one can reconsider this case.
!
      call fatal_error('get_cp1','SHOULD NOT BE CALLED WITH eos_ionization')
      cp1_=impossible
!
    endsubroutine get_cp1
!***********************************************************************
    subroutine get_cv1(cv1_)
!
!  22-dec-10/PJK: adapted from get_cp1
!
!  return the value of cv1 to outside modules
!
      real, intent(out) :: cv1_
!
!  for variable ionization, it doesn't make sense to calculate
!  just a single value of cv1, because it must depend on position.
!  Therefore, return impossible, so one can reconsider this case.
!
      call fatal_error('get_cv1','SHOULD NOT BE CALLED WITH eos_ionization')
      cv1_=impossible
!
    endsubroutine get_cv1
!***********************************************************************
    subroutine pressure_gradient_farray(f,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   17-nov-03/tobi: adapted from subroutine eoscalc
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: cs2,cp1tilde
      real, dimension(nx) :: lnrho,yH,lnTT
      real, dimension(nx) :: R,dlnTTdy,dRdy,temp
      real, dimension(nx) :: dlnPPdlnrho,fractions,fractions1
      real, dimension(nx) :: dlnPPdss,TT1
!
      lnrho=f(l1:l2,m,n,ilnrho)
      yH=f(l1:l2,m,n,iyH)
      lnTT=f(l1:l2,m,n,ilnTT)
      TT1=exp(-lnTT)
      fractions=(1+yH+xHe)
      fractions1=1/fractions
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
    subroutine temperature_laplacian(f,p)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   12-dec-05/tony: adapted from subroutine temperature_gradient
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      type (pencil_case) :: p
!
      call not_implemented('temperature_laplacian')
!
      p%del2lnTT=0.0
      call keep_compiler_quiet(f)
!
    endsubroutine temperature_laplacian
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
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(in), optional :: ee,pp,ss
!
      call not_implemented("eosperturb")
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(present(ee),present(pp),present(ss))
!
    endsubroutine eosperturb
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,ss,yH,lnTT,ee,pp,kapparho)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!
      use Sub
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(out), optional :: lnrho,ss
      real, dimension(psize), intent(out), optional :: yH,lnTT
      real, dimension(psize), intent(out), optional :: ee,pp,kapparho
      real, dimension(psize) :: lnrho_,ss_,yH_,lnTT_,TT_,fractions,exponent
!
!
      select case (psize)
!
      case (nx)
        lnrho_=f(l1:l2,m,n,ilnrho)
        ss_=f(l1:l2,m,n,iss)
        yH_=f(l1:l2,m,n,iyH)
        lnTT_=f(l1:l2,m,n,ilnTT)
!
      case (mx)
        lnrho_=f(:,m,n,ilnrho)
        ss_=f(:,m,n,iss)
        yH_=f(:,m,n,iyH)
        lnTT_=f(:,m,n,ilnTT)
!
      case default
        call stop_it("eoscalc: no such pencil size")
!
      end select
!
      TT_=exp(lnTT_)
      fractions=(1+yH_+xHe)
!
      if (present(lnrho)) lnrho=lnrho_
      if (present(ss)) ss=ss_
      if (present(yH)) yH=yH_
      if (present(lnTT)) lnTT=lnTT_
      if (present(ee)) ee=1.5*fractions*ss_ion*TT_+yH_*ee_ion
      if (present(pp)) pp=fractions*exp(lnrho_)*TT_*ss_ion
!
!  Hminus opacity
!
      if (present(kapparho)) then
!        lnchi=2*lnrho_-lnrho_e_+1.5*(lnTT_ion_-lnTT_) &
!             +TT_ion_/TT_+log(yH_+yMetals)+log(1-yH_+epsi)+lnchi0
        exponent = (2*lnrho_-lnrho_e_+1.5*(lnTT_ion_-lnTT_)+TT_ion_/TT_)
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
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
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
      real, dimension(nx), intent(out), optional :: ee,pp
      real, dimension(nx) :: lnrho_,ss_,yH_,lnTT_,TT_,TT1_,rho_,ee_,pp_
      real, dimension(nx) :: fractions,rhs,sqrtrhs
      integer :: i
!
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
        do i=1,nx
          call rtsafe(ilnrho_ss,lnrho_(i),ss_(i),yHmin,yHmax,yH_(i))
        enddo
        fractions=(1+yH_+xHe)
        lnTT_=(2.0/3.0)*((ss_/ss_ion+(1-yH_)*(log(1-yH_+epsi)-lnrho_H) &
                          +yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                          +xHe_term)/fractions+lnrho_-2.5)+lnTT_ion
        TT_=exp(lnTT_)
        rho_=exp(lnrho_)
        ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion
        pp_=fractions*rho_*TT_*ss_ion
!
      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
        TT1_=exp(-lnTT_)
        rhs=exp(lnrho_e-lnrho_+1.5*(lnTT_-lnTT_ion)-TT_ion*TT1_)
        rhs=max(rhs,tini)       ! avoid log(0.) below
        sqrtrhs=sqrt(rhs)
        yH_=2*sqrtrhs/(sqrtrhs+sqrt(4+rhs))
        fractions=(1+yH_+xHe)
        ss_=ss_ion*(fractions*(1.5*(lnTT_-lnTT_ion)-lnrho_+2.5) &
                   -yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                   -(1-yH_)*(log(1-yH_+epsi)-lnrho_H)-xHe_term)
        ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion
        pp_=(1+yH_+xHe)*exp(lnrho_)*TT_*ss_ion
!
      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
        yH_=yHmax
        yH_=0.5*min(ee_/ee_ion,yH_)
        do i=1,nx
          call rtsafe(ilnrho_ee,lnrho_(i),ee_(i),yHmin,yHmax*min(ee_(i)/ee_ion,1.0),yH_(i))
        enddo
        fractions=(1+yH_+xHe)
        TT_=(ee_-yH_*ee_ion)/(1.5*fractions*ss_ion)
        lnTT_=log(TT_)
        rho_=exp(lnrho_)
        ss_=ss_ion*(fractions*(1.5*(lnTT_-lnTT_ion)-lnrho_+2.5) &
                    -yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                    -(1-yH_)*(log(1-yH_+epsi)-lnrho_H)-xHe_term)
        pp_=fractions*rho_*TT_*ss_ion
!
      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
        yH_=0.5*yHmax
        do i=1,nx
          call rtsafe(ilnrho_pp,lnrho_(i),pp_(i),yHmin,yHmax,yH_(i))
        enddo
        fractions=(1+yH_+xHe)
        rho_=exp(lnrho_)
        TT_=pp_/(fractions*ss_ion*rho_)
        lnTT_=log(TT_)
        ss_=ss_ion*(fractions*(1.5*(lnTT_-lnTT_ion)-lnrho_+2.5) &
                   -yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                   -(1-yH_)*(log(1-yH_+epsi)-lnrho_H)-xHe_term)
        ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion
!
      case default
        call stop_it("eoscalc_pencil: I don't get what the independent variables are.")
!
      end select
!
      if (present(lnrho)) lnrho=lnrho_
      if (present(ss)) ss=ss_
      if (present(yH)) yH=yH_
      if (present(lnTT)) lnTT=lnTT_
      if (present(ee)) ee=ee_
      if (present(pp)) pp=pp_
!
    endsubroutine eoscalc_pencil
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
!
      use Mpicomm, only: stop_it
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp,cs2
      real :: lnrho_,ss_,yH_,lnTT_,TT_,TT1_,rho_,ee_,pp_,cs2_
      real :: fractions,rhs,sqrtrhs
!
      select case (ivars)
!
      case (ilnrho_ss)
        lnrho_=var1
        ss_=var2
        yH_=0.5*yHmax
        call rtsafe(ilnrho_ss,lnrho_,ss_,yHmin,yHmax,yH_)
        lnTT_=(ss_/ss_ion+(1-yH_)*(log(1-yH_+epsi)-lnrho_H) &
              +yH_*(2*log(yH_)-lnrho_e-lnrho_H)+xHe_term)/(1+yH_+xHe)
        lnTT_=(2.0/3.0)*(lnTT_+lnrho_-2.5)+lnTT_ion
!
        TT_=exp(lnTT_)
        rho_=exp(lnrho_)
        ee_=1.5*(1+yH_+xHe)*ss_ion*TT_+yH_*ee_ion
        pp_=(1+yH_+xHe)*rho_*TT_*ss_ion
        cs2_=impossible
!
      case (ilnrho_lnTT)
        lnrho_=var1
        lnTT_=var2
        TT_=exp(lnTT_)
        TT1_=1/TT_
        rhs=exp(lnrho_e-lnrho_+1.5*(lnTT_-lnTT_ion)-TT_ion*TT1_)
        rhs = max(rhs,tini)     ! avoid log(0.) below
        sqrtrhs=sqrt(rhs)
        yH_=2*sqrtrhs/(sqrtrhs+sqrt(4+rhs))
        fractions=(1+yH_+xHe)
        ss_=ss_ion*(fractions*(1.5*(lnTT_-lnTT_ion)-lnrho_+2.5) &
                   -yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                   -(1-yH_)*(log(1-yH_+epsi)-lnrho_H)-xHe_term)
        ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion
        pp_=(1+yH_+xHe)*exp(lnrho_)*TT_*ss_ion
        cs2_=impossible
!
      case (ilnrho_ee)
        lnrho_=var1
        ee_=var2
        yH_=0.5*min(ee_/ee_ion,yHmax)
        call rtsafe(ilnrho_ee,lnrho_,ee_,yHmin,yHmax*min(ee_/ee_ion,1.0),yH_)
        fractions=(1+yH_+xHe)
        TT_=(ee_-yH_*ee_ion)/(1.5*fractions*ss_ion)
        lnTT_=log(TT_)
        rho_=exp(lnrho_)
        ss_=ss_ion*(fractions*(1.5*(lnTT_-lnTT_ion)-lnrho_+2.5) &
                    -yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                    -(1-yH_)*(log(1-yH_+epsi)-lnrho_H)-xHe_term)
        pp_=fractions*rho_*TT_*ss_ion
        cs2_=impossible
!
      case (ilnrho_pp)
        lnrho_=var1
        pp_=var2
        yH_=0.5*yHmax
        call rtsafe(ilnrho_pp,lnrho_,pp_,yHmin,yHmax,yH_)
        fractions=(1+yH_+xHe)
        rho_=exp(lnrho_)
        TT_=pp_/(fractions*ss_ion*rho_)
        lnTT_=log(TT_)
        ss_=ss_ion*(fractions*(1.5*(lnTT_-lnTT_ion)-lnrho_+2.5) &
                   -yH_*(2*log(yH_)-lnrho_e-lnrho_H) &
                   -(1-yH_)*(log(1-yH_+epsi)-lnrho_H)-xHe_term)
        ee_=1.5*fractions*ss_ion*TT_+yH_*ee_ion
        cs2_=impossible
!
      case default
        call stop_it("eoscalc_point: I don't get what the independent variables are.")
!
      end select
!
      if (present(lnrho)) lnrho=lnrho_
      if (present(ss)) ss=ss_
      if (present(yH)) yH=yH_
      if (present(lnTT)) lnTT=lnTT_
      if (present(ee)) ee=ee_
      if (present(pp)) pp=pp_
      if (present(cs2)) cs2=cs2_
!
    endsubroutine eoscalc_point
!***********************************************************************
    subroutine read_eos_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=eos_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=eos_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_eos_init_pars
!***********************************************************************
    subroutine write_eos_init_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=eos_init_pars)
    endsubroutine write_eos_init_pars
!***********************************************************************
    subroutine read_eos_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=eos_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=eos_run_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_eos_run_pars
!***********************************************************************
    subroutine write_eos_run_pars(unit)
      integer, intent(in) :: unit
!
      write(unit,NML=eos_run_pars)
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
      real, dimension(mx) :: dyHold,dyH,yHlow,yHhigh,f,df
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
      f=lnrho_e-lnrho+1.5*lnTT_-TT1_+log(1-yH+epsi)-2*log(yH)
      dlnTT_=((2.0/3.0)*(-f-TT1_)-1)*fractions1
! wd: Need to add epsi at the end, as otherwise (yH-yHlow)*df  will
! wd: eventually yield an overflow.
! wd: Something like  sqrt(tini)  would probably also do instead of epsi,
! wd: but even epsi does not affect the auto-tests so far.
!      df=dlnTT_*(1.5+TT1_)-1/(1-yH+epsi)-2/yH
      df=dlnTT_*(1.5+TT1_)-1/(1-yH+epsi)-2/(yH+epsi)
!
      do i=1,maxit
        where (.not.found)
          where (      sign(1.,((yH-yHlow)*df-f)) &
                    == sign(1.,((yH-yHhigh)*df-f)) &
                  .or. abs(2*f) > abs(dyHold*df) )
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
            dyH=f/df
            ! Apply floor to dyH (necessary to avoid negative yH in samples
            ! /0d-tests/heating_ionize)
            dyH=min(dyH,yH-yHmin)
            dyH=max(dyH,yH-yHmax)    ! plausibly needed as well
            yH=yH-dyH
          endwhere
        endwhere
        where (abs(dyH)>yHacc*yH)
          fractions1=1/(1+yH+xHe)
          lnTT_=(2.0/3.0)*((ss/ss_ion+(1-yH)*(log(1-yH+epsi)-lnrho_H) &
                             +yH*(2*log(yH)-lnrho_e-lnrho_H) &
                             +xHe_term)*fractions1+lnrho-2.5)
          TT1_=exp(-lnTT_)
          f=lnrho_e-lnrho+1.5*lnTT_-TT1_+log(1-yH+epsi)-2*log(yH)
          dlnTT_=((2.0/3.0)*(-f-TT1_)-1)*fractions1
          df=dlnTT_*(1.5+TT1_)-1/(1-yH+epsi)-2/yH
          where (f<0)
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
      use Mpicomm, only: stop_it
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
        call stop_it("saha: I don't get what the independent variables are.")
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
    subroutine get_soundspeed(lnTT,cs2)
!
!  Calculate sound speed for given temperature
!
!  20-Oct-03/tobi: coded
!
      use Mpicomm
!
      real, intent(in)  :: lnTT
      real, intent(out) :: cs2
!
      call stop_it("get_soundspeed: with ionization, lnrho needs to be known here")
!
      call keep_compiler_quiet(lnTT,cs2)
!
    endsubroutine get_soundspeed
!***********************************************************************
    subroutine isothermal_entropy(f,T0)
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
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: T0
      real, dimension(nx) :: lnrho,ss,yH,K,sqrtK,yH_term,one_yH_term
!
      do n=n1,n2
      do m=m1,m2
!
        lnrho=f(l1:l2,m,n,ilnrho)
!
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
        ss=ss_ion*((1+yH+xHe)*(1.5*log(T0/TT_ion)-lnrho+2.5) &
                   -yH_term-one_yH_term-xHe_term)
!
        f(l1:l2,m,n,iss)=ss
!
      enddo
      enddo
!
    endsubroutine isothermal_entropy
!***********************************************************************
    subroutine isothermal_lnrho_ss(f,T0,rho0)
!
!  Isothermal stratification for lnrho and ss (for yH=0!)
!
!  Currently only implemented for ionization_fixed.
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, intent(in) :: T0,rho0
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(T0)
      call keep_compiler_quiet(rho0)
!
    endsubroutine isothermal_lnrho_ss
!***********************************************************************
     subroutine get_average_pressure(average_density,average_pressure)
!
!   01-dec-2009/piyali+dhrube: coded
!
      use Cdata
!
      real, intent(in):: average_density
      real, intent(out):: average_pressure
!
      call keep_compiler_quiet(average_density)
      call keep_compiler_quiet(average_pressure)
!
    endsubroutine get_average_pressure
!***********************************************************************
    subroutine bc_ss_flux(f,topbot)
!
!  constant flux boundary condition for entropy (called when bcz='c1')
!
!  23-jan-2002/wolf: coded
!  11-jun-2002/axel: moved into the entropy module
!   8-jul-2002/axel: split old bc_ss into two
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
      use Gravity
      use SharedVariables,only:get_shared_variable
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      logical, pointer :: lmultilayer, lheatc_chiconst
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my) :: tmp_xy,TT_xy,rho_xy,yH_xy
      integer :: i,ierr
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
!  Get the shared variables
!
      call get_shared_variable('hcond0',hcond0,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond0")
      call get_shared_variable('hcond1',hcond1,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting hcond1")
      call get_shared_variable('Fbot',Fbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting Fbot")
      call get_shared_variable('Ftop',Ftop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting Ftop")
      call get_shared_variable('FbotKbot',FbotKbot,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting FbotKbot")
      call get_shared_variable('FtopKtop',FtopKtop,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting FtopKtop")
      call get_shared_variable('chi',chi,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting chi")
      call get_shared_variable('lmultilayer',lmultilayer,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lmultilayer")
      call get_shared_variable('lheatc_chiconst',lheatc_chiconst,ierr)
      if (ierr/=0) call stop_it("bc_ss_flux: "//&
           "there was a problem when getting lheatc_chiconst")
!
      select case (topbot)
!
!  bottom boundary
!  ---------------
!
      case ('bot')
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
      case ('top')
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
          f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+ss_ion*(1+yH_xy+xHe)* &
              (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)-3*i*dz*tmp_xy)
        enddo
!
      case default
        print*,"bc_ss_flux: invalid argument"
        call stop_it("")
      endselect
!
    endsubroutine bc_ss_flux
!***********************************************************************
    subroutine bc_ss_flux_turb(f,topbot)
!
!  dummy routine
!
!   4-may-2009/axel: dummy routine
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
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
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_flux_turb_x
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
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_ss_temp_old: NOT IMPLEMENTED IN EOS_IONIZATION")
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
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_ss_temp_x: NOT IMPLEMENTED IN EOS_IONIZATION")
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
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_ss_temp_y: NOT IMPLEMENTED IN EOS_IONIZATION")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_ss_temp_z: NOT IMPLEMENTED IN EOS_IONIZATION")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_lnrho_temp_z(f,topbot)
!
!  boundary condition for density: constant temperature
!
!  19-aug-2005/tobi: distributed across ionization modules
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_lnrho_temp_z: NOT IMPLEMENTED IN EOS_IONIZATION")
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
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_lnrho_pressure_z: NOT IMPLEMENTED IN EOS_IONIZATION")
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
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_ss_temp2_z: NOT IMPLEMENTED IN EOS_IONIZATION")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_stemp_x(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_ss_stemp_x: NOT IMPLEMENTED IN EOS_IONIZATION")
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
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_ss_stemp_y: NOT IMPLEMENTED IN EOS_IONIZATION")
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
      use Mpicomm, only: stop_it
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,nghost) :: lnrho,ss,yH,lnTT,TT,K,sqrtK,yH_term,one_yH_term
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
      case ('bot')
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
      case ('top')
        do i=1,nghost
          f(:,:,n2+i,ilnTT) = f(:,:,n2-i,ilnTT)
        enddo
!
        lnrho=f(:,:,n2+1:mz,ilnrho)
        lnTT=f(:,:,n2+1:mz,ilnTT)
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
        f(:,:,n2+1:mz,iyH)=yH
        f(:,:,n2+1:mz,iss)=ss
!
      case default
        print*,"bc_ss_stemp_z: invalid argument"
        call stop_it("")
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
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_ss_a2stemp_x: NOT IMPLEMENTED IN EOS_IONIZATION")
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
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_ss_a2stemp_y: NOT IMPLEMENTED IN EOS_IONIZATION")
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
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_ss_a2stemp_z: NOT IMPLEMENTED IN EOS_IONIZATION")
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
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
!  The 'ce' boundary condition for entropy makes the energy constant at
!  the boundaries.
!  This assumes that the density is already set (ie density must register
!  first!)
!
      call stop_it("bc_ss_stemp_y: NOT IMPLEMENTED IN EOS_IONIZATION")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_stellar_surface: NOT IMPLEMENTED IN EOS_IDEALGAS")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_lnrho_cfb_r_iso: NOT IMPLEMENTED IN NOEOS")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_cfb_r_iso
!***********************************************************************
    subroutine bc_lnrho_hds_z_iso(f,topbot)
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_lnrho_hds_z_iso: NOT IMPLEMENTED IN NOEOS")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hds_z_iso
!***********************************************************************
    subroutine bc_lnrho_hdss_z_iso(f,topbot)
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
!
      call stop_it("bc_lnrho_hdss_z_iso: NOT IMPLEMENTED IN NOEOS")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine read_transport_data
!
       real, dimension (mx,my,mz,mfarray) :: f
!
       call keep_compiler_quiet(f)
!
    endsubroutine read_transport_data
!***********************************************************************
    subroutine write_thermodyn()
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine write_thermodyn
!***********************************************************************
    subroutine read_thermodyn(input_file)
!
      character (len=*), intent(in) :: input_file
!
      call keep_compiler_quiet(input_file)
!
    endsubroutine read_thermodyn
!***********************************************************************
    subroutine read_species(input_file)
!
      character (len=*) :: input_file
!
      call keep_compiler_quiet(input_file)
!
    endsubroutine read_species
!***********************************************************************
    subroutine find_species_index(species_name,ind_glob,ind_chem,found_specie)
!
      integer, intent(out) :: ind_glob
      integer, intent(inout) :: ind_chem
      character (len=*), intent(in) :: species_name
      logical, intent(out) :: found_specie
!
         call keep_compiler_quiet(ind_glob)
         call keep_compiler_quiet(ind_chem)
         call keep_compiler_quiet(species_name)
         call keep_compiler_quiet(found_specie)
!
     endsubroutine find_species_index
!***********************************************************************
     subroutine find_mass(element_name,MolMass)
!
      character (len=*), intent(in) :: element_name
      real, intent(out) :: MolMass
!
       call keep_compiler_quiet(element_name)
       call keep_compiler_quiet(MolMass)
!
     endsubroutine find_mass
!***********************************************************************
    subroutine read_Lewis

       real, dimension (mx,my,mz,mfarray) :: f

       call keep_compiler_quiet(f)

    endsubroutine read_Lewis
!***********************************************************************
endmodule EquationOfState
