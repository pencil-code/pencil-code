! $Id$
!
!  Routine for ideal gas with variable degree of ionization and hence
!  variable mean molecular weight. Here, the ionization fraction, yH,
!  is allocated as an additional auxiliary array in f.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: leos = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 1
!
! PENCILS PROVIDED ss; gss(3); ee; pp; lnTT; cs2; nabla_ad; glnTT(3); TT; TT1; gTT(3)
! PENCILS PROVIDED yH; del2ss; del2lnTT; del2TT; cv; cv1; cp; cp1; gamma; gamma_m1; gamma1
! PENCILS PROVIDED mu1; hlnTT(3,3); rho1gpp(3); delta; gradcp(3); del6lnTT
! PENCILS PROVIDED glnmumol(3); ppvap; csvap2; rho_anel
!
!***************************************************************
module EquationOfState
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include 'eos.h'
!
! integers specifying which independent variables to use in eoscalc
  integer, parameter :: ilnrho_ss=1,ilnrho_ee=2,ilnrho_pp=3,ilnrho_lnTT=4
  integer, parameter :: ilnrho_TT=5, irho_ss=7, irho_TT=10, ipp_ss=11
  integer, parameter :: ipp_cs2=12
  integer, parameter :: irho_eth=13, ilnrho_eth=14
  integer :: icv, ics, idelta, igamma, inabad
!  integer :: icp, icv, ics, idelta, igamma, inabad
  integer :: imass = 1
  !  secondary parameters calculated in initialize
  real :: mu1_0,Rgas
  real :: TT_ion,lnTT_ion,TT_ion_,lnTT_ion_
  real :: ss_ion,kappa0,pp_ion
  real :: rho_H,rho_e,rho_e_,rho_He
  real :: lnrho_H,lnrho_e,lnrho_e_,lnrho_He
!
  real :: xHe=0.1, yH_const=0.0, yMetals=0.0, tau_relax=1.0
  real :: chiH_eV=13.6, chiHminus_eV=0.754
  logical :: lrevise_chiH_eV=.false., lrevise_chiHminus_eV=.false.
  logical :: lconst_yH=.false., lHminus_opacity_correction=.false.
!
  real :: lnpp_bot=0.
  real :: ss_bot=0.
  real :: TTbot, TTtop
!
  real :: va2max_eos=huge1
  integer :: va2power_eos=5
  real, dimension (3) :: B_ext_eos=(/0.,0.,0./)
  real :: cs0=impossible, rho0=impossible, cp=impossible,cv=impossible
  real :: cs20=impossible, lnrho0=impossible
  logical :: lcalc_cp=.false.,lcalc_cp_full=.false.
  logical :: lss_as_aux=.false., lpp_as_aux=.false., lcs_as_aux=.false.
  logical :: lcp_as_aux=.false., lcv_as_aux=.false., lgamma_as_aux=.false.
  logical :: lnabad_as_aux=.false., ldelta_as_aux=.false.
  real :: gamma=5./3., gamma_m1=impossible, gamma1=impossible
!
! init parameters
!
  namelist /eos_init_pars/ xHe,lconst_yH,yH_const,yMetals,lnpp_bot,ss_bot, &
                           lrevise_chiH_eV, chiH_eV, &
                           lrevise_chiHminus_eV, chiHminus_eV, &
                           tau_relax,va2max_eos,va2power_eos,B_ext_eos, &
                           lss_as_aux,lpp_as_aux,lcp_as_aux,lcv_as_aux, &
                           lcs_as_aux,lgamma_as_aux,lnabad_as_aux, &
                           ldelta_as_aux, &
                           lHminus_opacity_correction, TTbot, TTtop
!
! run parameters
!
  namelist /eos_run_pars/ xHe,lconst_yH,yH_const,yMetals,lnpp_bot,ss_bot, &
                          lrevise_chiH_eV, chiH_eV, &
                          lrevise_chiHminus_eV, chiHminus_eV, &
                          tau_relax,va2max_eos,va2power_eos,B_ext_eos, &
                          lss_as_aux,lpp_as_aux,lcp_as_aux,lcv_as_aux, &
                          lcs_as_aux,lgamma_as_aux,lnabad_as_aux, &
                          ldelta_as_aux, &
                          lHminus_opacity_correction
!
  real :: cs2bot=impossible, cs2top=impossible
!
! Allocatable 3D-array for cp
!
  real, dimension (:,:,:), allocatable :: cp_full
!
  real, dimension(nchemspec,18) :: species_constants
!
  real :: Cp_const=impossible
  real :: Pr_number=0.7
  logical :: lpres_grad=.false.
!
  contains
!***********************************************************************
    subroutine register_eos
!
!  14-jun-03/axel: adapted from register_eos
!
      use FArrayManager
!
      leos_temperature_ionization=.true.
!
!  Set indices for auxiliary variables.
!
      call farray_register_auxiliary('yH',iyH)
!
!  Writing files for use with IDL
!
      if (naux < maux)  aux_var(aux_count)=',yH $'
      if (naux == maux) aux_var(aux_count)=',yH'
      aux_count=aux_count+1
      if (lroot) then
        write(15,*) 'yH = fltarr(mx,my,mz)*one'
      endif
!
!  Identify version number.
!
      if (lroot) call svn_id( &
           '$Id$')
!
    endsubroutine register_eos
!***********************************************************************
    subroutine initialize_eos
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  21-may-14/axel: adapted from eos_entropy
!
      use Sub, only: register_report_aux
      use SharedVariables, only: put_shared_variable
!
      if (lroot) print*,'initialize_eos: ENTER'
!
!  Useful constants for ionization
!  (Here we assume m_H = m_p = m_u, m_He = 4*m_u, and m_e << m_u)
!
      if (lrevise_chiH_eV) chiH=chiH_eV*eV
      if (lrevise_chiHminus_eV) chiH_=chiHminus_eV*eV
      if (lroot) print*,'initialize_eos: chiH=',chiH,chiH/eV
      if (lroot) print*,'initialize_eos: chiH_=',chiH_,chiH_/eV
!
      mu1_0 = 1/(1 + 4*xHe)
      Rgas = k_B/m_u
      TT_ion = chiH/k_B
      lnTT_ion = log(TT_ion)
      TT_ion_ = chiH_/k_B
      lnTT_ion_ = log(TT_ion_)
      rho_H = (1/mu1_0)*m_u*((m_u/hbar)*(chiH/hbar)/(2*pi))**(1.5)
      lnrho_H = log(rho_H)
      rho_e = (1/mu1_0)*m_u*((m_e/hbar)*(chiH/hbar)/(2*pi))**(1.5)
      lnrho_e = log(rho_e)
      rho_He = (1/mu1_0)*m_u*((4*m_u/hbar)*(chiH/hbar)/(2*pi))**(1.5)
      lnrho_He = log(rho_He)
      rho_e_ = (1/mu1_0)*m_u*((m_e/hbar)*(chiH_/hbar)/(2*pi))**(1.5)
      lnrho_e_ = log(rho_e_)
      kappa0 = sigmaH_*mu1_0/(4*m_u)
      pp_ion = Rgas*mu1_0*rho_e*TT_ion
!
!  pressure, cp, and cv as optional auxiliary variable
!
      if (lss_as_aux) call register_report_aux('sss',iss)
      if (lpp_as_aux) call register_report_aux('ppp',ipp)
      if (lcp_as_aux) call register_report_aux('cp',icp)
      if (lcv_as_aux) call register_report_aux('cv',icv)
      if (lcs_as_aux) call register_report_aux('cs',ics)
      if (ldelta_as_aux) call register_report_aux('delta',idelta)
      if (lgamma_as_aux) call register_report_aux('gamma',igamma)
      if (lnabad_as_aux) call register_report_aux('nabad',inabad)

      if (.not.ldensity) then
        call put_shared_variable('rho0',rho0,caller='initialize_eos')
        call put_shared_variable('lnrho0',lnrho0)
      else
        call put_shared_variable('TTtop',TTtop,caller='initialize_eos')
      endif
!
!  write scale non-free constants to file; to be read by idl
!
      if (lroot) then
        open (1,file=trim(datadir)//'/pc_constants.pro',position="append")
        write (1,*) 'mu1_0=',mu1_0
        write (1,*) 'Rgas=',Rgas
        write (1,*) 'TT_ion=',TT_ion
        write (1,*) 'lnTT_ion=',lnTT_ion
        write (1,*) 'TT_ion_=',TT_ion_
        write (1,*) 'lnTT_ion_=',lnTT_ion_
        write (1,*) 'lnrho_H=',lnrho_H
        write (1,*) 'lnrho_e=',lnrho_e
        write (1,*) 'lnrho_He=',lnrho_He
        write (1,*) 'lnrho_e_=',lnrho_e_
        write (1,*) 'kappa0=',kappa0
        close (1)
      endif
!
    endsubroutine initialize_eos
!***********************************************************************
    subroutine select_eos_variable(variable,findex)
!
!  dummy (but to be changed)
!
      character (len=*), intent(in) :: variable
      integer, intent(in) :: findex
!
      call keep_compiler_quiet(variable)
      call keep_compiler_quiet(findex)
!
    endsubroutine select_eos_variable
!***********************************************************************
    subroutine pencil_criteria_eos
!
!  All pencils that the EquationOfState module depends on are specified here.
!
!  21-may-14/axel: adapted from eos_ionization
!
!  EOS is a pencil provider but evolves nothing so it is unlokely that
!  it will require any pencils for it's own use.
!
!  pp pencil if lpp_as_aux
!
      if (lss_as_aux) lpenc_requested(i_ss)=.true.
      if (lpp_as_aux) lpenc_requested(i_pp)=.true.
      if (lcp_as_aux) lpenc_requested(i_cp)=.true.
      if (lcv_as_aux) lpenc_requested(i_cv)=.true.
      if (lcs_as_aux) lpenc_requested(i_cs2)=.true.
      if (ldelta_as_aux) lpenc_requested(i_delta)=.true.
      if (lgamma_as_aux) lpenc_requested(i_gamma)=.true.
      if (lnabad_as_aux) lpenc_requested(i_nabla_ad)=.true.
!
    endsubroutine pencil_criteria_eos
!***********************************************************************
    subroutine pencil_interdep_eos(lpencil_in)
!
!  dummy (but to be changed)
!
      logical, dimension(npencils) :: lpencil_in
!
      if (lpencil_in(i_cs2)) then
        lpencil_in(i_gamma)=.true.
        lpencil_in(i_rho1)=.true.
        lpencil_in(i_pp)=.true.
      endif
!
      if (lpencil_in(i_rho1gpp)) then
        lpencil_in(i_gamma1)=.true.
        lpencil_in(i_cs2)=.true.
        lpencil_in(i_glnrho)=.true.
        lpencil_in(i_delta)=.true.
        lpencil_in(i_glnTT)=.true.
      endif
!
      if (lpencil_in(i_nabla_ad)) then
        lpencil_in(i_mu1)=.true.
        lpencil_in(i_delta)=.true.
        lpencil_in(i_cp1)=.true.
      endif
!
      if (lpencil_in(i_gamma_m1)) lpencil_in(i_gamma)=.true.
!
      if (lpencil_in(i_gamma)) then
        lpencil_in(i_cp)=.true.
        lpencil_in(i_cv1)=.true.
      endif
!
      if (lpencil_in(i_gamma1)) then
        lpencil_in(i_cv)=.true.
        lpencil_in(i_cp1)=.true.
      endif
!
      if (lpencil_in(i_cv1)) lpencil_in(i_cv)=.true.
!
      if (lpencil_in(i_cv)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_TT1)=.true.
        lpencil_in(i_mu1)=.true.
      endif
!
      if (lpencil_in(i_cp1)) lpencil_in(i_cp)=.true.
!
      if (lpencil_in(i_cp)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_TT1)=.true.
        lpencil_in(i_mu1)=.true.
      endif
!
      if (lpencil_in(i_pp)) then
        lpencil_in(i_mu1)=.true.
        lpencil_in(i_rho)=.true.
        lpencil_in(i_TT)=.true.
      endif
!
      if (lpencil_in(i_mu1)) lpencil_in(i_yH)=.true.
!
      if (lpencil_in(i_ee)) then
        lpencil_in(i_mu1)=.true.
        lpencil_in(i_TT)=.true.
        lpencil_in(i_yH)=.true.
      endif
!
      if (lpencil_in(i_ss)) then
        lpencil_in(i_yH)=.true.
        lpencil_in(i_lnrho)=.true.
        lpencil_in(i_lnTT)=.true.
      endif
!
      if (lpencil_in(i_TT)) lpencil_in(i_lnTT)=.true.
      if (lpencil_in(i_TT1)) lpencil_in(i_TT)=.true.
!
      if (lpencil_in(i_gradcp)) lcalc_cp_full=.true.
!
      if (lpencil_in(i_gTT)) then
        lpencil_in(i_TT)=.true.
        lpencil_in(i_glnTT)=.true.
      endif
!
      if (lpencil_in(i_del2TT)) then
        lpencil_in(i_TT)=.true.
        lpencil_in(i_glnTT)=.true.
        lpencil_in(i_del2lnTT)=.true.
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
      real, dimension (mx,my,mz,mfarray),intent(IN)   :: f
      type (pencil_case),                intent(INOUT):: p
!
      call calc_pencils_eos_pencpar(f,p,lpencil)
!
    endsubroutine calc_pencils_eos_std
!***********************************************************************
    subroutine calc_pencils_eos_pencpar(f,p,lpenc_loc)
!
!  Calculate relevant eos pencils
!
!   9-oct-15/MR: added mask on pencil case as a parameter.
!
      use Sub, only: grad,del2,del6,g2ij
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
      logical, dimension(npencils) :: lpenc_loc
!
      real, dimension (nx) :: yH_term_cv,TT_term_cv
      real, dimension (nx) :: yH_term_cp,TT_term_cp
      real, dimension (nx) :: alpha1,tmp
      integer :: i
!
!  Temperature
!
      if (lpenc_loc(i_lnTT)) p%lnTT=f(l1:l2,m,n,ilnTT)
      if (lpenc_loc(i_TT)) p%TT=exp(p%lnTT)
      if (lpenc_loc(i_TT1)) p%TT1=1/p%TT
!
!
!  Temperature laplacian and gradient
!
      if (lpenc_loc(i_glnTT)) call grad(f,ilnTT,p%glnTT)
      if (lpenc_loc(i_hlnTT)) call g2ij(f,ilnTT,p%hlnTT)
      if (lpenc_loc(i_del2lnTT)) call del2(f,ilnTT,p%del2lnTT)
      if (lpenc_loc(i_del2TT)) then
        tmp=0.0
        do i=1,3
          tmp=tmp+p%glnTT(:,i)**2
        enddo
        p%del2TT=(p%del2lnTT+tmp)*p%TT
      endif
!
      if (lpenc_loc(i_del6lnTT)) call del6(f,ilnTT,p%del6lnTT)
      if (lpenc_loc(i_gTT)) then
        do i=1,3
          p%gTT(:,i) =p%TT * p%glnTT(:,i)
        enddo
      endif
!
!  Ionization fraction
!
      if (lpenc_loc(i_yH)) p%yH = f(l1:l2,m,n,iyH)
!
!  Mean molecular weight
!
      if (lpenc_loc(i_mu1)) p%mu1 = mu1_0*(1 + p%yH + xHe)
!
!  Pressure
!
      if (lpenc_loc(i_pp)) p%pp = Rgas*p%mu1*p%rho*p%TT
!
!  Common terms involving the ionization fraction
!
      if (lpenc_loc(i_cv)) then
        yH_term_cv = p%yH*(1-p%yH)/((2-p%yH)*(1+p%yH+xHe))
        TT_term_cv = 1.5 + p%TT1*TT_ion
      endif
!
      if (lpenc_loc(i_cp).or.lpenc_loc(i_delta)) then
        yH_term_cp = p%yH*(1-p%yH)/(2+xHe*(2-p%yH))
        TT_term_cp = 2.5 + p%TT1*TT_ion
      endif
!
!  Specific heat at constant volume (i.e. density)
!
      if (lpenc_loc(i_cv)) p%cv = Rgas*p%mu1*(1.5 + yH_term_cv*TT_term_cv**2)
      if (lpenc_loc(i_cv1)) p%cv1=1/p%cv
!
!  Specific heat at constant pressure
!
      if (lpenc_loc(i_cp)) then
        if (lcalc_cp_full) then
          p%cp = cp_full(l1:l2,m,n)
        else
          p%cp = Rgas*p%mu1*(2.5 + yH_term_cp*TT_term_cp**2)
        endif
      endif
      if (lpenc_loc(i_cp1)) p%cp1 = 1/p%cp
!
!  Gradient of the above
!
      if (lpenc_loc(i_gradcp)) call grad(cp_full,p%gradcp)
!
!  Polytropic index
!
      if (lpenc_loc(i_gamma)) p%gamma = p%cp*p%cv1
      if (lpenc_loc(i_gamma1)) p%gamma1 = p%cv*p%cp1
      if (lpenc_loc(i_gamma_m1)) p%gamma_m1 = p%gamma - 1
!
!  For the definition of delta, see Kippenhahn & Weigert
!
      if (lpenc_loc(i_delta)) p%delta = 1 + yH_term_cp*TT_term_cp
!
!  Sound speed
!
      if (lpenc_loc(i_cs2)) then
        alpha1 = (2+xHe*(2-p%yH))/((2-p%yH)*(1+p%yH+xHe))
        p%cs2 = p%gamma*p%rho1*p%pp*alpha1
      endif
!
!  Adiabatic temperature gradient
!
      if (lpenc_loc(i_nabla_ad)) p%nabla_ad = Rgas*p%mu1*p%delta*p%cp1
!
!  Logarithmic pressure gradient
!
      if (lpenc_loc(i_rho1gpp)) then
        do i=1,3
          p%rho1gpp(:,i) = p%gamma1*p%cs2*(p%glnrho(:,i)+p%delta*p%glnTT(:,i))
        enddo
      endif
!
!  Energy per unit mass
!
      if (lpenc_loc(i_ee)) p%ee = 1.5*Rgas*p%mu1*p%TT + p%yH*Rgas*mu1_0*TT_ion
!
!  Entropy per unit mass
!  The contributions from each particle species contain the mixing entropy
!
      if (lpenc_loc(i_ss)) then
        tmp = 2.5 - 1.5*(lnTT_ion-p%lnTT) - p%lnrho
        where (p%yH < 1) ! Neutral Hydrogen
          p%ss = (1-p%yH)*(tmp + lnrho_H - log(1-p%yH))
        endwhere
        where (p%yH > 0) ! Electrons and ionized Hydrogen
          p%ss = p%ss + p%yH*(tmp + lnrho_H - log(p%yH))
          p%ss = p%ss + p%yH*(tmp + lnrho_e - log(p%yH))
        endwhere
        if (xHe > 0) then ! Helium
          p%ss = p%ss + xHe*(tmp + lnrho_He - log(xHe))
        endif
        p%ss = Rgas*mu1_0*p%ss
      endif
!
      if (lpenc_loc(i_gss)) then
        call fatal_error('gss',"SHOULDN'T BE CALLED WITH eos_temperature_...")
      endif
!
      if (lpenc_loc(i_del2ss)) then
        call fatal_error('del2ss',"SHOULDN'T BE CALLED WITH eos_temperature_...")
      endif
!
      if (lpenc_loc(i_glnmumol)) p%glnmumol(:,:)=0.
!
!  pressure and cp as optional auxiliary variables
!
      if (lss_as_aux) f(l1:l2,m,n,iss)=p%ss
      if (lpp_as_aux) f(l1:l2,m,n,ipp)=p%pp
      if (lcp_as_aux) f(l1:l2,m,n,icp)=p%cp
      if (lcv_as_aux) f(l1:l2,m,n,icv)=p%cv
      if (lcs_as_aux) f(l1:l2,m,n,ics)=sqrt(p%cs2)
      if (ldelta_as_aux) f(l1:l2,m,n,idelta)=p%delta
      if (lgamma_as_aux) f(l1:l2,m,n,igamma)=p%gamma
      if (lnabad_as_aux) f(l1:l2,m,n,inabad)=p%nabla_ad
!
    endsubroutine calc_pencils_eos_pencpar
!***********************************************************************
    subroutine getmu(f,mu_tmp)
!
!  Calculate average particle mass in the gas relative to
!
!   12-aug-03/tony: implemented dummy
!
      real, dimension (mx,my,mz,mfarray), optional :: f
      real, optional, intent(out) :: mu_tmp
!
      mu_tmp=0.
      call keep_compiler_quiet(present(f))
!
    endsubroutine getmu
!***********************************************************************
    subroutine rprint_eos(lreset,lwrite)
!
      logical :: lreset
      logical, optional :: lwrite
!
      call keep_compiler_quiet(lreset)
      call keep_compiler_quiet(present(lwrite))
!
    endsubroutine rprint_eos
!***********************************************************************
    subroutine ioninit(f)
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine ioninit
!***********************************************************************
    subroutine ioncalc(f)
!
!   calculate degree of ionization and temperature
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      real, dimension (mx) :: yH,rho1,TT1,rhs,sqrtrhs
      real, dimension (mx) :: mu1,yH_term_cp,TT_term_cp
!
      if (.not.allocated(cp_full)) allocate (cp_full(mx,my,mz))
!
      if (lconst_yH) then
!
        f(:,:,:,iyH) = yH_const
!
      else
!
        do n=1,mz
        do m=1,my
!
          rho1 = exp(-f(:,m,n,ilnrho))
          TT1 = exp(-f(:,m,n,ilnTT))
!
          rhs = rho_e*rho1*(TT1*TT_ion)**(-1.5)*exp(-TT_ion*TT1)
          sqrtrhs = sqrt(rhs)
          yH = 2*sqrtrhs/(sqrtrhs+sqrt(4+rhs))
!
          f(:,m,n,iyH) = yH
!
          if (lcalc_cp_full) then
            mu1 = mu1_0*(1 + yH + xHe)
            yH_term_cp = yH*(1-yH)/(2+xHe*(2-yH))
            TT_term_cp = 2.5 + TT_ion*TT1
            cp_full(:,m,n) = Rgas*mu1*(2.5 + yH_term_cp*TT_term_cp**2)
          endif
!
        enddo
        enddo
!
      endif
!
    endsubroutine ioncalc
!***********************************************************************
    subroutine eosperturb(f,psize,ee,pp,ss)
!
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(in), optional :: ee,pp,ss
!
      call not_implemented("eosperturb")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(psize)
      call keep_compiler_quiet(present(ee))
      call keep_compiler_quiet(present(pp))
      call keep_compiler_quiet(present(ss))
!
    endsubroutine eosperturb
!***********************************************************************
    subroutine getdensity(EE,TT,yH,rho)
!
      real, intent(in) :: EE,TT,yH
      real, intent(inout) :: rho
!
      call fatal_error('getdensity','SHOULD NOT BE CALLED WITH'// &
          ' eos_temperature_ionization')
      rho = 0.
!
      call keep_compiler_quiet(yH)
      call keep_compiler_quiet(EE)
      call keep_compiler_quiet(TT)
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
    subroutine getpressure(pp_tmp,TT_tmp,rho_tmp,mu1_tmp)
!
     real, dimension (nx), intent(out) :: pp_tmp
     real, dimension (nx), intent(in)  :: TT_tmp,rho_tmp,mu1_tmp
!
     call fatal_error('getpressure','Should not be called with noeos.')
!
     call keep_compiler_quiet(pp_tmp)
     call keep_compiler_quiet(TT_tmp)
     call keep_compiler_quiet(rho_tmp)
     call keep_compiler_quiet(mu1_tmp)
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
      call fatal_error('get_cp1',"SHOULDN'T BE CALLED WITH eos_temperature_...")
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
      call fatal_error('get_cv1',"SHOULDN'T BE CALLED WITH eos_temperature_...")
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
!   02-apr-04/tony: implemented dummy
!
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(out) :: cs2,cp1tilde
!
!
      cs2=impossible
      cp1tilde=impossible
      call fatal_error('pressure_gradient_farray','SHOULD NOT BE CALLED WITH'// &
          ' eos_temperature_ionization')
!
      call keep_compiler_quiet(f)
!
    endsubroutine pressure_gradient_farray
!***********************************************************************
    subroutine pressure_gradient_point(lnrho,ss,cs2,cp1tilde)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   02-apr-04/tony: implemented dummy
!
      real, intent(in) :: lnrho,ss
      real, intent(out) :: cs2,cp1tilde
!
      cs2=impossible
      cp1tilde=impossible
      call fatal_error('pressure_gradient_point','SHOULD NOT BE CALLED WITH'// &
          ' eos_temperature_ionization')
!
      call keep_compiler_quiet(lnrho)
      call keep_compiler_quiet(ss)
!
    endsubroutine pressure_gradient_point
!***********************************************************************
    subroutine temperature_gradient(f,glnrho,gss,glnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   02-apr-04/tony: implemented dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(in) :: glnrho,gss
      real, dimension(nx,3), intent(out) :: glnTT
!
      glnTT=0.
      call fatal_error('temperature_gradient','SHOULD NOT BE CALLED WITH'// &
          ' eos_temperature_ionization')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(glnrho)
      call keep_compiler_quiet(gss)
!
    endsubroutine temperature_gradient
!***********************************************************************
    subroutine temperature_laplacian(f,del2lnrho,del2ss,del2lnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally glnPP and glnTT
!   gP/rho=cs2*(glnrho+cp1tilde*gss)
!
!   12-dec-05/tony: adapted from subroutine temperature_gradient
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx), intent(in) :: del2lnrho,del2ss
      real, dimension(nx), intent(out) :: del2lnTT
!
      call fatal_error('temperature_laplacian','SHOULD NOT BE CALLED WITH'// &
          ' eos_temperature_ionization')
!
      del2lnTT=0.
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(del2lnrho)
      call keep_compiler_quiet(del2ss)
!
    endsubroutine temperature_laplacian
!***********************************************************************
    subroutine temperature_hessian(f,hlnrho,hss,hlnTT)
!
!   Calculate thermodynamical quantities, cs2 and cp1tilde
!   and optionally hlnPP and hlnTT
!   hP/rho=cs2*(hlnrho+cp1tilde*hss)
!
!   13-may-04/tony: adapted from idealgas dummy
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, dimension(nx,3), intent(in) :: hlnrho,hss
      real, dimension(nx,3) :: hlnTT
!
      call fatal_error('temperature_hessian','now I do not believe you'// &
                           ' intended to call this!')
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(hlnrho)
      call keep_compiler_quiet(hss)
      call keep_compiler_quiet(hlnTT)
!
    endsubroutine temperature_hessian
!***********************************************************************
    subroutine eoscalc_farray(f,psize,lnrho,ss,yH,mu1,lnTT,ee,pp,cs2,kapparho)
!
!   Calculate thermodynamical quantities
!
!   04-apr-06/tobi: Adapted for this EOS module
!   27-jan-11/MR: caught zero in calculation of kapparho
!
      use Mpicomm, only: stop_it
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      integer, intent(in) :: psize
      real, dimension(psize), intent(out), optional :: lnrho,ss
      real, dimension(psize), intent(out), optional :: yH,lnTT,mu1
      real, dimension(psize), intent(out), optional :: ee,pp,kapparho
      real, dimension(psize), optional :: cs2
!
      real, dimension(psize) :: lnrho_,lnTT_,yH_,mu1_
      real, dimension(psize) :: TT1,tmp,tmpy,tmpy1
!
      select case (psize)
!
      case (nx)
        lnrho_=f(l1:l2,m,n,ilnrho)
        lnTT_=f(l1:l2,m,n,ilnTT)
        yH_=f(l1:l2,m,n,iyH)
!
      case (mx)
        lnrho_=f(:,m,n,ilnrho)
        lnTT_=f(:,m,n,ilnTT)
        yH_=f(:,m,n,iyH)
!
      case default
        call stop_it("eoscalc: no such pencil size")
!
      end select
!
      if (present(lnrho)) lnrho=lnrho_
!
      if (present(lnTT)) lnTT=lnTT_
!
      if (present(yH)) yH = yH_
!
      if (present(mu1).or.present(ss).or.present(ee).or.present(pp)) then
        mu1 = mu1_0*(1 + yH_ + xHe)
      endif
!
      if (present(ee)) ee = 1.5*Rgas*mu1*exp(lnTT_) + yH_*Rgas*mu1_0*TT_ion
!
      if (present(pp)) pp = Rgas*mu1*exp(lnrho_+lnTT_)
!
      if (present(ss)) then
        tmp = 2.5 - 1.5*(lnTT_ion-lnTT_) - lnrho_
        where (yH_ < 1) ! Neutral Hydrogen
          ss = (1-yH_)*(tmp + lnrho_H - log(1-yH_))
        endwhere
        where (yH_ > 0) ! Electrons and ionized Hydrogen
          ss = ss + yH_*(tmp + lnrho_H - log(yH_))
          ss = ss + yH_*(tmp + lnrho_e - log(yH_))
        endwhere
        if (xHe > 0) then ! Helium
          ss = ss + xHe*(tmp + lnrho_He - log(xHe))
        endif
        ss = Rgas*mu1_0*ss
      endif
!
      if (present(cs2)) &
        call fatal_error('eoscalc_farray','calculation of cs2 not implemented')
!
      if (present(kapparho)) then
!
!  Note: the factor of 2 in front of lnrho_ is because we compute kappa*rho.
!
        TT1 = exp(-lnTT_)
        tmp = 2*lnrho_-lnrho_e_+1.5*(lnTT_ion_-lnTT_)+TT_ion_*TT1
        tmpy = yH_+yMetals
!
        where ( tmpy==0. )              ! assumes that tmpy >= 0
          kapparho=0.
        elsewhere
          kapparho = (1-yH_)*kappa0*exp(min(tmp,log(huge1))+alog(tmpy))
        endwhere
        if (lHminus_opacity_correction) then
          mu1_ = mu1_0*(1+yH_+xHe)
          tmpy1 = 1./(1+yH_)
          kapparho = kapparho*(1+4*xHe)*mu1_*tmpy1
        endif
      endif
!
    endsubroutine eoscalc_farray
!***********************************************************************
    subroutine eoscalc_point(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
!   Calculate thermodynamical quantities
!
!   02-apr-04/tony: implemented dummy
!
      integer, intent(in) :: ivars
      real, intent(in) :: var1,var2
      real, intent(out), optional :: lnrho,ss
      real, intent(out), optional :: yH,lnTT
      real, intent(out), optional :: ee,pp
!
!
      call fatal_error('eoscalc_point','SHOULD NOT BE CALLED WITH'// &
          ' eos_temperature_ionization')
!
      if (present(lnrho)) lnrho=0.
      if (present(ss)) ss=0.
      if (present(yH)) yH=impossible
      if (present(lnTT)) lnTT=0.
      if (present(ee)) ee=0.
      if (present(pp)) pp=0.
!
      call keep_compiler_quiet(ivars)
      call keep_compiler_quiet(var1)
      call keep_compiler_quiet(var2)
!
    endsubroutine eoscalc_point
!***********************************************************************
    subroutine eoscalc_pencil(ivars,var1,var2,lnrho,ss,yH,lnTT,ee,pp)
!
!   Calculate thermodynamical quantities
!
!   2-feb-03/axel: simple example coded
!   13-jun-03/tobi: the ionization fraction as part of the f-array
!                   now needs to be given as an argument as input
!   17-nov-03/tobi: moved calculation of cs2 and cp1tilde to
!                   subroutine pressure_gradient
!
      integer, intent(in) :: ivars
      real, dimension(nx), intent(in) :: var1,var2
      real, dimension(nx), intent(out), optional :: lnrho,ss
      real, dimension(nx), intent(out), optional :: yH,lnTT
      real, dimension(nx), intent(out), optional :: ee,pp
!
      call fatal_error('eoscalc_pencil','SHOULD NOT BE CALLED WITH'// &
          ' eos_temperature_ionization')
!
      if (present(lnrho)) lnrho=0.
      if (present(ss)) ss=0.
      if (present(yH)) yH=impossible
      if (present(lnTT)) lnTT=0.
      if (present(ee)) ee=0.
      if (present(pp)) pp=0.
!
      call keep_compiler_quiet(ivars)
      call keep_compiler_quiet(var1)
      call keep_compiler_quiet(var2)
!
    endsubroutine eoscalc_pencil
!***********************************************************************
    subroutine get_soundspeed(TT,cs2)
!
!  Calculate sound speed for given temperature
!
!  02-apr-04/tony: dummy coded
!
      real, intent(in)  :: TT
      real, intent(out) :: cs2
!
      call not_implemented('get_soundspeed')
      cs2=0.
!
      call keep_compiler_quiet(TT)
!
    endsubroutine get_soundspeed
!***********************************************************************
    subroutine units_eos
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
    subroutine get_slices_eos(f,slices)
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_eos
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
      real, dimension(nx) :: lnrho,ss
      real :: ss_offset=0.
!
!  if T0 is different from unity, we interpret
!  ss_offset = ln(T0)/gamma as an additive offset of ss
!
      if (T0/=1.) ss_offset=alog(T0)/gamma
!
      do n=n1,n2
      do m=m1,m2
        lnrho=f(l1:l2,m,n,ilnrho)
        ss=-gamma_m1*(lnrho-lnrho0)/gamma
          !+ other terms for sound speed not equal to cs_0
        f(l1:l2,m,n,iss)=ss+ss_offset
      enddo
      enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      cs2bot=cs20
      cs2top=cs20
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
!   01-dec-2009/piyali+dhrube: coded
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
      use Mpicomm,only:stop_it
!
      real, pointer :: Fbot,Ftop,FtopKtop,FbotKbot,hcond0,hcond1,chi
      logical, pointer :: lmultilayer, lheatc_chiconst
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      logical, optional :: lone_sided
      real, dimension (size(f,1),size(f,2)) :: tmp_xy,cs2_xy,rho_xy
      integer :: i
!
      if (ldebug) print*,'bc_ss_flux: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c1' boundary condition (constant heat flux) for entropy.
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!  ===============
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
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n1,ilnrho)-lnrho0)+gamma*f(:,:,n1,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!
        if (lheatc_chiconst) then
          tmp_xy=Fbot/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FbotKbot/cs2_xy
        endif
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n1-i,iss)=f(:,:,n1+i,iss)+gamma_m1/gamma* &
              (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)+dz2_bound(-i)*tmp_xy)
        enddo
!
!  top boundary
!  ============
!
      case ('top')
        if (lmultilayer) then
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0*hcond1
        else
          if (headtt) print*,'bc_ss_flux: Ftop,hcond=',Ftop,hcond0
        endif
!
!  calculate Ftop/(K*cs2)
!
        rho_xy=exp(f(:,:,n2,ilnrho))
        cs2_xy=cs20*exp(gamma_m1*(f(:,:,n2,ilnrho)-lnrho0)+gamma*f(:,:,n2,iss))
!
!  check whether we have chi=constant at bottom, in which case
!  we have the nonconstant rho_xy*chi in tmp_xy.
!
        if (lheatc_chiconst) then
          tmp_xy=Ftop/(rho_xy*chi*cs2_xy)
        else
          tmp_xy=FtopKtop/cs2_xy
        endif
!
!  enforce ds/dz + gamma_m1/gamma*dlnrho/dz = - gamma_m1/gamma*Fbot/(K*cs2)
!
        do i=1,nghost
          f(:,:,n2+i,iss)=f(:,:,n2-i,iss)+gamma_m1/gamma* &
              (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)-dz2_bound(i)*tmp_xy)
        enddo
      case default
        call fatal_error('bc_ss_flux','invalid argument')
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
      real, dimension (:,:,:,:) :: f
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
!   31-may-2010/axel: dummy routine
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
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
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
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
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
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
      character (len=3) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,1),size(f,2)) :: tmp_xy
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp_old: ENTER - cs20,cs0=',cs20,cs0
!
!  Do the `c2' boundary condition (fixed temperature/sound speed) for entropy.
!  This assumes that the density is already set (ie density must register
!  first!)
!  tmp_xy = s(x,y) on the boundary.
!  gamma*s/cp = [ln(cs2/cs20)-(gamma-1)ln(rho/rho0)]
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if ((bcz12(ilnrho,1) /= 'a2') .and. (bcz12(ilnrho,1) /= 'a3')) &
          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3.')
        if (ldebug) print*, &
                'bc_ss_temp_old: set bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) &
              print*,'bc_ss_temp_old: cannot have cs2bot<=0'
        tmp_xy = (-gamma_m1*(f(:,:,n1,ilnrho)-lnrho0) &
             + alog(cs2bot/cs20)) / gamma
        f(:,:,n1,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n1-i,iss) = 2*tmp_xy - f(:,:,n1+i,iss)
        enddo
!
!  top boundary
!
      case ('top')
        if ((bcz12(ilnrho,2) /= 'a2') .and. (bcz12(ilnrho,2) /= 'a3')) &
          call fatal_error('bc_ss_temp_old','Inconsistent boundary conditions 3.')
        if (ldebug) print*, &
                   'bc_ss_temp_old: set top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                   'bc_ss_temp_old: cannot have cs2top<=0'
  !     if (bcz12(ilnrho,1) /= 'a2') &
  !          call fatal_error(bc_ss_temp_old','Inconsistent boundary conditions 4.')
        tmp_xy = (-gamma_m1*(f(:,:,n2,ilnrho)-lnrho0) &
                 + alog(cs2top/cs20)) / gamma
        f(:,:,n2,iss) = tmp_xy
        do i=1,nghost
          f(:,:,n2+i,iss) = 2*tmp_xy - f(:,:,n2-i,iss)
        enddo
      case default
        call fatal_error('bc_ss_temp_old','invalid argument')
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp_x: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_x: set x bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_x: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(l1,:,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(l1,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l1-i,:,:,iss) = -f(l1+i,:,:,iss) + tmp &
               - gamma_m1/gamma*(f(l1+i,:,:,ilnrho)+f(l1-i,:,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                       'bc_ss_temp_x: set x top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*, &
                       'bc_ss_temp_x: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(l2,:,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(l2,:,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(l2+i,:,:,iss) = -f(l2-i,:,:,iss) + tmp &
               - gamma_m1/gamma*(f(l2-i,:,:,ilnrho)+f(l2+i,:,:,ilnrho)-2*lnrho0)
        enddo
!
      case default
        call fatal_error('bc_ss_temp_x','invalid argument')
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp_y: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (ldebug) print*, &
                   'bc_ss_temp_y: set y bottom temperature - cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp_y: cannot have cs2bot<=0'
        tmp = 2/gamma*alog(cs2bot/cs20)
        f(:,m1,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,m1,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m1-i,:,iss) = -f(:,m1+i,:,iss) + tmp &
               - gamma_m1/gamma*(f(:,m1+i,:,ilnrho)+f(:,m1-i,:,ilnrho)-2*lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                     'bc_ss_temp_y: set y top temperature - cs2top=',cs2top
        if (cs2top<=0.) print*, &
                     'bc_ss_temp_y: cannot have cs2top<=0'
        tmp = 2/gamma*alog(cs2top/cs20)
        f(:,m2,:,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,m2,:,ilnrho)-lnrho0)
        do i=1,nghost
          f(:,m2+i,:,iss) = -f(:,m2-i,:,iss) + tmp &
               - gamma_m1/gamma*(f(:,m2-i,:,ilnrho)+f(:,m2+i,:,ilnrho)-2*lnrho0)
        enddo
!
      case default
        call fatal_error('bc_ss_temp_y','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp_y
!***********************************************************************
    subroutine bc_ss_temp_z(f,topbot,lone_sided)
!
!  boundary condition for entropy: constant temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!  13-sep-2016/axel: added TTbot, TTtop, to make this work
!
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      logical, optional :: lone_sided
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp_z: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
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
        if (ldebug) print*, &
                   'bc_ss_temp_z: set z bottom temperature: TTbot=',TTbot
        if (TTbot<=0.) print*, &
                   'bc_ss_temp_z: cannot have TTbot<=0'
        f(:,:,n1,ilnTT) = log(TTbot)
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                     'bc_ss_temp_z: set z top temperature: TTtop=',TTtop
        if (TTtop<=0.) print*,'bc_ss_temp_z: cannot have TTtop<=0'
        f(:,:,n2,ilnTT) = log(TTtop)
!
      case default
        call fatal_error('bc_ss_temp_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp_z
!***********************************************************************
    subroutine bc_lnrho_temp_z(f,topbot)
!
!  boundary condition for lnrho *and* ss: constant temperature
!
!  27-sep-2002/axel: coded
!  19-aug-2005/tobi: distributed across ionization modules
!
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'bc_lnrho_temp_z: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
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
        if (ldebug) print*, &
                 'bc_lnrho_temp_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0. .and. lroot) print*, &
                 'bc_lnrho_temp_z: cannot have cs2bot<=0'
        tmp = 2/gamma*log(cs2bot/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n1,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,:,n1,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n1-i,iss) = 2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2bot
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=-gravz/cs2bot
        do i=1,nghost
          f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) +f(:,:,n1+i,iss) &
                                                  -f(:,:,n1-i,iss) +dz2_bound(-i)*tmp
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                    'bc_lnrho_temp_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0. .and. lroot) print*, &
                    'bc_lnrho_temp_z: cannot have cs2top<=0'
        tmp = 2/gamma*log(cs2top/cs20)
!
!  set boundary value for entropy, then extrapolate ghost pts by antisymmetry
!
        f(:,:,n2,iss) = 0.5*tmp - gamma_m1/gamma*(f(:,:,n2,ilnrho)-lnrho0)
        do i=1,nghost; f(:,:,n2+i,iss) = 2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density in the ghost zones so that dlnrho/dz + ds/dz = gz/cs2top
!  for the time being, we don't worry about lnrho0 (assuming that it is 0)
!
        tmp=gravz/cs2top
        do i=1,nghost
          f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) +f(:,:,n2-i,iss) &
                                                  -f(:,:,n2+i,iss) +dz2_bound(i)*tmp
        enddo
!
      case default
        call fatal_error('bc_lnrho_temp_z','invalid argument')
      endselect
!
    endsubroutine bc_lnrho_temp_z
!***********************************************************************
    subroutine bc_lnrho_pressure_z(f,topbot)
!
!  boundary condition for lnrho: constant pressure
!
!   4-apr-2003/axel: coded
!   1-may-2003/axel: added the same for top boundary
!  19-aug-2005/tobi: distributed across ionization modules
!
      use Gravity, only: lnrho_bot,lnrho_top,ss_bot,ss_top
      use DensityMethods, only: putlnrho
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'bc_lnrho_pressure_z: cs20,cs0=',cs20,cs0
!
!  Constant pressure, i.e. antisymmetric
!  This assumes that the entropy is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is processor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('top')
        if (ldebug) print*,'bc_lnrho_pressure_z: lnrho_top,ss_top=',lnrho_top,ss_top
!
!  fix entropy if inflow (uz>0); otherwise leave s unchanged
!  afterwards set s antisymmetrically about boundary value
!
        if (lentropy) then
!         do m=m1,m2
!         do l=l1,l2
!           if (f(l,m,n1,iuz)>=0) then
!             f(l,m,n1,iss)=ss_bot
!           else
!             f(l,m,n1,iss)=f(l,m,n1+1,iss)
!           endif
!         enddo
!         enddo
          f(:,:,n2,iss)=ss_top
          do i=1,nghost; f(:,:,n2+i,iss)=2*f(:,:,n2,iss)-f(:,:,n2-i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n2,ilnrho)=lnrho_top+ss_top-f(:,:,n2,iss)
        else
          call putlnrho(f(:,:,n2,ilnrho),lnrho_top)
        endif
!
!  make density antisymmetric about boundary
!  another possibility might be to enforce hydrostatics
!  ie to set dlnrho/dz=-g/cs^2, assuming zero entropy gradient
!
        do i=1,nghost
          f(:,:,n2+i,ilnrho)=2*f(:,:,n2,ilnrho)-f(:,:,n2-i,ilnrho)
        enddo
!
!  top boundary
!
      case ('bot')
        if (ldebug) print*,'bc_lnrho_pressure_z: lnrho_bot,ss_bot=',lnrho_bot,ss_bot
!
!  fix entropy if inflow (uz>0); otherwise leave s unchanged
!  afterwards set s antisymmetrically about boundary value
!
        if (lentropy) then
!         do m=m1,m2
!         do l=l1,l2
!           if (f(l,m,n1,iuz)>=0) then
!             f(l,m,n1,iss)=ss_bot
!           else
!             f(l,m,n1,iss)=f(l,m,n1+1,iss)
!           endif
!         enddo
!         enddo
          f(:,:,n1,iss)=ss_bot
          do i=1,nghost; f(:,:,n1-i,iss)=2*f(:,:,n1,iss)-f(:,:,n1+i,iss); enddo
!
!  set density value such that pressure is constant at the bottom
!
          f(:,:,n1,ilnrho)=lnrho_bot+ss_bot-f(:,:,n1,iss)
        else
          call putlnrho(f(:,:,n1,ilnrho),lnrho_bot)
        endif
!
!  make density antisymmetric about boundary
!  another possibility might be to enforce hydrostatics
!  ie to set dlnrho/dz=-g/cs^2, assuming zero entropy gradient
!
        do i=1,nghost
          f(:,:,n1-i,ilnrho)=2*f(:,:,n1,ilnrho)-f(:,:,n1+i,ilnrho)
        enddo
!
      case default
        call fatal_error('bc_lnrho_pressure_z','invalid argument')
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real :: tmp
      integer :: i
!
      if (ldebug) print*,'bc_ss_temp2_z: cs20,cs0=',cs20,cs0
!
!  Constant temperature/sound speed for entropy, i.e. antisymmetric
!  ln(cs2) relative to cs2top/cs2bot.
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
        if (ldebug) print*, &
                   'bc_ss_temp2_z: set z bottom temperature: cs2bot=',cs2bot
        if (cs2bot<=0.) print*, &
                   'bc_ss_temp2_z: cannot have cs2bot<=0'
        tmp = 1/gamma*alog(cs2bot/cs20)
        do i=0,nghost
          f(:,:,n1-i,iss) = tmp &
               - gamma_m1/gamma*(f(:,:,n1-i,ilnrho)-lnrho0)
        enddo
!
!  top boundary
!
      case ('top')
        if (ldebug) print*, &
                     'bc_ss_temp2_z: set z top temperature: cs2top=',cs2top
        if (cs2top<=0.) print*,'bc_ss_temp2_z: cannot have cs2top<=0'
        tmp = 1/gamma*alog(cs2top/cs20)
        do i=0,nghost
          f(:,:,n2+i,iss) = tmp &
               - gamma_m1/gamma*(f(:,:,n2+i,ilnrho)-lnrho0)
        enddo
      case default
        call fatal_error('bc_ss_temp2_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_temp2_z
!***********************************************************************
    subroutine bc_ss_temp3_z(f,topbot)
!
!  31-jan-2013/axel: coded to impose cs2bot and dcs2bot at bottom
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call fatal_error('bc_ss_temp3_z', &
          'not implemented in eos_temperature_ionization.f90')
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_stemp_x: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (cs2bot<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2bot<=0'
        do i=1,nghost
          f(l1-i,:,:,iss) = f(l1+i,:,:,iss) &
               + gamma_m1/gamma*(f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho))
        enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0.) print*, &
                        'bc_ss_stemp_x: cannot have cs2top<=0'
        do i=1,nghost
          f(l2+i,:,:,iss) = f(l2-i,:,:,iss) &
               + gamma_m1/gamma*(f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho))
        enddo
!
      case default
        call fatal_error('bc_ss_stemp_x','invalid argument')
      endselect
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
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_stemp_y: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
      case ('bot')
        if (cs2bot<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2bot<=0'
        do i=1,nghost
          f(:,m1-i,:,iss) = f(:,m1+i,:,iss) &
               + gamma_m1/gamma*(f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho))
        enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0.) print*, &
                       'bc_ss_stemp_y: cannot have cs2top<=0'
        do i=1,nghost
          f(:,m2+i,:,iss) = f(:,m2-i,:,iss) &
               + gamma_m1/gamma*(f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho))
        enddo
!
      case default
        call fatal_error('bc_ss_stemp_y','invalid argument')
      endselect
!
    endsubroutine bc_ss_stemp_y
!***********************************************************************
    subroutine bc_ss_stemp_z(f,topbot)
!
!  boundary condition for entropy: symmetric temperature
!
!   3-aug-2002/wolf: coded
!  26-aug-2003/tony: distributed across ionization modules
!
      use Gravity
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_stemp_z: cs20,cs0=',cs20,cs0
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
          if (cs2bot<=0.) print*, &
                                  'bc_ss_stemp_z: cannot have cs2bot<=0'
          do i=1,nghost
             f(:,:,n1-i,iss) = f(:,:,n1+i,iss) &
                  + gamma_m1/gamma*(f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho))
          enddo
!
!  top boundary
!
      case ('top')
        if (cs2top<=0.) print*, &
                 'bc_ss_stemp_z: cannot have cs2top<=0'
         do i=1,nghost
           f(:,:,n2+i,iss) = f(:,:,n2-i,iss) &
                + gamma_m1/gamma*(f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho))
         enddo
      case default
        call fatal_error('bc_ss_stemp_z','invalid argument')
      endselect
!
    endsubroutine bc_ss_stemp_z
!***********************************************************************
    subroutine bc_ss_a2stemp_x(f,topbot)
!
!  boundary condition for entropy: adopt minimum value for entropy in ghost
!  zone, which satisfies either asymmetric temperature: vanishing 2nd deriv
!  or symmetric temperature to handle shock profiles.
!  Effectively caps temperature in ghost zones while otherwise fluctuating
!  with the outward flow.
!
!  22-sep-2010/fred: adapted from bc_ss_stemp_z
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_a2stemp_x: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
      select case (topbot)
!
!  bottom boundary
!
        case ('bot')
          if (cs2bot<=0.) print*, &
              'bc_ss_a2stemp_x: cannot have cs2bot<=0'
          do i=1,nghost
            f(l1-i,:,:,iss) = min( &
                2*f(l1+1-i,:,:,iss)-f(l1+2-i,:,:,iss)+gamma_m1/gamma* &
                (2*f(l1+1-i,:,:,ilnrho)-f(l1+2-i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)), &
                f(l1+i,:,:,iss)+gamma_m1/gamma* &
                (f(l1+i,:,:,ilnrho)-f(l1-i,:,:,ilnrho)))
          enddo
!
!  top boundary
!
        case ('top')
          if (cs2top<=0.) print*, &
              'bc_ss_a2stemp_x: cannot have cs2top<=0'
          do i=1,nghost
            f(l2+i,:,:,iss) = min( &
                2*f(l2-1+i,:,:,iss)-f(l2+2-i,:,:,iss)+gamma_m1/gamma* &
                (2*f(l2-1+i,:,:,ilnrho)-f(l2+2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho)), &
                f(l2+i,:,:,iss)+gamma_m1/gamma* &
                (f(l2-i,:,:,ilnrho)-f(l2+i,:,:,ilnrho)))
          enddo
!
        case default
          call fatal_error('bc_ss_a2stemp_x','invalid argument')
      endselect
!
    endsubroutine bc_ss_a2stemp_x
!***********************************************************************
    subroutine bc_ss_a2stemp_y(f,topbot)
!
!  boundary condition for entropy: adopt minimum value for entropy in ghost
!  zone, which satisfies either asymmetric temperature: vanishing 2nd deriv
!  or symmetric temperature to handle shock profiles.
!  Effectively caps temperature in ghost zones while otherwise fluctuating
!  with the outward flow.
!
!  22-sep-2010/fred: adapted from bc_ss_stemp_y
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_a2stemp_y: cs20,cs0=',cs20,cs0
!
!  Symmetric temperature/sound speed for entropy.
!  This assumes that the density is already set (ie density _must_ register
!  first!)
!
!  check whether we want to do top or bottom (this is precessor dependent)
!
!
      select case (topbot)
!
!  bottom boundary
!
        case ('bot')
          if (cs2bot<=0.) print*, &
              'bc_ss_a2stemp_y: cannot have cs2bot<=0'
          do i=1,nghost
            f(:,m1-i,:,iss) = min( &
                2*f(:,m1+1-i,:,iss)-f(:,m1+2-i,:,iss)+gamma_m1/gamma* &
                (2*f(:,m1+1-i,:,ilnrho)-f(:,m1+2-i,:,ilnrho)-f(:,m1-i,:,ilnrho)), &
                f(:,m1-i,:,iss)+gamma_m1/gamma* &
                (f(:,m1+i,:,ilnrho)-f(:,m1-i,:,ilnrho)))
          enddo
!
!  top boundary
!
        case ('top')
          if (cs2top<=0.) print*, &
              'bc_ss_a2stemp_y: cannot have cs2top<=0'
          do i=1,nghost
            f(:,m2+i,:,iss) = min( &
                2*f(:,m2-1+i,:,iss)-f(:,m2-2+i,:,iss)+gamma_m1/gamma* &
                (2*f(:,m2-1+i,:,ilnrho)-f(:,m2-2+i,:,ilnrho)-f(:,m2+i,:,ilnrho)), &
                f(:,m2-i,:,iss)+gamma_m1/gamma* &
                (f(:,m2-i,:,ilnrho)-f(:,m2+i,:,ilnrho)))
          enddo
!
        case default
          call fatal_error('bc_ss_a2stemp_y','invalid argument')
      endselect
!
    endsubroutine bc_ss_a2stemp_y
!***********************************************************************
    subroutine bc_ss_a2stemp_z(f,topbot)
!
!  boundary condition for entropy: adopt minimum value for entropy in ghost
!  zone, which satisfies either asymmetric temperature: vanishing 2nd deriv
!  or symmetric temperature to handle shock profiles.
!  Effectively caps temperature in ghost zones while otherwise fluctuating
!  with the outward flow.
!
!  22-sep-2010/fred: adapted from bc_ss_stemp_z
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      integer :: i
!
      if (ldebug) print*,'bc_ss_a2stemp_z: cs20,cs0=',cs20,cs0
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
          if (cs2bot<=0.) print*, &
              'bc_ss_a2stemp_z: cannot have cs2bot<=0'
          do i=1,nghost
            f(:,:,n1-i,iss) = min( &
                2*f(:,:,n1+1-i,iss)-f(:,:,n1+2-i,iss) + gamma_m1/gamma* &
                (2*f(:,:,n1+1-i,ilnrho)-f(:,:,n1+2-i,ilnrho)-f(:,:,n1-i,ilnrho)), &
                f(:,:,n1+i,iss)+gamma_m1/gamma* &
                (f(:,:,n1+i,ilnrho)-f(:,:,n1-i,ilnrho)))
          enddo
!
!  top boundary
!
        case ('top')
          if (cs2top<=0.) print*, &
              'bc_ss_a2stemp_z: cannot have cs2top<=0'
          do i=1,nghost
            f(:,:,n2+i,iss) = min( &
                2*f(:,:,n2-1+i,iss)-f(:,:,n2-2+i,iss) + gamma_m1/gamma* &
                (2*f(:,:,n2-1+i,ilnrho)-f(:,:,n2-2+i,ilnrho)-f(:,:,n2+i,ilnrho)), &
                f(:,:,n2-i,iss)+gamma_m1/gamma* &
                (f(:,:,n2-i,ilnrho)-f(:,:,n2+i,ilnrho)))
          enddo
        case default
          call fatal_error('bc_ss_a2stemp_z','invalid argument')
      endselect
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
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
      real, dimension (size(f,1),size(f,2)) :: cs2_2d
      integer :: i
!
!  The 'ce' boundary condition for entropy makes the energy constant at
!  the boundaries.
!  This assumes that the density is already set (ie density must register
!  first!)
!
    select case (topbot)
!
! Bottom boundary
!
    case ('bot')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma_m1*f(:,:,n1,ilnrho)+gamma*f(:,:,n1,iss))
      do i=1,nghost
         f(:,:,n1-i,iss)=1./gamma*(-gamma_m1*f(:,:,n1-i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
!
! Top boundary
!
    case ('top')
      !  Set cs2 (temperature) in the ghost points to the value on
      !  the boundary
      !
      cs2_2d=cs20*exp(gamma_m1*f(:,:,n2,ilnrho)+gamma*f(:,:,n2,iss))
      do i=1,nghost
         f(:,:,n2+i,iss)=1./gamma*(-gamma_m1*f(:,:,n2+i,ilnrho)-log(cs20)&
              +log(cs2_2d))
      enddo
    case default
      call fatal_error('bc_ss_energy','invalid argument')
    endselect
!
    endsubroutine bc_ss_energy
!***********************************************************************
    subroutine bc_stellar_surface(f,topbot)
!
!  Boundary condition for density.
!
!  We make both the lower and the upper boundary hydrostatic.
!    rho1*(dpp/dz) = gravz
!
!  Additionally, we make the lower boundary isentropic
!    (dss/dz) = 0
!  and the upper boundary isothermal
!    (dlnTT/dz) = 0
!
!  At the moment, this is probably only useful for solar surface convection.
!
!
!  11-May-2006/tobi: coded
!  16-May-2006/tobi: isentropic lower boundary
!
      use Gravity, only: gravz,gravz_profile,reduced_top
      use DensityMethods, only: getlnrho
!
      real, dimension (:,:,:,:), intent (inout) :: f
      character (len=3), intent (in) :: topbot
!
      real, dimension (size(f,1),size(f,2)) :: lnrho,lnTT,TT1
      real, dimension (size(f,1),size(f,2)) :: rhs,sqrtrhs,yH
      real, dimension (size(f,1),size(f,2)) :: mu1,rho1pp
      real, dimension (size(f,1),size(f,2)) :: yH_term_cv,yH_term_cp
      real, dimension (size(f,1),size(f,2)) :: TT_term_cv,TT_term_cp
      real, dimension (size(f,1),size(f,2)) :: alpha,delta
      real, dimension (size(f,1),size(f,2)) :: cv,cp,cs2,nabla_ad
      real, dimension (size(f,1),size(f,2)) :: dlnrhodz,dlnTTdz
      real :: fac
      integer :: i
!
      select case (topbot)
!
!  Bottom boundary
!
      case ('bot')
!
!  Boundary condition for density and temperature
!
        if (bcz12(ilnTT,1)/='StS'.and.bcz12(ilnTT,1)/='') then
          call fatal_error("bc_stellar_surface", &
                           "This boundary condition for density also sets "// &
                           "temperature. We therfore require "// &
                           "bcz12(ilnTT,1)='StS' or bcz12(ilnTT,1)=''")
        endif
!
!  Get variables from f-array
!
        call getlnrho(f(:,:,n1,ilnrho),lnrho)
        lnTT = f(:,:,n1,ilnTT)
        TT1 = exp(-lnTT)
!
!  Hydrogen ionization fraction
!
        rhs = exp(lnrho_e-lnrho + 1.5*(lnTT-lnTT_ion) - TT_ion*TT1)
        sqrtrhs = sqrt(rhs)
        yH = 2*sqrtrhs/(sqrtrhs+sqrt(4+rhs))
!
!  Inverse mean molecular weight
!
        mu1 = mu1_0*(1 + yH + xHe)
!
!  Pressure over density
!
        rho1pp = Rgas*mu1/TT1
!
!  Abreviations
!
        yH_term_cv = yH*(1-yH)/((2-yH)*(1+yH+xHe))
        TT_term_cv = 1.5 + TT_ion*TT1
!
        yH_term_cp = yH*(1-yH)/(2+xHe*(2-yH))
        TT_term_cp = 2.5 + TT_ion*TT1
!
!  Specific heats in units of Rgas/mu
!
        cv = 1.5 + yH_term_cv*TT_term_cv**2
        cp = 2.5 + yH_term_cp*TT_term_cp**2
!
!  See Kippenhahn & Weigert
!
        alpha = ((2-yH)*(1+yH+xHe))/(2+xHe*(2-yH))
        delta = 1 + yH_term_cp*TT_term_cp
!
!  Speed of sound
!
        cs2 = cp*rho1pp/(alpha*cv)
!
!  Adiabatic pressure gradient
!
        nabla_ad = delta/cp
!
!  z-derivatives of density and temperature on the boundary
!
        dlnrhodz = gravz/cs2
        dlnTTdz = (nabla_ad/rho1pp)*gravz
!
!  Fill ghost zones accordingly
!
        do i=1,nghost
          f(:,:,n1-i,ilnrho) = f(:,:,n1+i,ilnrho) - dz2_bound(-i)*dlnrhodz
          f(:,:,n1-i,ilnTT)  = f(:,:,n1+i,ilnTT)  - dz2_bound(-i)*dlnTTdz
        enddo
!
!  Top boundary
!
      case ('top')
!
!  Boundary condition for density, temperature, and vector potential
!
        if (bcz12(ilnTT,2)/='StS'.and.bcz12(ilnTT,2)/='') then
          call fatal_error("bc_stellar_surface", &
                           "This boundary condition for density also sets "// &
                           "temperature. We therfore require "// &
                           "bcz2(ilnTT)='StS' or bcz2(ilnTT)=''")
        endif
!
!  Get variables from f-array
!
        call getlnrho(f(:,:,n2,ilnrho),lnrho)
        lnTT = f(:,:,n2,ilnTT)
        TT1 = exp(-lnTT)
!
!  `Effective' gravitational acceleration (geff = gravz - rho^-1*dz1ppm)
!
        if (gravz_profile=='reduced_top') then
          fac = reduced_top
        else
          fac = 1.0
        endif
!
!  Hydrogen ionization fraction
!
        rhs = exp(lnrho_e-lnrho + 1.5*(lnTT-lnTT_ion) - TT_ion*TT1)
        sqrtrhs = sqrt(rhs)
        yH = 2*sqrtrhs/(sqrtrhs+sqrt(4+rhs))
!
!  Inverse mean molecular weight
!
        mu1 = mu1_0*(1 + yH + xHe)
!
!  Pressure over density
!
        rho1pp = Rgas*mu1/TT1
!
!  Pressure derivative
!
        alpha = ((2-yH)*(1+yH+xHe))/(2+xHe*(2-yH))
!
!  z-derivatives of density on the boundary
!
        dlnrhodz = fac*gravz*alpha/rho1pp
!
!  Fill ghost zones accordingly
!
        do i=1,nghost
          f(:,:,n2+i,ilnrho) = f(:,:,n2-i,ilnrho) + dz2_bound(i)*dlnrhodz
          f(:,:,n2+i,ilnTT) = f(:,:,n2-i,ilnTT)
        enddo
!
      case default
!
      endselect
!
    endsubroutine bc_stellar_surface
!***********************************************************************
    subroutine bc_lnrho_cfb_r_iso(f,topbot)
!
      use Mpicomm, only: stop_it
!
      character (len=3) :: topbot
      real, dimension (:,:,:,:) :: f
!
      call stop_it("bc_lnrho_cfb_r_iso: NOT IMPLEMENTED in "// &
          "eos_temperature_ionization")
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
      real, dimension (:,:,:,:) :: f
!
      call stop_it("bc_lnrho_hds_z_iso: NOT IMPLEMENTED in "// &
          "eos_temperature_ionization")
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
      real, dimension (:,:,:,:) :: f
!
      call stop_it("bc_lnrho_hdss_z_iso: NOT IMPLEMENTED in "// &
          "eos_temperature_ionization")
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(topbot)
!
    endsubroutine bc_lnrho_hdss_z_iso
!***********************************************************************
    subroutine write_thermodyn
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
    subroutine get_stratz(z, rho0z, dlnrho0dz, eth0z)
!
!  Get background stratification in z direction.
!
!  13-oct-14/ccyang: dummy
!
      real, dimension(:), intent(in) :: z
      real, dimension(:), intent(out), optional :: rho0z, dlnrho0dz, eth0z
!
      call fatal_error('get_stratz', 'Stratification for this EOS is not implemented. ')
!
      call keep_compiler_quiet(z)
      if (present(rho0z)) call keep_compiler_quiet(rho0z)
      if (present(dlnrho0dz)) call keep_compiler_quiet(dlnrho0dz)
      if (present(eth0z)) call keep_compiler_quiet(eth0z)
!
    endsubroutine get_stratz
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr

    integer, parameter :: n_pars=1
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(cs20,p_par(1))

    endsubroutine pushpars2c
!***********************************************************************
endmodule EquationOfState
