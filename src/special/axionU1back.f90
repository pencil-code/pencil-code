! $Id$
!
!  Solve axion-U(1) inflation model for many wavenumber values.
!  The different wavenumber correspond to 1D "meshpoints".
!
!  12-sep-25/ramkishor: adapted from axionSU2back.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 8
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
  include '../record_types.h'
!
!
! Declare index of variables
!
  integer :: iaxi_AL=0, iaxi_ALdot=0, iaxi_AR=0, iaxi_ARdot=0
  integer :: iaxi_imAL=0, iaxi_imALdot=0, iaxi_imAR=0, iaxi_imARdot=0
  integer :: iaxi_lna=0, iaxi_phi=0, iaxi_phidot=0, Ndivt=40, ikdata=100
!
  ! input parameters
  real :: a, k0=1.06e-4, dk=1.06e-4, ascale_ini=1.
  real :: alpf=50.
  real :: H=1.04e-6
  real :: H_init
  real :: mpl2=1., Hdot=0.
  real :: m_phi=1.06e-6, phi_ini=0.7, dphi_ini=-1.687e-7
  real :: alpha=0.1, m_alpha=3.285e-11, n_alpha=1.5
  real :: phi_0=0.014, lambda_gmssm=1.0, alpha1_gmssm=6.5e-11, const_sigma=0.0

  real :: lnkmin0=-1., lnkmin0_dummy, lnkmax0, dlnk
  real :: nmin0=-1., nmax0=3., horizon_factor=0., sgn=1.
  real :: edotb_sum, rhob, rhoe
  real, dimension (nx) :: dt1_special, lnk
  real, dimension (nx) :: edotb, e2, b2
  logical :: lbackreact=.false., lhubble=.true., lquant_filter=.true.
  logical :: ldo_adjust_krange=.false., lswap_sign=.false.
  logical :: lwrite_krange=.true., lwrite_backreact=.true.
  logical :: llnk_spacing_adjustable=.false., llnk_spacing=.false.
  character(len=50) :: init_axionU1back='standard'
  character (len=labellen) :: V_choice='quadratic'

! initial parameters
  namelist /special_init_pars/ &
    k0, dk, alpf, H, Ndivt, llnk_spacing_adjustable, llnk_spacing, &
    nmin0, nmax0, ldo_adjust_krange, lswap_sign, m_phi, &
    lhubble, V_choice, phi_ini, dphi_ini, alpha, m_alpha, n_alpha, mpl2, &
    const_sigma, phi_0, alpha1_gmssm, lambda_gmssm, lquant_filter

  ! run parameters
  namelist /special_run_pars/ &
    k0, dk, alpf, H, lbackreact, Ndivt, llnk_spacing_adjustable, llnk_spacing, &
    nmin0, nmax0, ldo_adjust_krange, lswap_sign, lwrite_krange, lwrite_backreact, &
    lhubble, ikdata, const_sigma, lquant_filter
!
  ! k array
  real, dimension (nx) :: k
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_a =0 ! DIAG_DOC: $a$
  integer :: idiag_phi   =0 ! DIAG_DOC: $phi$
  integer :: idiag_phidot   =0 ! DIAG_DOC: $phidot$
  integer :: idiag_H =0 ! DIAG_DOC: $Hubble_parameter$
  integer :: idiag_AL =0 ! DIAG_DOC: $\dot\psi$
  integer :: idiag_ALdot=0 ! DIAG_DOC: $\ddot\psi$
  integer :: idiag_AR     =0 ! DIAG_DOC: $T_R$
  integer :: idiag_ARdot     =0 ! DIAG_DOC: $T_L$
  integer :: idiag_imAR=0 ! DIAG_DOC: $\Im T_R$
  integer :: idiag_imAL  =0 ! DIAG_DOC: $T_R^{\rm anal}$
  integer :: idiag_imARdot=0 ! DIAG_DOC: $\Im T_R$
  integer :: idiag_imALdot  =0
  integer :: idiag_k0=0   ! DIAG_DOC: $k0$
  integer :: idiag_dk=0   ! DIAG_DOC: $dk$
  integer :: idiag_rhoe     =0 ! DIAG_DOC: $rho_e$
  integer :: idiag_rhob     =0 ! DIAG_DOC: $rho_B$
  integer :: idiag_edotb     =0 ! DIAG_DOC: $edotb$

!
! z averaged diagnostics given in zaver.in
!
  integer :: enum_v_choice = 0
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialized (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Set iaxionU1back to consecutive numbers
!
      call farray_register_pde('axi_AR'     ,iaxi_AR)
      call farray_register_pde('axi_ARdot'  ,iaxi_ARdot)
      call farray_register_pde('axi_AL'   ,iaxi_AL)
      call farray_register_pde('axi_ALdot',iaxi_ALdot)
!
!  Possibility to evolve the Hubble parameter (in cosmic time)
!
      if (lhubble) then
        call farray_register_ode('axi_lna'   ,iaxi_lna)
        call farray_register_ode('axi_phi'   ,iaxi_phi)
        call farray_register_ode('axi_phidot',iaxi_phidot)
      endif
!
!
!  The imaginary parts.
!
      call farray_register_pde('axi_imAR'   ,iaxi_imAR)
      call farray_register_pde('axi_imARdot',iaxi_imARdot)
      call farray_register_pde('axi_imAL'    ,iaxi_imAL)
      call farray_register_pde('axi_imALdot' ,iaxi_imALdot)
!
!
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: lnH, lna, a
      real :: kmax=2., lnkmax, lnk0=1.
      real :: phi, phidot, tmp
      real :: V, beta, alpha_gmssm
      integer :: ik
!
!  Initialize any module variables which are parameter dependent
!
!
!  Compute lnkmin0 and lnkmax0. Even for a linear k-range, dlnk
!  is needed to determine the output of k-range and grand etc.
!
      a=ascale_ini
      if (lhubble) then
        phidot=dphi_ini
        select case (V_choice)
          case ('alpha_attractors')
            beta=sqrt(2./(3.*alpha))
            V=alpha*m_alpha*(tanh(beta*phi_ini/2)**2)**n_alpha
          case ('quadratic') ; V=.5*(m_phi*phi_ini)**2
          case ('gmssm')
            tmp=phi_ini/phi_0
            alpha_gmssm=1-alpha1_gmssm
            V=lambda_gmssm**4*(0.5*tmp**2-(alpha_gmssm*onethird)*tmp**6+&
              (alpha_gmssm/10.)*tmp**10)
          case default
            call fatal_error("initialize_special: No such V_choice: ", trim(V_choice))
        endselect
        H=sqrt(8.*pi*onethird*(1/mpl2)*(.5*phidot**2+V))
      endif

      lna=alog(a)
      lnH=alog(H)
      if (lstart) then
        lnkmin0=nmin0+lnH+lna
        lnkmax0=nmax0+lnH+lna
      endif
      if (nxgrid==1) then
        dlnk=1.
      else
        dlnk=(nmax0-nmin0)/(nxgrid-1)
      endif
!
!  Initialize lnkmin0 and lnkmax0
!
      if (llnk_spacing_adjustable) then
        if (ldo_adjust_krange) then
          do ik=1,nx
            lnk(ik)=lnkmin0+dlnk*(ik-1+ipx*nx)
            k(ik)=exp(lnk(ik))
          enddo
        else
          do ik=1,nx
            lnk(ik)=nmin0+lnH+lna+dlnk*(ik-1+ipx*nx)
            k(ik)=exp(lnk(ik))
          enddo
          if (ip<10) print*,'iproc,lnk=',iproc,lnk
        endif
      elseif (llnk_spacing) then
        do ik=1,nx
          lnk(ik)=nmin0+lnH+lna+dlnk*(ik-1+ipx*nx)
          k(ik)=exp(lnk(ik))
        enddo
        lnkmin0_dummy=lnkmin0
        if (ip<10) print*,'iproc,lnk=',iproc,lnk
      else
        do ik=1,nx
          k(ik)=k0+dk*(ik-1+ipx*nx)
        enddo
        lnk=impossible
        lnkmin0_dummy=nmin0+lnH+lna
      endif
!
      call keep_compiler_quiet(f)
      H_init = H
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: AL, ALdot, AR, ARdot
      real, dimension (nx) :: imAL, imALdot, imAR, imARdot
      real :: V, Vprime, lnt, lnH, lna, a, beta, alpha_gmssm, tmp
      real :: kmax=2., lnkmax, lnk0=1.
      integer :: ik
!
      intent(inout) :: f
!
!  Initial condition; depends on k, which is here set to x.
!
      select case (init_axionU1back)
        case ('nothing'); if (lroot) print*,'nothing'
        case ('standard')
          a=ascale_ini
          if (lhubble) then
            f_ode(iaxi_lna)=alog(ascale_ini)
            f_ode(iaxi_phi)=phi_ini
            f_ode(iaxi_phidot)=dphi_ini
            select case (V_choice)
              case ('alpha_attractors')
                beta=sqrt(2./(3.*alpha))
                V=alpha*m_alpha*(tanh(beta*phi_ini/2)**2)**n_alpha
              case ('quadratic')
                V=.5*(m_phi*phi_ini)**2
                Vprime=m_phi**2*phi_ini
              case ('gmssm')
                tmp=phi_ini/phi_0
                alpha_gmssm=1-alpha1_gmssm
                V=lambda_gmssm**4*(0.5*tmp**2-(alpha_gmssm*onethird)*tmp**6+&
                  (alpha_gmssm/10.)*tmp**10)
              case default
                call fatal_error("init_special: No such V_choice: ", trim(V_choice))
            endselect
            H=sqrt(8.*pi*onethird*(1/mpl2)*(.5*f_ode(iaxi_phidot)**2+V))
          else
            a=exp(H*t)
          endif
!
!  need k
!
          AL=(ascale_ini/sqrt(2.*k))*cos(k/(ascale_ini*H))
          AR=(ascale_ini/sqrt(2.*k))*cos(k/(ascale_ini*H))
          ALdot=(k/sqrt(2.*k))*sin(k/(ascale_ini*H))
          ARdot=(k/sqrt(2.*k))*sin(k/(ascale_ini*H))
          imAL=(ascale_ini/sqrt(2.*k))*sin(k/(ascale_ini*H))
          imAR=(ascale_ini/sqrt(2.*k))*sin(k/(ascale_ini*H))
          imALdot=(-k/sqrt(2.*k))*cos(k/(ascale_ini*H))
          imARdot=(-k/sqrt(2.*k))*cos(k/(ascale_ini*H))
!
!  PDE variables
!
          do n=n1,n2
          do m=m1,m2
            f(l1:l2,m,n,iaxi_AL)=AL
            f(l1:l2,m,n,iaxi_ALdot)=ALdot
            f(l1:l2,m,n,iaxi_AR)=AR
            f(l1:l2,m,n,iaxi_ARdot)=ARdot
            f(l1:l2,m,n,iaxi_imAL)=imAL
            f(l1:l2,m,n,iaxi_imALdot)=imALdot
            f(l1:l2,m,n,iaxi_imAR)=imAR
            f(l1:l2,m,n,iaxi_imARdot)=imARdot
          enddo
          enddo
        case default
          call fatal_error("init_special","no such init_axionU1back: "//trim(init_axionU1back))
      endselect
!
1000  format(a,2x,i3,1p,80e14.6)
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
 
!  All pencils that this special module depends on are specified here.
!
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Most basic pencils should come first, as others may depend on them.
!
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine get_Hubble_and_scale_factor(f_ode)

      real, dimension(n_odevars), intent(IN) :: f_ode
      real :: phi, phidot, V, beta, Vprime, phiddot, tmp, alpha_gmssm
!
      if (lhubble) then
        phi=f_ode(iaxi_phi)
        phidot=f_ode(iaxi_phidot)
        select case (V_choice)
          case ('alpha_attractors')
            beta=sqrt(2./(3.*alpha))
            V=alpha*m_alpha*(tanh(beta*phi/2)**2)**n_alpha
            Vprime=alpha*m_alpha*n_alpha*beta*tanh(beta*phi/2)*(1./cosh(beta*phi/2)**2)
            case ('quadratic')
              V=.5*m_phi**2*phi**2
              Vprime=m_phi**2*phi
            case ('gmssm')
              tmp=phi/phi_0
              alpha_gmssm=1-alpha1_gmssm
              V=lambda_gmssm**4*(0.5*tmp**2-(alpha_gmssm*onethird)*tmp**6+&
                (alpha_gmssm/10.)*tmp**10)
              Vprime=lambda_gmssm**4*(tmp-2*alpha_gmssm*tmp**5+&
                alpha_gmssm*tmp**9)/phi_0
          case default
            call fatal_error("dspecial_dt: No such V_choice: ", trim(V_choice))
        endselect
        a=exp(f_ode(iaxi_lna))
        if (lbackreact) then
          H=sqrt(8.*pi*onethird*(1/mpl2)*(.5*phidot**2+V+rhoe+rhob))
        else
          H=sqrt(8.*pi*onethird*(1/mpl2)*(.5*phidot**2+V))
        endif
        phiddot=-3.*H*phidot-Vprime
      else
        H=H_init
        a=1.0
      endif
    endsubroutine get_Hubble_and_scale_factor
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  calculate right hand side of ONE OR MORE extra coupled PDEs
!  along the 'current' Pencil, i.e. f(l1:l2,m,n) where
!  m,n are global variables looped over in equ.f90
!
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df array.  NEVER reset it to zero.
!
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!
      use General, only: random_number_wrapper
      use Diagnostics
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: AL, ALdot, AR, ARdot, ALddot, ARddot
      real, dimension (nx) :: imAL , imALdot, imAR, imARdot, imALddot, imARddot
      real :: phi, phidot
      integer :: ik
      type (pencil_case) :: p
!
      intent(in) :: p
      intent(inout) :: f
      intent(inout) :: df
!
      call keep_compiler_quiet(p)
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
!  Set the 8 variable
!
      AL=f(l1:l2,m,n,iaxi_AL)
      ALdot=f(l1:l2,m,n,iaxi_ALdot)
      AR=f(l1:l2,m,n,iaxi_AR)
      ARdot=f(l1:l2,m,n,iaxi_ARdot)
      imAL=f(l1:l2,m,n,iaxi_imAL)
      imALdot=f(l1:l2,m,n,iaxi_imALdot)
      imAR=f(l1:l2,m,n,iaxi_imAR)
      imARdot=f(l1:l2,m,n,iaxi_imARdot)
!
      call get_Hubble_and_scale_factor(f_ode)
      phidot=f_ode(iaxi_phidot)
!
!
!  same with cosmic time; should also write in terms of psiddot and TRddot
!  But this is not currently done because of mixed terms.
!
      ALddot=-(const_sigma*a**(0.68)+H)*ALdot-(k/a)*(k/a-alpf*phidot)*AL
      ARddot=-(const_sigma*a**(0.68)+H)*ARdot-(k/a)*(k/a+alpf*phidot)*AR
      df(l1:l2,m,n,iaxi_AL)=df(l1:l2,m,n,iaxi_AL)+ALdot
      df(l1:l2,m,n,iaxi_AR )=df(l1:l2,m,n,iaxi_AR )+ARdot
      df(l1:l2,m,n,iaxi_ALdot)=df(l1:l2,m,n,iaxi_ALdot)+ALddot
      df(l1:l2,m,n,iaxi_ARdot)=df(l1:l2,m,n,iaxi_ARdot)+ARddot
      imALddot=-(const_sigma*a**(0.68)+H)*imALdot-(k/a)*(k/a-alpf*phidot)*imAL
      imARddot=-(const_sigma*a**(0.68)+H)*imARdot-(k/a)*(k/a+alpf*phidot)*imAR
      df(l1:l2,m,n,iaxi_imAL   )=df(l1:l2,m,n,iaxi_imAL   )+imALdot
      df(l1:l2,m,n,iaxi_imALdot)=df(l1:l2,m,n,iaxi_imALdot)+imALddot
      df(l1:l2,m,n,iaxi_imAR    )=df(l1:l2,m,n,iaxi_imAR    )+imARdot
      df(l1:l2,m,n,iaxi_imARdot )=df(l1:l2,m,n,iaxi_imARdot )+imARddot
!
!
!
      if (lfirst.and.ldt) then
        if (V_choice=='gmssm') then
          dt1_special = Ndivt*max(lambda_gmssm**2/phi_0, maxval(k)/a)
        else
          dt1_special = Ndivt*max(m_phi, maxval(k)/a)
        endif
        dt1_max=max(dt1_max,dt1_special)
      endif
!
!  diagnostics
!
      call calc_diagnostics_special(f,p)
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine calc_diagnostics_special(f,p)
!
!
      use Diagnostics

      real, dimension (nx) :: psi_anal, psidot_anal, TR_anal, TRdot_anal
      real, dimension(mx,my,mz,mfarray) :: f
      real, parameter :: fact=1.
      type(pencil_case) :: p
      if (ldiagnos) then
        call sum_mn_name(sqrt(f(l1:l2,m,n,iaxi_AL)**2),idiag_AL)
        call sum_mn_name(sqrt(f(l1:l2,m,n,iaxi_ALdot)**2),idiag_ALdot)
        call sum_mn_name(sqrt(f(l1:l2,m,n,iaxi_AR)**2),idiag_AR)
        call sum_mn_name(sqrt(f(l1:l2,m,n,iaxi_ARdot)**2),idiag_ARdot)
        call sum_mn_name(sqrt(f(l1:l2,m,n,iaxi_imAL)**2),idiag_imAL)
        call sum_mn_name(sqrt(f(l1:l2,m,n,iaxi_imAR)**2),idiag_imAR)
        call sum_mn_name(sqrt(f(l1:l2,m,n,iaxi_imALdot)**2),idiag_imALdot)
        call sum_mn_name(sqrt(f(l1:l2,m,n,iaxi_imARdot)**2),idiag_imARdot)
        call save_name(k0,idiag_k0)
        call save_name(dk,idiag_dk)
        call save_name(rhoe, idiag_rhoe)
        call save_name(rhob, idiag_rhob)
        call save_name(edotb_sum, idiag_edotb)
      endif
!
!
    endsubroutine calc_diagnostics_special
!***********************************************************************
    subroutine calc_ode_dt(f_ode,phi,phidot,phiddot,H)
!
      use General, only: random_number_wrapper
      use Mpicomm
      use Sub

      real, dimension(:), intent(IN) :: f_ode
      real, intent(OUT) :: phiddot, phidot, H, phi
      real :: beta, V, Vprime, alpha_gmssm, tmp
      real :: fact=1., sign_swap=1.
      
!
!  Set the all variable
!
      if (iaxi_phi > 0) phi = f_ode(iaxi_phi)
      if (iaxi_phidot > 0) phidot = f_ode(iaxi_phidot)
!
      call get_Hubble_and_scale_factor(f_ode)
      if (lhubble) then
        select case (V_choice)
          case ('alpha_attractors')
            beta=sqrt(2./(3.*alpha))
            V=alpha*m_alpha*(tanh(beta*phi/2)**2)**n_alpha
            Vprime=alpha*m_alpha*n_alpha*beta*tanh(beta*phi/2)*(1./cosh(beta*phi/2)**2)
          case ('quadratic')
            V=.5*m_phi**2*phi**2
            Vprime=m_phi**2*phi
          case ('gmssm')
            tmp=phi/phi_0
            alpha_gmssm=1-alpha1_gmssm
            V=lambda_gmssm**4*(0.5*tmp**2-(alpha_gmssm*onethird)*tmp**6+&
              (alpha_gmssm/10.)*tmp**10)
            Vprime=lambda_gmssm**4*(tmp-2*alpha_gmssm*tmp**5+&
                   alpha_gmssm*tmp**9)/phi_0
          case default
            call fatal_error("special_dspecial_dt_ode: No such V_choice: ", trim(V_choice))
        endselect
        if (lbackreact) then
          H=sqrt(8.*pi*onethird*(1/mpl2)*(.5*phidot**2+V+rhoe+rhob))
        else
          H=sqrt(8.*pi*onethird*(1/mpl2)*(.5*phidot**2+V))
        endif
        phiddot=-3.*H*phidot-Vprime
      endif
      if (lbackreact) then
        phiddot  =phiddot + alpf*edotb_sum/a**4
      endif
    endsubroutine calc_ode_dt
!***********************************************************************
    subroutine dspecial_dt_ode
!
!  Calculate right hand side(s) of one or more extra coupled ODEs.
!  Due to the multi-step Runge Kutta timestepping used one MUST always
!  add to the present contents of the df_ode array.  NEVER reset it to zero.
!
      real ::  phiddot, phidot
!
!  identify module and boundary conditions
!
      if (lgpu) then
        call read_sums_from_GPU
      endif
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
      call calc_ode_dt(f_ode,phi,phidot,phiddot,H)
!
!  Choice whether or not we want to update the background
!
!  Possibility to evolve the Hubble parameter (in cosmic time)
!
      if (lhubble) then
        df_ode(iaxi_lna)=df_ode(iaxi_lna)+H
        df_ode(iaxi_phi)=df_ode(iaxi_phi)+phidot
        df_ode(iaxi_phidot)=df_ode(iaxi_phidot)+phiddot
      endif
!
!  diagnostics
!
      if (ldiagnos .and. .not. lgpu) then
        call calc_ode_diagnostics_special(f_ode)
      endif
!
    endsubroutine dspecial_dt_ode
!***********************************************************************
    subroutine  calc_ode_diagnostics_special(f_ode)
!
      use Diagnostics
      real, dimension(n_odevars), intent(IN) :: f_ode
      real :: phiddot

      if (lhubble) then
        call get_Hubble_and_scale_factor(f_ode)
        call save_name(a     ,idiag_a)
        call save_name(f_ode(iaxi_phi)   ,idiag_phi)
        call save_name(f_ode(iaxi_phidot),idiag_phidot)
        call save_name(H     ,idiag_H)
      endif
    endsubroutine  calc_ode_diagnostics_special 
!***********************************************************************
    subroutine read_special_init_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(LEN=*), intent(out) :: iomsg
      integer :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
      !call keep_compiler_quiet(f)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iomsg)
!
      use File_io, only: parallel_unit
!
      character(LEN=*), intent(out) :: iomsg
      integer :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat, IOMSG=iomsg)
      if (iostat==0) iomsg=""
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine input_persist_special_id(id,done)
!
      use IO, only: read_persist
!
      integer :: id
      logical :: done
!
!      print*,'ram_persist',id
      select case (id)
        case (id_record_SPECIAL_LNKMIN0)
          done=read_persist ('SPECIAL_LNKMIN0', lnkmin0)
          if (lroot .and. .not. done) print *, 'input_persist_special: ', lnkmin0
      endselect
!
    endsubroutine input_persist_special_id
!***********************************************************************
    subroutine input_persist_special
!
!  Read in the persistent special variables.
!
      use IO, only: read_persist
!
      logical :: error
!
      error = read_persist ('SPECIAL_LNKMIN0', lnkmin0)
      if (lroot .and. .not. error) print *, 'input_persist_special: lnkmin0: ', lnkmin0
!
    endsubroutine input_persist_special
!***********************************************************************
    logical function output_persistent_special()
!
      use IO, only: write_persist
!
      if (ip<=6.and.lroot .and. (lnkmin0>=0.)) print *,'output_persistent_special: ',lnkmin0
!
!  write details
!
      output_persistent_special = .true.
!
      if (write_persist ('SPECIAL_LNKMIN0', id_record_SPECIAL_LNKMIN0, lnkmin0)) return
!
      output_persistent_special = .false.
!
    endfunction output_persistent_special
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  
    use Mpicomm, only: mpireduce_sum
    real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
    if (it == 1 .or. mod(it,ikdata)==0 .and. lfirst .and. llnk_spacing) then
      open (1, file=trim(directory_snap)//'/AR_and_AL.dat', form='formatted', position='append')
      write(1,*) t, lnk, f(l1:l2,m1,n1,iaxi_AR), f(l1:l2,m1,n1,iaxi_imAR), &
                               f(l1:l2,m1,n1,iaxi_AL), f(l1:l2,m1,n1,iaxi_imAL)
      close(1)
      open (1, file=trim(directory_snap)//'/ARdot_and_ALdot.dat', form='formatted', position='append')
      write(1,*) t, lnk, f(l1:l2,m1,n1,iaxi_ARdot), f(l1:l2,m1,n1,iaxi_imARdot), &
                               f(l1:l2,m1,n1,iaxi_ALdot), f(l1:l2,m1,n1,iaxi_imALdot)
      close(1)
    endif
    n=n1
    m=m1
    call calc_integrand(f)
    call mpireduce_sum(sum(edotb),edotb_sum,1)
    call mpireduce_sum(sum(e2),rhoe,1)
    call mpireduce_sum(sum(b2),rhob,1)
    if (mod(it,ikdata)==0 .and. lfirst .and. llnk_spacing) then
      open (1, file=trim(directory_snap)//'/backreact.dat', form='formatted', position='append')
      write(1,*) t, lnk, edotb, e2, b2
      close(1)
    endif
!            
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine calc_integrand(f)
!
      real, dimension (mx,my,mz,mfarray), intent(IN) :: f
      real, dimension (nx) :: AR, ARdot, imAR, imARdot, AReff2, ARARdoteff, imARARdoteff
      real, dimension (nx) :: ARdoteff2, imAReff2, imARdoteff2
      real, dimension (nx) :: AL, ALdot, imAL, imALdot, ALALdoteff, imALALdoteff
      real, dimension (nx) :: ALeff2,  ALdoteff2, imALeff2, imALdoteff2
      real, dimension (nx) :: psi, psidot, impsi , impsidot

      call get_Hubble_and_scale_factor(f_ode)
      AR   =f(l1:l2,m,n,iaxi_AR)
      ARdot=f(l1:l2,m,n,iaxi_ARdot)
      AL   =f(l1:l2,m,n,iaxi_AL)
      ALdot=f(l1:l2,m,n,iaxi_ALdot)
      imAR   =f(l1:l2,m,n,iaxi_imAR)
      imARdot=f(l1:l2,m,n,iaxi_imARdot)
      imAL   =f(l1:l2,m,n,iaxi_imAL)
      imALdot=f(l1:l2,m,n,iaxi_imALdot)
      AReff2=AR*AR
      ALeff2=AL*AL
      ARdoteff2=ARdot*ARdot
      ALdoteff2=ALdot*ALdot
      imAReff2=imAR*imAR
      imALeff2=imAL*imAL
      imARdoteff2=imARdot*imARdot
      imALdoteff2=imALdot*imALdot
      ARARdoteff=AR*ARdot
      imARARdoteff=imAR*imARdot
      ALALdoteff=AL*ALdot
      imALALdoteff=imAL*imALdot
      if (lquant_filter) then
        if (horizon_factor==0.) then
          if (headt.and.lfirst) print*,'horizon_factor=',horizon_factor
          where (AReff2<1./(2.*a*H))
            AReff2=0.
            ALeff2=0.
            imAReff2=0.
            imALeff2=0.
            ALdoteff2=0.
            ARdoteff2=0.
            imARdoteff2=0.
            imALdoteff2=0.
            ARARdoteff=0.
            imARARdoteff=0.
            ALALdoteff=0.
            imALALdoteff=0.
          endwhere
        elseif (horizon_factor>0.) then
          where (k>(a*H*horizon_factor))
            AReff2=0.
            ALeff2=0.
            imAReff2=0.
            imALeff2=0.
            ALdoteff2=0.
            ARdoteff2=0.
            imARdoteff2=0.
            imALdoteff2=0.
            ARARdoteff=0.
            imARARdoteff=0.
            ALALdoteff=0.
            imALALdoteff=0.
          endwhere
        endif
      endif
      if (llnk_spacing_adjustable .or. llnk_spacing) then
        b2=0.5*(4.*pi*k**5*dlnk)*a**(-4)*(AReff2+imAReff2+ALeff2+imALeff2)/twopi**3
        e2=0.5*(4.*pi*k**3*dlnk)*a**(-2)*(ARdoteff2+imARdoteff2+ALdoteff2+imALdoteff2)/twopi**3
        edotb=(4.*pi*k**4*dlnk)*a*(ARARdoteff+imARARdoteff-ALALdoteff-imALALdoteff)/twopi**3
      else
        b2=0.5*(4.*pi*k**4*dk)*(AReff2+imAReff2+ALeff2+imALeff2)/twopi**3
        e2=0.5*(4.*pi*k**2*dk)*a**2*(ARdoteff2+imARdoteff2+ALdoteff2+imALdoteff2)/twopi**3
        edotb=(4.*pi*k**3*dk)*a*(ARARdoteff+imARARdoteff-ALALdoteff-imALALdoteff)/twopi**3
      endif
    endsubroutine calc_integrand
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
      use Sub
!
!   SAMPLE IMPLEMENTATION
!
      integer :: iname, inamexy, inamexz
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
        idiag_AL=0;      idiag_ALdot=0;   idiag_imAL=0
        idiag_AR=0;      idiag_ARdot=0;   idiag_imAR=0
        idiag_imARdot=0; idiag_imALdot=0; idiag_k0=0;   
        idiag_dk=0;      idiag_rhoe=0
        idiag_rhob=0;    idiag_edotb=0
        if (lhubble) then
          idiag_a=0; idiag_phi=0; idiag_phidot=0; idiag_H=0
        endif
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'AL' ,idiag_AL)
        call parse_name(iname,cname(iname),cform(iname),'ALdot',idiag_ALdot)
        call parse_name(iname,cname(iname),cform(iname),'imAL' ,idiag_imAL)
        call parse_name(iname,cname(iname),cform(iname),'AR' ,idiag_AR)
        call parse_name(iname,cname(iname),cform(iname),'ARdot' ,idiag_ARdot)
        call parse_name(iname,cname(iname),cform(iname),'imAR' ,idiag_imAR)
        call parse_name(iname,cname(iname),cform(iname),'imARdot' ,idiag_imARdot)
        call parse_name(iname,cname(iname),cform(iname),'imALdot' ,idiag_imALdot)
        call parse_name(iname,cname(iname),cform(iname),'k0' ,idiag_k0)
        call parse_name(iname,cname(iname),cform(iname),'dk' ,idiag_dk)
        call parse_name(iname,cname(iname),cform(iname),'rhoe' ,idiag_rhoe)
        call parse_name(iname,cname(iname),cform(iname),'rhob' ,idiag_rhob)
        call parse_name(iname,cname(iname),cform(iname),'edotb' ,idiag_edotb)
        if (lhubble) then
          call parse_name(iname,cname(iname),cform(iname),'a' ,idiag_a)
          call parse_name(iname,cname(iname),cform(iname),'phi' ,idiag_phi)
          call parse_name(iname,cname(iname),cform(iname),'phidot' ,idiag_phidot)
          call parse_name(iname,cname(iname),cform(iname),'H' ,idiag_H)
        endif
      enddo
!
!  Check for those quantities for which we want z-averages.
!
!
    endsubroutine rprint_special
!***********************************************************************
! Subroutines below needed only for GPUs, if you do not care about GPUs don't worry about them
!***********************************************************************
    subroutine read_sums_from_GPU
      use GPU, only: get_gpu_reduced_vars
      real, dimension(12) :: tmp
      call get_gpu_reduced_vars(tmp)
      edotb_sum = tmp(1)
      rhoe      = tmp(2)
      rhob      = tmp(3)
    endsubroutine read_sums_from_GPU
!***********************************************************************
    subroutine pushpars2c(p_par)

    use Syscalls, only: copy_addr
    use General , only: string_to_enum

    integer, parameter :: n_pars=100
    integer(KIND=ikind8), dimension(n_pars) :: p_par

    call copy_addr(iaxi_al,p_par(6)) ! int
    call copy_addr(iaxi_aldot,p_par(7)) ! int
    call copy_addr(iaxi_imal,p_par(8)) ! int
    call copy_addr(iaxi_ar,p_par(9)) ! int
    call copy_addr(iaxi_ardot,p_par(10)) ! int
    call copy_addr(iaxi_imar,p_par(11)) ! int
    call copy_addr(iaxi_lna,p_par(22)) ! int
    call copy_addr(iaxi_phi,p_par(23)) ! int
    call copy_addr(iaxi_phidot,p_par(24)) ! int
    call copy_addr(ndivt,p_par(5)) ! int
    call copy_addr(alpf,p_par(25))
    call copy_addr(mpl2,p_par(29))
    call copy_addr(alpha,p_par(34))
    call copy_addr(m_alpha,p_par(35))
    call copy_addr(n_alpha,p_par(36))
    call copy_addr(lhubble,p_par(46)) ! bool
    call copy_addr(k,p_par(47)) ! (nx)
    call string_to_enum(enum_v_choice,v_choice)
    call copy_addr(enum_v_choice,p_par(48)) ! int
    call copy_addr(h_init,p_par(49))
    call copy_addr(llnk_spacing_adjustable,p_par(51)) ! bool
    call copy_addr(llnk_spacing,p_par(52)) ! bool
    call copy_addr(dlnk,p_par(53))
    call copy_addr(dk,p_par(54))
    call copy_addr(iaxi_imaldot,p_par(55)) ! int
    call copy_addr(iaxi_imardot,p_par(56)) ! int
    call copy_addr(m_phi,p_par(57))
    call copy_addr(phi_0,p_par(58))
    call copy_addr(lambda_gmssm,p_par(59))
    call copy_addr(alpha1_gmssm,p_par(60))
    call copy_addr(const_sigma,p_par(61))
    call copy_addr(lbackreact,p_par(64)) ! bool
    call copy_addr(lquant_filter,p_par(65)) ! bool
    call copy_addr(horizon_factor,p_par(65)) ! bool


    endsubroutine pushpars2c
!***********************************************************************
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
