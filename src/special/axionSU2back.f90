! $Id$
!
!  Solve the lucky droplet model for many realizations.
!  The different realizations correspond to "meshpoints".
!  To add the contributions for each step, we use the usual
!  time step in the Pencil Code, so t is just the step, and
!  the accumulated collision times (after 125 steps or so)
!  for all realizations at the same time are the values in
!  the f-array.
!
!  16-apr-20/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 4
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
!
! Declare index of variables
!
  integer :: iaxi_Q=0, iaxi_Qdot=0, iaxi_chi=0, iaxi_chidot=0, Ndivt=100
  integer :: iaxi_psi=0, iaxi_psidot=0, iaxi_TR=0, iaxi_TRdot=0
  integer :: iaxi_psiL=0, iaxi_psiLdot=0, iaxi_TL=0, iaxi_TLdot=0
  integer :: iaxi_impsi=0, iaxi_impsidot=0, iaxi_imTR=0, iaxi_imTRdot=0
  integer :: iaxi_impsiL=0, iaxi_impsiLdot=0, iaxi_imTL=0, iaxi_imTLdot=0
!
  ! input parameters
  real :: a, k0=1e-2, dk=1e-2, ascale_ini=1.
  real :: fdecay=.003, g=1.11e-2, lam=500., mu=1.5e-4
  real :: Q0=3e-4, Qdot0=0., chi_prefactor=.49, chidot0=0., H=1.04e-6
  real :: Mpl2=1., Hdot=0., lamf, Hscript
  real, dimension (nx) :: grand, grant, dgrant
  real, dimension (nx) :: xmask_axion
  real, dimension (2) :: axion_sum_range=(/0.,1./)
  integer, dimension (nx) :: kindex_array
  real :: grand_sum, grant_sum, dgrant_sum
  real :: TRdoteff2km_sum, TRdoteff2m_sum, TReff2km_sum, TReff2m_sum
  real :: TLdoteff2km_sum, TLdoteff2m_sum, TLeff2km_sum, TLeff2m_sum
  real :: TRpsim_sum, TRpsikm_sum, TRpsidotm_sum, TRdotpsim_sum
  real :: sbackreact_Q=1., sbackreact_chi=1., tback=1e6, dtback=1e6
  real :: lnkmin0, lnkmin0_dummy, lnkmax0, dlnk
  real :: nmin0=-1., nmax0=3., horizon_factor=0.
  real, dimension (nx) :: dt1_special, lnk
  logical :: lbackreact=.false., lwith_eps=.true., lupdate_background=.true.
  logical :: ldo_adjust_krange=.true., lswap_sign=.false.
  logical :: lconf_time=.false., lanalytic=.false., lvariable_k=.false.
  logical :: llnk_spacing_adjustable=.false., llnk_spacing=.false.
  logical :: lim_psi_TR=.false., lleft_psiL_TL=.false., lkeep_mQ_const=.false.
  character(len=50) :: init_axionSU2back='standard'
  namelist /special_init_pars/ &
    k0, dk, fdecay, g, lam, mu, Q0, Qdot0, chi_prefactor, chidot0, H, &
    lconf_time, Ndivt, lanalytic, lvariable_k, axion_sum_range, &
    llnk_spacing_adjustable, llnk_spacing, lim_psi_TR, lleft_psiL_TL, &
    nmin0, nmax0, ldo_adjust_krange, lswap_sign
!
  ! run parameters
  namelist /special_run_pars/ &
    k0, dk, fdecay, g, lam, mu, H, lwith_eps, lupdate_background, &
    lbackreact, sbackreact_Q, sbackreact_chi, tback, dtback, lconf_time, &
    Ndivt, lanalytic, lvariable_k, llnk_spacing_adjustable, llnk_spacing, &
    nmin0, nmax0, horizon_factor, axion_sum_range, lkeep_mQ_const, &
    ldo_adjust_krange, lswap_sign
!
  ! k array
  real, dimension (nx) :: k
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_Q   =0 ! DIAG_DOC: $Q$
  integer :: idiag_Qdot=0 ! DIAG_DOC: $\dot{Q}$
  integer :: idiag_Qddot=0! DIAG_DOC: $\ddot{Q}$
  integer :: idiag_chi =0 ! DIAG_DOC: $\chi$
  integer :: idiag_chidot =0 ! DIAG_DOC: $\dot{\chi}$
  integer :: idiag_chiddot =0 ! DIAG_DOC: $\ddot{\chi}$
  integer :: idiag_psi     =0 ! DIAG_DOC: $\psi$
  integer :: idiag_psiL    =0 ! DIAG_DOC: $\psi_L$
  integer :: idiag_psidot =0 ! DIAG_DOC: $\dot\psi$
  integer :: idiag_psiddot=0 ! DIAG_DOC: $\ddot\psi$
  integer :: idiag_TR     =0 ! DIAG_DOC: $T_R$
  integer :: idiag_TL     =0 ! DIAG_DOC: $T_L$
  integer :: idiag_TRdot =0 ! DIAG_DOC: $\dot T_R$
  integer :: idiag_TRddot=0 ! DIAG_DOC: $\ddot T_R$
  integer :: idiag_imTR=0 ! DIAG_DOC: $\Im T_R$
  integer :: idiag_psi_anal =0 ! DIAG_DOC: $\psi^{\rm anal}$
  integer :: idiag_TR_anal  =0 ! DIAG_DOC: $T_R^{\rm anal}$
  integer :: idiag_TReff2m  =0 ! DIAG_DOC: $|T_R|^2_{\rm eff}$
  integer :: idiag_TReff2km  =0 ! DIAG_DOC: $k|T_R|^2_{\rm eff}$
  integer :: idiag_TRdoteff2m  =0 ! DIAG_DOC: $|T_R\dot{T}_R|_{\rm eff}$
  integer :: idiag_TRdoteff2km  =0 ! DIAG_DOC: $k|T_R\dot{T}_R|_{\rm eff}$
  integer :: idiag_TRpsim  =0 ! DIAG_DOC: $\langle T_R^* \psi\rangle$
  integer :: idiag_TRpsikm  =0 ! DIAG_DOC: $\langle T_R^* \psi (k/a)\rangle$
  integer :: idiag_TRpsidotm  =0 ! DIAG_DOC: $\langle T_R^* \psi'\rangle$
  integer :: idiag_TRdotpsim  =0 ! DIAG_DOC: $\langle T_R^*' \psi\rangle$
  integer :: idiag_TLeff2m  =0 ! DIAG_DOC: $|T_R|^2_{\rm eff}$
  integer :: idiag_TLeff2km  =0 ! DIAG_DOC: $k|T_R|^2_{\rm eff}$
  integer :: idiag_TLdoteff2m  =0 ! DIAG_DOC: $|T_R\dot{T}_R|_{\rm eff}$
  integer :: idiag_TLdoteff2km  =0 ! DIAG_DOC: $k|T_R\dot{T}_R|_{\rm eff}$
  integer :: idiag_dgrant_up=0 ! DIAG_DOC: ${\cal T}^\chi$
  integer :: idiag_grand2=0 ! DIAG_DOC: ${\cal T}^Q$ (test)
  integer :: idiag_dgrant=0 ! DIAG_DOC: $\dot{\cal T}^\chi$
  integer :: idiag_fact=0   ! DIAG_DOC: $\Theta(t)$
  integer :: idiag_k0=0   ! DIAG_DOC: $k0$
  integer :: idiag_dk=0   ! DIAG_DOC: $dk$
!
! z averaged diagnostics given in zaver.in
!
  integer :: idiag_grandxy=0      ! ZAVG_DOC: $\left< {\cal T}^Q\right>_{z}$
  integer :: idiag_grantxy=0      ! ZAVG_DOC: $\left< {\cal T}^\chi\right>_{z}$
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  19-feb-2019/axel: coded
!
      use FArrayManager
      use SharedVariables, only: put_shared_variable
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Set iaxionSU2back to consecutive numbers
!
      call farray_register_ode('axi_Q'     ,iaxi_Q)
      call farray_register_ode('axi_Qdot'  ,iaxi_Qdot)
      call farray_register_ode('axi_chi'   ,iaxi_chi)
      call farray_register_ode('axi_chidot',iaxi_chidot)
!
      call farray_register_pde('axi_psi'   ,iaxi_psi)
      call farray_register_pde('axi_psidot',iaxi_psidot)
      call farray_register_pde('axi_TR'    ,iaxi_TR)
      call farray_register_pde('axi_TRdot' ,iaxi_TRdot)
!
!  Possibility of including the imaginary parts.
!
      if (lim_psi_TR) then
        call farray_register_pde('axi_impsi'   ,iaxi_impsi)
        call farray_register_pde('axi_impsidot',iaxi_impsidot)
        call farray_register_pde('axi_imTR'    ,iaxi_imTR)
        call farray_register_pde('axi_imTRdot' ,iaxi_imTRdot)
      endif
!
!  Possibility of including left-handed modes.
!
      if (lleft_psiL_TL) then
        call farray_register_pde('axi_psiL'   ,iaxi_psiL)
        call farray_register_pde('axi_psiLdot',iaxi_psiLdot)
        call farray_register_pde('axi_TL'     ,iaxi_TL)
        call farray_register_pde('axi_TLdot'  ,iaxi_TLdot)
        call farray_register_pde('axi_impsiL'   ,iaxi_impsiL)
        call farray_register_pde('axi_impsiLdot',iaxi_impsiLdot)
        call farray_register_pde('axi_imTL'     ,iaxi_imTL)
        call farray_register_pde('axi_imTLdot'  ,iaxi_imTLdot)
      endif
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  19-feb-2019/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: lnH, lna, a
      real :: kmax=2., lnkmax, lnk0=1.
      integer :: ik
!
      lamf=lam/fdecay
!
!  Compute lnkmin0 and lnkmax0. Even for a linear k-range, dlnk
!  is needed to determine the output of k-range and grand etc.
!
      if (lconf_time) then
        a=-1./(H*t)
      else
        a=exp(H*t)
      endif
      lna=alog(a)
      lnH=alog(H)
      lnkmin0=nmin0+lnH+lna
      lnkmax0=nmax0+lnH+lna
      dlnk=(lnkmax0-lnkmin0)/(ncpus*nx-1)
!
!  Initialize lnkmin0 and lnkmax0
!
      if (llnk_spacing_adjustable .and. .not.lstart) then
        do ik=1,nx
          lnk(ik)=lnkmin0+dlnk*(ik-1+iproc*nx)
          k(ik)=exp(lnk(ik))
        enddo
        kindex_array=nint((lnk-lnkmin0)/dlnk)
      elseif (llnk_spacing) then
        do ik=1,nx
          lnk(ik)=lnkmin0+dlnk*(ik-1+iproc*nx)
          k(ik)=exp(lnk(ik))
        enddo
        lnkmin0_dummy=lnkmin0
        kindex_array=nint((lnk-lnkmin0)/dlnk)
      else
        do ik=1,nx
          k(ik)=k0+dk*(ik-1+iproc*nx)
        enddo
        lnk=impossible
        lnkmin0_dummy=nmin0+lnH+lna
        kindex_array=nint((k-k0)/dk)
      endif
!
!  Compute mask for error diagnostics
!
      if (l1 == l2) then
        xmask_axion = 1.
      else
        where (      kindex_array >= nint(mx*axion_sum_range(1)-1) &
               .and. kindex_array <= nint(mx*axion_sum_range(2)-1))
          xmask_axion = 1.
        elsewhere
          xmask_axion = 0.
        endwhere
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!   2-dec-2022/axel: coded
!
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: psi, psidot, TR, TRdot
      real, dimension (nx) :: impsi, impsidot, imTR, imTRdot
      real :: chi0, Uprime0
      real :: lnt, lnH, lna, a
      real :: kmax=2., lnkmax, lnk0=1.
      integer :: ik
!
      intent(inout) :: f
!
!  Initialize any module variables which are parameter dependent
!  Compute lnkmin0. Reset tstart if conformal time.
!
      if (lconf_time) then
        tstart=-1./(ascale_ini*H)
        t=tstart
        a=-1./(H*t)
      else
        a=exp(H*t)
      endif
      lna=alog(a)
      lnH=alog(H)
      lnkmin0=nmin0+lnH+lna
      lnkmax0=nmax0+lnH+lna
!
!  Different k prescriptions
!
      if (llnk_spacing_adjustable) then
        dlnk=(lnkmax0-lnkmin0)/(ncpus*nx-1)
        do ik=1,nx
          lnk(ik)=lnkmin0+dlnk*(ik-1+iproc*nx)
          k(ik)=exp(lnk(ik))
        enddo
        if (ip<10) print*,'iproc,lnk=',iproc,lnk
        kindex_array=nint((lnk-lnkmin0)/dlnk)
      elseif (llnk_spacing) then
        dlnk=(lnkmax0-lnkmin0)/(ncpus*nx-1)
        do ik=1,nx
          lnk(ik)=lnkmin0+dlnk*(ik-1+iproc*nx)
          k(ik)=exp(lnk(ik))
        enddo
        if (ip<10) print*,'iproc,lnk=',iproc,lnk
        kindex_array=nint((lnk-lnkmin0)/dlnk)
      else
        do ik=1,nx
          k(ik)=k0+dk*(ik-1+iproc*nx)
        enddo
        kindex_array=nint((k-k0)/dk)
      endif
!
!  Initial condition; depends on k, which is here set to x.
!
      select case (init_axionSU2back)
        case ('nothing'); if (lroot) print*,'nothing'
        case ('standard')
          if (lconf_time) then
            if (ip<10) print*,'k=',k
            psi=(1./sqrt(2.*k))*cos(-k*t)
            psidot=(k/sqrt(2.*k))*sin(-k*t)
            TR=(1./sqrt(2.*k))*cos(-k*t)
            TRdot=(k/sqrt(2.*k))*sin(-k*t)
            if (lim_psi_TR) then
              impsi=(1./sqrt(2.*k))*sin(-k*t)
              impsidot=(-k/sqrt(2.*k))*cos(-k*t)
              imTR=(1./sqrt(2.*k))*sin(-k*t)
              imTRdot=(-k/sqrt(2.*k))*cos(-k*t)
            endif
          else
            if (ip<10) print*,'k=',k
            a=exp(H*t)
            psi=(ascale_ini/sqrt(2.*k))*cos(k/(ascale_ini*H))
            psidot=(k/sqrt(2.*k))*sin(k/(ascale_ini*H))
            TR=(ascale_ini/sqrt(2.*k))*cos(k/(ascale_ini*H))
            TRdot=(k/sqrt(2.*k))*sin(k/(ascale_ini*H))
          endif
!
!  ODE variables (exist only on root processor)
!
          if (lroot) then
            chi0=chi_prefactor*pi*fdecay
            Uprime0=-mu**4/fdecay*sin(chi0/fdecay)
            Q0=(-Uprime0/(3.*g*lamf*H))**onethird
            f_ode(iaxi_Q)=Q0
            f_ode(iaxi_Qdot)=Qdot0
            f_ode(iaxi_chi)=chi0
            f_ode(iaxi_chidot)=chidot0
          endif
!
!  PDE variables
!
          do n=n1,n2
          do m=m1,m2
            f(l1:l2,m,n,iaxi_psi)=psi
            f(l1:l2,m,n,iaxi_psidot)=psidot
            f(l1:l2,m,n,iaxi_TR)=TR
            f(l1:l2,m,n,iaxi_TRdot)=TRdot
            if (lim_psi_TR) then
              f(l1:l2,m,n,iaxi_impsi)=impsi
              f(l1:l2,m,n,iaxi_impsidot)=impsidot
              f(l1:l2,m,n,iaxi_imTR)=imTR
              f(l1:l2,m,n,iaxi_imTRdot)=imTRdot
            endif
            if (lleft_psiL_TL) then
              f(l1:l2,m,n,iaxi_psiL     )=psi
              f(l1:l2,m,n,iaxi_psiLdot  )=psidot
              f(l1:l2,m,n,iaxi_TL       )=TR
              f(l1:l2,m,n,iaxi_TLdot    )=TRdot
              f(l1:l2,m,n,iaxi_impsiL   )=impsi
              f(l1:l2,m,n,iaxi_impsiLdot)=impsidot
              f(l1:l2,m,n,iaxi_imTL     )=imTR
              f(l1:l2,m,n,iaxi_imTLdot  )=imTRdot
            endif
          enddo
          enddo 
          write(6,1000) 'iproc,TR=',iproc,TR
!
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_axionSU2back: No such value for init_axionSU2back: ', trim(init_axionSU2back)
          call stop_it("")
      endselect
!
1000  format(a,2x,i3,1p80e14.6)
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
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
!   2-dec-2022/axel: coded
!   1-sep-2023/axel: implemented lwith_eps with conformal time.
!
      use General, only: random_number_wrapper
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: psi , psidot , psiddot , TR, TRdot, TRddot
      real, dimension (nx) :: psiL, psiLdot, psiLddot, TL, TLdot, TLddot
      real, dimension (nx) :: psi_anal, psidot_anal, TR_anal, TRdot_anal
      real, dimension (nx) :: impsi , impsidot , impsiddot , imTR, imTRdot, imTRddot
      real, dimension (nx) :: impsiL, impsiLdot, impsiLddot, imTL, imTLdot, imTLddot
      real, dimension (nx) :: Uprime, mQ, xi, epsQE, epsQB
      real :: Q, Qdot, chi, chidot
      real :: fact=1., sign_swap=1.
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
      psi=f(l1:l2,m,n,iaxi_psi)
      psidot=f(l1:l2,m,n,iaxi_psidot)
      TR=f(l1:l2,m,n,iaxi_TR)
      TRdot=f(l1:l2,m,n,iaxi_TRdot)
!
      if (lim_psi_TR) then
        impsi=f(l1:l2,m,n,iaxi_impsi)
        impsidot=f(l1:l2,m,n,iaxi_impsidot)
        imTR=f(l1:l2,m,n,iaxi_imTR)
        imTRdot=f(l1:l2,m,n,iaxi_imTRdot)
      endif
!
      if (lleft_psiL_TL) then
        psiL     =f(l1:l2,m,n,iaxi_psiL)
        psiLdot  =f(l1:l2,m,n,iaxi_psiLdot)
        TL       =f(l1:l2,m,n,iaxi_TL)
        TLdot    =f(l1:l2,m,n,iaxi_TLdot)
        impsiL   =f(l1:l2,m,n,iaxi_impsiL)
        impsiLdot=f(l1:l2,m,n,iaxi_impsiLdot)
        imTL     =f(l1:l2,m,n,iaxi_imTL)
        imTLdot  =f(l1:l2,m,n,iaxi_imTLdot)
      endif
!
!  make ODE variables available (should exist on all processors)
!
      Q=f_ode(iaxi_Q)
      Qdot=f_ode(iaxi_Qdot)
      chi=f_ode(iaxi_chi)
      chidot=f_ode(iaxi_chidot)
!
!  Possibility of keeping mQ constant, i,e., we keep mQ=g*Q0/H
!  Need to have Q on all processors.
!
      if (lkeep_mQ_const) then
        mQ=g*Q0/H
      else
        mQ=g*Q/H
      endif
!
!  Set other parameters
!
      Uprime=-mu**4/fdecay*sin(chi/fdecay)
      if (lconf_time) then
        a=-1./(H*t)
        Hscript=a*H
        xi=lamf*chidot*(-0.5*t)
        epsQE=(Qdot/a+H*Q)**2/(Mpl2*H**2)
      else
        a=exp(H*t)
        xi=lamf*chidot/(2.*H)
        epsQE=(Qdot+H*Q)**2/(Mpl2*H**2)
      endif
      epsQB=g**2*Q**4/(Mpl2*H**2)
!
!  analytical solution
!
      if (lconf_time) then
        psi_anal=(1./sqrt(2.*k))*cos(-k*t)
        psidot_anal=(k/sqrt(2.*k))*sin(-k*t)
        TR_anal=(1./sqrt(2.*k))*cos(-k*t)
        TRdot_anal=(k/sqrt(2.*k))*sin(-k*t)
      else
        psi_anal=(1./sqrt(2.*k))*cos(k/(a*H))
        psidot_anal=(k/sqrt(2.*k))*sin(k/(a*H))
        TR_anal=(1./sqrt(2.*k))*cos(k/(a*H))
        TRdot_anal=(k/sqrt(2.*k))*sin(k/(a*H))
      endif
!
!  perturbation: when lwith_eps=T, we replace Q -> sqrt(epsQE)
!  and mQ*Q -> sqrt(epsQB).
!
      if (lconf_time) then
        if (lwith_eps) then
          psiddot=-(k**2-2./t**2)*psi &
            +(2.*sqrt(epsQE)/t)*TRdot+((2.*sqrt(epsQB)*(mQ+k*t))/t**2)*TR
          TRddot=-(k**2+(2.*(mQ*xi+k*t*(mQ+xi)))/t**2)*TR-(2.*sqrt(epsQE))/t*psidot &
            +(2.*sqrt(epsQE))/t**2*psi+(2.*sqrt(epsQB))/t**2*(mQ+k*t)*psi
          if (lim_psi_TR) then
            impsiddot=-(k**2-2./t**2)*impsi &
              +(2.*sqrt(epsQE)/t)*imTRdot+((2.*sqrt(epsQB)*(mQ+k*t))/t**2)*imTR
            imTRddot=-(k**2+(2.*(mQ*xi+k*t*(mQ+xi)))/t**2)*imTR-(2.*sqrt(epsQE))/t*impsidot &
              +(2.*sqrt(epsQE))/t**2*impsi+(2.*sqrt(epsQB))/t**2*(mQ+k*t)*impsi
          endif
!
!  Left-handed modes (conformal, with epsilon formulation)
!
          if (lleft_psiL_TL) then
            psiLddot=-(k**2-2./t**2)*psiL &
              +(2.*sqrt(epsQE)/t)*TLdot+((2.*sqrt(epsQB)*(mQ-k*t))/t**2)*TL
            TLddot=-(k**2+(2.*(mQ*xi-k*t*(mQ+xi)))/t**2)*TL-(2.*sqrt(epsQE))/t*psiLdot &
              +(2.*sqrt(epsQE))/t**2*psiL+(2.*sqrt(epsQB))/t**2*(mQ-k*t)*psiL
            impsiLddot=-(k**2-2./t**2)*impsiL &
              +(2.*sqrt(epsQE)/t)*imTLdot+((2.*sqrt(epsQB)*(mQ-k*t))/t**2)*imTL
            imTLddot=-(k**2+(2.*(mQ*xi-k*t*(mQ+xi)))/t**2)*imTL-(2.*sqrt(epsQE))/t*impsiLdot &
              +(2.*sqrt(epsQE))/t**2*impsiL+(2.*sqrt(epsQB))/t**2*(mQ-k*t)*impsiL
          endif
        else
          if (lswap_sign) then
            sign_swap=-1.
          else
            sign_swap=+1.
          endif
          psiddot=-(k**2-2.*(1.-Q**2*(mQ**2-1.))/t**2)*psi &
            +sign_swap*(2.*Q/t)*TRdot+((2.*mQ*Q*(mQ+k*t))/t**2)*TR
          TRddot=-(k**2+(2.*(mQ*xi+k*t*(mQ+xi)))/t**2)*TR-sign_swap*(2.*Q)/t*psidot &
            +sign_swap*(2.*Q)/t**2*psi+(2.*mQ*Q)/t**2*(mQ+k*t)*psi
          if (lim_psi_TR) then
            impsiddot=-(k**2-2.*(1-Q**2*(mQ**2-1.))/t**2)*impsi &
              +sign_swap*(2.*Q/t)*imTRdot+((2.*mQ*Q*(mQ+k*t))/t**2)*imTR
            imTRddot=-(k**2+(2.*(mQ*xi+k*t*(mQ+xi)))/t**2)*imTR-sign_swap*(2.*Q)/t*impsidot &
              +sign_swap*(2.*Q)/t**2*impsi+(2.*mQ*Q)/t**2*(mQ+k*t)*impsi
          endif
!
!  Left-handed modes
!
          if (lleft_psiL_TL) then
            psiLddot=-(k**2-2.*(1.-Q**2*(mQ**2-1.))/t**2)*psiL &
              +(2.*Q/t)*TLdot+((2.*mQ*Q*(mQ-k*t))/t**2)*TL
            TLddot=-(k**2+(2.*(mQ*xi-k*t*(mQ+xi)))/t**2)*TL-(2.*Q)/t*psiLdot &
              +(2.*Q)/t**2*psiL+(2.*mQ*Q)/t**2*(mQ-k*t)*psiL
            impsiLddot=-(k**2-2.*(1-Q**2*(mQ**2-1.))/t**2)*impsiL &
              +(2.*Q/t)*imTLdot+((2.*mQ*Q*(mQ-k*t))/t**2)*imTL
            imTLddot=-(k**2+(2.*(mQ*xi-k*t*(mQ+xi)))/t**2)*imTL-(2.*Q)/t*impsiLdot &
              +(2.*Q)/t**2*impsiL+(2.*mQ*Q)/t**2*(mQ-k*t)*impsiL
          endif
        endif
!
!  When lanalytic, overwrite corresponding points with zero.
!
        if (lanalytic) then
!          where (k>(a*H*4.5))
          where (k>(a*H*20.0))
            df(l1:l2,m,n,iaxi_psi)=0.
            df(l1:l2,m,n,iaxi_psidot)=0.
            df(l1:l2,m,n,iaxi_TR)=0.
            df(l1:l2,m,n,iaxi_TRdot)=0.
            f(l1:l2,m,n,iaxi_psi)=psi_anal
            f(l1:l2,m,n,iaxi_psidot)=psidot_anal
            f(l1:l2,m,n,iaxi_TR)=TR_anal
            f(l1:l2,m,n,iaxi_TRdot)=TRdot_anal
          elsewhere
            df(l1:l2,m,n,iaxi_psi)=df(l1:l2,m,n,iaxi_psi)+psidot
            df(l1:l2,m,n,iaxi_psidot)=df(l1:l2,m,n,iaxi_psidot)+psiddot
            df(l1:l2,m,n,iaxi_TR)=df(l1:l2,m,n,iaxi_TR)+TRdot
            df(l1:l2,m,n,iaxi_TRdot)=df(l1:l2,m,n,iaxi_TRdot)+TRddot
          endwhere
        else
          df(l1:l2,m,n,iaxi_psi   )=df(l1:l2,m,n,iaxi_psi   )+psidot
          df(l1:l2,m,n,iaxi_psidot)=df(l1:l2,m,n,iaxi_psidot)+psiddot
          df(l1:l2,m,n,iaxi_TR    )=df(l1:l2,m,n,iaxi_TR    )+TRdot
          df(l1:l2,m,n,iaxi_TRdot )=df(l1:l2,m,n,iaxi_TRdot )+TRddot
          if (lim_psi_TR) then
            df(l1:l2,m,n,iaxi_impsi   )=df(l1:l2,m,n,iaxi_impsi   )+impsidot
            df(l1:l2,m,n,iaxi_impsidot)=df(l1:l2,m,n,iaxi_impsidot)+impsiddot
            df(l1:l2,m,n,iaxi_imTR    )=df(l1:l2,m,n,iaxi_imTR    )+imTRdot
            df(l1:l2,m,n,iaxi_imTRdot )=df(l1:l2,m,n,iaxi_imTRdot )+imTRddot
          endif
!
!  left-handed modes (to be checked, not finalized)
!
          if (lleft_psiL_TL) then
            df(l1:l2,m,n,iaxi_psiL     )=df(l1:l2,m,n,iaxi_psiL     )+psiLdot
            df(l1:l2,m,n,iaxi_psiLdot  )=df(l1:l2,m,n,iaxi_psiLdot  )+psiLddot
            df(l1:l2,m,n,iaxi_TL       )=df(l1:l2,m,n,iaxi_TL       )+TLdot
            df(l1:l2,m,n,iaxi_TLdot    )=df(l1:l2,m,n,iaxi_TLdot    )+TLddot
            df(l1:l2,m,n,iaxi_impsiL   )=df(l1:l2,m,n,iaxi_impsiL   )+impsidot
            df(l1:l2,m,n,iaxi_impsiLdot)=df(l1:l2,m,n,iaxi_impsiLdot)+impsiddot
            df(l1:l2,m,n,iaxi_imTL     )=df(l1:l2,m,n,iaxi_imTL     )+imTLdot
            df(l1:l2,m,n,iaxi_imTLdot  )=df(l1:l2,m,n,iaxi_imTLdot  )+imTLddot
          endif
        endif
      else
!
!  same with cosmic time; should also write in terms of psiddot and TRddot
!  But this is not currently done because of mixed terms.
!
        if (lanalytic) then
!          where (k>(a*H*4.5))
          where (k>(a*H*20.0))
            df(l1:l2,m,n,iaxi_psi)=0.
            df(l1:l2,m,n,iaxi_psidot)=0.
            df(l1:l2,m,n,iaxi_TR)=0.
            df(l1:l2,m,n,iaxi_TRdot)=0.
            f(l1:l2,m,n,iaxi_psi)=psi_anal
            f(l1:l2,m,n,iaxi_psidot)=psidot_anal
            f(l1:l2,m,n,iaxi_TR)=TR_anal
            f(l1:l2,m,n,iaxi_TRdot)=TRdot_anal
          elsewhere
            df(l1:l2,m,n,iaxi_psi   )=df(l1:l2,m,n,iaxi_psi   )+psidot
            df(l1:l2,m,n,iaxi_psidot)=df(l1:l2,m,n,iaxi_psidot) &
              -H*psidot-(k**2/a**2-2.*H**2+2.*Q**2*H**2*(mQ**2-1.))*psi &
              -2.*H*Q*TRdot+2.*mQ*Q*H**2*(mQ-k/(a*H))*TR
            df(l1:l2,m,n,iaxi_TR)=df(l1:l2,m,n,iaxi_TR)+TRdot
            df(l1:l2,m,n,iaxi_TRdot)=df(l1:l2,m,n,iaxi_TRdot) &
              -H*TRdot-(k**2/a**2+2.*H**2*(mQ*xi-k/(a*H)*(mQ+xi)))*TR+2.*H*Q*psidot &
              +2.*Q*H**2*psi+2.*mQ*Q*H**2*(mQ-k/(a*H))*psi
          endwhere
        else
          df(l1:l2,m,n,iaxi_psi)=df(l1:l2,m,n,iaxi_psi)+psidot
          df(l1:l2,m,n,iaxi_TR )=df(l1:l2,m,n,iaxi_TR )+TRdot
          if (lwith_eps) then
            psiddot=-H*psidot-(k**2/a**2-2.*H**2)*psi-2.*H*sqrt(epsQE)*TRdot &
              +2.*H**2*sqrt(epsQB)*(mQ-k/(a*H))*TR
          else
            psiddot=-H*psidot-(k**2/a**2-2.*H**2+2.*Q**2*H**2*(mQ**2-1.))*psi &
              -2.*H*Q*TRdot+2.*mQ*Q*H**2*(mQ-k/(a*H))*TR
          endif
          if (lwith_eps) then
            TRddot=-H*TRdot-(k**2/a**2+2.*H**2*(mQ*xi-k/(a*H)*(mQ+xi)))*TR &
              +2.*H*sqrt(epsQE)*psidot &
              +2.*H**2*(sqrt(epsQB)*(mQ-k/(a*H))+sqrt(epsQE))*psi
          else
            TRddot=-H*TRdot-(k**2/a**2+2.*H**2*(mQ*xi-k/(a*H)*(mQ+xi)))*TR+2.*H*Q*psidot &
              +2.*Q*H**2*psi+2.*mQ*Q*H**2*(mQ-k/(a*H))*psi
          endif
          df(l1:l2,m,n,iaxi_psidot)=df(l1:l2,m,n,iaxi_psidot)+psiddot
          df(l1:l2,m,n,iaxi_TRdot)=df(l1:l2,m,n,iaxi_TRdot)+TRddot
        endif
      endif
!
      if (lfirst.and.ldt.and.lconf_time) then
        dt1_special = Ndivt*abs(Hscript)
        dt1_max=max(dt1_max,dt1_special)
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        call sum_mn_name(psi_anal,idiag_psi_anal)
        call sum_mn_name(TR_anal,idiag_TR_anal)
        call sum_mn_name(psi ,idiag_psi)
        call sum_mn_name(psiL,idiag_psiL)
        call sum_mn_name(psidot,idiag_psidot)
        call sum_mn_name(psiddot,idiag_psiddot)
        call sum_mn_name(TR,idiag_TR)
        call sum_mn_name(TL,idiag_TL)
        call sum_mn_name(TRdot,idiag_TRdot)
        call sum_mn_name(TRddot,idiag_TRddot)
        call sum_mn_name(imTR,idiag_imTR)
        call save_name(TReff2m_sum,idiag_TReff2m)
        call save_name(TReff2km_sum,idiag_TReff2km)
        call save_name(TRdoteff2m_sum,idiag_TRdoteff2m)
        call save_name(TRdoteff2km_sum,idiag_TRdoteff2km)
        call save_name(TRpsim_sum,idiag_TRpsim)
        call save_name(TRpsikm_sum,idiag_TRpsikm)
        call save_name(TRpsidotm_sum,idiag_TRpsidotm)
        call save_name(TRdotpsim_sum,idiag_TRdotpsim)
        call save_name(TLeff2m_sum,idiag_TLeff2m)
        call save_name(TLeff2km_sum,idiag_TLeff2km)
        call save_name(TLdoteff2m_sum,idiag_TLdoteff2m)
        call save_name(TLdoteff2km_sum,idiag_TLdoteff2km)
        call save_name(grand_sum,idiag_grand2)
        call save_name(dgrant_sum,idiag_dgrant)
        call sum_mn_name(dgrant*xmask_axion,idiag_dgrant_up,lplain=.true.)
        call save_name(fact,idiag_fact)
        call save_name(k0,idiag_k0)
        call save_name(dk,idiag_dk)
      endif
!
      if (l2davgfirst) then
        if (idiag_grandxy/=0)   call zsum_mn_name_xy(grand,idiag_grandxy)
        if (idiag_grantxy/=0)   call zsum_mn_name_xy(grant,idiag_grantxy)
      endif
!
1000  format(a,2x,i3,1p80e14.6)
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine dspecial_dt_ode
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
!   2-dec-2022/axel: coded
!   1-sep-2023/axel: implemented lwith_eps with conformal time.
!
      use General, only: random_number_wrapper
      use Diagnostics
      use Mpicomm
      use Sub
!
      real :: Q, Qdot, Qddot, chi, chidot, chiddot
      real :: Uprime, mQ, xi
      real :: fact=1., sign_swap=1.
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
!  Set the all variable
!
      Q=f_ode(iaxi_Q)
      Qdot=f_ode(iaxi_Qdot)
      chi=f_ode(iaxi_chi)
      chidot=f_ode(iaxi_chidot)
!
!  Possibility of keeping mQ constant, i,e., we keep mQ=g*Q0/H
!
      if (lkeep_mQ_const) then
        mQ=g*Q0/H
      else
        mQ=g*Q/H
      endif
!
!  Set other parameters
!
      Uprime=-mu**4/fdecay*sin(chi/fdecay)
      if (lconf_time) then
        a=-1./(H*t)
        Hscript=a*H
        xi=lamf*chidot*(-0.5*t)
      else
        a=exp(H*t)
        xi=lamf*chidot/(2.*H)
      endif
!
!  background
!
      if (lconf_time) then
        Qddot=g*lamf*a*chidot*Q**2-2.*Hscript*Qdot-(Hdot*a**2+2*Hscript**2)*Q-2.*g**2*a**2*Q**3
        chiddot=-3.*g*lamf*a*Q**2*(Qdot+Hscript*Q)-2.*Hscript*chidot-a**2*Uprime
      else
        Qddot=g*lamf*chidot*Q**2-3.*H*Qdot-(Hdot+2*H**2)*Q-2.*g**2*Q**3
        chiddot=-3.*g*lamf*Q**2*(Qdot+H*Q)-3.*H*chidot-Uprime
      endif
!
!  Optionally, include backreaction
!
      if (lbackreact) then
!
!  compute factor to switch on bachreaction gradually:
!
        if (tback/=0.) then
          fact=cubic_step(real(t),tback,dtback)
        else
          fact=1.
        endif
!
!  apply factor to switch on bachreaction gradually:
!
        if (lconf_time) then
          Qddot  =Qddot  -sbackreact_Q  *fact*a**2 *grand_sum
          chiddot=chiddot-sbackreact_chi*fact*a**2*dgrant_sum
        else
          Qddot  =Qddot  -sbackreact_Q  *fact *grand_sum
          chiddot=chiddot-sbackreact_chi*fact*dgrant_sum
        endif
      endif
!
!  Choice whether or not we want to update the background
!
      if (lupdate_background) then
        df_ode(iaxi_Q)     =df_ode(iaxi_Q     )+Qdot
        df_ode(iaxi_chi)   =df_ode(iaxi_chi   )+chidot
        df_ode(iaxi_Qdot)  =df_ode(iaxi_Qdot  )+Qddot
        df_ode(iaxi_chidot)=df_ode(iaxi_chidot)+chiddot
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        call save_name(Q   ,idiag_Q)
        call save_name(Qdot,idiag_Qdot)
        call save_name(Qddot,idiag_Qddot)
        call save_name(chi   ,idiag_chi)
        call save_name(chidot,idiag_chidot)
        call save_name(chiddot,idiag_chiddot)
      endif
!
    endsubroutine dspecial_dt_ode
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
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
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_run_pars, IOSTAT=iostat)
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
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  13-may-18/axel: added remove_mean_value for hij and gij
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
!  07-aug-17/axel: coded

      use Mpicomm, only: mpireduce_sum, mpibcast
!
      real, dimension (mx,my,mz,mfarray) :: tmp
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: TR, TRdot, imTR, imTRdot, TReff2, TRdoteff2
      real, dimension (nx) :: TL, TLdot, imTL, imTLdot, TLeff2, TLdoteff2
      real, dimension (nx) :: psi, psidot, impsi , impsidot
      real, dimension (nx) :: TRdoteff2km, TRdoteff2m, TReff2km, TReff2m
      real, dimension (nx) :: TRpsi , TRpsik , TRpsidot , TRdotpsi
      real, dimension (nx) :: TRpsim, TRpsikm, TRpsidotm, TRdotpsim
      real, dimension (nx) :: TLdoteff2km, TLdoteff2m, TLeff2km, TLeff2m
      real, dimension (nx) :: tmp_psi , tmp_psidot , tmp_TR, tmp_TRdot
      real, dimension (nx) :: tmp_psiL, tmp_psiLdot, tmp_TL, tmp_TLdot
      real, dimension (nx) :: tmp_impsi, tmp_impsidot, tmp_imTR, tmp_imTRdot
      real, dimension (nx) :: tmp_impsiL, tmp_impsiLdot, tmp_imTL, tmp_imTLdot
      real :: Q, Qdot, chi, chidot
      real :: mQ, xi
      real :: lnt, lnH, lna, a, lnkmin, lnkmax
      integer :: ik, nswitch
!
!  make ODE variables available (should exist on all processors)
!
      Q=f_ode(iaxi_Q)
      Qdot=f_ode(iaxi_Qdot)
      chidot=f_ode(iaxi_chidot)
!
!  Possibility of keeping mQ constant
!
      if (lkeep_mQ_const) then
        mQ=g*Q0/H
      else
        mQ=g*Q/H
      endif
!
!  For conformal time, there is a 1/a factor in Qdot/a+H
!
      if (lconf_time) then
        a=-1./(H*t)
        xi=lamf*chidot*(-0.5*t)
      else
        a=exp(H*t)
        xi=lamf*chidot/(2.*H)
      endif
!
!  decide about revising the k array
!  a=exp(N), N=H*t (=lna).
!
      if (llnk_spacing_adjustable .and. lfirst) then
        lna=alog(a)
        lnH=alog(H)
        lnkmin=nmin0+lnH+lna
        lnkmax=nmax0+lnH+lna
        if (lnkmin >= (lnkmin0+dlnk)) then
          nswitch=int((lnkmin-lnkmin0)/dlnk)
          if (ip<10) print*,'nswitch: ',a, lnkmin0, nswitch
          if (nswitch==0) call fatal_error('special_after_boundary','nswitch must not be zero')
          if (nswitch>1) call fatal_error('special_after_boundary','nswitch must not exceed 1')
!
!  calculate new k array (because nswitch=1)
!
          if (ldo_adjust_krange) then
            do ik=1,nx
              lnk(ik)=lnkmin0+dlnk*(ik-1+iproc*nx+nswitch)
              k(ik)=exp(lnk(ik))
            enddo
            kindex_array=nint((lnk-lnkmin0)/dlnk)
!
!  move f array. Compute new initial point for the last entry.
!
            tmp(l1:l2,:,:,iaxi_psi:iaxi_TRdot)=f(l1+nswitch:l2+nswitch,:,:,iaxi_psi:iaxi_TRdot)
            f(l1:l2,:,:,iaxi_psi:iaxi_TRdot)=tmp(l1:l2,:,:,iaxi_psi:iaxi_TRdot)
!
            if (lim_psi_TR) then
              tmp(l1:l2,:,:,iaxi_impsi:iaxi_imTRdot)=f(l1+nswitch:l2+nswitch,:,:,iaxi_impsi:iaxi_imTRdot)
              f(l1:l2,:,:,iaxi_impsi:iaxi_imTRdot)=tmp(l1:l2,:,:,iaxi_impsi:iaxi_imTRdot)
            endif
!
!  move f array. Compute new initial point for the last entry.
!  Compute still for the full array (but only on the last processor),
!  but only use the last point below. (TR is therefore no longer ok after this.)
!
            if (ipx==nprocx-1) then
              tmp_psi   =(1./sqrt(2.*k))*cos(-k*t)
              tmp_psidot= (k/sqrt(2.*k))*sin(-k*t)
              tmp_TR    =(1./sqrt(2.*k))*cos(-k*t)
              tmp_TRdot = (k/sqrt(2.*k))*sin(-k*t)
              n=n1
              m=m1
              f(l2,m,n,iaxi_psi)   =tmp_psi(nx)
              f(l2,m,n,iaxi_psidot)=tmp_psidot(nx)
              f(l2,m,n,iaxi_TR)    =tmp_TR(nx)
              f(l2,m,n,iaxi_TRdot) =tmp_TRdot(nx)
              if (lim_psi_TR) then
                tmp_impsi=(1./sqrt(2.*k))*sin(-k*t)
                tmp_impsidot=(-k/sqrt(2.*k))*cos(-k*t)
                tmp_imTR=(1./sqrt(2.*k))*sin(-k*t)
                tmp_imTRdot=(-k/sqrt(2.*k))*cos(-k*t)
                f(l2,m,n,iaxi_impsi)   =tmp_impsi(nx)
                f(l2,m,n,iaxi_impsidot)=tmp_impsidot(nx)
                f(l2,m,n,iaxi_imTR)    =tmp_imTR(nx)
                f(l2,m,n,iaxi_imTRdot) =tmp_imTRdot(nx)
              endif
!
!  Same for left-handed modes
!
              if (lleft_psiL_TL) then
                tmp_psiL   =(1./sqrt(2.*k))*cos(-k*t)
                tmp_psiLdot= (k/sqrt(2.*k))*sin(-k*t)
                tmp_TL    =(1./sqrt(2.*k))*cos(-k*t)
                tmp_TLdot = (k/sqrt(2.*k))*sin(-k*t)
                n=n1
                m=m1
                f(l2,m,n,iaxi_psiL)   =tmp_psiL(nx)
                f(l2,m,n,iaxi_psiLdot)=tmp_psiLdot(nx)
                f(l2,m,n,iaxi_TL)    =tmp_TL(nx)
                f(l2,m,n,iaxi_TLdot) =tmp_TLdot(nx)
                tmp_impsiL=(1./sqrt(2.*k))*sin(-k*t)
                tmp_impsiLdot=(-k/sqrt(2.*k))*cos(-k*t)
                tmp_imTL=(1./sqrt(2.*k))*sin(-k*t)
                tmp_imTLdot=(-k/sqrt(2.*k))*cos(-k*t)
                f(l2,m,n,iaxi_impsiL)   =tmp_impsiL(nx)
                f(l2,m,n,iaxi_impsiLdot)=tmp_impsiLdot(nx)
                f(l2,m,n,iaxi_imTL)    =tmp_imTL(nx)
                f(l2,m,n,iaxi_imTLdot) =tmp_imTLdot(nx)
              endif
            endif
          endif
!
!  output
!
          open (1, file=trim(directory_snap)//'/krange.dat', form='formatted', position='append')
          write(1,*) t, lnk, f(l1:l2,m1,n1,iaxi_psi), f(l1:l2,m1,n1,iaxi_impsi), &
                             f(l1:l2,m1,n1,iaxi_TR), f(l1:l2,m1,n1,iaxi_imTR)
          close(1)
!
!  output for left-handed modes
!
          if (lleft_psiL_TL) then
          open (1, file=trim(directory_snap)//'/krange_left.dat', form='formatted', position='append')
          write(1,*) t, lnk, f(l1:l2,m1,n1,iaxi_psiL), f(l1:l2,m1,n1,iaxi_impsiL), &
                             f(l1:l2,m1,n1,iaxi_TL),   f(l1:l2,m1,n1,iaxi_imTL)
          close(1)
          endif
!
!  reset lnkmin0
!
          lnkmin0=lnkmin0+dlnk
        else
          nswitch=0
        endif
!
!  for llnk_spacing=T, we still want output at the same times as in
!  the adjustable case, so her nswitch means just "output", but no switch
!
      elseif (lfirst) then
        lna=alog(a)
        lnH=alog(H)
        lnkmin=nmin0+lnH+lna
        if (lnkmin >= (lnkmin0_dummy+dlnk)) then
          nswitch=int((lnkmin-lnkmin0_dummy)/dlnk)
          if (nswitch==0) call fatal_error('special_after_boundary','nswitch must not be zero')
          if (nswitch>1) call fatal_error('special_after_boundary','nswitch must not exceed 1')
          open (1, file=trim(directory_snap)//'/krange.dat', form='formatted', position='append')
          write(1,*) t, lnk, f(l1:l2,m1,n1,iaxi_psi)
          close(1)
!
!  reset lnkmin0_dummy (but now it is a dummy)
!
          lnkmin0_dummy=lnkmin0_dummy+dlnk
        else
          nswitch=0
        endif
      endif
!
!  Now set TR, TRdot, and imaginary parts, after they have been updated.
!
      psi   =f(l1:l2,m,n,iaxi_psi)
      psidot=f(l1:l2,m,n,iaxi_psidot)
      TR   =f(l1:l2,m1,n1,iaxi_TR)
      TRdot=f(l1:l2,m1,n1,iaxi_TRdot)
      if (lim_psi_TR) then
        impsi   =f(l1:l2,m,n,iaxi_impsi)
        impsidot=f(l1:l2,m,n,iaxi_impsidot)
        imTR   =f(l1:l2,m1,n1,iaxi_imTR)
        imTRdot=f(l1:l2,m1,n1,iaxi_imTRdot)
      endif
!
!  integrand (for diagnostics)
!  Here, "eff" refers to either just the real part squared, or the modulus
!
      TReff2=TR**2
      TRdoteff2=TR*TRdot
      TRpsi=TR*psi
      TRpsidot=TR*psidot
      TRdotpsi=TRdot*psi
      if (lim_psi_TR) then
        !TReff2=TR**2+imTR**2
        !TRdoteff2=TR*TRdot+imTR*imTRdot
        TReff2=TReff2+imTR**2
        TRdoteff2=TRdoteff2+imTR*imTRdot
        TRpsi=TRpsi+imTR*impsi
        TRpsidot=TRpsidot+imTR*impsidot
        TRdotpsi=TRdotpsi+imTRdot*impsi
      endif
!
!  Apply horizon factor
!
      if (horizon_factor==0.) then
        if (headt.and.lfirst) print*,'horizon_factor=',horizon_factor
        where (TReff2<1./(2.*a*H))
          TReff2=0.
          TRdoteff2=0.
!         TReff=(1./sqrt(2.*k))*cos(-k*t)
!         TRdoteff=(k/sqrt(2.*k))*sin(-k*t)
          TRpsi=0.
          TRpsidot=0.
          TRdotpsi=0.
        endwhere
      elseif (horizon_factor>0.) then
        where (k>(a*H*horizon_factor))
          TReff2=0.
          TRdoteff2=0.
          TRpsi=0.
          TRpsidot=0.
          TRdotpsi=0.
        endwhere
      endif
!
!  Same for left-handed modes
!
      if (lleft_psiL_TL) then
        TL   =f(l1:l2,m1,n1,iaxi_TL)
        TLdot=f(l1:l2,m1,n1,iaxi_TLdot)
        imTL   =f(l1:l2,m1,n1,iaxi_imTL)
        imTLdot=f(l1:l2,m1,n1,iaxi_imTLdot)
        TLeff2=TL**2
        TLdoteff2=TL*TLdot
        TLeff2=TLeff2+imTL**2
        TLdoteff2=TLdoteff2+imTL*imTLdot
        if (horizon_factor==0.) then
          if (headt.and.lfirst) print*,'horizon_factor=',horizon_factor
          where (TLeff2<1./(2.*a*H))
            TLeff2=0.
            TLdoteff2=0.
          endwhere
        elseif (horizon_factor>0.) then
          where (k>(a*H*horizon_factor))
            TLeff2=0.
            TLdoteff2=0.
          endwhere
        endif
      endif
!
!  Calculation of grand and grant is different for logarithmic and
!  uniform spacings.
!
      if (llnk_spacing_adjustable .or. llnk_spacing) then
        TRpsim=(4.*pi*k**3*dlnk)*TRpsi
        TRpsikm=(4.*pi*k**3*dlnk)*TRpsi*(k/a)
        TRpsidotm=(4.*pi*k**3*dlnk)*TRpsidot
        TRdotpsim=(4.*pi*k**3*dlnk)*TRdotpsi
        TRdoteff2km=(4.*pi*k**3*dlnk)*TRdoteff2*(k/a)
        TRdoteff2m=(4.*pi*k**3*dlnk)*TRdoteff2
        TReff2km=(4.*pi*k**3*dlnk)*TReff2*(k/a)
        TReff2m=(4.*pi*k**3*dlnk)*TReff2
        grand=(4.*pi*k**3*dlnk)*(xi*H-k/a)*TReff2*(+   g/(3.*a**2))/twopi**3
        grant=(4.*pi*k**3*dlnk)*(mQ*H-k/a)*TReff2*(-lamf/(2.*a**2))/twopi**3
        if (lleft_psiL_TL) then
          TLdoteff2km=(4.*pi*k**3*dlnk)*TLdoteff2*(k/a)
          TLdoteff2m=(4.*pi*k**3*dlnk)*TLdoteff2
          TLeff2km=(4.*pi*k**3*dlnk)*TLeff2*(k/a)
          TLeff2m=(4.*pi*k**3*dlnk)*TLeff2
          grand=grand+(4.*pi*k**3*dlnk)*(xi*H-k/a)*TLeff2*(+   g/(3.*a**2))/twopi**3
          grant=grant+(4.*pi*k**3*dlnk)*(mQ*H-k/a)*TLeff2*(-lamf/(2.*a**2))/twopi**3
        endif
        if (lconf_time) then
          dgrant=(4.*pi*k**3*dlnk)*(-lamf/(2.*a**3))*( &
          (a*mQ*H**2+g*Qdot)*TReff2+(mQ*H-k/a)*2*TRdoteff2 &
          )/twopi**3
          if (lleft_psiL_TL) then
            dgrant=dgrant+(4.*pi*k**3*dlnk)*(-lamf/(2.*a**3))*( &
            (a*mQ*H**2+g*Qdot)*TLeff2+(mQ*H-k/a)*2*TLdoteff2 &
            )/twopi**3
          endif
        else
          dgrant=(4.*pi*k**3*dlnk)*(-lamf/(2.*a**3))*( &
          (a*mQ*H**2+a*g*Qdot)*TReff2+(a*mQ*H-k)*2*TRdoteff2 &
          )/twopi**3
          if (lleft_psiL_TL) then
            dgrant=dgrant+(4.*pi*k**3*dlnk)*(-lamf/(2.*a**3))*( &
            (a*mQ*H**2+a*g*Qdot)*TLeff2+(a*mQ*H-k)*2*TLdoteff2 &
            )/twopi**3
          endif
        endif
      else
        TRpsidotm=(4.*pi*k**2*dk)*TRpsidot
        TRdotpsim=(4.*pi*k**2*dk)*TRdotpsi
        TRdoteff2km=(4.*pi*k**2*dk)*TRdoteff2*(k/a)
        TRdoteff2m=(4.*pi*k**2*dk)*TRdoteff2
        TReff2km=(4.*pi*k**2*dk)*TReff2*(k/a)
        TReff2m=(4.*pi*k**2*dk)*TReff2
        grand=(4.*pi*k**2*dk)*(xi*H-k/a)*TReff2*(+   g/(3.*a**2))/twopi**3
        grant=(4.*pi*k**2*dk)*(mQ*H-k/a)*TReff2*(-lamf/(2.*a**2))/twopi**3
        if (lconf_time) then
          dgrant=(4.*pi*k**2*dk)*(-lamf/(2.*a**3))*( &
          (a*mQ*H**2+g*Qdot)*TReff2+(mQ*H-k/a)*2*TRdoteff2 &
          )/twopi**3
        else 
          dgrant=(4.*pi*k**2*dk)*(-lamf/(2.*a**3))*( &
          (a*mQ*H**2+a*g*Qdot)*TReff2+(a*mQ*H-k)*2*TRdoteff2 &
          )/twopi**3
        endif
      endif
!
!  output of integrand
!
      if ((llnk_spacing_adjustable.or.llnk_spacing) .and. lfirst) then
        if (nswitch>0) then
          open (1, file=trim(directory_snap)//'/backreact.dat', form='formatted', position='append')
          write(1,*) t, lnk, grand, dgrant
          close(1)
        endif
      elseif (lfirst) then
        if (nswitch>0) then
          open (1, file=trim(directory_snap)//'/backreact.dat', form='formatted', position='append')
          write(1,*) t, k, grand, dgrant
          close(1)
        endif
      endif
!
!  Compute the sum over all processors.
!  But result is only needed on root processor.
!
      call mpireduce_sum(sum(grand),grand_sum,1)
      call mpireduce_sum(sum(grant),grant_sum,1)
      call mpireduce_sum(sum(dgrant),dgrant_sum,1)
!
!  These 8 lines are only needed for diagnostics and could be escaped.
!
      call mpireduce_sum(sum(TRpsim),TRpsim_sum,1)
      call mpireduce_sum(sum(TRpsikm),TRpsikm_sum,1)
      call mpireduce_sum(sum(TRpsidotm),TRpsidotm_sum,1)
      call mpireduce_sum(sum(TRdotpsim),TRdotpsim_sum,1)
!
      call mpireduce_sum(sum(TRdoteff2km),TRdoteff2km_sum,1)
      call mpireduce_sum(sum(TRdoteff2m),TRdoteff2m_sum,1)
      call mpireduce_sum(sum(TReff2km),TReff2km_sum,1)
      call mpireduce_sum(sum(TReff2m),TReff2m_sum,1)
!
!  Same for left-handed modes
!
      call mpireduce_sum(sum(TLdoteff2km),TLdoteff2km_sum,1)
      call mpireduce_sum(sum(TLdoteff2m),TLdoteff2m_sum,1)
      call mpireduce_sum(sum(TLeff2km),TLeff2km_sum,1)
      call mpireduce_sum(sum(TLeff2m),TLeff2m_sum,1)
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!  19-feb-2019/axel: coded
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
        idiag_Q=0; idiag_Qdot=0; idiag_Qddot=0; idiag_chi=0; idiag_chidot=0; idiag_chiddot=0
        idiag_psi=0; idiag_psiL=0; idiag_psidot=0; idiag_psiddot=0
        idiag_TR=0; idiag_TL=0; idiag_TRdot=0; idiag_TRddot=0; idiag_psi_anal=0; idiag_TR_anal=0
        idiag_TReff2m=0; idiag_TReff2km=0; idiag_TRdoteff2m=0; idiag_TRdoteff2km=0
        idiag_TRpsim=0; idiag_TRpsikm=0; idiag_TRpsidotm=0; idiag_TRdotpsim=0
        idiag_TLeff2m=0; idiag_TLeff2km=0; idiag_TLdoteff2m=0; idiag_TLdoteff2km=0
        idiag_grand2=0; idiag_dgrant=0; idiag_dgrant_up=0; idiag_fact=0
        idiag_grandxy=0; idiag_grantxy=0; idiag_k0=0; idiag_dk=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'Q' ,idiag_Q)
        call parse_name(iname,cname(iname),cform(iname),'Qdot' ,idiag_Qdot)
        call parse_name(iname,cname(iname),cform(iname),'Qddot' ,idiag_Qddot)
        call parse_name(iname,cname(iname),cform(iname),'chi' ,idiag_chi)
        call parse_name(iname,cname(iname),cform(iname),'chidot' ,idiag_chidot)
        call parse_name(iname,cname(iname),cform(iname),'chiddot' ,idiag_chiddot)
        call parse_name(iname,cname(iname),cform(iname),'psi' ,idiag_psi)
        call parse_name(iname,cname(iname),cform(iname),'psiL',idiag_psiL)
        call parse_name(iname,cname(iname),cform(iname),'psidot' ,idiag_psidot)
        call parse_name(iname,cname(iname),cform(iname),'psiddot' ,idiag_psiddot)
        call parse_name(iname,cname(iname),cform(iname),'TR' ,idiag_TR)
        call parse_name(iname,cname(iname),cform(iname),'TL' ,idiag_TL)
        call parse_name(iname,cname(iname),cform(iname),'TRdot' ,idiag_TRdot)
        call parse_name(iname,cname(iname),cform(iname),'TRddot' ,idiag_TRddot)
        call parse_name(iname,cname(iname),cform(iname),'imTR' ,idiag_imTR)
        call parse_name(iname,cname(iname),cform(iname),'psi_anal' ,idiag_psi_anal)
        call parse_name(iname,cname(iname),cform(iname),'TR_anal' ,idiag_TR_anal)
        call parse_name(iname,cname(iname),cform(iname),'TReff2m' ,idiag_TReff2m)
        call parse_name(iname,cname(iname),cform(iname),'TReff2km' ,idiag_TReff2km)
        call parse_name(iname,cname(iname),cform(iname),'TRdoteff2m' ,idiag_TRdoteff2m)
        call parse_name(iname,cname(iname),cform(iname),'TRdoteff2km' ,idiag_TRdoteff2km)
        call parse_name(iname,cname(iname),cform(iname),'TRpsim' ,idiag_TRpsim)
        call parse_name(iname,cname(iname),cform(iname),'TRpsikm' ,idiag_TRpsikm)
        call parse_name(iname,cname(iname),cform(iname),'TRpsidotm' ,idiag_TRpsidotm)
        call parse_name(iname,cname(iname),cform(iname),'TRdotpsim' ,idiag_TRdotpsim)
        call parse_name(iname,cname(iname),cform(iname),'TLeff2m' ,idiag_TLeff2m)
        call parse_name(iname,cname(iname),cform(iname),'TLeff2km' ,idiag_TLeff2km)
        call parse_name(iname,cname(iname),cform(iname),'TLdoteff2m' ,idiag_TLdoteff2m)
        call parse_name(iname,cname(iname),cform(iname),'TLdoteff2km' ,idiag_TLdoteff2km)
        call parse_name(iname,cname(iname),cform(iname),'dgrant_up' ,idiag_dgrant_up)
        call parse_name(iname,cname(iname),cform(iname),'grand2' ,idiag_grand2)
        call parse_name(iname,cname(iname),cform(iname),'dgrant' ,idiag_dgrant)
        call parse_name(iname,cname(iname),cform(iname),'fact' ,idiag_fact)
        call parse_name(iname,cname(iname),cform(iname),'k0' ,idiag_k0)
        call parse_name(iname,cname(iname),cform(iname),'dk' ,idiag_dk)
      enddo
!
!  Check for those quantities for which we want z-averages.
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'grandxy',idiag_grandxy)
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'grantxy',idiag_grantxy)
      enddo
!
    endsubroutine rprint_special
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
