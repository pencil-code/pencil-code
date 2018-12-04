! $Id$
!
!  This module provide a way for users to specify custom
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc.
!
!  The module provides a set of standard hooks into the Pencil-Code and
!  currently allows the following customizations:
!
!  Description                                     | Relevant function
!  call
!  ---------------------------------------------------------------------------
!  Special variable registration                   | register_special
!    (pre parameter read)                          |
!  Special variable initialization                 | initialize_special
!    (post parameter read)                         |
!  Special variable finalization                   | finalize_special
!    (deallocation, etc.)                          |
!                                                  |
!  Special initial condition                       | init_special
!   this is called last so may be used to modify   |
!   the mvar variables declared by this module     |
!   or optionally modify any of the other f array  |
!   variables.  The latter, however, should be     |
!   avoided where ever possible.                   |
!                                                  |
!  Special term in the mass (density) equation     |
!  special_calc_density
!  Special term in the momentum (hydro) equation   | special_calc_hydro
!  Special term in the energy equation             | special_calc_energy
!  Special term in the induction (magnetic)        |
!  special_calc_magnetic
!     equation                                     |
!                                                  |
!  Special equation                                | dspecial_dt
!    NOT IMPLEMENTED FULLY YET - HOOKS NOT PLACED INTO THE PENCIL-CODE
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED muS; mu5; gmuS(3); gmu5(3)
! PENCILS PROVIDED ugmu5; ugmuS; del2mu5; del2muS
!***************************************************************
!
! HOW TO USE THIS FILE
! --------------------
!
! Change the line above to
!   lspecial = .true.
! to enable use of special hooks.
!
! The rest of this file may be used as a template for your own
! special module.  Lines which are double commented are intended
! as examples of code.  Simply fill out the prototypes for the
! features you want to use.
!
! Save the file with a meaningful name, eg. geo_kws.f90 and place
! it in the $PENCIL_HOME/src/special directory.  This path has
! been created to allow users ot optionally check their contributions
! in to the Pencil-Code SVN repository.  This may be useful if you
! are working on/using the additional physics with somebodyelse or
! may require some assistance from one of the main Pencil-Code team.
!
! To use your additional physics code edit the Makefile.local in
! the src directory under the run directory in which you wish to
! use your additional physics.  Add a line with all the module
! selections to say something like:
!
!   SPECIAL=special/geo_kws
!
! Where geo_kws it replaced by the filename of your new module
! upto and not including the .f90
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! Declare index of new variables in f array (if any).
!
   real :: amplmuS=0., kx_muS=0., ky_muS=0., kz_muS=0., phase_muS=0.
   real :: amplmu5=0., kx_mu5=0., ky_mu5=0., kz_mu5=0., phase_mu5=0.
   real :: diffmu5, diffmuS, lambda5, mu5_const=0., gammaf5, Cw=0.
   real :: muS_const=0., coef_muS=0., coef_mu5=0.
   real :: meanmu5=0., flucmu5=0.
   real, dimension (nx,3) :: aatest, bbtest
   real, dimension (nx,3,3) :: aijtest
   real, pointer :: eta
   real :: cdtchiral=1.
   real, dimension (nx) :: dt1_mu5_1, dt1_mu5_2, dt1_mu5_3, dt1_mu5_4
   real, dimension (nx) :: dt1_muS_1, dt1_muS_2, dt1_bb_1, dt1_special
   real, dimension (nx) :: uxbj
   integer :: imu5, imuS
   logical :: lmuS=.false., lCVE=.false.
   logical :: lmu5adv=.true., lmuSadv=.true.
!
  character (len=labellen) :: initspecial='nothing'
!
  namelist /special_init_pars/ &
      initspecial, mu5_const, &
      lmuS, lCVE, lmu5adv, lmuSadv, muS_const, coef_muS, &
      amplmuS, kx_muS, ky_muS, kz_muS, phase_muS, &
      amplmu5, kx_mu5, ky_mu5, kz_mu5, phase_mu5, &
      coef_muS, coef_mu5
!
  namelist /special_run_pars/ &
      diffmu5, diffmuS, lambda5, cdtchiral, gammaf5, &
      coef_muS, coef_mu5, Cw, lmuS, lCVE, lmu5adv
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_muSm=0      ! DIAG_DOC: $\left<\mu_S\right>$
  integer :: idiag_muSrms=0    ! DIAG_DOC: $\left<\mu_S^2\right>^{1/2}$
  integer :: idiag_mu5m=0      ! DIAG_DOC: $\left<\mu_5\right>$
  integer :: idiag_mu5rms=0    ! DIAG_DOC: $\left<\mu_5^2\right>^{1/2}$
  integer :: idiag_gmu5rms=0   ! DIAG_DOC: $\left<(\nabla\mu_5)^2\right>^{1/2}$     
  integer :: idiag_gmu5mx=0    ! DIAG_DOC: $\left<\nabla\mu_5\right>_x$   
  integer :: idiag_gmu5my=0    ! DIAG_DOC: $\left<\nabla\mu_5\right>_y$       
  integer :: idiag_gmu5mz=0    ! DIAG_DOC: $\left<\nabla\mu_5\right>_z$   
  integer :: idiag_bgmu5rms=0  ! DIAG_DOC: $\left<(\Bv\cdot\nabla\mu_5)^2\right>^{1/2}$ 
  integer :: idiag_bgmuSrms=0  ! DIAG_DOC: $\left<(\Bv\cdot\nabla\mu_S)^2\right>^{1/2}$ 
  integer :: idiag_mu5bjm=0    ! DIAG_DOC: $\left<\mu_5 ((\nabla\times\Bv)\cdot\Bv) \right>$
  integer :: idiag_mu5bjrms=0  ! DIAG_DOC: $\left<(\mu_5 ((\nabla\times\Bv)\cdot\Bv))^2 \right>^{1/2}$
  integer :: idiag_oogmu5rms=0  
  integer :: idiag_oogmuSrms=0  
  integer :: idiag_dt_mu5_1=0  ! DIAG_DOC: $\mathrm{min}(\mu_5/\Bv^2) \delta x/(\lambda \eta)$ 
  integer :: idiag_dt_mu5_2=0  ! DIAG_DOC: $(\lambda \eta \mathrm{min}(\Bv^2))^{-1}$ 
  integer :: idiag_dt_mu5_3=0  ! DIAG_DOC: $\delta x^2/D_5$   
  integer :: idiag_dt_mu5_4=0  ! DIAG_DOC: $\delta x^2/D_5$   
  integer :: idiag_dt_muS_1=0  ! DIAG_DOC: $\mathrm{min}(\mu_5/\Bv^2) \delta x/(\lambda \eta)$ 
  integer :: idiag_dt_muS_2=0  ! DIAG_DOC: $(\lambda \eta \mathrm{min}(\Bv^2))^{-1}$ 
  integer :: idiag_dt_bb_1=0   ! DIAG_DOC: $\delta x /(\eta \mathrm{max}(\mu_5))$ 
  integer :: idiag_dt_chiral=0 ! DIAG_DOC: total time-step contribution from chiral MHD
  integer :: idiag_mu5bxm=0    ! DIAG_DOC: $\left<\mu_5B_x\right>$
  integer :: idiag_mu5b2m=0    ! DIAG_DOC: $\left<\mu_5B^2\right>$
  integer :: idiag_jxm = 0     ! DIAG_DOC: $\langle J_x\rangle$
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      use FArrayManager, only: farray_register_pde
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_pde('mu5',imu5)
!
      if (lmuS) then
        call farray_register_pde('muS',imuS)
      endif
!
!!      call farray_register_auxiliary('specaux',ispecaux)
!!      call
!farray_register_auxiliary('specaux',ispecaux,communicated=.true.)
!
    endsubroutine register_special
!***********************************************************************
    subroutine register_particles_special(npvar)
!
!  Set up indices for particle variables in special modules.
!
!  4-jan-14/tony: coded
!
      integer :: npvar
!
      if (lroot) call svn_id( &
           "$Id$")
      call keep_compiler_quiet(npvar)
!
!
!!      iqp=npvar+1
!!      npvar=npvar+1
!
    endsubroutine register_particles_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use SharedVariables, only : get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ierr
!
      call keep_compiler_quiet(f)
!
      if (lmagnetic.and.lrun) then
          call get_shared_variable('eta',eta,ierr)
          if (ierr/=0) call fatal_error("initialize_special: ", &
              "cannot get shared var eta")
      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine finalize_special(f)
!
!  Called right before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine finalize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Initcond
!
      real, dimension (mx,my,mz,mfarray) :: f,df
!
      intent(inout) :: f
!
!  initial conditions
!
      select case (initspecial)
!
        case ('nothing'); if (lroot) print*,'init_special: nothing'
!
        case ('zero')
          f(:,:,:,imu5) = 0.
          if (lmuS) f(:,:,:,imuS) = 0.
!
        case ('const')
          f(:,:,:,imu5) = mu5_const
          if (lmuS) f(:,:,:,imuS) = muS_const
!
        case ('sinwave-phase')
          call sinwave_phase(f,imu5,amplmu5,kx_mu5,ky_mu5,kz_mu5,phase_mu5)
          if (lmuS) call sinwave_phase(f,imuS,amplmuS,kx_muS,ky_muS,kz_muS,phase_muS)
!
        case ('mu5const-muSsin')
          f(:,:,:,imu5) = mu5_const
          if (lmuS) call sinwave_phase(f,imuS,amplmuS,kx_muS,ky_muS,kz_muS,phase_muS)

!
        case default
          call fatal_error("init_special: No such value for initspecial:" &
              ,trim(initspecial))
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (lmuS) then
        lpenc_requested(i_muS)=.true.
        lpenc_requested(i_gmuS)=.true.
        lpenc_requested(i_ugmuS)=.true.
        if (diffmuS/=0.) lpenc_requested(i_del2muS)=.true.
      endif
      lpenc_requested(i_mu5)=.true.
      lpenc_requested(i_gmu5)=.true.
      lpenc_requested(i_ugmu5)=.true.
      if (ldt) lpenc_requested(i_rho1)=.true.
!      lpenc_requested(i_jjij)=.true.
      if (diffmu5/=0.) lpenc_requested(i_del2mu5)=.true.
      if (lhydro.or.lhydro_kinematic) lpenc_requested(i_uu)=.true.
      if (lmagnetic) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_b2)=.true.
      endif
!      if (lmagnetic) lpenc_requested(i_jij)=.true.
!      if (lmagnetic.and.lhydro) lpenc_requested(i_ub)=.true.
      if (lmagnetic.and.lhydro) lpenc_requested(i_jb)=.true.
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified
!  here.
!
!  18-07-06/tony: coded
!
      logical, dimension(npencils), intent(inout) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      use Sub, only: del2, dot2_mn, del2v_etc, grad, dot, u_dot_grad, gij
      use Sub, only: multsv, curl, curl_mn
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
      if (lmuS) then
        if (lpencil(i_muS)) p%muS=f(l1:l2,m,n,imuS)
        if (lpencil(i_gmuS)) call grad(f,imuS,p%gmuS)
        if (lpencil(i_ugmuS)) call dot(p%uu,p%gmuS,p%ugmuS)
        if (lpencil(i_del2muS)) call del2(f,imuS,p%del2muS)
      endif
      if (lpencil(i_mu5)) p%mu5=f(l1:l2,m,n,imu5)
      if (lpencil(i_gmu5)) call grad(f,imu5,p%gmu5)
      if (lpencil(i_ugmu5)) call dot(p%uu,p%gmu5,p%ugmu5)
      if (lpencil(i_del2mu5)) call del2(f,imu5,p%del2mu5)
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
!  Several precalculated Pencils of information are passed for
!  efficiency.
!
!  06-oct-03/tony: coded
!  29-sep-18/axel: included ldiffus_mu5_1_old and modified diffus_mu5_1
!  25-noc-18/jenny: included muS terms in timestep calculation
!  25-noc-18/jenny: added diffusion term to muS equation
!
      use Sub, only: multsv, dot_mn, dot2_mn, dot_mn_vm_trans, dot, curl_mn, gij
      use Diagnostics, only: sum_mn_name, max_mn_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      real, dimension (nx) :: bgmuS, bgmu5, EB, uujj, bbjj, gmu52, bdotgmuS, bdotgmu5
      real, dimension (nx) :: muSmu5, oobb, oogmuS, oogmu5
      real, dimension (nx,3) :: mu5bb, muSmu5oo
      real, parameter :: alpha_fine_structure=1./137.
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!!    if (headtt) call identify_bcs('mu5',imu5)
!
!  Compute E.B
!
      EB=eta*(p%jb-p%mu5*p%b2)
!
!  Evolution of mu5
!
      df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) &
        +diffmu5*p%del2mu5+lambda5*EB-gammaf5*p%mu5 
      if (lmu5adv) then
       df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) - p%ugmu5 
      endif
!  Contributions to timestep from mu5 equation
      dt1_mu5_1 = lambda5*eta*p%b2
      dt1_mu5_2 = diffmu5*dxyz_2
      if (lmuS) then
        dt1_mu5_3 = p%muS*coef_muS*sqrt(p%b2)
      endif
      dt1_mu5_4 = gammaf5
!
!  Evolution of muS
!
      if (lmuS) then
        muSmu5 = p%muS*p%mu5
        call dot(p%bb,p%gmu5,bdotgmu5)
        call dot(p%bb,p%gmuS,bdotgmuS)
        df(l1:l2,m,n,imuS) = df(l1:l2,m,n,imuS) &
          + diffmuS*p%del2muS - coef_muS*bdotgmu5
        if (lmuSadv) then
          df(l1:l2,m,n,imuS) = df(l1:l2,m,n,imuS) - p%ugmuS
        endif
        df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) &
          -coef_mu5*bdotgmuS  
        if (lCVE) then   
          call dot(p%oo,p%bb,oobb)
          call dot(p%oo,p%gmuS,oogmuS)
          call dot(p%oo,p%gmu5,oogmu5)
          df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) - lambda5*eta*muSmu5*oobb &
            -2.*Cw*(p%muS*oogmuS+p%mu5*oogmu5)
        endif
!  Contributions to timestep from muS equation
        dt1_muS_1 = p%mu5*coef_muS*sqrt(p%b2)
        dt1_muS_2 = diffmuS*dxyz_2
      endif
!                          
!  Additions to evolution of bb
!
      if (lmagnetic) then
        call multsv(p%mu5,p%bb,mu5bb)
        df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + eta*mu5bb 
        if (lCVE) then   
          call multsv(muSmu5,p%oo,muSmu5oo)
          df(l1:l2,m,n,iax:iaz) = df(l1:l2,m,n,iax:iaz) + eta*muSmu5oo
        endif
      endif
!  Contributions to timestep from bb equation
      dt1_bb_1 = eta*p%mu5*sqrt(dxyz_2)
!
!  Additions to evolution of uu
!
      if (lhydro) then
        df(l1:l2,m,n,iux:iuz) = df(l1:l2,m,n,iux:iuz)
      endif  
!
!  Additions to the test-field equations
      if (ltestfield) then
        aatest=f(l1:l2,m,n,iaxtest:iaztest)
        call gij(f,iaxtest,aijtest,1)
        call curl_mn(aijtest,bbtest,aatest)
        df(l1:l2,m,n,iaxtest:iaztest) = df(l1:l2,m,n,iaxtest:iaztest) &
                                        + eta*meanmu5*bbtest
      endif  
!
!  Todal contribution to the timestep
!
      if (lfirst.and.ldt) then
        if (lmuS) then
          dt1_special = cdtchiral*max(dt1_mu5_1, dt1_mu5_2, &
                          dt1_mu5_3, dt1_mu5_4, dt1_bb_1, &
                          dt1_muS_1, dt1_muS_2) 
        else
          dt1_special = cdtchiral*max(dt1_mu5_1, dt1_mu5_2, &
                          dt1_mu5_3, dt1_mu5_4, dt1_bb_1)
        endif
        maxdiffus=max(maxdiffus,dt1_special)
      endif
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_muSm/=0) call sum_mn_name(p%muS,idiag_muSm)
        if (idiag_muSrms/=0) call sum_mn_name(p%muS**2,idiag_muSrms,lsqrt=.true.)
        if (idiag_mu5m/=0) call sum_mn_name(p%mu5,idiag_mu5m)
        if (idiag_mu5rms/=0) call sum_mn_name(p%mu5**2,idiag_mu5rms,lsqrt=.true.)
        if (idiag_gmu5rms/=0) then
          call dot2_mn(p%gmu5,gmu52)
          call sum_mn_name(gmu52,idiag_gmu5rms,lsqrt=.true.)
        endif
        if (idiag_gmu5mx/=0) call sum_mn_name(p%gmu5(:,1),idiag_gmu5mx)
        if (idiag_gmu5my/=0) call sum_mn_name(p%gmu5(:,2),idiag_gmu5my)
        if (idiag_gmu5mz/=0) call sum_mn_name(p%gmu5(:,3),idiag_gmu5mz)
        if (idiag_bgmu5rms/=0) then
          call dot_mn(p%bb,p%gmu5,bgmu5)
          call sum_mn_name(bgmu5**2,idiag_bgmu5rms,lsqrt=.true.)
        endif
        if (idiag_bgmuSrms/=0) then
          call dot_mn(p%bb,p%gmuS,bgmuS)
          call sum_mn_name(bgmuS**2,idiag_bgmuSrms,lsqrt=.true.)
        endif
        if (idiag_mu5bjm/=0) then
          call dot_mn(p%bb,p%jj,bbjj)
          call sum_mn_name(p%mu5*bbjj,idiag_mu5bjm)
        endif
        if (idiag_mu5bjrms/=0) then
          call dot_mn(p%bb,p%jj,bbjj)
          call sum_mn_name((p%mu5*bbjj)**2,idiag_mu5bjrms,lsqrt=.true.)
        endif
        if (idiag_oogmu5rms/=0) call sum_mn_name(oogmu5**2,idiag_oogmu5rms,lsqrt=.true.)
        if (idiag_oogmuSrms/=0) call sum_mn_name(oogmuS**2,idiag_oogmuSrms,lsqrt=.true.)
        if (idiag_dt_mu5_1/=0) call max_mn_name(-(1./dt1_mu5_1),idiag_dt_mu5_1,lneg=.true.)
        if (idiag_dt_mu5_2/=0) call max_mn_name(-(1./dt1_mu5_2),idiag_dt_mu5_2,lneg=.true.)
        if (idiag_dt_mu5_3/=0) call max_mn_name(-(1./dt1_mu5_3),idiag_dt_mu5_3,lneg=.true.)
        if (idiag_dt_mu5_4/=0) call max_mn_name(-(1./dt1_mu5_4),idiag_dt_mu5_4,lneg=.true.)
        if (idiag_dt_bb_1/=0) call max_mn_name(-(1./dt1_bb_1),idiag_dt_bb_1,lneg=.true.)
        if (idiag_dt_muS_1/=0) call max_mn_name(-(1./dt1_muS_1),idiag_dt_muS_1,lneg=.true.)
        if (idiag_dt_muS_2/=0) call max_mn_name(-(1./dt1_muS_2),idiag_dt_muS_2,lneg=.true.)
        if (idiag_dt_chiral/=0) call max_mn_name(-(1./dt1_special),idiag_dt_chiral,lneg=.true.)
        if (idiag_mu5bxm/=0) call sum_mn_name(p%mu5*p%bb(:,1),idiag_mu5bxm)
        if (idiag_mu5b2m/=0) call sum_mn_name(p%mu5*p%b2,idiag_mu5b2m)
        if (idiag_jxm /= 0) then
          lpenc_diagnos(i_jj) = .true.
          call sum_mn_name(p%jj(:,1), idiag_jxm)
        endif
     endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=special_init_pars, IOSTAT=iostat)
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
      call keep_compiler_quiet(unit)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset
      logical, optional :: lwrite
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='muS'.or.cnamev=='mu5') cformv='DEFINED'
      endif
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_muSm=0; idiag_muSrms=0
        idiag_mu5m=0; idiag_mu5rms=0; 
        idiag_bgmu5rms=0; idiag_bgmuSrms=0;
        idiag_mu5bjm=0; idiag_mu5bjrms=0; idiag_gmu5rms=0;
        idiag_gmu5mx=0; idiag_gmu5my=0; idiag_gmu5mz=0;
        idiag_dt_chiral=0; idiag_dt_bb_1=0;
        idiag_dt_mu5_1=0; idiag_dt_mu5_2=0; idiag_dt_mu5_3=0;
        idiag_dt_mu5_4=0; idiag_dt_muS_1=0; idiag_dt_muS_2=0;
        idiag_jxm=0; idiag_oogmuSrms=0; idiag_oogmu5rms=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'muSm',idiag_muSm)
        call parse_name(iname,cname(iname),cform(iname),'muSrms',idiag_muSrms)
        call parse_name(iname,cname(iname),cform(iname),'mu5m',idiag_mu5m)
        call parse_name(iname,cname(iname),cform(iname),'mu5rms',idiag_mu5rms)
        call parse_name(iname,cname(iname),cform(iname),'gmu5rms',idiag_gmu5rms)
        call parse_name(iname,cname(iname),cform(iname),'gmu5mx',idiag_gmu5mx)
        call parse_name(iname,cname(iname),cform(iname),'gmu5my',idiag_gmu5my)
        call parse_name(iname,cname(iname),cform(iname),'gmu5mz',idiag_gmu5mz)
        call parse_name(iname,cname(iname),cform(iname),'bgmu5rms',idiag_bgmu5rms)
        call parse_name(iname,cname(iname),cform(iname),'bgmuSrms',idiag_bgmuSrms)
        call parse_name(iname,cname(iname),cform(iname),'mu5bjm',idiag_mu5bjm)
        call parse_name(iname,cname(iname),cform(iname),'mu5bjrms',idiag_mu5bjrms)
        call parse_name(iname,cname(iname),cform(iname),'oogmuSrms',idiag_oogmuSrms)
        call parse_name(iname,cname(iname),cform(iname),'oogmu5rms',idiag_oogmu5rms)
        call parse_name(iname,cname(iname),cform(iname),'dt_mu5_1',idiag_dt_mu5_1)
        call parse_name(iname,cname(iname),cform(iname),'dt_mu5_2',idiag_dt_mu5_2)
        call parse_name(iname,cname(iname),cform(iname),'dt_mu5_3',idiag_dt_mu5_3)
        call parse_name(iname,cname(iname),cform(iname),'dt_mu5_4',idiag_dt_mu5_4)
        call parse_name(iname,cname(iname),cform(iname),'dt_muS_1',idiag_dt_muS_1)
        call parse_name(iname,cname(iname),cform(iname),'dt_muS_2',idiag_dt_muS_2)
        call parse_name(iname,cname(iname),cform(iname),'dt_bb_1',idiag_dt_bb_1)
        call parse_name(iname,cname(iname),cform(iname),'dt_chiral',idiag_dt_chiral)
        call parse_name(iname,cname(iname),cform(iname),'mu5bxm',idiag_mu5bxm)
        call parse_name(iname,cname(iname),cform(iname),'mu5b2m',idiag_mu5b2m)
        call parse_name(iname, cname(iname), cform(iname), 'jxm', idiag_jxm)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
!
!   1-oct-18/axel: adapted from sample
!
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices.
!
      select case (trim(slices%name))
        case ('muS'); call assign_slices_scal(slices,f,imuS)
        case ('mu5'); call assign_slices_scal(slices,f,imu5)
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Calculate meanmu5
!
!  11-oct-15/jenny: coded
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: fact, meanmu5_tmp, nw1
      intent(inout) :: f
!
!  compute meanmu5
!
      meanmu5=0.
      do n=n1,n2; do m=m1,m2
        meanmu5=meanmu5+sum(f(l1:l2,m,n,imu5))
!        print*, "sum(f(l1:l2,m,n,imu5))", sum(f(l1:l2,m,n,imu5))
      enddo; enddo
!
!  communicate and divide by all mesh meshpoints
!
     if (nprocxy>1) then
   !    call mpiallreduce_sum(meanmu5,meanmu5_tmp,(/nx,ny,nz/))
       call mpiallreduce_sum(meanmu5,meanmu5_tmp)
     endif
!     fact=1./ncpus
!
! number of grid points
      nw1=1./(nxgrid*nygrid*nzgrid)
!
      meanmu5=nw1*meanmu5_tmp
!      flucmu5=p%mu5-meanmu5
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine special_calc_hydro(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  momentum equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_hydro
!***********************************************************************
    subroutine special_calc_density(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  continuity equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_density
!***********************************************************************
    subroutine special_calc_dustdensity(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  continuity equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_dustdensity
!***********************************************************************
    subroutine special_calc_energy(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  energy equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ient) = df(l1:l2,m,n,ient) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_energy
!***********************************************************************
    subroutine special_calc_magnetic(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  induction equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_magnetic
!***********************************************************************
    subroutine special_calc_pscalar(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  passive scalar equation.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array.
!
!  15-jun-09/anders: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,ilncc) = df(l1:l2,m,n,ilncc) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_pscalar
!***********************************************************************
    subroutine special_particles_bfre_bdary(f,fp,ineargrid)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (:,:), intent(in) :: fp
      integer, dimension(:,:) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine special_particles_bfre_bdary
!***********************************************************************
    subroutine special_calc_particles(f,fp,ineargrid)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (:,:), intent(in) :: fp
      integer, dimension(:,:) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine special_calc_particles
!***********************************************************************
    subroutine special_calc_chemistry(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  induction equation.
!
!
!  15-sep-10/natalia: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(inout) :: df
      type (pencil_case), intent(in) :: p
!!
!!  SAMPLE IMPLEMENTATION (remember one must ALWAYS add to df).
!!
!!  df(l1:l2,m,n,iux) = df(l1:l2,m,n,iux) + SOME NEW TERM
!!  df(l1:l2,m,n,iuy) = df(l1:l2,m,n,iuy) + SOME NEW TERM
!!  df(l1:l2,m,n,iuz) = df(l1:l2,m,n,iuz) + SOME NEW TERM
!!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine special_calc_chemistry
!***********************************************************************
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_boundconds(f,bc)
!
!  Some precalculated pencils of data are passed in for efficiency,
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      type (boundary_condition), intent(in) :: bc
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(bc)
!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_,llast)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      logical, intent(in) :: llast
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(dt_)
      call keep_compiler_quiet(llast)
!
    endsubroutine  special_after_timestep
!***********************************************************************
    subroutine set_init_parameters(Ntot,dsize,init_distr,init_distr2)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(ndustspec) :: dsize,init_distr,init_distr2
      real :: Ntot
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(dsize,init_distr,init_distr2)
      call keep_compiler_quiet(Ntot)
!
    endsubroutine  set_init_parameters
!***********************************************************************
  subroutine initialize_mult_special
          
  endsubroutine initialize_mult_special

!***********************************************************************
  subroutine finalize_mult_special
        

  endsubroutine finalize_mult_special
!*********************************************************************** 
endmodule Special
