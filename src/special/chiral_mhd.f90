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
! PENCILS PROVIDED del4mu5; del4muS; del6mu5; del6muS
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
   real :: diffmu5=0., diffmuS=0., diffmuSmax=0.
   real :: diffmu5_hyper2=0., diffmuS_hyper2=0.
   real :: diffmu5_hyper3=0., diffmuS_hyper3=0.
   real :: mu5_const=0., gammaf5=0.
   real :: gammaf5_input=0., t1_gammaf5=0., t2_gammaf5=0.
   real :: muS_const=0., coef_muS=0., coef_mu5=0., Cw=0.
   real :: meanmu5=0., flucmu5=0., meanB2=0., Brms=0.
   real :: initpower_mu5=0., cutoff_mu5=0.
   real :: initpower_muS=0., cutoff_muS=0.
   real :: kgaussian_mu5=0.,kpeak_mu5=0.
   real :: kgaussian_muS=0.,kpeak_muS=0.
   real :: radius_mu5=0., sigma_mu5=0., rescale_mu5=1.
   real :: diffmu5_=0., diffmu5_tdep_t0=1., diffmu5_tdep_toffset=1., diffmu5_tdep_exponent=0.
   real, dimension (nx,3) :: aatest, bbtest
   real, dimension (nx,3,3) :: aijtest
   real, pointer :: eta
   real :: cdtchiral=0.4
   real, dimension (nx) :: dt1_lambda5, dt1_D5, dt1_gammaf5
   real, dimension (nx) :: dt1_CMW, dt1_Dmu, dt1_vmu, dt1_special
   real, dimension (nx) :: dt1_CVE1, dt1_CVE2 
   real, dimension (nx) :: uxbj
   integer :: imu5, imuS
   logical :: lmuS=.false., lCVE=.false.
   logical :: ldiffmu5_hyper2_simplified=.false.
   logical :: ldiffmuS_hyper2_simplified=.false.
   logical :: ldiffmu5_hyper3_simplified=.false.
   logical :: ldiffmuS_hyper3_simplified=.false.
   logical :: lremove_mean_mu5=.false., lremove_mean_muS=.false.
   logical :: lmu5adv=.true., lmuSadv=.true.
   logical :: lmu5divu_term=.true., lmuSdivu_term=.true.
   logical :: ldt_chiral_mhd=.true., ldiffmu5_tdep=.false.
   logical :: reinitialize_mu5=.false.
   logical :: lupw_mu5=.false., lupw_muS=.false.
   real :: widthmu5=0.5, widthmuS=0.5
!
  character (len=labellen) :: initspecial='nothing'
  character (len=labellen) :: gammaf5_tdep='const'
!
  namelist /special_init_pars/ &
      initspecial, mu5_const, cdtchiral, &
      lmuS, lCVE, lmu5adv, lmuSadv, muS_const, &
      amplmuS, kx_muS, ky_muS, kz_muS, phase_muS, &
      amplmu5, kx_mu5, ky_mu5, kz_mu5, phase_mu5, &
      coef_muS, coef_mu5, initpower_mu5, cutoff_mu5, &
      initpower_muS, cutoff_muS, lremove_mean_mu5, lremove_mean_muS, &
      kgaussian_mu5, kpeak_mu5, kgaussian_muS, kpeak_muS, &
      radius_mu5, sigma_mu5, widthmu5, widthmuS
!
  namelist /special_run_pars/ &
      diffmu5, diffmuS, diffmuSmax, diffmuSmax, ldt_chiral_mhd, &
      initspecial, mu5_const, &
      diffmu5, diffmuS, diffmuSmax, diffmuSmax, &
      lambda5, cdtchiral, gammaf5, diffmu5_hyper2, diffmuS_hyper2, &
      ldiffmu5_hyper2_simplified, ldiffmuS_hyper2_simplified, &
      diffmu5_hyper3, diffmuS_hyper3, &
      ldiffmu5_hyper3_simplified, ldiffmuS_hyper3_simplified, &
      coef_muS, coef_mu5, Cw, lmuS, lCVE, lmu5adv, &
      lmu5divu_term, lmuSdivu_term, &
      reinitialize_mu5, rescale_mu5, gammaf5_tdep, t1_gammaf5, t2_gammaf5, &
      ldiffmu5_tdep, diffmu5_tdep_toffset, &
      diffmu5_tdep_t0, diffmu5_tdep_exponent, &
      lupw_mu5, lupw_muS
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_muSm=0      ! DIAG_DOC: $\left<\mu_S\right>$
  integer :: idiag_muSrms=0    ! DIAG_DOC: $\left<\mu_S^2\right>^{1/2}$
  integer :: idiag_muSmax=0    ! DIAG_DOC: $\max\mu$
  integer :: idiag_mu5m=0      ! DIAG_DOC: $\left<\mu_5\right>$
  integer :: idiag_mu51m=0     ! DIAG_DOC: $\left<|\mu_5|\right>$
  integer :: idiag_mu53m=0     ! DIAG_DOC: $\left<\mu_5^3\right>$
  integer :: idiag_mu54m=0     ! DIAG_DOC: $\left<\mu_5^4\right>$
  integer :: idiag_mu5rms=0    ! DIAG_DOC: $\left<\mu_5^2\right>^{1/2}$
  integer :: idiag_mu5min=0    ! DIAG_DOC: $\min\mu_5$
  integer :: idiag_mu5max=0    ! DIAG_DOC: $\max\mu_5$
  integer :: idiag_mu5abs=0    ! DIAG_DOC: $\max|\mu_5|$
  integer :: idiag_gamf5m=0    ! DIAG_DOC: $\left<\Gamma_5\right>$
  integer :: idiag_gmu5rms=0   ! DIAG_DOC: $\left<(\nabla\mu_5)^2\right>^{1/2}$     
  integer :: idiag_gmuSrms=0   ! DIAG_DOC: $\left<(\nabla\mu_S)^2\right>^{1/2}$     
  integer :: idiag_gmu5mx=0    ! DIAG_DOC: $\left<\nabla\mu_5\right>_x$   
  integer :: idiag_gmu5my=0    ! DIAG_DOC: $\left<\nabla\mu_5\right>_y$       
  integer :: idiag_gmu5mz=0    ! DIAG_DOC: $\left<\nabla\mu_5\right>_z$   
  integer :: idiag_bgmu5rms=0  ! DIAG_DOC: $\left<(\Bv\cdot\nabla\mu_5)^2\right>^{1/2}$ 
  integer :: idiag_bgmuSrms=0  ! DIAG_DOC: $\left<(\Bv\cdot\nabla\mu_S)^2\right>^{1/2}$ 
  integer :: idiag_mu5bjm=0    ! DIAG_DOC: $\left<\mu_5 ((\nabla\times\Bv)\cdot\Bv) \right>$
  integer :: idiag_mu5bjrms=0  ! DIAG_DOC: $\left<(\mu_5 ((\nabla\times\Bv)\cdot\Bv))^2 \right>^{1/2}$
  integer :: idiag_oogmu5rms=0  
  integer :: idiag_oogmuSrms=0  
  integer :: idiag_dt_lambda5=0! DIAG_DOC: $\mathrm{min}(\mu_5/\Bv^2) \delta x/(\lambda \eta)$ 
  integer :: idiag_dt_D5=0     ! DIAG_DOC: $(\lambda \eta \mathrm{min}(\Bv^2))^{-1}$ 
  integer :: idiag_dt_gammaf5=0! DIAG_DOC: $1/\Gamma_\mathrm{f}$   
  integer :: idiag_dt_CMW=0    ! DIAG_DOC: $\delta x/((C_\mu C_5)^{1/2} \mathrm{max}(|\Bv|))$ 
  integer :: idiag_dt_Dmu=0    ! DIAG_DOC: $(\lambda \eta \mathrm{min}(\Bv^2))^{-1}$ 
  integer :: idiag_dt_vmu=0    ! DIAG_DOC: $\delta x /(\eta \mathrm{max}(|\mu_5 |))$ 
  integer :: idiag_dt_chiral=0 ! DIAG_DOC: total time-step contribution from chiral MHD
  integer :: idiag_mu5bxm=0    ! DIAG_DOC: $\left<\mu_5B_x\right>$
  integer :: idiag_mu5b2m=0    ! DIAG_DOC: $\left<\mu_5B^2\right>$
  integer :: idiag_mu5jbm=0    ! DIAG_DOC: $\left<\mu_5\Jv\cdot\Bv\right>$
  integer :: idiag_jxm = 0     ! DIAG_DOC: $\langle J_x\rangle$
  integer :: idiag_Dmu5_tdep=0 ! DIAG_DOC: $D(t)$
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
      ispecialvar=imu5
!
      if (lmuS) then
        call farray_register_pde('muS',imuS)
        ispecialvar2=imuS
      endif
!
!!      call farray_register_auxiliary('specaux',ispecaux)
!!      call
!      farray_register_auxiliary('specaux',ispecaux,communicated=.true.)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use SharedVariables, only : get_shared_variable
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: ierr
!
!  set gammaf5_input to input value (which was gammaf5)
!
      gammaf5_input=gammaf5
!
!  Reinitialize GW field using a small selection of perturbations
!  that were mostly also available as initial conditions.
!
      if (reinitialize_mu5) then
        select case (initspecial)
!
        case ('rescale')
          f(:,:,:,imu5)=rescale_mu5*f(:,:,:,imu5)
!
        case ('const')
          f(:,:,:,imu5) = mu5_const
          if (lmuS) f(:,:,:,imuS) = muS_const
!
        case default
        endselect
      endif
!
      call keep_compiler_quiet(f)
!
!  give eta out as shared_variable
!
      if (lmagnetic.and.lrun) then
        call get_shared_variable('eta',eta,ierr)
        if (ierr/=0) call fatal_error("initialize_special: ", &
            "cannot get shared var eta")
      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Initcond
      use Sub, only: remove_mean_value, remove_mean, blob
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
        case ('double_sin')
          do n=n1,n2; do m=m1,m2
             f(l1:l2,m,n,imu5)=amplmu5*( sin(kx_mu5*x(l1:l2)) + sin(2.*kx_mu5*x(l1:l2)) )
          enddo; enddo
          if (lmuS) f(:,:,:,imuS) = muS_const
!
        case ('triple_sin')
          do n=n1,n2; do m=m1,m2
             f(l1:l2,m,n,imu5)=amplmu5*( sin(kx_mu5*x(l1:l2)) + sin(2.*kx_mu5*x(l1:l2))  &
              + sin(4.*kx_mu5*x(l1:l2)) )
          enddo; enddo
          if (lmuS) f(:,:,:,imuS) = muS_const
!
        case ('blob')
          call blob(amplmu5,f,imu5,radius_mu5,0.,0.,0.)

        case ('plusminus_sphere')
          call plusminus_sphere(amplmu5,f,imu5,radius_mu5,sigma_mu5,0.,0.,0.)

        case ('gaussian-noise')
          call gaunoise(amplmu5,f,imu5)
          if (lmuS) call gaunoise(amplmuS,f,imuS)
!
        case ('const_sinwave-phase')
            call sinwave_phase(f,imu5,amplmu5,kx_mu5,ky_mu5,kz_mu5,phase_mu5)
            f(:,:,:,imu5) = f(:,:,:,imu5) + mu5_const
          if (lmuS) then
            call sinwave_phase(f,imuS,amplmuS,kx_muS,ky_muS,kz_muS,phase_muS)
            f(:,:,:,imuS) = f(:,:,:,imuS) + muS_const
          endif
!
        case ('mu5const-muSsin')
          f(:,:,:,imu5) = mu5_const
          if (lmuS) call sinwave_phase(f,imuS,amplmuS,kx_muS,ky_muS,kz_muS,phase_muS)
!
        case ('hatwave-x')
          call hatwave(amplmu5,f,imu5,widthmu5,kx=kx_mu5)
          if (lmuS) call hatwave(amplmuS,f,imuS,widthmuS,kx=kx_muS)
!
        case ('power_randomphase')
          call power_randomphase(amplmu5,initpower_mu5,kgaussian_mu5,kpeak_mu5,cutoff_mu5,&
            f,imu5,imu5,lscale_tobox=.false.)
          if(lremove_mean_mu5) call remove_mean(f,imu5)
          if (lmuS) then
            call power_randomphase(amplmuS,initpower_muS,kgaussian_muS,kpeak_muS,cutoff_muS,&
              f,imuS,imuS,lscale_tobox=.false.)
            if(lremove_mean_muS) call remove_mean(f,imuS)
          endif
        case default
          call fatal_error("init_special: No such value for initspecial:" &
              ,trim(initspecial))
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine plusminus_sphere(ampl,f,i,radius,sigma,xblob,yblob,zblob)
!
!  Single sphere with negative outside and positive inside.
!
!  24-feb-21/axel: adapted from blob
!
      integer :: i
      real, dimension (mx,my,mz,mfarray) :: f
      real, optional :: xblob,yblob,zblob
      real :: ampl,radius,sigma,x01=0.,y01=0.,z01=0.,fact
!
!  Single  blob.
!
      if (present(xblob)) x01=xblob
      if (present(yblob)) y01=yblob
      if (present(zblob)) z01=zblob
      if (ampl==0) then
        if (lroot) print*,'ampl=0 in plusminus_sphere'
      else
        if (lroot.and.ip<14) print*,'plusminus_sphere: variable i,ampl=',i,ampl
        fact=1./sigma**2
        f(:,:,:,i)=f(:,:,:,i)-ampl*tanh(fact*(( &
           spread(spread((x-x01)**2,2,my),3,mz) &
          +spread(spread((y-y01)**2,1,mx),3,mz) &
          +spread(spread((z-z01)**2,1,mx),2,my))-radius**2))
      endif
!
    endsubroutine plusminus_sphere
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
      if (diffmu5/=0.) lpenc_requested(i_del2mu5)=.true.
      if (ldiffmu5_hyper2_simplified) lpenc_requested(i_del4mu5)=.true.
      if (ldiffmuS_hyper2_simplified) lpenc_requested(i_del4muS)=.true.
      if (ldiffmu5_hyper3_simplified) lpenc_requested(i_del6mu5)=.true.
      if (ldiffmuS_hyper3_simplified) lpenc_requested(i_del6muS)=.true.
      if (lhydro.or.lhydro_kinematic) then
         lpenc_requested(i_uu)=.true.
         if (lCVE) lpenc_requested(i_oo)=.true.
      endif
      if (lmu5divu_term) lpenc_requested(i_divu)=.true.
      if (lmuSdivu_term) lpenc_requested(i_divu)=.true.
      if (lmagnetic) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_b2)=.true.
      endif
      if (lmagnetic.and.lhydro) lpenc_requested(i_jb)=.true.
!
!  diagnostic pencils
!
      if (idiag_mu5jbm/=0) lpenc_diagnos(i_jb)=.true.
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
      use Sub, only: del2, del4, del6, dot2_mn, del2v_etc, grad, dot
      use Sub, only: u_dot_grad, gij, multsv, curl, curl_mn
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
        if (lpencil(i_ugmuS)) then
          call u_dot_grad(f,imuS,p%gmuS,p%uu,p%ugmuS,UPWIND=lupw_muS)
        endif
        if (lpencil(i_del2muS)) call del2(f,imuS,p%del2muS)
        if (lpencil(i_del4muS)) call del4(f,imuS,p%del4muS)
        if (lpencil(i_del6muS)) call del6(f,imuS,p%del6muS)
      endif
      if (lpencil(i_mu5)) p%mu5=f(l1:l2,m,n,imu5)
      if (lpencil(i_gmu5)) call grad(f,imu5,p%gmu5)
      if (lpencil(i_ugmu5)) then
        call u_dot_grad(f,imu5,p%gmu5,p%uu,p%ugmu5,UPWIND=lupw_mu5)
      endif
      if (lpencil(i_del2mu5)) call del2(f,imu5,p%del2mu5)
      if (lpencil(i_del4mu5)) call del4(f,imu5,p%del4mu5)
      if (lpencil(i_del6mu5)) call del6(f,imu5,p%del6mu5)
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
!  11-jun-20/jenny: added hyperdiffusion (hyper2)
!  28-aug-21/jenny: added hyperdiffusion (hyper3)
!
      use Sub, only: multsv, dot_mn, dot2_mn, dot_mn_vm_trans, dot, curl_mn, gij
      use Diagnostics, only: sum_mn_name, max_mn_name, save_name
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      real, dimension (nx) :: bgmuS, bgmu5, EB, uujj, bbjj, gmu52, bdotgmuS, bdotgmu5
      real, dimension (nx) :: muSmu5, oobb, oogmuS, oogmu5, gmuS2, unity=1.
      real, dimension (nx,3) :: mu5bb, muSmu5oo
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
          +lambda5*EB-gammaf5*p%mu5
!
!  Different diffusion operators.
!
      if (ldiffmu5_hyper2_simplified) then
         df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) &
            -diffmu5_hyper2*p%del4mu5
      else if (ldiffmu5_hyper3_simplified) then
         df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) &
            +diffmu5_hyper3*p%del6mu5
      else
        df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) + diffmu5_*p%del2mu5 
      endif
! 
      if (lmu5adv) then
        df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) - p%ugmu5 
      endif
!
!  Set lmu5divu_term=T to obey total chirality conservation in the compressible case.
!  This is not the default and was only used since Brandenburg (2021, ApJ 911, 110).
!
      if (lmu5divu_term) then
        df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) - p%mu5*p%divu
      endif
!
      if (lCVE) then
        call dot(p%oo,p%gmu5,oogmu5)
        df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) &
          -2.*Cw*p%mu5*oogmu5
      endif
!
!  Contributions to timestep from mu5 equation
      dt1_lambda5 = lambda5*eta*p%b2
      if (lmuS) then
         dt1_CVE2 = p%muS*lambda5*eta*p%b2
      endif
      dt1_D5 = diffmu5*dxyz_2
!      if (lmuS) then
!        dt1_mu5_3 = p%muS*coef_mu5*sqrt(p%b2)
!      endif
      dt1_gammaf5 = gammaf5
!
!  Evolution of muS
!
      if (lmuS) then
        muSmu5 = p%muS*p%mu5
        call dot(p%bb,p%gmu5,bdotgmu5)
        call dot(p%bb,p%gmuS,bdotgmuS)
        df(l1:l2,m,n,imuS) = df(l1:l2,m,n,imuS) &
           -coef_muS*bdotgmu5
! 
        if (lmuSdivu_term) then
          df(l1:l2,m,n,imuS) = df(l1:l2,m,n,imuS) - p%muS*p%divu
        endif
!
        if (ldiffmuS_hyper2_simplified) then
           df(l1:l2,m,n,imuS) = df(l1:l2,m,n,imuS) &
              -diffmuS_hyper2*p%del4muS
        else if (ldiffmuS_hyper3_simplified) then
           df(l1:l2,m,n,imuS) = df(l1:l2,m,n,imuS) &
              +diffmuS_hyper3*p%del6muS
        else
           df(l1:l2,m,n,imuS) = df(l1:l2,m,n,imuS) &
              +diffmuS*p%del2muS
        endif
! 
        if (lmuSadv) then
          df(l1:l2,m,n,imuS) = df(l1:l2,m,n,imuS) - p%ugmuS
        endif
        df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) &
          -coef_mu5*bdotgmuS  
        if (lCVE) then   
          call dot(p%oo,p%bb,oobb)
          call dot(p%oo,p%gmuS,oogmuS)
!          call dot(p%oo,p%gmu5,oogmu5)
          df(l1:l2,m,n,imu5) = df(l1:l2,m,n,imu5) - lambda5*eta*muSmu5*oobb &
            -2.*Cw*p%muS*oogmuS
        endif
!  Contributions to timestep from muS equation
        dt1_CMW = sqrt(coef_mu5*coef_muS)*sqrt(p%b2)*sqrt(dxyz_2)
        dt1_Dmu = diffmuS*dxyz_2
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
      dt1_vmu = eta*p%mu5*sqrt(dxyz_2)
      if (lCVE) then 
         dt1_CVE1 = eta*p%muS*p%mu5*sqrt(dxyz_2)
      endif
!
!  Additions to the test-field equations
!
      if (ltestfield) then
        aatest=f(l1:l2,m,n,iaxtest:iaztest)
        call gij(f,iaxtest,aijtest,1)
        call curl_mn(aijtest,bbtest,aatest)
        df(l1:l2,m,n,iaxtest:iaztest) = df(l1:l2,m,n,iaxtest:iaztest) &
                                        + eta*meanmu5*bbtest
      endif  
!
!  Total contribution to the timestep
!
      if (lfirst.and.ldt.and.ldt_chiral_mhd) then
        if (lmuS) then
          dt1_special = cdtchiral*max(dt1_lambda5, dt1_D5, &
                          dt1_gammaf5, dt1_vmu, &
                          dt1_CVE1, dt1_CVE2, &
                          dt1_CMW, dt1_Dmu) 
        else
          dt1_special = cdtchiral*max(dt1_lambda5, dt1_D5, &
                          dt1_gammaf5, dt1_vmu)
        endif
        dt1_max=max(dt1_max,dt1_special)  
      endif
!
!  Diagnostics
!
      if (ldiagnos) then
        if (idiag_muSm/=0) call sum_mn_name(p%muS,idiag_muSm)
        if (idiag_muSrms/=0) call sum_mn_name(p%muS**2,idiag_muSrms,lsqrt=.true.)
        if (idiag_muSmax/=0) call max_mn_name(p%muS,idiag_muSmax)
        if (idiag_mu5m/=0) call sum_mn_name(p%mu5,idiag_mu5m)
        if (idiag_mu51m/=0) call sum_mn_name(sqrt(p%mu5**2),idiag_mu51m)
        if (idiag_mu53m/=0) call sum_mn_name(p%mu5**3,idiag_mu53m)
        if (idiag_mu54m/=0) call sum_mn_name(p%mu5**4,idiag_mu54m)
        if (idiag_mu5rms/=0) call sum_mn_name(p%mu5**2,idiag_mu5rms,lsqrt=.true.)
        if (idiag_mu5min/=0) call max_mn_name(-p%mu5,idiag_mu5min,lneg=.true.)
        if (idiag_mu5max/=0) call max_mn_name(p%mu5,idiag_mu5max)
        if (idiag_mu5abs/=0) call max_mn_name(abs(p%mu5),idiag_mu5abs)
        if (idiag_gamf5m/=0) call save_name(gammaf5,idiag_gamf5m)
        if (idiag_gmu5rms/=0) then
          call dot2_mn(p%gmu5,gmu52)
          call sum_mn_name(gmu52,idiag_gmu5rms,lsqrt=.true.)
        endif
        if (idiag_gmuSrms/=0) then
          call dot2_mn(p%gmuS,gmuS2)
          call sum_mn_name(gmuS2,idiag_gmuSrms,lsqrt=.true.)
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
        if (idiag_dt_lambda5/=0) call max_mn_name(-(1./dt1_lambda5),idiag_dt_lambda5,lneg=.true.)
        if (idiag_dt_D5/=0) call max_mn_name(-(1./dt1_D5),idiag_dt_D5,lneg=.true.)
        if (idiag_dt_gammaf5/=0) call max_mn_name(-(1./dt1_gammaf5),idiag_dt_gammaf5,lneg=.true.)
        if (idiag_dt_vmu/=0) call max_mn_name(-(1./dt1_vmu),idiag_dt_vmu,lneg=.true.)
        if (idiag_dt_CMW/=0) call max_mn_name(-(1./dt1_CMW),idiag_dt_CMW,lneg=.true.)
        if (idiag_dt_Dmu/=0) call max_mn_name(-(1./dt1_Dmu),idiag_dt_Dmu,lneg=.true.)
        if (idiag_dt_chiral/=0) call max_mn_name(-(1./dt1_special),idiag_dt_chiral,lneg=.true.)
        if (idiag_mu5bxm/=0) call sum_mn_name(p%mu5*p%bb(:,1),idiag_mu5bxm)
        if (idiag_mu5b2m/=0) call sum_mn_name(p%mu5*p%b2,idiag_mu5b2m)
        if (idiag_mu5jbm/=0) call sum_mn_name(p%mu5*p%jb,idiag_mu5jbm)
        if (idiag_Dmu5_tdep/=0) call max_mn_name(unity*diffmu5_,idiag_Dmu5_tdep)
!
!AB: shouldn't this pencil be requested in pencil_criteria_special?
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
      write(unit, NML=special_run_pars)
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
      logical :: lreset,lwrite
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
        idiag_muSm=0; idiag_muSrms=0; idiag_muSmax=0;
        idiag_mu5m=0; idiag_mu51m=0; idiag_mu53m=0; idiag_mu54m=0; idiag_mu5rms=0;
        idiag_mu5min=0; idiag_mu5max=0; idiag_mu5abs=0;
        idiag_gamf5m=0; idiag_gmu5rms=0; idiag_gmuSrms=0; 
        idiag_bgmu5rms=0; idiag_bgmuSrms=0;
        idiag_mu5bjm=0; idiag_mu5bjrms=0;
        idiag_gmu5mx=0; idiag_gmu5my=0; idiag_gmu5mz=0;
        idiag_dt_chiral=0; idiag_dt_vmu=0;
        idiag_dt_lambda5=0; idiag_dt_D5=0;
        idiag_dt_gammaf5=0; idiag_dt_CMW=0; idiag_dt_Dmu=0;
        idiag_jxm=0; idiag_Dmu5_tdep=0; idiag_oogmuSrms=0; idiag_oogmu5rms=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'muSm',idiag_muSm)
        call parse_name(iname,cname(iname),cform(iname),'muSrms',idiag_muSrms)
        call parse_name(iname,cname(iname),cform(iname),'muSmax',idiag_muSmax)
        call parse_name(iname,cname(iname),cform(iname),'mu5m',idiag_mu5m)
        call parse_name(iname,cname(iname),cform(iname),'mu51m',idiag_mu51m)
        call parse_name(iname,cname(iname),cform(iname),'mu53m',idiag_mu53m)
        call parse_name(iname,cname(iname),cform(iname),'mu54m',idiag_mu54m)
        call parse_name(iname,cname(iname),cform(iname),'mu5rms',idiag_mu5rms)
        call parse_name(iname,cname(iname),cform(iname),'mu5min',idiag_mu5min)
        call parse_name(iname,cname(iname),cform(iname),'mu5max',idiag_mu5max)
        call parse_name(iname,cname(iname),cform(iname),'mu5abs',idiag_mu5abs)
        call parse_name(iname,cname(iname),cform(iname),'gamf5m',idiag_gamf5m)
        call parse_name(iname,cname(iname),cform(iname),'gmu5rms',idiag_gmu5rms)
        call parse_name(iname,cname(iname),cform(iname),'gmuSrms',idiag_gmuSrms)
        call parse_name(iname,cname(iname),cform(iname),'gmu5mx',idiag_gmu5mx)
        call parse_name(iname,cname(iname),cform(iname),'gmu5my',idiag_gmu5my)
        call parse_name(iname,cname(iname),cform(iname),'gmu5mz',idiag_gmu5mz)
        call parse_name(iname,cname(iname),cform(iname),'bgmu5rms',idiag_bgmu5rms)
        call parse_name(iname,cname(iname),cform(iname),'bgmuSrms',idiag_bgmuSrms)
        call parse_name(iname,cname(iname),cform(iname),'mu5bjm',idiag_mu5bjm)
        call parse_name(iname,cname(iname),cform(iname),'mu5bjrms',idiag_mu5bjrms)
        call parse_name(iname,cname(iname),cform(iname),'oogmuSrms',idiag_oogmuSrms)
        call parse_name(iname,cname(iname),cform(iname),'oogmu5rms',idiag_oogmu5rms)
        call parse_name(iname,cname(iname),cform(iname),'dt_lambda5',idiag_dt_lambda5)
        call parse_name(iname,cname(iname),cform(iname),'dt_D5',idiag_dt_D5)
        call parse_name(iname,cname(iname),cform(iname),'dt_gammaf5',idiag_dt_gammaf5)
        call parse_name(iname,cname(iname),cform(iname),'dt_CMW',idiag_dt_CMW)
        call parse_name(iname,cname(iname),cform(iname),'dt_Dmu',idiag_dt_Dmu)
        call parse_name(iname,cname(iname),cform(iname),'dt_vmu',idiag_dt_vmu)
        call parse_name(iname,cname(iname),cform(iname),'dt_chiral',idiag_dt_chiral)
        call parse_name(iname,cname(iname),cform(iname),'mu5bxm',idiag_mu5bxm)
        call parse_name(iname,cname(iname),cform(iname),'mu5b2m',idiag_mu5b2m)
        call parse_name(iname,cname(iname),cform(iname),'mu5jbm',idiag_mu5jbm)
        call parse_name(iname,cname(iname),cform(iname),'jxm', idiag_jxm)
        call parse_name(iname,cname(iname),cform(iname),'Dmu5_tdep',idiag_Dmu5_tdep)
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
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  22-aug-21/axel: temporal profile for gammaf5
!
      use Sub, only: remove_mean_value
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
!  Choice of gammaf5_tdep profiles.
!
      select case (gammaf5_tdep)
        case ('const')
          gammaf5=gammaf5_input
          if (headtt.or.ldebug) print*,'gammaf5=const=',gammaf5
!
!  Time-dependent profile for sigma.
!
        case ('step')
          if (t<=t1_gammaf5) then
            gammaf5=0.
          elseif (t<=t2_gammaf5) then
            gammaf5=gammaf5_input
          else
            gammaf5=0.
          endif
!
!  Default.
!
        case default
          call fatal_error("daa_dt: No such value for gammaf5_tdep:",trim(gammaf5_tdep))
      endselect
!
!  The option ldiffmu5_tdep=T allows for a time-dependent diffusivity.
!  This implementation is equivalent to that in viscosity and magnetic if
!  lresi_nu_tdep_t0_norm=T and lresi_eta_tdep=T (which is not the default).
!
      if (ldiffmu5_tdep) then
        diffmu5_=diffmu5*max(real(t-diffmu5_tdep_toffset)/diffmu5_tdep_t0,1.)**diffmu5_tdep_exponent
      else
        diffmu5_=diffmu5
      endif
!
    endsubroutine special_before_boundary
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Calculate meanmu5 and meanB2
!
!  11-oct-15/jenny: coded
!
      use Mpicomm, only: mpiallreduce_sum
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: meanmu5_tmp !, meanB2_tmp
      intent(inout) :: f
!
!  compute meanmu5 and meanB2
!
      meanmu5=0.
!      meanB2=0.
      do n=n1,n2; do m=m1,m2
        meanmu5=meanmu5+sum(f(l1:l2,m,n,imu5))
!        print*, "sum(f(l1:l2,m,n,imu5))", sum(f(l1:l2,m,n,imu5))
      enddo; enddo
!      meanB2=meanB2+sum(p%b2)
!
!  communicate and divide by all mesh meshpoints
!
     if (nprocxy>1) then
   !    call mpiallreduce_sum(meanmu5,meanmu5_tmp,(/nx,ny,nz/))
       call mpiallreduce_sum(meanmu5,meanmu5_tmp)
!       call mpiallreduce_sum(meanB2,meanB2_tmp)
     else
       meanmu5_tmp=meanmu5
     endif
!
! number of grid points
      meanmu5=meanmu5_tmp/nwgrid
!      meanB2=nw1*meanB2_tmp
!      flucmu5=p%mu5-meanmu5
!
    endsubroutine special_after_boundary
!***********************************************************************
!********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING        *************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!*********************************************************************** 
endmodule Special
