! $Id$
!
!  This module supplies a kinematic velocity field.
!  Most of the content of this module was moved with revision r12019
!  by Dhrubaditya Mitra on 5-nov-09 away from nohydro.f90
!  To inspect the revision history of the file before that time,
!  check out nohydro.f90 prior to or at revision r12018.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module.
!
! CPARAM logical, parameter :: lhydro = .false.
! CPARAM logical, parameter :: lhydro_kinematic = .true.
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED oo(3); ou; uij(3,3); uu(3); u2; sij(3,3)
! PENCILS PROVIDED der6u(3)
! PENCILS PROVIDED divu; uij5(3,3); graddivu(3)
!***********************************************************************
module Hydro
!
  use Cparam
  use Cdata
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include 'record_types.h'
  include 'hydro.h'
!
  real, dimension (mz,3) :: uumz=0.
  real, dimension (mz,3) :: uumzg=0.,guumz=0.
  real, dimension (mx,my,3) :: uumxy=0.
!
  real :: u_out_kep=0.0
  real :: tphase_kinflow=-1.,phase1=0., phase2=0., tsforce=0.
  real ::  dtforce=impossible
  real, dimension(3) :: location,location_fixed=(/0.,0.,0./)
  logical :: lcalc_uumean=.false.,lupw_uu=.false.
  logical :: lcalc_uumeanxy=.false.
!
  real, allocatable, dimension (:,:) :: KS_k,KS_A,KS_B !or through whole field for each wavenumber?
  real, allocatable, dimension (:) :: KS_omega !or through whole field for each wavenumber?
  integer :: KS_modes = 25
  real, allocatable, dimension (:) :: Zl,dZldr,Pl,dPldtheta
  real :: ampl_fcont_uu=1.
  logical :: lforcing_cont_uu=.false., lrandom_location=.false., lwrite_random_location=.false.
  real, dimension(nx) :: ck_r,ck_rsqr
!
!  init parameters
!  (none)
!
!  run parameters
!
  character (len=labellen) :: kinematic_flow='none'
  real :: ABC_A=1.0, ABC_B=1.0, ABC_C=1.0
  real :: wind_amp=0.,wind_rmin=impossible,wind_step_width=0.
  real :: circ_amp=0.,circ_rmax=0.,circ_step_width=0.
  real :: kx_uukin=1., ky_uukin=1., kz_uukin=1.
  real :: cx_uukin=0., cy_uukin=0., cz_uukin=0.
  real :: phasez_uukin=0.
  real :: radial_shear=0.,uphi_at_rzero=0.,uphi_at_rmax=0.,uphi_rmax=1.,&
          uphi_step_width=0.
  real :: gcs_rzero=0.,gcs_psizero=0.
  real :: kinflow_ck_Balpha=0.
  real :: eps_kinflow=0., omega_kinflow=0., ampl_kinflow=1.
  integer :: kinflow_ck_ell=0.
  character (len=labellen) :: wind_profile='none'
  logical, target :: lpressuregradient_gas=.false.
!
  namelist /hydro_run_pars/ &
      kinematic_flow,wind_amp,wind_profile,wind_rmin,wind_step_width, &
      circ_rmax,circ_step_width,circ_amp, ABC_A,ABC_B,ABC_C, &
      ampl_kinflow, &
      kx_uukin,ky_uukin,kz_uukin, &
      cx_uukin,cy_uukin,cz_uukin, &
      phasez_uukin, &
      lrandom_location,lwrite_random_location,location_fixed,dtforce, &
      radial_shear,uphi_at_rzero,uphi_rmax,uphi_step_width,gcs_rzero, &
      gcs_psizero,kinflow_ck_Balpha,kinflow_ck_ell, &
      eps_kinflow,omega_kinflow,ampl_kinflow
!
  integer :: idiag_u2m=0,idiag_um2=0,idiag_oum=0,idiag_o2m=0
  integer :: idiag_uxpt=0,idiag_uypt=0,idiag_uzpt=0
  integer :: idiag_dtu=0,idiag_urms=0,idiag_umax=0,idiag_uzrms=0
  integer :: idiag_uzmax=0,idiag_orms=0,idiag_omax=0
  integer :: idiag_ux2m=0,idiag_uy2m=0,idiag_uz2m=0
  integer :: idiag_uxuym=0,idiag_uxuzm=0,idiag_uyuzm=0,idiag_oumphi=0
  integer :: idiag_ruxm=0,idiag_ruym=0,idiag_ruzm=0,idiag_rumax=0
  integer :: idiag_uxmz=0,idiag_uymz=0,idiag_uzmz=0,idiag_umx=0
  integer :: idiag_umy=0,idiag_umz=0,idiag_uxmxy=0,idiag_uymxy=0,idiag_uzmxy=0
  integer :: idiag_Marms=0,idiag_Mamax=0,idiag_divu2m=0,idiag_epsK=0
  integer :: idiag_urmphi=0,idiag_upmphi=0,idiag_uzmphi=0,idiag_u2mphi=0
  integer :: idiag_phase1=0,idiag_phase2=0
  integer :: idiag_ekintot=0,idiag_ekin=0
  integer :: idiag_divum=0
!
  contains
!***********************************************************************
    subroutine register_hydro()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!   6-nov-01/wolf: coded
!
      use Mpicomm, only: lroot
      use SharedVariables, only: put_shared_variable
!
      integer :: ierr
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
          "$Id$")
!
      call put_shared_variable('lpressuregradient_gas',&
          lpressuregradient_gas,ierr)
      if (ierr/=0) call fatal_error('register_hydro',&
          'there was a problem sharing lpressuregradient_gas')
!
    endsubroutine register_hydro
!***********************************************************************
    subroutine initialize_hydro(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  24-nov-02/tony: coded
!
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      kinflow=kinematic_flow
      if (kinflow=='KS') then
        call periodic_KS_setup(-5./3.) !Kolmogorov spec. periodic KS
        !call random_isotropic_KS_setup(-5./3.,1.,(nxgrid)/2.) !old form
        !call random_isotropic_KS_setup_test !Test KS model code with 3 specific modes.
        elseif (kinflow=='ck') then
          call init_ck
      endif
!
!  Register an extra aux slot for uu if requested (so uu is written
!  to snapshots and can be easily analyzed later). For this to work you
!  must reserve enough auxiliary workspace by setting, for example,
!     ! MAUX CONTRIBUTION 3
!  in the beginning of your src/cparam.local file, *before* setting
!  ncpus, nprocy, etc.
!
!  After a reload, we need to rewrite index.pro, but the auxiliary
!  arrays are already allocated and must not be allocated again.
!
      if (lkinflow_as_aux) then
        if (iuu==0) then
          call farray_register_auxiliary('uu',iuu,vector=3)
          iux=iuu
          iuy=iuu+1
          iuz=iuu+2
        endif
! set the initial velocity to zero
        f(:,:,:,iux:iuz) = 0.
        if (iuu/=0.and.lroot) then
          print*, 'initialize_velocity: iuu = ', iuu
          open(3,file=trim(datadir)//'/index.pro', POSITION='append')
          write(3,*) 'iuu=',iuu
          write(3,*) 'iux=',iux
          write(3,*) 'iuy=',iuy
          write(3,*) 'iuz=',iuz
          close(3)
        endif
      endif
!
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_hydro
!***********************************************************************
    subroutine init_uu(f)
!
!  Initialise uu.
!
!   7-jun-02/axel: adapted from hydro
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_uu
!***********************************************************************
    subroutine pencil_criteria_hydro()
!
!  All pencils that the Hydro module depends on are specified here.
!
!  20-nov-04/anders: coded
!   1-jul-09/axel: added more for kinflow
!
!  pencils for kinflow
!
!DM: The following line with kinflow can be later removed and the variable
!DM: kinematic_flow replaced by kinflow.
!
      kinflow=kinematic_flow
      if (kinflow/='') then
        lpenc_requested(i_uu)=.true.
        if (kinflow=='eddy') then
          lpenc_requested(i_rcyl_mn)=.true.
          lpenc_requested(i_rcyl_mn1)=.true.
        endif
      endif
!
!  Diagnostic pencils.
!
      if (idiag_urms/=0 .or. idiag_umax/=0 .or. idiag_u2m/=0 .or. &
          idiag_um2/=0) lpenc_diagnos(i_u2)=.true.
      if (idiag_oum/=0) lpenc_diagnos(i_ou)=.true.
      if (idiag_divum/=0) lpenc_diagnos(i_divu)=.true.
!
    endsubroutine pencil_criteria_hydro
!***********************************************************************
    subroutine pencil_interdep_hydro(lpencil_in)
!
!  Interdependency among pencils from the Hydro module is specified here
!
!  20-nov-04/anders: coded
!
      logical, dimension (npencils) :: lpencil_in
!
      if (lpencil_in(i_uglnrho)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_ugrho)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_grho)=.true.
      endif
      if (lpencil_in(i_uij5glnrho)) then
        lpencil_in(i_uij5)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_u2)) lpencil_in(i_uu)=.true.
! oo
      if (lpencil_in(i_ou)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_oo)=.true.
      endif
!
    endsubroutine pencil_interdep_hydro
!***********************************************************************
    subroutine calc_pencils_hydro(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   08-nov-04/tony: coded
!
      use Diagnostics
      use General
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension(nx) :: kdotxwt, cos_kdotxwt, sin_kdotxwt
      real, dimension(nx) :: local_Omega
      real, dimension(nx) :: wind_prof,div_uprof,der6_uprof
      real, dimension(nx) :: div_vel_prof
      real, dimension(nx) :: vel_prof
      real, dimension(nx) :: tmp_mn, cos1_mn, cos2_mn
      real, dimension(nx) :: rone, argx
      real :: fac, fac2, argy, argz, cxt, cyt, czt, omt
      real :: fpara, dfpara, ecost, esint, epst, sin2t, cos2t
      real :: sqrt2, sqrt21k1, eps1=1., WW=0.25, k21
      real :: Balpha
      real :: theta,theta1
      integer :: modeN, ell
!
      intent(in) :: f
      intent(inout) :: p
!
!  Choose from a list of different flow profiles.
!  Begin with a 
!
      select case (kinematic_flow)
!
!constant flow in the x direction.
!
      case ('const-x')
        if (headtt) print*,'const-x'
        if (lpencil(i_uu)) then
          p%uu(:,1)=ampl_kinflow*cos(omega_kinflow*t)*exp(eps_kinflow*t)
          p%uu(:,2)=0.
          p%uu(:,3)=0.
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!constant flow in the
!  (ampl_kinflow_x,ampl_kinflow_y,ampl_kinflow_z) direction. 
!
      case ('const-xyz')
        if (headtt) print*,'const-xyz'
        if (lpencil(i_uu)) then
          p%uu(:,1)=ampl_kinflow_x*cos(omega_kinflow*t)*exp(eps_kinflow*t)
          p%uu(:,2)=ampl_kinflow_y*cos(omega_kinflow*t)*exp(eps_kinflow*t)
          p%uu(:,3)=ampl_kinflow_z*cos(omega_kinflow*t)*exp(eps_kinflow*t)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  ABC-flow
!
      case ('ABC') 
        if (headtt) print*,'ABC flow'
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=ABC_A*sin(kz_uukin*z(n))    +ABC_C*cos(ky_uukin*y(m))
          p%uu(:,2)=ABC_B*sin(kx_uukin*x(l1:l2))+ABC_A*cos(kz_uukin*z(n))
          p%uu(:,3)=ABC_C*sin(ky_uukin*y(m))    +ABC_B*cos(kx_uukin*x(l1:l2))
        endif
! divu
        if (lpencil(i_divu)) p%divu=0.
!
!  nocosine or Archontis flow
!
      case ('nocos') 
          if (headtt) print*,'nocosine or Archontis flow'
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=ABC_A*sin(kz_uukin*z(n))
          p%uu(:,2)=ABC_B*sin(kx_uukin*x(l1:l2))
          p%uu(:,3)=ABC_C*sin(ky_uukin*y(m))
        endif
! divu
        if (lpencil(i_divu)) p%divu=0.
!
!  Gen-Roberts flow (negative helicity)
!
      case ('roberts') 
        if (headtt) print*,'Glen Roberts flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
! uu
        sqrt2=sqrt(2.)
        if (lpencil(i_uu)) then
          eps1=1.-eps_kinflow
          p%uu(:,1)=+eps1*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,2)=-eps1*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,3)=sqrt2*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
        endif
! divu
        if (lpencil(i_divu)) & 
            p%divu= (kx_uukin-ky_uukin)*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
! uij
        if (lpencil(i_uij)) then
          p%uij(:,1,1)=+eps1*kx_uukin*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uij(:,1,2)=-eps1*ky_uukin*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uij(:,1,3)=+0.
          p%uij(:,2,1)=+eps1*kx_uukin*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uij(:,2,2)=-eps1*ky_uukin*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uij(:,2,3)=+0.
          p%uij(:,3,1)=sqrt2*kx_uukin*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uij(:,3,2)=sqrt2*ky_uukin*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uij(:,3,3)=+0.
        endif
!
! Chandrasekhar-Kendall Flow
!
      case ('ck') 
        if (headtt) print*,'Chandrasekhar-Kendall flow'
! uu
        ell=kinflow_ck_ell
        Balpha=kinflow_ck_Balpha
        ck_r = x(l1:l2)
        ck_rsqr = x(l1:l2)*x(l1:l2)
        if (lpencil(i_uu)) then
          p%uu(:,1)=ampl_kinflow*Pl(m)*(  &
              (ell*(ell+1)/(Balpha*ck_rsqr))*Zl(l1:l2) &
              -(2./(Balpha*ck_r))*dZldr(l1:l2)  )
          p%uu(:,2)=ampl_kinflow*( &
              dZldr(l1:l2)/Balpha- Zl(l1:l2)/(Balpha*ck_r) &
              )*dPldtheta(m)
          p%uu(:,3)=-ampl_kinflow*Zl(l1:l2)*dPldtheta(m)
        endif
! divu
        if (lpencil(i_divu)) p%divu= 0.
!
!  Glen-Roberts flow (positive helicity)
!
      case ('poshel-roberts') 
        if (headtt) print*,'Pos Helicity Roberts flow; eps1=',eps1
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  Glen-Roberts flow (positive helicity), alternative version
!
      case ('varhel-roberts')
        if (headtt) print*,'Pos Helicity Roberts flow; eps1=',eps1
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)*eps1
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  1-D Glen-Roberts flow (positive helicity, no y-dependence)
!
      case ('poshel-roberts-1d') 
        if (headtt) print*,'Pos Helicity Roberts flow; kx_uukin=',kx_uukin
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*eps1
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*sqrt(2.)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  Glen-Roberts flow (x-direction, positive helicity)
!  x -> y
!  y -> z
!  z -> x
!
      case ('xdir-roberts')
        if (headtt) print*,'x-dir Roberts flow; ky_uukin,kz_uukin=',ky_uukin,kz_uukin
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
! uu
        if (lpencil(i_uu)) then
          p%uu(:,2)=-fac*cos(ky_uukin*y(m))*sin(kz_uukin*z(n))*eps1
          p%uu(:,3)=+fac*sin(ky_uukin*y(m))*cos(kz_uukin*z(n))*eps1
          p%uu(:,1)=+fac*cos(ky_uukin*y(m))*cos(kz_uukin*z(n))*sqrt(2.)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  z-dependent Roberts flow (positive helicity)
!
      case ('zdep-roberts') 
        if (headtt) print*,'z-dependent Roberts flow; kx,ky=',kx_uukin,ky_uukin
        fpara=ampl_kinflow*(quintic_step(z(n),-1.+eps_kinflow,eps_kinflow) &
            -quintic_step(z(n),+1.-eps_kinflow,eps_kinflow))
        dfpara=ampl_kinflow*(quintic_der_step(z(n),-1.+eps_kinflow,eps_kinflow)&
            -quintic_der_step(z(n),+1.-eps_kinflow,eps_kinflow))
!
        sqrt2=sqrt(2.)
        sqrt21k1=1./(sqrt2*kx_uukin)
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m)) &
              -dfpara*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt21k1
          p%uu(:,2)=+sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)) &
              -dfpara*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*sqrt21k1
          p%uu(:,3)=+fpara*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt2
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with cosinusoidal helicity variation (version 1)
!
      case ('IncohRoberts1') 
        if (headtt) print*,'Roberts flow with cosinusoidal helicity;',&
            ' kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        eps1=(1.-eps_kinflow)*cos(omega_kinflow*t)
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with cosinusoidal helicity variation (version 2)
!
      case ('IncohRoberts2') 
        if (headtt) print*,'Roberts flow with cosinusoidal helicity;',&
            'kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        eps1=(1.-eps_kinflow)*cos(omega_kinflow*t)
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)*eps1
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with helicity variation and shear (version 1)
!  Must use shear module and set eps_kinflow equal to shear
!
      case ('ShearRoberts1') 
        if (headtt) print*,'Roberts flow with cosinusoidal helicity;',&
            'kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        ky_uukin=1.
        kx_uukin=ky_uukin*(mod(.5-eps_kinflow*t,1.D0)-.5)
        if (ip==11.and.m==4.and.n==4) write(21,*) t,kx_uukin
        eps1=cos(omega_kinflow*t)
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with helicity variation and shear (version 2)
!  Must use shear module and set eps_kinflow equal to shear
!
      case ('ShearRoberts2')
        if (headtt) print*,'Roberts flow with cosinusoidal helicity;',& 
            'kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        ky_uukin=1.
        kx_uukin=ky_uukin*(mod(.5-eps_kinflow*t,1.D0)-.5)
        if (ip==11.and.m==4.and.n==4) write(21,*) t,kx_uukin
        eps1=cos(omega_kinflow*t)
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)*eps1
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  Taylor-Green flow
!
      case ('TG') 
        if (headtt) print*,'Taylor-Green flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=+2.*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*cos(kz_uukin*z(n))
          p%uu(:,2)=-2.*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*cos(kz_uukin*z(n))
          p%uu(:,3)=+0.
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  Galloway-Proctor flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!
      case ('Galloway-Proctor') 
        if (headtt) print*,'Galloway-Proctor flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*sin(ky_uukin*y(m)    +esint)
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)
          p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+ecost)+cos(ky_uukin*y(m)+esint))
        endif
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Galloway-Proctor flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!
      case ('Galloway-Proctor-nohel') 
        if (headtt) print*,'nonhelical Galloway-Proctor flow;',& 
            'kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow*sqrt(1.5)
        fac2=ampl_kinflow*sqrt(6.)
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=+fac*cos(ky_uukin*y(m)    +esint)
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)
          p%uu(:,3)=-fac2*(sin(kx_uukin*x(l1:l2)+ecost)*cos(ky_uukin*y(m)+esint))
        endif
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Otani flow, U=curl(psi*zz) + psi*zz, where
!  psi = 2*cos^2t * cosx - 2*csin2t * cosy
!
      case ('Otani') 
        if (headtt) print*,'Otani flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=2.*ampl_kinflow
        sin2t=sin(omega_kinflow*t)**2
        cos2t=cos(omega_kinflow*t)**2
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=fac*sin2t*sin(ky_uukin*y(m))
          p%uu(:,2)=fac*cos2t*sin(kx_uukin*x(l1:l2))
          p%uu(:,3)=fac*(cos2t*cos(kx_uukin*x(l1:l2))-sin2t*cos(ky_uukin*y(m)))
        endif
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=p%uu
!
!  Tilgner flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!
      case ('Tilgner') 
        if (headtt) print*,'Tilgner flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow*sqrt(2.)
        epst=eps_kinflow*t
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac* sin(ky_uukin*y(m)         )
          p%uu(:,2)=+fac* sin(kx_uukin*(x(l1:l2)+epst))
          p%uu(:,3)=-fac*(cos(kx_uukin*(x(l1:l2)+epst))+cos(ky_uukin*y(m)))
        endif
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Tilgner flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!  Here, W in Tilgner's equation is chosen to be 0.25.
!
      case ('Tilgner-orig') 
        if (headtt) print*,'original Tilgner flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        epst=eps_kinflow*t
        kx_uukin=2.*pi
        ky_uukin=2.*pi
        sqrt2=sqrt(2.)
        WW=0.25
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=+fac*sqrt2*sin(kx_uukin*(x(l1:l2)-epst))*cos(ky_uukin*y(m))
          p%uu(:,2)=-fac*sqrt2*cos(kx_uukin*(x(l1:l2)-epst))*sin(ky_uukin*y(m))
          p%uu(:,3)=+fac*2.*WW*sin(kx_uukin*(x(l1:l2)-epst))*sin(ky_uukin*y(m))
        endif
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Galloway-Proctor flow with random temporal phase
!
      case ('Galloway-Proctor-RandomTemporalPhase') 
        if (headtt) print*,'GP-RandomTemporalPhase; kx,ky=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        if (t>tphase_kinflow) then
          call random_number_wrapper(fran1)
          tphase_kinflow=t+dtphase_kinflow
          phase1=pi*(2*fran1(1)-1.)
          phase2=pi*(2*fran1(2)-1.)
        endif
        ecost=eps_kinflow*cos(omega_kinflow*t+phase1)
        esint=eps_kinflow*sin(omega_kinflow*t+phase2)
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*sin(ky_uukin*y(m)    +esint)
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)
          p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+ecost)+cos(ky_uukin*y(m)+esint))
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  Galloway-Proctor flow with random phase
!
      case ('Galloway-Proctor-RandomPhase') 
        if (headtt) print*,'Galloway-Proctor-RandomPhase; kx,ky=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        if (t>tphase_kinflow) then
          call random_number_wrapper(fran1)
          tphase_kinflow=t+dtphase_kinflow
          phase1=eps_kinflow*pi*(2*fran1(1)-1.)
          phase2=eps_kinflow*pi*(2*fran1(2)-1.)
        endif
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*sin(ky_uukin*y(m)    +phase1)*ky_uukin
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+phase2)*kx_uukin
          p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+phase2)+cos(ky_uukin*y(m)+phase1))
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  original Galloway-Proctor flow
!
      case ('Galloway-Proctor-orig') 
        if (headtt) print*,'Galloway-Proctor-orig flow; kx_uukin=',kx_uukin
        fac=sqrt(1.5)*ampl_kinflow
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=+fac*cos(ky_uukin*y(m)    +esint)*ky_uukin
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)*kx_uukin
          p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+ecost)+sin(ky_uukin*y(m)+esint))
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  Potential flow, u=gradphi, with phi=coskx*X cosky*Y coskz*Z,
!  and X=x-ct, Y=y-ct, Z=z-ct.
!
      case ('potential') 
        if (headtt) print*,'potential; ampl_kinflow,omega_kinflow=',& 
            ampl_kinflow,omega_kinflow
        if (headtt) print*,'potential; ki_uukin=',kx_uukin,ky_uukin,kz_uukin
        fac=ampl_kinflow
        cxt=cx_uukin*t
        cyt=cy_uukin*t
        czt=cz_uukin*t
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*kx_uukin*&
              sin(kx_uukin*(x(l1:l2)-cxt))*cos(ky_uukin*(y(m)-cyt))*&
              cos(kz_uukin*(z(n)-czt-phasez_uukin))
          p%uu(:,2)=-fac*ky_uukin*&
              cos(kx_uukin*(x(l1:l2)-cxt))*sin(ky_uukin*(y(m)-cyt))*&
              cos(kz_uukin*(z(n)-czt-phasez_uukin))
          p%uu(:,3)=-fac*kz_uukin*&
              cos(kx_uukin*(x(l1:l2)-cxt))*cos(ky_uukin*(y(m)-cyt))*&
              sin(kz_uukin*(z(n)-czt-phasez_uukin))
        endif
        if (lpencil(i_divu)) p%divu=-fac*(kx_uukin**2+ky_uukin**2+kz_uukin**2) &
            *cos(kx_uukin*x(l1:l2)-cxt)*cos(ky_uukin*y(m)-cyt)*cos(kz_uukin*z(n)-czt)
!
!  2nd Potential flow, u=gradphi, with phi=cos(kx*X+ky*Y+kz*Z),
!  and X=x-ct, Y=y-ct, Z=z-ct.
!
      case ('potential2') 
        if (headtt) print*,'2nd potential; ampl_kinflow,omega_kinflow=',&
            ampl_kinflow,omega_kinflow
        if (headtt) print*,'2nd potential; ki_uukin=',kx_uukin,ky_uukin,kz_uukin
        fac=ampl_kinflow
        cxt=cx_uukin*t
        cyt=cy_uukin*t
        czt=cz_uukin*t
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*kx_uukin*&
              sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)+phasez_uukin)
          p%uu(:,2)=-fac*ky_uukin*&
              sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)+phasez_uukin)
          p%uu(:,3)=-fac*kz_uukin*&
            sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)+phasez_uukin)
        endif
        if (lpencil(i_divu)) p%divu=-fac*(kx_uukin**2+ky_uukin**2+kz_uukin**2) &
            *cos(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)+phasez_uukin)
!
!  1-D Potential flow, u=gradphi, with phi=cos(kx*X+ky*Y+kz*Z),
!  and X=x-ct, Y=y-ct, Z=z-ct.
!
      case ('potentialz') 
        if (headtt) print*,'1-D potential; ampl_kinflow,omega_kinflow=',&
            ampl_kinflow,omega_kinflow
        if (headtt) print*,'1-D potential; ki_uukin=',kx_uukin,ky_uukin,kz_uukin
        fac=ampl_kinflow
        omt=omega_kinflow*t
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*kx_uukin*&
              sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)-omt)
          p%uu(:,2)=-fac*ky_uukin*&
              sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)-omt)
          p%uu(:,3)=-fac*kz_uukin*&
              sin(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)-omt)
        endif
        if (lpencil(i_divu)) p%divu=-fac*(kx_uukin**2+ky_uukin**2+kz_uukin**2) &
            *cos(kx_uukin*x(l1:l2)+ky_uukin*y(m)+kz_uukin*z(n)-omt)
!
!  Potential random flow, u=gradphi, with phi=cos(x-x0)*cosy*cosz;
!  assume kx_uukin=ky_uukin=kz_uukin.
!
      case ('potential_random') 
        if (headtt) print*,'potential_random; kx_uukin,ampl_kinflow=',ampl_kinflow
        if (headtt) print*,'potential_random; kx_uukin=',kx_uukin,ky_uukin,kz_uukin
        fac=ampl_kinflow
        argx=kx_uukin*(x(l1:l2)-location(1))
        argy=ky_uukin*(y(m)-location(2))
        argz=kz_uukin*(z(n)-location(3))
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*kx_uukin*sin(argx)*cos(argy)*cos(argz)
          p%uu(:,2)=-fac*ky_uukin*cos(argx)*sin(argy)*cos(argz)
          p%uu(:,3)=-fac*kz_uukin*cos(argx)*cos(argy)*sin(argz)
        endif
        if (lpencil(i_divu)) p%divu=fac
!
!  Convection rolls
!  Stream function: psi_y = cos(kx*x) * cos(kz*z)
!
      case ('rolls')
        if (headtt) print*,'Convection rolls; kx_kinflow,kz_uukin=',kx_kinflow,kz_kinflow
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=+ampl_kinflow*kz_kinflow*cos(kx_kinflow*x(l1:l2))*sin(kz_kinflow*z(n))
          p%uu(:,2)=+0.
          p%uu(:,3)=-ampl_kinflow*kx_kinflow*sin(kx_kinflow*x(l1:l2))*cos(kz_kinflow*z(n))
        endif
! divu
        if (lpencil(i_divu)) p%divu=0.
!
!  Convection rolls
!  Stream function: psi_y = sin(kx*x) * sin(kz*z)
!
      case ('rolls2')
        if (headtt) print*,'Convection rolls2; kx_kinflow,kz_uukin=',kx_kinflow,kz_kinflow
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-ampl_kinflow*kz_kinflow*sin(kx_kinflow*x(l1:l2))*cos(kz_kinflow*z(n))
          p%uu(:,2)=+0.
          p%uu(:,3)=+ampl_kinflow*kx_kinflow*cos(kx_kinflow*x(l1:l2))*sin(kz_kinflow*z(n))
        endif
! divu
        if (lpencil(i_divu)) p%divu=0.
!
!  Twist (Yousef & Brandenburg 2003)
!
      case ('Twist') 
        if (headtt) print*,'Twist flow; eps_kinflow,kx=',eps_kinflow,kx_uukin
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=eps_kinflow*z(n)*sin(kx_uukin*x(l1:l2))
          p%uu(:,3)=1.+cos(kx_uukin*x(l1:l2))
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  Eddy (Brandenburg & Zweibel 1994)
!
      case ('eddy') 
        if (headtt) print*,'eddy flow; eps_kinflow,kx=',eps_kinflow,kx_uukin
        cos1_mn=max(cos(.5*pi*p%rcyl_mn),0.)
        cos2_mn=cos1_mn**2
        tmp_mn=-.5*pi*p%rcyl_mn1*sin(.5*pi*p%rcyl_mn)*ampl_kinflow* &
            4.*cos2_mn*cos1_mn
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=+tmp_mn*y(m)
          p%uu(:,2)=-tmp_mn*x(l1:l2)
          p%uu(:,3)=eps_kinflow*ampl_kinflow*cos2_mn**2
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  Shearing wave
!
      case ('ShearingWave') 
        if (headtt) print*,'ShearingWave flow; Sshear,eps_kinflow=',Sshear,eps_kinflow
        ky_uukin=1.
        kx_uukin=-ky_uukin*Sshear*t
        k21=1./(kx_uukin**2+ky_uukin**2)
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-ky_uukin*k21*cos(kx_uukin*x(l1:l2)+ky_uukin*y(m))
          p%uu(:,2)=+kx_uukin*k21*cos(kx_uukin*x(l1:l2)+ky_uukin*y(m))
          p%uu(:,3)=eps_kinflow*cos(kx_uukin*x(l1:l2)+ky_uukin*y(m))
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  Helical Shearing wave
!
      case ('HelicalShearingWave') 
        if (headtt) print*,'ShearingWave flow; Sshear,eps_kinflow=',Sshear,eps_kinflow
        ky_uukin=1.
        kx_uukin=-ky_uukin*Sshear*t
        k21=1./(kx_uukin**2+ky_uukin**2)
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
          p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
          p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
! step function along z
!
      case ('zstep') 
        if (headtt) print*,'wind:step function along z'
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=wind_amp*step_scalar(z(n),wind_rmin,wind_step_width)
        endif
!
! uniform radial shear with a cutoff at rmax
! uniform radial shear with saturation.
!
      case ('rshear-sat') 
        if (headtt) print*,'radial shear saturated'
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
!          omega0=uphi_at_rmax/uphi_rmax
!          shear = arctanh(omega0)/(uphi_rmax-x(l1))
          p%uu(:,3)=ampl_kinflow*x(l1:l2)*sinth(m)*tanh(10.*(x(l1:l2)-x(l1)))
!*tanh(10.*(x(l1:l2)-x(l1)))
!          write(*,*)'DM',m,n,p%uu(:,3)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  U_phi = r*sin(theta)*Omega(r) with Omega(r)=1-r
!
      case ('(1-x)*x*siny')
        if (headtt) print*,'Omega=Omega0*(1-r), Omega0=',ampl_kinflow
        local_Omega=ampl_kinflow*(1.-x(l1:l2))
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=local_Omega*x(l1:l2)*sinth(m)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  U_phi adapted from
!  Ghizaru-Charbonneau-Smolarkiewicz (ApJL 715:L133-L137, 2010)
!
      case ('gcs') 
        if (headtt) print*,'gcs:gcs_rzero ',gcs_rzero
        fac=ampl_kinflow
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          local_Omega=fac*exp(-((x(l1:l2)-xyz1(1))/gcs_rzero)**2- &
              ((pi/2-y(m))/gcs_psizero**2))
          p%uu(:,3)= local_Omega*x(l1:l2)*sinth(m)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  Vertical wind
!
      case ('vertical-wind') 
        if (headtt) print*,'vertical-wind along z'
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=wind_amp*z(n)
        endif
        if (lpencil(i_divu)) p%divu=wind_amp
!
!  Vertical wind that goes to zero for z < 0.
!
      case ('vertical-wind-surface') 
        if (headtt) print*,'vertical-wind along z'
! uu
        if (lpencil(i_uu)) then
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=wind_amp*z(n)*0.5*(1.+erfunc((z(n)/wind_step_width)))
        endif
        if (lpencil(i_divu)) p%divu=wind_amp
!
! Radial wind
!
      case ('radial-wind') 
        if (headtt) print*,'Radial wind'
        select case (wind_profile)
        case ('none'); wind_prof=0.;div_uprof=0.
        case ('constant'); wind_prof=1.;div_uprof=0.
        case ('radial-step')
          wind_prof=step(x(l1:l2),wind_rmin,wind_step_width)
          div_uprof=der_step(x(l1:l2),wind_rmin,wind_step_width)
!          der6_uprof=der6_step(x(l1:l2),wind_rmin,wind_step_width)
          der6_uprof=0.
        case default;
          call fatal_error('hydro_kinematic:kinflow = radial wind', &
              'no such wind profile. ')
        endselect
!
        if (lpencil(i_uu)) then
          p%uu(:,1)=wind_amp*wind_prof
          p%uu(:,2)=0.
          p%uu(:,3)=0.
        endif
        p%divu=wind_amp*div_uprof
        p%der6u(:,1)=wind_amp*der6_uprof
        p%der6u(:,2)= 0.
        p%der6u(:,3)= 0.
!
!  meridional circulation; psi=.5*(x-x0)*(x-x1)*(y-y0)*(y-y1), so
!  ux=+dpsi/dy=+(x-x0)*(x-x1)*y
!  uy=-dpsi/dx=-x*(y-y0)*(y-y1)
!
      case ('circ_cartesian') 
        if (headtt) print*,'just circulation'
        if (lpencil(i_uu)) then
          p%uu(:,1)=+(x(l1:l2)-xyz0(1))*(x(l1:l2)-xyz1(1))*y(m)
          p%uu(:,2)=-x(l1:l2)*(y(m)-xyz0(2))*(y(m)-xyz1(2))
          p%uu(:,3)=0.
        endif
        p%divu=0.
        p%der6u(:,1)=wind_amp*der6_uprof
        p%der6u(:,2)= 0.
        p%der6u(:,3)= 0.
!
!  meridional circulation
!
      case ('circ_spherical') 
        if (headtt) print*,'just circulation (spherical)'
        if (lpencil(i_uu)) then
          p%uu(:,1)=+(x(l1:l2)-xyz0(1))*(x(l1:l2)-xyz1(1))*y(m)
          p%uu(:,2)=-x(l1:l2)*(y(m)-xyz0(2))*(y(m)-xyz1(2))
          p%uu(:,3)=0.
        endif
        p%divu=0.
        p%der6u(:,1)=wind_amp*der6_uprof
        p%der6u(:,2)= 0.
        p%der6u(:,3)= 0.
!
! Radial wind+circulation
!
      case ('radial-wind+circ') 
        if (headtt) print*,'Radial wind and circulation'
        select case (wind_profile)
        case ('none'); wind_prof=0.;div_uprof=0.
        case ('constant'); wind_prof=1.;div_uprof=0.
        case ('radial-step')
          wind_prof=step(x(l1:l2),wind_rmin,wind_step_width)
          div_uprof=der_step(x(l1:l2),wind_rmin,wind_step_width)
          der6_uprof=0.
        case default;
          call fatal_error('hydro_kinematic: kinflow= radial_wind-circ',&
              'no such wind profile. ')
        endselect
        rone=xyz0(1)
        theta=y(m)
        theta1=xyz0(2)
        if (lpencil(i_uu)) then
          vel_prof=circ_amp*(1+stepdown(x(l1:l2),circ_rmax,circ_step_width))
          div_vel_prof=-der_step(x(l1:l2),circ_rmax,circ_step_width)
          p%uu(:,1)=vel_prof*(r1_mn**2)*(sin1th(m))*(&
              2*sin(theta-theta1)*cos(theta-theta1)*cos(theta)&
              -sin(theta)*sin(theta-theta1)**2)*&
              (x(l1:l2)-1.)*(x(l1:l2)-rone)**2 + &
              wind_amp*wind_prof
          p%uu(:,2)=-vel_prof*r1_mn*sin1th(m)*(&
              cos(theta)*sin(theta-theta1)**2)*&
              (x(l1:l2)-rone)*(3*x(l1:l2)-rone-2.)
          p%uu(:,3)=0.
        endif
        p%divu=div_uprof*wind_amp + &
            div_vel_prof*(r1_mn**2)*(sin1th(m))*(&
            2*sin(theta-theta1)*cos(theta-theta1)*cos(theta)&
            -sin(theta)*sin(theta-theta1)**2)*&
            (x(l1:l2)-1.)*(x(l1:l2)-rone)**2
        p%der6u(:,1)=wind_amp*der6_uprof
        p%der6u(:,2)= 0.
        p%der6u(:,3)= 0.
!
!  Radial shear + radial wind + circulation.
!
      case ('rshear-sat+rwind+circ') 
        if (headtt) print*,'radial shear, wind, circulation: not complete yet'
!
!  First set wind and cirulation.
!
        select case (wind_profile)
        case ('none'); wind_prof=0.;div_uprof=0.
        case ('constant'); wind_prof=1.;div_uprof=0.
!
        case ('radial-step')
          wind_prof=step(x(l1:l2),wind_rmin,wind_step_width)
          div_uprof=der_step(x(l1:l2),wind_rmin,wind_step_width)
          der6_uprof=0.
        case default;
          call fatal_error('hydro_kinematic: kinflow= radial-shear+rwind+rcirc',&
              'no such wind profile. ')
        endselect
        rone=xyz0(1)
        theta=y(m)
        theta1=xyz0(2)
        if (lpencil(i_uu)) then
          vel_prof=circ_amp*(1+stepdown(x(l1:l2),circ_rmax,circ_step_width))
          div_vel_prof=-der_step(x(l1:l2),circ_rmax,circ_step_width)
          p%uu(:,1)=vel_prof*(r1_mn**2)*(sin1th(m))*(&
              2*sin(theta-theta1)*cos(theta-theta1)*cos(theta)&
              -sin(theta)*sin(theta-theta1)**2)*&
              (x(l1:l2)-1.)*(x(l1:l2)-rone)**2 + &
              wind_amp*wind_prof
          p%uu(:,2)=-vel_prof*r1_mn*sin1th(m)*(&
              cos(theta)*sin(theta-theta1)**2)*&
              (x(l1:l2)-rone)*(3*x(l1:l2)-rone-2.)
          p%uu(:,3)=x(l1:l2)*sinth(m)*tanh(10.*(x(l1:l2)-x(l1)))
        endif
        p%divu=div_uprof*wind_amp + &
            div_vel_prof*(r1_mn**2)*(sin1th(m))*(&
            2*sin(theta-theta1)*cos(theta-theta1)*cos(theta)&
            -sin(theta)*sin(theta-theta1)**2)*&
            (x(l1:l2)-1.)*(x(l1:l2)-rone)**2
        p%der6u(:,1)=wind_amp*der6_uprof
        p%der6u(:,2)= 0.
        p%der6u(:,3)= 0.
!
!  KS-flow
!
      case ('KS') 
        p%uu=0.
        do modeN=1,KS_modes  ! sum over KS_modes modes
          kdotxwt=KS_k(1,modeN)*x(l1:l2)+(KS_k(2,modeN)*y(m)+KS_k(3,modeN)*z(n))+KS_omega(modeN)*t
          cos_kdotxwt=cos(kdotxwt) ;  sin_kdotxwt=sin(kdotxwt)
          if (lpencil(i_uu)) then
            p%uu(:,1) = p%uu(:,1) + cos_kdotxwt*KS_A(1,modeN) + sin_kdotxwt*KS_B(1,modeN)
            p%uu(:,2) = p%uu(:,2) + cos_kdotxwt*KS_A(2,modeN) + sin_kdotxwt*KS_B(2,modeN)
            p%uu(:,3) = p%uu(:,3) + cos_kdotxwt*KS_A(3,modeN) + sin_kdotxwt*KS_B(3,modeN)
          endif
        enddo
        if (lpencil(i_divu))  p%divu = 0.
!
! no kinematic flow.
!
      case ('none')
        if (headtt) print*,'kinflow = none'
        if (lpencil(i_uu)) p%uu=0.
! divu
        if (lpencil(i_divu)) p%divu=0.
      case default;
        call fatal_error('hydro_kinematic:', 'kinflow not found')
      end select
!
! kinflows end here
!
! u2
!
      if (lpencil(i_u2)) call dot2_mn(p%uu,p%u2)
      if (idiag_ekin/=0 .or. idiag_ekintot/=0) then
        lpenc_diagnos(i_rho)=.true.
        lpenc_diagnos(i_u2)=.true.
      endif
! ou
      if (lpencil(i_ou)) call dot_mn(p%oo,p%uu,p%ou)
!
!  Calculate maxima and rms values for diagnostic purposes
!
      if (ldiagnos) then
        if (idiag_urms/=0)  call sum_mn_name(p%u2,idiag_urms,lsqrt=.true.)
        if (idiag_umax/=0)  call max_mn_name(p%u2,idiag_umax,lsqrt=.true.)
        if (idiag_uzrms/=0) &
            call sum_mn_name(p%uu(:,3)**2,idiag_uzrms,lsqrt=.true.)
        if (idiag_uzmax/=0) &
            call max_mn_name(p%uu(:,3)**2,idiag_uzmax,lsqrt=.true.)
        if (idiag_u2m/=0)   call sum_mn_name(p%u2,idiag_u2m)
        if (idiag_um2/=0)   call max_mn_name(p%u2,idiag_um2)
!
        if (idiag_ekin/=0)  call sum_mn_name(.5*p%rho*p%u2,idiag_ekin)
        if (idiag_ekintot/=0) &
            call integrate_mn_name(.5*p%rho*p%u2,idiag_ekintot)
      endif
      if (idiag_divum/=0)  call sum_mn_name(p%divu,idiag_divum)
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_pencils_hydro
!***********************************************************************
    subroutine hydro_before_boundary(f)
!
!  Dummy routine
!
!   16-dec-10/bing: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine hydro_before_boundary
!***********************************************************************
    subroutine duu_dt(f,df,p)
!
!  velocity evolution, dummy routine
!
!   7-jun-02/axel: adapted from hydro
!
      use Diagnostics
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension (nx,3) :: uu
!
      intent(in)  :: df,p
      intent(out) :: f
!
!  uu/dx for timestep (if kinflow is set)
!
      if (kinflow/='') then
        if (lfirst.and.ldt) then
          uu=p%uu
           if (lspherical_coords) then
            advec_uu=abs(uu(:,1))*dx_1(l1:l2)+ &
                abs(uu(:,2))*dy_1(  m  )*r1_mn+ &
                abs(uu(:,3))*dz_1(  n  )*r1_mn*sin1th(m)
          elseif (lcylindrical_coords) then
            advec_uu=abs(uu(:,1))*dx_1(l1:l2)+ &
                abs(uu(:,2))*dy_1(  m  )*rcyl_mn1+ &
                abs(uu(:,3))*dz_1(  n  )
          else
            advec_uu=abs(uu(:,1))*dx_1(l1:l2)+ &
                abs(uu(:,2))*dy_1(  m  )+ &
                abs(uu(:,3))*dz_1(  n  )
          endif
        endif
      endif
      if (headtt.or.ldebug) print*, 'duu_dt: max(advec_uu) =', maxval(advec_uu)
!
!  Store uu in auxiliary variable if requested.
!  Just neccessary immediately before writing snapshots, but how would we
!  know we are?
!
     if (lkinflow_as_aux) f(l1:l2,m,n,iux:iuz)=p%uu
!
!  Calculate maxima and rms values for diagnostic purposes.
!
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'duu_dt: diagnostics ...'
        if (idiag_oum/=0) call sum_mn_name(p%ou,idiag_oum)
!
!  Kinetic field components at one point (=pt).
!
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_uxpt/=0) call save_name(p%uu(lpoint-nghost,1),idiag_uxpt)
          if (idiag_uypt/=0) call save_name(p%uu(lpoint-nghost,2),idiag_uypt)
          if (idiag_uzpt/=0) call save_name(p%uu(lpoint-nghost,3),idiag_uzpt)
          if (idiag_phase1/=0) call save_name(phase1,idiag_phase1)
          if (idiag_phase2/=0) call save_name(phase2,idiag_phase2)
        endif
      endif
!
      call keep_compiler_quiet(df)
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine time_integrals_hydro(f,p)
!
!   1-jul-08/axel: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f,p
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
!
    endsubroutine time_integrals_hydro
!***********************************************************************
    subroutine traceless_strain(uij,divu,sij,uu)
!
!  Calculates traceless rate-of-strain tensor sij from derivative tensor uij
!  and divergence divu within each pencil;
!  curvilinear co-ordinates require optional velocity argument uu
!
!  16-oct-09/MR: dummy
!
    real, dimension(nx,3,3)         :: uij, sij
    real, dimension(nx)             :: divu
    real, dimension(nx,3), optional :: uu
!
    intent(in) :: uij, divu, sij
!
    call keep_compiler_quiet(uij,sij)
    call keep_compiler_quiet(divu)
    call keep_compiler_quiet(present(uu))
!
    endsubroutine traceless_strain
!***********************************************************************
   subroutine coriolis_cartesian(df,uu,velind)
!
!  Coriolis terms for cartesian geometry.
!
!  30-oct-09/MR: outsourced, parameter velind added
!  checked to be an equivalent change by auot-test conv-slab-noequi, mdwarf
!
      real, dimension (mx,my,mz,mvar), intent(out) :: df
      real, dimension (nx,3),          intent(in)  :: uu
      integer,                         intent(in)  :: velind
!
      call keep_compiler_quiet(df)
      call keep_compiler_quiet(uu)
      call keep_compiler_quiet(velind)
!
   endsubroutine coriolis_cartesian
!***********************************************************************
    subroutine calc_lhydro_pars(f)
!
!  Dummy routine.
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lhydro_pars
!***********************************************************************
    subroutine random_isotropic_KS_setup(initpower,kmin,kmax)
!
!  Produces random, isotropic field from energy spectrum following the
!  KS method (Malik and Vassilicos, 1999.)
!
!  More to do; unsatisfactory so far - at least for a steep power-law
!  energy spectrum.
!
!  27-may-05/tony: modified from snod's KS hydro initial
!  03-feb-06/weezy: Attempted rewrite to guarantee periodicity of
!                    KS modes.
!
      use Sub, only: cross, dot2
      use General, only: random_number_wrapper
!
      integer :: modeN
!
      real, dimension (3) :: k_unit
      real, dimension (3) :: ee,e1,e2
      real, dimension (6) :: r
      real :: initpower,kmin,kmax
      real, dimension(KS_modes) :: k,dk,energy,ps
      real :: theta,phi,alpha,beta
      real :: ex,ey,ez,norm,a
!
      allocate(KS_k(3,KS_modes))
      allocate(KS_A(3,KS_modes))
      allocate(KS_B(3,KS_modes))
      allocate(KS_omega(KS_modes))
!
      kmin=2.*pi      !/(1.0*Lxyz(1))
      kmax=128.*pi    !nx*pi
      a=(kmax/kmin)**(1./(KS_modes-1.))
!
!  Loop over all modes.
!
      do modeN=1,KS_modes
!
!  Pick wavenumber.
!
        k=kmin*(a**(modeN-1.))
!
!  Pick 4 random angles for each mode.
!
        call random_number_wrapper(r);
        theta=pi*(r(1) - 0.)
        phi=pi*(2*r(2) - 0.)
        alpha=pi*(2*r(3) - 0.)
        beta=pi*(2*r(4) - 0.)
!
!  Make a random unit vector by rotating fixed vector to random position
!  (alternatively make a random transformation matrix for each k).
!
        k_unit(1)=sin(theta)*cos(phi)
        k_unit(2)=sin(theta)*sin(phi)
        k_unit(3)=cos(theta)
!
        energy=(((k/kmin)**2. +1.)**(-11./6.))*(k**2.) &
            *exp(-0.5*(k/kmax)**2.)
!
!  Make a vector KS_k of length k from the unit vector for each mode.
!
        KS_k(:,modeN)=k*k_unit(:)
        KS_omega(:)=sqrt(energy(:)*(k(:)**3.))
!
!  Construct basis for plane having rr normal to it
!  (bit of code from forcing to construct x', y').
!
      if ((k_unit(2)==0).and.(k_unit(3)==0)) then
          ex=0.; ey=1.; ez=0.
        else
          ex=1.; ey=0.; ez=0.
        endif
        ee = (/ex, ey, ez/)
!
        call cross(k_unit(:),ee,e1)
!  e1: unit vector perp. to KS_k
        call dot2(e1,norm); e1=e1/sqrt(norm)
        call cross(k_unit(:),e1,e2)
!  e2: unit vector perp. to KS_k, e1
        call dot2(e2,norm); e2=e2/sqrt(norm)
!
!  Make two random unit vectors KS_B and KS_A in the constructed plane.
!
        KS_A(:,modeN) = cos(alpha)*e1 + sin(alpha)*e2
        KS_B(:,modeN) = cos(beta)*e1  + sin(beta)*e2
!
!  Make sure dk is set.
!
        call error('random_isotropic_KS_setup', 'Using uninitialized dk')
        dk=0.                     ! to make compiler happy
!
        ps=sqrt(2.*energy*dk)   !/3.0)
!
!  Give KS_A and KS_B length ps.
!
        KS_A(:,modeN)=ps*KS_A(:,modeN)
        KS_B(:,modeN)=ps*KS_B(:,modeN)
!
      enddo
!
!  Form RA = RA x k_unit and RB = RB x k_unit.
!  Note: cannot reuse same vector for input and output.
!
      do modeN=1,KS_modes
        call cross(KS_A(:,modeN),k_unit(:),KS_A(:,modeN))
        call cross(KS_B(:,modeN),k_unit(:),KS_B(:,modeN))
      enddo
!
      call keep_compiler_quiet(initpower)
!
    endsubroutine random_isotropic_KS_setup
!***********************************************************************
    subroutine random_isotropic_KS_setup_test
!
!  Produces random, isotropic field from energy spectrum following the
!  KS method (Malik and Vassilicos, 1999.)
!  This test case only uses 3 very specific modes (useful for comparison
!  with Louise's kinematic dynamo code.
!
!  03-feb-06/weezy: modified from random_isotropic_KS_setup
!
      use Sub, only: cross
      use General, only: random_number_wrapper
!
      integer :: modeN
!
      real, dimension (3,KS_modes) :: k_unit
      real, dimension(KS_modes) :: k,dk,energy,ps
      real :: initpower,kmin,kmax
!
      allocate(KS_k(3,KS_modes))
      allocate(KS_A(3,KS_modes))
      allocate(KS_B(3,KS_modes))
      allocate(KS_omega(KS_modes))
!
      initpower=-5./3.
      kmin=10.88279619
      kmax=23.50952672
!
      KS_k(1,1)=2.00*pi
      KS_k(2,1)=-2.00*pi
      KS_k(3,1)=2.00*pi
!
      KS_k(1,2)=-4.00*pi
      KS_k(2,2)=0.00*pi
      KS_k(3,2)=2.00*pi
!
      KS_k(1,3)=4.00*pi
      KS_k(2,3)=2.00*pi
      KS_k(3,3)=-6.00*pi
!
      KS_k(1,1)=+1; KS_k(2,1)=-1; KS_k(3,1)=1
      KS_k(1,2)=+0; KS_k(2,2)=-2; KS_k(3,2)=1
      KS_k(1,3)=+0; KS_k(2,3)=-0; KS_k(3,3)=1
!
      k(1)=kmin
      k(2)=14.04962946
      k(3)=kmax
!
      do modeN=1,KS_modes
        k_unit(:,modeN)=KS_k(:,modeN)/k(modeN)
      enddo
!
      kmax=k(KS_modes)
      kmin=k(1)
!
      do modeN=1,KS_modes
        if (modeN==1) dk(modeN)=(k(modeN+1)-k(modeN))/2.
        if (modeN>1.and.modeN<KS_modes) &
            dk(modeN)=(k(modeN+1)-k(modeN-1))/2.
        if (modeN==KS_modes) dk(modeN)=(k(modeN)-k(modeN-1))/2.
      enddo
!
      do modeN=1,KS_modes
         energy(modeN)=((k(modeN)**2 +1.)**(-11./6.))*(k(modeN)**2) &
             *exp(-0.5*(k(modeN)/kmax)**2)
      enddo
!
      ps=sqrt(2.*energy*dk)
!
      KS_A(1,1)=1.00/sqrt(2.00)
      KS_A(2,1)=-1.00/sqrt(2.00)
      KS_A(3,1)=0.00
!
      KS_A(1,2)=1.00/sqrt(3.00)
      KS_A(2,2)=1.00/sqrt(3.00)
      KS_A(3,2)=-1.00/sqrt(3.00)
!
      KS_A(1,3)=-1.00/2.00
      KS_A(2,3)=-1.00/2.00
      KS_A(3,3)=1.00/sqrt(2.00)
!
      KS_B(1,3)=1.00/sqrt(2.00)
      KS_B(2,3)=-1.00/sqrt(2.00)
      KS_B(3,3)=0.00
!
      KS_B(1,1)=1.00/sqrt(3.00)
      KS_B(2,1)=1.00/sqrt(3.00)
      KS_B(3,1)=-1.00/sqrt(3.00)
!
      KS_B(1,2)=-1.00/2.00
      KS_B(2,2)=-1.00/2.00
      KS_B(3,2)=1.00/sqrt(2.00)
!
      do modeN=1,KS_modes
        KS_A(:,modeN)=ps(modeN)*KS_A(:,modeN)
        KS_B(:,modeN)=ps(modeN)*KS_B(:,modeN)
      enddo
!
!  Form RA = RA x k_unit and RB = RB x k_unit.
!
       do modeN=1,KS_modes
         call cross(KS_A(:,modeN),k_unit(:,modeN),KS_A(:,modeN))
         call cross(KS_B(:,modeN),k_unit(:,modeN),KS_B(:,modeN))
       enddo
!
    endsubroutine random_isotropic_KS_setup_test
!***********************************************************************
    subroutine periodic_KS_setup(initpower)
!
!  Produces random, isotropic field from energy spectrum following the
!  KS method, however this setup produces periodic velocity field
!  (assuming box at least (-pi,pi)).
!
!  22-jun-2010/abag coded
!
      use Sub
      use General
!
      real, intent(in) :: initpower
      real, allocatable, dimension(:,:) :: unit_k,k,A,B,orderK
      real, allocatable, dimension(:) :: kk,delk,energy,omega,klengths
      real, allocatable, dimension(:) :: ampA(:), ampB(:)
      real, dimension(3) :: angle,dir_in
      real :: k_option(3,10000),mkunit(10000)
      real :: arg, unity
      real :: bubble, max_box
      real :: j(3),l(3),newa(3),newa2(3)
      integer ::i,s1,num,direction(3)
      logical :: ne
!
!  For the calculation of the velocity field we require
!  KS_k, KS_A/KS_B, KS_omega.
!
      allocate(KS_k(3,KS_modes), KS_A(3,KS_modes),KS_B(3,KS_modes),KS_omega(KS_modes))
!
!  The rest are dummy arrays used to get above.
!
      allocate(unit_k(3,KS_modes), k(3,KS_modes), A(3,KS_modes), B(3,KS_modes))
      allocate(orderk(3,KS_modes), kk(KS_modes), delk(KS_modes), energy(KS_modes))
      allocate(omega(KS_modes), klengths(KS_modes), ampA(KS_modes), ampB(KS_modes))
!
!  Dummy variable.
!
      num=1
!
!  Space wavenumbers out.
!
      bubble=1.
!
!  Needs adaptingf for 2D runs.
!
      max_box=min(nxgrid,nygrid,nzgrid)
!
      if (mod(max_box,4.)/=0) print*, 'warning will not be periodic'
!
      print*, 'calculating KS wavenumbers'
!
      do i=1,10000
        call random_number_wrapper(angle)
        if ((angle(1)-0.0 < epsilon(0.0)) .or. &
            (angle(2)-0.0 < epsilon(0.0)) .or. &
            (angle(3)-0.0 < epsilon(0.0))) then
          call random_number_wrapper(angle)
        endif
!
!  Need 4 meshpoints to resolve a wave.
!
        angle=floor((max_box/4)*angle)
        call random_number_wrapper(dir_in)
        direction=nint(dir_in)
        direction=2*direction -1  !positive or negative directions
!
!  A possible orientation provided we havent's already got this length.
!
        k_option(1,i)=direction(1)*2.*pi*angle(1)
        k_option(2,i)=direction(2)*2.*pi*angle(2)
        k_option(3,i)=direction(3)*2.*pi*angle(3)
!
        if (i==1)then
          k_option(1,i)=2.*pi
          k_option(2,i)=0.
          k_option(3,i)=0.
        endif
!
!  Find the length of the current k_option vector.
!
        mkunit(i)=sqrt((k_option(1,i)**2)+(k_option(2,i)**2)+(k_option(3,i)**2))
!
        if (i==1.and.mkunit(i)>0.)then
          k(:,num)=k_option(:,i)
          klengths(num)=mkunit(i)
        endif
!
!  Now we check that the current length is unique (hasn't come before).
!
        if (i>1.and.num<KS_modes)then
          do s1=i-1,1,-1
            if (mkunit(i)>0.0.and. &
                mkunit(i)<=(mkunit(s1)-bubble).or. &
                mkunit(i)>=(mkunit(s1)+bubble)) then
              ne=.true.
            else
              ne=.false.
              exit
            endif
!
!  If length of current k_option is new......
!
            if (s1==1.and.ne) then
              num=num+1
!
!  Load current k_option into k that we keep.
!
              k(:,num)=k_option(:,i)
!
!  Store the length also.
!
              klengths(num)=mkunit(i)
            endif
          enddo
        endif
        if (i==10000.and.num<KS_modes) print*,"Haven't got",KS_modes,"modes!!!!"
      enddo
!
!  The 1 means ascending order.
!
      call KS_order(klengths,KS_modes,1,kk)
!
      do i=1,KS_modes
        do s1=1,KS_modes
          if (kk(i)==klengths(s1))then
            orderK(:,i)=k(:,s1)
          endif
        enddo
      enddo
!
      k=orderK
      do i=1,KS_modes
        unit_k(:,i)=k(:,i)/kk(i)
      enddo
!
      do i=1,KS_modes
!
!  Now we find delk as defined in Malik & Vassilicos' paper.
!
        if (i==1) delk(i)=(kk(i+1)-kk(i))/2.0
        if (i==KS_modes) delk(i)=(kk(i)-kk(i-1))/2.0
        if (i>1.and.i<KS_modes) delk(i)=(kk(i+1)-kk(i-1))/2.0
      enddo
!
!  Now find A&B that are perpendicular to each of our N wave-vectors.
!
      do i=1,KS_modes
!
!  Define "energy" - here we want k^{initpower} in the inertial range.
!
        energy(i)=1.0+(kk(i))**2
        energy(i)=(kk(i)**2)*(energy(i)**((initpower-2.)/2.))
        energy(i)=energy(i)*exp(-0.5*(kk(i)/kk(KS_modes))**2)
!
!  Set the lengths of A& B as defined in Malik & Vassilicos.
!
        ampA(i)=sqrt(2.0*energy(i)*delk(i))
        ampB(i)=ampA(i)
!
        call random_number(newa)
        call random_number(newa2)
        newa=2.0*newa -1.0
        newa2=2.0*newa2 -1.0
        j=newa
        l=newa2
        j=j/(sqrt(sum(newa**2)))
        l=l/(sqrt(sum(newa2**2)))
!
!  Now take the vector product of k with j (->A) and with l (->B).
!
        A(1,i)=(j(2)*k(3,i))-(j(3)*k(2,i))
        A(2,i)=(j(3)*k(1,i))-(j(1)*k(3,i))
        A(3,i)=(j(1)*k(2,i))-(j(2)*k(1,i))
        unity=sqrt((A(1,i)**2)+(A(2,i)**2)+(A(3,i)**2))
        A(:,i)=A(:,i)/unity
        B(1,i)=(l(2)*k(3,i))-(l(3)*k(2,i))
        B(2,i)=(l(3)*k(1,i))-(l(1)*k(3,i))
        B(3,i)=(l(1)*k(2,i))-(l(2)*k(1,i))
        unity=sqrt((B(1,i)**2)+(B(2,i)**2)+(B(3,i)**2))
        B(:,i)=B(:,i)/unity
!
!  Now that we have our unit A's & B's we multiply them by the amplitudes
!  we defined earlier, to create the spectrum.
!
        A(:,i)=ampA(i)*A(:,i)
        B(:,i)=ampB(i)*B(:,i)
      enddo
!
      do i=1,KS_modes
!
!  These are used to define omega - the unsteadyness frequency (co-eff of t).
!
        arg=energy(i)*(kk(i)**3)
        if (arg>0.0)omega(i)=sqrt(arg)
        if (arg==0.0)omega(i)=0.0
      enddo
!
      do i=1,KS_modes
        call cross(A(:,i),unit_k(:,i),KS_A(:,i))
        call cross(B(:,i),unit_k(:,i),KS_B(:,i))
      enddo
!
      KS_omega(:)=omega(:)
      KS_k=k
!
!  Tidy up.
!
      deallocate(unit_k,orderk,kk,delk,energy,omega,klengths,ampA,ampB)
!
    endsubroutine periodic_KS_setup
!***********************************************************************
    subroutine KS_order(ad_, i_N, i_ord, B)
!
!  Bubble sort algorithm.
!
!  22-feb-10/abag: coded
!
      integer, intent(in) :: i_N, i_ord
      real, intent(in)  :: ad_(i_N)
      real, intent(out) :: B(i_N)
      real :: c=1.
      integer :: n
!
      B = ad_
!
      do while(c>0.0)
        c = -1.
        do n = 1, i_N-1
          if (   (i_ord==1 .and. B(n)>B(n+1)) &
            .or. (i_ord==2 .and. B(n)<B(n+1)) ) then
            c        = B(n)
            B(n)   = B(n+1)
            B(n+1) = c
            c = 1.
          endif
        enddo
      enddo
!
    endsubroutine KS_order
!***********************************************************************
    subroutine input_persistent_hydro(id,lun,done)
!
!  Read in the stored time of the next random phase calculation.
!
!  12-apr-08/axel: adapted from input_persistent_forcing
!
      integer :: id,lun
      logical :: done
!
      if (id==id_record_NOHYDRO_TPHASE) then
        read (lun) tphase_kinflow
        done=.true.
      elseif (id==id_record_NOHYDRO_PHASE1) then
        read (lun) phase1
        done=.true.
      elseif (id==id_record_NOHYDRO_PHASE2) then
        read (lun) phase2
        done=.true.
      elseif (id==id_record_NOHYDRO_LOCATION) then
        read (lun) location
        done=.true.
      elseif (id==id_record_NOHYDRO_TSFORCE) then
        read (lun) tsforce
        done=.true.
      endif
!
      if (lroot) print*,'input_persistent_hydro: ',tphase_kinflow
!
    endsubroutine input_persistent_hydro
!***********************************************************************
    subroutine output_persistent_hydro(lun)
!
!  Writes out the time of the next random phase calculation.
!
!  12-apr-08/axel: adapted from output_persistent_forcing
!
      integer :: lun
!
      if (lroot.and.ip<14) then
        if (tphase_kinflow>=0.) &
            print*,'output_persistent_hydro: ',tphase_kinflow
      endif
!
!  Write details.
!
      write (lun) id_record_NOHYDRO_TPHASE
      write (lun) tphase_kinflow
      write (lun) id_record_NOHYDRO_PHASE1
      write (lun) phase1
      write (lun) id_record_NOHYDRO_PHASE2
      write (lun) phase2
      write (lun) id_record_NOHYDRO_LOCATION
      write (lun) location
      write (lun) id_record_NOHYDRO_TSFORCE
      write (lun) tsforce
!
    endsubroutine output_persistent_hydro
!***********************************************************************
    subroutine read_hydro_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_hydro_init_pars
!***********************************************************************
    subroutine write_hydro_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_hydro_init_pars
!***********************************************************************
    subroutine read_hydro_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=hydro_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=hydro_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_hydro_run_pars
!***********************************************************************
    subroutine write_hydro_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=hydro_run_pars)
!
    endsubroutine write_hydro_run_pars
!***********************************************************************
    subroutine rprint_hydro(lreset,lwrite)
!
!  Reads and registers print parameters relevant for hydro part.
!
!   8-jun-02/axel: adapted from hydro
!
      use Diagnostics, only: parse_name
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  Reset everything in case of reset
!  (this needs to be consistent with what is defined above!).
!
      if (lreset) then
        idiag_u2m=0; idiag_um2=0; idiag_oum=0; idiag_o2m=0
        idiag_uxpt=0; idiag_uypt=0; idiag_uzpt=0; idiag_dtu=0
        idiag_urms=0; idiag_umax=0; idiag_uzrms=0; idiag_uzmax=0
        idiag_phase1=0; idiag_phase2=0
        idiag_orms=0; idiag_omax=0; idiag_oumphi=0
        idiag_ruxm=0; idiag_ruym=0; idiag_ruzm=0; idiag_rumax=0
        idiag_ux2m=0; idiag_uy2m=0; idiag_uz2m=0
        idiag_uxuym=0; idiag_uxuzm=0; idiag_uyuzm=0
        idiag_umx=0; idiag_umy=0; idiag_umz=0
        idiag_Marms=0; idiag_Mamax=0; idiag_divu2m=0; idiag_epsK=0
        idiag_urmphi=0; idiag_upmphi=0; idiag_uzmphi=0; idiag_u2mphi=0
        idiag_ekin=0; idiag_ekintot=0
        idiag_divum=0
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if (lroot.and.ip<14) print*,'rprint_hydro: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'ekin',idiag_ekin)
        call parse_name(iname,cname(iname),cform(iname),'ekintot',idiag_ekintot)
        call parse_name(iname,cname(iname),cform(iname),'u2m',idiag_u2m)
        call parse_name(iname,cname(iname),cform(iname),'um2',idiag_um2)
        call parse_name(iname,cname(iname),cform(iname),'o2m',idiag_o2m)
        call parse_name(iname,cname(iname),cform(iname),'oum',idiag_oum)
        call parse_name(iname,cname(iname),cform(iname),'dtu',idiag_dtu)
        call parse_name(iname,cname(iname),cform(iname),'urms',idiag_urms)
        call parse_name(iname,cname(iname),cform(iname),'umax',idiag_umax)
        call parse_name(iname,cname(iname),cform(iname),'uzrms',idiag_uzrms)
        call parse_name(iname,cname(iname),cform(iname),'uzmax',idiag_uzmax)
        call parse_name(iname,cname(iname),cform(iname),'ux2m',idiag_ux2m)
        call parse_name(iname,cname(iname),cform(iname),'uy2m',idiag_uy2m)
        call parse_name(iname,cname(iname),cform(iname),'uz2m',idiag_uz2m)
        call parse_name(iname,cname(iname),cform(iname),'uxuym',idiag_uxuym)
        call parse_name(iname,cname(iname),cform(iname),'uxuzm',idiag_uxuzm)
        call parse_name(iname,cname(iname),cform(iname),'uyuzm',idiag_uyuzm)
        call parse_name(iname,cname(iname),cform(iname),'orms',idiag_orms)
        call parse_name(iname,cname(iname),cform(iname),'omax',idiag_omax)
        call parse_name(iname,cname(iname),cform(iname),'divum',idiag_divum)
        call parse_name(iname,cname(iname),cform(iname),'ruxm',idiag_ruxm)
        call parse_name(iname,cname(iname),cform(iname),'ruym',idiag_ruym)
        call parse_name(iname,cname(iname),cform(iname),'ruzm',idiag_ruzm)
        call parse_name(iname,cname(iname),cform(iname),'rumax',idiag_rumax)
        call parse_name(iname,cname(iname),cform(iname),'umx',idiag_umx)
        call parse_name(iname,cname(iname),cform(iname),'umy',idiag_umy)
        call parse_name(iname,cname(iname),cform(iname),'umz',idiag_umz)
        call parse_name(iname,cname(iname),cform(iname),'Marms',idiag_Marms)
        call parse_name(iname,cname(iname),cform(iname),'Mamax',idiag_Mamax)
        call parse_name(iname,cname(iname),cform(iname),'divu2m',idiag_divu2m)
        call parse_name(iname,cname(iname),cform(iname),'epsK',idiag_epsK)
        call parse_name(iname,cname(iname),cform(iname),'uxpt',idiag_uxpt)
        call parse_name(iname,cname(iname),cform(iname),'uypt',idiag_uypt)
        call parse_name(iname,cname(iname),cform(iname),'uzpt',idiag_uzpt)
        call parse_name(iname,cname(iname),cform(iname),'phase1',idiag_phase1)
        call parse_name(iname,cname(iname),cform(iname),'phase2',idiag_phase2)
      enddo
!
!  Write column where which variable is stored.
!
      if (lwr) then
        write(3,*) 'nname=',nname
        write(3,*) 'iuu=',iuu
        write(3,*) 'iux=',iux
        write(3,*) 'iuy=',iuy
        write(3,*) 'iuz=',iuz
      endif
!
    endsubroutine rprint_hydro
!***********************************************************************
    subroutine get_slices_hydro(f,slices)
!
!  Write slices for animation of Hydro variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_hydro
!***********************************************************************
    subroutine calc_mflow
!
!  Dummy routine.
!
!  19-jul-03/axel: adapted from hydro
!
    endsubroutine calc_mflow
!***********************************************************************
    subroutine remove_mean_momenta(f)
!
!  Dummy routine.
!
!  32-nov-06/tobi: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine remove_mean_momenta
!***********************************************************************
    subroutine impose_velocity_ceiling(f)
!
!  13-aug-2007/anders: dummy
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine impose_velocity_ceiling
!***********************************************************************
    subroutine init_ck
!
!  8-sep-2009/dhruba: coded
!
      integer :: l,m
      real :: Balpha,jl,jlp1,jlm1,LPl,LPlm1
      integer :: ell
!
      print*, 'Initializing variables for Chandrasekhar-Kendall flow'
      print*, 'Allocating..'
      allocate(Zl(mx),dZldr(mx))
      allocate(Pl(my),dPldtheta(my))
      print*, 'Allocation done'
      ell=kinflow_ck_ell
      Balpha=kinflow_ck_Balpha
      print*, 'ell=,alpha=',ell,Balpha
!
      do l=1,mx
        call sp_besselj_l(jl,ell,Balpha*x(l))
        call sp_besselj_l(jlp1,ell+1,Balpha*x(l))
        call sp_besselj_l(jlm1,ell-1,Balpha*x(l))
        Zl(l) = jl
        dZldr(l) = ell*jlm1-(ell+1)*jlp1
      enddo
!
      do m=1,my
        call legendre_pl(LPl,ell,y(m))
        call legendre_pl(LPlm1,ell-1,y(m))
        Pl(m) = Lpl
        dPldtheta(m) = -(1/sin(y(m)))*ell*(LPlm1-LPl)
      enddo
!
    endsubroutine init_ck
!***********************************************************************
    subroutine hydro_clean_up
!
!  Deallocate the variables allocated in nohydro.
!
!  8-sep-2009/dhruba: coded
!
      print*, 'Deallocating some nohydro variables ...'
      if (kinflow=='ck') then
        deallocate(Zl,dZldr)
        deallocate(Pl,dPldtheta)
      elseif (kinflow=='KS') then
        deallocate(KS_k)
        deallocate(KS_A)
        deallocate(KS_B)
        deallocate(KS_omega)
      endif
!
      print*, 'Done.'
!
    endsubroutine hydro_clean_up
!***********************************************************************
    subroutine kinematic_random_phase
!
!  Get a random phase to be used for the whole kinematic velocity field.
!
!  16-feb-2010/dhruba: coded
!
      use General, only: random_number_wrapper
!
      real, dimension(3) :: fran
!
!  Generate random numbers.
!
      if (t>tsforce) then
        if (lrandom_location) then
          call random_number_wrapper(fran)
          location=fran*Lxyz+xyz0
        else
          location=location_fixed
        endif
!
!  Writing down the location.
!
        if (lroot .and. lwrite_random_location) then
          open(1,file=trim(datadir)//'/random_location.dat',status='unknown',position='append')
          write(1,'(4f14.7)') t, location
          close (1)
        endif
!
!  Update next tsforce.
!
        tsforce=t+dtforce
        if (ip<=6) print*,'kinematic_random_phase: location=',location
      endif
!
    endsubroutine kinematic_random_phase
!***********************************************************************
endmodule Hydro
