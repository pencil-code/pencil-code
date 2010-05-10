! $Id$
!
!  This module supplies a kinematic velocity field.
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lhydro = .false.
! CPARAM logical, parameter :: lhydro_kinematic = .true.
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED oo(3); ou; uij(3,3); uu(3); u2; sij(3,3)
! PENCILS PROVIDED der6u(3)
! PENCILS PROVIDED divu; uij5(3,3); graddivu(3)
!************************************************************************
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
!
  real :: u_out_kep=0.0
  real :: tphase_kinflow=-1.,phase1=0., phase2=0., tsforce=0., dtforce=impossible
  real, dimension(3) :: location,location_fixed=(/0.,0.,0./)
  logical :: lpressuregradient_gas=.true.,lcalc_uumean=.false.,lupw_uu=.false.
!
  real, allocatable, dimension (:,:) :: KS_k,KS_A,KS_B !or through whole field for each wavenumber?
  real, allocatable, dimension (:) :: KS_omega !or through whole field for each wavenumber?
  integer :: KS_modes = 3
  real, allocatable, dimension (:) :: Zl,dZldr,Pl,dPldtheta
  real :: ampl_fcont_uu=1.
  logical :: lforcing_cont_uu=.false., lrandom_location=.false., lwrite_random_location=.false.
!
!
!init parameters
!
!
!run parameters
!
  character (len=labellen) :: kinematic_flow='none'
  real :: ABC_A=1.0, ABC_B=1.0, ABC_C=1.0
  real :: wind_amp=0.,wind_rmin=impossible,wind_step_width=0.
  real :: circ_amp=0.,circ_rmax=0.,circ_step_width=0.
  real :: kx_uukin=1., ky_uukin=1., kz_uukin=1.
  real :: radial_shear=0.,uphi_at_rzero=0.,uphi_at_rmax=0.,uphi_rmax=1.,&
          uphi_step_width=0.
  character (len=labellen) :: wind_profile='none'
  namelist /hydro_run_pars/ &
    kinematic_flow,wind_amp,wind_profile,wind_rmin,wind_step_width, &
    circ_rmax,circ_step_width,circ_amp, &
    ABC_A,ABC_B,ABC_C, &
    ampl_kinflow,kx_uukin,ky_uukin,kz_uukin, &
    lrandom_location,lwrite_random_location,location_fixed,dtforce,&
    radial_shear,uphi_at_rzero,uphi_rmax,uphi_step_width
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
!
  contains
!***********************************************************************
    subroutine register_hydro()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use Mpicomm, only: lroot
      use SharedVariables
!
      integer :: ierr
!
!  Identify version number (generated automatically by SVN).
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Share lpressuregradient_gas so Entropy module knows whether to apply
!  pressure gradient or not.
!
      call put_shared_variable('lpressuregradient_gas',lpressuregradient_gas,ierr)
      if (ierr/=0) call fatal_error('register_hydro','there was a problem sharing lpressuregradient_gas')
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
!        call random_isotropic_KS_setup(-5./3.,1.,(nxgrid)/2.)
!
!  Use constant values for testing KS model code with 3
!  specific modes.
!
        call random_isotropic_KS_setup_test
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
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_hydro
!***********************************************************************
    subroutine init_uu(f)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Hydro module, if there was one.
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
! DM
! The following line can be later removed and the variable kinematic_flow replaced
! by kinflow.
      kinflow=kinematic_flow
      if (kinflow/='') then
        lpenc_requested(i_uu)=.true.
        if (kinflow=='eddy') then
          lpenc_requested(i_rcyl_mn)=.true.
          lpenc_requested(i_rcyl_mn1)=.true.
        endif
      endif
!
!  disgnostic pencils
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
!ajwm May be overkill... Perhaps only needed for certain kinflow?
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
      real, dimension(nx) :: wind_prof,div_uprof,der6_uprof
      real, dimension(nx) :: div_vel_prof
      real, dimension(nx) :: vel_prof
      real, dimension(nx) :: tmp_mn, cos1_mn, cos2_mn
      real, dimension(nx) :: rone, argx
      real :: fac, fac2, argy, argz
      real :: fpara, dfpara, ecost, esint, epst, sin2t, cos2t
      integer :: modeN
      real :: sqrt2, sqrt21k1, eps1=1., WW=0.25, k21
      integer :: ell
      real :: Balpha
      real :: theta,theta1
!      real :: omega0,shear
!
      intent(in) :: f
      intent(inout) :: p
!
!  choose from a list of different flow profiles.
!  Begin with a constant flow in the x direction.
!
      if (kinflow=='const-x') then
        if (lpencil(i_uu)) then
          if (headtt) print*,'const-x'
          p%uu(:,1)=ampl_kinflow*cos(omega_kinflow*t)*exp(eps_kinflow*t)
          p%uu(:,2)=0.
          p%uu(:,3)=0.
        endif
        if (lpencil(i_divu)) p%divu=0.
!  Begin with a constant flow in the (ampl_kinflow_x,ampl_kinflow_y,ampl_kinflow_z) direction.
!
      elseif (kinflow=='const-xyz') then
        if (lpencil(i_uu)) then
          if (headtt) print*,'const-xyz'
          p%uu(:,1)=ampl_kinflow_x*cos(omega_kinflow*t)*exp(eps_kinflow*t)
          p%uu(:,2)=ampl_kinflow_y*cos(omega_kinflow*t)*exp(eps_kinflow*t)
          p%uu(:,3)=ampl_kinflow_z*cos(omega_kinflow*t)*exp(eps_kinflow*t)
        endif
        if (lpencil(i_divu)) p%divu=0.
!
!  choose from a list of different flow profiles
!  ABC-flow
!
      elseif (kinflow=='ABC') then
! uu
        if (lpencil(i_uu)) then
          if (headtt) print*,'ABC flow'
          p%uu(:,1)=ABC_A*sin(kz_uukin*z(n))    +ABC_C*cos(ky_uukin*y(m))
          p%uu(:,2)=ABC_B*sin(kx_uukin*x(l1:l2))+ABC_A*cos(kz_uukin*z(n))
          p%uu(:,3)=ABC_C*sin(ky_uukin*y(m))    +ABC_B*cos(kx_uukin*x(l1:l2))
        endif
! divu
        if (lpencil(i_divu)) p%divu=0.
!
!  nocosine or Archontis flow
!
      elseif (kinflow=='nocos') then
! uu
        if (lpencil(i_uu)) then
          if (headtt) print*,'nocosine or Archontis flow'
          p%uu(:,1)=ABC_A*sin(kz_uukin*z(n))
          p%uu(:,2)=ABC_B*sin(kx_uukin*x(l1:l2))
          p%uu(:,3)=ABC_C*sin(ky_uukin*y(m))
        endif
! divu
        if (lpencil(i_divu)) p%divu=0.
!
!  Gen-Roberts flow (negative helicity)
!
      elseif (kinflow=='roberts') then
! uu
        sqrt2=sqrt(2.)
        if (lpencil(i_uu)) then
          if (headtt) print*,'Glen Roberts flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
          eps1=1.-eps_kinflow
          p%uu(:,1)=+eps1*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
          p%uu(:,2)=-eps1*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
          p%uu(:,3)=sqrt2*sin(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
        endif
! divu
        if (lpencil(i_divu)) p%divu= (kx_uukin-ky_uukin)*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
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
      elseif (kinflow=='ck') then
! uu
        ell=kinflow_ck_ell
        Balpha=kinflow_ck_Balpha
        if (lpencil(i_uu)) then
          if (headtt) print*,'Chandrasekhar-Kendall flow'
          p%uu(:,1)=ampl_kinflow*(ell*(ell+1)/Balpha*x(l1:l2))*Zl(l1:l2)*Pl(m)
          p%uu(:,2)=ampl_kinflow*(1./Balpha*x(l1:l2))*(2*x(l1:l2)*Zl(l1:l2)+&
                         dZldr(l1:l2)*Balpha*x(l1:l2)**2)*dPldtheta(m)
          p%uu(:,3)=-ampl_kinflow*Zl(l1:l2)*dPldtheta(m)
        endif
! divu
        if (lpencil(i_divu)) p%divu= 0.
!
!  Glen-Roberts flow (positive helicity)
!
      elseif (kinflow=='poshel-roberts') then
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
        if (headtt) print*,'Pos Helicity Roberts flow; eps1=',eps1
        p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
        p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        if (lpencil(i_divu)) p%divu=0.
!
!  1-D Glen-Roberts flow (positive helicity, no y-dependence)
!
      elseif (kinflow=='poshel-roberts-1d') then
        if (headtt) print*,'Pos Helicity Roberts flow; kx_uukin=',kx_uukin
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
        p%uu(:,1)=0.
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*eps1
        p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*sqrt(2.)
        if (lpencil(i_divu)) p%divu=0.
!
!  Glen-Roberts flow (x-direction, positive helicity)
!  x -> y
!  y -> z
!  z -> x
!
      elseif (kinflow=='xdir-roberts') then
        if (headtt) print*,'x-dir Roberts flow; ky_uukin,kz_uukin=',ky_uukin,kz_uukin
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
        p%uu(:,2)=-fac*cos(ky_uukin*y(m))*sin(kz_uukin*z(n))*eps1
        p%uu(:,3)=+fac*sin(ky_uukin*y(m))*cos(kz_uukin*z(n))*eps1
        p%uu(:,1)=+fac*cos(ky_uukin*y(m))*cos(kz_uukin*z(n))*sqrt(2.)
        if (lpencil(i_divu)) p%divu=0.
!
!  z-dependent Roberts flow (positive helicity)
!
      elseif (kinflow=='zdep-roberts') then
        if (headtt) print*,'z-dependent Roberts flow; kx,ky=',kx_uukin,ky_uukin
        fpara=ampl_kinflow*(quintic_step(z(n),-1.+eps_kinflow,eps_kinflow) &
                           -quintic_step(z(n),+1.-eps_kinflow,eps_kinflow))
        dfpara=ampl_kinflow*(quintic_der_step(z(n),-1.+eps_kinflow,eps_kinflow)&
                            -quintic_der_step(z(n),+1.-eps_kinflow,eps_kinflow))
!
        sqrt2=sqrt(2.)
        sqrt21k1=1./(sqrt2*kx_uukin)
!
        p%uu(:,1)=-cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m)) &
           -dfpara*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt21k1
        p%uu(:,2)=+sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m)) &
           -dfpara*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*sqrt21k1
        p%uu(:,3)=+fpara*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt2
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with cosinusoidal helicity variation (version 1)
!
      elseif (kinflow=='IncohRoberts1') then
        if (headtt) print*,'Roberts flow with cosinusoidal helicity; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        eps1=(1.-eps_kinflow)*cos(omega_kinflow*t)
        p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
        p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with cosinusoidal helicity variation (version 2)
!
      elseif (kinflow=='IncohRoberts2') then
        if (headtt) print*,'Roberts flow with cosinusoidal helicity; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        eps1=(1.-eps_kinflow)*cos(omega_kinflow*t)
        p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
        p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)*eps1
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with helicity variation and shear (version 1)
!  Must use shear module and set eps_kinflow equal to shear
!
      elseif (kinflow=='ShearRoberts1') then
        if (headtt) print*,'Roberts flow with cosinusoidal helicity; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        ky_uukin=1.
        kx_uukin=ky_uukin*(mod(.5-eps_kinflow*t,1.D0)-.5)
if (ip.eq.11.and.m==4.and.n==4) write(21,*) t,kx_uukin
        eps1=cos(omega_kinflow*t)
        p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
        p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with helicity variation and shear (version 2)
!  Must use shear module and set eps_kinflow equal to shear
!
      elseif (kinflow=='ShearRoberts1') then
        if (headtt) print*,'Roberts flow with cosinusoidal helicity; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        ky_uukin=1.
        kx_uukin=ky_uukin*(mod(.5-eps_kinflow*t,1.D0)-.5)
if (ip.eq.11.and.m==4.and.n==4) write(21,*) t,kx_uukin
        eps1=cos(omega_kinflow*t)
        p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))
        p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)*eps1
        if (lpencil(i_divu)) p%divu=0.
!
!  Taylor-Green flow
!
      elseif (kinflow=='TG') then
        if (headtt) print*,'Taylor-Green flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        p%uu(:,1)=+2.*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*cos(kz_uukin*z(n))
        p%uu(:,2)=-2.*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*cos(kz_uukin*z(n))
        p%uu(:,3)=+0.
        if (lpencil(i_divu)) p%divu=0.
!
!  Galloway-Proctor flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!
      elseif (kinflow=='Galloway-Proctor') then
        if (headtt) print*,'Galloway-Proctor flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
        p%uu(:,1)=-fac*sin(ky_uukin*y(m)    +esint)
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)
        p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+ecost)+cos(ky_uukin*y(m)+esint))
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Galloway-Proctor flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!
      elseif (kinflow=='Galloway-Proctor-nohel') then
        if (headtt) print*,'nonhelical Galloway-Proctor flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow*sqrt(1.5)
        fac2=ampl_kinflow*sqrt(6.)
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
        p%uu(:,1)=+fac*cos(ky_uukin*y(m)    +esint)
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)
        p%uu(:,3)=-fac2*(sin(kx_uukin*x(l1:l2)+ecost)*cos(ky_uukin*y(m)+esint))
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Otani flow, U=curl(psi*zz) + psi*zz, where
!  psi = 2*cos^2t * cosx - 2*csin2t * cosy
!
      elseif (kinflow=='Otani') then
        if (headtt) print*,'Otani flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=2.*ampl_kinflow
        sin2t=sin(omega_kinflow*t)**2
        cos2t=cos(omega_kinflow*t)**2
        p%uu(:,1)=fac*sin2t*sin(ky_uukin*y(m))
        p%uu(:,2)=fac*cos2t*sin(kx_uukin*x(l1:l2))
        p%uu(:,3)=fac*(cos2t*cos(kx_uukin*x(l1:l2))-sin2t*cos(ky_uukin*y(m)))
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=p%uu
!
!  Tilgner flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!
      elseif (kinflow=='Tilgner') then
        if (headtt) print*,'Tilgner flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow*sqrt(2.)
        epst=eps_kinflow*t
kx_uukin=2.*pi
ky_uukin=2.*pi
        p%uu(:,1)=-fac* sin(ky_uukin*y(m)         )
        p%uu(:,2)=+fac* sin(kx_uukin*(x(l1:l2)+epst))
        p%uu(:,3)=-fac*(cos(kx_uukin*(x(l1:l2)+epst))+cos(ky_uukin*y(m)))
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Tilgner flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kx_uukin=ky_uukin
!  Here, W in Tilgner's equation is chosen to be 0.25.
!
      elseif (kinflow=='Tilgner-orig') then
        if (headtt) print*,'original Tilgner flow; kx_uukin,ky_uukin=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        epst=eps_kinflow*t
        kx_uukin=2.*pi
        ky_uukin=2.*pi
        sqrt2=sqrt(2.)
        WW=0.25
        p%uu(:,1)=+fac*sqrt2*sin(kx_uukin*(x(l1:l2)-epst))*cos(ky_uukin*y(m))
        p%uu(:,2)=-fac*sqrt2*cos(kx_uukin*(x(l1:l2)-epst))*sin(ky_uukin*y(m))
        p%uu(:,3)=+fac*2.*WW*sin(kx_uukin*(x(l1:l2)-epst))*sin(ky_uukin*y(m))
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kx_uukin*p%uu
!
!  Galloway-Proctor flow with random temporal phase
!
      elseif (kinflow=='Galloway-Proctor-RandomTemporalPhase') then
        if (headtt) print*,'GP-RandomTemporalPhase; kx,ky=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        if (t.gt.tphase_kinflow) then
          call random_number_wrapper(fran1)
          tphase_kinflow=t+dtphase_kinflow
          phase1=pi*(2*fran1(1)-1.)
          phase2=pi*(2*fran1(2)-1.)
        endif
        ecost=eps_kinflow*cos(omega_kinflow*t+phase1)
        esint=eps_kinflow*sin(omega_kinflow*t+phase2)
        p%uu(:,1)=-fac*sin(ky_uukin*y(m)    +esint)
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)
        p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+ecost)+cos(ky_uukin*y(m)+esint))
        if (lpencil(i_divu)) p%divu=0.
!
!  Galloway-Proctor flow with random phase
!
      elseif (kinflow=='Galloway-Proctor-RandomPhase') then
        if (headtt) print*,'Galloway-Proctor-RandomPhase; kx,ky=',kx_uukin,ky_uukin
        fac=ampl_kinflow
        if (t.gt.tphase_kinflow) then
          call random_number_wrapper(fran1)
          tphase_kinflow=t+dtphase_kinflow
          phase1=eps_kinflow*pi*(2*fran1(1)-1.)
          phase2=eps_kinflow*pi*(2*fran1(2)-1.)
        endif
        p%uu(:,1)=-fac*sin(ky_uukin*y(m)    +phase1)*ky_uukin
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+phase2)*kx_uukin
        p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+phase2)+cos(ky_uukin*y(m)+phase1))
        if (lpencil(i_divu)) p%divu=0.
!
!  original Galloway-Proctor flow
!
      elseif (kinflow=='Galloway-Proctor-orig') then
        if (headtt) print*,'Galloway-Proctor-orig flow; kx_uukin=',kx_uukin
        fac=sqrt(1.5)*ampl_kinflow
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
        p%uu(:,1)=+fac*cos(ky_uukin*y(m)    +esint)*ky_uukin
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2)+ecost)*kx_uukin
        p%uu(:,3)=-fac*(cos(kx_uukin*x(l1:l2)+ecost)+sin(ky_uukin*y(m)+esint))
        if (lpencil(i_divu)) p%divu=0.
!
!
!potential flow, u=gradphi, with phi=cosx*cosy*cosz
!  assume kx_uukin=ky_uukin=kz_uukin
!
      elseif (kinflow=='potential') then
        fac=ampl_kinflow
        if (headtt) print*,'potential; kx_uukin,ampl_kinflow=',ampl_kinflow
        if (headtt) print*,'potential; kx_uukin=',kx_uukin,ky_uukin,kz_uukin
        p%uu(:,1)=-fac*kx_uukin*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*cos(kz_uukin*z(n))
        p%uu(:,2)=-fac*ky_uukin*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*cos(kz_uukin*z(n))
        p%uu(:,3)=-fac*kz_uukin*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sin(kz_uukin*z(n))
        if (lpencil(i_divu)) p%divu=-fac*(kx_uukin**2+ky_uukin**2*kz_uukin**2) &
          *cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*cos(kz_uukin*z(n))
!
!potential random flow, u=gradphi, with phi=cos(x-x0)*cosy*cosz
!  assume kx_uukin=ky_uukin=kz_uukin
!
!
      elseif (kinflow=='potential_random') then
        fac=ampl_kinflow
        if (headtt) print*,'potential_random; kx_uukin,ampl_kinflow=',ampl_kinflow
        if (headtt) print*,'potential_random; kx_uukin=',kx_uukin,ky_uukin,kz_uukin
        argx=kx_uukin*(x(l1:l2)-location(1))
        argy=ky_uukin*(y(m)-location(2))
        argz=kz_uukin*(z(n)-location(3))
        p%uu(:,1)=-fac*kx_uukin*sin(argx)*cos(argy)*cos(argz)
        p%uu(:,2)=-fac*ky_uukin*cos(argx)*sin(argy)*cos(argz)
        p%uu(:,3)=-fac*kz_uukin*cos(argx)*cos(argy)*sin(argz)
        if (lpencil(i_divu)) p%divu=fac
!
!  Convection rolls
!  Stream function: psi_y = cos(kx*x) * cos(kz*z)
!
      elseif (kinflow=='rolls') then
        if (headtt) print*,'Convection rolls; kx_kinflow,kz_uukin=',kx_kinflow,kz_kinflow
        p%uu(:,1)=ampl_kinflow*kz_kinflow*cos(kx_kinflow*x(l1:l2))*sin(kz_kinflow*z(n))
        p%uu(:,2)=+0.
        p%uu(:,3)=ampl_kinflow*kx_kinflow*sin(kx_kinflow*x(l1:l2))*cos(kz_kinflow*z(n))
! divu
        if (lpencil(i_divu)) p%divu=0.
!
!  Twist (Yousef & Brandenburg 2003)
!
      elseif (kinflow=='Twist') then
        if (headtt) print*,'Twist flow; eps_kinflow,kx=',eps_kinflow,kx_uukin
        p%uu(:,1)=0.
        p%uu(:,2)=eps_kinflow*z(n)*sin(kx_uukin*x(l1:l2))
        p%uu(:,3)=1.+cos(kx_uukin*x(l1:l2))
        if (lpencil(i_divu)) p%divu=0.
!
!  Eddy (Brandenburg & Zweibel 1994)
!
      elseif (kinflow=='eddy') then
        if (headtt) print*,'eddy flow; eps_kinflow,kx=',eps_kinflow,kx_uukin
        cos1_mn=max(cos(.5*pi*p%rcyl_mn),0.)
        cos2_mn=cos1_mn**2
        tmp_mn=-.5*pi*p%rcyl_mn1*sin(.5*pi*p%rcyl_mn)*ampl_kinflow* &
          4.*cos2_mn*cos1_mn
        p%uu(:,1)=+tmp_mn*y(m)
        p%uu(:,2)=-tmp_mn*x(l1:l2)
        p%uu(:,3)=eps_kinflow*ampl_kinflow*cos2_mn**2
        if (lpencil(i_divu)) p%divu=0.
!
!  Shearing wave
!
      elseif (kinflow=='ShearingWave') then
        if (headtt) print*,'ShearingWave flow; Sshear,eps_kinflow=',Sshear,eps_kinflow
        ky_uukin=1.
        kx_uukin=-ky_uukin*Sshear*t
        k21=1./(kx_uukin**2+ky_uukin**2)
        p%uu(:,1)=-ky_uukin*k21*cos(kx_uukin*x(l1:l2)+ky_uukin*y(m))
        p%uu(:,2)=+kx_uukin*k21*cos(kx_uukin*x(l1:l2)+ky_uukin*y(m))
        p%uu(:,3)=eps_kinflow*cos(kx_uukin*x(l1:l2)+ky_uukin*y(m))
        if (lpencil(i_divu)) p%divu=0.
!
!  Helical Shearing wave
!
      elseif (kinflow=='HelicalShearingWave') then
        if (headtt) print*,'ShearingWave flow; Sshear,eps_kinflow=',Sshear,eps_kinflow
        ky_uukin=1.
        kx_uukin=-ky_uukin*Sshear*t
        k21=1./(kx_uukin**2+ky_uukin**2)
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
        p%uu(:,1)=-fac*cos(kx_uukin*x(l1:l2))*sin(ky_uukin*y(m))*eps1
        p%uu(:,2)=+fac*sin(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*eps1
        p%uu(:,3)=+fac*cos(kx_uukin*x(l1:l2))*cos(ky_uukin*y(m))*sqrt(2.)
        if (lpencil(i_divu)) p%divu=0.
!
! step function along z
!
      elseif (kinflow=='zstep') then
        if (lpencil(i_uu)) then
          if (headtt) print*,'wind:step function along z'
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=wind_amp*step_scalar(z(n),wind_rmin,wind_step_width)
        endif
!
! uniform radial shear with a cutoff at rmax
! uniform radial shear with saturation.
!
      elseif (kinflow=='rshear-sat') then
        if (lpencil(i_uu)) then
          if (headtt) print*,'radial shear saturated'
          p%uu(:,1)=0.
          p%uu(:,2)=0.
!          omega0=uphi_at_rmax/uphi_rmax
!          shear = arctanh(omega0)/(uphi_rmax-x(l1))
          p%uu(:,3)=x(l1:l2)*sinth(m)*tanh(10.*(x(l1:l2)-x(l1)))
!              (1+stepdown(x(l1:l2),uphi_rmax,uphi_step_width))
        endif
!
!  Vertical wind
!
      elseif (kinflow=='vertical-wind') then
        if (lpencil(i_uu)) then
          if (headtt) print*,'vertical-wind along z'
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=wind_amp*z(n)
        endif
        if (lpencil(i_divu)) p%divu=wind_amp
!
!  Vertical wind that goes to zero for z < 0.
!
      elseif (kinflow=='vertical-wind-surface') then
        if (lpencil(i_uu)) then
          if (headtt) print*,'vertical-wind along z'
          p%uu(:,1)=0.
          p%uu(:,2)=0.
          p%uu(:,3)=wind_amp*z(n)*0.5*(1.+erfunc((z(n)/wind_step_width)))
        endif
        if (lpencil(i_divu)) p%divu=wind_amp
!
! Radial wind
!
      elseif (kinflow=='radial-wind') then
        select case (wind_profile)
        case ('none'); wind_prof=0.;div_uprof=0.
        case ('constant'); wind_prof=1.;div_uprof=0.
!
        case ('radial-step')
          wind_prof=step(x(l1:l2),wind_rmin,wind_step_width)
          div_uprof=der_step(x(l1:l2),wind_rmin,wind_step_width)
!          der6_uprof=der6_step(x(l1:l2),wind_rmin,wind_step_width)
          der6_uprof=0.
        case default;
          call fatal_error('hydro_kinematic:kinflow = radial wind', 'no such wind profile. ')
        endselect
!
        if (lpencil(i_uu)) then
          if (headtt) print*,'Radial wind'
          p%uu(:,1)=wind_amp*wind_prof
          p%uu(:,2)=0.
          p%uu(:,3)=0.
        endif
        p%divu=wind_amp*div_uprof
        p%der6u(:,1)=wind_amp*der6_uprof
        p%der6u(:,2)= 0.
        p%der6u(:,3)= 0.
!
! Radial wind+circulation
!
      elseif (kinflow=='radial-wind+circ') then
        select case (wind_profile)
        case ('none'); wind_prof=0.;div_uprof=0.
        case ('constant'); wind_prof=1.;div_uprof=0.
!
        case ('radial-step')
          wind_prof=step(x(l1:l2),wind_rmin,wind_step_width)
          div_uprof=der_step(x(l1:l2),wind_rmin,wind_step_width)
!          der6_uprof=der6_step(x(l1:l2),wind_rmin,wind_step_width)
          der6_uprof=0.
        case default;
          call fatal_error('hydro_kinematic: kinflow= radial_wind-circ', 'no such wind profile. ')
        endselect
!
        rone=xyz0(1)
        theta=y(m)
        theta1=xyz0(2)
        if (lpencil(i_uu)) then
          if (headtt) print*,'Radial wind and circulation: not complete yet'
          vel_prof=circ_amp*(1+stepdown(x(l1:l2),circ_rmax,circ_step_width))
          div_vel_prof=-der_step(x(l1:l2),circ_rmax,circ_step_width)
!
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
! radial shear + radial wind + circulation
!
      elseif (kinflow=='rshear-sat+rwind+circ') then
!
! first set wind and cirulation
!
        select case (wind_profile)
        case ('none'); wind_prof=0.;div_uprof=0.
        case ('constant'); wind_prof=1.;div_uprof=0.
!
        case ('radial-step')
          wind_prof=step(x(l1:l2),wind_rmin,wind_step_width)
          div_uprof=der_step(x(l1:l2),wind_rmin,wind_step_width)
!          der6_uprof=der6_step(x(l1:l2),wind_rmin,wind_step_width)
          der6_uprof=0.
        case default;
          call fatal_error('hydro_kinematic: kinflow= radial-shear+rwind+rcirc', 'no such wind profile. ')
        endselect
!
        rone=xyz0(1)
        theta=y(m)
        theta1=xyz0(2)
        if (lpencil(i_uu)) then
          if (headtt) print*,'radial shear, wind, circulation: not complete yet'
          vel_prof=circ_amp*(1+stepdown(x(l1:l2),circ_rmax,circ_step_width))
          div_vel_prof=-der_step(x(l1:l2),circ_rmax,circ_step_width)
!
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
      elseif (kinflow=='KS') then
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
      else
! uu
        if (lpencil(i_uu)) then
          if (headtt) print*,'uu=0'
          p%uu=0.
        endif
! divu
        if (lpencil(i_divu)) p%divu=0.
      endif
! u2
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
!
      intent(in)  :: df,p
      intent(out) :: f
!
!  uu/dx for timestep (if kinflow is set)
!
      if (kinflow/='') then
        if (lfirst.and.ldt) advec_uu=abs(p%uu(:,1))*dx_1(l1:l2)+ &
                                     abs(p%uu(:,2))*dy_1(  m  )+ &
                                     abs(p%uu(:,3))*dz_1(  n  )
      endif
      if (headtt.or.ldebug) print*, 'duu_dt: max(advec_uu) =', maxval(advec_uu)
!
!  Store uu in auxiliary variable if requested.
!  Just neccessary immediately before writing snapshots, but how would we
!  know we are?
!
     if (lkinflow_as_aux) f(l1:l2,m,n,iux:iuz)=p%uu
!
!  Calculate maxima and rms values for diagnostic purposes
!
      if (ldiagnos) then
        if (headtt.or.ldebug) print*,'duu_dt: diagnostics ...'
        if (idiag_oum/=0) call sum_mn_name(p%ou,idiag_oum)
        !if (idiag_orms/=0) call sum_mn_name(p%o2,idiag_orms,lsqrt=.true.)
        !if (idiag_omax/=0) call max_mn_name(p%o2,idiag_omax,lsqrt=.true.)
!
!  kinetic field components at one point (=pt)
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
!  coriolis terms for cartesian geometry
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
!  dummy routine
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lhydro_pars
!***********************************************************************
    subroutine random_isotropic_KS_setup_tony(initpower,kmin,kmax)
!
!   produces random, isotropic field from energy spectrum following the
!   KS method (Malik and Vassilicos, 1999.)
!
!   more to do; unsatisfactory so far - at least for a steep power-law
!   energy spectrum
!
!   27-may-05/tony: modified from snod's KS hydro initial
!   03-feb-06/weezy: Tony's code doesn't appear to have the
!                    correct periodicity.
!                    renamed from random_isotropic_KS_setup
!
    use Sub, only: cross
    use General, only: random_number_wrapper
!
    integer :: modeN
!
    real, dimension (3) :: k_unit
    real, dimension (3) :: e1,e2
    real,dimension (6) :: r
    real,dimension (3) ::j,l  !get rid of this - these replace ee,ee1
    real :: initpower,kmin,kmax
    real, dimension(KS_modes) :: k,dk,energy,ps
    real :: theta,phi,alpha,beta
    real :: a,mkunit
    real :: newthet,newphi  !get rid of this line if there's no change
!
    allocate(KS_k(3,KS_modes))
    allocate(KS_A(3,KS_modes))
    allocate(KS_B(3,KS_modes))
    allocate(KS_omega(KS_modes))
!
!    minlen=Lxyz(1)/(nx-1)
!    kmax=2.*pi/minlen
!    KS_modes=int(0.5*(nx-1))
!    hh=Lxyz(1)/(nx-1)
!    pta=(nx)**(1.0/(nx-1))
!    do modeN=1,KS_modes
!       ggt=(kkmax-kkmin)/(KS_modes-1)
!       ggt=(kkmax/kkmin)**(1./(KS_modes-1))
!        k(modeN)=kmin+(ggt*(modeN-1))
!        k(modeN)=(modeN+3)*2*pi/Lxyz(1)
!       k(modeN)=kkmin*(ggt**(modeN-1)
!    enddo
!
!    do modeN=1,KS_modes
!       if (modeN.eq.1)delk(modeN)=(k(modeN+1)-K(modeN))
!       if (modeN.eq.KS_modes)delk(modeN)=(k(modeN)-k(modeN-1))
!       if (modeN.gt.1.and.modeN.lt.KS_modes)delk(modeN)=(k(modeN+1)-k(modeN-2))/2.0
!    enddo
!          mk=(k2*k2)*((1.0 + (k2/(bk_min*bk_min)))**(0.5*initpower-2.0))
!
!  set kmin
!
       kmin=2.*pi      !/(1.0*Lxyz(1))
!       kmin=kmin*2.*pi
       kmax=128.*pi    !nx*pi
       a=(kmax/kmin)**(1./(KS_modes-1.))
!
!
    do modeN=1,KS_modes
!
!  pick wavenumber
!
!       k=modeN*kmin
      k=kmin*(a**(modeN-1.))
!
!  calculate dk
!
!       print *,kmin,kmax,k
!       dk=1.0*kmin
!
      if (modeN==1)&
              dk=kmin*(a-1.)/2.
      if (modeN.gt.1.and.modeN.lt.KS_modes) &
              dk=(a**(modeN-2.))*kmin*((a**2.) -1.)/2.
      if (modeN==KS_modes) &
              dk=(a**(KS_modes -2.))*kmin*(a -1.)/2.
!
       call random_number_wrapper(r)
       theta=r(1)*pi
       phi=r(2)*2.0*pi
       alpha=r(3)*pi
       beta=r(4)*2.0*pi
       newthet=r(5)*pi
       newphi=r(6)*2.0*pi
!
       k_unit(1)=sin(theta)*cos(phi)
       k_unit(2)=sin(theta)*sin(phi)
       k_unit(3)=cos(theta)
!
       j(1)=sin(alpha)*cos(beta)
       j(2)=sin(alpha)*sin(beta)
       j(3)=cos(alpha)
!
       l(1)=sin(newthet)*cos(newphi)
       l(2)=sin(newthet)*sin(newphi)
       l(3)=cos(newthet)
!
       KS_k(:,modeN)=k*k_unit(:)
!
       call cross(KS_k(:,modeN),j,e1)
       call cross(KS_k(:,modeN),l,e2)
!
!  Make e1 & e2 unit vectors so that we can later make them
!  the correct lengths
!
       mkunit=sqrt(e1(1)**2+e1(2)**2+e1(3)**2)
       e1=e1/mkunit
!
       mkunit=sqrt(e2(1)**2+e2(2)**2+e2(3)**2)
       e2=e2/mkunit
!
!        energy=(((k/1.)**2. +1.)**(-11./6.))*(k**2.) &
!                            *exp(-0.5*(k/kmax)**2.)
!  The energy above is how this code has it. i
!  I've changed the divisor of k.
       energy=(((k/kmin)**2. +1.)**(-11./6.))*(k**2.) &
                       *exp(-0.5*(k/kmax)**2.)
       energy=1.*energy
       ps=sqrt(2.*energy*dk)   !/3.0)
!
       KS_A(:,modeN) = ps*e1
       KS_B(:,modeN) = ps*e2
!
    enddo
!
!   form RA = RA x k_unit and RB = RB x k_unit
!
    do modeN=1,KS_modes
      call cross(KS_A(:,modeN),k_unit(:),KS_A(:,modeN))
      call cross(KS_B(:,modeN),k_unit(:),KS_B(:,modeN))
    enddo
!
    call keep_compiler_quiet(initpower)
!
    endsubroutine random_isotropic_KS_setup_tony
!***********************************************************************
    subroutine random_isotropic_KS_setup(initpower,kmin,kmax)
!
!   produces random, isotropic field from energy spectrum following the
!   KS method (Malik and Vassilicos, 1999.)
!
!   more to do; unsatisfactory so far - at least for a steep power-law
!   energy spectrum
!
!   27-may-05/tony: modified from snod's KS hydro initial
!   03-feb-06/weezy: Attempted rewrite to guarantee periodicity of
!                    KS modes.
!
    use Sub, only: cross, dot2
    use General, only: random_number_wrapper
!
    integer :: modeN
!
    real, dimension (3) :: k_unit
    real, dimension (3) :: ee,e1,e2
!    real, dimension (4) :: r
    real,dimension (6) :: r
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
!    kmin=kmin*2.*pi
    kmax=128.*pi    !nx*pi
    a=(kmax/kmin)**(1./(KS_modes-1.))
!
!
    do modeN=1,KS_modes
!
!  pick wavenumber
!
!      k=modeN*kmin
      k=kmin*(a**(modeN-1.))
!
!weezy need to investigate if this is still needed
!weezy !
!weezy !  calculate dk
!weezy !
!weezy       print *,kmin,kmax,k
!weezy       dk=1.0*kmin
!
!weezy       if (modeN==1)dk=kmin*(a-1.)/2.
!weezy       if (modeN.gt.1.and.modeN.lt.KS_modes)dk=(a**(modeN-2.))*kmin*((a**2.) -1.)/2.
!weezy       if (modeN==KS_modes)dk=(a**(KS_modes -2.))*kmin*(a -1.)/2.
!
!
!  pick 4 random angles for each mode
!
      call random_number_wrapper(r);
      theta=pi*(r(1) - 0.)
      phi=pi*(2*r(2) - 0.)
      alpha=pi*(2*r(3) - 0.)
      beta=pi*(2*r(4) - 0.)
!
!  random phase?
!      call random_number_wrapper(r); gamma(modeN)=pi*(2*r - 0.)
!
!  make a random unit vector by rotating fixed vector to random position
!  (alternatively make a random transformation matrix for each k)
!
      k_unit(1)=sin(theta)*cos(phi)
      k_unit(2)=sin(theta)*sin(phi)
      k_unit(3)=cos(theta)
!
      energy=(((k/kmin)**2. +1.)**(-11./6.))*(k**2.) &
                       *exp(-0.5*(k/kmax)**2.)
!      energy=(((k/1.)**2. +1.)**(-11./6.))*(k**2.) &
!                       *exp(-0.5*(k/kmax)**2.)
!
!  make a vector KS_k of length k from the unit vector for each mode
!
      KS_k(:,modeN)=k*k_unit(:)
!      KS_omega(modeN)=k**(2./3.)
      KS_omega(:)=sqrt(energy(:)*(k(:)**3.))
!
!  construct basis for plane having rr normal to it
!  (bit of code from forcing to construct x', y')
!
      if ((k_unit(2).eq.0).and.(k_unit(3).eq.0)) then
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
!  make two random unit vectors KS_B and KS_A in the constructed plane
!
      KS_A(:,modeN) = cos(alpha)*e1 + sin(alpha)*e2
      KS_B(:,modeN) = cos(beta)*e1  + sin(beta)*e2
!
!  define the power spectrum (ps=sqrt(2.*power_spectrum(k)*delta_k/3.))
!
!      ps=(k**(initpower/2.))*sqrt(dk*2./3.)
!  The factor of 2 just after the sqrt may need to be 2./3.
!
!
!  With the `weezey' stuff above commented out, dk is currently used, but
!  never set, so we better abort
!
      call error('random_isotropic_KS_setup', 'Using uninitialized dk')
      dk=0.                     ! to make compiler happy
!
      ps=sqrt(2.*energy*dk)   !/3.0)
!
!  give KS_A and KS_B length ps
!
      KS_A(:,modeN)=ps*KS_A(:,modeN)
      KS_B(:,modeN)=ps*KS_B(:,modeN)
!
    enddo
!
!  form RA = RA x k_unit and RB = RB x k_unit
!  Note: cannot reuse same vector for input and output
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
!   produces random, isotropic field from energy spectrum following the
!   KS method (Malik and Vassilicos, 1999.)
!   This test case only uses 3 very specific modes (useful for comparison
!   with Louise's kinematic dynamo code.
!
!   03-feb-06/weezy: modified from random_isotropic_KS_setup
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
!-----------------------------
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
!-----------------------------
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
      if (modeN.gt.1.and.modeN.lt.KS_modes) &
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
!   form RA = RA x k_unit and RB = RB x k_unit
!
!
     do modeN=1,KS_modes
       call cross(KS_A(:,modeN),k_unit(:,modeN),KS_A(:,modeN))
       call cross(KS_B(:,modeN),k_unit(:,modeN),KS_B(:,modeN))
     enddo
!
    endsubroutine random_isotropic_KS_setup_test
!***********************************************************************
   ! subroutine random_isotropic_KS_setup_abag
!
!  ! produces random, isotropic field from energy spectrum following the
!  ! KS method, however this setup produces periodic velocity field
!  ! (assuming box (-pi,pi))
!
!  ! 28-mar-08/abag coded
!
   ! use Sub
   ! use General
   ! implicit none
   ! real,allocatable,dimension(:,:) :: unit_k,k,A,B,orderK
   ! real,allocatable,dimension(:) :: kk,delk,energy,omega,klengths
   ! real,dimension(3) :: angle,dir_in,u
   ! real :: k_option(3,10000),mkunit(10000)
   ! real :: arg
   ! real :: turn1,turnN
   ! integer ::i,s1,num,direction(3)
   ! logical :: ne
!
   ! allocate(KS_k(3,KS_modes))
   ! allocate(KS_A(3,KS_modes))
   ! allocate(KS_B(3,KS_modes))
   ! allocate(unit_k(3,KS_modes))
   ! allocate(k(3,KS_modes))
   ! allocate(A(3,KS_modes))
   ! allocate(B(3,KS_modes))
   ! allocate(orderk(3,KS_modes))
   ! allocate(KS_omega(KS_modes))
   ! allocate(kk(KS_modes))
   ! allocate(delk(KS_modes))
   ! allocate(energy(KS_modes))
   ! allocate(omega(KS_modes))
   ! allocate(klengths(KS_modes))
   ! num=1
   ! do i=1,10000
   !  call random_number(angle)
   !  if ((angle(1)-0.0 < epsilon(0.0)) .or. &
   !     (angle(2)-0.0 < epsilon(0.0)) .or. &
   !     (angle(3)-0.0 < epsilon(0.0))) then
   !     call random_number(angle)
   !  endif
   !  angle=floor(9.*angle)
   !  call random_number(dir_in)
   !  direction=nint(dir_in)
   !  direction=2*direction -1  !positive or negative directions
   !
   !  k_option(1,i)=direction(1)*angle(1)!a possible orientation
   !  k_option(2,i)=direction(2)*angle(2)   !provided we haven't
   !  k_option(3,i)=direction(3)*angle(3)  !already got this length
!
   !  !find the length of the current k_option vector
   !  mkunit(i)=dsqrt((k_option(1,i)**2)+(k_option(2,i)**2)+(k_option(3,i)**2))
!
   !  if (i==1.and.mkunit(i).gt.0.)then
   !    k(:,num)=k_option(:,i)
   !    klengths(num)=mkunit(i)
   !  endif
!
   !  !now we check that the current length is unique (hasn't come before)
   !  if (i.gt.1.and.num.lt.KS_modes)then
   !    do s1=i-1,1,-1
   !      if (mkunit(i).gt.0.0D0.and.mkunit(i) /= mkunit(s1))then
   !        ne=.true.
   !      else
   !        ne=.false.
   !        exit
   !      endif
   !      if (s1==1.and.ne)then !i.e. if length of current k_option is new......
   !        num=num+1
   !        k(:,num)=k_option(:,i) !load current k_option into k that we keep
   !        klengths(num)=mkunit(i)  ! store the length also
   !      endif
   !    enddo
   !   endif
   !   if (i==10000.and.num.lt.KS_modes)print*,"Haven't got",KS_modes,"modes!!!!"
   ! enddo
   ! do i=1,KS_modes
   !    do s1=1,KS_modes
   !       if (kk(i)==klengths(s1))then
   !          orderK(:,i)=k(:,s1)
   !       endif
   !    enddo
   ! enddo
   ! k=orderK
   ! do i=1,KS_modes
   !   unit_k(:,i)=k(:,i)/kk(i)
   ! enddo
   ! do i=1,N
   ! !now we find delk as defined in Malik & Vassilicos' paper
   !    if (i==1)delk(i)=(kk(i+1)-kk(i))/2.0D0
   !    if (i==KS_modes)delk(i)=(kk(i)-kk(i-1))/2.0D0
   !    if (i.gt.1.and.i.lt.KS_modes)delk(i)=(kk(i+1)-kk(i-1))/2.0D0
   ! enddo
   ! endsubroutine random_isotropic_KS_setup_abag
!***********************************************************************
    subroutine input_persistent_hydro(id,lun,done)
!
!  Read in the stored time of the next random phase calculation
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
      if (lroot) print*,'input_persistent_hydro: ',tphase_kinflow
!
    endsubroutine input_persistent_hydro
!***********************************************************************
    subroutine output_persistent_hydro(lun)
!
!  Writes out the time of the next random phase calculation
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
!  write details
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
!  reads and registers print parameters relevant for hydro part
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
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
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
!  write column where which hydro variable is stored
!
      if (lwr) then
        write(3,*) 'i_ekin=',idiag_ekin
        write(3,*) 'i_ekintot=',idiag_ekintot
        write(3,*) 'i_u2m=',idiag_u2m
        write(3,*) 'i_um2=',idiag_um2
        write(3,*) 'i_o2m=',idiag_o2m
        write(3,*) 'i_oum=',idiag_oum
        write(3,*) 'i_dtu=',idiag_dtu
        write(3,*) 'i_urms=',idiag_urms
        write(3,*) 'i_umax=',idiag_umax
        write(3,*) 'i_uzrms=',idiag_uzrms
        write(3,*) 'i_uzmax=',idiag_uzmax
        write(3,*) 'i_ux2m=',idiag_ux2m
        write(3,*) 'i_uy2m=',idiag_uy2m
        write(3,*) 'i_uz2m=',idiag_uz2m
        write(3,*) 'i_uxuym=',idiag_uxuym
        write(3,*) 'i_uxuzm=',idiag_uxuzm
        write(3,*) 'i_uyuzm=',idiag_uyuzm
        write(3,*) 'i_orms=',idiag_orms
        write(3,*) 'i_omax=',idiag_omax
        write(3,*) 'i_ruxm=',idiag_ruxm
        write(3,*) 'i_ruym=',idiag_ruym
        write(3,*) 'i_ruzm=',idiag_ruzm
        write(3,*) 'i_rumax=',idiag_rumax
        write(3,*) 'i_umx=',idiag_umx
        write(3,*) 'i_umy=',idiag_umy
        write(3,*) 'i_umz=',idiag_umz
        write(3,*) 'i_Marms=',idiag_Marms
        write(3,*) 'i_Mamax=',idiag_Mamax
        write(3,*) 'i_divu2m=',idiag_divu2m
        write(3,*) 'i_epsK=',idiag_epsK
        write(3,*) 'i_uxpt=',idiag_uxpt
        write(3,*) 'i_uypt=',idiag_uypt
        write(3,*) 'i_uzpt=',idiag_uzpt
        write(3,*) 'i_uxmz=',idiag_uxmz
        write(3,*) 'i_uymz=',idiag_uymz
        write(3,*) 'i_uzmz=',idiag_uzmz
        write(3,*) 'i_uxmxy=',idiag_uxmxy
        write(3,*) 'i_uymxy=',idiag_uymxy
        write(3,*) 'i_uzmxy=',idiag_uzmxy
        write(3,*) 'i_urmphi=',idiag_urmphi
        write(3,*) 'i_upmphi=',idiag_upmphi
        write(3,*) 'i_uzmphi=',idiag_uzmphi
        write(3,*) 'i_u2mphi=',idiag_u2mphi
        write(3,*) 'i_oumphi=',idiag_oumphi
        write(3,*) 'i_phase1=',idiag_phase1
        write(3,*) 'i_phase2=',idiag_phase2
        write(3,*) 'i_divum=',idiag_divum
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
!  dummy routine
!
!  19-jul-03/axel: adapted from hydro
!
    endsubroutine calc_mflow
!***********************************************************************
    subroutine remove_mean_momenta(f)
!
!  dummy routine
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
!  Deallocate the variables allocated in nohydro
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
      print*, 'Done.'
!
    endsubroutine hydro_clean_up
!*******************************************************************
    subroutine kinematic_random_phase
!
!  Get a random phase to be used for the whole kinematic velocity field.
!
!  16-feb-2010/dhruba: coded
!
      use General, only: random_number_wrapper
      real, dimension(3) :: fran
!
!  generate random numbers
!
      if (t>tsforce) then
        if (lrandom_location) then
          call random_number_wrapper(fran)
          location=fran*Lxyz+xyz0
        else
          location=location_fixed
        endif
!
!  writing down the location
!
        if (lroot .and. lwrite_random_location) then
          open(1,file=trim(datadir)//'/random_location.dat',status='unknown',position='append')
            write(1,'(4f14.7)') t, location
          close (1)
        endif
!
!  update next tsforce
!
        tsforce=t+dtforce
        if (ip<=6) print*,'kinematic_random_phase: location=',location
      endif
!
    endsubroutine kinematic_random_phase
!*******************************************************************
endmodule Hydro
