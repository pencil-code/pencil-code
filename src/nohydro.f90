! $Id$

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lhydro = .false.
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED oo(3); ou; uij(3,3); uu(3); u2; sij(3,3)
! PENCILS PROVIDED divu; uij5(3,3)
!
!***************************************************************

module Hydro

  use Cparam
  use Cdata, only: huge1
  use Messages
  use Sub, only: keep_compiler_quiet

  implicit none

  include 'record_types.h'
  include 'hydro.h'

  real, dimension (nz,3) :: uumz=0.
  real :: u_out_kep=0.0
  real :: tphase_kinflow=-1.,phase1=0., phase2=0.
  logical :: lpressuregradient_gas=.true.,lcalc_uumean=.false.

  real, allocatable, dimension (:,:) :: KS_k,KS_A,KS_B !or through whole field for each wavenumber?
  real, allocatable, dimension (:) :: KS_omega !or through whole field for each wavenumber?
  integer :: KS_modes = 3
  !namelist /hydro_init_pars/ dummyuu
  !namelist /hydro_run_pars/  dummyuu

  ! other variables (needs to be consistent with reset list below)
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
  integer :: idiag_ekintot=0, idiag_ekin=0

  contains

!***********************************************************************
    subroutine register_hydro()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  6-nov-01/wolf: coded
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
      use SharedVariables
      use Sub
!
      logical, save :: first=.true.
      integer :: ierr
!
      if (.not. first) call stop_it('register_hydro called twice')
      first = .false.
!
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
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
      use Cdata
      use FArrayManager
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      if (kinflow=='KS') then
!        call random_isotropic_KS_setup(-5./3.,1.,(nxgrid)/2.)
!
!  Use constant values for testing KS model code with 3
!  specific modes.
!
        call random_isotropic_KS_setup_test
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
      if (ip == 0) print*,f,lstarting  !(keep compiler quiet)
!
    endsubroutine initialize_hydro
!***********************************************************************
    subroutine read_hydro_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) call keep_compiler_quiet(iostat)
      call keep_compiler_quiet(unit)

    endsubroutine read_hydro_init_pars
!***********************************************************************
    subroutine write_hydro_init_pars(unit)
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_hydro_init_pars
!***********************************************************************
    subroutine read_hydro_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) call keep_compiler_quiet(iostat)
      call keep_compiler_quiet(unit)

    endsubroutine read_hydro_run_pars
!***********************************************************************
    subroutine write_hydro_run_pars(unit)
      integer, intent(in) :: unit

      call keep_compiler_quiet(unit)

    endsubroutine write_hydro_run_pars
!***********************************************************************
    subroutine init_uu(f,xx,yy,zz)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Hydro module, if there was one.
!
!   7-jun-02/axel: adapted from hydro
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(xx,yy,zz)

    endsubroutine init_uu
!***********************************************************************
    subroutine pencil_criteria_hydro()
!
!  All pencils that the Hydro module depends on are specified here.
!
!  20-11-04/anders: coded
!
      use Cdata
!
      if (kinflow/='') lpenc_requested(i_uu)=.true.
      if (idiag_urms/=0 .or. idiag_umax/=0 .or. idiag_u2m/=0 .or. &
          idiag_um2/=0) lpenc_diagnos(i_u2)=.true.
      if (idiag_oum/=0) lpenc_diagnos(i_ou)=.true.
!
    endsubroutine pencil_criteria_hydro
!***********************************************************************
    subroutine pencil_interdep_hydro(lpencil_in)
!
!  Interdependency among pencils from the Hydro module is specified here
!
!  20-11-04/anders: coded
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
      use Cdata
      use Sub, only: dot_mn,dot2_mn, sum_mn_name, max_mn_name, &
        integrate_mn_name, quintic_step, quintic_der_step
      use Magnetic, only: ABC_A,ABC_B,ABC_C,kx_aa,ky_aa,kz_aa
      use General
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      real, dimension(nx) :: kdotxwt,cos_kdotxwt,sin_kdotxwt
      real :: kkx_aa,kky_aa,kkz_aa, fac, fpara, dfpara, ecost, esint, epst
      integer :: modeN
      real :: sqrt2, sqrt21k1, eps1=1., WW=0.25, k21
!
      intent(in) :: f
      intent(inout) :: p
!
!  introduce new symbol for the first argument of kx_aa, etc, arrays
!
      kkx_aa=kx_aa(1)
      kky_aa=ky_aa(1)
      kkz_aa=kz_aa(1)
!
!  choose from a list of different flow profiles.
!  Begin with a constant flow in the x direction.
!
      if (kinflow=='const-x') then
        if (lpencil(i_uu)) then
          if (headtt) print*,'const-x'
          p%uu(:,1)=1.
          p%uu(:,2)=0.
          p%uu(:,3)=0.
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
          p%uu(:,1)=ABC_A*sin(kkz_aa*z(n))    +ABC_C*cos(kky_aa*y(m))
          p%uu(:,2)=ABC_B*sin(kkx_aa*x(l1:l2))+ABC_A*cos(kkz_aa*z(n))
          p%uu(:,3)=ABC_C*sin(kky_aa*y(m))    +ABC_B*cos(kkx_aa*x(l1:l2))
        endif
! divu
        if (lpencil(i_divu)) p%divu=0.
!
!  Gen-Roberts flow (negative helicity)
!
      elseif (kinflow=='roberts') then
! uu
        if (lpencil(i_uu)) then
          if (headtt) print*,'Glen Roberts flow; kx_aa,ky_aa=',kkx_aa,kky_aa
          eps1=1.-eps_kinflow
          p%uu(:,1)=+sin(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*eps1
          p%uu(:,2)=-cos(kkx_aa*x(l1:l2))*sin(kky_aa*y(m))*eps1
          p%uu(:,3)=+sin(kkx_aa*x(l1:l2))*sin(kky_aa*y(m))*sqrt(2.)
        endif
! divu
        if (lpencil(i_divu)) p%divu= (kkx_aa-kky_aa)*cos(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))
!
!  Glen-Roberts flow (positive helicity)
!
      elseif (kinflow=='poshel-roberts') then
        if (headtt) print*,'Pos Helicity Roberts flow; kx_aa,ky_aa=',kkx_aa,kky_aa
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
        p%uu(:,1)=-fac*cos(kkx_aa*x(l1:l2))*sin(kky_aa*y(m))*eps1
        p%uu(:,2)=+fac*sin(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*eps1
        p%uu(:,3)=+fac*cos(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*sqrt(2.)
        if (lpencil(i_divu)) p%divu=0.
!
!  Glen-Roberts flow (x-direction, positive helicity)
!  x -> y
!  y -> z
!  z -> x
!
      elseif (kinflow=='xdir-roberts') then
        if (headtt) print*,'x-dir Roberts flow; ky_aa,kz_aa=',kky_aa,kkz_aa
        fac=ampl_kinflow
        eps1=1.-eps_kinflow
        p%uu(:,2)=-fac*cos(kky_aa*y(m))*sin(kkz_aa*z(n))*eps1
        p%uu(:,3)=+fac*sin(kky_aa*y(m))*cos(kkz_aa*z(n))*eps1
        p%uu(:,1)=+fac*cos(kky_aa*y(m))*cos(kkz_aa*z(n))*sqrt(2.)
        if (lpencil(i_divu)) p%divu=0.
!
!  z-dependent Roberts flow (positive helicity)
!
      elseif (kinflow=='zdep-roberts') then
        if (headtt) print*,'z-dependent Roberts flow; kx,ky=',kkx_aa,kky_aa
        fpara=ampl_kinflow*(quintic_step(z(n),-1.+eps_kinflow,eps_kinflow) &
                           -quintic_step(z(n),+1.-eps_kinflow,eps_kinflow))
        dfpara=ampl_kinflow*(quintic_der_step(z(n),-1.+eps_kinflow,eps_kinflow)&
                            -quintic_der_step(z(n),+1.-eps_kinflow,eps_kinflow))
!
        sqrt2=sqrt(2.)
        sqrt21k1=1./(sqrt2*kx_aa(1))
!
        p%uu(:,1)=-cos(kkx_aa*x(l1:l2))*sin(kky_aa*y(m)) &
           -dfpara*sin(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*sqrt21k1
        p%uu(:,2)=+sin(kkx_aa*x(l1:l2))*cos(kky_aa*y(m)) &
           -dfpara*cos(kkx_aa*x(l1:l2))*sin(kky_aa*y(m))*sqrt21k1
        p%uu(:,3)=+fpara*cos(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*sqrt2
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with cosinusoidal helicity variation (version 1)
!
      elseif (kinflow=='IncohRoberts1') then
        if (headtt) print*,'Roberts flow with cosinusoidal helicity; kx_aa,ky_aa=',kkx_aa,kky_aa
        fac=ampl_kinflow
        eps1=(1.-eps_kinflow)*cos(omega_kinflow*t)
        p%uu(:,1)=-fac*cos(kkx_aa*x(l1:l2))*sin(kky_aa*y(m))*eps1
        p%uu(:,2)=+fac*sin(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*eps1
        p%uu(:,3)=+fac*cos(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*sqrt(2.)
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with cosinusoidal helicity variation (version 2)
!
      elseif (kinflow=='IncohRoberts2') then
        if (headtt) print*,'Roberts flow with cosinusoidal helicity; kx_aa,ky_aa=',kkx_aa,kky_aa
        fac=ampl_kinflow
        eps1=(1.-eps_kinflow)*cos(omega_kinflow*t)
        p%uu(:,1)=-fac*cos(kkx_aa*x(l1:l2))*sin(kky_aa*y(m))
        p%uu(:,2)=+fac*sin(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))
        p%uu(:,3)=+fac*cos(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*sqrt(2.)*eps1
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with helicity variation and shear (version 1)
!  Must use shear module and set eps_kinflow equal to shear
!
      elseif (kinflow=='ShearRoberts1') then
        if (headtt) print*,'Roberts flow with cosinusoidal helicity; kx_aa,ky_aa=',kkx_aa,kky_aa
        fac=ampl_kinflow
        kky_aa=1.
        kkx_aa=kky_aa*(mod(.5-eps_kinflow*t,1.)-.5)
if (ip.eq.11.and.m==4.and.n==4) write(21,*) t,kkx_aa
        eps1=cos(omega_kinflow*t)
        p%uu(:,1)=-fac*cos(kkx_aa*x(l1:l2))*sin(kky_aa*y(m))*eps1
        p%uu(:,2)=+fac*sin(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*eps1
        p%uu(:,3)=+fac*cos(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*sqrt(2.)
        if (lpencil(i_divu)) p%divu=0.
!
!  "Incoherent" Roberts flow with helicity variation and shear (version 2)
!  Must use shear module and set eps_kinflow equal to shear
!
      elseif (kinflow=='ShearRoberts1') then
        if (headtt) print*,'Roberts flow with cosinusoidal helicity; kx_aa,ky_aa=',kkx_aa,kky_aa
        fac=ampl_kinflow
        kky_aa=1.
        kkx_aa=kky_aa*(mod(.5-eps_kinflow*t,1.)-.5)
if (ip.eq.11.and.m==4.and.n==4) write(21,*) t,kkx_aa
        eps1=cos(omega_kinflow*t)
        p%uu(:,1)=-fac*cos(kkx_aa*x(l1:l2))*sin(kky_aa*y(m))
        p%uu(:,2)=+fac*sin(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))
        p%uu(:,3)=+fac*cos(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*sqrt(2.)*eps1
        if (lpencil(i_divu)) p%divu=0.
!
!  Taylor-Green flow
!
      elseif (kinflow=='TG') then
        if (headtt) print*,'Taylor-Green flow; kx_aa,ky_aa=',kkx_aa,kky_aa
        p%uu(:,1)=+2.*sin(kkx_aa*x(l1:l2))*cos(kky_aa*y(m))*cos(kkz_aa*z(n))
        p%uu(:,2)=-2.*cos(kkx_aa*x(l1:l2))*sin(kky_aa*y(m))*cos(kkz_aa*z(n))
        p%uu(:,3)=+0.
        if (lpencil(i_divu)) p%divu=0.
!
!  Galloway-Proctor flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kkx_aa=kky_aa
!
      elseif (kinflow=='Galloway-Proctor') then
        if (headtt) print*,'Galloway-Proctor flow; kx_aa,ky_aa=',kkx_aa,kky_aa
        fac=ampl_kinflow
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
        p%uu(:,1)=-fac*sin(kky_aa*y(m)    +esint)
        p%uu(:,2)=+fac*sin(kkx_aa*x(l1:l2)+ecost)
        p%uu(:,3)=-fac*(cos(kkx_aa*x(l1:l2)+ecost)+cos(kky_aa*y(m)+esint))
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kkx_aa*p%uu
!
!  Tilgner flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kkx_aa=kky_aa
!
      elseif (kinflow=='Tilgner') then
        if (headtt) print*,'Tilgner flow; kx_aa,ky_aa=',kkx_aa,kky_aa
        fac=ampl_kinflow*sqrt(2.)
        epst=eps_kinflow*t
kkx_aa=2.*pi
kky_aa=2.*pi
        p%uu(:,1)=-fac* sin(kky_aa*y(m)         )
        p%uu(:,2)=+fac* sin(kkx_aa*(x(l1:l2)+epst))
        p%uu(:,3)=-fac*(cos(kkx_aa*(x(l1:l2)+epst))+cos(kky_aa*y(m)))
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kkx_aa*p%uu
!
!  Tilgner flow, U=-z x grad(psi) - z k psi, where
!  psi = U0/kH * (cosX+cosY), so U = U0 * (-sinY, sinX, -cosX-cosY).
!  This makes sense only for kkx_aa=kky_aa
!  Here, W in Tilgner's equation is chosen to be 0.25.
!
      elseif (kinflow=='Tilgner-orig') then
        if (headtt) print*,'original Tilgner flow; kx_aa,ky_aa=',kkx_aa,kky_aa
        fac=ampl_kinflow
        epst=eps_kinflow*t
        kkx_aa=2.*pi
        kky_aa=2.*pi
        sqrt2=sqrt(2.)
        WW=0.25
        p%uu(:,1)=+fac*sqrt2*sin(kkx_aa*(x(l1:l2)-epst))*cos(kky_aa*y(m))
        p%uu(:,2)=-fac*sqrt2*cos(kkx_aa*(x(l1:l2)-epst))*sin(kky_aa*y(m))
        p%uu(:,3)=+fac*2.*WW*sin(kkx_aa*(x(l1:l2)-epst))*sin(kky_aa*y(m))
        if (lpencil(i_divu)) p%divu=0.
        if (lpencil(i_oo)) p%oo=-kkx_aa*p%uu
!
!  Galloway-Proctor flow with random temporal phase
!
      elseif (kinflow=='Galloway-Proctor-RandomTemporalPhase') then
        if (headtt) print*,'GP-RandomTemporalPhase; kx,ky=',kkx_aa,kky_aa
        fac=ampl_kinflow
        if (t.gt.tphase_kinflow) then
          call random_number_wrapper(fran1)
          tphase_kinflow=t+dtphase_kinflow
          phase1=pi*(2*fran1(1)-1.)
          phase2=pi*(2*fran1(2)-1.)
        endif
        ecost=eps_kinflow*cos(omega_kinflow*t+phase1)
        esint=eps_kinflow*sin(omega_kinflow*t+phase2)
        p%uu(:,1)=-fac*sin(kky_aa*y(m)    +esint)
        p%uu(:,2)=+fac*sin(kkx_aa*x(l1:l2)+ecost)
        p%uu(:,3)=-fac*(cos(kkx_aa*x(l1:l2)+ecost)+cos(kky_aa*y(m)+esint))
        if (lpencil(i_divu)) p%divu=0.
!
!  Galloway-Proctor flow with random phase
!
      elseif (kinflow=='Galloway-Proctor-RandomPhase') then
        if (headtt) print*,'Galloway-Proctor-RandomPhase; kx,ky=',kkx_aa,kky_aa
        fac=ampl_kinflow
        if (t.gt.tphase_kinflow) then
          call random_number_wrapper(fran1)
          tphase_kinflow=t+dtphase_kinflow
          phase1=eps_kinflow*pi*(2*fran1(1)-1.)
          phase2=eps_kinflow*pi*(2*fran1(2)-1.)
        endif
        p%uu(:,1)=-fac*sin(kky_aa*y(m)    +phase1)*kky_aa
        p%uu(:,2)=+fac*sin(kkx_aa*x(l1:l2)+phase2)*kkx_aa
        p%uu(:,3)=-fac*(cos(kkx_aa*x(l1:l2)+phase2)+cos(kky_aa*y(m)+phase1))
        if (lpencil(i_divu)) p%divu=0.
!
!  original Galloway-Proctor flow
!
      elseif (kinflow=='Galloway-Proctor-orig') then
        if (headtt) print*,'Galloway-Proctor-orig flow; kx_aa=',kkx_aa
        fac=sqrt(1.5)*ampl_kinflow
        ecost=eps_kinflow*cos(omega_kinflow*t)
        esint=eps_kinflow*sin(omega_kinflow*t)
        p%uu(:,1)=+fac*cos(kky_aa*y(m)    +esint)*kky_aa
        p%uu(:,2)=+fac*sin(kkx_aa*x(l1:l2)+ecost)*kkx_aa
        p%uu(:,3)=-fac*(cos(kkx_aa*x(l1:l2)+ecost)+sin(kky_aa*y(m)+esint))
        if (lpencil(i_divu)) p%divu=0.
!
!  Convection rolls
!  Stream function: psi_y = cos(kx*x) * cos(kz*z)
!
      elseif (kinflow=='rolls') then
        if (headtt) print*,'Convection rolls; kx_aa,kz_aa=',kkx_aa,kkz_aa
        p%uu(:,1)=kkz_aa*cos(kkx_aa*x(l1:l2))*sin(kkz_aa*z(n))
        p%uu(:,2)=+0.
        p%uu(:,3)=kkx_aa*sin(kkx_aa*x(l1:l2))*cos(kkz_aa*z(n))
! divu
        if (lpencil(i_divu)) p%divu=0.
!
!  Twist (Yousef & Brandenburg 2003)
!
      elseif (kinflow=='Twist') then
        if (headtt) print*,'Twist flow; eps_kinflow,kx=',eps_kinflow,kkx_aa
        p%uu(:,1)=0.
        p%uu(:,2)=eps_kinflow*z(n)*sin(kkx_aa*x(l1:l2))
        p%uu(:,3)=1.+cos(kkx_aa*x(l1:l2))
        if (lpencil(i_divu)) p%divu=0.
!
!  Shearing wave
!
      elseif (kinflow=='ShearingWave') then
        if (headtt) print*,'ShearingWave flow; eps_kinflow=',eps_kinflow
        kky_aa=1.
        kkx_aa=kky_aa*eps_kinflow*(10.-t)
        k21=1./(kkx_aa**2+kky_aa**2)
        p%uu(:,1)=-kky_aa*k21*cos(kkx_aa*x(l1:l2)+kky_aa*y(m))
        p%uu(:,2)=+kkx_aa*k21*cos(kkx_aa*x(l1:l2)+kky_aa*y(m))
        p%uu(:,3)=            cos(kkx_aa*x(l1:l2)+kky_aa*y(m))
        if (lpencil(i_divu)) p%divu=0.
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
      use Cdata
      use Sub
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
      call keep_compiler_quiet(f,df)
!
    endsubroutine duu_dt
!***********************************************************************
    subroutine time_integrals_hydro(f,p)
!
!   1-jul-08/axel: dummy
!
      use Cdata
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
    use Cdata, only: pi,Lxyz
    use Sub
    use General
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
!       if(modeN.eq.1)delk(modeN)=(k(modeN+1)-K(modeN))
!       if(modeN.eq.KS_modes)delk(modeN)=(k(modeN)-k(modeN-1))
!       if(modeN.gt.1.and.modeN.lt.KS_modes)delk(modeN)=(k(modeN+1)-k(modeN-2))/2.0
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
      if(modeN==1)&
              dk=kmin*(a-1.)/2.
      if(modeN.gt.1.and.modeN.lt.KS_modes) &
              dk=(a**(modeN-2.))*kmin*((a**2.) -1.)/2.
      if(modeN==KS_modes) &
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
    use Cdata, only: pi,Lxyz
    use Sub
    use General
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

!weezy       if(modeN==1)dk=kmin*(a-1.)/2.
!weezy       if(modeN.gt.1.and.modeN.lt.KS_modes)dk=(a**(modeN-2.))*kmin*((a**2.) -1.)/2.
!weezy       if(modeN==KS_modes)dk=(a**(KS_modes -2.))*kmin*(a -1.)/2.

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
      if((k_unit(2).eq.0).and.(k_unit(3).eq.0)) then
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
!  With the `weezey' stuff above commented out, dk is currently used, but
!  never set, so we better abort
!
      call error('random_isotropic_KS_setup', 'Using uninitialized dk')
      dk=0.                     ! to make compiler happy

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
    use Cdata, only: pi,Lxyz
    use Sub
    use General
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
    end do
!
    kmax=k(KS_modes)
    kmin=k(1)
!
    do modeN=1,KS_modes
      if(modeN==1) dk(modeN)=(k(modeN+1)-k(modeN))/2.
      if(modeN.gt.1.and.modeN.lt.KS_modes) &
                dk(modeN)=(k(modeN+1)-k(modeN-1))/2.
      if(modeN==KS_modes) dk(modeN)=(k(modeN)-k(modeN-1))/2.
    end do
!
    do modeN=1,KS_modes
       energy(modeN)=((k(modeN)**2 +1.)**(-11./6.))*(k(modeN)**2) &
                         *exp(-0.5*(k(modeN)/kmax)**2)
    end do
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
    end do
!
!   form RA = RA x k_unit and RB = RB x k_unit
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
   ! use Cdata, only: pi,Lxyz
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
   !  if((angle(1)-0.0 < epsilon(0.0)) .or. &
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

   !  !find the length of the current k_option vector
   !  mkunit(i)=dsqrt((k_option(1,i)**2)+(k_option(2,i)**2)+(k_option(3,i)**2))

   !  if(i==1.and.mkunit(i).gt.0.)then 
   !    k(:,num)=k_option(:,i)
   !    klengths(num)=mkunit(i)
   !  endif

   !  !now we check that the current length is unique (hasn't come before)
   !  if(i.gt.1.and.num.lt.KS_modes)then
   !    do s1=i-1,1,-1
   !      if(mkunit(i).gt.0.0D0.and.mkunit(i) /= mkunit(s1))then
   !        ne=.true.
   !      else
   !        ne=.false.
   !        exit
   !      endif
   !      if(s1==1.and.ne)then !i.e. if length of current k_option is new...... 
   !        num=num+1
   !        k(:,num)=k_option(:,i) !load current k_option into k that we keep
   !        klengths(num)=mkunit(i)  ! store the length also
   !      endif
   !    end do
   !   endif
   !   if(i==10000.and.num.lt.KS_modes)print*,"Haven't got",KS_modes,"modes!!!!"
   ! end do
   ! do i=1,KS_modes
   !    do s1=1,KS_modes
   !       if(kk(i)==klengths(s1))then
   !          orderK(:,i)=k(:,s1)
   !       endif
   !    end do
   ! end do
   ! k=orderK
   ! do i=1,KS_modes
   !   unit_k(:,i)=k(:,i)/kk(i)
   ! end do
   ! do i=1,N
   ! !now we find delk as defined in Malik & Vassilicos' paper
   !    if(i==1)delk(i)=(kk(i+1)-kk(i))/2.0D0 
   !    if(i==KS_modes)delk(i)=(kk(i)-kk(i-1))/2.0D0 
   !    if(i.gt.1.and.i.lt.KS_modes)delk(i)=(kk(i+1)-kk(i-1))/2.0D0  
   ! end do
   ! endsubroutine random_isotropic_KS_setup_abag
!***********************************************************************
    subroutine input_persistent_hydro(id,lun,done)
!
!  Read in the stored time of the next random phase calculation
!
!  12-apr-08/axel: adapted from input_persistent_forcing
!
      use Cdata, only: lroot
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
      use Cdata, only: lroot
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
!
    endsubroutine output_persistent_hydro
!***********************************************************************
    subroutine rprint_hydro(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   8-jun-02/axel: adapted from hydro
!
      use Cdata
      use Sub
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
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_hydro: run through parse list'
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
!  Write slices for animation of hydro variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices)
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

      call keep_compiler_quiet(f)

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
endmodule Hydro
