! $Id: density.f90,v 1.341 2007-08-19 17:33:33 wlyra Exp $

!  This module is used both for the initial condition and during run time.
!  It contains dlnrho_dt and init_lnrho, among other auxiliary routines.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: ldensity = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED lnrho,rho,rho1,glnrho,grho,uglnrho,ugrho
! PENCILS PROVIDED glnrho2,del2lnrho,del2rho
! PENCILS PROVIDED del6lnrho,del6rho,hlnrho,sglnrho,uij5glnrho
!
!***************************************************************

module Density

  use Cparam
  use Cdata
  use Messages
  use EquationOfState, only: cs0,cs20,lnrho0,rho0, &
                             gamma,gamma1,cs2top,cs2bot, &
                             mpoly,beta_glnrho_global,&
                             get_cp1,get_ptlaw

  use Special

  implicit none

  include 'density.h'

  real :: ampllnrho=0.,widthlnrho=.1,rho_left=1.,rho_right=1.
  real :: cdiffrho=0.,diffrho=0.,diffrho_hyper3=0.,diffrho_shock=0.
  real :: lnrho_const=0., rho_const=1.
  real :: amplrho=0, phase_lnrho=0.0
  real :: radius_lnrho=.5,kx_lnrho=1.,ky_lnrho=1.,kz_lnrho=1.
  real :: eps_planet=.5
  real :: q_ell=5., hh0=0.
  real :: xblob=0., yblob=0., zblob=0.
  real :: co1_ss=0.,co2_ss=0.,Sigma1=150.
  real :: lnrho_int=0.,lnrho_ext=0.,damplnrho_int=0.,damplnrho_ext=0.
  real :: wdamp=0.,plaw=0, density_floor=-1.0
  real, dimension(3) :: diffrho_hyper3_aniso=0.
  integer, parameter :: ndiff_max=4
  logical :: lmass_source=.false.,lcontinuity_gas=.true.
  logical :: lupw_lnrho=.false.,lupw_rho=.false.
  logical :: ldiff_normal=.false.,ldiff_hyper3=.false.,ldiff_shock=.false.
  logical :: ldiff_hyper3lnrho=.false.,ldiff_hyper3_aniso=.false.
  logical :: lfreeze_lnrhoint=.false.,lfreeze_lnrhoext=.false.
  logical :: lfreeze_lnrhosqu=.false.

  character (len=labellen), dimension(ninit) :: initlnrho='nothing'
  character (len=labellen) :: strati_type='lnrho_ss',initlnrho2='nothing'
  character (len=labellen), dimension(ndiff_max) :: idiff=''
  character (len=labellen) :: borderlnrho='nothing'
  character (len=5) :: iinit_str
  complex :: coeflnrho=0.

  integer :: iglobal_gg=0
  integer :: iglobal_rho=0

  namelist /density_init_pars/ &
       ampllnrho,initlnrho,initlnrho2,widthlnrho,    &
       rho_left,rho_right,lnrho_const,rho_const,cs2bot,cs2top, &
       radius_lnrho,eps_planet,xblob,yblob,zblob,    &
       b_ell,q_ell,hh0,rbound,                       &
       mpoly,strati_type,beta_glnrho_global,         &
       kx_lnrho,ky_lnrho,kz_lnrho,amplrho,phase_lnrho,coeflnrho, &
       co1_ss,co2_ss,Sigma1,idiff,ldensity_nolog,    &
       wdamp,plaw,lcontinuity_gas,density_floor

  namelist /density_run_pars/ &
       cdiffrho,diffrho,diffrho_hyper3,diffrho_shock,   &
       cs2bot,cs2top,lupw_lnrho,lupw_rho,idiff,lmass_source,     &
       lnrho_int,lnrho_ext,damplnrho_int,damplnrho_ext, &
       wdamp,lfreeze_lnrhoint,lfreeze_lnrhoext,         &
       lnrho_const,plaw,lcontinuity_gas,borderlnrho,    &
       diffrho_hyper3_aniso,lfreeze_lnrhosqu,density_floor

  ! diagnostic variables (need to be consistent with reset list below)
  integer :: idiag_rhom=0       ! DIAG_DOC: $\left<\varrho\right>$
                                ! DIAG_DOC:   \quad(mean density)
  integer :: idiag_rho2m=0      ! DIAG_DOC:
  integer :: idiag_lnrho2m=0    ! DIAG_DOC:
  integer :: idiag_rhomin=0     ! DIAG_DOC:
  integer :: idiag_rhomax=0     ! DIAG_DOC:
  integer :: idiag_uglnrhom=0   ! DIAG_DOC:
  integer :: idiag_lnrhomphi=0  ! DIAG_DOC:
  integer :: idiag_rhomphi=0    ! DIAG_DOC:
  integer :: idiag_dtd=0        ! DIAG_DOC:
  integer :: idiag_rhomz=0      ! DIAG_DOC:
  integer :: idiag_rhomy=0      ! DIAG_DOC:
  integer :: idiag_rhomx=0      ! DIAG_DOC:
  integer :: idiag_rhomxy=0     ! DIAG_DOC:
  integer :: idiag_rhomxz=0     ! DIAG_DOC:
  integer :: idiag_rhomr=0      ! DIAG_DOC:
  integer :: idiag_totmass=0    ! DIAG_DOC:

  contains

!***********************************************************************
    subroutine register_density()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrho; increase nvar accordingly.
!
!   4-jun-02/axel: adapted from hydro
!
      use Sub
      use FArrayManager
!
      logical, save :: first=.true.
!      integer, target :: tmp_ilnrho
!
      if (.not. first) call fatal_error('register_density','module registration called twice')
      first = .false.

      call farray_register_pde('lnrho',ilnrho)
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: density.f90,v 1.341 2007-08-19 17:33:33 wlyra Exp $")
!
    endsubroutine register_density
!***********************************************************************
    subroutine initialize_density(f,lstarting)
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  For compatibility with other applications, we keep the possibility
!  of giving diffrho units of dxmin*cs0, but cs0 is not well defined general
!
!  24-nov-02/tony: coded
!  31-aug-03/axel: normally, diffrho should be given in absolute units
!
!

      use CData, only: lfreeze_varext,lfreeze_varint,lreloading,ilnrho,lfreeze_varsquare
      use FArrayManager
      use EquationOfState, only: select_eos_variable
      use Gravity, only: lnumerical_equilibrium
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      integer :: i
      logical :: lnothing
!
!  initialize cs2cool to cs20
!  (currently disabled, because it causes problems with mdarf auto-test)
!     cs2cool=cs20
!
!
      if(diffrho==0.) then
!
!  Made to work by adding diffrho + cdiffrho to the rprint reset list.
!
        diffrho=cdiffrho*dxmin*cs0
      endif
!
!  Turn off continuity equation term for 0-D runs.
!
      if (nxgrid*nygrid*nzgrid==1) then
        lcontinuity_gas=.false.
        print*, 'initialize_density: 0-D run, turned off continity equation'
      endif
!
!  Initialize dust diffusion
!
      ldiff_normal=.false.
      ldiff_shock=.false.
      ldiff_hyper3=.false.
      ldiff_hyper3lnrho=.false.
      ldiff_hyper3_aniso=.false.
!
      lnothing=.false.
!
      do i=1,ndiff_max
        select case (idiff(i))
        case ('normal')
          if (lroot) print*,'diffusion: div(D*grad(rho))'
          ldiff_normal=.true.
        case ('hyper3')
          if (lroot) print*,'diffusion: (d^6/dx^6+d^6/dy^6+d^6/dz^6)rho'
          ldiff_hyper3=.true.
        case ('hyper3lnrho')
          if (lroot) print*,'diffusion: (d^6/dx^6+d^6/dy^6+d^6/dz^6)lnrho'
          ldiff_hyper3lnrho=.true.
       case ('hyper3_aniso')
          if (lroot) print*,'diffusion: (Dx*d^6/dx^6 + Dy*d^6/dy^6 + Dz*d^6/dz^6)rho'
          ldiff_hyper3_aniso=.true.
        case ('shock')
          if (lroot) print*,'diffusion: shock diffusion'
          ldiff_shock=.true.
        case ('')
          if (lroot .and. (.not. lnothing)) print*,'diffusion: nothing'
        case default
          write(unit=errormsg,fmt=*) 'initialize_density: ', &
              'No such value for idiff(',i,'): ', trim(idiff(i))
          call fatal_error('initialize_density',errormsg)
        endselect
        lnothing=.true.
      enddo
!
!  If we're timestepping, die or warn if the the diffusion coefficient that
!  corresponds to the chosen diffusion type is not set.
!
      if (lrun) then
        if (ldiff_normal.and.diffrho==0.0) &
            call warning('initialize_density', &
            'Diffusion coefficient diffrho is zero!')
        if ( (ldiff_hyper3 .or. ldiff_hyper3lnrho) &
            .and. diffrho_hyper3==0.0) &
            call fatal_error('initialize_density', &
            'Diffusion coefficient diffrho_hyper3 is zero!')
        if ( (ldiff_hyper3_aniso) .and.  &
             ((diffrho_hyper3_aniso(1)==0. .and. nxgrid/=1 ).or. &
              (diffrho_hyper3_aniso(2)==0. .and. nygrid/=1 ).or. &
              (diffrho_hyper3_aniso(3)==0. .and. nzgrid/=1 )) ) &
            call fatal_error('initialize_density', &
            'A diffusion coefficient of diffrho_hyper3 is zero!')
        if (ldiff_shock .and. diffrho_shock==0.0) &
            call fatal_error('initialize_density', &
            'diffusion coefficient diffrho_shock is zero!')
      endif
!
!  Hyperdiffusion only works with (not log) density. One must either use
!  ldensity_nolog=T or work with GLOBAL =   global_nolog_density.
!
      if ((ldiff_hyper3.or.ldiff_hyper3_aniso) .and. (.not. ldensity_nolog)) then
        if (lroot) print*,"initialize_density: Creating global array for rho to use hyperdiffusion"
        call farray_register_global('rho',iglobal_rho)
      endif
!
      if (lfreeze_lnrhoint) lfreeze_varint(ilnrho)    = .true.
      if (lfreeze_lnrhoext) lfreeze_varext(ilnrho)    = .true.
      if (lfreeze_lnrhosqu) lfreeze_varsquare(ilnrho) = .true.
!
! Tell the equation of state that we're here and what f variable we use
!
        if (ldensity_nolog) then
          call select_eos_variable('rho',ilnrho)
        else
          call select_eos_variable('lnrho',ilnrho)
        endif
!
        if (lnumerical_equilibrium) then
           if (lroot) print*,'initializing global gravity in density'
           call farray_register_global('gg',iglobal_gg,vector=3)
        endif
!
      if (NO_WARN) print*,f,lstarting  !(to keep compiler quiet)
!
    endsubroutine initialize_density
!***********************************************************************
    subroutine read_density_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=density_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=density_init_pars,ERR=99)
      endif


99    return
    endsubroutine read_density_init_pars
!***********************************************************************
    subroutine write_density_init_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=density_init_pars)

    endsubroutine write_density_init_pars
!***********************************************************************
    subroutine read_density_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat

      if (present(iostat)) then
        read(unit,NML=density_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=density_run_pars,ERR=99)
      endif

99    return
    endsubroutine read_density_run_pars
!***********************************************************************
    subroutine write_density_run_pars(unit)
      integer, intent(in) :: unit

      write(unit,NML=density_run_pars)

    endsubroutine write_density_run_pars
!***********************************************************************
    subroutine init_lnrho(f,xx,yy,zz)
!
!  initialise lnrho; called from start.f90
!
!  7-nov-01/wolf: coded
! 28-jun-02/axel: added isothermal
! 15-oct-03/dave: added spherical shell (kws)
!
      use General, only: chn,complex_phase
      use Gravity, only: zref,z1,z2,gravz,nu_epicycle,potential, &
                          lnumerical_equilibrium
      use Selfgravity,only: rhs_poisson_const
      use Initcond
      use Initcond_spec
      use IO
      use Mpicomm
      use Sub
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,tmp,pot,prof
      real :: lnrhoint,cs2int,pot0,lnrho_left,lnrho_right
      real :: pot_ext,lnrho_ext,cs2_ext,tmp1,k_j2
      real :: zbot,ztop,haut,TT
      real, dimension (nx) :: r_mn,lnrho,lnTT,ss
      logical :: lnothing
      complex :: omega_jeans

!
!  define bottom and top height
!
      zbot=xyz0(3)
      ztop=xyz0(3)+Lxyz(3)
!
!  Set default values for sound speed at top and bottom.
!  These may be updated in one of the following initialization routines.
!
      cs2top=cs20; cs2bot=cs20
!
!  different initializations of lnrho (called from start).
!  If initrho does't match, f=0 is assumed (default).
!
      lnrho0      = log(rho0)
      lnrho_left  = log(rho_left)
      lnrho_right = log(rho_right)

      lnothing=.true.

      do iinit=1,ninit

      if (initlnrho(iinit)/='nothing') then

      lnothing=.false.

      call chn(iinit,iinit_str)

      select case(initlnrho(iinit))

      case('zero', '0'); f(:,:,:,ilnrho)=0.
      case('const_lnrho'); f(:,:,:,ilnrho)=lnrho_const
      case('const_rho'); f(:,:,:,ilnrho)=log(rho_const)
      case('constant'); f(:,:,:,ilnrho)=log(rho_left)
      case('mode'); call modes(ampllnrho,coeflnrho,f,ilnrho,kx_lnrho,ky_lnrho,kz_lnrho,xx,yy,zz)
      case('blob'); call blob(ampllnrho,f,ilnrho,radius_lnrho,xblob,yblob,zblob)
      case('isothermal'); call isothermal_density(f)
      case('local-isothermal'); call local_isothermal_density(f)  
      case('stratification'); call stratification(f,strati_type)
      case('polytropic_simple'); call polytropic_simple(f)
      case('hydrostatic-z', '1'); print*, &
                              'init_lnrho: use polytropic_simple instead!'
      case('xjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'x')
      case('yjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'y')
      case('zjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'z')
      case('soundwave-x'); call soundwave(ampllnrho,f,ilnrho,kx=kx_lnrho)
      case('soundwave-y'); call soundwave(ampllnrho,f,ilnrho,ky=ky_lnrho)
      case('soundwave-z'); call soundwave(ampllnrho,f,ilnrho,kz=kz_lnrho)
      case('sinwave-phase'); call sinwave_phase(f,ilnrho,ampllnrho,kx_lnrho,ky_lnrho,kz_lnrho,phase_lnrho)
      case('coswave-phase'); call coswave_phase(f,ilnrho,ampllnrho,kx_lnrho,ky_lnrho,kz_lnrho,phase_lnrho)
      case('sinwave-x'); call sinwave(ampllnrho,f,ilnrho,kx=kx_lnrho)
      case('sinwave-y'); call sinwave(ampllnrho,f,ilnrho,ky=ky_lnrho)
      case('sinwave-z'); call sinwave(ampllnrho,f,ilnrho,kz=kz_lnrho)
      case('coswave-x'); call coswave(ampllnrho,f,ilnrho,kx=kx_lnrho)
      case('coswave-y'); call coswave(ampllnrho,f,ilnrho,ky=ky_lnrho)
      case('coswave-z'); call coswave(ampllnrho,f,ilnrho,kz=kz_lnrho)
      case('sinx_siny_sinz'); call sinx_siny_sinz(ampllnrho,f,ilnrho,kx_lnrho,ky_lnrho,kz_lnrho)
      case('corona'); call corona_init(f)
      case('gaussian3d'); call gaussian3d(ampllnrho,f,ilnrho,xx,yy,zz,radius_lnrho)
      case('gaussian-z')
        do n=n1,n2
          f(:,:,n,ilnrho) = f(:,:,n,ilnrho) + &
              alog(exp(f(:,:,n,ilnrho))+ampllnrho*(exp(-z(n)**2/(2*radius_lnrho**2))))
        enddo
      case('gaussian-noise')
        If (lnrho_left /= 0.) f(:,:,:,ilnrho)=lnrho_left
        call gaunoise(ampllnrho,f,ilnrho,ilnrho)
      case('gaussian-noise-x')
        !
        !  noise, but just x-dependent
        !
        call gaunoise(ampllnrho,f,ilnrho,ilnrho)
        f(:,:,:,ilnrho)=spread(spread(f(:,4,4,ilnrho),2,my),3,mz) !(watch 1-D)

      case('rho-jump-z', '2')
        !
        !  density jump (for shocks)
        !
        if (lroot) print*, &
               'init_lnrho: density jump; rho_left,right=',rho_left,rho_right
        if (lroot) print*, &
               'init_lnrho: density jump; widthlnrho=',widthlnrho
        prof=.5*(1.+tanh(zz/widthlnrho))
        f(:,:,:,ilnrho)=log(rho_left)+log(rho_left/rho_right)*prof

      case ('hydrostatic-z-2', '3')
        !
        !  hydrostatic density stratification for isentropic atmosphere
        !
        if (lgravz) then
          if (lroot) print*,'init_lnrho: vertical density stratification'
          !        f(:,:,:,ilnrho)=-zz
          ! isentropic case:
          !        zmax = -cs20/gamma1/gravz
          !        print*, 'zmax = ', zmax
          !        f(:,:,:,ilnrho) = 1./gamma1 * log(abs(1-zz/zmax))
          ! linear entropy gradient;
          f(:,:,:,ilnrho) = -grads0*zz &
                            + 1./gamma1*log( 1 + gamma1*gravz/grads0/cs20 &
                                                  *(1-exp(-grads0*zz)) )
        endif

      case ('hydrostatic-r')
        !
        !  hydrostatic radial density stratification for isentropic (or
        !  isothermal) sphere
        !
        if (lgravr) then
          if (lroot) print*, &
               'init_lnrho: radial density stratification (assumes s=const)'

          call potential(xx,yy,zz,POT=pot) ! gravity potential
          call potential(R=r_ref,POT=pot0)
          call output(trim(directory)//'/pot.dat',pot,1)
          !
          ! rho0, cs0, pot0 are the values at r=r_ref
          !
          if (gamma /= 1) then  ! isentropic
            f(:,:,:,ilnrho) = lnrho0 &
                              + log(1 - gamma1*(pot-pot0)/cs20) / gamma1
          else                  ! isothermal
            f(:,:,:,ilnrho) = lnrho0 - (pot-pot0)/cs20
          endif

          !
          ! the following sets gravity gg in order to achieve numerical
          ! exact equilibrium at t=0
          !

          if (lnumerical_equilibrium) call numerical_equilibrium(f)

        endif

      case ('sph_isoth')
        if (lgravr) then
          if (lroot) print*, 'init_lnrho: isothermal sphere'
          haut=cs20/gamma
          TT=cs20/gamma1
          lnTT=spread(alog(TT),1,nx)
          do n=n1,n2
          do m=m1,m2
            r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
            f(l1:l2,m,n,ilnrho)=lnrho0-r_mn/haut
            lnrho=f(l1:l2,m,n,ilnrho)
            call eoscalc(ilnrho_lnTT,lnrho,lnTT,ss=ss)
            f(l1:l2,m,n,iss)=ss
          enddo
          enddo
        endif

      case ('isentropic-star')
        !
        !  isentropic/isothermal hydrostatic sphere"
        !    ss  = 0       for r<R,
        !    cs2 = const   for r>R
        !
        !  Only makes sense if both initlnrho=initss='isentropic-star'
        !
        if (lgravr) then
          if (lroot) print*, &
               'init_lnrho: isentropic star with isothermal atmosphere'
          call potential(xx,yy,zz,POT=pot,POT0=pot0) ! gravity potential
          call output(trim(directory)//'/pot.dat',pot,1)
          !
          ! rho0, cs0, pot0 are the values in the centre
          !
          if (gamma /= 1) then
            ! Note:
            ! (a) `where' is expensive, but this is only done at
            !     initialization.
            ! (b) Comparing pot with pot_ext instead of r with r_ext will
            !     only work if grav_r<=0 everywhere -- but that seems
            !     reasonable.
            call potential(R=r_ext,POT=pot_ext) ! get pot_ext=pot(r_ext)
            ! Do consistency check before taking the log() of a potentially
            ! negative number
            tmp1 = 1 - gamma1*(pot_ext-pot0)/cs20
            if (tmp1 <= 0.) then
              if (lroot) then
                print*, 'BAD IDEA: Trying to calculate log(', tmp1, ')'
                print*, '  for r_ext -- need to increase cs20?'
              endif
              call error('init_lnrho', 'Imaginary density values')
            endif
            lnrho_ext = lnrho0 + log(tmp1) / gamma1
            cs2_ext   = cs20*tmp1
            ! Adjust for given cs2cool (if given) or set cs2cool (otherwise)
            if (cs2cool/=0) then
              lnrho_ext = lnrho_ext - log(cs2cool/cs2_ext)
            else
              cs2cool   = cs2_ext
            endif
            ! Add temperature and entropy jump (such that pressure
            ! remains continuous) if cs2cool was specified in start.in:
            ! where (sqrt(xx**2+yy**2+zz**2) <= r_ext) ! isentropic for r<r_ext
            where (pot <= pot_ext) ! isentropic for r<r_ext
              f(:,:,:,ilnrho) = lnrho0 &
                                + log(1 - gamma1*(pot-pot0)/cs20) / gamma1
            elsewhere           ! isothermal for r>r_ext
              f(:,:,:,ilnrho) = lnrho_ext - gamma*(pot-pot_ext)/cs2cool
            endwhere
          else                  ! gamma=1 --> simply isothermal (I guess [wd])
            f(:,:,:,ilnrho) = lnrho0 - (pot-pot0)/cs20
          endif
        endif

      case ('piecew-poly', '4')
        !
        !  piecewise polytropic for stellar convection models
        !
        if (lroot) print*, &
           'init_lnrho: piecewise polytropic vertical stratification (lnrho)'
        ! top region

        cs2int = cs0**2
        lnrhoint = lnrho0
        f(:,:,:,ilnrho) = lnrho0 ! just in case
        call polytropic_lnrho_z(f,mpoly2,zz,tmp,zref,z2,ztop+Lz, &
                            isothtop,cs2int,lnrhoint)
          ! unstable layer
        call polytropic_lnrho_z(f,mpoly0,zz,tmp,z2,z1,z2,0,cs2int,lnrhoint)
          ! stable layer
        call polytropic_lnrho_z(f,mpoly1,zz,tmp,z1,z0,z1,0,cs2int,lnrhoint)
        !
        !  calculate cs2bot and cs2top for run.x (boundary conditions)
        !
        cs2bot = cs2int + gamma/gamma1*gravz/(mpoly2+1)*(zbot-z0  )
        if (isothtop /= 0) then
          cs2top = cs20
        else
          cs2top = cs20 + gamma/gamma1*gravz/(mpoly0+1)*(ztop-zref)
        endif

      case ('piecew-disc', '41')
        !
        !  piecewise polytropic for accretion discs
        !
        if (lroot) print*, &
             'init_lnrho: piecewise polytropic disc stratification (lnrho)'
        ! bottom region

        cs2int = cs0**2
        lnrhoint = lnrho0
        f(:,:,:,ilnrho) = lnrho0 ! just in case
        call polytropic_lnrho_disc(f,mpoly1,zz,tmp,zref,z1,z1, &
                            0,cs2int,lnrhoint)
          ! unstable layer
        call polytropic_lnrho_disc(f,mpoly0,zz,tmp,z1,z2,z2,0,cs2int,lnrhoint)
          ! stable layer (top)
        call polytropic_lnrho_disc(f,mpoly2,zz,tmp,z2,ztop,ztop, &
                                   isothtop,cs2int,lnrhoint)
        !
        !  calculate cs2bot and cs2top for run.x (boundary conditions)
        !
        !cs2bot = cs2int + gamma/gamma1*gravz*nu_epicycle**2/(mpoly2+1)* &
        !         (zbot**2-z0**2)/2.
        cs2bot = cs20
        if (isothtop /= 0) then
          cs2top = cs20
        else
          cs2top = cs20 + gamma/gamma1*gravz*nu_epicycle**2/(mpoly0+1)* &
                  (ztop**2-zref**2)/2.
        endif

      case ('polytropic', '5')
        !
        !  polytropic stratification
        !  cs0, rho0 and ss0=0 refer to height z=zref
        !
        if (lroot) print*, &
                  'init_lnrho: polytropic vertical stratification (lnrho)'
        !
        cs2int = cs20
        lnrhoint = lnrho0
        f(:,:,:,ilnrho) = lnrho0 ! just in case
        ! only one layer
        call polytropic_lnrho_z(f,mpoly0,zz,tmp,zref,z0,z0+2*Lz, &
             0,cs2int,lnrhoint)
        !
        !  calculate cs2bot and cs2top for run.x (boundary conditions)
        !
        cs2bot = cs20 + gamma*gravz/(mpoly0+1)*(zbot-zref)
        cs2top = cs20 + gamma*gravz/(mpoly0+1)*(ztop-zref)

      case('sound-wave', '11')
        !
        !  sound wave (should be consistent with hydro module)
        !
        if (lroot) print*,'init_lnrho: x-wave in lnrho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=lnrho_const+ampllnrho*sin(kx_lnrho*xx)

      case('sound-wave-exp')
        !
        !  sound wave (should be consistent with hydro module)
        !
        if (lroot) print*,'init_lnrho: x-wave in rho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=log(rho_const+amplrho*sin(kx_lnrho*xx))

      case('sound-wave2')
        !
        !  sound wave (should be consistent with hydro module)
        !
        if (lroot) print*,'init_lnrho: x-wave in lnrho; ampllnrho=',ampllnrho
        f(:,:,:,ilnrho)=lnrho_const+ampllnrho*cos(kx_lnrho*xx)

      case('shock-tube', '13')
        !
        !  shock tube test (should be consistent with hydro module)
        !
        call information('init_lnrho','polytopic standing shock')
        prof=.5*(1.+tanh(xx/widthlnrho))
        f(:,:,:,ilnrho)=log(rho_left)+(log(rho_right)-log(rho_left))*prof

      case('sin-xy')
        !
        !  sin profile in x and y
        !
        call information('init_lnrho','lnrho=sin(x)*sin(y)')
        f(:,:,:,ilnrho) = &
             log(rho0) + ampllnrho*sin(kx_lnrho*xx)*sin(ky_lnrho*yy)

      case('sin-xy-rho')
        !
        !  sin profile in x and y, but in rho, not ln(rho)
        !
        call information('init_lnrho','rho=sin(x)*sin(y)')
        f(:,:,:,ilnrho) = &
             log(rho0*(1+ampllnrho*sin(kx_lnrho*xx)*sin(ky_lnrho*yy)))

      case('linear')
        !
        !  linear profile in kk.xxx
        !
        call information('init_lnrho','linear profile')
        f(:,:,:,ilnrho) = log(rho0) &
             + ampllnrho*(kx_lnrho*xx+ky_lnrho*yy+kz_lnrho*zz) &
               / sqrt(kx_lnrho**2+ky_lnrho**2+kz_lnrho**2)

      case('planet')
        !
        !  planet solution of Goodman, Narayan & Goldreich (1987)
        !  (Simple 3-D)
        !
        call planet(rbound,f,xx,yy,zz,eps_planet,radius_lnrho, &
            gamma,cs20,rho0,widthlnrho,hh0)

      case('planet_hc')
        !
        !  planet solution of Goodman, Narayan & Goldreich (1987)
        !  (3-D with hot corona)
        !
        call planet_hc(amplrho,f,xx,yy,zz,eps_planet,radius_lnrho, &
            gamma,cs20,rho0,widthlnrho)

      case('Ferriere'); call information('init_lnrho','Ferriere set in entropy')

      case('geo-kws')
      !
      ! radial hydrostatic profile in shell region only
      !
        call information('init_lnrho','kws hydrostatic in spherical shell region')
        call shell_lnrho(f)

      case('geo-kws-constant-T','geo-benchmark')
      !
      ! radial hydrostatic profile throughout box, which is consistent
      ! with constant temperature in exterior regions, and gives continuous
      ! density at shell boundaries
      !
        call information('init_lnrho','kws hydrostatic in spherical shell and exterior')
        call shell_lnrho(f)

     case ('step_xz')
        call fatal_error('init_lnrho','neutron_star initial condition is now in the special/neutron_star.f90 code')

     case('jeans-wave-x')
! soundwave + self gravity
        omega_jeans = sqrt(cmplx(cs20*kx_lnrho**2 - rhs_poisson_const*rho0,0.))/(rho0*kx_lnrho)
        print*,'Re(omega_jeans), Im(omega_jeans), Abs(omega_jeans)',&
          real(omega_jeans),aimag(omega_jeans),abs(omega_jeans)

        f(:,:,:,ilnrho) = lnrho_const + &
          ampllnrho*sin(kx_lnrho*xx+phase_lnrho)
        f(:,:,:,iux) = f(:,:,:,iux) &
             + abs(omega_jeans*ampllnrho) * &
             sin(kx_lnrho*xx+phase_lnrho+complex_phase(omega_jeans*ampllnrho))

     case('jeans-wave-oblique')
! soundwave + self gravity
        k_j2 = kx_lnrho**2 + ky_lnrho**2 + kz_lnrho**2
        omega_jeans = sqrt(cmplx(cs20*k_j2 - rhs_poisson_const*rho0,0.))/ &
            (rho0*sqrt(k_j2))
        print*,'Re(omega_jeans), Im(omega_jeans), Abs(omega_jeans)',&
          real(omega_jeans),aimag(omega_jeans),abs(omega_jeans)

        f(:,:,:,ilnrho) = lnrho_const + &
          ampllnrho*sin(kx_lnrho*xx + ky_lnrho*yy + kz_lnrho*zz)
        if (kx_lnrho .ne. 0) &
             f(:,:,:,iux) = f(:,:,:,iux) &
             + abs(omega_jeans*ampllnrho) * &
             sin(kx_lnrho*xx+complex_phase(omega_jeans*ampllnrho))
        if (ky_lnrho .ne. 0) &
             f(:,:,:,iuy) = f(:,:,:,iuy) &
             + abs(omega_jeans*ampllnrho) * &
             sin(ky_lnrho*yy+complex_phase(omega_jeans*ampllnrho))
        if (kz_lnrho .ne. 0) &
             f(:,:,:,iuz) = f(:,:,:,iuz) &
             + abs(omega_jeans*ampllnrho) * &
             sin(kz_lnrho*zz+complex_phase(omega_jeans*ampllnrho))

     case('toomre-wave-x')
! soundwave + self gravity + (differential) rotation
        omega_jeans = sqrt(cmplx(cs20*kx_lnrho**2 + Omega**2 - rhs_poisson_const*rho0,0.))/(rho0*kx_lnrho)

        print*,'Re(omega_jeans), Im(omega_jeans), Abs(omega_jeans)',&
          real(omega_jeans),aimag(omega_jeans),abs(omega_jeans)

        f(:,:,:,ilnrho) = lnrho_const + &
          ampllnrho*sin(kx_lnrho*xx)
        f(:,:,:,iux) = f(:,:,:,iux) &
             + abs(omega_jeans*ampllnrho) * &
             sin(kx_lnrho*xx+complex_phase(omega_jeans*ampllnrho))
        f(:,:,:,iuy) = f(:,:,:,iuy) &
             + abs(ampllnrho*cmplx(0,-0.5*Omega/(kx_lnrho*rho0))) * &
             sin(kx_lnrho*xx+complex_phase(ampllnrho*cmplx(0,-0.5*Omega/(kx_lnrho*rho0))))

      case('cylind-poly')
        call information('init_lnrho',' cylind-poly')
        call cylind_poly(f)

      case('compressive-shwave')
        ! should be consistent with density 
        f(:,:,:,ilnrho) = log(rho_const + f(:,:,:,ilnrho))
      case default
        !
        !  Catch unknown values
        !
        write(unit=errormsg,fmt=*) 'No such value for initlnrho(' &
                          //trim(iinit_str)//'): ',trim(initlnrho(iinit))
        call fatal_error('init_lnrho',errormsg)

      endselect

      if (lroot) print*,'init_lnrho: initlnrho(' &
                        //trim(iinit_str)//') = ',trim(initlnrho(iinit))

      endif

      enddo

      if (lnothing.and.lroot) print*,'init_lnrho: nothing'
!
!  check that cs2bot,cs2top are ok
!  for runs with ionization or fixed ionization, don't print them
!
      if (leos_ionization .or. leos_fixed_ionization) then
        cs2top=impossible
        cs2bot=impossible
      else
        if(lroot) print*,'init_lnrho: cs2bot,cs2top=',cs2bot,cs2top
      endif
!
!  Add some structures to the lnrho initialized above
!  Tobi: This is only kept for backwards compatibility
!
      select case(initlnrho2)

      case('nothing')

        if (lroot.and.ip<=5) print*,"init_lnrho: initlnrho2='nothing'"

      case('addblob')

        if (lroot) print*,'init_lnrho: add blob'
        if (lroot) print*,'initlnrho2 is deprecated. '// &
                          "use initlnrho='initcond1','initcond2' instead"
        f(:,:,:,ilnrho)=f(:,:,:,ilnrho) &
          +ampllnrho*exp(-(xx**2+yy**2+zz**2)/radius_lnrho**2)





      case default

        write(unit=errormsg,fmt=*) 'init_lnrho: No such value for initlnrho2: ', &
                          trim(initlnrho2)
        call fatal_error('init_lnrho',errormsg)

      endselect
!
!  If unlogarithmic density considered, take exp of lnrho resulting from
!  initlnrho
!
      if (ldensity_nolog) f(:,:,:,ilnrho)=exp(f(:,:,:,ilnrho))
!
!  sanity check
!
      if (notanumber(f(l1:l2,m1:m2,n1:n2,ilnrho))) then
        call error('init_lnrho', 'Imaginary density values')
      endif
!
    endsubroutine init_lnrho
!***********************************************************************
    subroutine polytropic_lnrho_z( &
         f,mpoly,zz,tmp,zint,zbot,zblend,isoth,cs2int,lnrhoint)
!
!  Implement a polytropic profile in ss above zbot. If this routine is
!  called several times (for a piecewise polytropic atmosphere), on needs
!  to call it from top to bottom.
!
!  zint    -- height of (previous) interface, where cs2int and lnrhoint
!             are set
!  zbot    -- z at bottom of layer
!  zblend  -- smoothly blend (with width whcond) previous ss (for z>zblend)
!             with new profile (for z<zblend)
!  isoth   -- flag for isothermal stratification;
!  lnrhoin -- value of lnrho at the interface, i.e. at the zint on entry,
!             at the zbot on exit
!  cs2int  -- same for cs2
!
      use Sub, only: step
      use Gravity, only: gravz
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: tmp,p,zz
      real, dimension (mz) :: stp
      real :: mpoly,zint,zbot,zblend,beta1,cs2int,lnrhoint
      integer :: isoth
!
      intent(in)    :: mpoly,zz,zint,zbot,zblend,isoth
      intent(out)   :: f,tmp
      intent(inout) :: cs2int,lnrhoint
!
      ! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
      if (isoth /= 0) then ! isothermal layer
        beta1 = 0.
        tmp = gamma*gravz/cs2int*(zz-zint)
        lnrhoint = lnrhoint + gamma*gravz/cs2int*(zbot-zint)
      else
        beta1 = gamma*gravz/(mpoly+1)
        tmp = 1 + beta1*(zz-zint)/cs2int
        ! Abort if args of log() are negative
        if (any((tmp <= 0.) .and. (zz <= zblend))) then
          call fatal_error('polytropic_lnrho_z', &
              'Imaginary density values -- your z_inf is too low.')
        endif
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = lnrhoint + mpoly*log(tmp)
        lnrhoint = lnrhoint + mpoly*log(1 + beta1*(zbot-zint)/cs2int)
      endif
      cs2int = cs2int + beta1*(zbot-zint) ! cs2 at layer interface (bottom)
      !
      ! smoothly blend the old value (above zblend) and the new one (below
      ! zblend) for the two regions:
      !
      stp = step(z,zblend,widthlnrho)
      p = spread(spread(stp,1,mx),2,my)
      f(:,:,:,ilnrho) = p*f(:,:,:,ilnrho) + (1-p)*tmp
!
    endsubroutine polytropic_lnrho_z
!***********************************************************************
    subroutine polytropic_lnrho_disc( &
         f,mpoly,zz,tmp,zint,zbot,zblend,isoth,cs2int,lnrhoint)
!
!  Implement a polytropic profile in a disc. If this routine is
!  called several times (for a piecewise polytropic atmosphere), on needs
!  to call it from bottom (middle of disc) upwards.
!
!  zint    -- height of (previous) interface, where cs2int and lnrhoint
!             are set
!  zbot    -- z at top of layer (name analogous with polytropic_lnrho_z)
!  zblend  -- smoothly blend (with width whcond) previous ss (for z>zblend)
!             with new profile (for z<zblend)
!  isoth   -- flag for isothermal stratification;
!  lnrhoint -- value of lnrho at the interface, i.e. at the zint on entry,
!             at the zbot on exit
!  cs2int  -- same for cs2
!
!  24-jun-03/ulf:  coded
!
      use Sub, only: step
      use Gravity, only: gravz, nu_epicycle
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz) :: tmp,p,zz
      real, dimension (mz) :: stp
      real :: mpoly,zint,zbot,zblend,beta1,cs2int,lnrhoint,nu_epicycle2
      integer :: isoth
!
      ! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
      nu_epicycle2 = nu_epicycle**2
      if (isoth /= 0) then ! isothermal layer
        beta1 = 0.
        tmp = gamma*gravz*nu_epicycle2/cs2int*(zz**2-zint**2)/2.
        lnrhoint = lnrhoint + gamma*gravz*nu_epicycle2/cs2int* &
                   (zbot**2-zint**2)/2.
      else
        beta1 = gamma*gravz*nu_epicycle2/(mpoly+1)
        tmp = 1 + beta1*(zz**2-zint**2)/cs2int/2.
        ! Abort if args of log() are negative
        if (any((tmp <= 0.) .and. (zz <= zblend))) then
          call fatal_error('polytropic_lnrho_disc', &
              'Imaginary density values -- your z_inf is too low.')
        endif
        tmp = max(tmp,epsi)  ! ensure arg to log is positive
        tmp = lnrhoint + mpoly*log(tmp)
        lnrhoint = lnrhoint + mpoly*log(1 + beta1*(zbot**2-zint**2)/cs2int/2.)
      endif
      cs2int = cs2int + beta1*(zbot**2-zint**2)/2. ! cs2 at layer interface (bottom)
      !
      ! smoothly blend the old value (above zblend) and the new one (below
      ! zblend) for the two regions:
      !
      stp = step(z,zblend,widthlnrho)
      p = spread(spread(stp,1,mx),2,my)
      f(:,:,:,ilnrho) = p*f(:,:,:,ilnrho) + (1-p)*tmp
!
    endsubroutine polytropic_lnrho_disc
!***********************************************************************
    subroutine shell_lnrho(f)
!
!  Initialize density based on specified radial profile in
!  a spherical shell
!
!  22-oct-03/dave -- coded
!
      use Gravity, only: g0,potential
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: pot, r_mn
      real :: beta1,lnrho_int,lnrho_ext,pot_int,pot_ext
!
      beta1=g0/(mpoly+1)*gamma/gamma1  ! gamma1/gamma=R_{*} (for cp=1)
!
!     densities at shell boundaries
      lnrho_int=lnrho0+mpoly*log(1+beta1*(r_ext/r_int-1.))
      lnrho_ext=lnrho0
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
!
        r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
!
        ! in the fluid shell
        where (r_mn < r_ext .AND. r_mn > r_int) f(l1:l2,m,n,ilnrho)=lnrho0+mpoly*log(1+beta1*(r_ext/r_mn-1.))
        ! outside the fluid shell
        if (initlnrho(1)=='geo-kws') then
          where (r_mn >= r_ext) f(l1:l2,m,n,ilnrho)=lnrho_ext
          where (r_mn <= r_int) f(l1:l2,m,n,ilnrho)=lnrho_int
        elseif (initlnrho(1)=='geo-kws-constant-T'.or.initlnrho(1)=='geo-benchmark') then
          call potential(R=r_int,POT=pot_int)
          call potential(R=r_ext,POT=pot_ext)
          call potential(RMN=r_mn,POT=pot)
          ! gamma/gamma1=1/R_{*} (for cp=1)
          where (r_mn >= r_ext) f(l1:l2,m,n,ilnrho)=lnrho_ext+(pot_ext-pot)*exp(-lnrho_ext/mpoly)*gamma/gamma1
          where (r_mn <= r_int) f(l1:l2,m,n,ilnrho)=lnrho_int+(pot_int-pot)*exp(-lnrho_int/mpoly)*gamma/gamma1
        endif
      enddo
!
    endsubroutine shell_lnrho
!***********************************************************************
    subroutine numerical_equilibrium(f)
!
!  sets gravity gg in order to achieve an numerical exact equilbrium
!  at t=0. This is only valid for the polytropic case, i.e.
!
!    (1/rho) grad(P) = cs20 (rho/rho0)^(gamma-2) grad(rho)
!
      use Sub, only: grad
      use IO

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: lnrho,cs2
      real, dimension (nx,3) :: glnrho
      real, dimension (nx,3) :: gg_mn
      integer :: i,j

      do m=m1,m2
      do n=n1,n2

        lnrho=f(l1:l2,m,n,ilnrho)
        cs2=cs20*exp(gamma1*(lnrho-lnrho0))
        call grad(f,ilnrho,glnrho)
        do j=1,3
          gg_mn(:,j)=cs2*glnrho(:,j)
        enddo
        f(l1:l2,m,n,iglobal_gg:iglobal_gg+2)=gg_mn

      enddo
      enddo

    endsubroutine numerical_equilibrium
!***********************************************************************
    subroutine pencil_criteria_density()
!
!  All pencils that the Density module depends on are specified here.
!
!  19-11-04/anders: coded
!
      use Cdata
!
      if (ldensity_nolog) lpenc_requested(i_rho)=.true.
      if (lcontinuity_gas) then
        lpenc_requested(i_divu)=.true.
        if (ldensity_nolog) then
          lpenc_requested(i_ugrho)=.true.
        else
          lpenc_requested(i_uglnrho)=.true.
        endif
      endif
      if (ldiff_shock) then
        lpenc_requested(i_shock)=.true.
        lpenc_requested(i_gshock)=.true.
        if (ldensity_nolog) then
          lpenc_requested(i_grho)=.true.
          lpenc_requested(i_del2rho)=.true.
        else
          lpenc_requested(i_glnrho)=.true.
          lpenc_requested(i_glnrho2)=.true.
          lpenc_requested(i_del2lnrho)=.true.
        endif
      endif
      if (ldiff_normal) then
        if (ldensity_nolog) then
          lpenc_requested(i_del2rho)=.true.
        else
          lpenc_requested(i_glnrho2)=.true.
          lpenc_requested(i_del2lnrho)=.true.
        endif
      endif
      if (ldiff_hyper3.or.ldiff_hyper3_aniso) lpenc_requested(i_del6rho)=.true.
      if (ldiff_hyper3.and..not.ldensity_nolog) lpenc_requested(i_rho)=.true.
      if (ldiff_hyper3lnrho) lpenc_requested(i_del6lnrho)=.true.
!
      if (lmass_source) lpenc_requested(i_rcyl_mn)=.true.
      if (borderlnrho=='stratification') then 
        lpenc_requested(i_cs2)=.true.
        lpenc_requested(i_r_mn)=.true.
      endif
!
      lpenc_diagnos2d(i_lnrho)=.true.
      lpenc_diagnos2d(i_rho)=.true.
!
      if (idiag_rhom/=0 .or. idiag_rhomz/=0 .or. idiag_rhomy/=0 .or. &
           idiag_rhomx/=0 .or. idiag_rho2m/=0 .or. idiag_rhomin/=0 .or. &
           idiag_rhomax/=0 .or. idiag_rhomxy/=0 .or. idiag_rhomxz/=0 .or. &
           idiag_totmass/=0) &
           lpenc_diagnos(i_rho)=.true.
      if (idiag_lnrho2m/=0) lpenc_diagnos(i_lnrho)=.true.
      if (idiag_uglnrhom/=0) lpenc_diagnos(i_uglnrho)=.true.
!
    endsubroutine pencil_criteria_density
!***********************************************************************
    subroutine pencil_interdep_density(lpencil_in)
!
!  Interdependency among pencils from the Density module is specified here.
!
!  19-11-04/anders: coded
!
      logical, dimension(npencils) :: lpencil_in
!
      if (ldensity_nolog) then
        if (lpencil_in(i_rho1)) lpencil_in(i_rho)=.true.
      else
        if (lpencil_in(i_rho)) lpencil_in(i_rho1)=.true.
      endif
      if (lpencil_in(i_uglnrho)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_ugrho)) then
        lpencil_in(i_uu)=.true.
        lpencil_in(i_grho)=.true.
      endif
      if (lpencil_in(i_glnrho2)) lpencil_in(i_glnrho)=.true.
      if (lpencil_in(i_sglnrho)) then
        lpencil_in(i_sij)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
      if (lpencil_in(i_uij5glnrho)) then
        lpencil_in(i_uij5)=.true.
        lpencil_in(i_glnrho)=.true.
      endif
!  The pencils glnrho and grho come in a bundle.
      if (lpencil_in(i_glnrho) .and. lpencil_in(i_grho)) then
        if (ldensity_nolog) then
          lpencil_in(i_grho)=.false.
        else
          lpencil_in(i_glnrho)=.false.
        endif
      endif
!
    endsubroutine pencil_interdep_density
!***********************************************************************
    subroutine calc_pencils_density(f,p)
!
!  Calculate Density pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  19-11-04/anders: coded
!
      use Sub, only: grad,dot,dot2,u_dot_grad,del2,del6,multmv,g2ij
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      integer :: i, mm, nn
!
      intent(inout) :: f,p
! lnrho
      if (lpencil(i_lnrho)) then
        if (ldensity_nolog) then
          p%lnrho=log(f(l1:l2,m,n,ilnrho))
        else
          p%lnrho=f(l1:l2,m,n,ilnrho)
        endif
      endif
! rho1 and rho
      if (ldensity_nolog) then
        if (lpencil(i_rho)) p%rho=f(l1:l2,m,n,ilnrho)
        if (lpencil(i_rho1)) p%rho1=1.0/p%rho
      else
        if (lpencil(i_rho1)) p%rho1=exp(-f(l1:l2,m,n,ilnrho))
        if (lpencil(i_rho)) p%rho=1.0/p%rho1
      endif
! glnrho and grho
      if (lpencil(i_glnrho).or.lpencil(i_grho)) then
        if (ldensity_nolog) then
          call grad(f,ilnrho,p%grho)
          if (lpencil(i_glnrho)) then
            do i=1,3
              p%glnrho(:,i)=p%grho(:,i)/p%rho
            enddo
          endif
        else
          call grad(f,ilnrho,p%glnrho)
          if (lpencil(i_grho)) then
            do i=1,3
              p%grho(:,i)=p%rho*p%glnrho(:,i)
            enddo
          endif
        endif
      endif
! uglnrho
      if (lpencil(i_uglnrho)) then
        if (ldensity_nolog) then
          call dot(p%uu,p%glnrho,p%uglnrho)
        else
          if (lupw_rho) call stop_it("calc_pencils_density: you switched "//&
               "lupw_rho instead of lupw_lnrho")
          call u_dot_grad(f,ilnrho,p%glnrho,p%uu,p%uglnrho,UPWIND=lupw_lnrho)
        endif
      endif
! ugrho
      if (lpencil(i_ugrho)) then
        if (ldensity_nolog) then
          if (lupw_lnrho) call stop_it("calc_pencils_density: you switched "//&
               "lupw_lnrho instead of lupw_rho")
          call u_dot_grad(f,ilnrho,p%grho,p%uu,p%ugrho,UPWIND=lupw_rho)
        else
          call dot(p%uu,p%grho,p%ugrho)
        endif
      endif
! glnrho2
      if (lpencil(i_glnrho2)) call dot2(p%glnrho,p%glnrho2)
! del2lnrho
      if (lpencil(i_del2lnrho)) then
        if (ldensity_nolog) then
          if (headtt) then
            call fatal_error('calc_pencils_density', &
                             'del2lnrho not available for non-logarithmic mass density')
          endif
        else
          call del2(f,ilnrho,p%del2lnrho)
        endif
      endif
! del2rho
      if (lpencil(i_del2rho)) then
        if (ldensity_nolog) then
          call del2(f,ilnrho,p%del2rho)
        else
          if (headtt) then
            call fatal_error('calc_pencils_density', &
                             'del2rho not available for logarithmic mass density')
          endif
        endif
      endif
! del6rho
      if (lpencil(i_del6rho)) then
        if (ldensity_nolog) then
          call del6(f,ilnrho,p%del6rho)
        else
          if (lfirstpoint) then
            !
            ! Fill global rho array using the ilnrho data
!ajwm This won't work unless earlt_finalize is used... ?
            if (iglobal_rho/=0) then
              do mm=1,my; do nn=1,mz
                f(:,mm,nn,iglobal_rho) = exp(f(:,mm,nn,ilnrho))
              enddo; enddo
            else
              call fatal_error("calc_pencils_density",&
                       "A global rho slot must be available for calculating del6rho from lnrho")
            endif
          endif
          if (iglobal_rho/=0) call del6(f,iglobal_rho,p%del6rho)
        endif
      endif
! del6lnrho
      if (lpencil(i_del6lnrho)) then
        if (ldensity_nolog) then
          if (headtt) then
            call fatal_error('calc_pencils_density', &
                             'del6lnrho not available for non-logarithmic mass density')
          endif
        else
          call del6(f,ilnrho,p%del6lnrho)
        endif
      endif
! hlnrho
      if (lpencil(i_hlnrho)) then
        if (ldensity_nolog) then
          if (headtt) then
            call fatal_error('calc_pencils_density', &
                             'hlnrho not available for non-logarithmic mass density')
          endif
        else
          call g2ij(f,ilnrho,p%hlnrho)
        endif
      endif
! sglnrho
      if (lpencil(i_sglnrho)) call multmv(p%sij,p%glnrho,p%sglnrho)
! uij5glnrho
      if (lpencil(i_uij5glnrho)) call multmv(p%uij5,p%glnrho,p%uij5glnrho)
!
    endsubroutine calc_pencils_density
!***********************************************************************
    subroutine dlnrho_dt(f,df,p)
!
!  continuity equation
!  calculate dlnrho/dt = - u.gradlnrho - divu
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Mpicomm, only: stop_it
      use Special, only: special_calc_density
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: fdiff, gshockglnrho, gshockgrho, tmp
!
      intent(in)  :: f,p
      intent(out) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dlnrho_dt: SOLVE dlnrho_dt'
      if (headtt) call identify_bcs('lnrho',ilnrho)
!
!  continuity equation
!
      if (lcontinuity_gas) then
        if (ldensity_nolog) then
          df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) - p%ugrho - p%rho*p%divu
        else
          df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) - p%uglnrho - p%divu
        endif
      endif
!
!  mass sources and sinks
!
      if (lmass_source) call mass_source(f,df,p)
!
!  Mass diffusion
!
      fdiff=0.0

      if (ldiff_normal) then  ! Normal diffusion operator
        if (ldensity_nolog) then
          fdiff = fdiff + diffrho*p%del2rho
        else
          fdiff = fdiff + diffrho*(p%del2lnrho+p%glnrho2)
        endif
        if (lfirst.and.ldt) diffus_diffrho=diffus_diffrho+diffrho*dxyz_2
        if (headtt) print*,'dlnrho_dt: diffrho=', diffrho
      endif
!
!  Hyper diffusion
!
      if (ldiff_hyper3) then
        if (ldensity_nolog) then
          fdiff = fdiff + diffrho_hyper3*p%del6rho
        else
          fdiff = fdiff + 1/p%rho*diffrho_hyper3*p%del6rho
        endif
        if (lfirst.and.ldt) diffus_diffrho=diffus_diffrho+diffrho_hyper3*dxyz_6
        if (headtt) print*,'dlnrho_dt: diffrho_hyper3=', diffrho_hyper3
      endif
!
      if (ldiff_hyper3_aniso) then
         if (ldensity_nolog) then
            call del6fj(f,diffrho_hyper3_aniso,ilnrho,tmp)
            fdiff = fdiff + tmp
            if (lfirst.and.ldt) diffus_diffrho=diffus_diffrho+&
                 diffrho_hyper3_aniso(1)*dx_1(l1:l2)**6 + &
                 diffrho_hyper3_aniso(2)*dy_1(m)**6 + &
                 diffrho_hyper3_aniso(3)*dz_1(n)**6
            if (headtt) &
                 print*,'dlnrho_dt: diffrho_hyper3=(Dx,Dy,Dz)=',diffrho_hyper3_aniso
         else
            call stop_it("anisotropic hyperdiffusion not implemented for lnrho")
         endif
      endif
!
      if (ldiff_hyper3lnrho) then
        if (.not. ldensity_nolog) then
          fdiff = fdiff + diffrho_hyper3*p%del6lnrho
        endif
        if (lfirst.and.ldt) diffus_diffrho=diffus_diffrho+diffrho_hyper3*dxyz_6
        if (headtt) print*,'dlnrho_dt: diffrho_hyper3=', diffrho_hyper3
      endif
!
!  Shock diffusion
!
      if (ldiff_shock) then
        if (ldensity_nolog) then
          call dot_mn(p%gshock,p%grho,gshockgrho)
          df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + &
              diffrho_shock*p%shock*p%del2rho + diffrho_shock*gshockgrho
        else
          call dot_mn(p%gshock,p%glnrho,gshockglnrho)
          df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + &
              diffrho_shock*p%shock*(p%del2lnrho+p%glnrho2) + &
              diffrho_shock*gshockglnrho
        endif
        if (lfirst.and.ldt) &
            diffus_diffrho=diffus_diffrho+diffrho_shock*p%shock*dxyz_2
        if (headtt) print*,'dlnrho_dt: diffrho_shock=', diffrho_shock
      endif
!
!  Add diffusion term to continuity equation
!
      df(l1:l2,m,n,ilnrho) = df(l1:l2,m,n,ilnrho) + fdiff
      if (headtt.or.ldebug) &
          print*,'dlnrho_dt: max(diffus_diffrho) =', maxval(diffus_diffrho)
!
!
!
!ajwm  Cannot alter special interface!!
      if (lspecial) call special_calc_density(f,df,p)
!
!  Apply border profile
!
      if (lborder_profiles) call set_border_density(f,df,p)
!
!  2d-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        if (idiag_lnrhomphi/=0) call phisum_mn_name_rz(p%lnrho,idiag_lnrhomphi)
        if (idiag_rhomphi/=0)   call phisum_mn_name_rz(p%rho,idiag_rhomphi)
        if (idiag_rhomxz/=0)    call ysum_mn_name_xz(p%rho,idiag_rhomxz)
        if (idiag_rhomxy/=0)    call zsum_mn_name_xy(p%rho,idiag_rhomxy)
      endif
!
!  1d-averages. Happens at every it1d timesteps, NOT at every it1
!
      if (l1ddiagnos) then
         if (idiag_rhomr/=0)    call phizsum_mn_name_r(p%rho,idiag_rhomr)
         if (idiag_rhomz/=0)    call xysum_mn_name_z(p%rho,idiag_rhomz)
         if (idiag_rhomx/=0)    call yzsum_mn_name_x(p%rho,idiag_rhomx)
         if (idiag_rhomy/=0)    call xzsum_mn_name_y(p%rho,idiag_rhomy)
      endif
!
!  Calculate density diagnostics
!
      if (ldiagnos) then
        if (idiag_rhom/=0)     call sum_mn_name(p%rho,idiag_rhom)
        if (idiag_totmass/=0)  call integrate_mn_name(p%rho,idiag_totmass)
        if (idiag_rhomin/=0) &
            call max_mn_name(-p%rho,idiag_rhomin,lneg=.true.)
        if (idiag_rhomax/=0)   call max_mn_name(p%rho,idiag_rhomax)
        if (idiag_rho2m/=0)    call sum_mn_name(p%rho**2,idiag_rho2m)
        if (idiag_lnrho2m/=0)  call sum_mn_name(p%lnrho**2,idiag_lnrho2m)
        if (idiag_uglnrhom/=0) call sum_mn_name(p%uglnrho,idiag_uglnrhom)
        if (idiag_dtd/=0) &
            call max_mn_name(diffus_diffrho/cdtv,idiag_dtd,l_dt=.true.)
      endif
!
    endsubroutine dlnrho_dt
!***********************************************************************
    subroutine set_border_density(f,df,p)
!
!  Calculates the driving term for the border profile
!  of the lnrho variable.
!
!  28-jul-06/wlad: coded
!
      use BorderProfiles,  only: border_driving
      use EquationOfState, only: cs0,cs20
      use Sub,             only: power_law
      use Mpicomm,         only: stop_it
      use Gravity,         only: potential
      use FArrayManager
!
      real, dimension(mx,my,mz,mfarray) :: f
      real, dimension(mx,my,mz,mvar) :: df
      real, dimension(nx) :: f_target!,lnrhomid,pot,tmp1,tmp2
      real :: lnrhomid,pot,tmp1,tmp2
      type (pencil_case)  :: p
      integer            :: i

      select case(borderlnrho)

      case('zero','0')
        if (plaw.ne.0) call stop_it("borderlnrho: density is not flat but "//&
             "you are calling zero border")
        if (ldensity_nolog) then  
          f_target=0.
        else
          f_target=1.
        endif
      case('constant')
        if (plaw.ne.0) call stop_it("borderlnrho: density is not flat but "//&
             "you are calling constant border")
        if (ldensity_nolog) then 
          f_target=rho_const
        else
          f_target=lnrho_const
        endif
      case('power-law')
        if (plaw.eq.0) call stop_it("borderlnrho: no need to call a power-"//&
             "law border for a flat density profile")
        if (ldensity_nolog) then
          call power_law(rho_const,p%rcyl_mn,plaw,f_target)
!          f_target=rho_const*p%rcyl_mn1**(plaw)
        else
          f_target=lnrho_const - plaw*log(p%rcyl_mn)
        endif
      case('stratification')

        do i=1,nx
          if ( ((p%rcyl_mn(i).ge.r_int).and.(p%rcyl_mn(i).le.r_int+2*wborder_int)).or.&
               ((p%rcyl_mn(i).ge.r_ext-2*wborder_ext).and.(p%rcyl_mn(i).le.r_ext))) then
            
            lnrhomid=log(rho0) - plaw*log(p%rcyl_mn(i))
            call potential(R=p%r_mn(i),POT=tmp1)
            call potential(R=p%rcyl_mn(i),POT=tmp2)
            pot=-gamma*(tmp1-tmp2)/p%cs2(i)

            f_target(i)=lnrhomid+pot

            if (ldensity_nolog) &
                 f_target(i)=exp(f_target(i))
          endif
        enddo

     case('nothing')
        if (lroot.and.ip<=5) &
             print*,"set_border_lnrho: borderlnrho='nothing'"
      case default
        write(unit=errormsg,fmt=*) &
             'set_border_lnrho: No such value for borderlnrho: ', &
             trim(borderlnrho)
        call fatal_error('set_border_lnrho',errormsg)
      endselect
!
      if (borderlnrho /= 'nothing') then
        call border_driving(f,df,p,f_target,ilnrho)
      endif
!
    endsubroutine set_border_density
!***********************************************************************
    subroutine rprint_density(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Sub
!
      logical :: lreset
      logical, optional :: lwrite
!
      integer :: iname, inamex, inamey, inamez, inamexy, inamexz, irz, inamer
      logical :: lwr
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_rhom=0; idiag_rho2m=0; idiag_lnrho2m=0; idiag_uglnrhom=0
        idiag_rhomin=0; idiag_rhomax=0; idiag_dtd=0
        idiag_lnrhomphi=0; idiag_rhomphi=0
        idiag_rhomz=0; idiag_rhomy=0; idiag_rhomx=0 
        idiag_rhomxy=0; idiag_rhomr=0; idiag_totmass=0
        idiag_rhomxz=0
        cdiffrho=0.
        diffrho=0.
      endif
!
!  iname runs through all possible names that may be listed in print.in
!
      if(lroot.and.ip<14) print*,'rprint_density: run through parse list'
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'rhom',idiag_rhom)
        call parse_name(iname,cname(iname),cform(iname),'rho2m',idiag_rho2m)
        call parse_name(iname,cname(iname),cform(iname),'rhomin',idiag_rhomin)
        call parse_name(iname,cname(iname),cform(iname),'rhomax',idiag_rhomax)
        call parse_name(iname,cname(iname),cform(iname),'lnrho2m',idiag_lnrho2m)
        call parse_name(iname,cname(iname),cform(iname),'uglnrhom',idiag_uglnrhom)
        call parse_name(iname,cname(iname),cform(iname),'dtd',idiag_dtd)
        call parse_name(iname,cname(iname),cform(iname),'totmass',idiag_totmass)
      enddo
!
!  check for those quantities for which we want xy-averages
!
      do inamez=1,nnamez
        call parse_name(inamez,cnamez(inamez),cformz(inamez),'rhomz',idiag_rhomz)
      enddo
!
!  check for those quantities for which we want xz-averages
!
      do inamey=1,nnamey
        call parse_name(inamey,cnamey(inamey),cformy(inamey),'rhomy',idiag_rhomy)
      enddo
!
!  check for those quantities for which we want yz-averages
!
      do inamex=1,nnamex
        call parse_name(inamex,cnamex(inamex),cformx(inamex),'rhomx',idiag_rhomx)
      enddo
!
!  check for those quantities for which we want phiz-averages
!
      do inamer=1,nnamer
        call parse_name(inamer,cnamer(inamer),cformr(inamer),'rhomr',idiag_rhomr)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do inamexz=1,nnamexz
        call parse_name(inamexz,cnamexz(inamexz),cformxz(inamexz),'rhomxz',idiag_rhomxz)
      enddo
!
!  check for those quantities for which we want z-averages
!
      do inamexy=1,nnamexy
        call parse_name(inamexy,cnamexy(inamexy),cformxy(inamexy),'rhomxy',idiag_rhomxy)
      enddo
!
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),&
            'lnrhomphi',idiag_lnrhomphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rhomphi',idiag_rhomphi)
      enddo
!
!  write column where which density variable is stored
!
      if (lwr) then
        write(3,*) 'i_rhom=',idiag_rhom
        write(3,*) 'i_rho2m=',idiag_rho2m
        write(3,*) 'i_rhomin=',idiag_rhomin
        write(3,*) 'i_rhomax=',idiag_rhomax
        write(3,*) 'i_lnrho2m=',idiag_lnrho2m
        write(3,*) 'i_uglnrhom=',idiag_uglnrhom
        write(3,*) 'i_rhomz=',idiag_rhomz
        write(3,*) 'i_rhomy=',idiag_rhomy
        write(3,*) 'i_rhomx=',idiag_rhomx
        write(3,*) 'i_rhomxy=',idiag_rhomxy
        write(3,*) 'i_rhomxz=',idiag_rhomxz
        write(3,*) 'nname=',nname
        write(3,*) 'ilnrho=',ilnrho
        write(3,*) 'i_lnrhomphi=',idiag_lnrhomphi
        write(3,*) 'i_rhomphi=',idiag_rhomphi
        write(3,*) 'i_rhomr=',idiag_rhomr
        write(3,*) 'i_dtd=',idiag_dtd
        write(3,*) 'i_totmass=',idiag_totmass
      endif
!
    endsubroutine rprint_density
!***********************************************************************
!  Here comes a collection of different density stratification routines
!***********************************************************************
    subroutine isothermal_density(f)
!
!  Isothermal stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!  When entropy is present, this module also initializes entropy.
!
!  Sound speed (and hence Temperature), and density (at infinity) are
!  initialised to their respective reference values:
!           sound speed: cs^2_0            from start.in
!           density: rho0 = exp(lnrho0)
!
!   8-jul-02/axel: incorporated/adapted from init_lnrho
!  11-jul-02/axel: fixed sign; should be tmp=gamma*pot/cs20
!  02-apr-03/tony: made entropy explicit rather than using tmp/-gamma
!  11-jun-03/tony: moved entropy initialisation to separate routine
!                  to allow isothermal condition for arbitrary density
!
      use Gravity
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: pot,tmp
      real :: cp1
!
!  Stratification depends on the gravity potential
!
      if (lroot) print*,'isothermal_density: isothermal stratification'
      if (gamma/=1.0) then
        if ((.not. lentropy) .and. (.not. ltemperature)) & 
          call fatal_error('isothermal_density','for gamma/=1.0, you need entropy or temperature!');
      endif
!
      call get_cp1(cp1)
      do n=n1,n2
        do m=m1,m2
          call potential(x(l1:l2),y(m),z(n),pot=pot)
          tmp=-gamma*pot/cs20
          f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + lnrho0 + tmp
          if (lentropy) f(l1:l2,m,n,iss) = f(l1:l2,m,n,iss) &
               -gamma1*(f(l1:l2,m,n,ilnrho)-lnrho0)/gamma
          if (ltemperature) f(l1:l2,m,n,ilnTT)=log(cs20*cp1/gamma1)
        enddo
      enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  The values calculated here may be revised in the entropy module.
!
      cs2bot=cs20
      cs2top=cs20
!
    endsubroutine isothermal_density
!***********************************************************************
    subroutine local_isothermal_density(f)
!                                                                   
!  Stratification depends on the gravity potential, which in turn   
!  varies with radius. This reproduces the initial condition of the 
!  locally isothermal approximation in which temperature is a power 
!  law of radial distance
!
!  18-apr-07/wlad : coded
!
      use FArrayManager
      use Mpicomm
      use Initcond,only:set_thermodynamical_quantities
      use Gravity, only:potential
      use Sub,     only:power_law,grad,get_radial_distance
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx,3) :: glnrho
      real, dimension (nx)   :: pot,tmp1,tmp2,corr,gslnrho
      real, dimension (nx)   :: rr_sph,rr,rr_cyl,lnrhomid
      real                   :: ptlaw,g0_
      integer, pointer       :: iglobal_cs2,iglobal_glnTT
      integer                :: i
      logical                :: lheader
!
      if (lroot) print*,'isothermal_density: local isothermal stratification'
      if (lroot) print*,'Radial stratification with power law=',plaw
!
      call get_ptlaw(ptlaw)
      call set_thermodynamical_quantities(f,iglobal_cs2,iglobal_glnTT,ptlaw)
!
      do n=n1,n2
        do m=m1,m2
          lheader=lroot.and.(m==m1).and.(n==n1)
          call get_radial_distance(rr)
          if (lcartesian_coords) then 
            rr_cyl=sqrt(x(l1:l2)**2+y(m)**2)
            rr_sph=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          elseif (lcylindrical_coords) then
            rr_cyl=rr
            rr_sph=sqrt(rr**2+z(n)**2)
          endif

          !lnrhomid=log(rho0)-plaw*log(rr_cyl)
          lnrhomid=log(rho0)-.5*plaw*log(rr_cyl**2+rsmooth**2)
          if (.not.lcylindrical_gravity) then 
            if (lheader) &
                 print*,'Adding vertical stratification with scale height h/r=',cs0
!
!  The subroutine "potential" yields the whole gradient.
!  I want the function that partially derived in 
!  z gives g0/r^3 * z. This is NOT -g0/r
!  The second call takes care of normalizing it 
!  i.e., there should be no correction at midplane
!
            if (.not.lcartesian_coords) &
                 call stop_it("fix the potential for cylindrical coords")
            call potential(x(l1:l2),y(m),z(n),POT=tmp1,RMN=rr_sph)
            call potential(x(l1:l2),y(m),z(n),POT=tmp2,RMN=rr_cyl)
!
            pot=-(tmp1-tmp2)/f(l1:l2,m,n,iglobal_cs2)
            if (ltemperature) pot=gamma*pot
          else 
            pot=0.
          endif
          f(l1:l2,m,n,ilnrho) = lnrhomid+pot
        enddo
      enddo
! 
!  Correct for density gradient term in the centrifugal force. The
!  temperature gradient was already corrected for when setting the
!  thermodynamical quantities.   
! 
!  Had to split in two loops because I need the (log)density to 
!  be fully coded before taking its gradient in y or z
!      
      do m=m1,m2
        do n=n1,n2
!
          call grad(f,ilnrho,glnrho)
          call get_radial_distance(rr)
          if (lcartesian_coords) then 
            rr_cyl=sqrt(x(l1:l2)**2+y(m)**2)
            rr_sph=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          elseif (lcylindrical_coords) then
            rr_cyl=rr
            rr_sph=sqrt(rr**2+z(n)**2)
          endif
          !gs= gx*cos + gy*sin
          if (lcartesian_coords) then
            gslnrho=(glnrho(:,1)*x(l1:l2) + glnrho(:,2)*y(m))/rr_cyl
          else
            gslnrho=glnrho(:,1)
          endif
          corr=gslnrho*f(l1:l2,m,n,iglobal_cs2)  
          if (ltemperature) corr=corr/gamma
!
          if (lcartesian_coords) then
            tmp1=(f(l1:l2,m,n,iux)**2+f(l1:l2,m,n,iuy)**2)/rr_cyl**2
          elseif (lcylindrical_coords) then
            tmp1=(f(l1:l2,m,n,iuy)/rr_cyl)**2
          endif
!
          tmp2=tmp1 + corr/rr_cyl
!            
          do i=1,nx
            if (tmp2(i).lt.0.) then
              if (rr_cyl(i) .lt. r_int) then
                !it's inside the frozen zone, so 
                !just set tmp2 to zero and emit a warning
                tmp2(i)=0.
                if (ip<=10) call warning('local_isothermal_density',&
                     'the density gradient is too steep inside the frozen zone')
              else
                print*,'local_isothermal_density: the density gradient '
                print*,'is too steep at x,y,z=',x(i+l1-1),y(m),z(n)
                call stop_it("")
              endif
            endif
          enddo
          if (lcartesian_coords) then
            f(l1:l2,m,n,iux)=-sqrt(tmp2)*y(  m  )
            f(l1:l2,m,n,iuy)= sqrt(tmp2)*x(l1:l2)
          elseif (lcylindrical_coords) then
            f(l1:l2,m,n,iuy)= sqrt(tmp2)*rr_cyl
          endif
        enddo
      enddo
!
    endsubroutine local_isothermal_density
!***********************************************************************
    subroutine polytropic_simple(f)
!
!  Polytropic stratification (for lnrho and ss)
!  This routine should be independent of the gravity module used.
!
!  To maintain continuity with respect to the isothermal case,
!  one may want to specify cs20 (=1), and so zinfty is calculated from that.
!  On the other hand, for polytropic atmospheres it may be more
!  physical to specify zinfty (=1), ie the top of the polytropic atmosphere.
!  This is done if zinfty is different from 0.
!
!   8-jul-02/axel: incorporated/adapted from init_lnrho
!
      use Gravity, only: grav_profile,gravz,zinfty,zref,zgrav,  &
                             potential,nu_epicycle
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: pot,dlncs2,r_mn
      real :: ggamma,ztop,zbot,zref2,pot_ext,lnrho_ref,ptop,pbot
!
!  identifier
!
      if (lroot) print*,'polytropic_simple: mpoly=',mpoly
!
!  The following is specific only to cases with gravity in the z direction
!  zref is calculated such that rho=rho0 and cs2=cs20 at z=zref.
!  Note: gravz is normally negative!
!
      if (lgravz) then
        if (grav_profile=='const') then
          if(lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zref=zinfty-(mpoly+1.)*cs20/(-gamma*gravz)
        elseif (grav_profile=='const_zero') then
          if(lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zref=zinfty-(mpoly+1.)*cs20/(-gamma*gravz)
        elseif (grav_profile=='linear') then
          if(lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zref2=zinfty**2-(mpoly+1.)*cs20/(0.5*gamma*nu_epicycle**2)
          if(zref2<0) then
            if(lroot) print*,'polytropic_simple: zref**2<0 is not ok'
            zref2=0. !(and see what happens)
          endif
          zref=sqrt(zref2)
        else
          if(lroot) print*,'polytropic_simple: zref not prepared!'
        endif
        if (lroot) print*,'polytropic_simple: zref=',zref
!
!  check whether zinfty lies outside the domain (otherwise density
!  would vanish within the domain). At the moment we are not properly
!  testing the lower boundary on the case of a disc (commented out).
!
        ztop=xyz0(3)+Lxyz(3)
        zbot=xyz0(3)
        !-- if(zinfty<min(ztop,zgrav) .or. (-zinfty)>min(zbot,zgrav)) then
        if(zinfty<min(ztop,zgrav)) then
          if(lroot) print*,'polytropic_simple: domain too big; zinfty=',zinfty
          !call stop_it( &
          !         'polytropic_simply: rho and cs2 will vanish within domain')
        endif
      endif
!
!  stratification Gamma (upper case in the manual)
!
      ggamma=1.+1./mpoly
!
!  polytropic sphere with isothermal exterior
!  calculate potential at the stellar surface, pot_ext
!
      if (lgravr) then
        call potential(R=r_ext,POT=pot_ext)
        cs2top=-gamma/(mpoly+1.)*pot_ext
        lnrho_ref=mpoly*log(cs2top)-(mpoly+1.)
        print*,'polytropic_simple: pot_ext=',pot_ext
        do n=n1,n2
        do m=m1,m2
          r_mn=sqrt(x(l1:l2)**2+y(m)**2+z(n)**2)
          call potential(x(l1:l2),y(m),z(n),pot=pot)
!
!  density
!  these formulae assume lnrho0=0 and cs0=1
!
          where (r_mn > r_ext)
            f(l1:l2,m,n,ilnrho)=lnrho_ref-gamma*pot/cs2top
          elsewhere
            dlncs2=log(-gamma*pot/((mpoly+1.)*cs20))
            f(l1:l2,m,n,ilnrho)=lnrho0+mpoly*dlncs2
          endwhere
!
!  entropy
!
          if(lentropy) then
            where (r_mn > r_ext)
              f(l1:l2,m,n,iss)=-(1.-1./gamma)*f(l1:l2,m,n,ilnrho)+log(cs2top)/gamma
            elsewhere
              dlncs2=log(-gamma*pot/((mpoly+1.)*cs20))
              f(l1:l2,m,n,iss)=mpoly*(ggamma/gamma-1.)*dlncs2
            endwhere
          endif
        enddo
        enddo
      else
!
!  cartesian case with gravity in the z direction
!
        do n=n1,n2
        do m=m1,m2
          call potential(x(l1:l2),y(m),z(n),pot=pot)
          dlncs2=log(-gamma*pot/((mpoly+1.)*cs20))
          f(l1:l2,m,n,ilnrho)=lnrho0+mpoly*dlncs2
          if(lentropy) f(l1:l2,m,n,iss)=mpoly*(ggamma/gamma-1.)*dlncs2
!         if(ltemperature) f(l1:l2,m,n,ilnTT)=dlncs2-log(gamma1)
          if(ltemperature) f(l1:l2,m,n,ilnTT)=log(-gamma*pot/(mpoly+1.)/gamma1)
        enddo
        enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  In spherical geometry, ztop is z at the outer edge of the box,
!  so this calculation still makes sense.
!
        call potential(xyz0(1),xyz0(2),ztop,pot=ptop)
        cs2top=-gamma*ptop/(mpoly+1.)
!
!  In spherical geometry ztop should never be used.
!  Even in slab geometry ztop is not normally used.
!
        call potential(xyz0(1),xyz0(2),zbot,pot=pbot)
        cs2bot=-gamma*pbot/(mpoly+1.)
      endif
!
    endsubroutine polytropic_simple
!***********************************************************************
    subroutine mass_source(f,df,p)
!
!  add mass sources and sinks
!
!  28-apr-2005/axel: coded
!
      use Cdata
      use Sub, only: step
      use Gravity, only: lnrho_bot,lnrho_top,ss_bot,ss_top
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      real, dimension(nx) :: fint,fext,pdamp
!
      if(ldebug) print*,'mass_source: cs20,cs0=',cs20,cs0
!
!  cylindrical profile for inner cylinder
!
      pdamp=1-step(p%rcyl_mn,r_int,wdamp) ! inner damping profile
      fint=-damplnrho_int*pdamp*(f(l1:l2,m,n,ilnrho)-lnrho_int)
      df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+fint
!
!  cylindrical profile for outer cylinder
!
      pdamp=step(p%rcyl_mn,r_ext,wdamp) ! outer damping profile
      fext=-damplnrho_ext*pdamp*(f(l1:l2,m,n,ilnrho)-lnrho_ext)
      df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+fext
!     df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+f(l1:l2,m,n,iux)*fext
!     df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+f(l1:l2,m,n,iuy)*fext
!
    endsubroutine mass_source
!***********************************************************************
    subroutine get_plaw(plaw_)
!
      real :: plaw_
!
      plaw_=plaw
!
    endsubroutine get_plaw
!***********************************************************************
    subroutine cylind_poly(f)
!
!  Initialize density and entropy/temperature on a cylindrical grid
!  using a polytrope under a constant (radial) gravity
!
!  07-mar-07/dintrans: coded
!
      use Gravity, only: g0
      use EquationOfState
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (nx) :: temp,lnTT,lnrho,ss
      real :: beta1,temp_top,temp_bot
!
      mpoly0=mpoly ! to be consistent with the Fbot definition
      beta1=-g0/(mpoly+1)*gamma/gamma1  ! gamma1/gamma=R_{*} (for cp=1)
      temp_top=cs20/gamma1
      temp_bot=temp_top+beta1*(r_int-r_ext)
      cs2top=cs20
      cs2bot=gamma1*temp_bot
      if (lroot) &
        print*,'cylind_poly: beta1,temp_top,r_ext=',beta1,temp_top,r_ext
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
!
        temp=temp_top+beta1*(x(l1:l2)-r_ext)
        lnTT=alog(temp)
        f(l1:l2,m,n,ilnrho)=lnrho0+mpoly*lnTT-mpoly*log(temp_top)
        if (lentropy) then
          lnrho=f(l1:l2,m,n,ilnrho)
          call eoscalc(ilnrho_lnTT,lnrho,lnTT,ss=ss)
          f(l1:l2,m,n,iss)=ss
        endif
        if (ltemperature) f(l1:l2,m,n,ilnTT)=lnTT
      enddo
!
    endsubroutine cylind_poly
!***********************************************************************
    subroutine impose_density_floor(f)
!
!  Impose a minimum density by setting all lower densities to the minimum
!  value (density_floor). Useful for debugging purposes.
!
!  13-aug-2007/anders: implemented.
!
      use Cdata
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
!
      real :: density_floor_log
      logical, save :: lfirstcall=.true.
!
      if (density_floor>0.) then
        if (lfirstcall) then
          density_floor_log=alog(density_floor)
          lfirstcall=.false.
        endif
!
        if (ldensity_nolog) then
          where (f(:,:,:,ilnrho)<density_floor) &
              f(:,:,:,ilnrho)=density_floor
        else
          where (f(:,:,:,ilnrho)<density_floor_log) &
              f(:,:,:,ilnrho)=density_floor_log
        endif
      endif
!
    endsubroutine impose_density_floor
!***********************************************************************
endmodule Density
