! $Id: density.f90,v 1.218 2006-02-02 13:13:01 wlyra Exp $

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
                             mpoly,beta_glnrho_global

  implicit none

  include 'density.h'

  real :: ampllnrho=0.,widthlnrho=.1,rho_left=1.,rho_right=1.
  real :: cdiffrho=0.,diffrho=0.,diffrho_hyper3=0.,diffrho_shock=0.
  real :: lnrho_const=0., rho_const=1.
  real :: amplrho=0, phase_lnrho=0.0
  real :: radius_lnrho=.5,kx_lnrho=1.,ky_lnrho=1.,kz_lnrho=1.
  real :: eps_planet=.5
  real :: q_ell=5., hh0=0.
  real :: co1_ss=0.,co2_ss=0.,Sigma1=150.
  real :: lnrho_int=0.,lnrho_ext=0.,damplnrho_int=0.,damplnrho_ext=0.
  real :: wdamp=0.,plaw=0.
  integer, parameter :: ndiff_max=4
  logical :: lupw_lnrho=.false.,lmass_source=.false.,lcontinuity_gas=.true.
  logical :: ldiff_normal=.false.,ldiff_hyper3=.false.,ldiff_shock=.false.
  logical :: ldiff_hyper3lnrho=.false.
  logical :: lfreeze_lnrhoint=.false.,lfreeze_lnrhoext=.false.

  character (len=labellen), dimension(ninit) :: initlnrho='nothing'
  character (len=labellen) :: strati_type='lnrho_ss',initlnrho2='nothing'
  character (len=labellen), dimension(ndiff_max) :: idiff=''
  character (len=4) :: iinit_str
  complex :: coeflnrho=0.

  namelist /density_init_pars/ &
       ampllnrho,initlnrho,initlnrho2,widthlnrho,    &
       rho_left,rho_right,lnrho_const,cs2bot,cs2top, &
       radius_lnrho,eps_planet,                      &
       b_ell,q_ell,hh0,rbound,                       &
       mpoly,strati_type,beta_glnrho_global,         &
       kx_lnrho,ky_lnrho,kz_lnrho,amplrho,phase_lnrho,coeflnrho, &
       co1_ss,co2_ss,Sigma1,idiff,ldensity_nolog,    &
       wdamp,plaw,lcontinuity_gas

  namelist /density_run_pars/ &
       cdiffrho,diffrho,diffrho_hyper3,diffrho_shock,   &
       cs2bot,cs2top,lupw_lnrho,idiff,lmass_source,     &
       lnrho_int,lnrho_ext,damplnrho_int,damplnrho_ext, &
       wdamp,lfreeze_lnrhoint,lfreeze_lnrhoext,         &
       lnrho_const,plaw,lcontinuity_gas
  
  ! diagnostic variables (needs to be consistent with reset list below)
  integer :: idiag_rhom=0,idiag_rho2m=0,idiag_lnrho2m=0
  integer :: idiag_rhomin=0,idiag_rhomax=0
  integer :: idiag_lnrhomphi=0,idiag_rhomphi=0,idiag_dtd=0
  integer :: idiag_rhomz=0, idiag_rhomy=0, idiag_rhomx=0

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
!
      logical, save :: first=.true.
!
      if (.not. first) call fatal_error('register_density','module registration called twice')
      first = .false.
!
!ajwm      ldensity = .true.
!
      ilnrho = nvar+1           ! indix to access lam
      nvar = nvar+1             ! added 1 variable
!
      if ((ip<=8) .and. lroot) then
        print*, 'register_density: nvar = ', nvar
        print*, 'register_density: ilnrho = ', ilnrho
      endif
!
!  Put variable name in array
!
      varname(ilnrho) = 'lnrho'
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: density.f90,v 1.218 2006-02-02 13:13:01 wlyra Exp $")
!
      if (nvar > mvar) then
        if (lroot) write(0,*) 'nvar = ', nvar, ', mvar = ', mvar
        call fatal_error('register_density','nvar > mvar')
      endif
!
!  Writing files for use with IDL
!
      if (lroot) then
        if (maux == 0) then
          if (nvar < mvar) write(4,*) ',lnrho $'
          if (nvar == mvar) write(4,*) ',lnrho'
        else
          write(4,*) ',lnrho $'
        endif
        write(15,*) 'lnrho = fltarr(mx,my,mz)*one'
      endif
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

      use CData, only: lfreeze_varext,lfreeze_varint,ilnrho

      real, dimension (mx,my,mz,mvar+maux) :: f
      logical :: lstarting
!
      integer :: i
      logical :: lnothing
!
!   make sure all relevant parameters are set for spherical shell problems
!
      select case(initlnrho(1))
        case('geo-kws')
          if (lroot) print*,'initialize_density: set reference sound speed for spherical shell problem'
!ajwm Shouldn't be done here they should be set as input parameters in
!ajwm start.in.  There may be a problem at present with run since these
!ajwm vales will not get set consistantly.
!ajwm cs0 and lnrho0 are parameters of the EOS not of density
!        cs20=gamma
!        cs0=sqrt(cs20)
!        cp=5./2.
!ajwm Routine should really be called initialize_eos eventually
!ajwm important thing is the parameters of the Equation of State 
!ajwm have changed so we'd better recalculate everything consistently
!ajwm but that won't work either since 
!        call initialize_ionization()
      endselect
!
!  initialize cs2cool to cs20
!  (currently disabled, because it causes problems with mdarf auto-test)
!     cs2cool=cs20
!
!
      if(diffrho==0.) then
!
!  FIXME: This will not work with RELOAD if cdiffrho is changed in run.in
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
      if (ldiff_normal.and. diffrho==0.0) then
        call warning('initialize_density','Diffusion coefficient diffrho is zero!')
        ldiff_normal=.false.
      endif
      if ( (ldiff_hyper3 .or. ldiff_hyper3lnrho) .and. diffrho_hyper3==0.0) then
        call warning('initialize_density','Diffusion coefficient diffrho_hyper3 is zero!')
        ldiff_hyper3=.false.
      endif
      if (ldiff_shock .and. diffrho_shock==0.0) then
        call warning('initialize_density','diffusion coefficient diffrho_shock is zero!')
        ldiff_shock=.false.
      endif
!
      if (ldiff_hyper3 .and. (.not. ldensity_nolog) .and. &
          (.not. lglobal_nolog_density) ) then
         call fatal_error('initialize_density','must have '// &
             'global_nolog_density module for del6rho with '// &
             'logarithmic density')
      endif
!
      if (lfreeze_lnrhoint) lfreeze_varint(ilnrho) = .true.
      if (lfreeze_lnrhoext) lfreeze_varext(ilnrho) = .true.
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
      use General, only: chn
      use Global
      use Gravity, only: zref,z1,z2,gravz,nu_epicycle,potential, &
                          lnumerical_equilibrium
      use Initcond
      use Initcond_spec
      use IO
      use Mpicomm
      use Sub
      use EquationOfState
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz,tmp,pot,prof
      real :: lnrhoint,cs2int,pot0,lnrho_left,lnrho_right
      real :: pot_ext,lnrho_ext,cs2_ext,tmp1
      real :: zbot,ztop
      logical :: lnothing
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
      case('constant'); f(:,:,:,ilnrho)=log(rho_left)
      case('mode'); call modes(ampllnrho,coeflnrho,f,ilnrho,kx_lnrho,ky_lnrho,kz_lnrho,xx,yy,zz)
      case('blob'); call blob(ampllnrho,f,ilnrho,radius_lnrho,0.,0.,0.)
      case('isothermal'); call isothermal_density(f)
      case('stratification'); call stratification(f,strati_type)
      case('polytropic_simple'); call polytropic_simple(f)
      case('hydrostatic-z', '1'); print*, &
                              'init_lnrho: use polytropic_simple instead!'
      case('xjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'x')
      case('yjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'y')
      case('zjump'); call jump(f,ilnrho,lnrho_left,lnrho_right,widthlnrho,'z')
      case('soundwave-x'); call soundwave(ampllnrho,f,ilnrho,kx=1.)
      case('soundwave-y'); call soundwave(ampllnrho,f,ilnrho,ky=1.)
      case('soundwave-z'); call soundwave(ampllnrho,f,ilnrho,kz=1.)
      case('sinwave-phase'); call sinwave_phase(f,ilnrho,ampllnrho,kx_lnrho,ky_lnrho,kz_lnrho,phase_lnrho)
      case('sinwave-x'); call sinwave(ampllnrho,f,ilnrho,kx=kx_lnrho)
      case('sinwave-y'); call sinwave(ampllnrho,f,ilnrho,ky=ky_lnrho)
      case('sinwave-z'); call sinwave(ampllnrho,f,ilnrho,kz=kz_lnrho)
      case('coswave-x'); call coswave(ampllnrho,f,ilnrho,kx=kx_lnrho)
      case('coswave-y'); call coswave(ampllnrho,f,ilnrho,ky=ky_lnrho)
      case('coswave-z'); call coswave(ampllnrho,f,ilnrho,kz=kz_lnrho)
      case('corona'); call corona_init(f)
      case('gaussian3d'); call gaussian3d(ampllnrho,f,ilnrho,xx,yy,zz,radius_lnrho)
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
!          call initialize_gravity()     ! get coefficients cpot(1:5)
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

      case('solar_nebula')
         !
         !minimum mass solar nebula
         !
         if (lroot)  print*,'init_lnrho: initialize initial condition for planet building'
         call power_law(f,xx,yy,zz,lnrho_const,plaw)

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
        STOP "init_lnrho: Imaginary density values"
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: tmp,p,zz
      real, dimension (mz) :: stp
      real :: mpoly,zint,zbot,zblend,beta1,cs2int,lnrhoint
      integer :: isoth
!
      ! NB: beta1 is not dT/dz, but dcs2/dz = (gamma-1)c_p dT/dz
      if (isoth /= 0) then ! isothermal layer
        beta1 = 0.
        tmp = gamma*gravz/cs2int*(zz-zint)
        lnrhoint = lnrhoint + gamma*gravz/cs2int*(zbot-zint)
      else
        beta1 = gamma*gravz/(mpoly+1)
        tmp = 1 + beta1*(zz-zint)/cs2int
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      use Sub, only: calc_unitvects_sphere
!
      real, dimension (mx,my,mz,mvar+maux), intent(inout) :: f
      real, dimension (nx) :: pot
      real :: beta1,lnrho_int,lnrho_ext,pot_int,pot_ext
!
      beta1 = g0/(mpoly+1)
!
!     densities at shell boundaries
      lnrho_int=mpoly*log(1+beta1*(1/r_int-1))
      lnrho_ext=lnrho0
!
      do imn=1,ny*nz
        n=nn(imn)
        m=mm(imn)
!
        call calc_unitvects_sphere()
!
        ! in the fluid shell
        where (r_mn < r_ext .AND. r_mn > r_int) f(l1:l2,m,n,ilnrho)=mpoly*log(1+beta1*(1/r_mn-1))
        ! outside the fluid shell
        if (initlnrho(1)=='geo-kws') then
          where (r_mn >= r_ext) f(l1:l2,m,n,ilnrho)=lnrho_ext
          where (r_mn <= r_int) f(l1:l2,m,n,ilnrho)=lnrho_int
        elseif (initlnrho(1)=='geo-kws-constant-T'.or.initlnrho(1)=='geo-benchmark') then
          call potential(R=r_int,POT=pot_int)
          call potential(R=r_ext,POT=pot_ext)
          call potential(RMN=r_mn,POT=pot)
          where (r_mn >= r_ext) f(l1:l2,m,n,ilnrho)=lnrho_ext+(pot_ext-pot)*exp(-lnrho_ext/mpoly)
          where (r_mn <= r_int) f(l1:l2,m,n,ilnrho)=lnrho_int+(pot_int-pot)*exp(-lnrho_int/mpoly)
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
      use Global, only: set_global
      use IO

      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: lnrho,cs2
      real, dimension (nx,3) :: glnrho
      real, dimension (nx,3) :: gg_mn
      integer :: j
      
      do m=m1,m2
      do n=n1,n2

        lnrho=f(l1:l2,m,n,ilnrho)
        cs2=cs20*exp(gamma1*(lnrho-lnrho0))
        call grad(f,ilnrho,glnrho)
        do j=1,3
          gg_mn(:,j)=cs2*glnrho(:,j)
        enddo
        call set_global(gg_mn,m,n,'gg',nx)

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
        if (ldensity_nolog) lpenc_requested(i_ugrho)=.true.
        if (.not. ldensity_nolog) lpenc_requested(i_uglnrho)=.true.
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
      if ( ldiff_normal .and. .not. ldensity_nolog) then
        lpenc_requested(i_glnrho2)=.true.
        lpenc_requested(i_del2lnrho)=.true.
      endif
      if (ldiff_normal .and. ldensity_nolog) lpenc_requested(i_del2rho)=.true.
      if (ldiff_hyper3) lpenc_requested(i_del6rho)=.true.
      if (ldiff_hyper3 .and. .not. ldensity_nolog) lpenc_requested(i_rho)=.true.
      if (ldiff_hyper3lnrho) lpenc_requested(i_del6lnrho)=.true.
!
      lpenc_diagnos2d(i_lnrho)=.true.
      lpenc_diagnos2d(i_rho)=.true.
!
      if (idiag_rhom/=0 .or. idiag_rhomz/=0 .or. idiag_rhomy/=0 .or. &
           idiag_rhomx/=0 .or. idiag_rho2m/=0 .or. idiag_rhomin/=0 .or. &
           idiag_rhomax/=0) lpenc_diagnos(i_rho)=.true.
      if (idiag_lnrho2m/=0) lpenc_diagnos(i_lnrho)=.true.
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
      use Global, only: set_global,global_derivs
      use Sub
!      
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (pencil_case) :: p
!      
      integer :: i, mm, nn
!
      intent(in) :: f
      intent(inout) :: p
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
          call dot_mn(p%uu,p%glnrho,p%uglnrho)
        else
          call u_dot_gradf(f,ilnrho,p%glnrho,p%uu,p%uglnrho,UPWIND=lupw_lnrho)
        endif
      endif
! ugrho
      if (lpencil(i_ugrho)) call dot_mn(p%uu,p%grho,p%ugrho)
! glnrho2
      if (lpencil(i_glnrho2)) call dot2_mn(p%glnrho,p%glnrho2)
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
          if (lfirstpoint .and. lglobal_nolog_density) then
            do mm=1,my; do nn=1,mz
              call set_global(exp(f(:,mm,nn,ilnrho)),mm,nn,'rho',mx)
            enddo; enddo
          endif
          if (lglobal_nolog_density) call global_derivs(m,n,'rho',der6=p%del6rho) 
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
      if (lpencil(i_sglnrho)) call multmv_mn(p%sij,p%glnrho,p%sglnrho)
! uij5glnrho
      if (lpencil(i_uij5glnrho)) call multmv_mn(p%uij5,p%glnrho,p%uij5glnrho)
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
      use Special, only: special_calc_density
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
      logical :: ljunk=.true.
!      
      real, dimension (nx) :: fdiff, gshockglnrho, gshockgrho
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
      if (lmass_source) call mass_source(f,df)
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
      if (lspecial) call special_calc_density(f,df,p%uu,p%glnrho,p%divu,p%lnrho)
!
!  phi-averages
!  Note that this does not necessarily happen with ldiagnos=.true.
!
      if (l2davgfirst) then
        call phisum_mn_name_rz(p%lnrho,idiag_lnrhomphi)
        call phisum_mn_name_rz(p%rho,idiag_rhomphi)
      endif
!       
!  Calculate density diagnostics
!
      if (ldiagnos) then
        if (idiag_rhom/=0)    call sum_mn_name(p%rho,idiag_rhom)
        if (idiag_rhomin/=0) &
            call max_mn_name(-p%rho,idiag_rhomin,lneg=.true.)      
        if (idiag_rhomax/=0)  call max_mn_name(p%rho,idiag_rhomax)
        if (idiag_rho2m/=0)   call sum_mn_name(p%rho**2,idiag_rho2m)
        if (idiag_lnrho2m/=0) call sum_mn_name(p%lnrho**2,idiag_lnrho2m)
        if (idiag_rhomz/=0)   call xysum_mn_name_z(p%rho,idiag_rhomz)
        if (idiag_rhomx/=0)   call yzsum_mn_name_x(p%rho,idiag_rhomx)
        if (idiag_rhomy/=0)   call xzsum_mn_name_y(p%rho,idiag_rhomy)
        if (idiag_dtd/=0) &
            call max_mn_name(diffus_diffrho/cdtv,idiag_dtd,l_dt=.true.)
      endif
!
    endsubroutine dlnrho_dt
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
      integer :: iname,inamez,inamey,inamex,irz
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
        idiag_rhom=0; idiag_rho2m=0; idiag_lnrho2m=0
        idiag_rhomin=0; idiag_rhomax=0; idiag_dtd=0
        idiag_lnrhomphi=0; idiag_rhomphi=0
        idiag_rhomz=0; idiag_rhomy=0; idiag_rhomx=0
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
        call parse_name(iname,cname(iname),cform(iname),'dtd',idiag_dtd)
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
!  check for those quantities for which we want phi-averages
!
      do irz=1,nnamerz
        call parse_name(irz,cnamerz(irz),cformrz(irz),&
            'lnrhomphi',idiag_lnrhomphi)
        call parse_name(irz,cnamerz(irz),cformrz(irz),'rhomphi',idiag_rhomphi)
      enddo
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_rhom=',idiag_rhom
        write(3,*) 'i_rho2m=',idiag_rho2m
        write(3,*) 'i_rhomin=',idiag_rhomin
        write(3,*) 'i_rhomax=',idiag_rhomax
        write(3,*) 'i_lnrho2m=',idiag_lnrho2m
        write(3,*) 'i_rhomz=',idiag_rhomz
        write(3,*) 'i_rhomy=',idiag_rhomy
        write(3,*) 'i_rhomx=',idiag_rhomx
        write(3,*) 'nname=',nname
        write(3,*) 'ilnrho=',ilnrho
        write(3,*) 'i_lnrhomphi=',idiag_lnrhomphi
        write(3,*) 'i_rhomphi=',idiag_rhomphi
        write(3,*) 'i_dtd=',idiag_dtd
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: pot,tmp
!
!  Stratification depends on the gravity potential
!
      if (lroot) print*,'isothermal_density: isothermal stratification'
      if ( (.not. lentropy) .and. (gamma/=1.0) ) then
        call fatal_error('isothermal_density','for gamma/=1.0, you need entropy!');
      endif
      do n=n1,n2
        do m=m1,m2
          call potential(x(l1:l2),y(m),z(n),pot=pot)
          tmp=-gamma*pot/cs20
          f(l1:l2,m,n,ilnrho) = f(l1:l2,m,n,ilnrho) + lnrho0 + tmp

          if (lentropy) f(l1:l2,m,n,iss) = f(l1:l2,m,n,iss) &
               -gamma1*(f(l1:l2,m,n,ilnrho)-lnrho0)/gamma
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
                             potential
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (nx) :: pot,dlncs2,ptop,pbot,zero=0.
      real :: ggamma,ztop,zbot,zinfty2,pot_ext,lnrho_ref
!
!  identifier
!
      if (lroot) print*,'polytropic_simple: mpoly=',mpoly
!
!  The following is specific only to cases with gravity in the z direction
!  zinfty is calculated such that rho=rho0 and cs2=cs20 at z=zref.
!  Note: gravz is normally negative!
!
      if (lgravz) then
        if (grav_profile=='const') then
          if(lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zinfty=zref+(mpoly+1.)*cs20/(-gamma*gravz)
        elseif (grav_profile=='const_zero') then
          if(lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zinfty=zref+(mpoly+1.)*cs20/(-gamma*gravz)
        elseif (grav_profile=='linear') then
          if(lroot.and.gravz==0.) print*,'polytropic_simple: divide by gravz=0'
          zinfty2=zref**2+(mpoly+1.)*cs20/(-.5*gamma*gravz)
          if(zinfty2<0) then
            if(lroot) print*,'polytropic_simple: zinfty**2<0 is not ok'
            zinfty2=0. !(and see what happens)
          endif
          zinfty=sqrt(zinfty2)
        else
          if(lroot) print*,'polytropic_simple: zinfty not prepared!'
        endif
!
!  check whether zinfty lies outside the domain (otherwise density
!  would vanish within the domain). At the moment we are not properly
!  testing the lower boundary on the case of a disk (commented out).
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
        enddo
        enddo
!
!  cs2 values at top and bottom may be needed to boundary conditions.
!  In spherical geometry, ztop is z at the outer edge of the box,
!  so this calculation still makes sense.
!
        call potential(zero,0.,ztop,pot=ptop)
        cs2top=-gamma/(mpoly+1.)*ptop(1)
!
!  In spherical geometry ztop should never be used.
!  Even in slab geometry ztop is not normally used.
!
        call potential(zero,0.,zbot,pot=pbot)
        cs2bot=-gamma/(mpoly+1.)*pbot(1)
      endif
!
    endsubroutine polytropic_simple
!***********************************************************************
    subroutine mass_source(f,df)
!
!  add mass sources and sinks
!
!  28-apr-2005/axel: coded
!
      use Cdata
      use Sub, only: step
      use Gravity, only: lnrho_bot,lnrho_top,ss_bot,ss_top
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension(nx) :: fint,fext,pdamp
!
      if(ldebug) print*,'mass_source: cs20,cs0=',cs20,cs0
!
!  cylindrical profile for inner cylinder
!
      pdamp=1-step(rcyl_mn,r_int,wdamp) ! inner damping profile
      fint=-damplnrho_int*pdamp*(f(l1:l2,m,n,ilnrho)-lnrho_int)
      df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+fint
!
!  cylindrical profile for outer cylinder
!
      pdamp=step(rcyl_mn,r_ext,wdamp) ! outer damping profile
      fext=-damplnrho_ext*pdamp*(f(l1:l2,m,n,ilnrho)-lnrho_ext)
      df(l1:l2,m,n,ilnrho)=df(l1:l2,m,n,ilnrho)+fext
!     df(l1:l2,m,n,iux)=df(l1:l2,m,n,iux)+f(l1:l2,m,n,iux)*fext
!     df(l1:l2,m,n,iuy)=df(l1:l2,m,n,iuy)+f(l1:l2,m,n,iuy)*fext
!
    endsubroutine mass_source
!***********************************************************************
endmodule Density
