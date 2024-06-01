! $Id$
!
!  This module provide a way for users to specify custom
!  (i.e. not in the standard Pencil Code) physics, diagnostics etc.
!
!  The module provides a set of standard hooks into the Pencil-Code and
!  currently allows the following customizations:
!
!  Description                                     | Relevant function call
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
!  Special term in the mass (density) equation     | special_calc_density
!  Special term in the momentum (hydro) equation   | special_calc_hydro
!  Special term in the energy equation             | special_calc_energy
!  Special term in the induction (magnetic)        | special_calc_magnetic
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
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED infl_phi; infl_dphi; gphi(3)
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
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
!
! Declare index of new variables in f array (if any).
!
!  integer :: iinfl_phi=0, iinfl_dphi=0, iinfl_hubble=0, iinfl_lna=0, Ndiv=100
  integer :: iinfl_phi=0, iinfl_dphi=0, iinfl_lna=0, Ndiv=100
  real :: ncutoff_phi=1., infl_v=.1
  real :: axionmass=1.06e-6, axionmass2, ascale_ini=1.
  real :: phi0=.44, dphi0=-1e-5, c_light_axion=1., lambda_axion=0., eps=.01
  real :: amplphi=.1, ampldphi=.0, kx_phi=1., ky_phi=0., kz_phi=0., phase_phi=0., width=.1, offset=0.
  real :: initpower_phi=0.,  cutoff_phi=0.,  initpower2_phi=0.
  real :: initpower_dphi=0., cutoff_dphi=0., initpower2_dphi=0.
  real :: kgaussian_phi=0.,kpeak_phi=0., kgaussian_dphi=0., kpeak_dphi=0.
  real :: relhel_phi=0.
  real :: ddotam, a2rhopm, a2rhopm_all, a2rhom, a2rhom_all, edotbm, edotbm_all, a2rhophim, a2rhophim_all
  real :: lnascale, ascale, a2, a21, Hscript
  real :: Hscript0=0.
  real, target :: ddotam_all
  real, pointer :: alpf
  real, dimension (nx) :: dt1_special
  logical :: lbackreact_infl=.true., lem_backreact=.true., lzeroHubble=.false.
  logical :: lscale_tobox=.true.,ldt_backreact_infl=.true.
  logical :: lskip_projection_phi=.false., lvectorpotential=.false., lflrw=.false.
  logical, pointer :: lphi_hom
!
  character (len=labellen) :: Vprime_choice='quadratic', Hscript_choice='default'
  character (len=labellen), dimension(ninit) :: initspecial='nothing'
!
  namelist /special_init_pars/ &
      initspecial, phi0, dphi0, axionmass, eps, ascale_ini, &
      c_light_axion, lambda_axion, amplphi, ampldphi, &
      kx_phi, ky_phi, kz_phi, phase_phi, width, offset, &
      initpower_phi, initpower2_phi, cutoff_phi, kgaussian_phi, kpeak_phi, &
      initpower_dphi, initpower2_dphi, cutoff_dphi, kpeak_dphi, &
      ncutoff_phi, lscale_tobox, Hscript0, Hscript_choice, infl_v, lflrw
!
  namelist /special_run_pars/ &
      initspecial, phi0, dphi0, axionmass, eps, ascale_ini, &
      lbackreact_infl, lem_backreact, c_light_axion, lambda_axion, Vprime_choice, &
      lzeroHubble, ldt_backreact_infl, Ndiv, Hscript0, Hscript_choice, infl_v, &
      lflrw
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_phim=0      ! DIAG_DOC: $\left<\phi\right>$
  integer :: idiag_phi2m=0     ! DIAG_DOC: $\left<\phi^2\right>$
  integer :: idiag_phirms=0    ! DIAG_DOC: $\left<\phi^2\right>^{1/2}$
  integer :: idiag_dphim=0     ! DIAG_DOC: $\left<\phi'\right>$
  integer :: idiag_dphi2m=0    ! DIAG_DOC: $\left<(\phi')^2\right>$
  integer :: idiag_dphirms=0   ! DIAG_DOC: $\left<(\phi')^2\right>^{1/2}$
  integer :: idiag_Hubblem=0   ! DIAG_DOC: $\left<{\cal H}\right>$
  integer :: idiag_lnam=0      ! DIAG_DOC: $\left<\ln a\right>$
  integer :: idiag_ddotam=0    ! DIAG_DOC: $a''/a$
  integer :: idiag_a2rhopm=0   ! DIAG_DOC: $a^2 (rho+p)$
  integer :: idiag_a2rhom=0   ! DIAG_DOC: $a^2 rho$
  integer :: idiag_a2rhophim=0   ! DIAG_DOC: $a^2 rho$
!
  contains
!****************************************************************************
    subroutine register_special
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
  use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_pde('infl_phi',iinfl_phi)
      call farray_register_pde('infl_dphi',iinfl_dphi)
!
     if (lflrw) then
!     call farray_register_ode('infl_hubble',iinfl_hubble)
       call farray_register_ode('infl_lna',iinfl_lna)
     endif
!
!  for power spectra, it is convenient to use ispecialvar and
!
      ispecialvar=iinfl_phi
      ispecialvar2=iinfl_dphi
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      use SharedVariables, only: get_shared_variable, put_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      axionmass2=axionmass**2
!
      call put_shared_variable('ddotam',ddotam_all,caller='initialize_backreact_infl_ode')
!
      if (lmagnetic .and. lem_backreact) then
        call get_shared_variable('alpf',alpf)
        call get_shared_variable('lphi_hom',lphi_hom)
      else
        allocate(alpf)
        allocate(lphi_hom)
        alpf=0.
        lphi_hom=.false.
      endif
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Initcond, only: gaunoise, sinwave_phase, hat, power_randomphase_hel, power_randomphase
      use Mpicomm, only: mpibcast_real
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: Vpotential, Hubble_ini, infl_gam
      integer :: j
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      do j=1,ninit

        select case (initspecial(j))
          case ('nothing'); if (lroot) print*,'init_special: nothing'
          case ('phi=sinkx')
            f(:,:,:,iinfl_phi)=f(:,:,:,iinfl_phi) &
              +spread(spread(amplphi*sin(kx_phi*x),2,my),3,mz)
          case ('phi=tanhkx')
            f(:,:,:,iinfl_phi)=f(:,:,:,iinfl_phi) &
              +spread(spread(.5*amplphi*(1.+tanh(kx_phi*(x-offset))),2,my),3,mz)
          case ('phi=atan_exp_kx')
            infl_gam=1./sqrt(1.-infl_v**2)
            f(:,:,:,iinfl_phi)=f(:,:,:,iinfl_phi) &
              +spread(spread(4.*amplphi*atan(exp(infl_gam*kx_phi*(x-offset))),2,my),3,mz)
            f(:,:,:,iinfl_dphi)=f(:,:,:,iinfl_dphi)+spread(spread( &
              -4.*amplphi*kx_phi*infl_gam*infl_v*exp(infl_gam*kx_phi*(x-offset)) &
              /(exp(2.*infl_gam*kx_phi*(x-offset))+1.) &
              ,2,my),3,mz)
          case ('nophi')
            Vpotential=.5*axionmass2*phi0**2
            dphi0=0.
            tstart=-sqrt(3./(8.*pi))/(ascale_ini*sqrt(Vpotential))
            t=tstart
            Hubble_ini=sqrt(8.*pi/3.*(.5*dphi0**2+.5*axionmass2*phi0**2*ascale_ini**2))
            lnascale=log(ascale_ini)
            if (lroot .and. lflrw) then
              f_ode(iinfl_lna) =lnascale
!              f(iinfl_hubble) =Hubble_ini
            endif
          case ('default')
            Vpotential=.5*axionmass2*phi0**2
            dphi0=-ascale_ini*sqrt(2*eps/3.*Vpotential)
            tstart=-sqrt(3./(8.*pi))/(ascale_ini*sqrt(Vpotential))
            t=tstart
            Hubble_ini=sqrt(8.*pi/3.*(.5*dphi0**2+.5*axionmass2*phi0**2*ascale_ini**2))
            lnascale=log(ascale_ini)
            f(:,:,:,iinfl_phi)   =f(:,:,:,iinfl_phi)   +phi0
            f(:,:,:,iinfl_dphi)  =f(:,:,:,iinfl_dphi)  +dphi0
            if (lroot .and. lflrw) then
              f_ode(iinfl_lna)   =lnascale
              a2                 =exp(f_ode(iinfl_lna))**2
              Hscript            =Hubble_ini/exp(lnascale)
!              f(iinfl_hubble)   =Hscript
            endif
          case ('gaussian-noise')
            call gaunoise(amplphi,f,iinfl_phi)
          case ('sinwave-phase')
            !call sinwave_phase(f,iinfl_phi,amplphi,kx_phi,ky_phi,kz_phi,phase_phi)
            !f(:,:,:,iinfl_phi)=tanh(f(:,:,:,iinfl_phi)/width)
            call hat(amplphi,f,iinfl_phi,width,kx_phi,ky_phi,kz_phi)
            f(:,:,:,iinfl_phi)=f(:,:,:,iinfl_phi)+offset
          case ('phi_power_randomphase')
            call power_randomphase_hel(amplphi,initpower_phi,initpower2_phi, &
              cutoff_phi,ncutoff_phi,kpeak_phi,f,iinfl_phi,iinfl_phi, &
              relhel_phi,kgaussian_phi, lskip_projection_phi, lvectorpotential, &
              lscale_tobox, lpower_profile_file=.false.)
          case ('dphi_power_randomphase')
            call power_randomphase_hel(ampldphi,initpower_dphi,initpower2_dphi, &
              cutoff_dphi,ncutoff_phi,kpeak_dphi,f,iinfl_dphi,iinfl_dphi, &
              relhel_phi,kgaussian_phi, lskip_projection_phi, lvectorpotential, &
              lscale_tobox, lpower_profile_file=.false.)
  
          case default
            call fatal_error("init_special: No such initspecial: ", trim(initspecial(j)))
        endselect
      enddo
      call mpibcast_real(a2)
      call mpibcast_real(Hscript)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
!      if (lmagnetic .and. lbackreact_infl) lpenc_requested(i_infl_a21)=.true.
!
!  pencil for gradient of phi
!
      lpenc_requested(i_gphi)=.true.
!
!  Magnetic field needed for Maxwell stress
!
      if (lmagnetic) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_el)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Special pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!  24-nov-04/tony: coded
!
      use Sub, only: grad
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
! infl_phi
      if (lpencil(i_infl_phi)) p%infl_phi=f(l1:l2,m,n,iinfl_phi)
!
! infl_dphi
      if (lpencil(i_infl_dphi)) p%infl_dphi=f(l1:l2,m,n,iinfl_dphi)
!
! infl_gphi
      if (lpencil(i_gphi)) call grad(f,iinfl_phi,p%gphi)
!
    endsubroutine calc_pencils_special
!***********************************************************************
    subroutine dspecial_dt(f,df,p)
!
!  The entire module could be renamed to Klein-Gordon or Scalar field equation.
!  Calculate right hand side of ONE OR MORE extra coupled PDEs
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
!   2-nov-21/axel: first set of equations coded
!
      use Diagnostics, only: sum_mn_name, max_mn_name, save_name
      use Sub, only: dot_mn, del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: phi, dphi, Vprime
      real, dimension (nx) :: tmp, del2phi
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
      phi=f(l1:l2,m,n,iinfl_phi)
      dphi=f(l1:l2,m,n,iinfl_dphi)
!
!  Choice of prescription for Hscript
!
      select case (Hscript_choice)
        case ('default')
          Hscript=sqrt((8.*pi/3.)*a2rhom_all)
        case ('set')
          Hscript=Hscript0
          a2=1.
          a21=1./a2
        case default
          call fatal_error("dspecial_dt: No such Hscript_choice: ", trim(Hscript_choice))
      endselect

!  Possibility of turning off evolution of scale factor and Hubble parameter
!  By default, lzeroHubble=F, so we use the calculation from above.
!
      if (lzeroHubble) then
        a2=1.
        a21=1./a2
        Hscript=0.
      endif
!
!  Choice of different potentials.
!  For the 1-cos profile, -Vprime (on the rhs) enters with -sin().
!
      select case (Vprime_choice)
        case ('quadratic'); Vprime=axionmass2*phi
        case ('quartic'); Vprime=axionmass2*phi+(lambda_axion/6.)*phi**3
        case ('cos-profile'); Vprime=axionmass2*lambda_axion*sin(lambda_axion*phi)
        case default
          call fatal_error("dspecial_dt: No such Vprime_choice: ", trim(Vprime_choice))
      endselect
!
!  Update df.
!  dphi/dt = psi
!  dpsi/dt = - ...
!
        df(l1:l2,m,n,iinfl_phi)=df(l1:l2,m,n,iinfl_phi)+f(l1:l2,m,n,iinfl_dphi)
        df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)-2.*Hscript*dphi-a2*Vprime
!
!  speed of light term
!
        if (c_light_axion/=0.) then
          call del2(f,iinfl_phi,del2phi)
          df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+c_light_axion**2*del2phi
        endif
!
!  magnetic terms, add (alpf/a^2)*(E.B) to dphi'/dt equation
!
      if (lmagnetic .and. lem_backreact) then
        if (lphi_hom) then
          df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*edotbm_all*a21
        else
          call dot_mn(p%el,p%bb,tmp)
          df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*tmp*a21
        endif
      endif
!
!  Total contribution to the timestep
!
      if (lfirst.and.ldt.and.ldt_backreact_infl) then
        dt1_special = Ndiv*abs(Hscript)
        dt1_max=max(dt1_max,dt1_special) 
      endif
!
!  Diagnostics
!
      if (ldiagnos) then
        call sum_mn_name(phi,idiag_phim)
        if (idiag_phi2m/=0) call sum_mn_name(phi**2,idiag_phi2m)
        if (idiag_phirms/=0) call sum_mn_name(phi**2,idiag_phirms,lsqrt=.true.)
        call sum_mn_name(dphi,idiag_dphim)
        if (idiag_dphi2m/=0) call sum_mn_name(dphi**2,idiag_dphi2m)
        if (idiag_dphirms/=0) call sum_mn_name(dphi**2,idiag_dphirms,lsqrt=.true.)
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine dspecial_dt_ode
!
      use Diagnostics, only: save_name
!
      if (lflrw) then
        df_ode(iinfl_lna)=df_ode(iinfl_lna)+Hscript
      endif
!
!  Diagnostics
!
      if (ldiagnos) then
        call save_name(Hscript,idiag_Hubblem)
        call save_name(lnascale,idiag_lnam)
        call save_name(ddotam_all,idiag_ddotam)
        call save_name(a2rhopm_all,idiag_a2rhopm)
        call save_name(a2rhom_all,idiag_a2rhom)
        call save_name(a2rhophim_all,idiag_a2rhophim)
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
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_phim=0; idiag_phi2m=0; idiag_phirms=0
        idiag_dphim=0; idiag_dphi2m=0; idiag_dphirms=0
        idiag_Hubblem=0; idiag_lnam=0; idiag_ddotam=0
        idiag_a2rhopm=0; idiag_a2rhom=0; idiag_a2rhophim=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'phim',idiag_phim)
        call parse_name(iname,cname(iname),cform(iname),'phi2m',idiag_phi2m)
        call parse_name(iname,cname(iname),cform(iname),'phirms',idiag_phirms)
        call parse_name(iname,cname(iname),cform(iname),'dphim',idiag_dphim)
        call parse_name(iname,cname(iname),cform(iname),'dphi2m',idiag_dphi2m)
        call parse_name(iname,cname(iname),cform(iname),'dphirms',idiag_dphirms)
        call parse_name(iname,cname(iname),cform(iname),'Hubblem',idiag_Hubblem)
        call parse_name(iname,cname(iname),cform(iname),'lnam',idiag_lnam)
        call parse_name(iname,cname(iname),cform(iname),'ddotam',idiag_ddotam)
        call parse_name(iname,cname(iname),cform(iname),'a2rhopm',idiag_a2rhopm)
        call parse_name(iname,cname(iname),cform(iname),'a2rhom',idiag_a2rhom)
        call parse_name(iname,cname(iname),cform(iname),'a2rhophim',idiag_a2rhophim)
      enddo
!!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        call farray_index_append('idiag_SPECIAL_DIAGNOSTIC',idiag_SPECIAL_DIAGNOSTIC)
!!      endif
!!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      use Mpicomm, only: mpireduce_sum, mpiallreduce_sum, mpibcast_real
      use Sub, only: dot2_mn, grad, curl, dot_mn
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!      real, dimension (nx,3) :: el, bb, gphi
!      real, dimension (nx) :: e2, b2, gphi2, dphi, a21, a2rhop, a2rho
!      real, dimension (nx) :: ddota, phi, a2, Vpotential, edotb
!
!  if requested, calculate here <dphi**2+gphi**2+(4./3.)*(E^2+B^2)/a^2>
!
      if (lroot) then
        if(lflrw) then
          lnascale=f_ode(iinfl_lna)
          a2=exp(2*lnascale)
        else
          a2=1.
        endif
      endif
      a21=1./a2
      call mpibcast_real(a2)
      call mpibcast_real(a21)
!
      ddotam=0.; a2rhopm=0.; a2rhom=0.; edotbm=0; a2rhophim=0.;
      do n=n1,n2
      do m=m1,m2
        call prep_ode_right(f)
      enddo
      enddo
!
      a2rhopm=a2rhopm/nwgrid
      a2rhom=a2rhom/nwgrid
      a2rhophim=a2rhophim/nwgrid
      ddotam=(four_pi_over_three/nwgrid)*ddotam
      if (lphi_hom) then
!          edotbm=edotbm/nwgrid
          call mpireduce_sum(edotbm,edotbm_all)
      endif
!
      call mpireduce_sum(a2rhopm,a2rhopm_all)
      call mpiallreduce_sum(a2rhom,a2rhom_all)
      call mpireduce_sum(a2rhophim,a2rhophim_all)
      call mpireduce_sum(ddotam,ddotam_all)
!
      if (lroot .and. lflrw) then
        Hscript=(8.*pi/3.)*sqrt(a2rhom_all)
      endif
      call mpibcast_real(Hscript)
!
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine prep_ode_right(f)
!
      use Sub, only: dot2_mn, grad, curl, dot_mn
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx,3) :: el, bb, gphi
      real, dimension (nx) :: e2, b2, gphi2, dphi, a2rhop, a2rho
      real, dimension (nx) :: ddota, phi, Vpotential, edotb
!
!  if requested, calculate here <dphi**2+gphi**2+(4./3.)*(E^2+B^2)/a^2>
!
      phi=f(l1:l2,m,n,iinfl_phi)
      dphi=f(l1:l2,m,n,iinfl_dphi)
      call grad(f,iinfl_phi,gphi)    !MR: the ghost zones are not necessarily updated!!!
      call dot2_mn(gphi,gphi2)
      a2rhop=dphi**2+gphi2
      a2rho=0.5*(dphi**2+gphi2)
      a2rhophim=a2rhophim+sum(a2rho)

!
      if (iex/=0 .and. lem_backreact) then
        el=f(l1:l2,m,n,iex:iez)
        call curl(f,iaa,bb)          !MR: the ghost zones are not necessarily updated!!!
        call dot2_mn(bb,b2)
        call dot2_mn(el,e2)
        a2rhop=a2rhop+(.5*fourthird)*(e2+b2)*a21
        a2rho=a2rho+.5*(e2+b2)*a21
      endif
!
      a2rhopm=a2rhopm+sum(a2rhop)
!
!  Choice of different potentials
!
      select case (Vprime_choice)
        case ('quadratic')  ; Vpotential=.5*axionmass2*phi**2
        case ('quartic')    ; Vpotential=axionmass2*phi+(lambda_axion/6.)*phi**3  !(to be corrected)
        case ('cos-profile'); Vpotential=axionmass2*lambda_axion*sin(lambda_axion*phi)  !(to be corrected)
        case default
          call fatal_error("special_after_boundary: No such Vprime_choice: ",trim(Vprime_choice))
      endselect
!
!  compute ddotam = a"/a (needed for GW module)
!
      ddota=-dphi**2-gphi2+4.*a2*Vpotential
      ddotam=ddotam+sum(ddota)
      a2rho=a2rho+a2*Vpotential
      a2rhom=a2rhom+sum(a2rho)
      if (lmagnetic .and. lem_backreact .and. lphi_hom) then
        call dot_mn(el,bb,edotb)
        edotbm=edotbm+sum(edotb)
      endif
!
    endsubroutine prep_ode_right
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
