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
! MVAR CONTRIBUTION 4
! MAUX CONTRIBUTION 0
!
! PENCILS PROVIDED infl_phi; infl_dphi; infl_a2; infl_a21; gphi(3)
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
  integer :: iinfl_phi=0, iinfl_dphi=0, iinfl_hubble=0, iinfl_lna=0, Ndiv=100
  real :: ncutoff_phi=1.
  real :: axionmass=1.06e-6, axionmass2, ascale_ini=1.
  real :: phi0=.44, dphi0=-1e-5, c_light_axion=1., lambda_axion=0., eps=.01
  real :: amplphi=.1, ampldphi=.0, kx_phi=1., ky_phi=0., kz_phi=0., phase_phi=0., width=.1, offset=0.
  real :: initpower_phi=0.,  cutoff_phi=0.,  initpower2_phi=0.
  real :: initpower_dphi=0., cutoff_dphi=0., initpower2_dphi=0.
  real :: kgaussian_phi=0.,kpeak_phi=0., kgaussian_dphi=0., kpeak_dphi=0.
  real :: relhel_phi=0.
  real :: ddotam, a2rhopm, a2rhopm_all
  real, target :: ddotam_all
  real, pointer :: alpf
  real, dimension (nx) :: dt1_special
  logical :: lbackreact_infl=.true., lzeroHubble=.false.
  logical :: lscale_tobox=.true.,ldt_backreact_infl=.true.
  logical :: lskip_projection_phi=.false., lvectorpotential=.false.
!
  character (len=labellen) :: Vprime_choice='quadratic'
  character (len=labellen), dimension(ninit) :: initspecial='nothing'
!
  namelist /special_init_pars/ &
      initspecial, phi0, dphi0, axionmass, eps, ascale_ini, &
      c_light_axion, lambda_axion, amplphi, ampldphi, &
      kx_phi, ky_phi, kz_phi, phase_phi, width, offset, &
      initpower_phi, cutoff_phi, kgaussian_phi, kpeak_phi, &
      initpower_dphi, cutoff_dphi, kpeak_dphi, &
      ncutoff_phi, lscale_tobox
!
  namelist /special_run_pars/ &
      initspecial, phi0, dphi0, axionmass, eps, ascale_ini, &
      lbackreact_infl, c_light_axion, lambda_axion, Vprime_choice, &
      lzeroHubble, ldt_backreact_infl, Ndiv
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
      call farray_register_pde('infl_hubble',iinfl_hubble)
      call farray_register_pde('infl_lna',iinfl_lna)
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
      call put_shared_variable('ddotam',ddotam_all,caller='initialize_backreact_infl')
!
      if (lmagnetic .and. lbackreact_infl) call get_shared_variable('alpf',alpf)
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
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: Vpotential, Hubble_ini, lnascale
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
          case ('nophi')
            Vpotential=.5*axionmass2*phi0**2
            dphi0=0.
            tstart=-sqrt(3./(8.*pi))/(ascale_ini*sqrt(Vpotential))
            t=tstart
            Hubble_ini=sqrt(8.*pi/3.*(.5*dphi0**2+.5*axionmass2*phi0**2*ascale_ini**2))
            lnascale=log(ascale_ini)
            f(:,:,:,iinfl_hubble)=f(:,:,:,iinfl_hubble)+Hubble_ini
            f(:,:,:,iinfl_lna)=f(:,:,:,iinfl_lna)+lnascale
          case ('default')
            Vpotential=.5*axionmass2*phi0**2
            dphi0=-ascale_ini*sqrt(2*eps/3.*Vpotential)
            tstart=-sqrt(3./(8.*pi))/(ascale_ini*sqrt(Vpotential))
            t=tstart
            Hubble_ini=sqrt(8.*pi/3.*(.5*dphi0**2+.5*axionmass2*phi0**2*ascale_ini**2))
            lnascale=log(ascale_ini)
            f(:,:,:,iinfl_phi)   =f(:,:,:,iinfl_phi)   +phi0
            f(:,:,:,iinfl_dphi)  =f(:,:,:,iinfl_dphi)  +dphi0
            f(:,:,:,iinfl_hubble)=f(:,:,:,iinfl_hubble)+Hubble_ini
            f(:,:,:,iinfl_lna)   =f(:,:,:,iinfl_lna)   +lnascale
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
              cutoff_phi,ncutoff_phi,kpeak_dphi,f,iinfl_dphi,iinfl_dphi, &
              relhel_phi,kgaussian_phi, lskip_projection_phi, lvectorpotential, &
              lscale_tobox, lpower_profile_file=.false.)
  
          case default
            call fatal_error("init_special: No such initspecial: ", trim(initspecial(j)))
        endselect
      enddo
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (lmagnetic .and. lbackreact_infl) lpenc_requested(i_infl_a21)=.true.
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
! infl_a2
      if (lpencil(i_infl_a2)) p%infl_a2=exp(2.*f(l1:l2,m,n,iinfl_lna))
!
! infl_a21
      if (lpencil(i_infl_a21)) then
        if (lpencil(i_infl_a2)) then
          p%infl_a21=1./p%infl_a2
        else
          p%infl_a21=exp(-2.*f(l1:l2,m,n,iinfl_lna))
        endif
      endif
!
! infl_gphi
      if (lpencil(i_gphi)) call grad(f,iinfl_phi,p%gphi)
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
!   2-nov-21/axel: first set of equations coded
!
      use Diagnostics, only: sum_mn_name, max_mn_name, save_name
      use Sub, only: dot_mn, del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: phi, dphi, Hscript, lnascale, ascale, a2scale, a2rhop, Vprime
      real, dimension (nx) :: tmp, del2phi
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
      Hscript=f(l1:l2,m,n,iinfl_hubble)
      lnascale=f(l1:l2,m,n,iinfl_lna)
      ascale=exp(lnascale)
      a2scale=ascale**2
!     a2rhop=dphi**2
!     a2rhopm=<dphi**2+gphi**2+(4./3.)*a^*(E^2+B^2)
!
!  Possibility of turning off evolution of scale factor and Hubble parameter
!  By default, lzeroHubble=F, so we use the calculation from above.
!
      if (lzeroHubble) then
        a2scale=1.
        Hscript=0.
      endif
!
!  Choice of different potentials
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
        df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)-2.*Hscript*dphi-a2scale*Vprime
        df(l1:l2,m,n,iinfl_hubble)=df(l1:l2,m,n,iinfl_hubble)-4.*pi*a2rhopm_all+Hscript**2
        df(l1:l2,m,n,iinfl_lna)=df(l1:l2,m,n,iinfl_lna)+Hscript
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
      if (lmagnetic .and. lbackreact_infl) then
        call dot_mn(p%el,p%bb,tmp)
        df(l1:l2,m,n,iinfl_dphi)=df(l1:l2,m,n,iinfl_dphi)+alpf*p%infl_a21*tmp
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
        call sum_mn_name(Hscript,idiag_Hubblem)
        call sum_mn_name(lnascale,idiag_lnam)
        call save_name(ddotam_all,idiag_ddotam)
      endif

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
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_phim=0; idiag_phi2m=0; idiag_phirms=0
        idiag_dphim=0; idiag_dphi2m=0; idiag_dphirms=0
        idiag_Hubblem=0; idiag_lnam=0; idiag_ddotam=0
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
      use Mpicomm, only: mpiallreduce_sum
      use Sub, only: dot2_mn, grad, curl
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (nx,3) :: el, bb, gphi
      real, dimension (nx) :: e2, b2, gphi2, dphi, a21, a2rhop
      real, dimension (nx) :: ddota, phi, a2, Vpotential
!
!  if requested, calculate here <dphi**2+gphi**2+(4./3.)*(E^2+B^2)/a^2>
!
      ddotam=0.; a2rhopm=0.

      do n=n1,n2
      do m=m1,m2
        
        phi=f(l1:l2,m,n,iinfl_phi)
        dphi=f(l1:l2,m,n,iinfl_dphi)
        call grad(f,iinfl_phi,gphi)    !MR: the ghost zones are not necessarily updated!!!
        call dot2_mn(gphi,gphi2)
        a2rhop=dphi**2+gphi2

        if (iex/=0.and.lbackreact_infl) then
          a2=exp(2.*f(l1:l2,m,n,iinfl_lna))
          a21=1./a2
          el=f(l1:l2,m,n,iex:iez)
          call curl(f,iaa,bb)          !MR: the ghost zones are not necessarily updated!!!
          call dot2_mn(bb,b2)
          call dot2_mn(el,e2)
          a2rhop=a2rhop+(.5*fourthird)*(e2+b2)*a21
        endif

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
      enddo
      enddo
     
      a2rhopm=a2rhopm/nwgrid
      ddotam=(four_pi_over_three/nwgrid)*ddotam

      call mpiallreduce_sum(a2rhopm,a2rhopm_all)
      call mpiallreduce_sum(ddotam,ddotam_all)
!
    endsubroutine special_after_boundary
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
