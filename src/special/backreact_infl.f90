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
! PENCILS PROVIDED infl_phi; infl_dphi; infl_a2; infl_a21
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
!
! Declare index of new variables in f array (if any).
!
  integer :: ispecial=0
  real :: axionmass=1.06e-6, axionmass2, ascale_ini=1.
  real :: phi0=.44, dphi0=-1e-5, c_light_axion=0., lambda_axion=0.
  real :: amplphi=.1, kx_phi=1., ky_phi=0., kz_phi=0., phase_phi=0., width=.1, offset=0.
  real, pointer :: alpf
  logical :: lbackreact_infl=.true., lzeroHubble=.false.
!
  character (len=labellen) :: initspecial='nothing', Vprime_choice='quadratic'
!
  namelist /special_init_pars/ &
      initspecial, phi0, dphi0, axionmass, ascale_ini, &
      c_light_axion, lambda_axion, amplphi, &
      kx_phi, ky_phi, kz_phi, phase_phi, width, offset
!
  namelist /special_run_pars/ &
      initspecial, phi0, dphi0, axionmass, ascale_ini, &
      lbackreact_infl, c_light_axion, lambda_axion, Vprime_choice, &
      lzeroHubble
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_phim=0      ! DIAG_DOC: $\left<\phi\right>$
  integer :: idiag_dphim=0      ! DIAG_DOC: $\left<d\phi\right>$
  integer :: idiag_Hubblem=0      ! DIAG_DOC: $\left<{\cal H}\right>$
  integer :: idiag_lnam=0      ! DIAG_DOC: $\left<\ln a\right>$
!
  contains
!****************************************************************************
  subroutine initialize_mult_special
!
! Dummy routine.
!
  endsubroutine initialize_mult_special
!***********************************************************************
  subroutine finalize_mult_special
!
! Dummy routine.
!
  endsubroutine finalize_mult_special
!***********************************************************************
    subroutine register_special
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
! use Sub, only: farray_register_pde
! use Sub, only: register_report_aux
  use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_pde('special',ispecial,array=4)
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
      use SharedVariables, only: get_shared_variable
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      axionmass2=axionmass**2
!
      if (lmagnetic .and. lbackreact_infl) then
        call get_shared_variable('alpf', alpf, caller='initialize_backreact_infl')
      endif
!
      call keep_compiler_quiet(f)
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
      use Initcond, only: gaunoise, sinwave_phase, hat
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: Vpotential, eps=.01, Hubble_ini, lnascale
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      select case (initspecial)
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('nophi')
          Vpotential=.5*axionmass2*phi0**2
          dphi0=0.
          tstart=-sqrt(3./(8.*pi))/(ascale_ini*sqrt(Vpotential))
          t=tstart
          Hubble_ini=sqrt(8.*pi/3.*(.5*dphi0**2+.5*axionmass2*phi0**2*ascale_ini**2))
          lnascale=log(ascale_ini)
          f(:,:,:,ispecial+0)=0.
          f(:,:,:,ispecial+1)=0.
          f(:,:,:,ispecial+2)=Hubble_ini
          f(:,:,:,ispecial+3)=lnascale
        case ('default')
          Vpotential=.5*axionmass2*phi0**2
          dphi0=-ascale_ini*sqrt(2*eps/3.*Vpotential)
          tstart=-sqrt(3./(8.*pi))/(ascale_ini*sqrt(Vpotential))
          t=tstart
          Hubble_ini=sqrt(8.*pi/3.*(.5*dphi0**2+.5*axionmass2*phi0**2*ascale_ini**2))
          lnascale=log(ascale_ini)
          f(:,:,:,ispecial+0)=phi0
          f(:,:,:,ispecial+1)=dphi0
          f(:,:,:,ispecial+2)=Hubble_ini
          f(:,:,:,ispecial+3)=lnascale
        case ('gaussian-noise')
          call gaunoise(amplphi,f,ispecial+0)
          f(:,:,:,ispecial+1)=0.
          f(:,:,:,ispecial+2)=0.
          f(:,:,:,ispecial+3)=0.
        case ('sinwave-phase')
          !call sinwave_phase(f,ispecial+0,amplphi,kx_phi,ky_phi,kz_phi,phase_phi)
          !f(:,:,:,ispecial+0)=tanh(f(:,:,:,ispecial+0)/width)
          call hat(amplphi,f,ispecial+0,width,kx_phi,ky_phi,kz_phi)
          f(:,:,:,ispecial+0)=f(:,:,:,ispecial+0)+offset
          f(:,:,:,ispecial+1)=0.
          f(:,:,:,ispecial+2)=0.
          f(:,:,:,ispecial+3)=0.
        case default
          call fatal_error("init_special: No such value for initspecial:" &
              ,trim(initspecial))
      endselect
!
      call keep_compiler_quiet(f)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
!
      if (lmagnetic .and. lbackreact_infl) then
        lpenc_requested(i_infl_a21)=.true.
      endif
!
    endsubroutine pencil_criteria_special
!***********************************************************************
    subroutine pencil_interdep_special(lpencil_in)
!
!  Interdependency among pencils provided by this module are specified here.
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
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
! infl_phi
      if (lpencil(i_infl_phi)) p%infl_dphi=f(l1:l2,m,n,ispecial+0)
!
! infl_dphi
      if (lpencil(i_infl_dphi)) p%infl_dphi=f(l1:l2,m,n,ispecial+1)
!
! infl_a2
      if (lpencil(i_infl_a2)) p%infl_a2=exp(2.*f(l1:l2,m,n,ispecial+3))
!
! infl_a21
      if (lpencil(i_infl_a21)) p%infl_a21=exp(-2.*f(l1:l2,m,n,ispecial+3))
!
!  Magnetic field needed for Maxwell stress
!
      if (lmagnetic) then
        lpenc_requested(i_bb)=.true.
        lpenc_requested(i_el)=.true.
      endif
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
!       if (headtt) call identify_bcs('special',ispecial)
!
      phi=f(l1:l2,m,n,ispecial+0)
      dphi=f(l1:l2,m,n,ispecial+1)
      Hscript=f(l1:l2,m,n,ispecial+2)
      lnascale=f(l1:l2,m,n,ispecial+3)
      ascale=exp(lnascale)
      a2scale=ascale**2
      a2rhop=dphi**2
!
!  Possibility of turning off evolution of scale factor and Hubble parameter
!
      if (lzeroHubble) then
        a2scale=1.
        Hscript=0.
      endif
;
!  Choice of different potentials
!
      select case (Vprime_choice)
        case ('quadratic'); Vprime=axionmass2*phi
        case ('quartic'); Vprime=axionmass2*phi+(lambda_axion/6.)*phi**3
        case ('cos-profile'); Vprime=axionmass2*lambda_axion*sin(lambda_axion*phi)
        case default
          call fatal_error("init_special: No such value for initspecial:" &
              ,trim(initspecial))
      endselect
!
!  Update df.
!  dphi/dt = psi
!  dpsi/dt = - ...
!
        df(l1:l2,m,n,ispecial+0)=df(l1:l2,m,n,ispecial+0)+f(l1:l2,m,n,ispecial+1)
        df(l1:l2,m,n,ispecial+1)=df(l1:l2,m,n,ispecial+1)-2.*Hscript*dphi-a2scale*Vprime
        df(l1:l2,m,n,ispecial+2)=df(l1:l2,m,n,ispecial+2)-4.*pi*a2rhop+Hscript**2
        df(l1:l2,m,n,ispecial+3)=df(l1:l2,m,n,ispecial+3)+Hscript
!
!  speed of light term
!
        if (c_light_axion/=0.) then
          call del2(f,ispecial+0,del2phi)
          df(l1:l2,m,n,ispecial+1)=df(l1:l2,m,n,ispecial+1)+c_light_axion**2*del2phi
        endif
!
!  magnetic terms
!
      if (lmagnetic .and. lbackreact_infl) then
        call dot_mn(p%el,p%bb,tmp)
        df(l1:l2,m,n,ispecial+1)=df(l1:l2,m,n,ispecial+1)+alpf*p%infl_a21*tmp
      endif
!
!  Diagnostics
!
      if (ldiagnos) then
        if (idiag_phim/=0) call sum_mn_name(phi,idiag_phim)
        if (idiag_dphim/=0) call sum_mn_name(dphi,idiag_dphim)
        if (idiag_Hubblem/=0) call sum_mn_name(Hscript,idiag_Hubblem)
        if (idiag_lnam/=0) call sum_mn_name(lnascale,idiag_lnam)
      endif

      call keep_compiler_quiet(p)
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
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_phim=0; idiag_dphim=0; idiag_Hubblem=0; idiag_lnam=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'phim',idiag_phim)
        call parse_name(iname,cname(iname),cform(iname),'dphim',idiag_dphim)
        call parse_name(iname,cname(iname),cform(iname),'Hubblem',idiag_Hubblem)
        call parse_name(iname,cname(iname),cform(iname),'lnam',idiag_lnam)
      enddo
!!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        call farray_index_append('idiag_SPECIAL_DIAGNOSTIC',idiag_SPECIAL_DIAGNOSTIC)
!!      endif
!!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
!
!  26-jun-06/tony: dummy
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(slices%ready)
!
    endsubroutine get_slices_special
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
        subroutine special_calc_particles(f,df,fp,dfp,ineargrid)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
      real, dimension (mx,my,mz,mvar), intent(in) :: df
      real, dimension (:,:), intent(in) :: fp,dfp
      integer, dimension(:,:) :: ineargrid
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(fp,dfp)
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
    subroutine special_calc_spectra(f,spectrum,spectrumhel,lfirstcall,kind)

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spectrum,spectrumhel
      logical :: lfirstcall
      character(LEN=3) :: kind

      call keep_compiler_quiet(f)
      call keep_compiler_quiet(spectrum,spectrumhel)
      call keep_compiler_quiet(lfirstcall)
      call keep_compiler_quiet(kind)
  
    endsubroutine special_calc_spectra
!***********************************************************************
    subroutine special_calc_spectra_byte(f,spectrum,spectrumhel,lfirstcall,kind,len)

      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (:) :: spectrum,spectrumhel
      logical :: lfirstcall
      integer(KIND=ikind1), dimension(3) :: kind
      integer :: len

      !call keep_compiler_quiet(char(kind))
      call keep_compiler_quiet(len)
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(spectrum,spectrumhel)
      call keep_compiler_quiet(lfirstcall)
      
    endsubroutine special_calc_spectra_byte
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
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine special_after_boundary
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
    subroutine special_particles_after_dtsub(f, dtsub, fp, dfp, ineargrid)
!
!  Possibility to modify fp in the end of a sub-time-step.
!
!  28-aug-18/ccyang: coded
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, intent(in) :: dtsub
      real, dimension(:,:), intent(in) :: fp, dfp
      integer, dimension(:,:), intent(in) :: ineargrid
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(dtsub)
      call keep_compiler_quiet(fp)
      call keep_compiler_quiet(dfp)
      call keep_compiler_quiet(ineargrid)
!
    endsubroutine special_particles_after_dtsub
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
endmodule Special
