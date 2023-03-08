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
! PENCILS PROVIDED rrr, drrr, d2rrr, bet, mass, drho, grav
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
  integer :: irrr=0, ibet=0
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
      initspecial
!
  namelist /special_run_pars/ &
      initspecial
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_betm=0      ! DIAG_DOC: $\left<\beta\right>$
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
      call farray_register_pde('rrr',irrr)
      call farray_register_pde('bet',ibet)
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
!  pencil for ...
!
      lpenc_requested(i_rrr)=.true.
      lpenc_requested(i_drrr)=.true.
      lpenc_requested(i_d2rrr)=.true.
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
      use Deriv, only: der
!
      real, dimension (mx,my,mz,mfarray) :: f
      real :: cs2=onethird
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
! mass
      if (lpencil(i_mass)) p%mass=x(l1:l2)
! bet
      if (lpencil(i_bet)) p%bet=f(l1:l2,m,n,ibet)
! rrr
      if (lpencil(i_rrr)) p%rrr=f(l1:l2,m,n,irrr)
! drrr
      if (lpencil(i_drrr)) call der(f,irrr,p%drrr,1)
! d2rrr
      if (lpencil(i_d2rrr)) call der(f,irrr,p%drrr,2)
! drho
      if (lpencil(i_drho)) p%drho=-cs2*(p%d2rrr/(p%rrr*p%drrr)**2+2./p%rrr**3)/(4.*pi*p%drrr)
! drho
      if (lpencil(i_grav)) p%grav=-p%mass/abs(p%rrr*(p%rrr-2.*p%mass))

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
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
      df(l1:l2,m,n,irrr)=df(l1:l2,m,n,irrr)+p%bet/sqrt(1.+p%bet**2)
      df(l1:l2,m,n,ibet)=df(l1:l2,m,n,ibet)+p%grav-p%drho
!
!  Diagnostics
!
      if (ldiagnos) then
        call sum_mn_name(p%bet,idiag_betm)
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
        idiag_betm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'betm',idiag_betm)
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
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
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
