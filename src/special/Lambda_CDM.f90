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
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!! PENCILS PROVIDED infl_phi; infl_dphi; gphi(3)
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
! Declare index of new variables in f array (if any).
! ascale_type and nconformal are defined in cdata.
!
  integer :: iLCDM_lna=0, iLCDM_tph=0
  real :: Omega_Lam=.73, Omega_rad=1e-4, Omega_mat, Hubble0=0.072
  real :: lna, tph, redshift0=4500.
!
  namelist /special_init_pars/ &
      Omega_Lam, Omega_rad, Hubble0, redshift0, nconformal, ascale_type
!
  namelist /special_run_pars/ &
      Omega_Lam, Omega_rad, Hubble0, ascale_type
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_redshift=0 ! DIAG_DOC: redshift $z$
  integer :: idiag_Hubble=0   ! DIAG_DOC: $H(a)$
  integer :: idiag_ascale=0   ! DIAG_DOC: $a$
  integer :: idiag_lna=0      ! DIAG_DOC: $\ln a$
  integer :: idiag_tph=0      ! DIAG_DOC: $t_\mathrm{phys}$
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
      use SharedVariables, only: put_shared_variable
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_ode('LLCDM_lna',iLCDM_lna)
      call farray_register_ode('iLCDM_tph',iLCDM_tph)
!
    endsubroutine register_special
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
!  Adjust Omega_mat to conformally flat universe.
!
      Omega_mat=1.-(Omega_Lam+Omega_rad)
      ascale=exp(f_ode(iLCDM_lna))
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
      real :: Vpotential, tph_init
      integer :: j
!
      intent(inout) :: f
!
!  energy density of the charged particles
!
      ascale=1./(1.+redshift0)
      Hubble=Hubble0*sqrt(Omega_mat/ascale**3+Omega_Lam+Omega_rad/ascale**4)
!
!  initial condition for physical time.
!
      tph_init=ascale**2/(2.*Hubble0*Omega_rad**.5)
!
      f_ode(iLCDM_lna)=alog(ascale)
      f_ode(iLCDM_tph)=tph_init
!
      call mpibcast_real(ascale)
      call mpibcast_real(Hubble)
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  24-02-25/axel: adapted
!
!     if (ltemperature .and..not. lentropy .and..not. lthermal_energy) lpenc_requested(i_TT)=.true.
!
!  Magnetic field needed for Maxwell stress
!
!     if (lmagnetic) then
!       lpenc_requested(i_bb)=.true.
!       lpenc_requested(i_el)=.true.
!       if (lrho_chi .or. lnoncollinear_EB .or. lnoncollinear_EB_aver) lpenc_requested(i_e2)=.true.
!     endif
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
!     if (lpencil(i_infl_phi)) p%infl_phi=f(l1:l2,m,n,iinfl_phi)
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
!   6-oct-03/tony: coded
!   2-nov-21/axel: first set of equations coded
!
!     use Diagnostics, only: sum_mn_name, max_mn_name, save_name
!     use Sub, only: dot_mn, del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
!     real, dimension (nx) :: phi, dphi, Vprime
!     real, dimension (nx) :: tmp, del2phi
!     real :: tmp2
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine dspecial_dt_ode
!
      use Diagnostics, only: save_name
      use SharedVariables, only: get_shared_variable
!
      lna=f_ode(iLCDM_lna)
      tph=f_ode(iLCDM_tph)
!
!  dlna/dtph=H, and since dt=dtph/a^nconformal, we have
!  so dlna/dt=dlna/dtph*dtph/dt=H*dtph/dt=H*a^nconformal.
!
      df_ode(iLCDM_lna)=df_ode(iLCDM_lna)+ascale**nconformal*Hubble
      df_ode(iLCDM_tph)=df_ode(iLCDM_tph)+ascale**nconformal
!
!  Diagnostics
!
      if (ldiagnos) then
        call save_name(1./ascale-1.,idiag_redshift)
        call save_name(Hubble,idiag_Hubble)
        call save_name(ascale,idiag_ascale)
        call save_name(lna,idiag_lna)
        call save_name(tph,idiag_tph)
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
        idiag_redshift=0; idiag_Hubble=0; idiag_ascale=0
        idiag_lna=0; idiag_tph=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'redshift',idiag_redshift)
        call parse_name(iname,cname(iname),cform(iname),'Hubble',idiag_Hubble)
        call parse_name(iname,cname(iname),cform(iname),'ascale',idiag_ascale)
        call parse_name(iname,cname(iname),cform(iname),'lna',idiag_lna)
        call parse_name(iname,cname(iname),cform(iname),'tph',idiag_tph)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  06-jul-06/tony: coded
!
!     use Mpicomm, only: mpireduce_sum, mpiallreduce_sum, mpibcast_real
!     use Sub, only: dot2_mn, grad, curl, dot_mn
!
      real, dimension (mx,my,mz,mfarray), intent(in) :: f
!
!  Compute terms routinely used during this time substep.
!
      lna=f_ode(iLCDM_lna)
      ascale=exp(lna)
      sqrt_ascale=sqrt(ascale)
      Hubble=Hubble0*sqrt(Omega_mat/ascale**3+Omega_Lam+Omega_rad/ascale**4)
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
