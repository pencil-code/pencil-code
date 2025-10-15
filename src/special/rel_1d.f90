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
! Declare (for generation of rel_1d_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 2
!
! PENCILS PROVIDED bet
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
module rel_1d
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
  integer :: ieee=0, isss=0, ibss=0, ibet=0, ippp=0
  real, dimension (mx) :: ralp
  real :: amplspecial=.1, r0=3., width=.1
  real :: nu=0., alp=0.
  !real :: bet0=10*tini
  real :: bet0=1e-3
  logical :: lbet_as_aux=.false.
!
  character (len=labellen), dimension(ninit) :: initspecial='nothing'
!
  namelist /rel_1d_init_pars/ &
      initspecial, amplspecial, r0, width, alp, lbet_as_aux
!
  namelist /rel_1d_run_pars/ &
      nu, alp, lbet_as_aux
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_betm=0      ! DIAG_DOC: $\left<\beta\right>$
  integer :: idiag_betmax=0    ! DIAG_DOC: $\beta_{\max}$
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
      call farray_register_pde('eee',ieee)
      call farray_register_pde('sss',isss)
!
      call farray_register_auxiliary('bss',ibss)
      call farray_register_auxiliary('ppp',ippp)
!
      if (lbet_as_aux) call farray_register_auxiliary('bet',ibet)

    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      ralp=x**alp
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
      real :: Vpotential, eps=.01, Hubble_ini, lnascale
      integer :: j
!
      intent(inout) :: f
!
!  SAMPLE IMPLEMENTATION
!
      f(:,:,:,ieee)=0.
      f(:,:,:,isss)=0.
      do j=1,ninit
        select case (initspecial(j))
          case ('nothing'); if (lroot) print*,'init_special: nothing for j=',j
          case ('shock')
            f(:,m1,n1,ieee)=f(:,m1,n1,ieee)+((1.-amplspecial)*.5*(1.-tanh((x-r0)/width))+amplspecial)*ralp
          case default
            call fatal_error("init_special: No such value for initspecial:" &
                ,trim(initspecial(j)))
        endselect
      enddo
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special
!
!  All pencils that this special module depends on are specified here.
!
!  18-jul-06/tony: coded
!
      lpenc_requested(i_bet)=.true.
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
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (nx) :: rat
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
! bet
      if (lpencil(i_bet)) then
        p%bet=.75*f(l1:l2,m,n,isss)/f(l1:l2,m,n,ieee)
        where (p%bet >= bet0)
          rat=2.*f(l1:l2,m,n,ieee)/f(l1:l2,m,n,isss)
          p%bet=rat-sqrt(rat**2-3.)
        endwhere
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
      use Diagnostics, only: sum_mn_name, max_mn_name
      use Deriv, only: der
      use Sub, only: del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      real, dimension (nx) :: del2eee, del2sss, dbss, dsss, dppp, &
        diffus_special
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
      call del2(f,ieee,del2eee)
      call del2(f,isss,del2sss)
!
      call der(f,isss,dsss,1)
      call der(f,ibss,dbss,1)
      call der(f,ippp,dppp,1)

!  Update df
!
        df(l1:l2,m,n,ieee)=df(l1:l2,m,n,ieee)+nu*del2eee-dsss
        df(l1:l2,m,n,isss)=df(l1:l2,m,n,isss)+nu*del2sss-dbss-ralp(l1:l2)*dppp
!
!  For the timestep calculation, need maximum diffusion
!
      if (lfirst.and.ldt) then
        diffus_special=nu*dxyz_2
        if (headtt.or.ldebug) print*,'special: max(diffus_special) =', &
                                      maxval(diffus_special)
        maxdiffus=max(maxdiffus,diffus_special)
      endif
!
!  Diagnostics
!
      if (ldiagnos) then
        if (idiag_betm/=0) call sum_mn_name(p%bet,idiag_betm)
        if (idiag_betmax/=0) call max_mn_name(p%bet,idiag_betmax)
      endif

    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=rel_1d_init_pars, IOSTAT=iostat)
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=rel_1d_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(iostat)
!
      use File_io, only: parallel_unit
!
      integer, intent(out) :: iostat
!
      read(parallel_unit, NML=rel_1d_run_pars, IOSTAT=iostat)
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit, NML=rel_1d_run_pars)
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
        idiag_betm=0; idiag_betmax=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'betm',idiag_betm)
        call parse_name(iname,cname(iname),cform(iname),'betmax',idiag_betmax)
      enddo
!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='bet'.or.cnamev=='ppp') cformv='DEFINED'
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of Special variables.
!
!  26-jun-06/tony: dummy
!
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (slice_data) :: slices
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  bet
!
        case ('bet')
          if (lbet_as_aux) then
            call assign_slices_scal(slices,f,ibet)
          else
            call fatal_error('get_slices_special','bet not as aux')
          endif
!
!  ppp
!
        case ('ppp')
          call assign_slices_scal(slices,f,ippp)
!
      endselect
!
    endsubroutine get_slices_special
!***********************************************************************
    subroutine special_after_boundary(f)
!
!  Possibility to modify the f array after the boundaries are
!  communicated.
!
!  28-dec-21/axel: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      real, dimension (mx) :: bet, bet2, rat
!
      bet=.75*f(:,m1,n1,isss)/f(:,m1,n1,ieee)
      where (bet >= bet0)
        rat=2.*f(:,m1,n1,ieee)/f(:,m1,n1,isss)
        bet=rat-sqrt(rat**2-3.)
      endwhere
      bet2=bet**2
      f(:,m1,n1,ibss)=bet*f(:,m1,n1,isss)
      f(:,m1,n1,ippp)=f(:,m1,n1,ieee)*(1.-bet2)/(ralp*(3.+bet2))
      if (lbet_as_aux) f(:,m1,n1,ibet)=bet
!
    endsubroutine special_after_boundary
!***********************************************************************
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../rel_1d_dummies.inc'
!********************************************************************
endmodule rel_1d
