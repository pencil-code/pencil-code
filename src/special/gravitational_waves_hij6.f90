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
! MVAR CONTRIBUTION 12
!
! PENCILS PROVIDED stress_ij(6)
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
  use Initcond
  use General, only: keep_compiler_quiet
  use Messages, only: svn_id, fatal_error
!
  implicit none
!
  include '../special.h'
!
! Declare index of new variables in f array (if any).
!
  character (len=labellen) :: inithij='nothing'
  character (len=labellen) :: initgij='nothing'
  real :: amplhij=0., amplgij=0.
  real :: kx_hij=0., ky_hij=0., kz_hij=0.
  real :: kx_gij=0., ky_gij=0., kz_gij=0.
  real :: dummy=0., diffhh=0., diffgg=0.
  logical :: lno_transverse_part=.false., lsame_diffgg_as_hh=.true.
  real, dimension(3,3) :: ij_table
!
! input parameters
  namelist /special_init_pars/ &
    inithij, initgij, amplhij, amplgij, &
    kx_hij, ky_hij, kz_hij, &
    kx_gij, ky_gij, kz_gij
!
! run parameters
  namelist /special_run_pars/ &
    diffhh, diffgg, lsame_diffgg_as_hh
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_g22pt=0       ! DIAG_DOC: $g_{22}(x_1,y_1,z_1,t)$
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      use Sub, only: register_report_aux
      use FArrayManager
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_pde('hij',ihij,vector=6)
      call farray_register_pde('gij',igij,vector=6)
!
      if (lroot) call svn_id( &
           "$Id$")
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  Called after reading parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical, pointer :: lbb_as_comaux
!
!  Check whether diffgg=diffhh (which  is the default)
!
      if (lmagnetic) then
        call get_shared_variable('lbb_as_comaux',lbb_as_comaux)
        if (.not.lbb_as_comaux) call fatal_error('initialize_special', &
            'lbb_as_comaux needs to be T')
      endif
!
!  Check whether diffgg=diffhh (which  is the default)
!
      if (lsame_diffgg_as_hh) diffgg=diffhh
!
!  set index table
!
      ij_table(1,1)=1
      ij_table(2,2)=2
      ij_table(3,3)=3
      ij_table(1,2)=4
      ij_table(2,3)=5
      ij_table(3,1)=6
      ij_table(2,1)=4
      ij_table(3,2)=5
      ij_table(1,3)=6
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
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  initial condition for hij
!
      select case (inithij)
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('coswave-kx'); call coswave(amplhij,f,ihij,kx=kx_hij)
        case default
          call fatal_error("init_special: No such value for inithij:" &
              ,trim(inithij))
      endselect
!
!  initial condition for gij
!
      select case (initgij)
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('coswave-kx'); call coswave(amplgij,f,igij,kx=kx_gij)
        case default
          call fatal_error("init_special: No such value for initgij:" &
              ,trim(initgij))
      endselect
!
    endsubroutine init_special
!***********************************************************************
    subroutine pencil_criteria_special()
!
!  All pencils that this special module depends on are specified here.
!
!  18-07-06/tony: coded
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
      integer :: i, j, ij
!
!  Construct stress tensor
!
      p%stress_ij=0.0
      do j=1,3
      do i=1,j
        ij=ij_table(i,j)
        if (lhydro) p%stress_ij(:,ij)=p%stress_ij(:,ij) &
          +f(l1:l2,m,n,iux+i-1) &
          *f(l1:l2,m,n,iux+j-1)
        if (lmagnetic) p%stress_ij(:,ij)=p%stress_ij(:,ij) &
          +f(l1:l2,m,n,ibx+i-1) &
          *f(l1:l2,m,n,ibx+j-1)
      enddo
      enddo
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
!
      use Diagnostics
      use Sub, only: del2
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,6) :: del2hij, del2gij
      type (pencil_case) :: p
!
      integer :: ij,jhij,jgij
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  Identify module and boundary conditions.
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dspecial_dt'
!
      do ij=1,6
        jhij=ihij-1+ij
        jgij=igij-1+ij
        call del2(f,jhij,del2hij(:,ij))
        df(l1:l2,m,n,jhij)=df(l1:l2,m,n,jhij)+f(l1:l2,m,n,jgij)
        df(l1:l2,m,n,jgij)=df(l1:l2,m,n,jgij)+del2hij(:,ij)+p%stress_ij(:,ij)
        if (diffhh/=0.) then
          call del2(f,jhij,del2hij(:,ij))
          df(l1:l2,m,n,jhij)=df(l1:l2,m,n,jhij)+diffhh*del2hij(:,ij)
        endif
        if (diffgg/=0.) then
          call del2(f,jgij,del2gij(:,ij))
          df(l1:l2,m,n,jgij)=df(l1:l2,m,n,jgij)+diffgg*del2gij(:,ij)
        endif
      enddo
!
!  diagnostics
!
       if (ldiagnos) then
!        if (lroot.and.m==mpoint.and.n==npoint) then
!          if (idiag_g22pt/=0) call save_name(p%bb(lpoint-nghost,2),idiag_g22pt)
!        endif
       endif
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
    subroutine special_before_boundary(f)
!
!  Possibility to modify the f array before the boundaries are
!  communicated.
!
!  Some precalculated pencils of data are passed in for efficiency
!  others may be calculated directly from the f array
!
!  30-mar-17/axel: moved stuff from special_after_boundary to here
!
      !use Boundcond, only: zero_ghosts, update_ghosts
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
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
    endsubroutine special_after_boundary
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics
!
      integer :: iname
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!!!
!!!  reset everything in case of reset
!!!  (this needs to be consistent with what is defined above!)
!!!
      if (lreset) then
        idiag_g22pt=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'g22pt',idiag_g22pt)
      enddo
!!
!!!  write column where which magnetic variable is stored
!!      if (lwr) then
!!        write(3,*) 'i_SPECIAL_DIAGNOSTIC=',i_SPECIAL_DIAGNOSTIC
!!      endif
!!
    endsubroutine rprint_special
!***********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
