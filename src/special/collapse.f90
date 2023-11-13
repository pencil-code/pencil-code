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
! PENCILS PROVIDED rrr, drrr, d2rrr, bet, mass, dlnrho, grav
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
  real :: mass_blackhole=0., cs2=onethird
!
  character (len=labellen), dimension(ninit) :: initspecial='nothing'
!
  namelist /special_init_pars/ &
      initspecial, mass_blackhole, cs2
!
  namelist /special_run_pars/ &
      initspecial, mass_blackhole, cs2
!
! Diagnostic variables (needs to be consistent with reset list below).
!
  integer :: idiag_betm=0      ! DIAG_DOC: $\left<\beta\right>$
  integer :: idiag_massm=0      ! DIAG_DOC: $\left<m\right>$
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
      real, dimension(:), allocatable :: m13_file, r_file
      real :: dm13_file, m13, m13a, m13b, r, r1, r2
      integer :: j, im13_file, nm13_file, l
!
      intent(inout) :: f
!
!  Initial condition.
!
      do j=1,ninit
        select case (initspecial(j))
          case ('nothing'); if (lroot) print*,'init_special: nothing'
          case ('r_vs_m13')
            open(1,file='r_vs_m13.dat',status='old')
            read(1,*) nm13_file, dm13_file
          if (allocated(m13_file)) deallocate(m13_file, r_file)
          allocate(m13_file(nm13_file), r_file(nm13_file))
          do im13_file=1,nm13_file
            read(1,*) m13_file(im13_file), r_file(im13_file)
          enddo
          close(1)
          do l=l1,l2
            m13=x(l)**onethird
            im13_file=1+int(m13/dm13_file)
!
!  Check upper limit:
!
            if (im13_file+1 >= nm13_file) then
              print*,x(l)**onethird,im13_file
              call fatal_error("init_special","im13_file too large")
            endif
!
!  Linear interpolation:
!
            r1=r_file(im13_file)
            r2=r_file(im13_file+1)
            m13a=m13_file(im13_file)
            m13b=m13_file(im13_file+1)
            r=r1+(r2-r1)*(m13-m13a)/(m13b-m13a)
            if (mod(l,100)==0) print*,'interpol(sampled): ',m13a, m13, m13b, r1, r, r2
            f(l,m1,n1,irrr)=f(l,m1,n1,irrr)+r
          enddo
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
!  Request pencils.
!
      lpenc_requested(i_mass)=.true.
      lpenc_requested(i_bet)=.true.
      lpenc_requested(i_rrr)=.true.
      lpenc_requested(i_drrr)=.true.
      lpenc_requested(i_d2rrr)=.true.
      lpenc_requested(i_dlnrho)=.true.
      lpenc_requested(i_grav)=.true.
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
      use Deriv, only: der, der2
!
      real, dimension (mx,my,mz,mfarray) :: f
      type (pencil_case) :: p
!
      intent(in) :: f
      intent(inout) :: p
!
! mass
      if (lpencil(i_mass)) p%mass=x(l1:l2)+mass_blackhole
! bet
      if (lpencil(i_bet)) p%bet=f(l1:l2,m,n,ibet)
! rrr
      if (lpencil(i_rrr)) p%rrr=f(l1:l2,m,n,irrr)
! drrr
      if (lpencil(i_drrr)) call der(f,irrr,p%drrr,1)
! d2rrr
      if (lpencil(i_d2rrr)) call der2(f,irrr,p%d2rrr,1)
! dlnrho
      if (lpencil(i_dlnrho)) p%dlnrho=-(p%d2rrr/p%drrr**2+2./p%rrr)
! grav
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
      df(l1:l2,m,n,ibet)=df(l1:l2,m,n,ibet)+p%grav-cs2*p%dlnrho
!
!  Diagnostics
!
      if (ldiagnos) then
        call sum_mn_name(p%bet,idiag_betm)
        call sum_mn_name(p%mass,idiag_massm)
      endif

    endsubroutine dspecial_dt
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
      logical :: lreset,lwr,lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        idiag_betm=0; idiag_massm=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'betm',idiag_betm)
        call parse_name(iname,cname(iname),cform(iname),'massm',idiag_massm)
      enddo
!!
!  check for those quantities for which we want video slices
!
      if (lwrite_slices) then
        where(cnamev=='rrr') cformv='DEFINED'
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
      endif
!
    endsubroutine rprint_special
!***********************************************************************
    subroutine get_slices_special(f,slices)
!
!  Write slices for animation of electric potential
!
!  26-feb-07/axel: adapted from gross_pitaevskii
!
      use Slices_methods, only: assign_slices_scal
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      type (slice_data) :: slices
!
      integer :: inamev
!
!  Loop over slices
!
      select case (trim(slices%name))
!
!  Electric field.
!
      case ('rrr'); call assign_slices_scal(slices,f,irrr)
!
      endselect
!
    endsubroutine get_slices_special
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
                                         