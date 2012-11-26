! $Id: twist_inject.f90 19193 2012-06-30 12:55:46Z wdobler $
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
!  Special term in the entropy equation            | special_calc_entropy
!  Special term in the induction (magnetic)        | special_calc_magnetic
!     equation                                     |
!                                                  |
!  Special equation                                | dspecial_dt
!    NOT IMPLEMENTED FULLY YET - HOOKS NOT PLACED INTO THE PENCIL-CODE
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
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
!!  namelist /special_init_pars/ dummy
!
  real :: fring=1e-3,r0=0.2,tilt=0.0,width=0.02,&
          dIring=0.0,dposx=0.0,dtilt=0.0,Ilimit=0.15,poslimit=0.98
  real :: posz=0.0
  real, save :: posx,Iring
  logical :: lset_boundary_emf=.false.
  namelist /special_run_pars/ Iring,dIring,fring,r0,width,&
           posx,dposx,posz,tilt,dtilt,Ilimit,poslimit,lset_boundary_emf
!
! Declare index of new variables in f array (if any).
!
!!   integer :: ispecial=0
!!   integer :: ispecaux=0
!
!! Diagnostic variables (needs to be consistent with reset list below).
!
   integer :: idiag_posx=0,idiag_Iring=0
!
  contains
!***********************************************************************
    subroutine register_special()
!
!  Set up indices for variables in special modules.
!
!  6-oct-03/tony: coded
!
      if (lroot) call svn_id( &
           "$Id: nospecial.f90 19193 2012-06-30 12:55:46Z wdobler $")
!
!!      call farray_register_pde('special',ispecial)
!!      call farray_register_auxiliary('specaux',ispecaux)
!!      call farray_register_auxiliary('specaux',ispecaux,communicated=.true.)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  Called by start.f90 together with lstarting=.true.   and then
!  called by run.f90   together with lstarting=.false.  after reading
!  parameters, but before the time loop.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
      if (lfargo_advection) then
        print*,''
        print*,'Switch '
        print*,' SPECIAL = special/fargo'
        print*,'in src/Makefile.local if you want to use the fargo algorithm'
        print*,''
        call fatal_error('nospecial','initialize_special()')
      endif
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine finalize_special(f,lstarting)
!
!  Called by start.f90 together with lstarting=.true.   and then
!  called by run.f90   together with lstarting=.false.  before exiting.
!
!  14-aug-2011/Bourdin.KIS: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      logical, intent(in) :: lstarting
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
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
!!
!!  SAMPLE IMPLEMENTATION
!!
!!      select case (initspecial)
!!        case ('nothing'); if (lroot) print*,'init_special: nothing'
!!        case ('zero', '0'); f(:,:,:,iSPECIAL_VARIABLE_INDEX) = 0.
!!        case default
!!          call fatal_error("init_special: No such value for initspecial:" &
!!              ,trim(initspecial))
!!      endselect
!
      call keep_compiler_quiet(f)
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
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(p)
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
!!      if (headtt) call identify_bcs('special',ispecial)
!
!!
!! SAMPLE DIAGNOSTIC IMPLEMENTATION
!!
!!      if (ldiagnos) then
!!        if (idiag_SPECIAL_DIAGNOSTIC/=0) then
!!          call sum_mn_name(MATHEMATICAL EXPRESSION,idiag_SPECIAL_DIAGNOSTIC)
!!! see also integrate_mn_name
!!        endif
!!      endif
!
      call keep_compiler_quiet(f,df)
      call keep_compiler_quiet(p)
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      call keep_compiler_quiet(unit)
      call keep_compiler_quiet(present(iostat))
!
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
      call keep_compiler_quiet(unit)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99    return
!
    endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
!
      integer, intent(in) :: unit
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  Reads and registers print parameters relevant to special.
!
!  06-oct-03/tony: coded
!
      use Diagnostics
      integer :: iname
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
        idiag_posx=0
        idiag_Iring=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),&
            'posx',idiag_posx)
        call parse_name(iname,cname(iname),cform(iname),&
            'Iring',idiag_Iring)
      enddo
!  write column where which magnetic variable is stored
      if (lwr) then
        write(3,*) 'idiag_posx=',idiag_posx
        write(3,*) 'idiag_Iring=',idiag_Iring
      endif
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
    subroutine calc_lspecial_pars(f)
!
!  Dummy routine.
!
!  15-jan-08/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
      call keep_compiler_quiet(f)
!
    endsubroutine calc_lspecial_pars
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
    subroutine special_calc_entropy(f,df,p)
!
!  Calculate an additional 'special' term on the right hand side of the
!  entropy equation.
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
    endsubroutine special_calc_entropy
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
    subroutine special_calc_particles(fp)
!
!  Called before the loop, in case some particle value is needed
!  for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (:,:), intent(in) :: fp
!
      call keep_compiler_quiet(fp)
!
    endsubroutine special_calc_particles
!***********************************************************************
    subroutine special_calc_particles_nbody(fsp)
!
!  Called before the loop, in case some massive particles value
!  is needed for the special density/hydro/magnetic/entropy.
!
!  20-nov-08/wlad: coded
!
      real, dimension (:,:), intent(in) :: fsp
!
      call keep_compiler_quiet(fsp)
!
    endsubroutine special_calc_particles_nbody
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
    subroutine special_boundconds(f,bc)
!
!  Some precalculated pencils of data are passed in for efficiency,
!  others may be calculated directly from the f array.
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray), intent(inout) :: f
      character (len=3) :: topbot
      type (boundary_condition), intent(inout) :: bc
      integer :: j
!
      select case (bc%bcname)
      case ('nfc')
      j=bc%ivar
        select case (bc%location)
        case (iBC_X_TOP)
          topbot='top'
          call bc_nfc_x(f,topbot,j)
      bc%done=.true.
        case (iBC_X_BOT)
          topbot='bot'
          call bc_nfc_x(f,topbot,j)
      bc%done=.true.
        end select
      case ('go')
      j=bc%ivar
        select case (bc%location)
        case (iBC_X_TOP)
          topbot='top'
          call bc_go_x(f,topbot,j)
      bc%done=.true.
        case (iBC_X_BOT)
          topbot='bot'
          call bc_go_x(f,topbot,j)
      bc%done=.true.
        end select
      end select

!
    endsubroutine special_boundconds
!***********************************************************************
    subroutine special_after_timestep(f,df,dt_)
!
!  Possibility to modify the f and df after df is updated.
!  Used for the Fargo shift, for instance.
!
!  27-nov-08/wlad: coded
!
      use Deriv, only: der
      use Mpicomm, only: mpibcast_double
      use Diagnostics, only: save_name
!      use Sub, only: cross,curl_mn,gij,gij_etc
      real, dimension(mx,my,mz,mfarray), intent(inout) :: f
      real, dimension(mx,my,mz,mvar), intent(inout) :: df
      real, intent(in) :: dt_
!      real, dimension(nx,3,3) :: aij,bij
      real, dimension(nx) :: dfy,dfz
      real, dimension(3) :: tmpv,vv,uu
      real :: xi,xx0,yy0,zz0,xx1,yy1,zz1,dist,distxy,distyz,phi,rr,r1,&
              prof,ymid,zmid,umax
      real :: tmpx,tmpy,tmpz,posxold,Iringold
      logical :: lcorner,lcorner_y,lcorner_z
      integer :: l,k
!
!  IMPLEMENTATION OF INSERTION OF BIPOLES (Non-Potential part)
!  (Yeates, Mackay and van Ballegooijen 2008, Sol. Phys, 247, 103)
!
      ymid=0.5*(xyz0(2)+xyz1(2))
      zmid=0.5*(xyz0(3)+xyz1(3))
      if (lroot) then
        posxold=posx
        Iringold=Iring
        if (posxold.gt.poslimit) dposx=0.0
        if (Iringold.gt.Ilimit) dIring=0.0
        posx=posx+dposx*dt_
        Iring=Iring+dIring*dt_
        tilt=tilt+dtilt*dt_
      endif
      call mpibcast_double(posxold)
      call mpibcast_double(posx)
      call mpibcast_double(dposx)
      call mpibcast_double(Iringold)
      call mpibcast_double(Iring)
      call mpibcast_double(dIring)
      if (ldiagnos) then
        if (idiag_posx/=0) &
          call save_name(posx,idiag_posx)
        if (idiag_Iring/=0) &
          call save_name(Iring,idiag_Iring)
      endif
!
      do n=n1,n2
        do m=m1,m2
!
        if (lfirst_proc_x) then
          if (lcartesian_coords) then
            xi=((x(l)**2+z(n)**2)/2.+y(m)**2)/r0**2
            f(l1,m,n,iax) = fring*Iring*z(n)*exp(-2*xi)
            f(l1,m,n,iay) = fring*r0*exp(-xi)
            f(l1,m,n,iaz) = -fring*Iring*x(l)*exp(-2*xi)
          else if (lcylindrical_coords) then
            call fatal_error('special_after_timestep',&
            'Bipoles not coded for cylindrical coordinates')
          endif
!
!  Then set up the helical field
!
          if (lspherical_coords) then
              xx0=x(l1)*sinth(m)*cos(z(n))
              yy0=x(l1)*sinth(m)*sin(z(n))
              zz0=x(l1)*costh(m)
          ! Calculate D^(-1)*(xxx-disp)
              xx1=xx0-posxold
              yy1=cos(tilt*pi/180.0)*yy0+sin(tilt*pi/180.0)*(zz0-posz)
              zz1=-sin(tilt*pi/180.0)*yy0+cos(tilt*pi/180.0)*(zz0-posz)
              phi=atan2(yy1,xx1)
            if (dposx.ne.0.or.dtilt.ne.0) then
              call norm_ring(xx1,yy1,zz1,fring,Iringold,r0,width,posxold,tmpv,PROFILE='gaussian')
            ! calculate D*tmpv
              tmpx=tmpv(1)
              tmpy=cos(tilt*pi/180.0)*tmpv(2)-sin(tilt*pi/180.0)*tmpv(3)
              tmpz=sin(tilt*pi/180.0)*tmpv(2)+cos(tilt*pi/180.0)*tmpv(3)
              dist=sqrt(xx0**2+yy0**2+zz0**2)
              distxy=sqrt(xx0**2+yy0**2)
! Subtract original ring
              f(l1,m,n,iax)=f(l1,m,n,iax)-(xx0*tmpx/dist+yy0*tmpy/dist+zz0*tmpz/dist)
              f(l1,m,n,iay)= f(l1,m,n,iay)-(xx0*zz0*tmpx/(dist*distxy)+yy0*zz0*tmpy/(dist*distxy) &
                      -distxy*tmpz/dist)
              f(l1,m,n,iaz) = f(l1,m,n,iaz)-(-yy0*tmpx/distxy+xx0*tmpy/distxy)
! Add new ring
              xx1=xx0-posx
              phi=atan2(yy1,xx1)
              call norm_ring(xx1,yy1,zz1,fring,Iring,r0,width,posx,tmpv,PROFILE='gaussian')
            ! calculate D*tmpv
              tmpx=tmpv(1)
              tmpy=cos(tilt*pi/180.0)*tmpv(2)-sin(tilt*pi/180.0)*tmpv(3)
              tmpz=sin(tilt*pi/180.0)*tmpv(2)+cos(tilt*pi/180.0)*tmpv(3)
              f(l1,m,n,iax)=f(l1,m,n,iax)+(xx0*tmpx/dist+yy0*tmpy/dist+zz0*tmpz/dist)
              f(l1,m,n,iay)= f(l1,m,n,iay)+(xx0*zz0*tmpx/(dist*distxy)+yy0*zz0*tmpy/(dist*distxy) &
                      -distxy*tmpz/dist)
              f(l1,m,n,iaz) = f(l1,m,n,iaz)+(-yy0*tmpx/distxy+xx0*tmpy/distxy)
            endif
!
            if (dposx.ne.0) then
          ! Calculate D^(-1)*(xxx-disp)
              distyz=sqrt((sqrt(xx1**2+yy1**2)-r0)**2+zz1**2)
              if (distyz.lt.2*width) then
                vv(1)=dposx
                vv(2)=0.0
                vv(3)=0.0
              else
                vv=0.0
              endif
              tmpx=vv(1)
              tmpy=cos(tilt*pi/180.0)*vv(2)-sin(tilt*pi/180.0)*vv(3)
              tmpz=sin(tilt*pi/180.0)*vv(2)+cos(tilt*pi/180.0)*vv(3)
              f(l1,m,n,iux) = (xx0*tmpx/dist+yy0*tmpy/dist+zz0*tmpz/dist)
              f(l1,m,n,iuy) = (xx0*zz0*tmpx/(dist*distxy)+yy0*zz0*tmpy/(dist*distxy) &
                      -distxy*tmpz/dist)
              f(l1,m,n,iuz) = (-yy0*tmpx/distxy+xx0*tmpy/distxy)
            else if (dIring.eq.0.0.and.dposx.eq.0) then
              f(l1,m,n,iux:iuz)=0.0
            endif
!
          endif
        endif
          if (lset_boundary_emf) then
            lcorner=.false.
            lcorner_y=.false.
            lcorner_z=.false.
            if (llast_proc_y.and.m.eq.m2)  lcorner_y=.true.
            if (lfirst_proc_y.and.m.eq.m1) lcorner_y=.true.
            if (llast_proc_z.and.n.eq.n2)  lcorner_z=.true.
            if (lfirst_proc_z.and.n.eq.n1) lcorner_z=.true.
            if (lcorner_y.or.lcorner_z) lcorner=.true.
            if (lcorner) then
              do k=iay, iaz; call bc_nfc_x(f,'top',k); enddo
            endif
            call bc_emf_x(f,df,dt_,'top',iax)
            call bc_emf_x(f,df,dt_,'top',iay)
            call bc_emf_x(f,df,dt_,'top',iaz)
          endif
        enddo
      enddo
!      if (dIring.ne.0) then
!        call update_ghosts(f)
!        do n=n1,n2
!          do m=m1,m2
!            call der(f,iux,dfz,3)
!            call der(f,iux,dfy,2)
!!            
!! Normalize to unity.          
!!
!            f(l1,m,n,iuy)=dfz(1)
!            f(l1,m,n,iuz)=-dfy(1)
!          enddo
!        enddo
!        f(l1,m1:m2,n1:n2,iux)=0.0
!        call find_umax(f,umax)
!        f(l1,m1:m2,n1:n2,iuy)=dIring*f(l1,m1:m2,n1:n2,iuy)/umax
!        f(l1,m1:m2,n1:n2,iuz)=dIring*f(l1,m1:m2,n1:n2,iuz)/umax
!      endif

!
    endsubroutine  special_after_timestep
!***********************************************************************
    subroutine norm_ring(xx1,yy1,zz1,fring,Iring,r0,width,posx,vv,profile)
!
!  Generate vector potential for a flux ring of magnetic flux FRING,
!  current Iring (not correctly normalized), radius R0 and thickness
!  WIDTH in normal orientation (lying in the x-y plane, centred at (0,0,0)).
!
!   1-may-02/wolf: coded (see the manual, Section C.3)
!   7-jun-09/axel: added gaussian and constant (or box) profiles
!
      use Mpicomm, only: stop_it
      use Sub, only: erfunc
!
      real, dimension (3) :: vv
      real :: xx1,yy1,zz1,phi,tmp,cs,sn
      real :: fring,Iring,r0,width,posx,width2
      character (len=*) :: profile
!
!  magnetic ring, define r-R
!
      tmp = sqrt(xx1**2+yy1**2)-r0
!
!  choice of different profile functions
!
      select case (profile)
!
!  gaussian profile, exp(-.5*(x/w)^2)/(sqrt(2*pi)*eps),
!  so its derivative is .5*(1.+erf(-x/(sqrt(2)*eps))
!
      case ('gaussian')
        width2=50.0*width
        if (abs(zz1).le.width2) then
          if (sqrt(tmp**2+zz1**2).le.width2) then
            vv(3) = + fring * .5*(erfunc(sqrt(width2**2-zz1**2)/width/sqrt(2.))- &
                      erfunc(tmp/(sqrt(2.)*width))) &
                     *exp(-.5*(zz1/width)**2)/(sqrt(2.*pi)*width)
          else 
            if (tmp.lt.0) then
              vv(3) = + fring *(erfunc(sqrt(width2**2-zz1**2)/width/sqrt(2.))) &
                     *exp(-.5*(zz1/width)**2)/(sqrt(2.*pi)*width)
            else
              vv(3)=0.0
            endif         
          endif
        else
          vv(3)=0.0
        endif  
!
!  tanh profile, so the delta function is approximated by 1/cosh^2.
!  The name tanh is misleading, because the actual B frofile is
!  1./cosh^2, but this is harder to write.
!
      case ('tanh')
        vv(3) = - fring * 0.5*(1+tanh(tmp/width)) &
                          * 0.5/width/cosh(zz1/width)**2
!
!  constant profile, so the delta function is approximated by the function
!  delta(x) = 1/2w, if -w < x < w.
!
      case ('const')
        vv(3) = - fring * 0.5*(1.+max(-1.,min(tmp/width,1.))) &
                          * 0.25/width*(1.-sign(1.,abs(zz1)-width))
      case default
        call stop_it('norm_ring: No such fluxtube profile')
      endselect
!
!  current ring (to twist the B-lines)
!
      phi = atan2(yy1,xx1)
      tmp = sqrt((xx1-r0*cos(phi))**2 + (yy1-r0*sin(phi))**2+zz1**2)
      tmp = Iring*fring*exp(-0.5*(tmp/width)**2)/(sqrt(2.*pi)*width)  ! Now the A_phi component
      vv(1) = - tmp*sin(phi)
      vv(2) =   tmp*cos(phi)
!
    endsubroutine norm_ring
!***********************************************************************
    subroutine find_umax(f,umax)
!
!  Find the absolute maximum of the velocity.
!
!  19-aug-2011/ccyang: coded
!  13-sep-2012/piyali: Added this routine from hydro since
!  initial_condition cannot use
!  the hydro module because of Makefile.depend
!
      use Mpicomm, only: mpiallreduce_max
!
      real, dimension(mx,my,mz,mfarray), intent(in) :: f
      real, intent(out) :: umax
!
      real :: umax1
!
!  Find the maximum.
!
      umax1 = sqrt(maxval(f(l1,m1:m2,n1:n2,iux)**2 &
                        + f(l1,m1:m2,n1:n2,iuy)**2 &
                        + f(l1,m1:m2,n1:n2,iuz)**2))
      call mpiallreduce_max(umax1, umax)
!
    endsubroutine find_umax
!***********************************************************************
    subroutine bc_nfc_x(f,topbot,j)
!
!  Normal-field (or angry-hedgehog) boundary condition for spherical
!  coordinate system.
!  d_r(A_{\theta}) = -A_{\theta}/r  with A_r = 0 sets B_{r} to zero
!  in spherical coordinate system.
!  (compare with next subroutine sfree )
!
!  25-Aug-2012/piyali: adapted from bc_set_nfr_x
!
      character (len=bclen), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension(3,3) :: mat
      real, dimension(3) :: rhs
      real :: x2,x3
      integer, intent (in) :: j
      integer :: k
!
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        do k=1,nghost
          if (lfirst_proc_x) &
          f(l1-k,:,:,j)= f(l1+k,:,:,j)*(x(l1+k)/x(l1-k))
        enddo
!
     case ('top')               ! top boundary
       x2=x(l2+1)
       x3=x(l2+2)
       do m=m1, m2
         do n=n1,n2
           mat(1,1:3)=[2*x2*(60*dx+197*x3),x2*(60*dx+167*x3),-(6*x2*x3)]/ &
                       (7500*dx*x2+10020*dx*x3+24475*x2*x3+3600*dx**2)
!
           mat(2,1:3)=[-10*x3*(12*dx+x2), 120*x2*x3, x3*(12*dx+25*x2)]/ &
                       (1500*dx*x2+2004*dx*x3+4895*x2*x3+720*dx**2)
!
           mat(3,1:3)=[420*dx*x2+924*dx*x3+1259*x2*x3+720*dx**2,   &
                       -(9*x2*(60*dx+47*x3)), 9*x3*(12*dx+31*x2)]/ &
                       (1500*dx*x2+2004*dx*x3+4895*x2*x3+720*dx**2)
           if (llast_proc_x) then
!
             rhs(1)=-f(l2,m,n,j)*60*dx/x(l2)+45*f(l2-1,m,n,j)-9*f(l2-2,m,n,j)+f(l2-3,m,n,j)
             rhs(2)=f(l2,m,n,j)*80-30*f(l2-1,m,n,j)+8*f(l2-2,m,n,j)-f(l2-3,m,n,j)
             rhs(3)=-f(l2,m,n,j)*100+50*f(l2-1,m,n,j)-15*f(l2-2,m,n,j)+2*f(l2-3,m,n,j)
             f(l2+1:l2+3,m,n,j)=matmul(mat,rhs)
           endif
        enddo
      enddo

!
      case default
        print*, "bc_nfc_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_nfc_x
!***********************************************************************
    subroutine bc_emf_x(f,df,dt_,topbot,j)
      character (len=bclen), intent (in) :: topbot
      real, dimension (mx,my,mz,mfarray), intent (inout) :: f
      real, dimension (mx,my,mz,mvar), intent(in) :: df
      real, intent(in) :: dt_
      integer, intent (in) :: j
      integer :: i
      select case (topbot)
!
      case ('bot')               ! bottom boundary
        print*, "bc_emf_x: ", topbot, " should be 'top'"
      case ('top')               ! top boundary
        if (llast_proc_x) then
          do i=1,nghost
            f(l2+i,m,n,j)=f(l2+i,m,n,j)+df(l2,m,n,j)*dt_ 
          enddo
        endif
      
      case default
        print*, "bc_emf_x: ", topbot, " should be 'top'"
      endselect
    endsubroutine bc_emf_x
!***********************************************************************
    subroutine bc_go_x(f,topbot,j,lforce_ghost)
!
!  Outflow boundary conditions.
!
!  If the velocity vector points out of the box, the velocity boundary
!  condition is set to 's', otherwise it is set to 'a'.
!  If 'lforce_ghost' is true, the boundary and ghost cell values are forced
!  to not point inwards. Otherwise the boundary value is forced to be 0.
!
!  14-jun-2011/axel: adapted from bc_outflow_z
!  17-sep-2012/piyali: adapted from bc_outflow_x
!
      character (len=bclen) :: topbot
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
      logical, optional :: lforce_ghost
!
      integer :: i, iy, iz
      logical :: lforce
!
      lforce = .false.
      if (present (lforce_ghost)) lforce = lforce_ghost
!
      select case (topbot)
!
!  Bottom boundary.
!
      case ('bot')
        do iy=1,my; do iz=1,mz
          if (f(l1,iy,iz,j)<0.0) then  ! 's'
            do i=1,nghost; f(l1-i,iy,iz,j)=+f(l1+i,iy,iz,j); enddo
          else                         ! 'a'
            do i=1,nghost; f(l1-i,iy,iz,j)=-f(l1+i,iy,iz,j); enddo
            f(l1,iy,iz,j)=0.0
          endif
          if (lforce) then
            do i = 1, nghost
              if (f(l1-i,iy,iz,j) > 0.0) f(l1-i,iy,iz,j) = 0.0
            enddo
          endif
        enddo; enddo
!
!  Top boundary.
!
      case ('top')
        if (llast_proc_x) then
          do iy=1,my
          do iz=1,mz
            if (f(l2,iy,iz,j)>0.0) then
              do i=1,nghost
                f(l2+i,iy,iz,j)=+f(l2-i,iy,iz,j)
              enddo
            else
              do i=1,nghost 
                f(l2+i,iy,iz,j)=-f(l2-i,iy,iz,j)
              enddo
              f(l2,iy,iz,j)=0.0
              if (lforce) then
                do i = 1, nghost
                  if (f(l2+i,iy,iz,j) < 0.0) f(l2+i,iy,iz,j) = 0.0
                enddo
              endif
            endif
          enddo
          enddo
        endif
!
!  Default.
!
      case default
        print*, "bc_go_x: ", topbot, " should be 'top' or 'bot'"
!
      endselect
!
    endsubroutine bc_go_x


!***********************************************************************
!
endmodule Special
