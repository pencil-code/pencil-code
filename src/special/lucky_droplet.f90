! $Id$
!
!  Solve the lucky droplet model for many realizations.
!  The different realizations correspond to "meshpoints".
!  To add the contributions for each step, we use the usual
!  time step in the Pencil Code, so t is just the step, and
!  the accumulated collision times (after 125 steps or so)
!  for all realizations at the same time are the values in
!  the f-array.
!
!  16-apr-20/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of special_dummies.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 1
! MAUX CONTRIBUTION 0
!
!***************************************************************
!
module Special
!
  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages
!
  implicit none
!
  include '../special.h'
!
! Declare index of variables
!
   integer :: ispecial=0
!
  ! input parameters
  real :: gam_lucky=fourthird
  character(len=50) :: init_qq='zero'
  namelist /special_init_pars/ &
    gam_lucky
!
  ! run parameters
  logical :: lMFT=.false., lrA=.false., lrB=.false.
  namelist /special_run_pars/ &
    gam_lucky, lMFT, lrA, lrB
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_tt1m=0 ! DIAG_DOC: $\langle T \rangle$
  integer :: idiag_qq1m=0 ! DIAG_DOC: $\langle \ln T \rangle$
  integer :: idiag_qq2m=0 ! DIAG_DOC: $\langle \ln T^2 \rangle$
  integer :: idiag_qq3m=0 ! DIAG_DOC: $\langle \ln T^3 \rangle$
  integer :: idiag_qq4m=0 ! DIAG_DOC: $\langle \ln T^4 \rangle$
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  19-feb-2019/axel: coded
!
      use FArrayManager
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
!  Set ichemistry to consecutive numbers nvar+1, nvar+2, ..., nvar+nchemspec.
!
      call farray_register_pde('special',ispecial)
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  19-feb-2019/axel: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
!
!  Initialize any module variables which are parameter dependent
!
      call keep_compiler_quiet(f)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  16-apr-2020/axel: coded
!
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
!
      intent(inout) :: f
!
!  Initial condition; same for every population.
!
      select case (init_qq)
        case ('nothing'); if (lroot) print*,'init_qq: nothing'
        case ('zero'); f(:,:,:,ispecial)=0.
!
        case default
          !
          !  Catch unknown values
          !
          if (lroot) print*,'init_qq: No such value for init_qq: ', trim(init_qq)
          call stop_it("")
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
      logical, dimension(npencils) :: lpencil_in
!
      call keep_compiler_quiet(lpencil_in)
!
    endsubroutine pencil_interdep_special
!***********************************************************************
    subroutine calc_pencils_special(f,p)
!
!  Calculate Hydro pencils.
!  Most basic pencils should come first, as others may depend on them.
!
!   24-nov-04/tony: coded
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
!  several precalculated Pencils of information are passed if for
!  efficiency.
!
!  16-apr-2020/axel: coded
!
      use General, only: random_number_wrapper
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: rrr, tauk, lamk
      real, dimension (nx) :: tt, qq
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
      call keep_compiler_quiet(p)
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt'
!
!  Define state vector
!
      if (t==tstart) then
        tt=.0
        qq=.0
      else
        tt=f(l1:l2,m,n,ispecial)
        qq=alog(tt)
      endif
!
!  Selection of different combinations of rA and rB
!
      if (lrA.and.lrB) then
        lamk=(t**onethird+1.)**2*(t**twothird-1.)
      elseif (lrA.and..not.lrB) then
        lamk=(t**onethird+1.)**2*t**twothird
      elseif (.not.lrA.and.lrB) then
        lamk=t**twothird*(t**twothird-1.)
      else
        lamk=t**gam_lucky
      endif
!
!  Produce exponentially distributed random numbers,
!  but can also do mean-field theory as a test.
!
      if (lMFT) then
        tauk=lamk
      else
        call random_number_wrapper(rrr)
        tauk=-alog(rrr)/lamk
      endif
      df(l1:l2,m,n,ispecial)=df(l1:l2,m,n,ispecial)+tauk
!
!  diagnostics
!
      if (ldiagnos) then
        if (idiag_tt1m/=0) call sum_mn_name(tt   ,idiag_tt1m)
        if (idiag_qq1m/=0) call sum_mn_name(qq   ,idiag_qq1m)
        if (idiag_qq2m/=0) call sum_mn_name(qq**2,idiag_qq2m)
        if (idiag_qq3m/=0) call sum_mn_name(qq**3,idiag_qq3m)
        if (idiag_qq4m/=0) call sum_mn_name(qq**4,idiag_qq4m)
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
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!  19-feb-2019/axel: coded
!
      use Diagnostics
      use FArrayManager, only: farray_index_append
      use Sub
!
!   SAMPLE IMPLEMENTATION
!
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
        idiag_tt1m=0; idiag_qq1m=0; idiag_qq2m=0; idiag_qq3m=0; idiag_qq4m=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'tt1m',idiag_tt1m)
        call parse_name(iname,cname(iname),cform(iname),'qq1m',idiag_qq1m)
        call parse_name(iname,cname(iname),cform(iname),'qq2m',idiag_qq2m)
        call parse_name(iname,cname(iname),cform(iname),'qq3m',idiag_qq3m)
        call parse_name(iname,cname(iname),cform(iname),'qq4m',idiag_qq4m)
      enddo
!
    endsubroutine rprint_special
!***********************************************************************
!
!********************************************************************
!************        DO NOT DELETE THE FOLLOWING       **************
!********************************************************************
!**  This is an automatically generated include file that creates  **
!**  copies dummy routines from nospecial.f90 for any Special      **
!**  routines not implemented in this file                         **
!**                                                                **
    include '../special_dummies.inc'
!********************************************************************
endmodule Special
