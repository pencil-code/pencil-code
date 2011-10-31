! $Id$
!
!  Solve the Van der Pol oscillator equations
!
!  9-oct-09/tgastine: adapted from oscillation_0D.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 2
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special
!
  use Cdata
  use Cparam
  use Messages
  use Sub, only: keep_compiler_quiet
!
  implicit none
!
  include '../special.h'
!
! Declare index of variables
!
  integer :: ispecial=0,ispecial1=0,ispecial2=0
  real :: u1ini=1.5e-2,u2ini=0
  real :: tau=3.e-3, om1=3.8, finalamp=0.025
  real :: om_forc=3.8, amp_forc=0.1
  character(len=50) :: init='zero'
!
! input parameters
!
  namelist /special_init_pars/ init,u1ini,u2ini
!
! run parameters
!
  namelist /special_run_pars/ tau, finalamp, om1, om_forc, amp_forc
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_u1=0,idiag_u2=0
!
  contains
!
!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      use FArrayManager, only: farray_register_pde
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_pde('ispecial',ispecial,vector=2)
!
      ispecial1=ispecial
      ispecial2=ispecial+1
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f,lstarting)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mfarray) :: f
      logical :: lstarting
!
!  Initialize any module variables which are parameter dependent
!
      call keep_compiler_quiet(f)
      call keep_compiler_quiet(lstarting)
!
    endsubroutine initialize_special
!***********************************************************************
    subroutine init_special(f)
!
!  initialise special condition; called from start.f90
!  06-oct-2003/tony: coded
!
      use Mpicomm, only: stop_it
!
      real, dimension (mx,my,mz,mfarray) :: f
      intent(inout) :: f
!
!  initial condition
!
      select case (init)
        case ('nothing')
          if (lroot) print*,'init_special: nothing'
        case ('zero')
          f(:,:,:,ispecial1)=0.
        case ('set')
          f(:,:,:,ispecial1)=u1ini
          f(:,:,:,ispecial2)=u2ini
        case default
          if (lroot) print*,'init_special: No such value for init: ', trim(init)
          call stop_it("")
      endselect
!
    endsubroutine init_special
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
!   06-oct-03/tony: coded
!
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: u1,u2
      type (pencil_case) :: p
!
      intent(in) :: f,p
      intent(inout) :: df
!
!  identify module and boundary conditions
!
      if (headtt.or.ldebug) print*,'dspecial_dt: SOLVE dSPECIAL_dt', tau, finalamp
!
!  Solve the equations du1/dt = u2
!                      du2/dt = 2*tau*(1-u1**2/b**2)*u2 - om1**2*u1
!
        u1=f(l1:l2,m,n,ispecial1)
        u2=f(l1:l2,m,n,ispecial2)
!
        df(l1:l2,m,n,ispecial1)=df(l1:l2,m,n,ispecial1)+u2
        df(l1:l2,m,n,ispecial2)=df(l1:l2,m,n,ispecial2)+ &
                          2*tau*(1.-u1**2/finalamp**2)*u2 -om1**2*u1 + &
                          om1**2 * amp_forc*cos(om_forc*t)
!
!  diagnostics
!
      if (ldiagnos) then
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_u1/=0) call save_name(u1(lpoint-nghost),idiag_u1)
          if (idiag_u2/=0) call save_name(u2(lpoint-nghost),idiag_u2)
        endif
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
!  read namelist
!
      if (present(iostat)) then
        read(unit,NML=special_init_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_init_pars,ERR=99)
      endif
!
99    return
    endsubroutine read_special_init_pars
!***********************************************************************
    subroutine write_special_init_pars(unit)
!
      integer, intent(in) :: unit
!
!  write name list
!
      write(unit,NML=special_init_pars)
!
    endsubroutine write_special_init_pars
!***********************************************************************
    subroutine read_special_run_pars(unit,iostat)
!
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
!  write name list
!
      if (present(iostat)) then
        read(unit,NML=special_run_pars,ERR=99, IOSTAT=iostat)
      else
        read(unit,NML=special_run_pars,ERR=99)
      endif
!
99  endsubroutine read_special_run_pars
!***********************************************************************
    subroutine write_special_run_pars(unit)
      integer, intent(in) :: unit
!
!  write name list
!
      write(unit,NML=special_run_pars)
!
    endsubroutine write_special_run_pars
!***********************************************************************
    subroutine rprint_special(lreset,lwrite)
!
!  reads and registers print parameters relevant to special
!
!   06-oct-03/tony: coded
!
      use Diagnostics
      use Sub
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
        idiag_u1=0; idiag_u2=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'u1',idiag_u1)
        call parse_name(iname,cname(iname),cform(iname),'u2',idiag_u2)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_u1=',idiag_u1
        write(3,*) 'i_u2=',idiag_u2
      endif
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
    include 'special_dummies.inc'
!********************************************************************

endmodule Special

