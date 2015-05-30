! $Id$
!
!  Solve for a set of two ODEs, used to test time step
!
!  19-apr-09/axel: adapted from nospecial.f90
!
!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! CPARAM logical, parameter :: lspecial = .true.
!
! MVAR CONTRIBUTION 3
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Special

  use Cparam
  use Cdata
  use General, only: keep_compiler_quiet
  use Messages

  implicit none

  include '../special.h'

!
! Declare index of variables
!
   integer :: ispecial=0,ispecial1=0,ispecial2=0,ispecial3=0

  ! input parameters
  real :: bet,gam,rho,xxini,yyini,zzini
  real, dimension (ninit) :: ampl=0.
  character (len=labellen), dimension(ninit) :: init='nothing'
  namelist /special_init_pars/ &
    init,ampl,bet,gam,rho,xxini,yyini,zzini

  ! run parameters
  namelist /special_run_pars/ &
    bet,gam,rho
!
! other variables (needs to be consistent with reset list below)
!
  integer :: idiag_xx=0,idiag_yy=0,idiag_zz=0
!
  contains

!***********************************************************************
    subroutine register_special()
!
!  Configure pre-initialised (i.e. before parameter read) variables
!  which should be know to be able to evaluate
!
!  6-oct-03/tony: coded
!
      use FArrayManager
!
!  Identify CVS/SVN version information.
!
      if (lroot) call svn_id( &
           "$Id$")
!
      call farray_register_pde('spec_3vec',ispecial,vector=3)
!
      ispecial1=ispecial
      ispecial2=ispecial+1
      ispecial3=ispecial+2
!
    endsubroutine register_special
!***********************************************************************
    subroutine initialize_special(f)
!
!  called by run.f90 after reading parameters, but before the time loop
!
!  06-oct-03/tony: coded
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
!  06-oct-2003/tony: coded
!
      use Initcond
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      integer :: j
!
      intent(inout) :: f
!
!  initial condition
!
      do j=1,ninit
        select case (init(j))
        case ('zero'); f(:,:,:,ispecial1:ispecial3)=0.
        case ('nothing'); if (lroot) print*,'init_special: nothing'
        case ('gaussian-noise'); call gaunoise(ampl(j),f,ispecial1,ispecial3)
        case ('set')
          f(:,:,:,ispecial1)=f(:,:,:,ispecial1)+xxini
          f(:,:,:,ispecial2)=f(:,:,:,ispecial2)+yyini
          f(:,:,:,ispecial3)=f(:,:,:,ispecial3)+zzini
        case default
          if (lroot) print*,'init_special: No such value: ',trim(init(j))
          call stop_it("")
        endselect
      enddo
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
!   06-oct-03/tony: coded
!
      use Diagnostics
      use Mpicomm
      use Sub
!
      real, dimension (mx,my,mz,mfarray) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx) :: xx,yy,zz
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
!  Solve the Lorenz equations
!
        xx=f(l1:l2,m,n,ispecial1)
        yy=f(l1:l2,m,n,ispecial2)
        zz=f(l1:l2,m,n,ispecial3)
!
        df(l1:l2,m,n,ispecial1)=df(l1:l2,m,n,ispecial1)+gam*(yy-xx)
        df(l1:l2,m,n,ispecial2)=df(l1:l2,m,n,ispecial2)+rho*xx-yy-xx*zz
        df(l1:l2,m,n,ispecial3)=df(l1:l2,m,n,ispecial3)+xx*yy-bet*zz
!
!  diagnostics
!
      if (ldiagnos) then
        if (lroot.and.m==mpoint.and.n==npoint) then
          if (idiag_xx/=0) call save_name(xx(lpoint-nghost),idiag_xx)
          if (idiag_yy/=0) call save_name(yy(lpoint-nghost),idiag_yy)
          if (idiag_zz/=0) call save_name(zz(lpoint-nghost),idiag_zz)
        endif
      endif
!
    endsubroutine dspecial_dt
!***********************************************************************
    subroutine read_special_init_pars(iostat)
!
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include "../parallel_unit.h"
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
      use File_io, only: get_unit
!
      integer, intent(out) :: iostat
      include "../parallel_unit.h"
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
!   06-oct-03/tony: coded
!
      use Diagnostics
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
        idiag_xx=0; idiag_yy=0; idiag_zz=0
      endif
!
      do iname=1,nname
        call parse_name(iname,cname(iname),cform(iname),'xx',idiag_xx)
        call parse_name(iname,cname(iname),cform(iname),'yy',idiag_yy)
        call parse_name(iname,cname(iname),cform(iname),'zz',idiag_zz)
      enddo
!
!  write column where which variable is stored
!
      if (lwr) then
        write(3,*) 'i_xx=',idiag_xx
        write(3,*) 'i_yy=',idiag_yy
        write(3,*) 'i_zz=',idiag_zz
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
    include '../special_dummies.inc'
!********************************************************************

endmodule Special
