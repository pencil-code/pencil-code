!$Id: global_gg_hcond.f90,v 1.1 2006-11-13 14:38:56 dintrans Exp $

module Global

!
!  13-Nov-2006/bdintrans: coded
!  Global module for the gravity field gg, hcond and grad(log hcond).
!  To be used in connection with the 'lhcond_global' flag settled in
!  the 'entropy_run_pars' namelist
!
!
  use Cparam
  use Mpicomm

  implicit none

  include 'global.h'

  interface set_global
    module procedure set_global_vect
    module procedure set_global_scal
  endinterface

  interface set_global_point
    module procedure set_global_vect_point
    module procedure set_global_scal_point
  endinterface

  interface get_global
    module procedure get_global_vect
    module procedure get_global_scal
  endinterface


  interface get_global_point
    module procedure get_global_vect_point
    module procedure get_global_scal_point
  endinterface
!
!  we could (and should in some sense) initialize gg and halo; however,
!  with the Intel compiler this produces significantly larger
!  executables, compiles far longer, and (for 256^3) results in the
!  compiler message `Size of initialised array entity exceeds
!  implementation limit'
!
  real, dimension (mx,my,mz,3) :: gg,glhc
  real, dimension (mx,my,mz)   :: hcond

  contains

!***********************************************************************
    subroutine register_global()
!
!  Register Global module.
!
!  13-jun-05/anders: coded
!
      use Cdata, only: lglobal
!
      lglobal=.true.
!
    endsubroutine register_global
!***********************************************************************
    subroutine set_global_vect(var,m,n,label,length)
!
!  set (m,n)-pencil of the global vector variable identified by LABEL
!
!  18-jul-02/wolf coded
!
      use Mpicomm, only: stop_it
!
      integer :: length
      real, dimension(length,3) :: var
      integer :: m,n
      character (len=*) ::label
!
      if (length/=nx) then
        print*, 'set_global_vect: only nx pencils allowed!'
        call stop_it('set_global_vect')
      endif
!
      select case(label)

      case ('gg')
        gg(l1:l2,m,n,1:3) = var

      case ('glhc')
        glhc(l1:l2,m,n,1:3) = var

      case default
        if (lroot) print*, 'set_global_vect: No such value for label', trim(label)
        call stop_it('set_global_vect')

      endselect
!
    endsubroutine set_global_vect
!***********************************************************************
    subroutine set_global_scal(var,m,n,label,length)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!
!  18-jul-02/wolf coded
!
      use Mpicomm, only: stop_it
!
      integer :: length
      real, dimension(length) :: var
      integer :: m,n
      character (len=*) ::label
!
      if (length/=nx) then
        print*, 'set_global_scal: only nx pencils allowed!'
            call stop_it('set_global_scal')
      endif
!
      select case(label)

      case ('hcond')
        hcond(l1:l2,m,n) = var

      case default
        if (lroot) print*, 'set_global_scal: No such value for label', trim(label)
        call stop_it('set_global_scal')

      endselect
!
    endsubroutine set_global_scal
!***********************************************************************
    subroutine set_global_vect_point(var,l,m,n,label)
!
!  set point of global vector variable identified by LABEL
!
!  20-jun-05/anders: dummy
!
      real, dimension(3) :: var
      integer :: l,m,n
      character (len=*) :: label
!
      if (NO_WARN) print*, l, var(1), m, n, label ! keep compiler quiet
!
    endsubroutine set_global_vect_point
!***********************************************************************
    subroutine set_global_scal_point(var,l,m,n,label)
!
!  set point of global scalar variable identified by LABEL
!
!  20-jun-05/anders: dummy
!
      real :: var
      integer :: l,m,n
      character (len=*) :: label
!
      if (NO_WARN) print*, l, var, m, n, label ! keep compiler quiet
!
    endsubroutine set_global_scal_point
!***********************************************************************
    subroutine reset_global(label)
!
!  reset global variable identified by LABEL
!
!  20-jun-05/anders: dummy
!
      character (len=*) :: label
!
      if (ip == 0) print*, label ! keep compiler quiet
!
    endsubroutine reset_global
!***********************************************************************
    subroutine get_global_vect(var,m,n,label)
!
!  Get (m,n)-pencil of the global vector variable identified by LABEL.
!
!  18-jul-02/wolf coded
!
      real, dimension(nx,3) :: var
      integer :: m,n
      character (len=*) ::label
!
      select case(label)

      case ('gg')
        var = gg(l1:l2,m,n,1:3)

      case ('glhc')
        var = glhc(l1:l2,m,n,1:3)

      case default
        if (lroot) print*, 'get_global_vect: No such value for label', trim(label)
        call stop_it('get_global_vect')

      endselect
!
    endsubroutine get_global_vect
!***********************************************************************
    subroutine get_global_scal(var,m,n,label)
!
!  Get (m,n)-pencil of the global scalar variable identified by LABEL.
!
!  18-jul-02/wolf coded
!
      real, dimension(nx) :: var
      integer :: m,n
      character (len=*) ::label
!
      select case(label)

      case ('hcond')
        var = hcond(l1:l2,m,n)

      case default
        if (lroot) print*, 'get_global_scal: No such value for label', trim(label)
        call stop_it('get_global_scal')

      endselect
!
    endsubroutine get_global_scal
!***********************************************************************
    subroutine get_global_vect_point(var,l,m,n,label)
!
!  Get (l,m,n)-point of the global vector variable identified by LABEL.
!
!  15-sep-05/anders: adapted
!
      real, dimension(3) :: var
      integer :: l,m,n
      character (len=*) ::label
!
      select case(label)

      case ('gg')
        var = gg(l,m,n,1:3)

      case ('glhc')
        var = glhc(l,m,n,1:3)

      case default
        if (lroot) print*, &
            'get_global_vect_point: No such value for label', trim(label)
        call stop_it('get_global_vect_point')

      endselect
!
    endsubroutine get_global_vect_point
!***********************************************************************
    subroutine get_global_scal_point(var,l,m,n,label)
!
!  Get (l,m,n)-point of the global scalar variable identified by LABEL
!
!  15-sep-05/anders: adapted
!
      real :: var
      integer :: l,m,n
      character (len=*) ::label
!
      select case(label)

      case ('hcond')
        var = hcond(l,m,n)

      case default
        if (lroot) print*, &
            'get_global_scal_point: No such value for label', trim(label)
        call stop_it('get_global_scal_point')

      endselect
!
    endsubroutine get_global_scal_point
!***********************************************************************
    subroutine global_derivs(m,n,label,der6)
!
!  take any derivative of global scalar variable.
!
!  13-jun-05/anders: dummy
!
      real, dimension (nx), optional :: der6
      integer :: m,n
      character (len=*) ::label
!
      if (NO_WARN) print*, m, n, label
!
    endsubroutine global_derivs
!***********************************************************************
    subroutine wglobal()
!
!  write global variables
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use IO, only: output
!!      use IO
!
!  No need to read/write them as run.f90 will recalculate anyway
!
      if (ip <= 4) then
          call output(trim(directory)//'/gg.dat'  ,gg  ,3)
      endif
!
    endsubroutine wglobal
!***********************************************************************
    subroutine rglobal()
!
!  read global variables
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use IO, only: input
!!      use IO
!
!  No need to read/write them as run.f90 will recalculate anyway
!      call input(trim(directory)//'/halo.dat',halo,1,0)
!      call input(trim(directory)//'/gg.dat'  ,gg  ,3,0)
!
    endsubroutine rglobal
!***********************************************************************

endmodule Global
