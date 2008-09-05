! $Id$

module Global

!
!  A module container for additional variables which are globally needed
!  --- here this is spherical radius (for calculating gravity and heating
!  /cooling functions)
!  Put into a module, so one can easily switch.
!  NB: These variables use half as much memory as a new variable, so
!  keep their number at minimum.
!
  use Cparam

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

  real, dimension (mx,my,mz) :: rr

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
!  set (m-n)-pencil of global vector variable identified by LABEL
!  [dummy routine]
!
      use Mpicomm, only: stop_it
!
      integer :: length
      real, dimension(length,3) :: var
      integer :: m,n
      character (len=*) :: label
!
      call stop_it("GLOBAL_RR: set_global_vect not implemented")
      if (ip == 0) print*, var(1,1),m,n,label ! keep compiler quiet
!
    endsubroutine set_global_vect
!***********************************************************************
    subroutine set_global_scal(var,m,n,label,length)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!  [dummy routine]
!
      use Mpicomm, only: stop_it
!
      integer :: length
      real, dimension(length) :: var
      integer :: m,n
      character (len=*) :: label
!
      call stop_it("GLOBAL_RR: set_global_scal not implemented")
      if (ip == 0) print*, var, m, n, label ! keep compiler quiet
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
      if (NO_WARN) print*, l, var(1), m, n, label ! keep compiler quiet
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
      if (NO_WARN) print*, label ! keep compiler quiet
!
    endsubroutine reset_global
!***********************************************************************
    subroutine get_global_vect(var,m,n,label)
!
!  Get (m,n)-pencil of the global vector variable identified by LABEL.
!  [dummy routine]
!
      use Mpicomm, only: stop_it
!
      real, dimension(nx,3) :: var
      integer :: m,n
      character (len=*) :: label
!
      call stop_it("GLOBAL_RR: get_global_vect not implemented")
      if (ip == 0) print*, var(1,1),m,n,label ! keep compiler quiet
!
    endsubroutine get_global_vect
!***********************************************************************
    subroutine get_global_scal(var,m,n,label)
!
!  Get (m,n)-pencil of the global scalar variable identified by LABEL.
!  [dummy routine]
!
      use Mpicomm, only: stop_it
!
      real, dimension(nx) :: var
      integer :: m,n
      character (len=*) :: label
!
      call stop_it("GLOBAL_RR: get_global_scla not implemented")
      if (ip == 0) print*, var(1),m,n,label ! keep compiler quiet
!
    endsubroutine get_global_scal
!***********************************************************************
    subroutine get_global_vect_point(var,l,m,n,label)
!
!  Get (l,m,n)-point of the global vector variable identified by LABEL.
!  [dummy routine]
!
      use Mpicomm, only: stop_it
!
      real, dimension(3) :: var
      integer :: l,m,n
      character (len=*) :: label
!
      call stop_it("GLOBAL_RR: get_global_vect not implemented")
      if (NO_WARN) print*, var(1),l,m,n,label ! keep compiler quiet
!
    endsubroutine get_global_vect_point
!***********************************************************************
    subroutine get_global_scal_point(var,l,m,n,label)
!
!  Get (l,m,n)-pointof the global scalar variable identified by LABEL.
!  [dummy routine]
!
      use Mpicomm, only: stop_it
!
      real :: var
      integer :: l,m,n
      character (len=*) :: label
!
      call stop_it("GLOBAL_RR: get_global_scla not implemented")
      if (ip == 0) print*, var(1),l,m,n,label ! keep compiler quiet
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
      use Mpicomm
      use IO
!
      call output(trim(directory)//'/rr.dat',rr,1)
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
      use Mpicomm
      use IO
!
      call input(trim(directory)//'/rr.dat',rr,1,0)
!
    endsubroutine rglobal
!***********************************************************************

endmodule Global
