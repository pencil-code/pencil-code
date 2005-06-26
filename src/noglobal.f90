! $Id: noglobal.f90,v 1.6 2005-06-26 17:34:13 eos_merger_tony Exp $

module Global

!
!  A dummy module for the (lucky) case when we need no global variables.
!
  use Cparam

  implicit none

  include 'global.inc'

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

  contains

!***********************************************************************
    subroutine register_global()
!
!  Register Global module.
!
!  13-jun-05/anders: dummy
!
    endsubroutine register_global
!***********************************************************************
    subroutine set_global_vect(var,m,n,label,length)
!
!  set (m-n)-pencil of global vector variable identified by LABEL
!  [dummy routine]
!
!  18-jul-02/wolf coded
!
      integer :: length
      real, dimension(length,3) :: var
      integer :: m,n
      character (len=*) :: label
!
      if (NO_WARN) print*, length, var(1,1), m, n, label ! keep compiler quiet
!
    endsubroutine set_global_vect
!***********************************************************************
    subroutine set_global_scal(var,m,n,label,length)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!  [dummy routine]
!
!  18-jul-02/wolf coded
!
      integer :: length
      real, dimension(nx) :: var
      integer :: m,n
      character (len=*) :: label
!
      if (NO_WARN) print*, length, var(1), m, n, label ! keep compiler quiet
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
      if (NO_WARN) print*, label ! keep compiler quiet
!
    endsubroutine reset_global
!***********************************************************************
    subroutine get_global_vect(var,m,n,label)
!
!  get (m,n)-pencil of the global vector variable identified by LABEL
!  [dummy routine]
!
!  18-jul-02/wolf coded
!
!      use Cparam
!
      real, dimension(nx,3) :: var
      integer :: m,n
      character (len=*) :: label
!
      if (NO_WARN) print*, var(1,1),m,n,label ! keep compiler quiet
!
    endsubroutine get_global_vect
!***********************************************************************
    subroutine get_global_scal(var,m,n,label)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!  [dummy routine]
!
!  18-jul-02/wolf coded
!
!      use Cparam
!
      real, dimension(nx) :: var
      integer :: m,n
      character (len=*) :: label
!
      if (NO_WARN) print*, var(1),m,n,label ! keep compiler quiet
!
    endsubroutine get_global_scal
!***********************************************************************
    subroutine global_derivs(m,n,label,der6)
!
!  take any derivative of global scalar variable.
!
!  13-jun-05/anders: dummy
!
      use Cdata, only: lroot
      use General, only: warning
!
      real, dimension (nx), optional :: der6
      integer :: m,n
      character (len=*) ::label
!
      call warning('global_derivs')
      if (lroot) print*, 'global_derivs: not implemented'
!
      if (NO_WARN) print*, m, n, label
!
    endsubroutine global_derivs
!***********************************************************************
    subroutine wglobal()
!
!  write global variables
!
    endsubroutine wglobal
!***********************************************************************
    subroutine rglobal()
!
!  read global variables
!
    endsubroutine rglobal
!***********************************************************************

endmodule Global
