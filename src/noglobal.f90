! $Id: noglobal.f90,v 1.5 2004-06-10 00:47:18 dobler Exp $

module Global

!
!  A dummy module for the (lucky) case when we need no global variables.
!
  use Cparam

  implicit none

  interface set_global
    module procedure set_global_vect
    module procedure set_global_scal
  endinterface

  interface get_global
    module procedure get_global_vect
    module procedure get_global_scal
  endinterface

  contains

!***********************************************************************
    subroutine set_global_vect(var,m,n,label)
!
!  set (m-n)-pencil of global vector variable identified by LABEL
!  [dummy routine]
!
!  18-jul-02/wolf coded
!
      real, dimension(nx,3) :: var
      integer :: m,n
      character (len=*) :: label
!
      if (ip == 0) print*, var(1,1),m,n,label ! keep compiler quiet
!
    endsubroutine set_global_vect
!***********************************************************************
    subroutine set_global_scal(var,m,n,label)
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
      if (ip == 0) print*, var(1),m,n,label ! keep compiler quiet
!
    endsubroutine set_global_scal
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
      if (ip == 0) print*, var(1,1),m,n,label ! keep compiler quiet
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
      if (ip == 0) print*, var(1),m,n,label ! keep compiler quiet
!
    endsubroutine get_global_scal
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
