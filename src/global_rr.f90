! $Id: global_rr.f90,v 1.4 2004-06-10 00:46:49 dobler Exp $

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

  interface set_global
    module procedure set_global_vect
    module procedure set_global_scal
  endinterface

  interface get_global
    module procedure get_global_vect
    module procedure get_global_scal
  endinterface

  real, dimension (mx,my,mz) :: rr

  contains

!***********************************************************************
    subroutine set_global_vect(var,m,n,label)
!
!  set (m-n)-pencil of global vector variable identified by LABEL
!  [dummy routine]
!
      use Mpicomm, only: stop_it
!
      real, dimension(nx,3) :: var
      integer :: m,n
      character (len=*) :: label
!
      call stop_it("GLOBAL_RR: set_global_vect not implemented")
      if (ip == 0) print*, var(1,1),m,n,label ! keep compiler quiet
!
    endsubroutine set_global_vect
!***********************************************************************
    subroutine set_global_scal(var,m,n,label)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!  [dummy routine]
!
      use Mpicomm, only: stop_it
!
      real, dimension(nx) :: var
      integer :: m,n
      character (len=*) :: label
!
      call stop_it("GLOBAL_RR: set_global_scal not implemented")
      if (ip == 0) print*, var(1),m,n,label ! keep compiler quiet
!
    endsubroutine set_global_scal
!***********************************************************************
    subroutine get_global_vect(var,m,n,label)
!
!  get (m,n)-pencil of the global vector variable identified by LABEL
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
!  set (m,n)-pencil of the global scalar variable identified by LABEL
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
