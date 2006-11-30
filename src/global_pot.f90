! $Id: global_pot.f90,v 1.6 2006-11-30 09:03:35 dobler Exp $

module Global

!
!  A module container for additional variables which are globally needed
!  --- here this is the (negative) gravity potential.
!  Put into a module, so one can easily switch.
!  NB: These variables use half as much memory as a new variable, so
!  keep their number at minimum.
!
  use Cparam

  implicit none

  include 'global.h'

  real, dimension (mx,my,mz) :: m_pot

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
!  18-jul-02/wolf coded
!
      real, dimension(nx,3) :: var
      integer :: m,n
      character (len=labellen) ::label
!
      if (ip == 0) print*, var(1,1),m,n,label ! keep compiler quiet
!
    endsubroutine set_global_vect
!***********************************************************************
    subroutine set_global_scal(var,m,n,label)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!
!  18-jul-02/wolf coded
!
!      use Cparam
!
      real, dimension(nx) :: var
      integer :: m,n
      character (len=labellen) ::label
!
      if (ip == 0) print*, var(1),m,n,label ! keep compiler quiet
!
    endsubroutine set_global_scal
!***********************************************************************
    subroutine get_global_vect(var,m,n,label)
!
!  set (m,n)-pencil of the global vector variable identified by LABEL
!
!  18-jul-02/wolf coded
!
!      use Cparam
!
      real, dimension(nx,3) :: var
      integer :: m,n
      character (len=labellen) ::label
!
      if (ip == 0) print*, var(1,1),m,n,label ! keep compiler quiet
!
    endsubroutine get_global_vect
!***********************************************************************
    subroutine get_global_scal(var,m,n,label)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!
!  18-jul-02/wolf coded
!
!      use Cparam
!
      real, dimension(nx) :: var
      integer :: m,n
      character (len=labellen) ::label
!
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
      call output(trim(directory)//'/global.dat',m_pot,1)
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
      call input(trim(directory)//'/global.dat',m_pot,1,0)
!
    endsubroutine rglobal
!***********************************************************************

endmodule Global
