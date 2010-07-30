 ! $Id$

module Global

!
!  A module container and access functions (`methods') for additional
!  variables which are globally needed --- here this is the gravity field
!  gg, and an external magnetic field (such that B=curl(A)+B_ext_pot)
!

  use Cparam
  use Mpicomm

  implicit none

  interface set_global
    module procedure set_global_vect
    !module procedure set_global_scal
  endinterface

  interface get_global
    module procedure get_global_vect
    !module procedure get_global_scal
  endinterface

  real, dimension (mx,my,mz,3) :: gg
  real, dimension (mx,my,mz,3) :: B_ext_pot

  contains

!***********************************************************************
    subroutine set_global_vect(var,m,n,label)
!
!  set (m,n)-pencil of the global vector variable identified by LABEL
!
!  18-jul-02/wolf coded
!
!      use Cparam
!
      real, dimension(nx,3) :: var
      integer :: m,n
      character (len=*) ::label
!
      select case (label)

      case ('gg')
        gg(l1:l2,m,n,1:3) = var

      case ('B_ext_pot')
        B_ext_pot(l1:l2,m,n,1:3) = var

      case default
        if (lroot) print*, 'set_global_vect: No such value for label ' &
                         , trim(label)
        call stop_it('set_global_vect')

      endselect
!
    endsubroutine set_global_vect
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
      character (len=*) ::label
!
      select case (label)

      case ('gg')
        var = gg(l1:l2,m,n,1:3)

      case ('B_ext_pot')
        var = B_ext_pot(l1:l2,m,n,1:3)

      case default
        if (lroot) print*, 'get_global_vect: No such value for label ', &
                           trim(label)
        call stop_it('get_global_vect')

      endselect
!
    endsubroutine get_global_vect
!***********************************************************************
    subroutine wglobal()
!
!  write global variables
!
!  10-jan-02/wolf: coded
!
      use Cdata
      use IO
!
      call output(trim(directory)//'/gg.dat'  ,gg  ,3)
      call output(trim(directory)//'/B_ext_pot.dat'  ,B_ext_pot  ,3)
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
      use IO
!
      call input(trim(directory)//'/B_ext_pot.dat'  ,B_ext_pot  ,3,0)
      call input(trim(directory)//'/gg.dat'  ,gg  ,3,0)
!
    endsubroutine rglobal
!***********************************************************************

endmodule Global
