 ! $Id: global_gg_halo.f90,v 1.1 2002-07-18 23:17:50 dobler Exp $

module Global

!
!  A module container and access functions (`methods') for additional
!  variables which are globally needed --- here this is the gravity field
!  gg, and the halo profile for identifying a spherical halo region
!  outside a star.
!
  use Cparam
  use Mpicomm

  implicit none

  interface set_global
    module procedure set_global_vect
    module procedure set_global_scal
  endinterface

  interface get_global
    module procedure get_global_vect
    module procedure get_global_scal
  endinterface

  real, dimension (mx,my,mz,3) :: gg=1.
  real, dimension (mx,my,mz)   :: halo=0.

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
      select case(label)

      case ('gg')
        gg(l1:l2,m,n,1:3) = var

      case default
        if (lroot) print*, 'No such value for label', trim(label)
        call stop_it('SET_GLOBAL_VECT')

      endselect
!
    endsubroutine set_global_vect
!***********************************************************************
    subroutine set_global_scal(var,m,n,label)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!
!  18-jul-02/wolf coded
!
      real, dimension(nx) :: var
      integer :: m,n
      character (len=*) ::label
!
      select case(label)

      case ('halo')
        halo(l1:l2,m,n) = var

      case default
        if (lroot) print*, 'No such value for label', trim(label)
        call stop_it('SET_GLOBAL_VECT')

      endselect
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
      character (len=*) ::label
!
      select case(label)

      case ('gg')
        var = gg(l1:l2,m,n,1:3)

      case default
        if (lroot) print*, 'No such value for label', trim(label)
        call stop_it('SET_GLOBAL_VECT')

      endselect
!
    endsubroutine get_global_vect
!***********************************************************************
    subroutine get_global_scal(var,m,n,label)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!
!  18-jul-02/wolf coded
!
      real, dimension(nx) :: var
      integer :: m,n
      character (len=*) ::label
!
      select case(label)

      case ('halo')
        var = halo(l1:l2,m,n)

      case default
        if (lroot) print*, 'No such value for label', trim(label)
        call stop_it('SET_GLOBAL_VECT')

      endselect
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
      use IO
!
      call output(trim(directory)//'/halo.dat',halo,1)
      call output(trim(directory)//'/gg.dat'  ,gg  ,3)
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
      call input(trim(directory)//'/halo.dat',halo,1,0)
      call input(trim(directory)//'/gg.dat'  ,gg  ,3,0)
!
    endsubroutine rglobal
!***********************************************************************

endmodule Global
