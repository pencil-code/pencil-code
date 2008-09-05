 ! $Id$

module Global

!
!  A module container and access functions (`methods') for additional
!  variables which are globally needed --- here ee_ext, bb_ext, and jj_ext.
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

  real, dimension (mx,my,mz,3) :: eee_ext,bbb_ext,jjj_ext

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
      if (NO_WARN) print*,var,m,n,label
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
      select case(label)

      case ('ee_ext')
        var = eee_ext(l1:l2,m,n,:)

      case ('bb_ext')
        var = bbb_ext(l1:l2,m,n,:)

      case ('jj_ext')
        var = jjj_ext(l1:l2,m,n,:)

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
      call input(trim(directory)//'/ee_ext.dat',eee_ext,3,0)
      call input(trim(directory)//'/bb_ext.dat',bbb_ext,3,0)
!     call input(trim(directory)//'/jj_ext.dat',jjj_ext,3,0)
!
    endsubroutine rglobal
!***********************************************************************

endmodule Global
