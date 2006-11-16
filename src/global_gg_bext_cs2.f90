 ! $Id: global_gg_bext_cs2.f90,v 1.3 2006-11-16 19:58:18 mee Exp $

module Global

!
!  A module container and access functions (`methods') for additional
!  variables which are globally needed --- here this is the gravity field
!  gg, an external magnetic field (such that B=curl(A)+B_ext_pot), and
!  the local sound speed profile cs2 falling with r-1 for the global-disc
!  comparison project setup (i.e. extremely NOT general!!!!)
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
  real, dimension (mx,my,mz,3) :: gg
  real, dimension (mx,my,mz,3) :: B_ext_pot
  real, dimension (mx,my,mz) :: cs2

  contains

!***********************************************************************
    subroutine register_global()
!
!  Register Global module.
!
!  13-jun-05/anders: coded
!
      use Cdata, only: lglobal, lglobal_nolog_density
!
      lglobal=.true.
      lglobal_nolog_density=.true.
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

     case ('B_ext_pot')
        B_ext_pot(l1:l2,m,n,1:3) = var

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
!  13-jun-05/anders: adapted
!
      integer :: length
      real, dimension(length) :: var
      integer :: m,n
      character (len=*) ::label
!
      select case(label)
!
      case ('cs2')
         if (length==nx) then
            cs2(l1:l2,m,n) = var
         elseif (length==mx) then
            cs2(:,m,n) = var
         endif
!
      case default
        if (lroot) &
            print*, 'set_global_scal: No such value for label', trim(label)
        call stop_it('set_global_scal')
!
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
!
      case ('gg')
         var = gg(l1:l2,m,n,1:3)
         !
      case ('B_ext_pot')
         var = B_ext_pot(l1:l2,m,n,1:3)
!
      case default
         if (lroot) print*, 'get_global_vect: No such value for label', trim(label)
         call stop_it('get_global_vect')

      endselect
!
    endsubroutine get_global_vect
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
!
      case ('gg')
         var = gg(l,m,n,1:3)
!
      case ('B_ext_pot')
         var = B_ext_pot(l,m,n,1:3)
!
      case default
        if (lroot) print*, &
            'get_global_vect_point: No such value for label', trim(label)
        call stop_it('get_global_vect_point')
!
      endselect
!
    endsubroutine get_global_vect_point
!********************************************************************
    subroutine get_global_scal_point(var,l,m,n,label)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!
!  04-oct-05/tony: adapted
!
      real :: var
      integer :: l,m,n
      character (len=*) ::label
!
      select case(label)
!
      case ('cs2')
        var = cs2(l,m,n)
!
      case default
        if (lroot) print*, 'get_global_scal: No such value for label', trim(label)
        call stop_it('get_global_scal')
!
      endselect
!
    endsubroutine get_global_scal_point
!***********************************************************************
   subroutine get_global_scal(var,m,n,label)
!
!  set (m,n)-pencil of the global scalar variable identified by LABEL
!
!  13-jun-05/anders: adapted
!
      real, dimension(nx) :: var
      integer :: m,n
      character (len=*) ::label
!
      select case(label)
!
      case ('cs2')
         var = cs2(l1:l2,m,n)
!
      case default
        if (lroot) print*, 'get_global_scal: No such value for label', trim(label)
        call stop_it('get_global_scal')
!
      endselect
!
    endsubroutine get_global_scal
!***********************************************************************
    subroutine global_derivs(m,n,label,der6)
!
!  take any derivative of global scalar variable.
!
!  13-jun-05/anders: coded
!
      use Sub, only: del6_other
!
      real, dimension (nx), optional :: der6
      integer :: m,n
      character (len=*) ::label
!
      select case(label)
!
      case ('cs2')
        if (present(der6)) call del6_other(cs2,der6)
!
      case default
        if (lroot) &
            print*, 'global_derivs: No such value for label', trim(label)
        call stop_it('global_derivs')
!
      endselect
!
    endsubroutine global_derivs
!***********************************************************************
    subroutine wglobal()
!
!  write global variables
!
!  10-jan-02/wolf: coded
!
      use Cdata, only: directory
      use IO, only: output
!
      call output(trim(directory)//'/cs2.dat',cs2,1)
      call output(trim(directory)//'/gg.dat'  ,gg  ,3)
      call output(trim(directory)//'/B_ext_pot.dat'  ,B_ext_pot  ,3)
!
    endsubroutine wglobal
!***********************************************************************
    subroutine rglobal()
!
!  read global variables
!  10-jan-02/wolf: coded
!
      use Cdata
      use IO
!
      call input(trim(directory)//'/B_ext_pot.dat'  ,B_ext_pot  ,3,0)
      call input(trim(directory)//'/gg.dat'  ,gg  ,3,0)
      call input(trim(directory)//'/cs2.dat',cs2,1,0)
!
    endsubroutine rglobal
!***********************************************************************

endmodule Global
