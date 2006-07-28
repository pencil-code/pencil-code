 ! $Id: global_avgs.f90,v 1.6 2006-07-28 13:19:38 wlyra Exp $

module Global

!
  use Cparam
  use Mpicomm

  implicit none

  include 'global.h'

  interface set_global
    module procedure set_global_mvar
    module procedure set_global_vect
    module procedure set_global_scal
    module procedure set_global_coarse_vect
    module procedure set_global_coarse_scal
  endinterface

  interface set_global_point
    module procedure set_global_scal_point
    module procedure set_global_vect_point
  endinterface

  interface get_global
    module procedure get_global_mvar 
    module procedure get_global_vect
    module procedure get_global_scal
    module procedure get_global_coarse_vect
    module procedure get_global_coarse_scal
  endinterface

  interface get_global_point
     module procedure get_global_vect_point
     module procedure get_global_scal_point
  endinterface
!
!
  real, dimension (mx,my,mz,mvar) :: fborder
  real, dimension (mx,my,mz,3) :: gg
  real, dimension (mx,my,mz) :: rho,cs2,rhos
  real, dimension (mx,my,mz,3) :: bbs
  real, dimension (mx,my,mz,3) :: uus
!
  real, dimension (10,3) :: bavg_coarse,uavg_coarse
  real, dimension (10) :: rhoavg_coarse
!  
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
    subroutine set_global_mvar(var,m,n,j,label,length)
!
!  set (m,n)-pencil of the global vector variable identified by LABEL
!
!  27-jul-06/wlad coded
!
      use Mpicomm, only: stop_it
!
      integer :: length
      real, dimension(length) :: var
      integer :: m,n,j
      character (len=*) ::label
!
      if (length/=nx) then
        print*, 'set_global_mvar: only nx pencils allowed!'
        call stop_it('set_global_mvar')
      endif
!
      select case(label)
!
      case ('fborder')
         fborder(l1:l2,m,n,j) = var
!
      case default
         if (lroot) print*, 'set_global_vect: No such value for label', trim(label)
         call stop_it('set_global_mvar')
!
      endselect
!
    endsubroutine set_global_mvar
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
!
      case ('gg')
         gg(l1:l2,m,n,1:3) = var
!
      case ('bbs')
         bbs(l1:l2,m,n,1:3) = var
!
      case ('uus')
         uus(l1:l2,m,n,1:3) = var
!
      case default
         if (lroot) print*, 'set_global_vect: No such value for label', trim(label)
         call stop_it('set_global_vect')
!
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
      case ('rho')
         if (length==nx) then
            rho(l1:l2,m,n) = var
         elseif (length==mx) then
            rho(:,m,n) = var
         endif

      case ('cs2')
         if (length==nx) then
            cs2(l1:l2,m,n) = var
         elseif (length==mx) then
            cs2(:,m,n) = var
         endif

      case ('rhos')
         if (length==nx) then
            rhos(l1:l2,m,n) = var
         elseif (length==mx) then
            rhos(:,m,n) = var
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
    subroutine set_global_coarse_vect(var,label,length)
!
!  set global variable identified by LABEL
!
!  13-jun-05/anders: adapted
!
      integer :: length
      real, dimension(length,3) :: var
      character (len=*) ::label
!
      select case(label)
!
      case ('uavg')
         uavg_coarse(1:length,1:3) = var
      case ('bavg')
         bavg_coarse(1:length,1:3) = var
!
      case default
         if (lroot) &
              print*, 'set_global_coarse_vect: No such value for label', trim(label)
         call stop_it('set_global_coarse_vect')
!
      endselect
!
    endsubroutine set_global_coarse_vect
!***********************************************************************
    subroutine set_global_coarse_scal(var,label,length)
!
!  set global variable identified by LABEL
!
!  13-jun-05/anders: adapted
!
      integer :: length
      real, dimension(length) :: var
      character (len=*) ::label
!
      select case(label)
!
      case ('rhoavg')
         rhoavg_coarse(1:length) = var
!
      case default
         if (lroot) &
              print*, 'set_global_coarse_scal: No such value for label', trim(label)
         call stop_it('set_global_coarse_scal')
!
      endselect
!
    endsubroutine set_global_coarse_scal
!**********************************************************************
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
    subroutine get_global_mvar(var,m,n,j,label)
!
!  Get (m,n)-pencil of the global vector variable identified by LABEL.
!
!  27-jul-06/wolf coded
!
      real, dimension(nx) :: var
      integer :: m,n,j
      character (len=*) ::label
!
      select case(label)
!
      case ('fborder')
        var = fborder(l1:l2,m,n,j)
!
      case default
        if (lroot) print*, 'get_global_mvar: No such value for label', trim(label)
        call stop_it('get_global_mvar')
!
      endselect
!
    endsubroutine get_global_mvar
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
      case ('bbs')
        var = bbs(l1:l2,m,n,1:3)
!
      case ('uus')
        var = uus(l1:l2,m,n,1:3)
!
      case default
        if (lroot) print*, 'get_global_vect: No such value for label', trim(label)
        call stop_it('get_global_vect')
!
      endselect
!
    endsubroutine get_global_vect
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
      case ('rho')
         var = rho(l1:l2,m,n)
!
      case ('cs2')
         var = cs2(l1:l2,m,n)
!
      case ('rhos')
         var = rhos(l1:l2,m,n)
!
      case default
         if (lroot) print*, 'get_global_scal: No such value for label', trim(label)
         call stop_it('get_global_scal')
!
      endselect
!
    endsubroutine get_global_scal
!***********************************************************************
    subroutine get_global_coarse_vect(var,label,length)
!                                                                               
!  Get (m,n)-pencil of the global vector variable identified by LABEL.          
!                                                                               
!  18-jul-02/wolf coded                                                         
!                             
      integer :: length
      real, dimension(length,3) :: var
      character (len=*) ::label
!
      select case(label)
!
      case ('uavg')
        var = uavg_coarse
!
      case ('bavg')
        var = bavg_coarse
!
      case default
        if (lroot) print*, 'get_global_coarse_vect: No such value for label', trim(label)
        call stop_it('get_global_coarse_vect')
!
      endselect
!
    endsubroutine get_global_coarse_vect
!********************************************************************
    subroutine get_global_coarse_scal(var,label,length)
!
!  Get (m,n)-pencil of the global vector variable identified by LABEL.
!
!  18-jul-02/wolf: coded
!
      integer :: length
      real, dimension(length) :: var
      character (len=*) ::label
!
      select case(label)
!
      case ('rhoavg')
        var = rhoavg_coarse
!
      case default
        if (lroot) print*, 'get_global_coarse_scal: No such value for label', trim(label)
        call stop_it('get_global_coarse_scal')
!
      endselect
!
    endsubroutine get_global_coarse_scal
!*************************************************************************
    subroutine set_global_scal_point(var,l,m,n,label)
!
!  set point value of the global scalar variable identified by LABEL
!
!  20-jun-05/anders: adapted
!
      real :: var
      integer :: l,m,n
      character (len=*) ::label
!
      select case(label)
!
      case ('cs2')
         cs2(l,m,n) = cs2(l,m,n) + var
!
      case ('rho')
         rho(l,m,n) = rho(l,m,n) + var
!
      case ('rhos')
         rhos(l,m,n) = rhos(l,m,n) + var
!
      case default
         if (lroot) print*, &
              'set_global_scal_point: No such value for label=', trim(label)
         call stop_it('set_global_scal_point')
!
      endselect
!
    endsubroutine set_global_scal_point
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
!
      case ('cs2')
         var = cs2(l,m,n)
!
      case ('rho')
         var = rho(l,m,n)
!
      case ('rhos')
         var = rhos(l,m,n)
!
      case default
        if (lroot) print*, &
            'get_global_scal_point: No such value for label', trim(label)
        call stop_it('get_global_scal_point')
!
      endselect
!
    endsubroutine get_global_scal_point
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
      case ('uus')
        var = uus(l,m,n,1:3)
!
      case ('bbs')
        var = bbs(l,m,n,1:3)
!
      case default
        if (lroot) print*, &
            'get_global_vect_point: No such value for label', trim(label)
        call stop_it('get_global_vect_point')
!
      endselect
!
    endsubroutine get_global_vect_point
!*********************************************************************
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
      case ('rho')
        if (present(der6)) call del6_other(rho,der6)
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
      call output(trim(directory)//'/rho.dat',rho,1)
      call output(trim(directory)//'/cs2.dat',cs2,1)
      call output(trim(directory)//'/rhos.dat',rhos,1)
      call output(trim(directory)//'/gg.dat',gg,3)
      call output(trim(directory)//'/bbs.dat',bbs,3)
      call output(trim(directory)//'/uus.dat',uus,3)
      call output(trim(directory)//'/fborder.dat',fborder,mvar)
      call output(trim(directory)//'/uavg_coarse.dat',uavg_coarse,3,10)
      call output(trim(directory)//'/bavg_coarse.dat',bavg_coarse,3,10)
      call output(trim(directory)//'/rhoavg_coarse.dat',rhoavg_coarse,1,10)
!
    endsubroutine wglobal
!***********************************************************************
    subroutine rglobal()
!
!  read global variables
!
!  10-jan-02/wolf: coded
!
      use Cdata, only: directory
      use IO, only: input,input_coarse
!
      call input(trim(directory)//'/rho.dat',rho,1,0)
      call input(trim(directory)//'/cs2.dat',cs2,1,0)
      call input(trim(directory)//'/rhos.dat',rhos,1,0)
      call input(trim(directory)//'/gg.dat',gg,3,0)
      call input(trim(directory)//'/bbs.dat',bbs,3,0)
      call input(trim(directory)//'/uus.dat',uus,3,0)
      call input(trim(directory)//'/fborder.dat',fborder,mvar,0)
      call input_coarse(trim(directory)//'/uavg_coarse.dat',uavg_coarse,3,0,10)
      call input_coarse(trim(directory)//'/bavg_coarse.dat',bavg_coarse,3,0,10)
      call input_coarse(trim(directory)//'/rhoavg_coarse.dat',rhoavg_coarse,1,0,10)
!
    endsubroutine rglobal
!***********************************************************************

endmodule Global
