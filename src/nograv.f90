! $Id: nograv.f90,v 1.12 2002-07-08 19:28:14 brandenb Exp $

module Gravity

!
!  Dummy model: no gravity
!

  use Cparam

  implicit none

  real :: z1,z2,zref,gravz,zinfty  !(used by Entropy and Density)
  character (len=30) :: grav_profile='const'  !(used by Density)

  integer :: dummy              ! We cannot define empty namelists
  namelist /grav_init_pars/ dummy
  namelist /grav_run_pars/  dummy

  contains

!***********************************************************************
    subroutine register_grav()
!
!  initialise gravity flags
!
!  9-jan-02/wolf: coded
! 28-mar-02/axel: adapted from grav_z
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_grav called twice')
      first = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: nograv.f90,v 1.12 2002-07-08 19:28:14 brandenb Exp $")
!
      lgrav = .false.
      lgravz = .false.
      lgravr = .false.
!
    endsubroutine register_grav
!***********************************************************************
    subroutine init_grav(f,xx,yy,zz)
!
!  initialise gravity; called from start.f90
!   9-jan-02/wolf: coded
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
! Not doing anything (this might change if we decide to store gg)
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
    endsubroutine init_grav
!***********************************************************************
    subroutine duu_dt_grav(f,df)
!
!  add nothing to duu/dt
!
! 28-mar-02/axel: adapted from grav_z
!
      use Cdata
!      use Mpicomm
      use Sub
      use Slices
!
      real, dimension (mx,my,mz,mvar) :: f,df
!
      if(ip==0) print*,f,df  !(keep compiler quiet)
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential(xmn,ymn,zmn,pot,grav,rmn)
!
!  gravity potential
!  28-mar-02/axel: adapted from grav_z
!
      use Cdata, only: nx,lroot
!
      real, dimension (nx) :: xmn,pot
      real :: ymn,zmn
      real, optional, dimension (nx) :: rmn
      real, optional, dimension (nx,3) :: grav
!
      if (lroot) print*,'potential: should not have been called'
      pot = 0.
!
      if(ip==0) print*,xmn,ymn,zmn,rmn,grav
    endsubroutine potential
!***********************************************************************

endmodule Gravity
