! $Id: nograv.f90,v 1.4 2002-06-01 02:56:21 brandenb Exp $

module Gravity

!
!  Dummy model: no gravity
!

  use Cparam

  implicit none

  integer :: dummy              ! We cannot define empty namelists
  namelist /grav_z_init_pars/ dummy
  namelist /grav_z_run_pars/  dummy

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
           "$RCSfile: nograv.f90,v $", &
           "$Revision: 1.4 $", &
           "$Date: 2002-06-01 02:56:21 $")
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
    endsubroutine duu_dt_grav
!***********************************************************************
    subroutine potential(xmn,ymn,zmn,rmn, pot)
!
!  gravity potential
!  28-mar-02/axel: adapted from grav_z
!
      use Cdata, only: nx,ny,nz,gravz,lroot
      use Sub, only: poly
!
      real, dimension (nx,1,1) :: xmn,ymn,zmn,rmn, pot
!
      if (lroot) print*,'potential: should not have been called'
      pot = 0.
!
    endsubroutine potential
!***********************************************************************

endmodule Gravity
