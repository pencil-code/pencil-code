! $Id: noradiation.f90,v 1.1 2002-07-16 08:09:33 nilshau Exp $


module Radiation

  use Cparam

  implicit none

  integer :: dummyuu           ! We cannot define empty namelists
  namelist /radiation_init_pars/ dummyuu
  namelist /radiation_run_pars/  dummyuu

  contains

!***********************************************************************
    subroutine register_rad()
!
!  15-jul-2002/nils: dummy routine
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_rad called twice')
      first = .false.
!
      lradiation = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: noradiation.f90,v 1.1 2002-07-16 08:09:33 nilshau Exp $")
!
    endsubroutine register_rad
!***********************************************************************
    subroutine init_rad(f,xx,yy,zz)
!
!  initialise radiation; called from start.f90
!  15-jul-2002/nils: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz !(keep compiler quiet)
    endsubroutine init_rad
!***********************************************************************
   subroutine de_dt(f,df,rho1,divu,uu,uij,TT1)
!
!  15-jul-2002/nils: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu
      real, dimension (nx) :: rho1,TT1
      real, dimension (nx,3,3) :: uij
      real, dimension (nx) :: divu
!
      if(ip==0) print*,f,df,uu,rho1,TT1 !(keep compiler quiet)
    endsubroutine de_dt
!***********************************************************************
  end module Radiation
