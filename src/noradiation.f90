! $Id: noradiation.f90,v 1.2 2002-07-16 11:24:11 nilshau Exp $


module Radiation

  use Cparam

  implicit none

  integer :: dummyuu           ! We cannot define empty namelists
  namelist /radiation_init_pars/ dummyuu
  namelist /radiation_run_pars/  dummyuu

  ! other variables (needs to be consistent with reset list below)
  integer :: i_frms=0,i_fmax=0,i_Erms=0,i_Emax=0

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
           "$Id: noradiation.f90,v 1.2 2002-07-16 11:24:11 nilshau Exp $")
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
   subroutine de_dt(f,df,rho1,divu,uu,uij,TT1,gamma)
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
      real :: gamma
!
      if(ip==0) print*,f,df,uu,rho1,TT1 !(keep compiler quiet)
    endsubroutine de_dt
!*******************************************************************
    subroutine rprint_radiation(lreset)
!
!  reads and registers print parameters relevant for radiative part
!
!  16-jul-02/nils: adapted from rprint_hydro
!
      use Cdata
      use Sub
!
      integer :: iname
      logical :: lreset
!
!  write column where which radiative variable is stored
!
      write(3,*) 'i_frms=',i_frms
      write(3,*) 'i_fmax=',i_fmax
      write(3,*) 'i_Erms=',i_Erms
      write(3,*) 'i_Emax=',i_Emax
      write(3,*) 'nname=',nname
      write(3,*) 'ie=',ie
      write(3,*) 'ifx=',ifx
      write(3,*) 'ify=',ify
      write(3,*) 'ifz=',ifz
!
    endsubroutine rprint_radiation
!***********************************************************************
  end module Radiation
