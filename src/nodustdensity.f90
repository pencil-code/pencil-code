! $Id: nodustdensity.f90,v 1.3 2003-03-19 06:58:18 brandenb Exp $

!  This module is used both for the initial condition and during run time.
!  It contains dlnrhod_dt and init_lnrhod, among other auxiliary routines.

module Dustdensity

  use Cparam
  use Cdata

  implicit none

  integer :: dummy              ! We cannot define empty namelists
  namelist /dustdensity_init_pars/ dummy
  namelist /dustdensity_run_pars/  dummy

  ! diagnostic variables (needs to be consistent with reset list below)
  integer :: i_rhodm=0

  contains

!***********************************************************************
    subroutine register_dustdensity()
!
!  Initialise variables which should know that we solve the
!  compressible hydro equations: ilnrhod; increase nvar accordingly.
!
!  18-mar-03/axel: adapted from dustdensity
!
      use Mpicomm, only: stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_dustdensity called twice')
      first = .false.
!
      ldustdensity = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: nodustdensity.f90,v 1.3 2003-03-19 06:58:18 brandenb Exp $")
!
    endsubroutine register_dustdensity
!***********************************************************************
    subroutine initialize_dustdensity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-mar-03/axel: adapted from dustdensity
!
!  do nothing
!
    endsubroutine initialize_dustdensity
!***********************************************************************
    subroutine init_lnrhod(f,xx,yy,zz)
!
!  initialise lnrhod; called from start.f90
!
!  18-mar-03/axel: adapted from dustdensity
!
      use Mpicomm
      use Sub
      use IO
      use Global
      use Gravity
      use Initcond
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz ! keep compiler quiet
    endsubroutine init_lnrhod
!***********************************************************************
    subroutine dlnrhod_dt(f,df,uud,glnrhod,divud,lnrhod)
!
!  continuity equation
!  calculate dlnrhod/dt = - u.gradlnrhod - divud
!
!  18-mar-03/axel: adapted from dustdensity
!
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uud,glnrhod
      real, dimension (nx) :: lnrhod,divud
!
      if(ip==0) print*,f,df,uud,glnrhod,divud,lnrhod ! keep compiler quiet
    endsubroutine dlnrhod_dt
!***********************************************************************
    subroutine rprint_dustdensity(lreset)
!
!  reads and registers print parameters relevant for compressible part
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Sub
!
      logical :: lreset
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_rhodm=',i_rhodm
      write(3,*) 'nname=',nname
      write(3,*) 'ilnrhod=',ilnrhod
!
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_dustdensity
!***********************************************************************

endmodule Dustdensity
