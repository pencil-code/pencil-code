! $Id: nodustdensity.f90,v 1.20 2004-05-20 03:58:44 brandenb Exp $

!  This module is used both for the initial condition and during run time.
!  It contains dlnrhod_dt and init_lnrhod, among other auxiliary routines.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Dustdensity

  use Cparam
  use Cdata

  implicit none

  integer :: dummy              ! We cannot define empty namelists
  logical :: ldustnulling=.false.

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
           "$Id: nodustdensity.f90,v 1.20 2004-05-20 03:58:44 brandenb Exp $")
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
    subroutine init_nd(f,xx,yy,zz)
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz ! keep compiler quiet
    endsubroutine init_nd
!***********************************************************************
    subroutine dndmd_dt(f,df,rho1,TT1,cs2,uud,divud,cc,cc1,gnd)
!
!  continuity equation
!  calculate dlnrhod/dt = - u.gradlnrhod - divud
!
!  18-mar-03/axel: adapted from dustdensity
!
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu,uud,gnd
      real, dimension (nx) :: rho1,TT1,cs2,divud,cc,cc1
!
      if(ip==0) print*,f,df,rho1,uu,uud,divud,gnd ! keep compiler quiet
    endsubroutine dndmd_dt
!***********************************************************************
    subroutine redist_mdbins(f)
!
!  Redistribute dust number density and dust density in mass bins
!
!  4-may-2004/wolf: Adapted from dustdensity.f90
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if(ip==0) print*,f(1,1,1,1)
    endsubroutine redist_mdbins
!***********************************************************************
    subroutine null_dust_vars(f)
!
!  Dummy routine
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if(ip==0) print*,f(1,1,1,1)
    endsubroutine null_dust_vars
!***********************************************************************
    subroutine rprint_dustdensity(lreset,lwrite)
!
!  reads and registers print parameters relevant for compressible part
!
!   7-jun-02/axel: incoporated from subroutine pde
!
      use Sub
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_rhodm=',i_rhodm
        write(3,*) 'nname=',nname
        write(3,*) 'ind=',ind
      endif
!
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_dustdensity
!***********************************************************************

endmodule Dustdensity
