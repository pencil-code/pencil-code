! $Id: nodustvelocity.f90,v 1.8 2003-10-24 13:17:31 dobler Exp $


!  This module takes care of everything related to velocity

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Dustvelocity

  use Cparam

  implicit none

  integer :: dummy              ! We cannot define empty namelists
  namelist /dustvelocity_init_pars/ dummy
  namelist /dustvelocity_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ud2m=0,i_udm2=0,i_oudm=0,i_od2m=0
  integer :: i_udrms=0,i_udmax=0,i_odrms=0,i_odmax=0
  integer :: i_udxmz=0,i_udymz=0,i_udzmz=0,i_udmx=0,i_udmy=0,i_udmz=0
  integer :: i_udxmxy=0,i_udymxy=0,i_udzmxy=0
  integer :: i_divud2m=0,i_epsKd=0

  contains

!***********************************************************************
    subroutine register_dustvelocity()
!
!  Initialise variables which should know that we solve the hydro
!  equations: iuu, etc; increase nvar accordingly.
!
!  18-mar-03/axel: dummy routine
!
      use Cdata
      use Mpicomm, only: lroot,stop_it
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_dustvelocity called twice')
      first = .false.
!
      ldustvelocity = .false.
!
!  identify version number (generated automatically by CVS)
!
      if (lroot) call cvs_id( &
           "$Id: nodustvelocity.f90,v 1.8 2003-10-24 13:17:31 dobler Exp $")
!
    endsubroutine register_dustvelocity
!***********************************************************************
    subroutine initialize_dustvelocity()
!
!  Perform any post-parameter-read initialization i.e. calculate derived
!  parameters.
!
!  18-mar-03/axel: dummy routine
!
    endsubroutine initialize_dustvelocity
!***********************************************************************
    subroutine init_uud(f,xx,yy,zz)
!
!  initialise uu and lnrho; called from start.f90
!  Should be located in the Dustvelocity module, if there was one.
!
!  18-mar-03/axel: dummy routine
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz) :: xx,yy,zz
!
      if (ip==0) print*,f,xx,yy,zz  !(keep compiler quiet)
    endsubroutine init_uud
!***********************************************************************
    subroutine duud_dt(f,df,uu,uud,divud,ud2,udij)
!
!  velocity evolution
!  calculate dud/dt = - ud.gradud - 2Omega x ud + grav + Fvisc
!  no pressure gradient force for dust!
!
!  18-mar-03/axel: dummy routine
!
      use Cdata
      use Sub
      use IO
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3,3) :: udij
      real, dimension (nx,3) :: uu,uud
      real, dimension (nx) :: ud2,divud
!
      if(ip==0) print*,f,df,uu,uud,divud,ud2,udij  !(keep compiler quiet)
    endsubroutine duud_dt
!***********************************************************************
    subroutine rprint_dustvelocity(lreset,lwrite)
!
!  reads and registers print parameters relevant for hydro part
!
!   3-may-02/axel: coded
!  27-may-02/axel: added possibility to reset list
!
      use Cdata
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
        write(3,*) 'i_ud2m=',i_ud2m
        write(3,*) 'i_udm2=',i_udm2
        write(3,*) 'i_od2m=',i_od2m
        write(3,*) 'i_oudm=',i_oudm
        write(3,*) 'i_udrms=',i_udrms
        write(3,*) 'i_udmax=',i_udmax
        write(3,*) 'i_odrms=',i_odrms
        write(3,*) 'i_odmax=',i_odmax
        write(3,*) 'i_udmx=',i_udmx
        write(3,*) 'i_udmy=',i_udmy
        write(3,*) 'i_udmz=',i_udmz
        write(3,*) 'i_divud2m=',i_divud2m
        write(3,*) 'i_epsKd=',i_epsKd
        write(3,*) 'nname=',nname
        write(3,*) 'iuud=',iuud
        write(3,*) 'iudx=',iudx
        write(3,*) 'iudy=',iudy
        write(3,*) 'iudz=',iudz
        write(3,*) 'i_udxmz=',i_udxmz
        write(3,*) 'i_udymz=',i_udymz
        write(3,*) 'i_udzmz=',i_udzmz
        write(3,*) 'i_udxmxy=',i_udxmxy
        write(3,*) 'i_udymxy=',i_udymxy
        write(3,*) 'i_udzmxy=',i_udzmxy
      endif
!
      if(ip==0) print*,lreset  !(to keep compiler quiet)
    endsubroutine rprint_dustvelocity
!***********************************************************************

endmodule Dustvelocity
