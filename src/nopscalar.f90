! $Id: nopscalar.f90,v 1.1 2002-07-06 21:00:14 brandenb Exp $

!  This modules solves the passive scalar advection equation

module Pscalar

  use Cparam
  use Cdata

  implicit none

  integer :: ilncc=0

  integer :: dummy           ! We cannot define empty namelists
  namelist /pscalar_init_pars/ dummy
  namelist /pscalar_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: i_ccm=0,i_ccmax=0

  contains

!***********************************************************************
    subroutine register_lncc()
!
!  Initialise variables which should know that we solve for passive
!  scalar: ilncc; increase nvar accordingly
!
!  6-jul-02/axel: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_lncc called twice')
      first = .false.
!
      lpscalar = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: nopscalar.f90,v 1.1 2002-07-06 21:00:14 brandenb Exp $")
!
    endsubroutine register_lncc
!***********************************************************************
    subroutine init_lncc(f,xx,yy,zz)
!
!  initialise magnetic field; called from start.f90
!  AB: maybe we should here call different routines (such as rings)
!  AB: and others, instead of accummulating all this in a huge routine.
!  We have an init parameter (initlncc) to stear magnetic i.c. independently.
!
!   6-jul-02/axel: coded
!
      use Cdata
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_lncc
!***********************************************************************
    subroutine dlncc_dt(f,df,uu,glnrho)
!
!  magnetic field evolution
!
!  calculate dc/dt=-uu.gradlncc
!
!   6-jul-02/axel: coded
!
      use Sub
!
      real, dimension (mx,my,mz,mvar) :: f,df
      real, dimension (nx,3) :: uu,glnrho
!
      intent(in)  :: f,df,uu,glnrho
!
      if(ip==0) print*,f,df,uu,glnrho
    endsubroutine dlncc_dt
!***********************************************************************
    subroutine rprint_pscalar(lreset)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   6-jul-02/axel: coded
!
      use Sub
!
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_ccm=0; i_ccmax=0
      endif
!
!  write column where which magnetic variable is stored
!
      write(3,*) 'i_ccm=',i_ccm
      write(3,*) 'i_ccmax=',i_ccmax
      write(3,*) 'ilncc=',ilncc
!
    endsubroutine rprint_pscalar
!***********************************************************************

endmodule Pscalar
