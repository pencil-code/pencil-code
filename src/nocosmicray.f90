! $Id: nocosmicray.f90,v 1.4 2003-10-10 17:07:59 snod Exp $

!  This modules solves the passive scalar advection equation

module CosmicRay

  use Cparam
  use Cdata

  implicit none


  integer :: dummy           ! We cannot define empty namelists
  namelist /cosmicray_init_pars/ dummy
  namelist /cosmicray_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
!  integer :: i_rhocrm=0,i_crmax=0,i_lncrm=0,i_lncrmz=0

  contains

!***********************************************************************
    subroutine register_cosmicray()
!
!  Initialise variables which should know that we solve for active 
!  scalar: iecr - the cosmic ray energy density; increase nvar accordingly
!
!  09-oct-03/tony: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_cosmicray called twice')
      first = .false.
!
      lcosmicray = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: nocosmicray.f90,v 1.4 2003-10-10 17:07:59 snod Exp $")
!
    endsubroutine register_cosmicray
!***********************************************************************
    subroutine initialize_cosmicray(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  09-oct-03/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if(ip==0) print*,'f=',f
    endsubroutine initialize_cosmicray
!***********************************************************************
    subroutine init_ecr(f,xx,yy,zz)
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_ecr
!***********************************************************************
    subroutine decr_dt(f,df,uu,rho1,divu,bij,bb)
!
!  cosmic ray density evolution
!
!  DON'T calculate dcr/dt=????
!
!   09-oct-03/tony: coded
!
      use Sub
!
      real, intent(in), dimension (mx,my,mz,mvar+maux) :: f
      real, intent(inout), dimension (mx,my,mz,mvar) :: df
      real, intent(in), dimension (nx,3,3) :: bij
      real, intent(in), dimension (nx,3) :: uu,bb
      real, intent(in), dimension (nx) :: divu,rho1
!
      if(ip==0) print*,f,df,uu,rho1,divu,bij,bb
    endsubroutine decr_dt
!***********************************************************************
    subroutine rprint_cosmicray(lreset)
!
!  reads and registers print parameters relevant for cosmic rays
!
!   09-oct-03/tony: coded
!
      use Sub
!
      logical :: lreset
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
!        i_rhoccm=0; i_ccmax=0; i_lnccm=0; i_lnccmz=0
      endif
!
!  write column where which magnetic variable is stored
!
!      write(3,*) 'i_lnccmz=',i_lnccmz
      write(3,*) 'iecr=',iecr
!
    endsubroutine rprint_cosmicray
!***********************************************************************

endmodule CosmicRay
