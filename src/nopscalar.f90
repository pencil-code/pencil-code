! $Id: nopscalar.f90,v 1.9 2003-10-24 12:09:15 dobler Exp $

!  This modules solves the passive scalar advection equation

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Pscalar

  use Cparam
  use Cdata

  implicit none


  integer :: dummy           ! We cannot define empty namelists
  namelist /pscalar_init_pars/ dummy
  namelist /pscalar_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: i_rhoccm=0,i_ccmax=0,i_lnccm=0,i_lnccmz=0

  contains

!***********************************************************************
    subroutine register_pscalar()
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
           "$Id: nopscalar.f90,v 1.9 2003-10-24 12:09:15 dobler Exp $")
!
    endsubroutine register_pscalar
!***********************************************************************
    subroutine initialize_pscalar(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  24-nov-02/tony: coded
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
!  set to zero and then call the same initial condition
!  that was used in start.csh
!
      if(ip==0) print*,'f=',f
    endsubroutine initialize_pscalar
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
      real, dimension (mx,my,mz,mvar+maux) :: f
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
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu,glnrho
!
      intent(in)  :: f,df,uu,glnrho
!
      if(ip==0) print*,f,df,uu,glnrho
    endsubroutine dlncc_dt
!***********************************************************************
    subroutine rprint_pscalar(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!
!   6-jul-02/axel: coded
!
      use Sub
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=.true.
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_rhoccm=0; i_ccmax=0; i_lnccm=0; i_lnccmz=0
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_rhoccm=',i_rhoccm
        write(3,*) 'i_ccmax=',i_ccmax
        write(3,*) 'i_lnccm=',i_lnccm
        write(3,*) 'i_lnccmz=',i_lnccmz
        write(3,*) 'ilncc=',ilncc
      endif
!
    endsubroutine rprint_pscalar
!***********************************************************************
    subroutine calc_mpscalar
    endsubroutine calc_mpscalar
!***********************************************************************

endmodule Pscalar
