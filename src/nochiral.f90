! $Id: nochiral.f90,v 1.1 2004-05-29 06:33:31 brandenb Exp $

!  This modules solves two reactive scalar advection equations
!  This is used for modeling the spatial evolution of left and
!  right handed aminoacids.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Chiral

  use Cparam
  use Cdata

  implicit none

  integer :: dummy           ! We cannot define empty namelists
  namelist /chiral_init_pars/ dummy
  namelist /chiral_run_pars/  dummy

  ! other variables (needs to be consistent with reset list below)
  integer :: i_XX_chiralmax=0, i_YY_chiralmax=0

  contains

!***********************************************************************
    subroutine register_chiral()
!
!  Initialise variables which should know that we solve for passive
!  scalar: iXX_chiral and iYY_chiral; increase nvar accordingly
!
!  28-may-04/axel: adapted from pscalar
!
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_chiral called twice')
      first = .false.
!
      lchiral = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: nochiral.f90,v 1.1 2004-05-29 06:33:31 brandenb Exp $")
!
    endsubroutine register_chiral
!***********************************************************************
    subroutine initialize_chiral(f)
!
!  Perform any necessary post-parameter read initialization
!  Dummy routine
!
!  28-may-04/axel: adapted from pscalar
!
      real, dimension (mx,my,mz,mvar+maux) :: f
! 
!  set to zero and then call the same initial condition
!  that was used in start.csh
!   
      if(ip==0) print*,'f=',f
    endsubroutine initialize_chiral
!***********************************************************************
    subroutine init_chiral(f,xx,yy,zz)
!
!  initialise passive scalar field; called from start.f90
!
!  28-may-04/axel: adapted from pscalar
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz,prof
!
      if(ip==0) print*,xx,yy,zz !(prevent compiler warnings)
    endsubroutine init_chiral
!***********************************************************************
    subroutine dXY_chiral_dt(f,df,uu)
!
!  passive scalar evolution
!  calculate dc/dt=-uu.gcc + pscaler_diff*[del2cc + glnrho.gcc]
!
!  28-may-04/axel: adapted from pscalar
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu
!
      intent(in)  :: f,df,uu
!
      if(ip==0) print*,f,df,uu
    endsubroutine dXY_chiral_dt
!***********************************************************************
    subroutine rprint_chiral(lreset,lwrite)
!
!  reads and registers print parameters relevant for magnetic fields
!
!  28-may-04/axel: adapted from pscalar
!
      logical :: lreset,lwr
      logical, optional :: lwrite
!
      lwr = .false.
      if (present(lwrite)) lwr=lwrite
!
!  reset everything in case of reset
!  (this needs to be consistent with what is defined above!)
!
      if (lreset) then
        i_XX_chiralmax=0
        i_YY_chiralmax=0
      endif
!
!  write column where which magnetic variable is stored
!
      if (lwr) then
        write(3,*) 'i_XX_chiralmax=',i_XX_chiralmax
        write(3,*) 'i_YY_chiralmax=',i_YY_chiralmax
        write(3,*) 'iXX_chiral=',iXX_chiral
        write(3,*) 'iYY_chiral=',iYY_chiral
      endif
!
    endsubroutine rprint_chiral
!***********************************************************************

endmodule Chiral
