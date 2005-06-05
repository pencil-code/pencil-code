! $Id: notestfield.f90,v 1.1 2005-06-05 12:44:27 brandenb Exp $

!  This modules deals with all aspects of testfield fields; if no
!  testfield fields are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  testfield relevant subroutines listed in here.

!** AUTOMATIC CPARAM.INC GENERATION ****************************
! Declare (for generation of cparam.inc) the number of f array
! variables and auxiliary variables added by this module
!
! MVAR CONTRIBUTION 0
! MAUX CONTRIBUTION 0
!
!***************************************************************

module Testfield

  use Cparam

  implicit none

  real :: dummy=0.
  namelist /testfield_init_pars/ &
       dummy
  namelist /testfield_run_pars/ &
       dummy

  contains

!***********************************************************************
    subroutine register_testfield()
!
!  Dummy routine
!
    endsubroutine register_testfield
!***********************************************************************
    subroutine initialize_testfield(f)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
!
      if(ip==0) print*,f  !(to keep compiler quiet)
    endsubroutine initialize_testfield
!***********************************************************************
    subroutine init_aatest(f,xx,yy,zz)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz)      :: xx,yy,zz
!
      if(ip==0) print*,f,xx,yy,zz  !(to keep compiler quiet)
    endsubroutine init_aatest
!***********************************************************************
    subroutine daatest_dt(f,df,uu)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      real, dimension (nx,3) :: uu
!
      intent(in)     :: f,uu
      intent(inout)  :: df     
!
      if(ip==0) print*,f,df,uu  !(to keep compiler quiet)
    endsubroutine daatest_dt
!***********************************************************************
    subroutine rprint_testfield(lreset,lwrite)
!
!  Dummy routine
!
      use Cdata
!
      logical :: lreset
      logical, optional :: lwrite
!
      if(ip==0) print*,lreset,lwrite  !(to keep compiler quiet)
    endsubroutine rprint_testfield
!***********************************************************************

endmodule Testfield
