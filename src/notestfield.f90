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

  include 'testfield.h'

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
    subroutine pencil_criteria_testfield()
!
!   All pencils that the Testfield module depends on are specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
!
    endsubroutine pencil_criteria_testfield
!***********************************************************************
    subroutine pencil_interdep_testfield(lpencil_in)
!
!  Interdependency among pencils from the Testfield module is specified here.
!
!  26-jun-05/anders: adapted from magnetic
!
      use Cdata
!
      logical, dimension(npencils) :: lpencil_in
!
    endsubroutine pencil_interdep_testfield
!***********************************************************************
    subroutine read_testfield_init_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_testfield_init_pars
!***********************************************************************
    subroutine write_testfield_init_pars(unit)
      integer, intent(in) :: unit
!
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine write_testfield_init_pars
!***********************************************************************
    subroutine read_testfield_run_pars(unit,iostat)
      integer, intent(in) :: unit
      integer, intent(inout), optional :: iostat
!
      if (present(iostat).and.NO_WARN) print*,iostat
      if (NO_WARN) print*,unit !(keep compiler quiet)
!
    endsubroutine read_testfield_run_pars
!***********************************************************************
    subroutine write_testfield_run_pars(unit)
      integer, intent(in) :: unit

      if (NO_WARN) print*,unit !(keep compiler quiet)

    endsubroutine write_testfield_run_pars
!***********************************************************************
    subroutine daatest_dt(f,df,p)
!
!  Dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
      type (pencil_case) :: p
!
      intent(in)     :: f, p
      intent(inout)  :: df     
!
      if(ip==0) print*, f, df, p  !(to keep compiler quiet)
!        
    endsubroutine daatest_dt
!***********************************************************************
    subroutine calc_ltestfield_pars(f)
!
!  calculate <uxb>, which is needed when lsoca=.false.
!
!  29-jan-06/axel: dummy routine
!
      use Cdata
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      intent(in)     :: f
!
    endsubroutine calc_ltestfield_pars
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
