! $Id: noshear.f90,v 1.4 2003-06-16 04:41:11 brandenb Exp $

!  This modules deals with all aspects of shear; if no
!  shear are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  shear relevant subroutines listed in here.

module Shear

  use Sub
  use Cdata

  implicit none

  integer :: dummy              ! We cannot define empty namelists
  namelist /shear_init_pars/ &
       dummy

  namelist /shear_run_pars/ &
       dummy

  contains

!***********************************************************************
    subroutine register_shear()
!
!  Initialise variables
!
!  2-july-02/nils: coded
!
      use Mpicomm
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_shear called twice')
      first = .false.
!
      lshear = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: noshear.f90,v 1.4 2003-06-16 04:41:11 brandenb Exp $")
!
    endsubroutine register_shear
!***********************************************************************
    subroutine initialize_shear()
!
!  21-nov-02/tony: coded

! dummy

  endsubroutine initialize_shear
!***********************************************************************
    subroutine shearing(f,df)
!
!  Calculates the actuall shear terms
!
!  2-july-02/nils: coded
!
      use Cparam
      use Deriv
!
      real, dimension (mx,my,mz,mvar+maux) :: f
      real, dimension (mx,my,mz,mvar) :: df
!
      if(ip==0) print*,f,df !(to keep compiler quiet)
    end subroutine shearing
!***********************************************************************
    subroutine advance_shear
!
!  Dummy routine: deltay remains unchanged
!
! 18-aug-02/axel: incorporated from nompicomm.f90
!
      use Cdata
!
!  print identifier
!
      if (headtt.or.ldebug) print*,'advance_shear: deltay=const=',deltay
!
    end subroutine advance_shear
!***********************************************************************
  end module Shear
