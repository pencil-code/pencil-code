! $Id: noshear.f90,v 1.1 2002-07-04 10:10:55 nilshau Exp $

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
           "$Id: noshear.f90,v 1.1 2002-07-04 10:10:55 nilshau Exp $")
!
    endsubroutine register_shear
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
      real, dimension (mx,my,mz,mvar) :: f,df
!
      if(ip==0) print*,f,df
!
    end subroutine shearing
!***********************************************************************
  end module Shear
