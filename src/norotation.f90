! $Id: norotation.f90,v 1.1 2002-07-02 17:14:41 nilshau Exp $

!  This modules deals with all aspects of rotation; if no
!  rotation are invoked, a corresponding replacement dummy
!  routine is used instead which absorbs all the calls to the
!  magnetically relevant subroutines listed in here.

module Rotation

  use Sub
  use Cdata

  implicit none

  integer :: dummy              ! We cannot define empty namelists
  namelist /rotation_init_pars/ &
       dummy

  contains

!***********************************************************************
    subroutine register_rot()
!
!  Initialise variables
!
!  2-july-02/nils: coded
!
      use Mpicomm
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_rot called twice')
      first = .false.
!
      lrotation = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$Id: norotation.f90,v 1.1 2002-07-02 17:14:41 nilshau Exp $")
!
    endsubroutine register_rot
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
  end module Rotation
