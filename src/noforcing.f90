module Forcing

!! Dummy module for hydro/mhd without forcing

  use Cdata

  implicit none

  integer :: dummy              ! We cannot define empty namelists
  namelist /forcing_init_pars/ dummy
  namelist /forcing_run_pars/  dummy

  contains

!***********************************************************************
    subroutine register_forcing()
!
!  add forcing in timestep()
!  11-may-2002/wolf: coded
!
      use Cdata
      use Mpicomm
      use Sub
!
      logical, save :: first=.true.
!
      if (.not. first) call stop_it('register_forcing called twice')
      first = .false.
!
      lforcing = .false.
!
!  identify version number
!
      if (lroot) call cvs_id( &
           "$RCSfile: noforcing.f90,v $", &
           "$Revision: 1.2 $", &
           "$Date: 2002-05-11 12:18:48 $")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine addforce(df)
!
      use Cdata
!
!  add forcing in timestep()
!
      real, dimension (mx,my,mz,mvar) :: df
!
    endsubroutine addforce
!***********************************************************************
    subroutine forcing1
!
!  forcing function, using a set of precomputed wavevectors
!
    endsubroutine forcing1
!***********************************************************************
    subroutine forcing2
!
!  helical forcing function, using a set of precomputed wavevectors
!
    endsubroutine forcing2
!***********************************************************************

endmodule Forcing
