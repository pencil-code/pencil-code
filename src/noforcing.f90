! $Id: noforcing.f90,v 1.5 2002-06-25 14:58:47 dobler Exp $

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
           "$Id: noforcing.f90,v 1.5 2002-06-25 14:58:47 dobler Exp $")
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
      if(ip==1) print*,df !(to remove compiler warnings)
    endsubroutine addforce
!***********************************************************************

endmodule Forcing
