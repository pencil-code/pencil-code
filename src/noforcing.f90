! $Id: noforcing.f90,v 1.11 2003-04-10 06:58:24 brandenb Exp $

module Forcing

!! Dummy module for hydro/mhd without forcing

  use Cdata
  use General

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
           "$Id: noforcing.f90,v 1.11 2003-04-10 06:58:24 brandenb Exp $")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine initialize_forcing(lstart)
!
!  initialize random number generator in processor-dependent fashion
!  see comments in start.f90 for details
!
      use Cdata
!
      logical :: lstart
!
      if(ip==0) print*,'lstart=',lstart !(to keep compiler quiet)
    endsubroutine initialize_forcing
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
