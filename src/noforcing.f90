! $Id: noforcing.f90,v 1.8 2002-09-26 16:21:25 brandenb Exp $

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
           "$Id: noforcing.f90,v 1.8 2002-09-26 16:21:25 brandenb Exp $")
!
    endsubroutine register_forcing
!***********************************************************************
    subroutine forcing_run_hook()
!
!  initialize random number generator in processor-dependent fashion
!  see comments in start.f90 for details
!
      use Cdata
!
      if (lroot.and.ip<14) print*, 'initializing seed'
      call random_seed_wrapper(get=seed(1:nseed))
      seed(1) = -(10+iproc)    ! different random numbers on different CPUs
      call random_seed_wrapper(put=seed(1:nseed))
!
    endsubroutine forcing_run_hook
!***********************************************************************
    subroutine param_check_forcing
!
!  dummy routine
!
    endsubroutine param_check_forcing
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
